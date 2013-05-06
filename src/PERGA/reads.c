/*
 * reads.c
 *
 *  Created on: Dec 6, 2012
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"



/**
 * Build read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructReadset(readSet_t **readSet, char **readsFileNames, int readsFileNum, int reserveHashItemBlocksFlag)
{
	//printf("\n============= Begin to construct read set, please wait ... =============\n");

	if(readsFileFormatType==FILE_FORMAT_FASTA)
	{
		if(pairedMode==0)
		{
			if(constructReadsetBySEFasta(readSet, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the read set, Error!\n", __LINE__, __func__);
				return FAILED;
			}

			//outputReadseqInReadset("../output/reads.txt", *readSet);

		}else if(pairedMode==1)
		{
			if(constructReadsetByPEFastaSeparate(readSet, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the read set, Error!\n", __LINE__, __func__);
				return 1;
			}
		}else if(pairedMode==2)
		{
			if(constructReadsetByPEFastaInterleaved(readSet, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{
			printf("Invalid paired mode %d, error!\n", pairedMode);
			return FAILED;
		}
	}else if(readsFileFormatType==FILE_FORMAT_FASTQ)
	{
		if(pairedMode==0)
		{
			if(constructReadsetBySEFastq(readSet, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(pairedMode==1)
		{
			if(constructReadsetByPEFastqSeparate(readSet, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(pairedMode==2)
		{
			if(constructReadsetByPEFastqInterleaved(readSet, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else
		{
			printf("Invalid paired mode %d, error!\n", pairedMode);
			return FAILED;
		}
	}else
	{
		printf("Invalid file format that neither in fasta or fastq, error!\n");
		return FAILED;
	}

	// compute the maxReadLen
	if(computeMaxReadLenInReadset(*readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the maximal read length in read set, Error!\n", __LINE__, __func__);
		return FAILED;
	}

#if (DEBUG_PARA_PRINT==YES)
	printf("totalItemNumRead=%ld, totalItemNumReadseq=%ld, totalItemNumReadseqHashItem=%ld\n", (*readSet)->totalItemNumRead, (*readSet)->totalItemNumReadseq, (*readSet)->totalItemNumReadseqHashItem);
#endif

	if(reserveHashItemBlocksFlag==NO)
	{
		releaseHashItemReadset(*readSet);
	}

	//printf("============= End constructed read set. =============\n");

	return SUCCESSFUL;
}


/**
 * Build read set by single end files in fasta format.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructReadsetBySEFasta(readSet_t **readSet, char **readsFileNames, int readsFileNum)
{
	FILE* srcfp;
	int i, tmpFileID, percent;
	char ch, seq_data[MAX_READ_LEN_IN_BUF+1]; // read sequence
	int64_t tmpReadCount;

	struct timeval tpstart, tpend;
	double timeused_readset;
	gettimeofday(&tpstart, NULL);

	printf("Loading reads information ...\n");

	tmpReadCount = 0;
	percent = 0;

	if(initReadSet(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	pReadBlockTmp = (*readSet)->readBlockArr + (*readSet)->blocksNumRead - 1;
	pReadTmpDoing = pReadBlockTmp->readArr;

	pReadseqBlockTmp = (*readSet)->readseqBlockArr + (*readSet)->blocksNumReadseq - 1;
	pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;

	pReadseqHashItemBlockTmp = (*readSet)->readseqHashItemBlockArr + (*readSet)->blocksNumReadseqHashItem - 1;
	pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;


	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if(!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			return FAILED;
		}else
		{
			//printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		ch = fgetc(srcfp);
		while(!feof(srcfp))
		{
			while(ch!='\n') ch = fgetc(srcfp);

			i = 0;
			ch = fgetc(srcfp);
			while(ch!='>' && ch!=-1)
			{
				if(ch!='\n' && i<readLenCutOff)
					seq_data[i++] = ch;
				ch = fgetc(srcfp);
			}
			seq_data[i] = '\0';

			tmpReadCount ++;

			//add a read
			if(addReadToReadset(seq_data, NULL, i, *readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot add read to read set, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(tmpReadCount%1000000==0)
				printf("Sequences processed: %ld\n", tmpReadCount);
		}

		fclose(srcfp);
		srcfp = NULL;
	}

	printf("Sequences processed: %ld\n", tmpReadCount);

	// the last read block is empty, remove it
	if(pReadBlockTmp->itemNum==0)
	{
		free(pReadBlockTmp->readArr);
		(*readSet)->readBlockArr[(*readSet)->blocksNumRead-1].readArr = NULL;
		(*readSet)->blocksNumRead --;
	}

	// the last readseq block is empty, remove it
	if(pReadseqBlockTmp->rowsNum==0)
	{
		free(pReadseqBlockTmp->readseqArr);
		(*readSet)->readseqBlockArr[(*readSet)->blocksNumReadseq-1].readseqArr = NULL;
		(*readSet)->blocksNumReadseq --;
	}

	// the last readseq hash item block is empty, remove it
	if(pReadseqHashItemBlockTmp->itemNum==0)
	{
		free(pReadseqHashItemBlockTmp->readseqHashItemArr);
		(*readSet)->readseqHashItemBlockArr[(*readSet)->blocksNumReadseqHashItem-1].readseqHashItemArr = NULL;
		(*readSet)->blocksNumReadseqHashItem --;
	}

	gettimeofday(&tpend, NULL);
	timeused_readset = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Loading reads information used time: %.2f seconds.\n", timeused_readset);

	return SUCCESSFUL;
}

/**
 * Build read set by paired ends in separate fasta files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructReadsetByPEFastaSeparate(readSet_t **readSet, char **readsFileNames, int readsFileNum)
{
	readBuf_t *readBuf[2];
	FILE *fp1, *fp2;
	uint64_t i, j, tmpReadsNum[2], tmpFileID;
	int percent;
	int64_t tmpReadCount;

	struct timeval tpstart, tpend;
	double timeused_readset;
	gettimeofday(&tpstart, NULL);


	printf("Loading reads information ...\n");

	for(i=0; i<2; i++)
	{
		readBuf[i] = (readBuf_t*)calloc(MAX_READ_BUF_SIZE, sizeof(readBuf_t));
		if(readBuf[i]==NULL)
		{
			printf("line=%d, In %s(), can not allocate memory, Error.\n", __LINE__, __func__);
			return FAILED;
		}
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			readBuf[i][j].seq = (char *)malloc(sizeof(char)*(MAX_READ_LEN_IN_BUF+1));
			if(readBuf[i][j].seq==NULL)
			{
				printf("line=%d, In %s(), can not allocate memory, Error.\n", __LINE__, __func__);
				return FAILED;
			}
			readBuf[i][j].qual = NULL;
		}
	}

	// initialize read set
	if(initReadSet(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	pReadBlockTmp = (*readSet)->readBlockArr + (*readSet)->blocksNumRead - 1;
	pReadTmpDoing = pReadBlockTmp->readArr;

	pReadseqBlockTmp = (*readSet)->readseqBlockArr + (*readSet)->blocksNumReadseq - 1;
	pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;

	pReadseqHashItemBlockTmp = (*readSet)->readseqHashItemBlockArr + (*readSet)->blocksNumReadseqHashItem - 1;
	pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;

	percent = 0;
	tmpReadCount = 0;

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID+=2)
	{
		fp1 = fopen(readsFileNames[tmpFileID], "r");
		if(fp1==NULL)
		{
			printf("line=%d, In %s(), can not open file [ %s ], error.\n", __LINE__, __func__, readsFileNames[tmpFileID]);
			return FAILED;
		}else
		{
			//printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}
		fp2 = fopen(readsFileNames[tmpFileID+1], "r");
		if(fp2==NULL)
		{
			printf("line=%d, In %s(), can not open file [ %s ], error.\n", __LINE__, __func__, readsFileNames[tmpFileID+1]);
			return FAILED;
		}else
		{
			//printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID+1]);
		}

		while(1)
		{
			// check the end of files
			if(feof(fp1) && feof(fp2))
			{
				break;
			}else if(feof(fp1) || feof(fp2))
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// fill the reads to reads buffers
			if(fillReadsToBufFasta(fp1, readBuf[0], tmpReadsNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(fillReadsToBufFasta(fp2, readBuf[1], tmpReadsNum+1)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(tmpReadsNum[0]!=tmpReadsNum[1])
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			//...
			for(i=0; i<tmpReadsNum[0]; i++)
			{
				for(j=0; j<2; j++)
				{
					tmpReadCount ++;

					if(readBuf[j][i].len>readLenCutOff)
					{
						readBuf[j][i].seq[readLenCutOff] = '\0';
						readBuf[j][i].len = readLenCutOff;
					}

					//add a read
					if(addReadToReadset(readBuf[j][i].seq, NULL, readBuf[j][i].len, *readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read to read set, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(tmpReadCount%1000000==0)
						printf("Sequences processed: %ld\n", tmpReadCount);

				}
			}
		}

		fclose(fp1);
		fp1 = NULL;
		fclose(fp2);
		fp2 = NULL;
	}

	printf("Sequences processed: %ld\n", tmpReadCount);

	for(i=0; i<2; i++)
	{
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			free(readBuf[i][j].seq);
		}
		free(readBuf[i]);
	}

	// the last read block is empty, remove it
	if(pReadBlockTmp->itemNum==0)
	{
		free(pReadBlockTmp->readArr);
		(*readSet)->readBlockArr[(*readSet)->blocksNumRead-1].readArr = NULL;
		(*readSet)->blocksNumRead --;
	}

	// the last readseq block is empty, remove it
	if(pReadseqBlockTmp->rowsNum==0)
	{
		free(pReadseqBlockTmp->readseqArr);
		(*readSet)->readseqBlockArr[(*readSet)->blocksNumReadseq-1].readseqArr = NULL;
		(*readSet)->blocksNumReadseq --;
	}

	// the last readseq hash item block is empty, remove it
	if(pReadseqHashItemBlockTmp->itemNum==0)
	{
		free(pReadseqHashItemBlockTmp->readseqHashItemArr);
		(*readSet)->readseqHashItemBlockArr[(*readSet)->blocksNumReadseqHashItem-1].readseqHashItemArr = NULL;
		(*readSet)->blocksNumReadseqHashItem --;
	}


	gettimeofday(&tpend, NULL);
	timeused_readset = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Loading reads used time: %.2f seconds.\n", timeused_readset);

	return SUCCESSFUL;
}

/**
 * Build read set by paired ends in interleaved fasta files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructReadsetByPEFastaInterleaved(readSet_t **readSet, char **readsFileNames, int readsFileNum)
{
	FILE* srcfp;
	int i, tmpFileID, percent;
	char ch, seq_data[MAX_READ_LEN_IN_BUF+1]; // read sequence, read quality data which is encoded in ASCII
	int64_t tmpReadCount;

	struct timeval tpstart, tpend;
	double timeused_readset;
	gettimeofday(&tpstart, NULL);


	printf("Loading reads information ...\n");

	tmpReadCount = 0;
	percent = 0;

	// initialize read set
	if(initReadSet(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	pReadBlockTmp = (*readSet)->readBlockArr + (*readSet)->blocksNumRead - 1;
	pReadTmpDoing = pReadBlockTmp->readArr;

	pReadseqBlockTmp = (*readSet)->readseqBlockArr + (*readSet)->blocksNumReadseq - 1;
	pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;

	pReadseqHashItemBlockTmp = (*readSet)->readseqHashItemBlockArr + (*readSet)->blocksNumReadseqHashItem - 1;
	pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;


	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if(!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			return FAILED;
		}else
		{
			//printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		ch = fgetc(srcfp);
		while(!feof(srcfp))
		{
			while(ch!='\n') ch = fgetc(srcfp);

			i = 0;
			ch = fgetc(srcfp);
			while(ch!='>' && ch!=-1)
			{
				if(ch!='\n' && i<readLenCutOff)
					seq_data[i++] = ch;
				ch = fgetc(srcfp);
			}
			seq_data[i] = '\0';

			tmpReadCount ++;

			//add a read
			if(addReadToReadset(seq_data, NULL, i, *readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot add read to read set, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(tmpReadCount%1000000==0)
				printf("Sequences processed: %ld\n", tmpReadCount);

		}

		fclose(srcfp);
		srcfp = NULL;
	}

	printf("Sequences processed: %ld\n", tmpReadCount);

	// the last read block is empty, remove it
	if(pReadBlockTmp->itemNum==0)
	{
		free(pReadBlockTmp->readArr);
		(*readSet)->readBlockArr[(*readSet)->blocksNumRead-1].readArr = NULL;
		(*readSet)->blocksNumRead --;
	}

	// the last readseq block is empty, remove it
	if(pReadseqBlockTmp->rowsNum==0)
	{
		free(pReadseqBlockTmp->readseqArr);
		(*readSet)->readseqBlockArr[(*readSet)->blocksNumReadseq-1].readseqArr = NULL;
		(*readSet)->blocksNumReadseq --;
	}

	// the last readseq hash item block is empty, remove it
	if(pReadseqHashItemBlockTmp->itemNum==0)
	{
		free(pReadseqHashItemBlockTmp->readseqHashItemArr);
		(*readSet)->readseqHashItemBlockArr[(*readSet)->blocksNumReadseqHashItem-1].readseqHashItemArr = NULL;
		(*readSet)->blocksNumReadseqHashItem --;
	}

	gettimeofday(&tpend, NULL);
	timeused_readset = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Loading reads used time: %.2f seconds.\n", timeused_readset);

	return SUCCESSFUL;
}

/**
 * Build read set by single ends in fastq files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructReadsetBySEFastq(readSet_t **readSet, char **readsFileNames, int readsFileNum)
{
	FILE* srcfp;
	int i, tmpLen, line_index, tmpFileID, percent;
	char ch, seq_data[MAX_READ_LEN_IN_BUF+1], qual_data[MAX_READ_LEN_IN_BUF+1]; // read sequence, read quality data which is encoded in ASCII
	int64_t tmpReadCount;

	struct timeval tpstart, tpend;
	double timeused_readset;
	gettimeofday(&tpstart, NULL);


	printf("Loading reads information ...\n");

	tmpReadCount = 0;
	percent = 0;

	// initialize read set
	if(initReadSet(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	pReadBlockTmp = (*readSet)->readBlockArr + (*readSet)->blocksNumRead - 1;
	pReadTmpDoing = pReadBlockTmp->readArr;

	pReadseqBlockTmp = (*readSet)->readseqBlockArr + (*readSet)->blocksNumReadseq - 1;
	pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;

	pReadseqHashItemBlockTmp = (*readSet)->readseqHashItemBlockArr + (*readSet)->blocksNumReadseqHashItem - 1;
	pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;


	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if(!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			return FAILED;
		}else
		{
			//printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		line_index = 0;
		ch = fgetc(srcfp);
		while(!feof(srcfp))
		{
			if(line_index==0)  //the sequence name line
			{
				ch = fgetc(srcfp);
				while(ch!='\n' && ch!=-1)
				{
					ch = fgetc(srcfp);
				}
			}else if(line_index==1)  //the sequence line
			{
				tmpLen = 0;
				ch = fgetc(srcfp);
				while(ch!='\n' && ch!=-1)
				{
					if(tmpLen<readLenCutOff)
						seq_data[tmpLen++] = ch;
					ch = fgetc(srcfp);
				}
				seq_data[tmpLen] = '\0';
			}else if(line_index==2)  //the sequence name line
			{
				ch = fgetc(srcfp);
				while(ch!='\n' && ch!=-1)
				{
					ch = fgetc(srcfp);
				}
			}else
			{
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n' && ch!=-1)
				{
					if(i<readLenCutOff)
						qual_data[i++] = ch;
					ch = fgetc(srcfp);
				}
				qual_data[i] = '\0';
			}
			line_index++;

			if(line_index==4)  //the sequence is read finished, construct the read
			{
				tmpReadCount ++;

				//add a read
				if(addReadToReadset(seq_data, qual_data, tmpLen, *readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read to read set, error!\n", __LINE__, __func__);
					return FAILED;
				}

				if(tmpReadCount%1000000==0)
					printf("Sequences processed: %ld\n", tmpReadCount);

				line_index = 0;
			}
		}

		fclose(srcfp);
		srcfp = NULL;
	}

	printf("Sequences processed: %ld\n", tmpReadCount);

	// the last read block is empty, remove it
	if(pReadBlockTmp->itemNum==0)
	{
		free(pReadBlockTmp->readArr);
		(*readSet)->readBlockArr[(*readSet)->blocksNumRead-1].readArr = NULL;
		(*readSet)->blocksNumRead --;
	}

	// the last readseq block is empty, remove it
	if(pReadseqBlockTmp->rowsNum==0)
	{
		free(pReadseqBlockTmp->readseqArr);
		(*readSet)->readseqBlockArr[(*readSet)->blocksNumReadseq-1].readseqArr = NULL;
		(*readSet)->blocksNumReadseq --;
	}

	// the last readseq hash item block is empty, remove it
	if(pReadseqHashItemBlockTmp->itemNum==0)
	{
		free(pReadseqHashItemBlockTmp->readseqHashItemArr);
		(*readSet)->readseqHashItemBlockArr[(*readSet)->blocksNumReadseqHashItem-1].readseqHashItemArr = NULL;
		(*readSet)->blocksNumReadseqHashItem --;
	}


	gettimeofday(&tpend, NULL);
	timeused_readset = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Loading reads used time: %.2f seconds.\n", timeused_readset);

	return SUCCESSFUL;
}

/**
 * Build read set by paired ends in separated fastq files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructReadsetByPEFastqSeparate(readSet_t **readSet, char **readsFileNames, int readsFileNum)
{
	readBuf_t *readBuf[2];
	FILE *fp1, *fp2;
	uint64_t i, j, tmpReadsNum[2], tmpFileID;
	int percent;
	int64_t tmpReadCount;

	struct timeval tpstart, tpend;
	double timeused_readset;
	gettimeofday(&tpstart, NULL);


	printf("Loading reads information ...\n");

	for(i=0; i<2; i++)
	{
		readBuf[i] = (readBuf_t*)calloc(MAX_READ_BUF_SIZE, sizeof(readBuf_t));
		if(readBuf[i]==NULL)
		{
			printf("line=%d, In %s(), can not allocate memory, Error.\n", __LINE__,  __func__);
			return FAILED;
		}
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			readBuf[i][j].seq = (char *)malloc(sizeof(char)*(MAX_READ_LEN_IN_BUF+1));
			if(readBuf[i][j].seq==NULL)
			{
				printf("line=%d, In %s(), can not allocate memory, Error.\n", __LINE__,  __func__);
				return FAILED;
			}
			readBuf[i][j].qual = (char *)malloc(sizeof(char)*(MAX_READ_LEN_IN_BUF+1));
			if(readBuf[i][j].qual==NULL)
			{
				printf("line=%d, In %s(), can not allocate memory, Error.\n", __LINE__,  __func__);
				return FAILED;
			}
#if(DRAW_CURVE_FLAG==YES)
			readBuf[i][j].headname = (char *)malloc(sizeof(char)*(MAX_READ_LEN_IN_BUF+1));
			if(readBuf[i][j].qual==NULL)
			{
				printf("line=%d, In %s(), can not allocate memory, Error.\n", __LINE__,  __func__);
				return FAILED;
			}
#endif
		}
	}

	// initialize read set
	if(initReadSet(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	pReadBlockTmp = (*readSet)->readBlockArr + (*readSet)->blocksNumRead - 1;
	pReadTmpDoing = pReadBlockTmp->readArr;

	pReadseqBlockTmp = (*readSet)->readseqBlockArr + (*readSet)->blocksNumReadseq - 1;
	pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;

	pReadseqHashItemBlockTmp = (*readSet)->readseqHashItemBlockArr + (*readSet)->blocksNumReadseqHashItem - 1;
	pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;


	percent = 0;
	tmpReadCount = 0;

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID+=2)
	{
		fp1 = fopen(readsFileNames[tmpFileID], "r");
		if(fp1==NULL)
		{
			printf("In %s(), can not open file [ %s ], error.\n", __func__, readsFileNames[tmpFileID]);
			return FAILED;
		}else
		{
			//printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}
		fp2 = fopen(readsFileNames[tmpFileID+1], "r");
		if(fp2==NULL)
		{
			printf("In %s(), can not open file [ %s ], error.\n", __func__, readsFileNames[tmpFileID+1]);
			return FAILED;
		}else
		{
			//printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID+1]);
		}

		while(1)
		{
			// check the end of files
			if(feof(fp1) && feof(fp2))
			{
				break;
			}else if(feof(fp1) || feof(fp2))
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			// fille the reads to reads buffers
			if(fillReadsToBuf(fp1, readBuf[0], tmpReadsNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(fillReadsToBuf(fp2, readBuf[1], tmpReadsNum+1)==FAILED)
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}
			if(tmpReadsNum[0]!=tmpReadsNum[1])
			{
				printf("line=%d, In %s(), cannot fill the read buffer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			//...
			for(i=0; i<tmpReadsNum[0]; i++)
			{
				for(j=0; j<2; j++)
				{
					tmpReadCount ++;

					if(readBuf[j][i].len>readLenCutOff)
					{
						readBuf[j][i].seq[readLenCutOff] = '\0';
						readBuf[j][i].qual[readLenCutOff] = '\0';
						readBuf[j][i].len = readLenCutOff;
					}

					//add a read
					if(addReadToReadset(readBuf[j][i].seq, readBuf[j][i].qual, readBuf[j][i].len, *readSet)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read to read set, error!\n", __LINE__, __func__);
						return FAILED;
					}

					if(tmpReadCount%1000000==0)
						printf("Sequences processed: %ld\n", tmpReadCount);

				}
			}
		}

		fclose(fp1);
		fp1 = NULL;
		fclose(fp2);
		fp2 = NULL;
	}

	printf("Sequences processed: %ld\n", tmpReadCount);

	for(i=0; i<2; i++)
	{
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			free(readBuf[i][j].seq);
			free(readBuf[i][j].qual);
#if(DRAW_CURVE_FLAG==YES)
			free(readBuf[i][j].headname);
#endif
		}
		free(readBuf[i]);
	}

	// the last read block is empty, remove it
	if(pReadBlockTmp->itemNum==0)
	{
		free(pReadBlockTmp->readArr);
		(*readSet)->readBlockArr[(*readSet)->blocksNumRead-1].readArr = NULL;
		(*readSet)->blocksNumRead --;
	}

	// the last readseq block is empty, remove it
	if(pReadseqBlockTmp->rowsNum==0)
	{
		free(pReadseqBlockTmp->readseqArr);
		(*readSet)->readseqBlockArr[(*readSet)->blocksNumReadseq-1].readseqArr = NULL;
		(*readSet)->blocksNumReadseq --;
	}

	// the last readseq hash item block is empty, remove it
	if(pReadseqHashItemBlockTmp->itemNum==0)
	{
		free(pReadseqHashItemBlockTmp->readseqHashItemArr);
		(*readSet)->readseqHashItemBlockArr[(*readSet)->blocksNumReadseqHashItem-1].readseqHashItemArr = NULL;
		(*readSet)->blocksNumReadseqHashItem --;
	}

	gettimeofday(&tpend, NULL);
	timeused_readset = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Loading reads used time: %.2f seconds.\n", timeused_readset);

	return SUCCESSFUL;
}


/**
 * Build read set by paired ends in interleaved fastq files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructReadsetByPEFastqInterleaved(readSet_t **readSet, char **readsFileNames, int readsFileNum)
{
	FILE* srcfp;
	int i, tmpLen, line_index, tmpFileID, percent, headlen;
	char ch, seq_data[MAX_READ_LEN_IN_BUF+1], qual_data[MAX_READ_LEN_IN_BUF+1]; // read sequence, read quality data which is encoded in ASCII
	char headname[MAX_READ_LEN_IN_BUF+1];
	int64_t tmpReadCount;

	struct timeval tpstart, tpend;
	double timeused_readset;
	gettimeofday(&tpstart, NULL);


	printf("Loading reads information ...\n");


	tmpReadCount = 0;
	percent = 0;

	// initialize read set
	if(initReadSet(readSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	pReadBlockTmp = (*readSet)->readBlockArr + (*readSet)->blocksNumRead - 1;
	pReadTmpDoing = pReadBlockTmp->readArr;

	pReadseqBlockTmp = (*readSet)->readseqBlockArr + (*readSet)->blocksNumReadseq - 1;
	pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;

	pReadseqHashItemBlockTmp = (*readSet)->readseqHashItemBlockArr + (*readSet)->blocksNumReadseqHashItem - 1;
	pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;


	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if(!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			return FAILED;
		}else
		{
			//printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		line_index = 0;
		ch = fgetc(srcfp);
		while(!feof(srcfp))
		{
			if(line_index==0)  //the sequence name line
			{
				headlen = 0;
				ch = fgetc(srcfp);
				while(ch!='\n' && ch!=-1)
				{
#if(DRAW_CURVE_FLAG==YES)
					headname[ headlen ++ ] = ch;
#endif
					ch = fgetc(srcfp);
				}
			}else if(line_index==1)  //the sequence line
			{
				tmpLen = 0;
				ch = fgetc(srcfp);
				while(ch!='\n' && ch!=-1)
				{
					if(tmpLen<readLenCutOff)
						seq_data[tmpLen++] = ch;
					ch = fgetc(srcfp);
				}
				seq_data[tmpLen] = '\0';
			}else if(line_index==2)  //the sequence name line
			{
				ch = fgetc(srcfp);
				while(ch!='\n' && ch!=-1)
				{
					ch = fgetc(srcfp);
				}
			}else
			{
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n' && ch!=-1)
				{
					if(i<readLenCutOff)
						qual_data[i++] = ch;
					ch = fgetc(srcfp);
				}
				qual_data[i] = '\0';
			}
			line_index++;

			if(line_index==4)  //the sequence is read finished, construct the read
			{
				tmpReadCount ++;

#if(DRAW_CURVE_FLAG==YES)
				//add a read
				if(addReadToReadsetWithRefpos(seq_data, qual_data, tmpLen, headname, headlen, *readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read to read set, error!\n", __LINE__, __func__);
					return FAILED;
				}
#else
				//add a read
				if(addReadToReadset(seq_data, qual_data, tmpLen, *readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read to read set, error!\n", __LINE__, __func__);
					return FAILED;
				}
#endif

				if(tmpReadCount%1000000==0)
					printf("Sequences processed: %ld\n", tmpReadCount);

				line_index = 0;
			}
		}

		fclose(srcfp);
		srcfp = NULL;
	}

	printf("Sequences processed: %ld\n", tmpReadCount);

	// the last read block is empty, remove it
	if(pReadBlockTmp->itemNum==0)
	{
		free(pReadBlockTmp->readArr);
		(*readSet)->readBlockArr[(*readSet)->blocksNumRead-1].readArr = NULL;
		(*readSet)->blocksNumRead --;
	}

	// the last readseq block is empty, remove it
	if(pReadseqBlockTmp->rowsNum==0)
	{
		free(pReadseqBlockTmp->readseqArr);
		(*readSet)->readseqBlockArr[(*readSet)->blocksNumReadseq-1].readseqArr = NULL;
		(*readSet)->blocksNumReadseq --;
	}

	// the last readseq hash item block is empty, remove it
	if(pReadseqHashItemBlockTmp->itemNum==0)
	{
		free(pReadseqHashItemBlockTmp->readseqHashItemArr);
		(*readSet)->readseqHashItemBlockArr[(*readSet)->blocksNumReadseqHashItem-1].readseqHashItemArr = NULL;
		(*readSet)->blocksNumReadseqHashItem --;
	}

	gettimeofday(&tpend, NULL);
	timeused_readset = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Loading reads used time: %.2f seconds.\n", timeused_readset);

	return SUCCESSFUL;
}


/**
 * Add read to read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadToReadset(char *seq, char *qual_data, int32_t seqLen, readSet_t *readSet)
{
	int32_t entriesNum, validFlag, unknownBaseNum;
	uint64_t hashcode;
	readseqHashItem_t *pReadseqHashItem;

	if(seqLen>=kmerSize && getRatioBaseA(seq)<ARTIFACTS_BASE_A_THRESHOLD)
	{
		unknownBaseNum = getUnknownBaseNum(seq);
		if(unknownBaseNum==0)
			validFlag = YES;
		else
		{
			validFlag = NO;
			if(replaceUnknownBasesInReads(seq, unknownBaseNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot replace unknown bases in reads, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}else
		validFlag = NO;

	// check the quality values
	if(validFlag==YES && qual_data!=NULL)
	{
		if(calcAverQual5End(qual_data, seqLen)>=AVERAGE_QUAL_THRESHOLD_5End && singleQualSatisfied(qual_data)==YES && calcAverQual3End(qual_data, seqLen)>=AVERAGE_QUAL_THRESHOLD_3End)
			validFlag = YES;
		else
			validFlag = NO;
	}

	entriesNum = ((seqLen - 1) >> 5) + 1;

	if(validFlag==YES)
	{
		// set the valid pointer
		if(pReadseqBlockTmp->rowsNum + entriesNum > readSet->maxEntryNumReadseqBlock)
		{
			// add new readseq block
			if(addNewBlockReadseq(readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new readseq block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
			pReadseqBlockTmp = readSet->readseqBlockArr + readSet->blocksNumReadseq - 1;
			pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;
		}

		// generate the readseq
		if(generateReadseqInt(pReadseqTmpDoing, seq, seqLen, entriesNum)==FAILED)
		{
			printf("line=%d, In %s(), unknownBaseNum=%d, cannot generate the readseq integer sequence, error!\n", __LINE__, __func__, unknownBaseNum);
			return FAILED;
		}

		// get the hashcode and add the readseq
		hashcode = readseqHashInt(pReadseqTmpDoing, seqLen, entriesNum);
		pReadseqHashItem = getReadseqHashItemByHash(hashcode, pReadseqTmpDoing, seqLen, entriesNum, readSet);
		if(pReadseqHashItem)
		{ // set the readseq to the readseq blocks
			pReadTmpDoing->readseqBlockID = pReadseqHashItem->readseqBlockID;
			pReadTmpDoing->rowReadseqInBlock = pReadseqHashItem->rowReadseqInBlock;
		}
		else
		{ // add the new readseq item into readseq blocks and readseq hash item blocks
			// read block operations
			pReadTmpDoing->readseqBlockID = pReadseqBlockTmp->blockID;
			pReadTmpDoing->rowReadseqInBlock = pReadseqBlockTmp->rowsNum;

			// readseq hash item block operations
			// add new item in readseq hash item block
			pReadseqHashItemTmpDoing->readseqBlockID = pReadseqBlockTmp->blockID;
			pReadseqHashItemTmpDoing->rowReadseqInBlock = pReadseqBlockTmp->rowsNum;
			pReadseqHashItemTmpDoing->seqlen = seqLen;
			pReadseqHashItemTmpDoing->nextHashItemBlockID = readSet->readseqHashtable[hashcode].hashItemBlockID;
			pReadseqHashItemTmpDoing->nextItemRowHashItemBlock = readSet->readseqHashtable[hashcode].itemRowHashItemBlock;
			readSet->readseqHashtable[hashcode].hashItemBlockID = pReadseqHashItemBlockTmp->blockID;
			readSet->readseqHashtable[hashcode].itemRowHashItemBlock = pReadseqHashItemBlockTmp->itemNum;

			readSet->totalItemNumReadseqHashItem ++;
			pReadseqHashItemBlockTmp->itemNum ++;
			pReadseqHashItemTmpDoing ++;
			if(pReadseqHashItemBlockTmp->itemNum >= readSet->maxItemNumPerReadseqHashItemBlock)
			{
				// add new readseq hash item block
				if(addNewBlockReadseqHashItem(readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new readseq hash item block, Error!\n", __LINE__, __func__);
					return FAILED;
				}
				pReadseqHashItemBlockTmp = readSet->readseqHashItemBlockArr + readSet->blocksNumReadseqHashItem - 1;
				pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;
			}

			// readseq block operations
			pReadseqTmpDoing += entriesNum;
			readSet->totalItemNumReadseq ++;
			pReadseqBlockTmp->rowsNum += entriesNum;

			// set the valid pointer
			if(pReadseqBlockTmp->rowsNum >= readSet->maxEntryNumReadseqBlock)
			{
				// add new readseq block
				if(addNewBlockReadseq(readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new readseq block, Error!\n", __LINE__, __func__);
					return FAILED;
				}
				pReadseqBlockTmp = readSet->readseqBlockArr + readSet->blocksNumReadseq - 1;
				pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;
			}
		}
	}


	// add the read info
	pReadTmpDoing->seqlen = seqLen;
	pReadTmpDoing->nBaseNum = unknownBaseNum;
	pReadTmpDoing->validFlag = validFlag;
	pReadTmpDoing->modified = NO;
	pReadBlockTmp->itemNum ++;
	readSet->totalItemNumRead ++;
	if(validFlag==YES)
		readSet->totalValidItemNumRead ++;
	pReadTmpDoing++;

	if(pReadBlockTmp->itemNum >= readSet->maxItemNumPerReadBlock)
	{
		// add new read block
		if(addNewBlockRead(readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot add new read block, Error!\n", __LINE__, __func__);
			return FAILED;
		}
		pReadBlockTmp = readSet->readBlockArr + readSet->blocksNumRead - 1;
		pReadTmpDoing = pReadBlockTmp->readArr;
	}


	return SUCCESSFUL;
}

/**
 * Add read to read set with ref pos.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addReadToReadsetWithRefpos(char *seq, char *qual_data, int32_t seqLen, char *headname, int32_t headlen, readSet_t *readSet)
{
	int32_t entriesNum, validFlag, unknownBaseNum;
	int32_t strandPosArray[2];
	uint64_t hashcode;
	readseqHashItem_t *pReadseqHashItem;

	if(seqLen>=kmerSize && getRatioBaseA(seq)<ARTIFACTS_BASE_A_THRESHOLD)
	{
		unknownBaseNum = getUnknownBaseNum(seq);
		if(unknownBaseNum==0)
			validFlag = YES;
		else
		{
			validFlag = NO;
			if(replaceUnknownBasesInReads(seq, unknownBaseNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot replace unknown bases in reads, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}else
		validFlag = NO;

	// check the quality values
	if(validFlag==YES && qual_data!=NULL)
	{
		if(calcAverQual5End(qual_data, seqLen)>=AVERAGE_QUAL_THRESHOLD_5End && singleQualSatisfied(qual_data)==YES && calcAverQual3End(qual_data, seqLen)>=AVERAGE_QUAL_THRESHOLD_3End)
			validFlag = YES;
		else
			validFlag = NO;
	}

	if(extractRefPosFromHeadName(strandPosArray, headname, headlen)==FAILED)
	{
		printf("line=%d, In %s(), cannot extract ref pos info from read head name, error!\n", __LINE__, __func__);
		return FAILED;
	}

	entriesNum = ((seqLen - 1) >> 5) + 1;

	if(validFlag==YES)
	{
		// set the valid pointer
		if(pReadseqBlockTmp->rowsNum + entriesNum > readSet->maxEntryNumReadseqBlock)
		{
			// add new readseq block
			if(addNewBlockReadseq(readSet)==FAILED)
			{
				printf("line=%d, In %s(), cannot add new readseq block, Error!\n", __LINE__, __func__);
				return FAILED;
			}
			pReadseqBlockTmp = readSet->readseqBlockArr + readSet->blocksNumReadseq - 1;
			pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;
		}

		// generate the readseq
		if(generateReadseqInt(pReadseqTmpDoing, seq, seqLen, entriesNum)==FAILED)
		{
			printf("line=%d, In %s(), unknownBaseNum=%d, cannot generate the readseq integer sequence, error!\n", __LINE__, __func__, unknownBaseNum);
			return FAILED;
		}

		// get the hashcode and add the readseq
		hashcode = readseqHashInt(pReadseqTmpDoing, seqLen, entriesNum);
		pReadseqHashItem = getReadseqHashItemByHash(hashcode, pReadseqTmpDoing, seqLen, entriesNum, readSet);
		if(pReadseqHashItem)
		{ // set the readseq to the readseq blocks
			pReadTmpDoing->readseqBlockID = pReadseqHashItem->readseqBlockID;
			pReadTmpDoing->rowReadseqInBlock = pReadseqHashItem->rowReadseqInBlock;
		}
		else
		{ // add the new readseq item into readseq blocks and readseq hash item blocks
			// read block operations
			pReadTmpDoing->readseqBlockID = pReadseqBlockTmp->blockID;
			pReadTmpDoing->rowReadseqInBlock = pReadseqBlockTmp->rowsNum;

			// readseq hash item block operations
			// add new item in readseq hash item block
			pReadseqHashItemTmpDoing->readseqBlockID = pReadseqBlockTmp->blockID;
			pReadseqHashItemTmpDoing->rowReadseqInBlock = pReadseqBlockTmp->rowsNum;
			pReadseqHashItemTmpDoing->seqlen = seqLen;
			pReadseqHashItemTmpDoing->nextHashItemBlockID = readSet->readseqHashtable[hashcode].hashItemBlockID;
			pReadseqHashItemTmpDoing->nextItemRowHashItemBlock = readSet->readseqHashtable[hashcode].itemRowHashItemBlock;
			readSet->readseqHashtable[hashcode].hashItemBlockID = pReadseqHashItemBlockTmp->blockID;
			readSet->readseqHashtable[hashcode].itemRowHashItemBlock = pReadseqHashItemBlockTmp->itemNum;

			readSet->totalItemNumReadseqHashItem ++;
			pReadseqHashItemBlockTmp->itemNum ++;
			pReadseqHashItemTmpDoing ++;
			if(pReadseqHashItemBlockTmp->itemNum >= readSet->maxItemNumPerReadseqHashItemBlock)
			{
				// add new readseq hash item block
				if(addNewBlockReadseqHashItem(readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new readseq hash item block, Error!\n", __LINE__, __func__);
					return FAILED;
				}
				pReadseqHashItemBlockTmp = readSet->readseqHashItemBlockArr + readSet->blocksNumReadseqHashItem - 1;
				pReadseqHashItemTmpDoing = pReadseqHashItemBlockTmp->readseqHashItemArr;
			}

			// readseq block operations
			pReadseqTmpDoing += entriesNum;
			readSet->totalItemNumReadseq ++;
			pReadseqBlockTmp->rowsNum += entriesNum;

			// set the valid pointer
			if(pReadseqBlockTmp->rowsNum >= readSet->maxEntryNumReadseqBlock)
			{
				// add new readseq block
				if(addNewBlockReadseq(readSet)==FAILED)
				{
					printf("line=%d, In %s(), cannot add new readseq block, Error!\n", __LINE__, __func__);
					return FAILED;
				}
				pReadseqBlockTmp = readSet->readseqBlockArr + readSet->blocksNumReadseq - 1;
				pReadseqTmpDoing = pReadseqBlockTmp->readseqArr;
			}
		}
	}


	// add the read info
	pReadTmpDoing->seqlen = seqLen;
	pReadTmpDoing->nBaseNum = unknownBaseNum;
	pReadTmpDoing->validFlag = validFlag;
	pReadTmpDoing->modified = NO;
#if(DRAW_CURVE_FLAG==YES)
	pReadTmpDoing->refStrand = strandPosArray[0];
	pReadTmpDoing->refPos = strandPosArray[1];
#endif
	pReadBlockTmp->itemNum ++;
	readSet->totalItemNumRead ++;
	if(validFlag==YES)
		readSet->totalValidItemNumRead ++;
	pReadTmpDoing++;

	if(pReadBlockTmp->itemNum >= readSet->maxItemNumPerReadBlock)
	{
		// add new read block
		if(addNewBlockRead(readSet)==FAILED)
		{
			printf("line=%d, In %s(), cannot add new read block, Error!\n", __LINE__, __func__);
			return FAILED;
		}
		pReadBlockTmp = readSet->readBlockArr + readSet->blocksNumRead - 1;
		pReadTmpDoing = pReadBlockTmp->readArr;
	}

	return SUCCESSFUL;
}


//==========================================================================================================
//==========================================================================================================
//==========================================================================================================

/**
 * Initialize reads set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initReadSet(readSet_t **pReadSet)
{
	//allocate readSet node
	*pReadSet = (readSet_t *) calloc (1, sizeof(readSet_t));
	if((*pReadSet)==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize read blocks
	if(initReadBlockInReadset(*pReadSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize read blocks, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize readseq blocks
	if(initReadseqBlockInReadset(*pReadSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize read blocks, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize readseq hash table
	if(initReadseqHashtableInReadset(*pReadSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize readseq hash table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Initialize read blocks.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initReadBlockInReadset(readSet_t *pReadSet)
{
	pReadSet->bytesPerRead = sizeof(read_t);
	pReadSet->blocksNumRead = 0;
	pReadSet->maxBlocksNumRead = MAX_BLOCKS_NUM_READ;
	pReadSet->maxReadLen = 0;
	pReadSet->totalItemNumRead = 0;
	pReadSet->totalValidItemNumRead = 0;
	pReadSet->maxItemNumPerReadBlock = BLOCK_SIZE_PER_READ / pReadSet->bytesPerRead;

	// allocate read blocks
	pReadSet->readBlockArr = (readBlock_t *) malloc(pReadSet->maxBlocksNumRead * sizeof(readBlock_t));
	if( pReadSet->readBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate the readBlockArr for the k-mer hash table, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	pReadSet->blocksNumRead = 0;

	// add new read block
	if(addNewBlockRead(pReadSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot allocate memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}


/**
 * Add new read block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockRead(readSet_t *pReadSet)
{
	if(pReadSet->blocksNumRead>=pReadSet->maxBlocksNumRead)
	{
		pReadSet->readBlockArr = (readBlock_t *) realloc (pReadSet->readBlockArr, pReadSet->maxBlocksNumRead*2*sizeof(readBlock_t));
		if(pReadSet->readBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pReadSet->maxBlocksNumRead *= 2;
	}

	pReadSet->readBlockArr[pReadSet->blocksNumRead].readArr = (read_t *) calloc (pReadSet->maxItemNumPerReadBlock, pReadSet->bytesPerRead);
	if(pReadSet->readBlockArr[pReadSet->blocksNumRead].readArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	pReadSet->readBlockArr[pReadSet->blocksNumRead].blockID = pReadSet->blocksNumRead + 1;
	pReadSet->readBlockArr[pReadSet->blocksNumRead].itemNum = 0;

	pReadSet->blocksNumRead ++;

	return SUCCESSFUL;
}

/**
 * Initialize readseq blocks.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initReadseqBlockInReadset(readSet_t *pReadSet)
{
	pReadSet->blocksNumReadseq = 0;
	pReadSet->maxBlocksNumReadseq = MAX_BLOCKS_NUM_READ_SEQ;
	pReadSet->totalItemNumReadseq = 0;
	pReadSet->bytesPerEntryReadseq = sizeof(uint64_t);
	pReadSet->maxEntryNumReadseqBlock = BLOCK_SIZE_PER_READ_SEQ / pReadSet->bytesPerEntryReadseq;

	// allocate readseq blocks
	pReadSet->readseqBlockArr = (readseqBlock_t *) malloc(pReadSet->maxBlocksNumReadseq * sizeof(readseqBlock_t));
	if( pReadSet->readseqBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate the readseqBlockArr for the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}
	pReadSet->blocksNumReadseq = 0;

	// add new readseq block
	if(addNewBlockReadseq(pReadSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot add new readseq block for the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Add new readseq block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockReadseq(readSet_t *pReadSet)
{
	if(pReadSet->blocksNumReadseq>=pReadSet->maxBlocksNumReadseq)
	{
		pReadSet->readseqBlockArr = (readseqBlock_t *) realloc (pReadSet->readseqBlockArr, pReadSet->maxBlocksNumReadseq*2*sizeof(readseqBlock_t));
		if(pReadSet->readseqBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pReadSet->maxBlocksNumReadseq *= 2;
	}

	pReadSet->readseqBlockArr[pReadSet->blocksNumReadseq].readseqArr = (uint64_t *) calloc (pReadSet->maxEntryNumReadseqBlock, pReadSet->bytesPerEntryReadseq);
	if(pReadSet->readseqBlockArr[pReadSet->blocksNumReadseq].readseqArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	pReadSet->readseqBlockArr[pReadSet->blocksNumReadseq].blockID = pReadSet->blocksNumReadseq + 1;
	pReadSet->readseqBlockArr[pReadSet->blocksNumReadseq].rowsNum = 0;

	pReadSet->blocksNumReadseq ++;

	return SUCCESSFUL;
}

/**
 * Initialize k-mer hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initReadseqHashtableInReadset(readSet_t *pReadSet)
{
	pReadSet->hashTableSizeReadseq = hashTableSizeReadseq;

	pReadSet->readseqHashtable = (readseqHashBucket_t *) calloc(pReadSet->hashTableSizeReadseq, sizeof(readseqHashBucket_t));
	if(pReadSet->readseqHashtable==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// initialize readseq hash item block
	if(initReadseqHashItemBlockInGraph(pReadSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the readseq hash item blocks in read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Initialize readseq hash item blocks.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initReadseqHashItemBlockInGraph(readSet_t *pReadSet)
{
	pReadSet->blocksNumReadseqHashItem = 0;
	pReadSet->maxBlocksNumReadseqHashItem = MAX_BLOCKS_NUM_READ_SEQ_HASH;
	pReadSet->totalItemNumReadseqHashItem = 0;
	pReadSet->bytesPerReadseqHashItem = sizeof(readseqHashItem_t);
	pReadSet->maxItemNumPerReadseqHashItemBlock = BLOCK_SIZE_PER_READ_SEQ_HASH / pReadSet->bytesPerReadseqHashItem;

	// allocate readseq hash item blocks
	pReadSet->readseqHashItemBlockArr = (readseqHashItemBlock_t *) malloc(pReadSet->maxBlocksNumReadseqHashItem * sizeof(readseqHashItemBlock_t));
	if( pReadSet->readseqHashItemBlockArr == NULL )
	{
		printf("line=%d, In %s(), cannot allocate the readseq hash item blocks for the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}
	pReadSet->blocksNumReadseqHashItem = 0;

	// add new readseq hash item block
	if(addNewBlockReadseqHashItem(pReadSet)==FAILED)
	{
		printf("line=%d, In %s(), cannot add new readseq hash item block for the read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Add new readseq hash item block.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addNewBlockReadseqHashItem(readSet_t *pReadSet)
{
	if(pReadSet->blocksNumReadseqHashItem>=pReadSet->maxBlocksNumReadseqHashItem)
	{
		pReadSet->readseqHashItemBlockArr = (readseqHashItemBlock_t *) realloc (pReadSet->readseqHashItemBlockArr, pReadSet->maxBlocksNumReadseqHashItem*2*sizeof(readseqHashItemBlock_t));
		if(pReadSet->readseqHashItemBlockArr==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		pReadSet->maxBlocksNumReadseqHashItem *= 2;
	}

	pReadSet->readseqHashItemBlockArr[pReadSet->blocksNumReadseqHashItem].readseqHashItemArr = (readseqHashItem_t *) calloc (pReadSet->maxItemNumPerReadseqHashItemBlock, pReadSet->bytesPerReadseqHashItem);
	if(pReadSet->readseqHashItemBlockArr[pReadSet->blocksNumReadseqHashItem].readseqHashItemArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, Error!\n", __LINE__, __func__);
		return FAILED;
	}
	pReadSet->readseqHashItemBlockArr[pReadSet->blocksNumReadseqHashItem].blockID = pReadSet->blocksNumReadseqHashItem + 1;
	pReadSet->readseqHashItemBlockArr[pReadSet->blocksNumReadseqHashItem].itemNum = 0;

	pReadSet->blocksNumReadseqHashItem ++;

	return SUCCESSFUL;
}

/**
 * Release read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short releaseReadset(readSet_t *pReadSet)
{
	int32_t i;

	// release read blocks
	for(i=0; i<pReadSet->blocksNumRead; i++)
		free(pReadSet->readBlockArr[i].readArr);
	free(pReadSet->readBlockArr);

	// release readseq blocks
	for(i=0; i<pReadSet->blocksNumReadseq; i++)
		free(pReadSet->readseqBlockArr[i].readseqArr);
	free(pReadSet->readseqBlockArr);

	// release readseq hash blocks
	releaseHashItemReadset(pReadSet);

	free(pReadSet);

	return SUCCESSFUL;
}

/**
 * Release readseq hash blocks.
 */
void releaseHashItemReadset(readSet_t *readSet)
{
	int32_t i;

	// release readseq hash blocks
	if(readSet->readseqHashtable)
	{
		for(i=0; i<readSet->blocksNumReadseqHashItem; i++)
			free(readSet->readseqHashItemBlockArr[i].readseqHashItemArr);
		free(readSet->readseqHashItemBlockArr);
		readSet->readseqHashItemBlockArr = NULL;

		free(readSet->readseqHashtable);
		readSet->readseqHashtable = NULL;
	}
}

/**
 * Compute the readseq hash code.
 *  @return:
 *    The hash code.
 */
uint64_t readseqHashInt(uint64_t *seqInt, int32_t baseNum, int32_t entriesNum)
{
	//return *seqInt;

	uint64_t hashcode;
	int32_t i, j, lastEntryBaseNumTmp;

	lastEntryBaseNumTmp = baseNum % 32;

	hashcode = 5381;
	for(i=0; i<entriesNum-1; i++)
	{
		for(j=0; j<32; j++)
		{
			hashcode += (hashcode << 5) | ((seqInt[i] >> (62-2*j)) & 3);
		}
	}

	for(j=0; j<lastEntryBaseNumTmp; j++)
	{
		hashcode += (hashcode << 5) | ((seqInt[entriesNum-1] >> (2*(lastEntryBaseNumTmp-j-1))) & 3);
	}

	//return (hashcode & 0x7FFFFFFF) % hashTableSizeReadseq;
	return hashcode % hashTableSizeReadseq;
}

/**
 * Generate the kmer interger sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short generateReadseqInt(uint64_t *seqInt, char *seq, int32_t seqLen, int32_t entriesNum)
{
	int32_t i, j, baseInt;
	char *pseq;

	for(i=0; i<entriesNum; i++) seqInt[i] = 0;

	j = 0;
	i = 0;
	pseq = seq;
	while(*pseq)
	{
		switch(*pseq)
		{
			case 'A':
			case 'a': baseInt = 0; break;
			case 'C':
			case 'c': baseInt = 1; break;
			case 'G':
			case 'g': baseInt = 2; break;
			case 'T':
			case 't': baseInt = 3; break;
			default: printf("line=%d, In %s(), base %c, error!\n", __LINE__, __func__, *pseq); return FAILED;
		}

		seqInt[j] = (seqInt[j] << 2) | baseInt;
		i ++;

		if(i%32==0)
			j ++;

		pseq ++;
	}

	return SUCCESSFUL;
}

/**
 * Generate the kmer interger sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
inline readseqHashItem_t *getReadseqHashItemByHash(uint64_t hashvalue, uint64_t *readseqInt, int32_t seqLen, int32_t entriesNum, readSet_t *readSet)
{
	readseqHashItem_t *pReadseqHashItem;
	readseqHashBucket_t *pHashBucket;
	uint64_t *readseq;

	pHashBucket = readSet->readseqHashtable + hashvalue;

	if(pHashBucket->hashItemBlockID>0)
	{
		pReadseqHashItem = readSet->readseqHashItemBlockArr[pHashBucket->hashItemBlockID-1].readseqHashItemArr + pHashBucket->itemRowHashItemBlock;
		while(pReadseqHashItem)
		{
			if(pReadseqHashItem->seqlen==seqLen)
			{
				readseq = readSet->readseqBlockArr[pReadseqHashItem->readseqBlockID-1].readseqArr + pReadseqHashItem->rowReadseqInBlock;
				if(identicalReadseq(readseqInt, readseq, entriesNum)==YES)
					break;
			}

			if(pReadseqHashItem->nextHashItemBlockID>0)
				pReadseqHashItem = readSet->readseqHashItemBlockArr[pReadseqHashItem->nextHashItemBlockID-1].readseqHashItemArr + pReadseqHashItem->nextItemRowHashItemBlock;
			else
				pReadseqHashItem = NULL;
		}
	}
	else
		pReadseqHashItem = NULL;

	return pReadseqHashItem;
}

/**
 * Check whether the two sequence is identical.
 *  @return:
 *  	If identical, return YES; otherwise, return NO.
 */
short identicalReadseq(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2, int32_t entriesNum)
{
	int i;
	for(i=0; i<entriesNum; i++)
	{
		if(kmerSeqInt1[i] != kmerSeqInt2[i])
			return NO;
	}

	return YES;
}

/**
 * Get the read bases from integer.
 */
char *getReadBaseByInt(uint64_t *readseqInt, int32_t seqLen)
{
	int32_t i, j, k, baseInt, entriesNum, lastBaseNum;

	entriesNum = ((seqLen - 1) >> 5) + 1;
	lastBaseNum = seqLen % 32;

	k = 0;
	for(i=0; i<entriesNum-1; i++)
	{
		for(j=0; j<32; j++)
		{
			baseInt = (readseqInt[i] >> (62-2*j)) & 3;
			switch(baseInt)
			{
				case 0: baseSeq[k]='A'; break;
				case 1: baseSeq[k]='C'; break;
				case 2: baseSeq[k]='G'; break;
				case 3: baseSeq[k]='T'; break;
			}
			k ++;
		}
	}

	for(j=0; j<lastBaseNum; j++)
	{
		baseInt = (readseqInt[entriesNum-1] >> (2*(lastBaseNum-j-1))) & 3;
		switch(baseInt)
		{
			case 0: baseSeq[k]='A'; break;
			case 1: baseSeq[k]='C'; break;
			case 2: baseSeq[k]='G'; break;
			case 3: baseSeq[k]='T'; break;
		}
		k ++;
	}

	baseSeq[k] = '\0';

	return baseSeq;
}

/**
 * Get the reverse read bases from integer.
 */
char *getReverseReadBaseByInt(uint64_t *readseqInt, int32_t seqLen)
{
	uint64_t readseqIntRev[((seqLen-1)>>5)+1];
	if(getReverseReadseqInt(readseqIntRev, readseqInt, seqLen)==FAILED)
	{
		printf("line=%d, In %s(), cannot get reversed read bases, error!\n", __LINE__, __func__);
		return NULL;
	}

	return getReadBaseByInt(readseqIntRev, seqLen);
}

/**
 * Get the reversed read integer sequence of a read by its read integer sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReverseReadseqInt(uint64_t *readseqIntRev, uint64_t *readseqInt, int32_t seqLen)
{
	int i, j, baseInt;
	int entryIndexRev, baseNum, entriesNum, baseNumLastReadEntry;

	entriesNum = ((seqLen - 1) >> 5) + 1;
	baseNumLastReadEntry = ((seqLen - 1) % 32) + 1;

	for(i=0; i<entriesNum; i++) readseqIntRev[i] = 0;

	entryIndexRev = baseNum = 0;
	for(j=0; j<baseNumLastReadEntry; j++)
	{
		baseInt = (~(readseqInt[entriesNum-1] >> (j*2))) & 3;
		readseqIntRev[entryIndexRev] = (readseqIntRev[entryIndexRev] << 2) | baseInt;
		baseNum ++;
		if(baseNum==32)
		{
			entryIndexRev ++;
			baseNum = 0;
		}
	}
	for(i=entriesNum-2; i>=0; i--)
	{
		for(j=0; j<32; j++)
		{
			baseInt = (~(readseqInt[i] >> (2*j))) & 3;
			readseqIntRev[entryIndexRev] = (readseqIntRev[entryIndexRev] << 2) | baseInt;
			baseNum ++;
			if(baseNum==32)
			{
				entryIndexRev ++;
				baseNum = 0;
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Generate the kmer interger sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short replaceUnknownBasesInReads(char *seq, int32_t nBaseNum)
{
	char *pseq;
	int32_t nBaseNumTmp;

	nBaseNumTmp = 0;
	pseq = seq;
	while(*pseq)
	{
		if(*pseq=='N' || *pseq=='n' || *pseq=='.')
		{
			*pseq = UNKNOWN_BASE_REPLACE_CHAR;
			nBaseNumTmp ++;
			if(nBaseNumTmp>=nBaseNum)
				break;
		}
		pseq ++;
	}

	return SUCCESSFUL;
}


/**
 * Compute the maximal read length in read set.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short computeMaxReadLenInReadset(readSet_t *readSet)
{
	int64_t i, j, maxReadLenTmp;
	read_t *readArr;
	uint32_t itemNumInReadBlock;

	maxReadLenTmp = 0;
	for(i=0; i<readSet->blocksNumRead; i++)
	{
		readArr = readSet->readBlockArr[i].readArr;
		itemNumInReadBlock = readSet->readBlockArr[i].itemNum;
		for(j=0; j<itemNumInReadBlock; j++)
			if(maxReadLenTmp<readArr[j].seqlen)
				maxReadLenTmp = readArr[j].seqlen;
	}

	readSet->maxReadLen = maxReadLenTmp;

	if(maxReadLenTmp==0)
	{
		printf("line=%d, In %s(), cannot get the maximal read length in read set, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Extract reference position from read herad name.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short extractRefPosFromHeadName(int32_t *strandPosArray, char *headname, int32_t headlen)
{
	int32_t i, j, startPos, underlineNum, arrayRow;
	char tmpStr[2][20];

	startPos = 0;
	underlineNum = 0;
	for(i=headlen-1; i>=0; i--)
	{
		if(headname[i]=='_')
		{
			underlineNum ++;
			if(underlineNum>=2)
			{
				startPos = i + 1;
				break;
			}
		}
	}

	arrayRow = 0;
	j = 0;
	for(i=startPos; i<headlen; i++)
	{
		if(headname[i]=='_')
		{
			arrayRow ++;
			j = 0;
		}else
		{
			tmpStr[arrayRow][j++] = headname[i];
		}
	}

	strandPosArray[0] = atoi(tmpStr[0]);
	strandPosArray[1] = atoi(tmpStr[1]);

	if(strandPosArray[0]==1)
		strandPosArray[0] = STRAND_PLUS;
	else if(strandPosArray[0]==2)
		strandPosArray[0] = STRAND_MINUS;
	else
	{
		printf("line=%d, In %s(), cannot extract read reference strand, strand=%d, error!\n", __LINE__, __func__, strandPosArray[0]);
		return FAILED;
	}

	return SUCCESSFUL;
}