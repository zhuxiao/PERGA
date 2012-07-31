/*
 * graph.c
 *
 *  Created on: Jul 2, 2010
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build De Bruijn graph by single end file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraph(char *graphFileName, char **readsFileNames, int readsFileNum)
{
	printf("\n============= Begin to construct graph, please wait ... =============\n");

	if(readsFileFormatType==FILE_FORMAT_FASTA)
	{
		if(pairedMode==0)
		{
			if(constructGraphBySEFasta(graphFileName, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(pairedMode==1)
		{
			if(constructGraphByPEFastaSeparate(graphFile, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return 1;
			}
		}else if(pairedMode==2)
		{
			if(constructGraphByPEFastaInterleaved(graphFile, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return 1;
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
			if(constructGraphBySEFastq(graphFileName, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(pairedMode==1)
		{
			if(constructGraphByPEFastqSeparate(graphFile, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return 1;
			}
		}else if(pairedMode==2)
		{
			if(constructGraphByPEFastqInterleaved(graphFile, readsFileNames, readsFileNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
				return 1;
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

	printf("============= End constructed graph. =============\n");

	return SUCCESSFUL;
}

/**
 * Build De Bruijn graph by single end file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraphBySEFasta(char *graphFile, char **readsFileNames, int readsFileNum)
{
	struct timeval tpstart,tpend;
	gettimeofday(&tpstart,NULL);

	newKmerNum = 0;
	newRidposNum = 0;

	//统计kmer的数量
	graphtype *graph = ReadDataFilePreBySEFasta(readsFileNames, readsFileNum);
	if(graph)
	{
		//fill kmers
		graph = ReadDataFileBySEFasta(readsFileNames, readsFileNum, graph);
		if(graph==NULL)
		{
			printf("line=%d, In %s(), cannot construct the graph, error!\n", __LINE__, __func__);
			return FAILED;
		}

/*
		// ############################ Debug information ##############################
		if(checkGraph(graph)==FAILED)
		{
			printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
			return FAILED;
		}
		// ############################ Debug information ##############################
*/

		// output the graph to file
		if(outputGraphToFile(graphFile, graph)==FAILED)
		{
			printf("line=%d, In %s(), can not output graph to file. Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), cannot construct graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	releaseGraph(graph);  //释放de Bruijn图的内存

	printf("totalReadNum=%lu, validReadNum=%lu, validRatio=%.4f\n", totalReadNum, validReadNum, (double)validReadNum/totalReadNum);
	printf("newKmerNum=%ld, newRidposNum=%ld, folds=%f, nonempty=%f\n", newKmerNum, newRidposNum, (newKmerNum>0)?(double)newRidposNum/newKmerNum:0, (double)newKmerNum/hashTableSize);

	gettimeofday(&tpend, NULL);
	timeuse_deBruijn = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Graph Used Time: %f Seconds.\n", timeuse_deBruijn);

	return SUCCESSFUL;
}

/**
 * Build De Bruijn graph by single end file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraphBySEFastq(char *graphFile, char **readsFileNames, int readsFileNum)
{
	struct timeval tpstart,tpend;
	gettimeofday(&tpstart,NULL);

	newKmerNum = 0;
	newRidposNum = 0;

	//统计kmer的数量
	graphtype *graph = ReadDataFilePreBySEFastq(readsFileNames, readsFileNum);
	if(graph)
	{
		//填充kmer
		graph = ReadDataFileBySEFastq(readsFileNames, readsFileNum, graph);
		if(graph==NULL)
		{
			printf("line=%d, In %s(), cannot construct the graph, error!\n", __LINE__, __func__);
			return FAILED;
		}

/*
		// ############################ Debug information ##############################
		if(checkGraph(graph)==FAILED)
		{
			printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
			return FAILED;
		}
		// ############################ Debug information ##############################
*/

		// output the graph to file
		if(outputGraphToFile(graphFile, graph)==FAILED)
		{
			printf("line=%d, In %s(), can not output graph to file. Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), cannot construct graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	releaseGraph(graph);  //释放de Bruijn图的内存

	printf("totalReadNum=%lu, validReadNum=%lu, validRatio=%.4f\n", totalReadNum, validReadNum, (double)validReadNum/totalReadNum);
	printf("newKmerNum=%ld, newRidposNum=%ld, folds=%f, nonempty=%f\n", newKmerNum, newRidposNum, (newKmerNum>0)?(double)newRidposNum/newKmerNum:0, (double)newKmerNum/hashTableSize);

	gettimeofday(&tpend, NULL);
	timeuse_deBruijn = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Graph Used Time: %f Seconds.\n", timeuse_deBruijn);

	return SUCCESSFUL;
}

/**
 * Build De Bruijn graph by paired ends in separate files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraphByPEFastaSeparate(char *graphFile, char **readsFileNames, int readsFileNum)
{
	struct timeval tpstart,tpend;
	gettimeofday(&tpstart,NULL);

	newKmerNum = 0;
	newRidposNum = 0;

	//统计kmer的数量
	graphtype *graph = ReadDataFilePreByPEFastaSeparate(readsFileNames, readsFileNum);
	if(graph)
	{
		//填充kmer
		graph = ReadDataFileByPEFastaSeparate(readsFileNames, readsFileNum, graph);
		if(graph==NULL)
		{
			printf("line=%d, In %s(), cannot construct the graph, error!\n", __LINE__, __func__);
			return FAILED;
		}

/*
		// ############################ Debug information ##############################
		if(checkGraph(graph)==FAILED)
		{
			printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
			return FAILED;
		}
		// ############################ Debug information ##############################
*/

		// output the graph to file
		if(outputGraphToFile(graphFile, graph)==FAILED)
		{
			printf("line=%d, In %s(), can not output graph to file. Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), cannot construct graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	releaseGraph(graph);  //释放de Bruijn图的内存

	printf("totalReadNum=%lu, validReadNum=%lu, validRatio=%.4f\n", totalReadNum, validReadNum, (double)validReadNum/totalReadNum);
	printf("newKmerNum=%ld, newRidposNum=%ld, folds=%f, nonempty=%f\n", newKmerNum, newRidposNum, (newKmerNum>0)?(double)newRidposNum/newKmerNum:0, (double)newKmerNum/hashTableSize);

	gettimeofday(&tpend, NULL);
	timeuse_deBruijn = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Graph Used Time: %f Seconds.\n", timeuse_deBruijn);

	return SUCCESSFUL;
}

/**
 * Build De Bruijn graph by paired ends in separate files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraphByPEFastqSeparate(char *graphFile, char **readsFileNames, int readsFileNum)
{
	struct timeval tpstart,tpend;
	gettimeofday(&tpstart,NULL);

	newKmerNum = 0;
	newRidposNum = 0;

	//统计kmer的数量
	graphtype *graph = ReadDataFilePreByPEFastqSeparate(readsFileNames, readsFileNum);
	if(graph)
	{
		//填充kmer
		graph = ReadDataFileByPEFastqSeparate(readsFileNames, readsFileNum, graph);
		if(graph==NULL)
		{
			printf("line=%d, In %s(), cannot construct the graph, error!\n", __LINE__, __func__);
			return FAILED;
		}

/*
		// ############################ Debug information ##############################
		if(checkGraph(graph)==FAILED)
		{
			printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
			return FAILED;
		}
		// ############################ Debug information ##############################
*/

		// output the graph to file
		if(outputGraphToFile(graphFile, graph)==FAILED)
		{
			printf("line=%d, In %s(), can not output graph to file. Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), cannot construct graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	releaseGraph(graph);  //释放de Bruijn图的内存

	printf("totalReadNum=%lu, validReadNum=%lu, validRatio=%.4f\n", totalReadNum, validReadNum, (double)validReadNum/totalReadNum);
	printf("newKmerNum=%ld, newRidposNum=%ld, folds=%f, nonempty=%f\n", newKmerNum, newRidposNum, (newKmerNum>0)?(double)newRidposNum/newKmerNum:0, (double)newKmerNum/hashTableSize);

	gettimeofday(&tpend, NULL);
	timeuse_deBruijn = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Graph Used Time: %f Seconds.\n", timeuse_deBruijn);

	return SUCCESSFUL;
}

/**
 * Build De Bruijn graph by paired ends in interleaved single files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraphByPEFastaInterleaved(char *graphFile, char **readsFileNames, int readsFileNum)
{
	struct timeval tpstart,tpend;
	gettimeofday(&tpstart,NULL);

	newKmerNum = 0;
	newRidposNum = 0;

	//统计kmer的数量
	graphtype *graph = ReadDataFilePreByPEFastaInterleaved(readsFileNames, readsFileNum);
	if(graph)
	{
		//填充kmer
		graph = ReadDataFileByPEFastaInterleaved(readsFileNames, readsFileNum, graph);
		if(graph==NULL)
		{
			printf("line=%d, In %s(), cannot construct the graph, error!\n", __LINE__, __func__);
			return FAILED;
		}

/*
		// ############################ Debug information ##############################
		if(checkGraph(graph)==FAILED)
		{
			printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
			return FAILED;
		}
		// ############################ Debug information ##############################
*/

		// output the graph to file
		if(outputGraphToFile(graphFile, graph)==FAILED)
		{
			printf("line=%d, In %s(), can not output graph to file. Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), cannot construct graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	releaseGraph(graph);  //释放de Bruijn图的内存

	printf("totalReadNum=%lu, validReadNum=%lu, validRatio=%.4f\n", totalReadNum, validReadNum, (double)validReadNum/totalReadNum);
	printf("newKmerNum=%ld, newRidposNum=%ld, folds=%f, nonempty=%f\n", newKmerNum, newRidposNum, (newKmerNum>0)?(double)newRidposNum/newKmerNum:0, (double)newKmerNum/hashTableSize);

	gettimeofday(&tpend, NULL);
	timeuse_deBruijn = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Graph Used Time: %f Seconds.\n", timeuse_deBruijn);

	return SUCCESSFUL;
}

/**
 * Build De Bruijn graph by paired ends in interleaved single files.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short constructGraphByPEFastqInterleaved(char *graphFile, char **readsFileNames, int readsFileNum)
{
	struct timeval tpstart,tpend;
	gettimeofday(&tpstart,NULL);

	newKmerNum = 0;
	newRidposNum = 0;

	//统计kmer的数量
	graphtype *graph = ReadDataFilePreByPEFastqInterleaved(readsFileNames, readsFileNum);
	if(graph)
	{
		//填充kmer
		graph = ReadDataFileByPEFastqInterleaved(readsFileNames, readsFileNum, graph);
		if(graph==NULL)
		{
			printf("line=%d, In %s(), cannot construct the graph, error!\n", __LINE__, __func__);
			return FAILED;
		}

/*
		// ############################ Debug information ##############################
		if(checkGraph(graph)==FAILED)
		{
			printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
			return FAILED;
		}
		// ############################ Debug information ##############################
*/

		// output the graph to file
		if(outputGraphToFile(graphFile, graph)==FAILED)
		{
			printf("line=%d, In %s(), can not output graph to file. Error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), cannot construct graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	releaseGraph(graph);  //释放de Bruijn图的内存

	printf("totalReadNum=%lu, validReadNum=%lu, validRatio=%.4f\n", totalReadNum, validReadNum, (double)validReadNum/totalReadNum);
	printf("newKmerNum=%ld, newRidposNum=%ld, folds=%f, nonempty=%f\n", newKmerNum, newRidposNum, (newKmerNum>0)?(double)newRidposNum/newKmerNum:0, (double)newKmerNum/hashTableSize);

	gettimeofday(&tpend, NULL);
	timeuse_deBruijn = tpend.tv_sec-tpstart.tv_sec+ (double)(tpend.tv_usec-tpstart.tv_usec)/1000000;

	printf("Graph Used Time: %f Seconds.\n", timeuse_deBruijn);

	return SUCCESSFUL;
}

/**
 * Get the read length from fasta file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMinReadLenFromFastaFiles(int *readLenInFile, char **readFilesInput, int readFileNum)
{
	int i, tmpReadLen;

	*readLenInFile = INT_MAX;
	for(i=0; i<readFileNum; i++)
	{
		if(getReadLenFromFasta(&tmpReadLen, readFilesInput[i])==FAILED)
		{
			//printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(tmpReadLen<*readLenInFile)
			*readLenInFile = tmpReadLen;
	}

	return SUCCESSFUL;
}

/**
 * Get the read length from fastq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getMinReadLenFromFastqFiles(int *readLenInFile, char **readFilesInput, int readFileNum)
{
	int i, tmpReadLen;

	*readLenInFile = INT_MAX;
	for(i=0; i<readFileNum; i++)
	{
		if(getReadLenFromFastq(&tmpReadLen, readFilesInput[i])==FAILED)
		{
			//printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
			return FAILED;
		}

		if(tmpReadLen<*readLenInFile)
			*readLenInFile = tmpReadLen;
	}

	return SUCCESSFUL;
}

/**
 * Get the read length from fasta file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int getReadLenFromFasta(int *tmpReadLen, char *fastqFile)
{
	FILE *fpFasta;
	readBuf_t readBuf;
	char seq_data[5000];
	int i, returnCode, tmpLen, validFlag;

	fpFasta = fopen(fastqFile, "r");
	if(fpFasta==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fastqFile);
		return FAILED;
	}

	readBuf.seq = seq_data;

	tmpLen = 0;
	validFlag = YES;
	for(i=0; i<100; i++)
	{
		returnCode = getSingleReadFasta(fpFasta, &readBuf);
		if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}else if(returnCode==FAILED)
		{
			if(tmpLen>0 && tmpLen!=strlen(readBuf.seq))
			{
				validFlag = NO;
				break;
			}

			break;
		}

		if(i==0)
			tmpLen = strlen(readBuf.seq);
		else if(tmpLen!=strlen(readBuf.seq))
		{
			validFlag = NO;
			break;
		}
	}

	if(validFlag==YES)
		*tmpReadLen = tmpLen;
	else
	{
		*tmpReadLen = 0;

		fclose(fpFasta);
		fpFasta = NULL;

		printf("Reads in data sets should be in equal size.\n");
		return FAILED;
	}

	if(*tmpReadLen<kmerSize)
	{
		printf("readLen=%d < kmerSize=%d, error!\n", *tmpReadLen, kmerSize);
		return FAILED;
	}

	fclose(fpFasta);
	fpFasta = NULL;

	return SUCCESSFUL;
}

/**
 * Get the read length from fastq file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int getReadLenFromFastq(int *tmpReadLen, char *fastqFile)
{
	FILE *fpFastq;
	readBuf_t readBuf;
	char seq_data[5000], qual_data[5000];
	int i, returnCode, tmpLen, validFlag;

	fpFastq = fopen(fastqFile, "r");
	if(fpFastq==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fastqFile);
		return FAILED;
	}

	readBuf.seq = seq_data;
	readBuf.qual = qual_data;

	tmpLen = 0;
	validFlag = YES;
	for(i=0; i<100; i++)
	{
		returnCode = getSingleReadFastq(fpFastq, &readBuf);
		if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}else if(returnCode==FAILED)
		{
			if(tmpLen>0 && tmpLen!=strlen(readBuf.seq))
			{
				validFlag = NO;
				break;
			}

			break;
		}

		if(i==0)
			tmpLen = strlen(readBuf.seq);
		else if(tmpLen!=strlen(readBuf.seq))
		{
			validFlag = NO;
			break;
		}
	}

	if(validFlag==YES)
		*tmpReadLen = tmpLen;
	else
	{
		*tmpReadLen = 0;

		fclose(fpFastq);
		fpFastq = NULL;

		printf("Reads in data sets should be in equal size.\n");
		return FAILED;
	}

	if(*tmpReadLen<kmerSize)
	{
		printf("readLen=%d < kmerSize=%d, error!\n", *tmpReadLen, kmerSize);
		return FAILED;
	}

	fclose(fpFastq);
	fpFastq = NULL;

	return SUCCESSFUL;
}

/**
 * 统计de Bruijn图的信息.
 *
 *   @ return :
 *      成功, 返回graph指针; 失败, 返回NULL.
 */
graphtype *ReadDataFilePreBySEFasta(char **readsFileNames, int readsFileNum)
{
	FILE* srcfp;
	int i, line_index, tmpFileID;
	char ch, seq_data[5000]; // read sequence

	totalReadNum = 0;
	validReadNum = 0;

	graphtype *graph = initgraph(); //初始化de Bruijn图
	if(!graph)
	{
		printf("line=%d, In %s(), cannot initialize the graph, error!\n", __LINE__, __func__);
		return NULL;
	}

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if(!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		line_index = 0;
		ch = fgetc(srcfp);
		while(!feof(srcfp))
		{
			while(ch!='\n') ch = fgetc(srcfp);

			i = 0;
			ch = fgetc(srcfp);
			while(ch!='>' && ch!=-1)
			{
				if(ch!='\n' && i<readLen)
					seq_data[i++] = ch;
				ch = fgetc(srcfp);
			}
			seq_data[i] = '\0';

			totalReadNum ++;

			//printf("len = %d, %s\n", len, seq_data);
			if(containUnknownBase(seq_data)==NO && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
			{ //不包含‘N’
				//readnum++;
				validReadNum ++;
				newRidposNum += readLen - kmerSize + 1;

				//count the kmers
				if(addReadPre(seq_data, graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
					releaseGraph(graph);
					return NULL;
				}

			}

			line_index = 0;

		}

		fclose(srcfp);
		srcfp = NULL;
	}

	return graph;
}

/**
 * 统计de Bruijn图的信息.
 *
 *   @ return :
 *      成功, 返回graph指针; 失败, 返回NULL.
 */
graphtype *ReadDataFilePreBySEFastq(char **readsFileNames, int readsFileNum)
{
	FILE* srcfp;
	int i, line_index, tmpFileID;
	char ch, seq_data[5000], qual_data[5000]; // read sequence, read quality data which is encoded in ASCII

	totalReadNum = 0;
	validReadNum = 0;

	graphtype *graph = initgraph(); //初始化de Bruijn图
	if(!graph)
	{
		printf("line=%d, In %s(), cannot initialize the graph, error!\n", __LINE__, __func__);
		return NULL;
	}

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if(!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
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
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n')
				{
					if(i<readLen)
						seq_data[i++] = ch;
					ch = fgetc(srcfp);
				}
				seq_data[i] = '\0';
			}else if(line_index==2)  //the sequence name line
			{
				ch = fgetc(srcfp);
				while(ch!='\n')
				{
					ch = fgetc(srcfp);
				}
			}else
			{
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n'  && ch!=-1)
				{
					if(i<readLen)
						qual_data[i++] = ch;
					ch = fgetc(srcfp);
				}
				qual_data[i] = '\0';
			}
			line_index++;

			if(line_index==4)  //the sequence is read finished, construct the read
			{
				totalReadNum ++;

				//printf("len = %d, %s\n", len, seq_data);
				//if(contianUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
				if(containUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && qualSatisfied(qual_data)==YES && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
				{ //不包含‘N’, 并且相应碱基的平均质量大于15
					//readnum++;
					validReadNum ++;
					newRidposNum += readLen - kmerSize + 1;

					//count the kmers
					if(addReadPre(seq_data, graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
						releaseGraph(graph);
						return NULL;
					}

				}

				line_index = 0;
			}
		}

		fclose(srcfp);
		srcfp = NULL;
	}

	return graph;
}

graphtype *ReadDataFilePreByPEFastaSeparate(char **readsFileNames, int readsFileNum)
{
	readBuf_t *readBuf[2];
	FILE *fp1, *fp2;
	uint64_t i, j, tmpReadsNum[2], tmpFileID;
	char *seq_data;

	for(i=0; i<2; i++)
	{
		readBuf[i] = (readBuf_t*)calloc(MAX_READ_BUF_SIZE, sizeof(readBuf_t));
		if(readBuf[i]==NULL)
		{
			printf("In %s(), can not allocate memory, Error.\n", __func__);
			return NULL;
		}
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			readBuf[i][j].seq = (char *)malloc(sizeof(char)*(readLenInFile+1));
			if(readBuf[i][j].seq==NULL)
			{
				printf("In %s(), can not allocate memory, Error.\n", __func__);
				return NULL;
			}
		}
	}

	graphtype *graph = initgraph(); //初始化de Bruijn图
	if(!graph)
	{
		printf("In %s(), init graph error!\n", __func__);
		return NULL;
	}

	totalReadNum = 0;
	validReadNum = 0;

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID+=2)
	{
		fp1 = fopen(readsFileNames[tmpFileID], "r");
		if(fp1==NULL)
		{
			printf("In %s(), can not open file [ %s ], error.\n", __func__, readsFileNames[tmpFileID]);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}
		fp2 = fopen(readsFileNames[tmpFileID+1], "r");
		if(fp2==NULL)
		{
			printf("In %s(), can not open file [ %s ], error.\n", __func__, readsFileNames[tmpFileID+1]);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID+1]);
		}

		while(1)
		{
			// check the end of files
			if(feof(fp1) && feof(fp2))
			{
				break;
			}else if(feof(fp1) || feof(fp2))
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}

			// fille the reads to reads buffers
			if(fillReadsToBufFasta(fp1, readBuf[0], tmpReadsNum)==FAILED)
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}
			if(fillReadsToBufFasta(fp2, readBuf[1], tmpReadsNum+1)==FAILED)
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}
			if(tmpReadsNum[0]!=tmpReadsNum[1])
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}

			//...
			for(i=0; i<tmpReadsNum[0]; i++)
			{
				for(j=0; j<2; j++)
				{
					seq_data = readBuf[j][i].seq;

					if(readLen<readLenInFile)
						seq_data[readLen] = '\0';

					//printf("len = %d, %s\n", len, seq_data);
					if(containUnknownBase(seq_data)==NO  && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
					{ //不包含‘N’, 并且相应碱基的平均质量大于15
						//readnum++;
						validReadNum ++;
						newRidposNum += readLen - kmerSize + 1;

						//count the kmers
						if(addReadPre(seq_data, graph)==FAILED)
						{
							printf("line=%d, In %s(), cannot count kmer, error!\n", __LINE__, __func__);
							releaseGraph(graph);
							return NULL;
						}
					}
					totalReadNum ++;
				}
			}
		}

		fclose(fp1);
		fp1 = NULL;
		fclose(fp2);
		fp2 = NULL;
	}

	for(i=0; i<2; i++)
	{
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			free(readBuf[i][j].seq);
		}
		free(readBuf[i]);
	}

	return graph;
}

graphtype *ReadDataFilePreByPEFastqSeparate(char **readsFileNames, int readsFileNum)
{
	readBuf_t *readBuf[2];
	FILE *fp1, *fp2;
	uint64_t i, j, tmpReadsNum[2], tmpFileID;
	char *seq_data, *qual_data;

	for(i=0; i<2; i++)
	{
		readBuf[i] = (readBuf_t*)calloc(MAX_READ_BUF_SIZE, sizeof(readBuf_t));
		if(readBuf[i]==NULL)
		{
			printf("In %s(), can not allocate memory, Error.\n", __func__);
			return NULL;
		}
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			readBuf[i][j].seq = (char *)malloc(sizeof(char)*(readLenInFile+1));
			if(readBuf[i][j].seq==NULL)
			{
				printf("In %s(), can not allocate memory, Error.\n", __func__);
				return NULL;
			}
			readBuf[i][j].qual = (char *)malloc(sizeof(char)*(readLenInFile+1));
			if(readBuf[i][j].qual==NULL)
			{
				printf("In %s(), can not allocate memory, Error.\n", __func__);
				return NULL;
			}
		}
	}

	graphtype *graph = initgraph(); //初始化de Bruijn图
	if(!graph)
	{
		printf("In %s(), init graph error!\n", __func__);
		return NULL;
	}

	totalReadNum = 0;
	validReadNum = 0;

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID+=2)
	{
		fp1 = fopen(readsFileNames[tmpFileID], "r");
		if(fp1==NULL)
		{
			printf("In %s(), can not open file [ %s ], error.\n", __func__, readsFileNames[tmpFileID]);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}
		fp2 = fopen(readsFileNames[tmpFileID+1], "r");
		if(fp2==NULL)
		{
			printf("In %s(), can not open file [ %s ], error.\n", __func__, readsFileNames[tmpFileID+1]);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID+1]);
		}

		while(1)
		{
			// check the end of files
			if(feof(fp1) && feof(fp2))
			{
				break;
			}else if(feof(fp1) || feof(fp2))
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}

			// fille the reads to reads buffers
			if(fillReadsToBuf(fp1, readBuf[0], tmpReadsNum)==FAILED)
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}
			if(fillReadsToBuf(fp2, readBuf[1], tmpReadsNum+1)==FAILED)
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}
			if(tmpReadsNum[0]!=tmpReadsNum[1])
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}

			//...
			for(i=0; i<tmpReadsNum[0]; i++)
			{
				for(j=0; j<2; j++)
				{
					seq_data = readBuf[j][i].seq;
					qual_data = readBuf[j][i].qual;

					if(readLen<readLenInFile)
					{
						seq_data[readLen] = '\0';
						qual_data[readLen] = '\0';
					}

					//printf("len = %d, %s\n", len, seq_data);
					//if(contianUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
					if(containUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && qualSatisfied(qual_data)==YES && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
					{ //不包含‘N’, 并且相应碱基的平均质量大于15
						//readnum++;
						validReadNum ++;
						newRidposNum += readLen - kmerSize + 1;

						//count the kmers
						if(addReadPre(seq_data, graph)==FAILED)
						{
							printf("line=%d, In %s(), cannot count kmer, error!\n", __LINE__, __func__);
							releaseGraph(graph);
							return NULL;
						}
					}
					totalReadNum ++;
				}
			}
		}

		fclose(fp1);
		fp1 = NULL;
		fclose(fp2);
		fp2 = NULL;
	}

	for(i=0; i<2; i++)
	{
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			free(readBuf[i][j].seq);
			free(readBuf[i][j].qual);
		}
		free(readBuf[i]);
	}

	return graph;
}

graphtype *ReadDataFilePreByPEFastaInterleaved(char **readsFileNames, int readsFileNum)
{
	FILE* srcfp;
	int i, pairedID, tmpFileID;
	char ch, seq_data[2][5000]; // read sequence, read quality data which is encoded in ASCII

	totalReadNum = 0;
	validReadNum = 0;

	graphtype *graph = initgraph(); //初始化de Bruijn图
	if(!graph)
	{
		printf("line=%d, In %s(), cannot initialize the graph, error!\n", __LINE__, __func__);
		return NULL;
	}

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if(!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		pairedID = 0;
		ch = fgetc(srcfp);
		while(!feof(srcfp))
		{
			while(ch!='\n') ch = fgetc(srcfp);

			i = 0;
			ch = fgetc(srcfp);
			while(ch!='>' && ch!=-1)
			{
				if(ch!='\n' && i<readLen)
					seq_data[pairedID][i++] = ch;
				ch = fgetc(srcfp);
			}
			seq_data[pairedID][i] = '\0';

			totalReadNum ++;

			//printf("len = %d, %s\n", len, seq_data);
			if(containUnknownBase(seq_data[pairedID])==NO && getRatioBaseA(seq_data[pairedID])<ARTIFACTS_BASE_A_THRESHOLD)
			{ //不包含‘N’
				//readnum++;
				validReadNum ++;
				newRidposNum += readLen - kmerSize + 1;

				//count the kmers
				if(addReadPre(seq_data[pairedID], graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
					releaseGraph(graph);
					return NULL;
				}

			}

			if(pairedID==1)
				pairedID = 0;
			else
				pairedID ++;

		}

		fclose(srcfp);
		srcfp = NULL;
	}

	return graph;
}

graphtype *ReadDataFilePreByPEFastqInterleaved(char **readsFileNames, int readsFileNum)
{
	FILE* srcfp;
	int i, line_index, pairedID, tmpFileID;
	char ch, seq_data[2][5000], qual_data[2][5000]; // read sequence, read quality data which is encoded in ASCII

	totalReadNum = 0;
	validReadNum = 0;

	graphtype *graph = initgraph(); //初始化de Bruijn图
	if(!graph)
	{
		printf("line=%d, In %s(), cannot initialize the graph, error!\n", __LINE__, __func__);
		return NULL;
	}

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if(!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		pairedID = 0;
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
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n')
				{
					if(i<readLen)
						seq_data[pairedID][i++] = ch;
					ch = fgetc(srcfp);
				}
				seq_data[pairedID][i] = '\0';
			}else if(line_index==2)  //the sequence name line
			{
				ch = fgetc(srcfp);
				while(ch!='\n')
				{
					ch = fgetc(srcfp);
				}
			}else
			{
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n'  && ch!=-1)
				{
					if(i<readLen)
						qual_data[pairedID][i++] = ch;
					ch = fgetc(srcfp);
				}
				qual_data[pairedID][i] = '\0';
			}
			line_index++;

			if(line_index==4)  //the sequence is read finished, construct the read
			{
				totalReadNum ++;

				//printf("len = %d, %s\n", len, seq_data);
				//if(contianUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
				if(containUnknownBase(seq_data[pairedID])==NO && calcAverQual5End(qual_data[pairedID])>=AVERAGE_QUAL_THRESHOLD_5End && qualSatisfied(qual_data[pairedID])==YES && calcAverQual3End(qual_data[pairedID])>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data[pairedID])<ARTIFACTS_BASE_A_THRESHOLD)
				{ //不包含‘N’, 并且相应碱基的平均质量大于15
					//readnum++;
					validReadNum ++;
					newRidposNum += readLen - kmerSize + 1;

					//count the kmers
					if(addReadPre(seq_data[pairedID], graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot count the kmer, error!\n", __LINE__, __func__);
						releaseGraph(graph);
						return NULL;
					}

				}

				line_index = 0;

				if(pairedID==1)
					pairedID = 0;
				else
					pairedID ++;

			}
		}

		fclose(srcfp);
		srcfp = NULL;
	}

	return graph;
}

/**
 * Fill reads to reads buffer.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadsToBufFasta(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum)
{
	uint64_t i;
	int returnCode;

	*readsNum = 0;
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		returnCode = getSingleReadFasta(fpReads, pBuf+i);
		if(returnCode==FAILED)
			break;
		else if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*readsNum = i;

	return SUCCESSFUL;
}

/**
 * Fill reads to reads buffer.
 *  @return:
 *  	If succeed, return SUCCESSFUL; otherwise, return FAILED.
 */
short fillReadsToBuf(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum)
{
	uint64_t i;
	int returnCode;

	*readsNum = 0;
	for(i=0; i<MAX_READ_BUF_SIZE; i++)
	{
		returnCode = getSingleReadFastq(fpReads, pBuf+i);
		if(returnCode==FAILED)
			break;
		else if(returnCode==ERROR)
		{
			printf("line=%d, In %s(), cannot get single read, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	*readsNum = i;

	return SUCCESSFUL;
}

/**
 * get the single read from the given file.
 *  @return:
 *  	If succeed, return SUCCESSFUL;
 *  	else if the file end is reached, return FAILED;
 *  	otherwise, return ERROR.
 */
short getSingleReadFasta(FILE *fpPE, readBuf_t *pReadBuf)
{
	int i;
	//char qual_data[5000];	// read quality data which is encoded in ASCII
	char *readSeq;
	readSeq = pReadBuf->seq;

	char ch = fgetc(fpPE);
	if(feof(fpPE))// the file end is reached.
	{
		return FAILED;
	}

	while(ch!='\n') ch = fgetc(fpPE);

	while(!feof(fpPE))
	{
		i = 0;
		ch = fgetc(fpPE);
		while(ch!='>' && ch!=-1)
		{
			if(ch!='\n')
				readSeq[i++] = ch;
			ch = fgetc(fpPE);
		}
		readSeq[i] = '\0';

		if(ch=='>')  //a read is read finished
		{
			break;
		}
	}

	return SUCCESSFUL;
}

/**
 * get the single read from the given file.
 *  @return:
 *  	If succeed, return SUCCESSFUL;
 *  	else if the file end is reached, return FAILED;
 *  	otherwise, return ERROR.
 */
short getSingleReadFastq(FILE *fpPE, readBuf_t *pReadBuf)
{
	unsigned short i, line_index = 0;
	//char qual_data[5000];	// read quality data which is encoded in ASCII
	char *readSeq, *qual_data;
	readSeq = pReadBuf->seq;
	qual_data = pReadBuf->qual;

	char ch = fgetc(fpPE);
	if(feof(fpPE))// the file end is reached.
	{
		return FAILED;
	}

	while(!feof(fpPE))
	{
		if(line_index==0)  //the sequence name line
		{
			i = 0;
			ch = fgetc(fpPE);
			while(ch!='\n' && ch!=-1)
			{
				ch = fgetc(fpPE);
			}
		}else if(line_index==1)  //the sequence line
		{
			i = 0;
			ch = fgetc(fpPE);
			while(ch!='\n')
			{
				readSeq[i++] = ch;
				ch = fgetc(fpPE);
			}
			readSeq[i] = '\0';
		}else if(line_index==2)  //the sequence name line
		{
			ch = fgetc(fpPE);
			while(ch!='\n')
			{
				ch = fgetc(fpPE);
			}
		}else
		{
			i = 0;
			ch = fgetc(fpPE);
			while(ch!='\n'  && ch!=-1)
			{
				qual_data[i++] = ch;
				ch = fgetc(fpPE);
			}
			qual_data[i] = '\0';
		}
		line_index++;

		if(line_index==4)  //a read is read finished
		{
			break;
		}
	}


	// check the return code
	if(line_index==4)
	{
		return SUCCESSFUL;
	}else
	{
		return ERROR;
	}
}


/**
 * 检测read中是否有未知碱基.
 *  @ return:
 *  	如果有, 返回YES; 否则, 返回NO.
 */
short containUnknownBase(char *seq)
{
	char *p_str = seq;
	while(*p_str)
	{
		if(*p_str=='N' || *p_str=='.' || *p_str=='n')
			return YES;
		p_str++;
	}
	return NO;
}

short qualSatisfied(char *qual_data)
{
	char *p_qual = qual_data;
	int i;
	//for(i=0; i<READ_LEN-ERROR_REGION_LEN_3End; i++)
	for(i=0; i<readLen; i++)
	{
		if(*p_qual-33 < SINGLE_QUAL_THRESHOLD)
			return NO;
		p_qual ++;
	}
/*
	for(i=0; i<QUAL_BASE_NUM_5End; i++)
	{
		if(*p_qual-33 < SINGLE_QUAL_THRESHOLD)
			return NO;
		p_qual ++;
	}

	p_qual = qual_data + READ_LEN - QUAL_BASE_NUM_3End;
	for(i=0; i<QUAL_BASE_NUM_3End-ERROR_REGION_LEN_3End; i++)
	{
		if(*p_qual-33 < SINGLE_QUAL_THRESHOLD)
			return NO;
		p_qual ++;
	}
*/
	return YES;
}

/**
 * 计算read的3'端碱基的平均质量.
 */
float calcAverQual3End(char *qual_data)
{
	float sumQual = 0;
	char *p_qual = qual_data + readLen - qualityBaseNumEnd3;
	while(*p_qual)
	{
		sumQual += (*p_qual) - 33;
		p_qual ++;
	}
	//return sumQual/READ_LEN;
	return sumQual/qualityBaseNumEnd3;
}

/**
 * 计算read的5'端碱基的平均质量.
 */
float calcAverQual5End(char *qual_data)
{
	float sumQual = 0;
	char *p_qual = qual_data;
	int i;
	for(i=0; i<qualityBaseNumEnd5; i++)
	{
		sumQual += (*p_qual) - 33;
		p_qual ++;
	}

	return sumQual/qualityBaseNumEnd5;
}


/**
 * 统计read中A的百分比.
 * 	@ return:
 * 		返回该百分比.
 */
float getRatioBaseA(char *seq)
{
	float ratioA = 0, baseNum = 0;
	char *p_seq = seq;
	while(*p_seq)
	{
		if(*p_seq=='A')
		{
			ratioA ++;
		}
		baseNum++;
		p_seq++;
	}
	return ratioA/baseNum;
}

/**
 * Count the kmers of a read.
 *  @return:
 *  	If succeeds, return SUCESSFUL; otherwise, return FAILED.
 */
short addReadPre(char *seq, graphtype *graph)
{
	int i, j, baseInt;
	uint64_t hashcode;

	// generate the kmer integer sequence
	if(generateKmerSeqInt(kmerSeqInt, seq)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the hashcode and count the kmer occurrence
	hashcode = kmerhashInt(kmerSeqInt);
	if(countKmer(hashcode, kmerSeqInt, graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=kmerSize; i<readLen; i++)
	{
		// get the baseInt
		switch(seq[i])
		{
			case 'A':
			case 'a': baseInt = 0; break;
			case 'C':
			case 'c': baseInt = 1; break;
			case 'G':
			case 'g': baseInt = 2; break;
			case 'T':
			case 't': baseInt = 3; break;
			default: printf("line=%d, In %s(), base %c, error!\n", __LINE__, __func__, seq[i]); return FAILED;
		}

		// generate the kmer integer sequence
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				kmerSeqInt[j] = (kmerSeqInt[j] << 2) | (kmerSeqInt[j+1] >> 62);
			}

			kmerSeqInt[entriesPerKmer-2] = (kmerSeqInt[entriesPerKmer-2] << 2) | (kmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		kmerSeqInt[entriesPerKmer-1] = ((kmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

		// get the hashcode and count the kmer occurrence
		hashcode = kmerhashInt(kmerSeqInt);
		if(countKmer(hashcode, kmerSeqInt, graph)==FAILED)
		{
			printf("line=%d, In %s(), cannot count kmer occurrence, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Generate the kmer interger sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int generateKmerSeqInt(uint64_t *seqInt, char *seq)
{
	int i, j, baseInt;

	for(i=0; i<entriesPerKmer; i++) seqInt[i] = 0;

	j = 0;
	i = 0;
	while(i<kmerSize)
	{
		switch(seq[i])
		{
			case 'A':
			case 'a': baseInt = 0; break;
			case 'C':
			case 'c': baseInt = 1; break;
			case 'G':
			case 'g': baseInt = 2; break;
			case 'T':
			case 't': baseInt = 3; break;
			default: printf("line=%d, In %s(), base %c, error!\n", __LINE__, __func__, seq[i]); return FAILED;
		}

		seqInt[j] = (seqInt[j] << 2) | baseInt;
		i ++;

		if(i%32==0)
			j ++;
	}

	return SUCCESSFUL;
}

/**
 * fill kmers.
 *   @ return:
 *     成功, 返回graph指针; 失败, 返回NULL.
 */
graphtype *ReadDataFileBySEFasta(char **readsFileNames, int readsFileNum, graphtype *graph)
{
	int i, line_index, tmpFileID;
	int64_t readnum;  //read number
	char ch, seq_data[5000]; // read sequence, read quality data which is encoded in ASCII
	FILE* srcfp;

	readnum = 0;  //read number

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if (!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			releaseGraph(graph);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		line_index = 0;
		ch = fgetc(srcfp);
		while(!feof(srcfp))
		{
			while(ch!='\n') ch = fgetc(srcfp);

			i = 0;
			ch = fgetc(srcfp);
			while(ch!='>' && ch!=-1)
			{
				if(ch!='>' && i<readLen)
					seq_data[i++] = ch;
				ch = fgetc(srcfp);
			}
			seq_data[i] = '\0';


			readnum ++;

			//printf("len = %d, %s\n", len, seq_data);
			if(containUnknownBase(seq_data)==NO && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
			{ //不包含‘N’

				//add the read
				if(addRead(seq_data, readnum, graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot add a read to graph, error!\n", __LINE__, __func__);
					releaseGraph(graph);
					return NULL;
				}

			}

		}

		fclose(srcfp);
		srcfp = NULL;
	}

	return graph;
}

/**
 * fill kmers.
 *   @ return:
 *     成功, 返回graph指针; 失败, 返回NULL.
 */
graphtype *ReadDataFileBySEFastq(char **readsFileNames, int readsFileNum, graphtype *graph)
{
	int i, line_index, tmpFileID;
	int64_t readnum;  //read number
	char ch, seq_data[5000], qual_data[5000]; // read sequence, read quality data which is encoded in ASCII
	FILE* srcfp;

	readnum = 0;  //read number

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if (!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			releaseGraph(graph);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
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
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n')
				{
					if(i<readLen)
						seq_data[i++] = ch;
					ch = fgetc(srcfp);
				}
				seq_data[i] = '\0';
			}else if(line_index==2)  //the sequence name line
			{
				ch = fgetc(srcfp);
				while(ch!='\n')
				{
					ch = fgetc(srcfp);
				}
			}else
			{
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n'  && ch!=-1)
				{
					if(i<readLen)
						qual_data[i++] = ch;
					ch = fgetc(srcfp);
				}
				qual_data[i] = '\0';
			}
			line_index++;

			if(line_index==4)  //the sequence is read finished
			{
				readnum ++;

				//printf("len = %d, %s\n", len, seq_data);
				//if(contianUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
				if(containUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && qualSatisfied(qual_data)==YES && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
				{ //不包含‘N’, 并且相应碱基的平均质量大于15

					//add the read
					if(addRead(seq_data, readnum, graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot add a read to graph, error!\n", __LINE__, __func__);
						releaseGraph(graph);
						return NULL;
					}

				}

				line_index = 0;
			}
		}

		fclose(srcfp);
		srcfp = NULL;
	}

	return graph;
}

graphtype *ReadDataFileByPEFastaSeparate(char **readsFileNames, int readsFileNum, graphtype *graph)
{
	readBuf_t *readBuf[2];
	FILE *fp1, *fp2;
	uint64_t i, j, tmpReadsNum[2];
	char *seq_data;
	uint64_t readID, tmpFileID;

	for(i=0; i<2; i++)
	{
		readBuf[i] = (readBuf_t*)calloc(MAX_READ_BUF_SIZE, sizeof(readBuf_t));
		if(readBuf[i]==NULL)
		{
			printf("In %s(), can not allocate memory, Error.\n", __func__);
			return NULL;
		}
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			readBuf[i][j].seq = (char *)malloc(sizeof(char)*(readLenInFile+1));
			if(readBuf[i][j].seq==NULL)
			{
				printf("In %s(), can not allocate memory, Error.\n", __func__);
				return NULL;
			}
		}
	}

	readID = 1;
	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID+=2)
	{
		fp1 = fopen(readsFileNames[tmpFileID], "r");
		if(fp1==NULL)
		{
			printf("In %s(), can not open file [ %s ], Rrror.\n", __func__, readsFileNames[tmpFileID]);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}
		fp2 = fopen(readsFileNames[tmpFileID+1], "r");
		if(fp2==NULL)
		{
			printf("In %s(), can not open file [ %s ], Rrror.\n", __func__, readsFileNames[tmpFileID+1]);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID+1]);
		}

		while(1)
		{
			// check the end of files
			if(feof(fp1) && feof(fp2))
			{
				break;
			}else if(feof(fp1) || feof(fp2))
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}

			// fille the reads to reads buffers
			if(fillReadsToBufFasta(fp1, readBuf[0], tmpReadsNum)==FAILED)
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}
			if(fillReadsToBufFasta(fp2, readBuf[1], tmpReadsNum+1)==FAILED)
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}
			if(tmpReadsNum[0]!=tmpReadsNum[1])
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}

			//...
			for(i=0; i<tmpReadsNum[0]; i++)
			{
				for(j=0; j<2; j++)
				{
					seq_data = readBuf[j][i].seq;

					if(readLen<readLenInFile)
						seq_data[readLen] = '\0';

					//printf("len = %d, %s\n", len, seq_data);
					if(containUnknownBase(seq_data)==NO && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
					{ //不包含‘N’

						//add the read
						if(addRead(seq_data, readID, graph)==FAILED)
						{
							printf("line=%d, In %s(), cannot add read to graph, error!\n", __LINE__, __func__);
							releaseGraph(graph);
							return NULL;
						}

					}
					readID ++;
				}
			}
		}

		fclose(fp1);
		fp1 = NULL;
		fclose(fp2);
		fp2 = NULL;
	}

	for(i=0; i<2; i++)
	{
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
			free(readBuf[i][j].seq);
		free(readBuf[i]);
	}

	return graph;
}

graphtype *ReadDataFileByPEFastqSeparate(char **readsFileNames, int readsFileNum, graphtype *graph)
{
	readBuf_t *readBuf[2];
	FILE *fp1, *fp2;
	uint64_t i, j, tmpReadsNum[2];
	char *seq_data, *qual_data;
	uint64_t readID, tmpFileID;

	for(i=0; i<2; i++)
	{
		readBuf[i] = (readBuf_t*)calloc(MAX_READ_BUF_SIZE, sizeof(readBuf_t));
		if(readBuf[i]==NULL)
		{
			printf("In %s(), can not allocate memory, Error.\n", __func__);
			return NULL;
		}
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			readBuf[i][j].seq = (char *)malloc(sizeof(char)*(readLenInFile+1));
			if(readBuf[i][j].seq==NULL)
			{
				printf("In %s(), can not allocate memory, Error.\n", __func__);
				return NULL;
			}
			readBuf[i][j].qual = (char *)malloc(sizeof(char)*(readLenInFile+1));
			if(readBuf[i][j].qual==NULL)
			{
				printf("In %s(), can not allocate memory, Error.\n", __func__);
				return NULL;
			}
		}
	}

	readID = 1;
	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID+=2)
	{
		fp1 = fopen(readsFileNames[tmpFileID], "r");
		if(fp1==NULL)
		{
			printf("In %s(), can not open file [ %s ], Rrror.\n", __func__, readsFileNames[tmpFileID]);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}
		fp2 = fopen(readsFileNames[tmpFileID+1], "r");
		if(fp2==NULL)
		{
			printf("In %s(), can not open file [ %s ], Rrror.\n", __func__, readsFileNames[tmpFileID+1]);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID+1]);
		}

		while(1)
		{
			// check the end of files
			if(feof(fp1) && feof(fp2))
			{
				break;
			}else if(feof(fp1) || feof(fp2))
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}

			// fille the reads to reads buffers
			if(fillReadsToBuf(fp1, readBuf[0], tmpReadsNum)==FAILED)
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}
			if(fillReadsToBuf(fp2, readBuf[1], tmpReadsNum+1)==FAILED)
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}
			if(tmpReadsNum[0]!=tmpReadsNum[1])
			{
				printf("In %s(), cannot fill the read buffer, error!\n", __func__);
				return NULL;
			}

			//...
			for(i=0; i<tmpReadsNum[0]; i++)
			{
				for(j=0; j<2; j++)
				{
					seq_data = readBuf[j][i].seq;
					qual_data = readBuf[j][i].qual;

					if(readLen<readLenInFile)
					{
						seq_data[readLen] = '\0';
						qual_data[readLen] = '\0';
					}

					//printf("len = %d, %s\n", len, seq_data);
					//if(contianUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
					if(containUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && qualSatisfied(qual_data)==YES && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
					{ //不包含‘N’, 并且相应碱基的平均质量大于15

						//add the read
						if(addRead(seq_data, readID, graph)==FAILED)
						{
							printf("line=%d, In %s(), cannot add read to graph, error!\n", __LINE__, __func__);
							releaseGraph(graph);
							return NULL;
						}

					}
					readID ++;
				}
			}
		}

		fclose(fp1);
		fp1 = NULL;
		fclose(fp2);
		fp2 = NULL;
	}

	for(i=0; i<2; i++)
	{
		for(j=0; j<MAX_READ_BUF_SIZE; j++)
		{
			free(readBuf[i][j].seq);
			free(readBuf[i][j].qual);
		}
		free(readBuf[i]);
	}

	return graph;
}

/**
 * fill kmers.
 *   @ return:
 *     成功, 返回graph指针; 失败, 返回NULL.
 */
graphtype *ReadDataFileByPEFastaInterleaved(char **readsFileNames, int readsFileNum, graphtype *graph)
{
	int i, pairedID, tmpFileID;
	uint64_t readnum;  //read number
	char ch, seq_data[2][5000]; // read sequence, read quality data which is encoded in ASCII
	FILE* srcfp;

	readnum = 0;  //read number

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if (!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			releaseGraph(graph);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		pairedID = 0;
		ch = fgetc(srcfp);
		while(!feof(srcfp))
		{
			while(ch!='\n') ch = fgetc(srcfp);

			i = 0;
			ch = fgetc(srcfp);
			while(ch!='>' && ch!=-1)
			{
				if(ch!='\n' && i<readLen)
					seq_data[pairedID][i++] = ch;
				ch = fgetc(srcfp);
			}
			seq_data[pairedID][i] = '\0';

			readnum ++;

			//printf("len = %d, %s\n", len, seq_data);
			if(containUnknownBase(seq_data[pairedID])==NO && getRatioBaseA(seq_data[pairedID])<ARTIFACTS_BASE_A_THRESHOLD)
			{ //不包含‘N’, 并且相应碱基的平均质量大于15

				//add the read
				if(addRead(seq_data[pairedID], readnum, graph)==FAILED)
				{
					printf("line=%d, In %s(), cannot add a read to graph, error!\n", __LINE__, __func__);
					releaseGraph(graph);
					return NULL;
				}

			}

			if(pairedID==1)
				pairedID = 0;
			else
				pairedID ++;

		}

		fclose(srcfp);
		srcfp = NULL;
	}

	return graph;
}

/**
 * fill kmers.
 *   @ return:
 *     成功, 返回graph指针; 失败, 返回NULL.
 */
graphtype *ReadDataFileByPEFastqInterleaved(char **readsFileNames, int readsFileNum, graphtype *graph)
{
	int i, line_index, pairedID, tmpFileID;
	uint64_t readnum;  //read number
	char ch, seq_data[2][5000], qual_data[2][5000]; // read sequence, read quality data which is encoded in ASCII
	FILE* srcfp;

	readnum = 0;  //read number

	for(tmpFileID=0; tmpFileID<readsFileNum; tmpFileID++)
	{
		srcfp = fopen(readsFileNames[tmpFileID], "r");
		if (!srcfp)
		{
			printf("Data File [%s] cannot open! Error Code:%d\n", readsFileNames[tmpFileID], errno);
			releaseGraph(graph);
			return NULL;
		}else
		{
			printf("Data File [%s] opened ok!\n", readsFileNames[tmpFileID]);
		}

		pairedID = 0;
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
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n')
				{
					if(i<readLen)
						seq_data[pairedID][i++] = ch;
					ch = fgetc(srcfp);
				}
				seq_data[pairedID][i] = '\0';
			}else if(line_index==2)  //the sequence name line
			{
				ch = fgetc(srcfp);
				while(ch!='\n')
				{
					ch = fgetc(srcfp);
				}
			}else
			{
				i = 0;
				ch = fgetc(srcfp);
				while(ch!='\n'  && ch!=-1)
				{
					if(i<readLen)
						qual_data[pairedID][i++] = ch;
					ch = fgetc(srcfp);
				}
				qual_data[pairedID][i] = '\0';
			}
			line_index++;

			if(line_index==4)  //the sequence is read finished
			{
				readnum ++;

				//printf("len = %d, %s\n", len, seq_data);
				//if(contianUnknownBase(seq_data)==NO && calcAverQual5End(qual_data)>=AVERAGE_QUAL_THRESHOLD_5End && calcAverQual3End(qual_data)>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data)<ARTIFACTS_BASE_A_THRESHOLD)
				if(containUnknownBase(seq_data[pairedID])==NO && calcAverQual5End(qual_data[pairedID])>=AVERAGE_QUAL_THRESHOLD_5End && qualSatisfied(qual_data[pairedID])==YES && calcAverQual3End(qual_data[pairedID])>=AVERAGE_QUAL_THRESHOLD_3End && getRatioBaseA(seq_data[pairedID])<ARTIFACTS_BASE_A_THRESHOLD)
				{ //不包含‘N’, 并且相应碱基的平均质量大于15

					//add the read
					if(addRead(seq_data[pairedID], readnum, graph)==FAILED)
					{
						printf("line=%d, In %s(), cannot add a read to graph, error!\n", __LINE__, __func__);
						releaseGraph(graph);
						return NULL;
					}

				}

				if(pairedID==1)
					pairedID = 0;
				else
					pairedID ++;

				line_index = 0;
			}
		}

		fclose(srcfp);
		srcfp = NULL;
	}

	return graph;
}

/**
 * Add a read to graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addRead(char *seq, uint64_t rid, graphtype *graph)
{
	int i, j, baseInt, pos;
	uint64_t hashcode;

	// generate the kmer integer sequence
	if(generateKmerSeqInt(kmerSeqInt, seq)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the hashcode and count the kmer occurrence
	hashcode = kmerhashInt(kmerSeqInt);
	if(addKmer(hashcode, kmerSeqInt, rid, 1, graph)==FAILED)
	{
		printf("line=%d, In %s(), cannot add kmer for read %lu, error!\n", __LINE__, __func__, rid);
		return FAILED;
	}

	pos = 2;
	for(i=kmerSize; i<readLen; i++, pos++)
	{
		// get the baseInt
		switch(seq[i])
		{
			case 'A':
			case 'a': baseInt = 0; break;
			case 'C':
			case 'c': baseInt = 1; break;
			case 'G':
			case 'g': baseInt = 2; break;
			case 'T':
			case 't': baseInt = 3; break;
			default: printf("line=%d, In %s(), base %c, error!\n", __LINE__, __func__, seq[i]); return FAILED;
		}

		// generate the kmer integer sequence
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				kmerSeqInt[j] = (kmerSeqInt[j] << 2) | (kmerSeqInt[j+1] >> 62);
			}
			kmerSeqInt[entriesPerKmer-2] = (kmerSeqInt[entriesPerKmer-2] << 2) | (kmerSeqInt[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		kmerSeqInt[entriesPerKmer-1] = ((kmerSeqInt[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

		// get the hashcode and count the kmer occurrence
		hashcode = kmerhashInt(kmerSeqInt);
		if(addKmer(hashcode, kmerSeqInt, rid, pos, graph)==FAILED)
		{
			printf("line=%d, In %s(), cannot add kmer, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute the kmer hash code.
 *  @return:
 *    The hash code.
 */
uint64_t kmerhashInt(uint64_t *seqInt)
{
	//return *seqInt;

	uint64_t hashcode;
	int i, j;

	hashcode = 5381;
	for(i=0; i<entriesPerKmer-1; i++)
	{
		for(j=0; j<32; j++)
		{
			hashcode += (hashcode << 5) | ((seqInt[i] >> (62-2*j)) & 3);
		}
	}

	for(j=0; j<lastEntryBaseNum; j++)
	{
		hashcode += (hashcode << 5) | ((seqInt[entriesPerKmer-1] >> (2*(lastEntryBaseNum-j-1))) & 3);
	}

	//return (hashcode & 0x7FFFFFFF) % hashTableSize;
	return hashcode % hashTableSize;

}

/**
 * 初始化de Bruijn图.
 *
 *   @ return:
 *     成功, 返回成功标记; 失败, 返回失败标记.
 */
graphtype *initgraph()
{
	//countMemory = 0;
	graphtype *graph = (graphtype *) malloc(sizeof(graph));
	if(!graph)
	{
		printf("In initgraph(), cannot malloc the graph, Error!\n");
		return NULL;
	}
	//countMemory += sizeof(graphtype);

	graph->pkmers = ( kmertype**)calloc(hashTableSize, sizeof( kmertype* ));
	//countMemory += sizeof( kmertype* ) * TABLE_SIZE_DE_BRUIJN;
	if( graph->pkmers == NULL )
	{
		printf("In initgraph(), cannot malloc the graph, Error!\n");
		free(graph); graph = NULL;
		return NULL;
	}

	return graph;
}

/**
 * 统计kmer的数目。
 *
 *  当de Bruijn图中不存在有kmer时，需要开辟一个新的kmer节点，并赋给初值；
 *  当de Bruijn图中已经存在有kmer时，只需增加其arraysize字段。
 */
short countKmer(uint64_t hashcode, uint64_t *kmerSeqInt, graphtype *graph)
{
	uint64_t *tmpSeq;
	kmertype *kmer;

	//unsigned int hashcode = kmerhash(str);
	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{ // de Bruijn图中不存在该kmer, 开辟kmertype节点
		kmer = ( kmertype* )malloc( sizeof( kmertype ) );
		if(!kmer)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

		kmer -> arraysize = 1;
		kmer -> multiplicity = 0;
		kmer -> ppos = NULL;

		tmpSeq = (uint64_t*) malloc(entriesPerKmer * sizeof(uint64_t));
		if(tmpSeq==NULL)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(memcpy(tmpSeq, kmerSeqInt, entriesPerKmer * sizeof(uint64_t))==NULL)
		{
			printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		kmer -> kmerseq = tmpSeq;

		kmer->next = graph->pkmers[ hashcode ];
		graph->pkmers[ hashcode ] = kmer;

		//countMemory += sizeof( kmertype );
		newKmerNum ++;
	}else
	{ //该kmer已经存在, 则只需更新数量
		kmer ->arraysize ++;
	}

	return SUCCESSFUL;
}

/**
 * Add a kmer to graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short addKmer(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph)
{
	kmertype *kmer;

	//unsigned int hashcode = kmerhash(str);
	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(kmer) //kmer存在，添加位置信息
	{
		if(kmer->ppos==NULL)  //如果位置数组指针为空，则开辟数组
		{
			kmer->ppos = (ridpostype*)malloc(sizeof(ridpostype)*(kmer->arraysize));
			if(kmer->ppos==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			kmer->multiplicity = 0;
			//countMemory += sizeof( ridpostype )*(kmer->arraysize);
		}
		kmer->ppos[ kmer->multiplicity ].delsign = NO;
		kmer->ppos[ kmer->multiplicity ].pos = rpos;
		kmer->ppos[ kmer->multiplicity ].rid = rid;
		kmer->ppos[ kmer->multiplicity ].reserved = 0;
		kmer->multiplicity ++;

	}else
	{
		printf("line=%d, In %s(), can not get the kmer %s. Error and exit.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt));
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * 删除kmer中rid的位置记录, 只作标记.
 *    @ return:
 *      成功，返回SUCCESSFUL；失败，返回FAILED。
 */
short delKmer(char *str, uint64_t rid, unsigned short rpos, graphtype *graph)
{
	kmertype *kmer;
	ridpostype *ridpostable = NULL;
	uint64_t hashcode;

	// generate the kmer integer sequence
	if(generateKmerSeqInt(kmerSeqInt, str)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	hashcode = kmerhashInt(kmerSeqInt);
	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{ //kmer不存在，删除失败。
		//printf("line=%d, In %s(), can not delete the kmer %s of (%d,%d), it does not exist in the graph. Error!\n", __LINE__, __func__, str, rid, rpos);
		return FAILED;
	}

	ridpostable = kmer->ppos;
	if(!ridpostable)
	{ //位置数组不存在，删除失败。
		printf("line=%d, In %s(), can not delete the kmer %s of (%lu, %u), the position array does not exist in the graph. Error!\n", __LINE__, __func__, str, rid, rpos);
		return FAILED;
	}

	if(kmer->multiplicity==0)
	{ //kmer已经被全部删除, 则打印出错信息 (该信息实际说明该(rid, rpos)不存在)
		//printf("line=%d, In %s(), can not delete the kmer %s of (%lu,%u), all of its positions have been deleted.\n", __LINE__, __func__, str, rid, rpos);
		return FAILED;
	}

	int index = findDirectIndex(rid, rpos, ridpostable, kmer->arraysize);//二分查找rid_pos_table
	if(index>=0)
	{  //找到该rid_pos

		if(ridpostable[index].delsign==0)
		{ //该kmer未被删除, 则将其删除
			ridpostable[index].delsign = 1;
			kmer->multiplicity--;
		}else //该read已经被删除, 则打印出错信息
		{
			printf("line=%d, In %s(), can not delete the kmer %s (%lu, %u), it has been deleted in the graph.\n", __LINE__, __func__, str, rid, rpos);
			return FAILED;
		}

	}else
	{  //该rid_pos不存在
		//printf("line=%d, In %s(), can not delete the kmer %s s(%d, %d), it does not exist in the graph.\n", __LINE__, __func__, str, rid, rpos);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 *
 */
/**
 * 删除kmer中rid的位置记录, 只作标记.
 *    @ return:
 *      成功，返回SUCCESSFUL；失败，返回FAILED。
 */
short delKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, unsigned short rpos, graphtype *graph)
{
	ridpostype *ridpostable = NULL;
	kmertype *kmer;

	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{ //kmer不存在，删除失败。
		//printf("line=%d, In %s(), can not delete the kmer %s of (%d,%d), it does not exist in the graph. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	ridpostable = kmer->ppos;
	if(!ridpostable)
	{ //位置数组不存在，删除失败。
		printf("line=%d, In %s(), can not delete the kmer %s of (%lu,%u), the position array does not exist in the graph. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	if(kmer->multiplicity==0)
	{ //kmer已经被全部删除, 则打印出错信息 (该信息实际说明该(rid, rpos)不存在)
		//printf("line=%d, In %s(), can not delete the kmer %s of (%lu, %u), all of its positions have been deleted.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	int index = findDirectIndex(rid, rpos, ridpostable, kmer->arraysize);//二分查找rid_pos_table
	if(index>=0)
	{  //找到该rid_pos

		if(ridpostable[index].delsign==0)
		{ //该kmer未被删除, 则将其删除
			ridpostable[index].delsign = 1;
			kmer->multiplicity--;
		}else //该read已经被删除, 则打印出错信息
		{
			printf("line=%d, In %s(), can not delete the kmer %s of (%lu, %u), it has been deleted in the graph.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
			return FAILED;
		}

	}else
	{  //该rid_pos不存在
		//printf("line=%d, In %s(), can not delete the kmer %s of (%lu, %u), it does not exist in the graph.\n", __LINE__, __func__,  getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * 恢复kmer中rid的位置记录, 也即恢复删除标记.
 *    @ return:
 *      成功，返回SUCCESSFUL；失败，返回FAILED。
 */
short recoverKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph)
{
	ridpostype *ridpostable;
	kmertype *kmer;

	hashcode = kmerhashInt(kmerSeqInt);
	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);
	if(!kmer)
	{ //kmer不存在，删除失败。
		//printf("line=%d, In %s(), can not recover the kmer %s of (%lu,%u), it does not exist in the graph, kmer=NULL. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	ridpostable = kmer->ppos;
	if(!ridpostable)
	{ //位置数组不存在，删除失败。
		printf("line=%d, In %s(), can not recover the kmer %s of (%lu,%u), the position array does not exist in the graph. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	if(kmer->multiplicity==kmer->arraysize)
	{
		//printf("line=%d, In %s(), can not recover the kmer %s of (%lu,%u), multiplicity==arraysize. Error!\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	int index = findDirectIndex(rid, rpos, ridpostable, kmer->arraysize);//二分查找rid_pos_table
	if(index>=0)
	{  //找到该rid_pos
		if(ridpostable[index].delsign==0)  //该read已经被恢复, 不能重复恢复
		{
			printf("line=%d, In %s(), can not recover the kmer %s of (%lu, %u), it has been recovered in the graph.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
			return FAILED;
		}
		ridpostable[index].delsign = 0;
		kmer->multiplicity++;

	}else
	{  //该rid_pos不存在
		//printf("line=%d, In %s(), can not recover the kmer %s of (%lu,%u), it does not exist in the graph.\n", __LINE__, __func__, getKmerBaseByInt(kmerSeqInt), rid, rpos);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * 直接二分查找ridpostable表，返回下表，失败返回-1.
 */
int findDirectIndex(uint64_t rid, uint16_t rpos, ridpostype *ridpostable, int posNum)
{
	int index = -1;
	int startIndex = findStartIndex(rid, ridpostable, posNum);//二分查找rid_pos_table
	if(startIndex>=0)
	{
		index = getExectIndex(rid, rpos, startIndex, ridpostable, posNum);
	}
	return index;
}

/**
 * 在ridpostable表中查找rid的起始位置, 不存在返回-1.
 */
int findStartIndex(uint64_t rid, ridpostype *rid_pos_table, int posNum)
{
	int left, right, middle, existFlag = 0;
	left = 0;
	right = posNum - 1;
	while(left<=right)
	{
		middle = (left+right)/2;
		if(rid==rid_pos_table[middle].rid)
		{
			//return middle;
			existFlag = 1;
			break;
		}
		if(rid>rid_pos_table[middle].rid)
			left = middle+1;
		else
			right = middle-1;
	}

	if(existFlag) //判断是否存在rid位置
	{
		while(middle>0 && rid_pos_table[middle-1].rid==rid)  //找到第一个rid的位置
			middle--;
		return middle; //返回该位置
	}

	return -1;
}

/**
 * 从起始位置开始查找rid和rpos, 并返回该位置，如果不存在返回-1.
 */
int getExectIndex(uint64_t rid, uint16_t rpos, int startIndex, ridpostype *ridpostable, int posNum)
{
	int i = startIndex, exectIndex = -1;
	while(i<posNum && ridpostable[i].rid==rid)
	{
		if(ridpostable[i].pos==rpos)
		{
			exectIndex = i;
			break;
		}
		//ridpostable[i].reserved = 1;  //将reversed标志位置为1
		i++;
	}

	return exectIndex;
}

/**
 * 按照kmer碱基序列, 取得kmer.
 */
inline kmertype *getKmerByBase(char *str, graphtype *graph)
{
	uint64_t hashcode;
	kmertype *kmer;

	// generate the kmer integer sequence
	if(generateKmerSeqInt(kmerSeqInt, str)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return NULL;
	}

	hashcode = kmerhashInt(kmerSeqInt);
	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);

	return kmer;
}

inline kmertype *getKmer(uint64_t *kmerSeqInt, graphtype *graph)
{
	uint64_t hashcode;

	hashcode = kmerhashInt(kmerSeqInt);

	return getKmerByHash(hashcode, kmerSeqInt, graph);
}

inline kmertype *getKmerByHash(uint64_t hashvalue, uint64_t *kmerSeqInt, graphtype *graph)
{
	kmertype *kmer;

	kmer = graph->pkmers[hashvalue];
	while(kmer)
	{
		if(identicalKmerSeq(kmerSeqInt, kmer->kmerseq)==YES)
			break;

		kmer = kmer->next;
	}

	return kmer;
}

/**
 * Check whether the two sequence is identical.
 *  @return:
 *  	If identical, return YES; otherwise, return NO.
 */
int identicalKmerSeq(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2)
{
	int i;
	for(i=0; i<entriesPerKmer; i++)
	{
		if(kmerSeqInt1[i] != kmerSeqInt2[i])
			return NO;
	}

	return YES;
}

/**
 * Get the kmer bases from integer.
 */
char *getKmerBaseByInt(uint64_t *kmerSeqInt)
{
	int i, j, k, baseInt;

	k = 0;
	for(i=0; i<entriesPerKmer-1; i++)
	{
		for(j=0; j<32; j++)
		{
			baseInt = (kmerSeqInt[i] >> (62-2*j)) & 3;
			switch(baseInt)
			{
				case 0: kmerBaseSeq[k]='A'; break;
				case 1: kmerBaseSeq[k]='C'; break;
				case 2: kmerBaseSeq[k]='G'; break;
				case 3: kmerBaseSeq[k]='T'; break;
			}
			k ++;
		}
	}

	for(j=0; j<lastEntryBaseNum; j++)
	{
		baseInt = (kmerSeqInt[entriesPerKmer-1] >> (2*(lastEntryBaseNum-j-1))) & 3;
		switch(baseInt)
		{
			case 0: kmerBaseSeq[k]='A'; break;
			case 1: kmerBaseSeq[k]='C'; break;
			case 2: kmerBaseSeq[k]='G'; break;
			case 3: kmerBaseSeq[k]='T'; break;
		}
		k ++;
	}

	kmerBaseSeq[k] = '\0';

	return kmerBaseSeq;
}


/**
 * Get the kmer integer sequence of a kmer by its kmer integer sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int getReverseKmerseqInt(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt)
{
	int i, j, baseInt;
	int entryIndexRev, baseNum;

	for(i=0; i<entriesPerKmer; i++) kmerSeqIntRev[i] = 0;

	entryIndexRev = baseNum = 0;
	for(j=0; j<lastEntryBaseNum; j++)
	{
		baseInt = (~(kmerSeqInt[entriesPerKmer-1] >> (j*2))) & 3;
		kmerSeqIntRev[entryIndexRev] = (kmerSeqIntRev[entryIndexRev] << 2) | baseInt;
		baseNum ++;
		if(baseNum==32)
		{
			entryIndexRev ++;
			baseNum = 0;
		}
	}
	for(i=entriesPerKmer-2; i>=0; i--)
	{
		for(j=0; j<32; j++)
		{
			baseInt = (~(kmerSeqInt[i] >> (2*j))) & 3;
			kmerSeqIntRev[entryIndexRev] = (kmerSeqIntRev[entryIndexRev] << 2) | baseInt;
			baseNum ++;
			if(baseNum==32)
			{
				entryIndexRev ++;
				baseNum = 0;
			}
		}
	}

	// ########################### Debug information ########################
	if(kmerSize<32 && baseNum!=kmerSize)
	{
		printf("line=%d, In %s(), baseNum=%d != kmerSize=%d, error!\n", __LINE__, __func__, baseNum, kmerSize);
		return FAILED;
	}
	// ########################### Debug information ########################

	return SUCCESSFUL;
}

/**
 * The fast version of getReverseKmer().
 */
kmertype *getReverseKmer(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt, graphtype *graph)
{
	uint64_t hashcode;

	if(getReverseKmerseqInt(kmerSeqIntRev, kmerSeqInt)==FAILED)
	{
		printf("line=%d, In %s(), cannot get the reverse kmer, error!\n", __LINE__, __func__);
		exit(1);
	}

	hashcode = kmerhashInt(kmerSeqIntRev);

	return getKmerByHash(hashcode, kmerSeqIntRev, graph);
}

/**
 * 对一个kmer的碱基序列进行反向互补操作.
 */
//void reverse(char *str)
//{
//	if(!str) { printf("In reverse(): the tuple is NULL! Error!\n"); exit(1); }
//	unsigned short i = 0;
//	char tmp[KMER_SIZE+1];
//	strncpy(tmp, str, KMER_SIZE);
//
//	while(i<KMER_SIZE)
//	{
//		switch(tmp[i])
//		{
//			case 'A': str[KMER_SIZE-i-1] = 'T'; break;
//			case 'C': str[KMER_SIZE-i-1] = 'G'; break;
//			case 'G': str[KMER_SIZE-i-1] = 'C'; break;
//			case 'T': str[KMER_SIZE-i-1] = 'A'; break;
//			case 'a': str[KMER_SIZE-i-1] = 't'; break;
//			case 'c': str[KMER_SIZE-i-1] = 'g'; break;
//			case 'g': str[KMER_SIZE-i-1] = 'c'; break;
//			case 't': str[KMER_SIZE-i-1] = 'a'; break;
//			default:
//				printf("In reverse(), unknown base %c, Error!\n", tmp[i]);
//				exit(1);
//		}
//		i++;
//	}
//}


/**
 *  Output the graph to file.
 *   File format:
 *   	(1) arraySize1 for kmer array, arraySize for ridposArray, and readLen, kmerSize, hashtableSize;
 *   	(2) kmerArray;
 *   	(3) ridposArray.
 *   	(4) kmerseqArray.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int outputGraphToFile(char *graphFile, graphtype *graph)
{
	struct kmerArrayNode{
		uint64_t hashCode;
		uint64_t arraysize;
	};
	struct ridposArrayNode{
		uint64_t rid:48;
		uint64_t rpos:16;
	};

	uint64_t kmerArraySize, ridposArraySize, tmp[6];
	uint64_t i, j, rowsNumKmerArray, rowsNumRidposArray, itemNumKmerSeq;
	struct kmerArrayNode kmerArray;
	struct ridposArrayNode ridposArray;
	kmertype *kmer;
	ridpostype *ridpos;
	FILE *fpGraph;

	fpGraph = fopen(graphFile, "wb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), can not open the file [ %s ], Error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}

	kmerArraySize = ridposArraySize = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];
		while(kmer)
		{
			kmerArraySize ++;
			ridposArraySize += kmer->arraysize;

			kmer = kmer->next;
		}
	}

	tmp[0] = kmerArraySize;
	tmp[1] = ridposArraySize;
	tmp[2] = readLen;
	tmp[3] = kmerSize;
	tmp[4] = hashTableSize;
	tmp[5] = pairedMode;

	if(fwrite(tmp, sizeof(uint64_t), 6, fpGraph)!=6)
	{
		printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
		return FAILED;
	}

	rowsNumKmerArray = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];
		while(kmer)
		{
			kmerArray.hashCode = i;
			kmerArray.arraysize = kmer->arraysize;
			rowsNumKmerArray ++;

			if(fwrite(&kmerArray, sizeof(struct kmerArrayNode), 1, fpGraph)!=1)
			{
				printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
				return FAILED;
			}

			kmer = kmer->next;
		}
	}

	rowsNumRidposArray = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];
		while(kmer)
		{
			ridpos = kmer->ppos;
			for(j=0; j<kmer->arraysize; j++)
			{
				ridposArray.rid = ridpos[j].rid;
				ridposArray.rpos = ridpos[j].pos;
				rowsNumRidposArray ++;

				if(fwrite(&ridposArray, sizeof(struct ridposArrayNode), 1, fpGraph)!=1)
				{
					printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			kmer = kmer->next;
		}
	}

	itemNumKmerSeq = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];
		while(kmer)
		{
			if(fwrite(kmer->kmerseq, sizeof(uint64_t), entriesPerKmer, fpGraph)!=entriesPerKmer)
			{
				printf("line=%d, In %s(), fwrite Error!\n", __LINE__, __func__);
				return FAILED;
			}
			itemNumKmerSeq ++;

			kmer = kmer->next;
		}
	}

	// ############################ Debug information ##############################
	if(kmerArraySize!=rowsNumKmerArray)
	{
		printf("line=%d, In %s(),kmerArraySize=%lu != rowsNumKmerArray=%lu. Error\n", __LINE__, __func__, kmerArraySize, rowsNumKmerArray);
		return FAILED;
	}
	if(itemNumKmerSeq!=rowsNumKmerArray)
	{
		printf("line=%d, In %s(),itemNumKmerSeq=%lu != rowsNumKmerArray=%lu. Error\n", __LINE__, __func__, itemNumKmerSeq, rowsNumKmerArray);
		return FAILED;
	}
	if(ridposArraySize!=rowsNumRidposArray)
	{
		printf("line=%d, In %s(),ridposArraySize=%lu != rowsNumRidposArray=%lu. Error\n", __LINE__, __func__, ridposArraySize, rowsNumRidposArray);
		return FAILED;
	}
	// ############################ Debug information ##############################

	fclose(fpGraph);
	fpGraph = NULL;

	return SUCCESSFUL;
}

/**
 *  Load the graph to memory.
 *   File format:
 *   	(1) arraySize1 for kmer array, arraySize for ridposArray, and readLen, kmerSize, hashtableSize, pairedMode;
 *   	(2) kmerArray;
 *   	(3) ridposArray.
 *   	(4) kmerseqArray.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int loadGraph(graphtype **graph, char *graphFile)
{
	struct kmerArrayNode{
		uint64_t hashCode;
		uint64_t arraysize;
	};
	struct ridposArrayNode{
		uint64_t rid:48;
		uint64_t rpos:16;
	};

	FILE *fpGraph;
	uint64_t kmerArraySize, ridposArraySize, tmp[6], hashcode, arraysize;
	uint64_t preHashCode;
	kmertype *preKmer;
	uint64_t i, j;
	struct kmerArrayNode kmerArray;
	struct ridposArrayNode ridposArray;
	kmertype *kmer;
	ridpostype *ridpos;

	printf("Begin loading the graph ...\n");

	fpGraph = fopen(graphFile, "rb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, graphFile);
		return FAILED;
	}

	// get the item number of the kmer array and ridpos array, respectively
	if(fread(tmp, sizeof(uint64_t), 6, fpGraph)!=6)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	kmerArraySize = tmp[0];
	ridposArraySize = tmp[1];
	readLen = tmp[2];
	kmerSize = tmp[3];
	hashTableSize = tmp[4];
	pairedMode = tmp[5];

	// initialize the graph
	*graph = initgraph();
	if(*graph==NULL)
	{
		printf("line=%d, In %s(), init graph error!\n", __LINE__, __func__);
		return FAILED;
	}

	// fill the data
	preHashCode = -1;
	preKmer = NULL;
	for(i=0; i<kmerArraySize; i++)
	{
		if(fread(&kmerArray, sizeof(struct kmerArrayNode), 1, fpGraph)!=1)
		{
			printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
			return FAILED;
		}

		hashcode = kmerArray.hashCode;
		kmer = (kmertype*) malloc(sizeof(kmertype));
		if(!kmer)
		{
			printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
			return FAILED;
		}
		kmer->arraysize = kmer->multiplicity = kmerArray.arraysize;
		kmer->ppos = NULL;
		kmer->kmerseq = NULL;
		kmer->next = NULL;

		if(hashcode==preHashCode)
		{ // they are in the same list
			preKmer->next = kmer;
		}else
		{ // the first node in the list
			preHashCode = hashcode;
			(*graph)->pkmers[hashcode] = kmer;
		}
		//preHashCode = hashcode;
		preKmer = kmer;
	}

	for(i=0; i<hashTableSize; i++)
	{
		kmer = (*graph)->pkmers[i];
		while(kmer)
		{
			arraysize = kmer->arraysize;
			ridpos = (ridpostype*) malloc(arraysize * sizeof(ridpostype));
			if(!ridpos)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			for(j=0; j<arraysize; j++)
			{
				if(fread(&ridposArray, sizeof(struct ridposArrayNode), 1, fpGraph)!=1)
				{
					printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
					return FAILED;
				}

				ridpos[j].rid = ridposArray.rid;
				ridpos[j].pos = ridposArray.rpos;
				ridpos[j].delsign = NO;
				ridpos[j].reserved = 0;
			}
			kmer->ppos = ridpos;

			kmer = kmer->next;
		}
	}

	for(i=0; i<hashTableSize; i++)
	{
		kmer = (*graph)->pkmers[i];
		while(kmer)
		{
			kmer->kmerseq = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
			if(!kmer->kmerseq)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(fread(kmer->kmerseq, sizeof(uint64_t), entriesPerKmer, fpGraph)!=entriesPerKmer)
			{
				printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
				return FAILED;
			}

			kmer = kmer->next;
		}
	}

	fclose(fpGraph);
	fpGraph = NULL;

	printf("End loading the graph.\n");

	return SUCCESSFUL;
}

/**
 * Get the global parameters from graph.
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
short GlobalParasFromGraph(int *readLenPara, int *kmerSizePara, uint64_t *hashTableSizePara, int *pairedModePara, char *graphFileName)
{
	FILE *fpGraph;
	uint64_t tmp[6];

	fpGraph = fopen(graphFileName, "rb");
	if(fpGraph==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, graphFileName);
		return FAILED;
	}

	// get the item number of the kmer array and ridpos array, respectively
	if(fread(tmp, sizeof(uint64_t), 6, fpGraph)!=6)
	{
		printf("line=%d, In %s(), fread Error!\n", __LINE__, __func__);
		return FAILED;
	}

	*readLenPara = tmp[2];
	*kmerSizePara = tmp[3];
	*hashTableSizePara = tmp[4];
	*pairedModePara = tmp[5];

	fclose(fpGraph);
	fpGraph = NULL;

	return SUCCESSFUL;
}

/**
 * Reset the de Bruijn graph.
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
int resetGraph(graphtype *graph)
{
	int64_t i, j, arraysize;
	kmertype *kmer;
	ridpostype *ridposArray;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];
		while(kmer)
		{
			arraysize = kmer->arraysize;
			ridposArray = kmer->ppos;

			kmer->multiplicity = arraysize;
			for(j=0; j<arraysize; j++)
			{
				ridposArray[j].delsign = NO;
				ridposArray[j].reserved = NO;
			}

			kmer = kmer->next;
		}
	}

	return SUCCESSFUL;
}
