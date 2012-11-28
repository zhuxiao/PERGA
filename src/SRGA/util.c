/*
 * util.c
 *
 *  Created on: May 29, 2010
 *      Author: zhuxiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Output all contig nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContig(contigtype *contighead)
{
	contigtype *contig;
	successRead_t *ridposorientation = NULL;
	int i = 0, j = 0, num = 0;
	char tmp_base;
	contig = contighead;
	while(contig)
	{
		switch(contig->base)
		{
			case 0: tmp_base = 'A'; break;
			case 1: tmp_base = 'C'; break;
			case 2: tmp_base = 'G'; break;
			case 3: tmp_base = 'T'; break;
			default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contig->base); return FAILED;
		}

		printf("%d\t%c\t%d", contig->index, tmp_base, contig->ridposnum);
		ridposorientation = contig->pridposorientation;
		num = contig->ridposnum;
		for(i=0; i<num; i++)
		{
			printf("\t(%lu,%d,%d,%c)", ridposorientation[i].rid, ridposorientation[i].pos, ridposorientation[i].matchnum, ridposorientation[i].orientation);
		}
		printf("\n");
		contig = contig->next;
		j++;
	}
	printf("There are %d nodes.\n", j);

	return SUCCESSFUL;
}

/**
 * Output the contig nodes near 5' end nodeNum bases.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigEnd5(contigtype *contighead, int nodeNum)
{
	contigtype *contig;
	successRead_t *ridposorientation = NULL;
	int i, j, num;
	char tmp_base;

	i = 0;
	j = 0;
	num = 0;
	contig = contighead;
	while(contig)
	{
		switch(contig->base)
		{
			case 0: tmp_base = 'A'; break;
			case 1: tmp_base = 'C'; break;
			case 2: tmp_base = 'G'; break;
			case 3: tmp_base = 'T'; break;
			default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contig->base); return FAILED;
		}

		printf("%d\t%c\t%d", contig->index, tmp_base, contig->ridposnum);
		ridposorientation = contig->pridposorientation;
		num = contig->ridposnum;
		for(i=0; i<num; i++)
		{
			printf("\t(%lu,%d,%d,%c)", ridposorientation[i].rid, ridposorientation[i].pos, ridposorientation[i].matchnum, ridposorientation[i].orientation);
		}
		printf("\n");
		contig = contig->next;
		j++;
		if(j>=nodeNum)
			break;
	}
	printf("There are %d nodes.\n", j);

	return SUCCESSFUL;
}

/**
 * Output the contig nodes near 3' end nodeNum bases.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigEnd3(contigtype *contighead, int nodeNum)
{
	contigtype *contig;
	successRead_t *ridposorientation = NULL;
	int i, j, num, startID, totalNodeNum = 0;
	char tmp_base;

	contig = contighead;
	while(contig)
	{
		totalNodeNum ++;
		contig = contig->next;
	}
	startID = totalNodeNum - nodeNum + 1;

	contig = contighead;
	while(contig)
	{
		if(contig->index==startID)
			break;
		contig = contig->next;
	}

	i = 0;
	j = 0;
	num = 0;
	while(contig)
	{
		switch(contig->base)
		{
			case 0: tmp_base = 'A'; break;
			case 1: tmp_base = 'C'; break;
			case 2: tmp_base = 'G'; break;
			case 3: tmp_base = 'T'; break;
			default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contig->base); return FAILED;
		}

		printf("%d\t%c\t%d", contig->index, tmp_base, contig->ridposnum);
		ridposorientation = contig->pridposorientation;
		num = contig->ridposnum;
		for(i=0; i<num; i++)
		{
			printf("\t(%lu,%d,%d,%c)", ridposorientation[i].rid, ridposorientation[i].pos, ridposorientation[i].matchnum, ridposorientation[i].orientation);
		}
		printf("\n");
		contig = contig->next;
		j++;
	}
	printf("There are %d nodes.\n", j);

	return SUCCESSFUL;
}

/*
 * 输出未被删除的kmer及其pos.
 */
short outputUndelKmerpos(graphtype *graph)
{
	uint64_t i, count;
	int posNum, j;
	kmertype *kmer;
	ridpostype *ridpostable;

	printf("Begin checking undelKmerpos:\n");

	count = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];  //取得kmer
		if(kmer)
		{
			posNum = kmer->arraysize;
			ridpostable = kmer->ppos;
			if(kmer->multiplicity>0)
			{
				for(j=0; j<posNum; j++)
				{
					if(ridpostable->delsign==0)  //未被删除，则输出该rid和pos
					{
						//printf("(%d,%d) ", ridpostable->rid, ridpostable->pos);
						count++;
					}
					ridpostable++;
				}
				//printf("\n");
			}
		}
	}
	printf("checking undelKmerpos finished, the count=%lu\n", count);

	return SUCCESSFUL;
}

/*
 * Output the remained Kmer and its multiplicity and arraySize.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputRemainedKmers(graphtype *graph)
{
	FILE *fpRemainedKmer;
	int64_t i, count;
	int posNum, j;
	kmertype *kmer;
	ridpostype *ridpostable;

	printf("Begin outputting  remained k-mers:\n");

	fpRemainedKmer = fopen("../remainedKmers.txt", "w");
	if(fpRemainedKmer==NULL)
	{
		printf("line=%d, In %s(), cannot open file [../remainedKmers.txt], error!\n", __LINE__, __func__);
		return FAILED;
	}


	count = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];  //取得kmer
		while(kmer)
		{
			if(kmer->multiplicity>0)
				fprintf(fpRemainedKmer, "%ld\t\t%u\t%u\n", i, kmer->multiplicity, kmer->arraysize);

//			posNum = kmer->arraysize;
//			ridpostable = kmer->ppos;
//			if(kmer->multiplicity>0)
//			{
//				for(j=0; j<posNum; j++)
//				{
//					if(ridpostable->delsign==0)  //未被删除，则输出该rid和pos
//					{
//						//printf("(%d,%d) ", ridpostable->rid, ridpostable->pos);
//						count++;
//					}
//					ridpostable++;
//				}
//				//printf("\n");
//			}

			kmer = kmer->next;
		}
	}

	fclose(fpRemainedKmer);

	printf("checking undelKmerpos finished, the count=%lu\n", count);

	return SUCCESSFUL;
}

/**
 * Release the De Bruijn graph.
 */
short releaseGraph(graphtype *graph)
{
	uint64_t i;
	kmertype *kmer, *head;

	for(i=0; i<hashTableSize; i++)
	{
		head = graph->pkmers[i];  //取得kmer
		while(head)
		{
			kmer = head->next;
			free(head->ppos);
			free(head->kmerseq);
			free(head);
			head = kmer;
		}
		graph->pkmers[i] = NULL;
	}

	free(graph->pkmers);

	free(graph);

	return SUCCESSFUL;
}

/**
 * 输出决策表中的所有的reads.
 */
short outputReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum)
{
	assemblingreadtype *this_assemblingRead = decisionTable;
	int i = 0, lockedNum = 0, appearNum = 0;
	for(; i<readsNum; i++)
	{
		printf("rid=%lu, firstpos=%d, orientation=%c, status=%d, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, locked=%d, matedFlag=%d\n",
				(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
				this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
				this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->locked, this_assemblingRead->matedFlag);

		if(this_assemblingRead->locked==1)
			lockedNum++;

		if(this_assemblingRead->lastpos>0)
			appearNum ++;

		this_assemblingRead++;
	}
	printf("The readsNum=%d, appearNum=%d, locked=%d\n", readsNum, appearNum, lockedNum);
	return SUCCESSFUL;
}

short outputReadsInDecisionTableToFile(assemblingreadtype *decisionTable, int readsNum, int contigID, int contigIndex)
{
	FILE *fpReads;
	assemblingreadtype *this_assemblingRead = decisionTable;
	int i, lockedNum, appearNum, failedNum;
	char fileName[256];

	sprintf(fileName, "tmp/readsInDT_%d_%d.txt", contigID, contigIndex);

	fpReads = fopen(fileName, "w");
	if(fpReads==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, "../readsInDT.txt");
		return FAILED;
	}

	lockedNum = 0;
	appearNum = 0;
	for(i=0; i<readsNum; i++)
	{
		fprintf(fpReads, "rid=%lu, firstpos=%d, orientation=%c, status=%d, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, locked=%d, matedFlag=%d, status=%d\n",
				(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
				this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
				this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->locked, this_assemblingRead->matedFlag, this_assemblingRead->status);

		if(this_assemblingRead->locked==1)
			lockedNum++;

		if(this_assemblingRead->lastpos>0)
			appearNum ++;

		this_assemblingRead++;
	}
	fprintf(fpReads, "readsNum=%d, appearNum=%d, locked=%d\n", readsNum, appearNum, lockedNum);

	failedNum = 0;
	fprintf(fpReads, "\n\nfailed reads:\n");
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->status==FAILED_STATUS)
		{
			fprintf(fpReads, "rid=%lu, firstpos=%d, orientation=%c, status=%d, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, locked=%d, matedFlag=%d, status=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
					this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
					this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->locked, this_assemblingRead->matedFlag, this_assemblingRead->status);

			failedNum ++;
		}

		this_assemblingRead ++;
	}
	fprintf(fpReads, "failedNum=%d\n", failedNum);

	fclose(fpReads);
	fpReads = NULL;

	return SUCCESSFUL;
}

short outputFailedReadsInDecisionTable(assemblingreadtype *decisionTable, int itemNumDecisionTable, int contigID, int contigIndex)
{
	assemblingreadtype *this_assemblingRead = decisionTable;
	int i, failedNum;

	failedNum = 0;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		if(this_assemblingRead->status==FAILED_STATUS)
		{
			printf("\trid=%lu, firstpos=%d, orientation=%c, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, locked=%d, matedFlag=%d, status=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation,
					this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
					this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->locked, this_assemblingRead->matedFlag, this_assemblingRead->status);

			failedNum ++;
		}

		this_assemblingRead ++;
	}
	if(failedNum>0)
		printf("contigID=%d, contigIndex=%d, failedNum=%d\n", contigID, contigIndex, failedNum);

	return SUCCESSFUL;
}

/**
 * 输出决策表中的锁定的reads.
 */
short outputLockedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum)
{
	if(lockedReadsNum==0)
	{ //没有锁定的reads, 直接返回.
		return SUCCESSFUL;
	}
	assemblingreadtype *this_assemblingRead = decisionTable;
	printf("Output locked assemblingreads:\n");
	int i, lockedNum = 0;
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->locked==1)
		{ //如果该read被锁定, 输出该read
			printf("rid=%lu, firstpos=%d, orientation=%c, status=%d, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, matedFlag=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
					this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
					this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->matedFlag);
			lockedNum++;
		}

		this_assemblingRead++;
	}
	printf("readsNum=%d, LockedNum=%d\n", readsNum, lockedNum);
	return SUCCESSFUL;
}

/**
 * 输出决策表中的锁定的reads.
 */
short outputMatedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum)
{
	assemblingreadtype *this_assemblingRead = decisionTable;
	printf("Output mated reads in decision table:\n");
	int i, matedNum = 0;
	for(i=0; i<readsNum; i++)
	{
		if(this_assemblingRead->matedFlag==YES)
		{ //如果该read被锁定, 输出该read
			printf("rid=%lu, firstpos=%d, orientation=%c, status=%d, kmerappeartimes=%d, kmerunappearblocks=%d, kmerunappeartimes=%d, lastappearpos=%d, lastpos=%d, matedFlag=%d\n",
					(uint64_t)this_assemblingRead->rid, this_assemblingRead->firstpos, this_assemblingRead->orientation, this_assemblingRead->status,
					this_assemblingRead->kmerappeartimes, this_assemblingRead->kmerunappearblocks, this_assemblingRead->kmerunappeartimes,
					this_assemblingRead->lastappearpos,	this_assemblingRead->lastpos, this_assemblingRead->matedFlag);
			matedNum++;
		}

		this_assemblingRead ++;
	}
	printf("readsNum=%d, matedNum=%d\n", readsNum, matedNum);
	return SUCCESSFUL;
}

/**
 * 输出kmer中的内容.
 */
short outputKmer(graphtype *graph, int hashcode, uint64_t *kmerSeqInt)
{
	kmertype *kmer;

	kmer = getKmerByHash(hashcode, kmerSeqInt, graph);;
	if(!kmer)
	{
		printf("In outputKmer(), the kmer==NULL.\n");
		return FAILED;
	}
	printf("kmerseq=%s, multi=%u, arraysize=%u\n", getKmerBaseByInt(kmerSeqInt), kmer->multiplicity, kmer->arraysize);
	ridpostype *rid_pos = kmer->ppos;
	int i = 0, posNum = kmer->arraysize;
	for(; i<posNum; i++)  //输出ridpos表
	{
		printf("\trid=%lu, pos=%u, delsign=%u, reserved=%u\n", (uint64_t)rid_pos->rid, rid_pos->pos, rid_pos->delsign, rid_pos->reserved);
		rid_pos++;
	}
	return SUCCESSFUL;
}

/**
 * 输出ridpos表中的内容.
 */
short outputRidpos(ridpostype *ridpos, int posNum)
{
	printf("arraysize=%d\n", posNum);
	int i = 0;
	for(; i<posNum; i++)  //输出ridpos表
	{
		printf("\trid=%lu, pos=%u, delsign=%u, reserved=%u\n", (uint64_t)ridpos->rid, ridpos->pos, ridpos->delsign, ridpos->reserved);
		ridpos++;
	}
	return SUCCESSFUL;
}

/**
 * 输出成功的reads信息.
 */
void outputSuccessReads(successRead_t *successReadArray, int successReadNum)
{
	int i = 0;
	for(; i<successReadNum; i++)
	{
		printf("Success Reads[%d]: rid=%lu, pos=%u, matchnum=%u, orientation=%c\n",
				i+1, successReadArray[i].rid, successReadArray[i].pos, successReadArray[i].matchnum, successReadArray[i].orientation);
	}
	printf("There are %d success reads\n", i);
}


short checkGraph(graphtype *graph)
{
	printf("Begin checking De Bruijn graph, please wait ...\n");

	uint64_t i, j, hashcode;
	kmertype *kmer;

	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];
		while(kmer)
		{
			hashcode = kmerhashInt(kmer->kmerseq);
			if(hashcode!=i)
			{
				printf("line=%d, In %s(), hashcode=%lu, i=%lu, error!\n", __LINE__, __func__, hashcode, i);
				return FAILED;
			}

			if(kmer->arraysize!=kmer->multiplicity)
			{
				printf("line=%d, In %s(), arraysize=%u != multiplicity=%u, error!\n", __LINE__, __func__, kmer->arraysize, kmer->multiplicity);
				return FAILED;
			}

			for(j=0; j<kmer->arraysize-1; j++)
			{
				if(kmer->ppos[j].rid>kmer->ppos[j+1].rid)
				{
					printf("line=%d, In %s(), rID Error.\n", __LINE__, __func__);
					return FAILED;
				}
			}

			kmer = kmer->next;
		}
	}
	printf("End checking the graph, congratulations.\n");

	return SUCCESSFUL;
}


short outputContigToTmpFile(contigtype *contighead, int outFileType)
{
	FILE *fpContigTmp;
	contigtype *contig;
	successRead_t *ridposorientation = NULL;
	int i, contigNodeNum, num;
	char tmp_base, fileName[256];

	if(outFileType==BASE_TYPE_FASTA_CONTIG_FILE)
		strcpy(fileName, "../tmpContig.fa");
	else
		strcpy(fileName, "../tmpContig_hanging.fa");

	fpContigTmp = fopen(fileName, "w");
	if(fpContigTmp==NULL)
	{
		printf("line=%d, In %s(), can not open file [%s].\n", __LINE__, __func__, fileName);
		return FAILED;
	}

	contigNodeNum = 0;
	contig = contighead;
	while(contig)
	{
		contigNodeNum ++;
		contig = contig->next;
	}

	if(outFileType==BASE_TYPE_FASTA_CONTIG_FILE)
	{
		fprintf(fpContigTmp, ">%d length: %d\n", 1, contigNodeNum); //输出标题信息
		contig = contighead;
		while(contig)
		{
			switch(contig->base)
			{
				case 0: tmp_base = 'A'; break;
				case 1: tmp_base = 'C'; break;
				case 2: tmp_base = 'G'; break;
				case 3: tmp_base = 'T'; break;
				default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contig->base); return FAILED;
			}

			fprintf(fpContigTmp, "%c", tmp_base); //输出碱基
			contig = contig->next;
			//contigID++;
		}
		fprintf(fpContigTmp, "\n");
	}else
	{
		fprintf(fpContigTmp, ">%d length: %d\n", 1, contigNodeNum); //输出标题信息

		contig = contighead;
		while(contig)
		{
			switch(contig->base)
			{
				case 0: tmp_base = 'A'; break;
				case 1: tmp_base = 'C'; break;
				case 2: tmp_base = 'G'; break;
				case 3: tmp_base = 'T'; break;
				default: printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contig->base); return FAILED;
			}

			fprintf(fpContigTmp, "%d\t%c\t%d", contig->index, tmp_base, contig->ridposnum); //输出碱基
			ridposorientation = contig->pridposorientation;
			num = contig->ridposnum;
			for(i=0; i<num; i++)
			{
				fprintf(fpContigTmp, "\t(%lu,%d,%d,%c)",
						ridposorientation[i].rid,
						ridposorientation[i].pos,
						ridposorientation[i].matchnum,
						ridposorientation[i].orientation);
			}
			fprintf(fpContigTmp, "\n");

			contig = contig->next;
			//contigID++;
		}
	}

	fclose(fpContigTmp);
	fpContigTmp = NULL;

	return SUCCESSFUL;
}


short outputPEHashArray(PERead_t **PEHashArray)
{
	int i, totalReadNum;
	PERead_t *tmpRead;

	totalReadNum = 0;
	for(i=0; i<TABLE_SIZE_HASH_PE; i++)
	{
		if(PEHashArray[i])
		{
			tmpRead = PEHashArray[i];
			while(tmpRead)
			{
				printf("[%d]: rid=%lu, cpos=%u, orient=%u\n", i, tmpRead->rid, tmpRead->cpos, tmpRead->orient);
				totalReadNum ++;
				tmpRead = tmpRead->next;
			}
		}
	}

	printf("totalReadNum=%d\n", totalReadNum);

	return SUCCESSFUL;
}


short checkReadListArr(readList_t *readListArray, int64_t itemNumInReadListArray)
{
	int64_t i;

	for(i=0; i<itemNumInReadListArray-1; i++)
	{
		if(readListArray[i].readID >= readListArray[i+1].readID)
		{
			printf("line=%d, In %s(), the data was not correctly ordered, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}


short outputReadListToFile(char *readListFile, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray)
{
	uint64_t i, j, matchNum, readID;
	FILE *fpReadList;
	readPos_t *pReadPos;

	fpReadList = fopen(readListFile, "w");
	if(fpReadList==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readListFile);
		return FAILED;
	}

	for(i=0; i<itemNumInReadListArray; i++)
	{
		readID = readListArray[i].readID;
		matchNum = readListArray[i].matchNum;
		pReadPos = readPosArray + readListArray[i].firstRow;
		for(j=0; j<matchNum; j++)
		{
			fprintf(fpReadList, "[%lu]:\t%u\t%u\t%c\n", readID, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, pReadPos[j].orientation);
		}
	}

	fclose(fpReadList);
	fpReadList = NULL;

	return SUCCESSFUL;
}


short outputMatedReadsInReadListToFile(char *readListFile, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray)
{
	uint64_t i, j, matchNum, readID;
	FILE *fpReadList;
	readPos_t *pReadPos;

	fpReadList = fopen(readListFile, "w");
	if(fpReadList==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readListFile);
		return FAILED;
	}

	i = 0;
	while(i<itemNumInReadListArray-1)
	{
		if(readListArray[i].readID % 2 == 1 && readListArray[i].readID+1==readListArray[i+1].readID)
		{
			if(readPosArray[readListArray[i].firstRow].matchBaseNum==readLen && readPosArray[readListArray[i+1].firstRow].matchBaseNum==readLen)
			{
				readID = readListArray[i].readID;
				matchNum = readListArray[i].matchNum;
				pReadPos = readPosArray + readListArray[i].firstRow;
				for(j=0; j<matchNum; j++)
				{
					fprintf(fpReadList, "[%lu]:\t%u\t%u\t%c\n", readID, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, pReadPos[j].orientation);
				}

				readID = readListArray[i+1].readID;
				matchNum = readListArray[i+1].matchNum;
				pReadPos = readPosArray + readListArray[i+1].firstRow;
				for(j=0; j<matchNum; j++)
				{
					fprintf(fpReadList, "[%lu]:\t%u\t%u\t%c\n", readID, pReadPos[j].contigPos, pReadPos[j].matchBaseNum, pReadPos[j].orientation);
				}
			}

			i += 2;
		}else
		{
			i ++;
		}
	}

	fclose(fpReadList);
	fpReadList = NULL;

	return SUCCESSFUL;
}


short convertFragmentSizeFileInText(const char *fragmentSizeFile)
{
	int32_t fragmentSize;
	char fileName[256];
	FILE *fpFragSize_In, *fpFragSize_Out;

	strcpy(fileName, fragmentSizeFile);
	strcat(fileName, ".txt");

	fpFragSize_In = fopen(fragmentSizeFile, "rb");
	if(fpFragSize_In==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fragmentSizeFile);
		return FAILED;
	}

	fpFragSize_Out = fopen(fileName, "w");
	if(fpFragSize_Out==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fileName);
		return FAILED;
	}

	while(1)
	{
		if(fread(&fragmentSize, sizeof(int32_t), 1, fpFragSize_In)!=1)
		{
			if(feof(fpFragSize_In))
			{
				break;
			}else
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		fprintf(fpFragSize_Out, "%d\n", fragmentSize);
	}

	fclose(fpFragSize_Out);
	fpFragSize_Out = NULL;
	fclose(fpFragSize_In);
	fpFragSize_Out = NULL;

	return SUCCESSFUL;
}

short outputReadPosInGraph(int64_t readID, graphtype *graph)
{
	int64_t i, j;
	kmertype *kmer;
	ridpostype *ridpostable;
	int posNum;

	for(i=0; i<hashTableSize; i++)
	{
		kmer = graph->pkmers[i];
		while(kmer)
		{
			ridpostable = kmer->ppos;
			posNum = kmer->arraysize;
			for(j=0; j<posNum; j++)
			{
				if(ridpostable[j].rid==readID)
				{
					printf("[%ld]: rid=%lu, pos=%lu, kmerSeq=%s\n", j, (uint64_t)ridpostable[j].rid, (uint64_t)ridpostable[j].pos, getKmerBaseByInt(kmer->kmerseq));
				}
			}

			kmer = kmer->next;
		}
	}

	return SUCCESSFUL;
}

/**
 * Output items in navigation occurrence queue.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputNaviOccQueue(double *naviOccQueuePara, int itemNumNaviOccQueuePara, int frontRowNaviOccQueuePara)
{
	int i, j;

	j = frontRowNaviOccQueuePara;
	for(i=0; i<itemNumNaviOccQueuePara; i++)
	{
		printf("naviOccQueue[%d]: %.2f\n", i, naviOccQueuePara[j]);
		j = (j+1) % maxItemNumNaviOccQueue;
	}

	return SUCCESSFUL;
}
