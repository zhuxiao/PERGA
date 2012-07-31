/*
 * hashPE.c
 *
 *  Created on: Dec 1, 2011
 *      Author: xiao
 */

#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Estimate insertSize and Sdev.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int estimateInsertSizeAndSdev()
{
	printf("Begin estimating the insert size and standard deviation of paired end fragment library ...\n");

	estContigArr = (estContig_t *) calloc(MAX_NUM_EST_CONTIG, sizeof(estContig_t));
	if(estContigArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(PEGivenType>=INSERT_PE_GIVEN_TYPE)
		minContigLenEst = meanSizeInsert * MIN_CONTIG_LEN_EST_FACTOR;
	else
		minContigLenEst = MIN_CONTIG_LEN_EST;

	if(initPEHashParas()==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the PEHash table parameters, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// build estimated contigs
	estimateSuccessFlag = buildEstContigs(sampleContigsFile);
	if(estimateSuccessFlag!=FAILED)
	{
		//printf("line=%d, In %s(), can not build estimated contigs. Error!\n", __LINE__, __func__);
		//return FAILED;
		// estimation
		if(meanSizeInsertAndSdevEstimation(fragmentSizeFile, estContigArr, contigNumEstContigArr)==FAILED)
		{
			printf("line=%d, In %s(), cannot estimate the insert size and standard deviation of library fragments, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{

	}

	// reset De Bruijn graph
	if(resetGraph(deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot reset the De Bruijn graph, error!\n", __LINE__, __func__);
		return FAILED;
	}

	int i;
	for(i=0; i<contigNumEstContigArr; i++)
	{
		releaseContig(estContigArr[i].contighead);
	}
	free(estContigArr);
	estContigArr = NULL;

	printf("End estimating the insert size and standard deviation of paired end fragment library.\n");

	return SUCCESSFUL;
}

/**
 * Initialzie the PEHash table parameters.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int initPEHashParas()
{
	if(PEGivenType==INSERT_PE_GIVEN_TYPE)
	{
		standardDev = meanSizeInsert * DRAFT_SDEV_FACTOR;
	}

	if(PEGivenType>NONE_PE_GIVEN_TYPE)
	{
		standardDevFactor = SDEV_FACTOR;  ///////////////////////////////////
		shiftLenRound1 = readLen;
		minContigLenUsingPE = minMarginLenPEHash = (int)(meanSizeInsert - standardDevFactor * standardDev - (readLen-kmerSize+1) + 1);
		//if(minContigLenUsingPE<KMER_SIZE)
		if(minContigLenUsingPE<readLen+kmerSize)
		{
			//minContigLenUsingPE = minMarginLenPEHash = KMER_SIZE;
			minContigLenUsingPE = minMarginLenPEHash = readLen + kmerSize;
		}
		maxMarginLenPEHash = (int)(meanSizeInsert + standardDevFactor * standardDev - (readLen-kmerSize+1));
		if(maxMarginLenPEHash<=minContigLenUsingPE)
		{
			printf("line=%d, In %s(), maxMarginLenPEHash=%d <= minContigLenUsingPE=%d, incorrect given insert size, error!\n", __LINE__, __func__, maxMarginLenPEHash, minContigLenUsingPE);
			return FAILED;
		}
		maxRegLenPEHash = maxMarginLenPEHash - minContigLenUsingPE + 1;
		minRegLenUsingPE = maxRegLenPEHash * REG_LEN_PE_HASH_FACTOR;
	}

	// ######################### Debug information ##########################
//	if(PEGivenType==INSERT_PE_GIVEN_TYPE)
//	{
//		printf("ReadLen=%d, meanSizeInsert=%.2f\n", readLen, meanSizeInsert);
//		printf("The standard deviation is set to be %.2f by default.\n", standardDev);
//	}else if(PEGivenType==BOTH_PE_GIVEN_TYPE)
//		printf("readLen=%d, meanSizeInsert=%.2f, standardDev=%.2f\n", readLen, meanSizeInsert, standardDev);
	// ######################### Debug information ##########################

	return SUCCESSFUL;
}

/**
 * Update the PE hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int updatePEHashTable(int contigNodesNum, int assemblyRound)
{
	contigtype *tmpContig;
	int i, newContigIndex, ridposnum;
	successRead_t *tmpRidposOrient;
	contigtype *tmpRegLeftContig, *tmpRegRightContig;

	if(allowedUpdatePEHashArrFlag==NO)
	{
		allowedUpdatePEHashArrFlag = YES;
		return SUCCESSFUL;
	}

	if(contigNodesNum>=maxMarginLenPEHash)
	{ // slide the window by one base pair
		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			tmpRegLeftContig = shiftedRegLeftContig;
			tmpRegRightContig = shiftedRegRightContig;
		}else
		{ // the second round
			tmpRegLeftContig = hashRegLeftContig;
			tmpRegRightContig = hashRegRightContig;
		}

		// remove the reads hanging on the left margin
		ridposnum = tmpRegLeftContig->ridposnum;
		if(ridposnum>0)
		{
			tmpRidposOrient = tmpRegLeftContig->pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==readLen)
				{
					if(delReadfromPEHashtable(tmpRidposOrient[i].rid)==FAILED)
					{
						printf("line=%d, In %s(), cannot remove read %lu from PE hash table, error!\n", __LINE__, __func__, tmpRidposOrient[i].rid);
						return FAILED;
					}
				}
			}
		}

		// add new reads
		if(tmpRegRightContig->next==NULL)
		{
			printf("line=%d, In %s(), pContig->next==NULL, error!\n", __LINE__, __func__);
			return FAILED;
		}
		tmpRegRightContig = tmpRegRightContig->next;

		ridposnum = tmpRegRightContig->ridposnum;
		if(ridposnum>0)
		{
			tmpRidposOrient = tmpRegRightContig->pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==readLen)
				{
					if(addReadToPEHashtable(tmpRidposOrient+i, tmpRegRightContig->index, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, tmpRidposOrient[i].rid);
						return FAILED;
					}
				}
			}
		}

		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			shiftedRegLeftContig = shiftedRegLeftContig->next;
			shiftedRegRightContig = shiftedRegRightContig->next;
		}else
		{ // the second round
			hashRegLeftContig = hashRegLeftContig->next;
			hashRegRightContig = hashRegRightContig->next;
		}
	}else if(contigNodesNum>=minMarginLenPEHash)
	{
		if(contigNodesNum==minMarginLenPEHash)
		{ // initialize the PE hash table
			if(readsNumInPEHashArr>0)
			{
				if(cleanReadsFromPEHashtable()==FAILED)
				{
					printf("line=%d, In %s(), cannot clean PE hash table, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			readsNumInPEHashArr = 0;
			regLenPEHash = 1;
			if(assemblyRound==FIRST_ROUND_ASSEMBLY)
			{ // the first round
				newContigIndex = shiftLenRound1 + 1;
				tmpContig = contighead;
				while(tmpContig)
				{
					if(tmpContig->index==newContigIndex)
						break;
					tmpContig = tmpContig->next;
				}
				shiftedRegLeftContig = shiftedRegRightContig = tmpContig;
				validReadOrientPEHash = ORIENTATION_MINUS;

				// add reads to PE hash table
				ridposnum = shiftedRegRightContig->ridposnum;
				if(ridposnum>0)
				{
					tmpRidposOrient = shiftedRegRightContig->pridposorientation;
					for(i=0; i<ridposnum; i++)
					{
						if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==readLen)
						{
							if(addReadToPEHashtable(tmpRidposOrient+i, shiftedRegRightContig->index, assemblyRound)==FAILED)
							{
								printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, tmpRidposOrient[i].rid);
								return FAILED;
							}
						}
					}
				}
			}else
			{ // the second round
				hashRegLeftContig = hashRegRightContig = contighead;
				validReadOrientPEHash = ORIENTATION_PLUS;

				// add reads to PE hash table
				ridposnum = hashRegRightContig->ridposnum;
				if(ridposnum>0)
				{
					tmpRidposOrient = hashRegRightContig->pridposorientation;
					for(i=0; i<ridposnum; i++)
					{
						if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==readLen)
						{
							if(addReadToPEHashtable(tmpRidposOrient+i, hashRegRightContig->index, assemblyRound)==FAILED)
							{
								printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, tmpRidposOrient[i].rid);
								return FAILED;
							}
						}
					}
				}
			}
		}else
		{ // enlarge the window by one base pair
			if(assemblyRound==FIRST_ROUND_ASSEMBLY)
			{ // the first round
				tmpRegRightContig = shiftedRegRightContig;
			}else
			{ // the second round
				tmpRegRightContig = hashRegRightContig;
			}

			regLenPEHash ++;

			// add new reads
			if(tmpRegRightContig->next==NULL)
			{
				printf("line=%d, In %s(), pContig->next==NULL, error!\n", __LINE__, __func__);
				return FAILED;
			}
			tmpRegRightContig = tmpRegRightContig->next;

			ridposnum = tmpRegRightContig->ridposnum;
			if(ridposnum>0)
			{
				tmpRidposOrient = tmpRegRightContig->pridposorientation;
				for(i=0; i<ridposnum; i++)
				{
					if(tmpRidposOrient[i].orientation==validReadOrientPEHash && tmpRidposOrient[i].matchnum==readLen)
					{
						if(addReadToPEHashtable(tmpRidposOrient+i, tmpRegRightContig->index, assemblyRound)==FAILED)
						{
							printf("line=%d, In %s(), cannot add read %lu to PE hash table, error!\n", __LINE__, __func__, tmpRidposOrient[i].rid);
							return FAILED;
						}
					}
				}
			}

			if(assemblyRound==FIRST_ROUND_ASSEMBLY)
			{ // the first round
				shiftedRegRightContig = shiftedRegRightContig->next;
			}else
			{ // the second round
				hashRegRightContig = hashRegRightContig->next;
			}
		}
	}else
	{ // do nothing

	}

	return SUCCESSFUL;
}

/**
 * Get the paired read by the given readID.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int getReadFromPEHashtable(PERead_t **pRead, uint64_t readID)
{
	uint64_t hashcode;

	hashcode = readID & RID_LOW_BITS_MASK;
	if(PEHashArr[hashcode])
	{
		*pRead = PEHashArr[hashcode];
		while(*pRead)
		{
			if((*pRead)->rid==readID)
			{
				break;
			}
			*pRead = (*pRead)->next;
		}
	}else
	{
		*pRead = NULL;
	}

	return SUCCESSFUL;
}

/**
 * Add a read into PE hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int addReadToPEHashtable(successRead_t *ridposOrient, int contigPos, int assemblyRound)
{
	uint64_t hashcode;
	PERead_t *tmpRead;

	// ######################## Debug information ############################
	//if(ridposOrient->rid==4882354)
	//{
	//	printf("rid=%d, matchnum=%d, orient=%c, firstpos=%d\n", ridposOrient->rid, ridposOrient->matchnum, ridposOrient->orientation, ridposOrient->pos);
	//}
	// ######################## Debug information ############################

	hashcode = ridposOrient->rid & RID_LOW_BITS_MASK;

	tmpRead = (PERead_t*) malloc(sizeof(PERead_t));
	if(!tmpRead)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first round
		tmpRead->cpos = contigPos - ridposOrient->matchnum + 1;
		if(tmpRead->cpos<=0)
		{
			printf("line=%d, In %s(), cpos=%d, error!\n", __LINE__, __func__, tmpRead->cpos);
			return FAILED;
		}
	}else
	{ // the second round
		tmpRead->cpos = contigPos;
	}

	tmpRead->rid = ridposOrient->rid;
	tmpRead->orient = ORIENTATION_PLUS;
	//insert the new node from head
	tmpRead->next = PEHashArr[hashcode];
	PEHashArr[hashcode] = tmpRead;

	readsNumInPEHashArr ++;

	return SUCCESSFUL;
}

/**
 * Remove a read into PE hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int delReadfromPEHashtable(uint64_t readID)
{
	uint64_t hashcode;
	PERead_t *tmpRead, *preRead;

	// ######################## Debug information ############################
	//if(readID==4882354)
	//{
	//	printf("readID=%lu\n", readID);
	//}
	// ######################## Debug information ############################

	hashcode = readID & RID_LOW_BITS_MASK;
	tmpRead = PEHashArr[hashcode];
	if(!tmpRead)
	{
		printf("line=%d, In %s(), cannot get the read %lu from PE hash table, error!\n", __LINE__, __func__, readID);
		return FAILED;
	}

	preRead = NULL;
	while(tmpRead)
	{
		if(tmpRead->rid==readID)
			break;

		preRead = tmpRead;
		tmpRead = tmpRead->next;
	}

	if(!tmpRead)
	{
		printf("line=%d, In %s(), cannot get read %lu from PE hash table, error!\n", __LINE__, __func__, readID);
		return FAILED;
	}

	if(preRead==NULL)
	{
		PEHashArr[hashcode] = tmpRead->next;
	}else
	{
		preRead->next = tmpRead->next;
	}
	free(tmpRead);

	readsNumInPEHashArr --;

	return SUCCESSFUL;
}

/**
 * Clean the PE hash table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int cleanReadsFromPEHashtable()
{
	int i;
	PERead_t *tmpRead, *head;

	//###################### Debug information ########################
	//printf("line=%d, In %s(), before clean, readsNumInPEHashArr=%d\n", __LINE__, __func__, readsNumInPEHashArr);
	//###################### Debug information ########################

	for(i=0; i<TABLE_SIZE_HASH_PE; i++)
	{
		if(PEHashArr[i])
		{
			head = PEHashArr[i];
			while(head)
			{
				tmpRead = head->next;
				free(head);
				head = tmpRead;
				readsNumInPEHashArr --;
			}
			PEHashArr[i] = NULL;
		}
	}

	//###################### Debug information ########################
	//printf("line=%d, In %s(), after clean, readsNumInPEHashArr=%d\n", __LINE__, __func__, readsNumInPEHashArr);
	//###################### Debug information ########################

	regLenPEHash = 0;

	// ######################### Debug information ##############################
	if(readsNumInPEHashArr!=0)
	{
		printf("line=%d, In %s(), readsNumInPEHashArr=%d != 0, error!\n", __LINE__, __func__, readsNumInPEHashArr);
		return FAILED;
	}
	// ######################### Debug information ##############################

	return SUCCESSFUL;
}

/**
 * Initialize the PE hash table before second round assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int initPEHashtableSecondAssembly(contigtype *contighead, int contigNodesNum)
{
	int i, tmpLeftContigIndex, tmpRightContigIndex;
	contigtype *tmpContig;
	int ridposnum;
	successRead_t *tmpRidposorient;

	if(cleanReadsFromPEHashtable()==FAILED)
	{
		printf("line=%d, In %s, cannot clean PE hash table, error!\n", __LINE__, __func__);
		return ERROR;
	}

	//###################### Debug information ########################
	if(readsNumInPEHashArr!=0)
	{
		printf("line=%d, In %s(), readsNumInPEHashArr=%d != 0, error!\n", __LINE__, __func__, readsNumInPEHashArr);
		return FAILED;
	}
	//###################### Debug information ########################

	if(contigNodesNum>=maxMarginLenPEHash)
	{ // full length of the hash region
		tmpRightContigIndex = contigNodesNum - minMarginLenPEHash + 1;
		tmpLeftContigIndex = contigNodesNum - maxMarginLenPEHash + 1;

		regLenPEHash = tmpRightContigIndex - tmpLeftContigIndex + 1;
		if(regLenPEHash!=maxRegLenPEHash)
		{
			printf("line=%d, In %s(), regLenPEHash=%d, != maxRegLenPEHash=%d, error!\n", __LINE__, __func__, regLenPEHash, maxRegLenPEHash);
			return FAILED;
		}

		tmpContig = contighead;
		while(tmpContig)
		{
			if(tmpContig->index==tmpLeftContigIndex)
				break;
			tmpContig = tmpContig->next;
		}
		hashRegLeftContig = tmpContig;

		while(tmpContig)
		{
			if(tmpContig->index==tmpRightContigIndex)
				break;
			tmpContig = tmpContig->next;
		}
		hashRegRightContig = tmpContig;
		validReadOrientPEHash = ORIENTATION_PLUS;

		// add reads into PE hash table
		tmpContig = hashRegLeftContig;
		while(tmpContig)
		{
			ridposnum = tmpContig->ridposnum;
			tmpRidposorient = tmpContig->pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				if(tmpRidposorient[i].orientation==validReadOrientPEHash && tmpRidposorient[i].matchnum==readLen)
				{
					if(addReadToPEHashtable(tmpRidposorient+i, tmpContig->index, SECOND_ROUND_ASSEMBLY)==FAILED)
					{
						printf("In %s(), cannot add read %lu to PE hash table, error!\n", __func__, tmpRidposorient[i].rid);
						return FAILED;
					}
				}
			}

			if(tmpContig==hashRegRightContig)
				break;

			tmpContig = tmpContig->next;
		}

		//########################## Debug information ###########################
		if(tmpContig==NULL)
		{
			printf("line=%d, In %s(), PE hash region error!\n", __LINE__, __func__);
			return FAILED;
		}
		//########################## Debug information ###########################

	}else if(contigNodesNum>=minMarginLenPEHash)
	{ // some reads in the hash region
		tmpRightContigIndex = contigNodesNum - minMarginLenPEHash + 1;
		tmpLeftContigIndex = 1;

		regLenPEHash = tmpRightContigIndex - tmpLeftContigIndex + 1;
		if(regLenPEHash<=0)
		{
			printf("line=%d, In %s(), regLenPEHash=%d, error!\n", __LINE__, __func__, regLenPEHash);
			return FAILED;
		}

		hashRegLeftContig = contighead;

		tmpContig = contighead;
		while(tmpContig)
		{
			if(tmpContig->index==tmpRightContigIndex)
				break;
			tmpContig = tmpContig->next;
		}
		hashRegRightContig = tmpContig;
		validReadOrientPEHash = ORIENTATION_PLUS;

		// add reads into PE hash table
		tmpContig = hashRegLeftContig;
		while(tmpContig)
		{
			ridposnum = tmpContig->ridposnum;
			tmpRidposorient = tmpContig->pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				if(tmpRidposorient[i].orientation==validReadOrientPEHash && tmpRidposorient[i].matchnum==readLen)
				{
					if(addReadToPEHashtable(tmpRidposorient+i, tmpContig->index, SECOND_ROUND_ASSEMBLY)==FAILED)
					{
						printf("In %s(), cannot add read %lu to PE hash table, error!\n", __func__, tmpRidposorient[i].rid);
						return FAILED;
					}
				}
			}

			if(tmpContig==hashRegRightContig)
				break;

			tmpContig = tmpContig->next;
		}

		//########################## Debug information ###########################
		if(tmpContig==NULL)
		{
			printf("line=%d, In %s(), hash region error!\n", __LINE__, __func__);
			return FAILED;
		}
		//########################## Debug information ###########################

	}else
	{ // do nothing

	}

	return SUCCESSFUL;
}

// ============================== Estimation ================================
/**
 * Estimate the insert size and standard deviation of library fragments.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int meanSizeInsertAndSdevEstimation(const char *fragmentSizeFile, estContig_t *estContigArray, int contigNumEstContigArray)
{
	int i;
	FILE *fpFragmentSize;

	if(contigNumEstContigArray<=0)
	{
		printf("line=%d, In %s(), contigNumEstContigArray=%d, cannot estimate the mean and sdev. of fragments in paired ends library, error!\n", __LINE__, __func__, contigNumEstContigArray);
		return FAILED;
	}

	fpFragmentSize = fopen(fragmentSizeFile, "wb");
	if(fpFragmentSize==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fragmentSizeFile);
		return FAILED;
	}

	for(i=0; i<contigNumEstContigArray; i++)
	//for(i=19; i<22; i++)
	{
		// process single contig
		if(getPairedEndsFromSingleContig(fpFragmentSize, estContigArray+i)==FAILED)
		{
			printf("line=%d, In %s(), cannot get paired reads from contig: %d, error!\n", __LINE__, __func__, estContigArray[i].contigID);
			return FAILED;
		}
	}

	fclose(fpFragmentSize);
	fpFragmentSize = NULL;

/*
	// ################################# Debug information ###############################
	if(convertFragmentSizeFileInText(fragmentSizeFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot convert fragment size file in text mode, error!\n", __LINE__, __func__);
		return FAILED;
	}
	// ################################# Debug information ###############################
*/

	// reload the file and estimate the insert size and standard deviation
	if(computeInsertSizeAndSdev(&meanSizeInsertEst, &standardDevEst, fragmentSizeFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot compute the insert size and standard deviation, error!\n", __LINE__, __func__);
		return FAILED;
	}
	printf("Number of estimated contigs: %d\n", contigNumEstContigArray);

	if(PEGivenType>NONE_PE_GIVEN_TYPE)
	{
		oldMeanSizeInsert = meanSizeInsert;
		oldStandardDev = standardDev;
	}

	meanSizeInsert = meanSizeInsertEst;
	standardDev = standardDevEst;

	oldPEGivenType = PEGivenType;
	PEGivenType = BOTH_PE_GIVEN_TYPE;

	return SUCCESSFUL;
}

/**
 * Get the paired ends from single contig, and output them to a binary file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int getPairedEndsFromSingleContig(FILE *fpFragSize, estContig_t *estContig)
{
	// initialize the memory
	if(initMemGetPESingleContig(estContig)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize memory for getting paired ends of contig %d, error!\n", __LINE__, __func__, estContig->contigID);
		return FAILED;
	}

	// fill data
	if(fillDataReadPosTmpArr(readPosTmpArr, estContig->contighead)==FAILED)
	{
		printf("line=%d, In %s(), cannot fill data of readPosTmp array of contig %d, error!\n", __LINE__, __func__, estContig->contigID);
		return FAILED;
	}

	// sort the data in readPosTmp array
	if(radixSortReadPosTmpArr(readPosTmpArr, readPosTmpArrBuf, readsNumSingleContig)==FAILED)
	{
		printf("line=%d, In %s(), cannot sort the reads in readPosTmp Array, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Fill data of the readPosTmp array to read list (RL)
	if(fillDataReadList(readListArr, readPosArr, &itemNumInReadListArr, &itemNumInReadPosArr, readPosTmpArr, readsNumSingleContig)==FAILED)
	{
		printf("line=%d, In %s(), cannot convert the readPosTmp array to read list (RL), error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################### Debug information ###############################
	//if(outputMatedReadsInReadListToFile("../matedReads", readListArr, readPosArr, itemNumInReadListArr)==FAILED)
	//{
	//	printf("line=%d, In %s(), cannot mated reads in read list (RL) to file, error!\n", __LINE__, __func__);
	//	return FAILED;
	//}
	// ############################### Debug information ###############################

	// get the valid insert size, and output their fragment size to binary file
	if(outputInsertSizeToFile(fpFragSize, readListArr, readPosArr, itemNumInReadListArr)==FAILED)
	{
		printf("line=%d, In %s(), cannot output valid insert size to file, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// free memory
	freeMemGetPESingleContig();

	return SUCCESSFUL;
}

/**
 * Initialize the memory for getting paired ends of single contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int initMemGetPESingleContig(estContig_t *estContig)
{
	// get the total reads number in single contig
	if(getTotalReadsNumOfSingleContig(&readsNumSingleContig, estContig->contighead)==FAILED)
	{
		printf("line=%d, In %s(), cannot get total reads number of contig %d, error!\n", __LINE__, __func__, estContig->contigID);
		return FAILED;
	}
	if(readsNumSingleContig<=0)
	{
		printf("line=%d, In %s(), contigID=%d, readsNumSingleContig=%ld, error!\n", __LINE__, __func__, estContig->contigID, readsNumSingleContig);
		return FAILED;
	}

	// allocate memory
	readPosTmpArr = (readPosTemp_t*) malloc(readsNumSingleContig * sizeof(readPosTemp_t));
	if(readPosTmpArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}
	readPosTmpArrBuf = (readPosTemp_t*) malloc(readsNumSingleContig * sizeof(readPosTemp_t));
	if(readPosTmpArrBuf==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	readListArr = (readList_t *) calloc(readsNumSingleContig, sizeof(readList_t));
	if(readListArr==NULL)
	{
		printf("line=%d, In %s(), cannot callocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	readPosArr = (readPos_t *) calloc(readsNumSingleContig, sizeof(readPos_t));
	if(readPosArr==NULL)
	{
		printf("line=%d, In %s(), cannot callocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Free memory of getting paired ends of single contig.
 */
void freeMemGetPESingleContig()
{
	readsNumSingleContig = 0;

	free(readPosTmpArr);
	readPosTmpArr = NULL;
	free(readPosTmpArrBuf);
	readPosTmpArrBuf = NULL;

	free(readListArr);
	readListArr = NULL;
	free(readPosArr);
	readPosArr = NULL;
}

/**
 * Get total reads number of a contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int getTotalReadsNumOfSingleContig(int64_t *totalReadsNum, contigtype *contighead)
{
	contigtype *contig;

	*totalReadsNum = 0;
	contig = contighead;
	while(contig)
	{
		if(contig->ridposnum>0)
		{
			*totalReadsNum += contig->ridposnum;
		}
		contig = contig->next;
	}

	return SUCCESSFUL;
}

/**
 * Fill the data of readPosTmp Array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int fillDataReadPosTmpArr(readPosTemp_t *readPosTmpArray, contigtype *contighead)
{
	int64_t tmpReadsNum;
	int i, ridposnum;
	contigtype *contig;
	successRead_t *pReadsArray;

	tmpReadsNum = 0;
	contig = contighead;
	while(contig)
	{
		if(contig->ridposnum>0)
		{
			ridposnum = contig->ridposnum;
			pReadsArray = contig->pridposorientation;
			for(i=0; i<ridposnum; i++)
			{
				readPosTmpArray[tmpReadsNum].readID = pReadsArray[i].rid;
				readPosTmpArray[tmpReadsNum].contigPos = contig->index;
				readPosTmpArray[tmpReadsNum].orientation = pReadsArray[i].orientation;
				readPosTmpArray[tmpReadsNum].matchBaseNum = pReadsArray[i].matchnum;
				tmpReadsNum ++;

				// ################### Debug information ##################
				if(tmpReadsNum>readsNumSingleContig)
				{
					printf("line=%d, In %s(), tmpReadsNum=%ld > readsNumSingleContig=%ld, error!\n", __LINE__, __func__, tmpReadsNum, readsNumSingleContig);
					return FAILED;
				}
				// ################### Debug information ##################
			}
		}
		contig = contig->next;
	}

	// ################### Debug information ##################
	if(tmpReadsNum!=readsNumSingleContig)
	{
		printf("line=%d, In %s(), tmpReadsNum=%ld != readsNumSingleContig=%ld, error!\n", __LINE__, __func__, tmpReadsNum, readsNumSingleContig);
		return FAILED;
	}
	// ################### Debug information ##################

	return SUCCESSFUL;
}

/**
 * Sort the reads in readPosTmp array to ascending order according to their readID.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int radixSortReadPosTmpArr(readPosTemp_t *readPosTmpArray, readPosTemp_t *readPosTmpArrayBuf, int64_t readsNumInArray)
{
	struct partNode
	{
		uint32_t curItemNum;
		uint32_t totalItemNum;
		uint64_t firstRow;
	};

	int64_t i, step, total;
	readPosTemp_t *data, *buf;
	struct partNode *part;
	int partArrSize, stepBits, maxStepLen;
	unsigned int bitMask;
	unsigned int hashcode, firstRow, curItemNum;

	stepBits = 16;
	maxStepLen = 64;
	partArrSize = 1 << stepBits;
	bitMask = (1 << stepBits) - 1;

	part = (struct partNode *) malloc(partArrSize * sizeof(struct partNode));
	if(part==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// Begin to sort
	step = 0;
	while(step!=maxStepLen)
	{
		// set the data and buf
		if((step/stepBits)%2==1)
		{
			buf = readPosTmpArray;
			data = readPosTmpArrayBuf;
		}else
		{
			data = readPosTmpArray;
			buf = readPosTmpArrayBuf;
		}

		// count the number
		if(memset(part, 0, partArrSize * sizeof(struct partNode))==NULL)
		{
			printf("line=%d, In %s(), cannot reset memory, error!\n", __LINE__, __func__);
			free(part);
			return FAILED;
		}
		for(i=0; i<readsNumInArray; i++)
			part[ (data[i].readID >> step) & bitMask ].totalItemNum ++;

		// initialize the part array
		total = 0;
		for(i=0; i<partArrSize; i++)
		{
			part[i].firstRow = total;
			total += part[i].totalItemNum;
		}

		// copy the data to the right place
		for(i=0; i<readsNumInArray; i++)
		{
			hashcode = (data[i].readID >> step) & bitMask;
			firstRow = part[hashcode].firstRow;
			curItemNum = part[hashcode].curItemNum;
			if(memcpy(buf+firstRow+curItemNum, data+i, sizeof(readPosTemp_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				free(part);
				return FAILED;
			}
			part[hashcode].curItemNum ++;
		}

		step += stepBits;

		//######################## Debug information #######################
		for(i=0; i<partArrSize; i++)
		{
			if(part[i].curItemNum!=part[i].totalItemNum)
			{
				printf("line=%d, In %s(), in part[%ld], curItemNum=%u != totalItemNum=%u, error!\n", __LINE__, __func__, i, part[i].curItemNum, part[i].totalItemNum);
				free(part);
				return FAILED;
			}
		}
		//######################## Debug information #######################
	}

	free(part);

	return SUCCESSFUL;
}

/**
 * Fill data of the readPosTmp array to read list (RL).
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int fillDataReadList(readList_t *readListArray, readPos_t *readPosArray, int64_t *itemNumInReadListArray, int64_t *itemNumInReadPosArray, readPosTemp_t *readPosTmpArray, int64_t itemNumInReadPosTmpArray)
{
	int64_t i, j;
	int sameItemNum;

	*itemNumInReadListArray = *itemNumInReadPosArray = 0;

	i = 0;
	while(i<itemNumInReadPosTmpArray)
	{
		sameItemNum = 1;
		for(j=i+1; j<itemNumInReadPosTmpArray; j++)
		{
			if(readPosTmpArray[i].readID!=readPosTmpArray[j].readID)
				break;
			sameItemNum ++;
		}

		// fill data of the item
		readListArray[*itemNumInReadListArray].readID = readPosTmpArray[i].readID;
		readListArray[*itemNumInReadListArray].matchNum = sameItemNum;
		readListArray[*itemNumInReadListArray].firstRow = *itemNumInReadPosArray;
		(*itemNumInReadListArray) ++;

		for(j=0; j<sameItemNum; j++)
		{
			readPosArray[(*itemNumInReadPosArray)+j].contigPos = readPosTmpArray[i+j].contigPos;
			readPosArray[(*itemNumInReadPosArray)+j].matchBaseNum = readPosTmpArray[i+j].matchBaseNum;
			readPosArray[(*itemNumInReadPosArray)+j].orientation = readPosTmpArray[i+j].orientation;
		}
		(*itemNumInReadPosArray) += sameItemNum;

		i += sameItemNum;
	}

	// ###################### Debug information ########################
	if((*itemNumInReadPosArray)!=itemNumInReadPosTmpArray)
	{
		printf("line=%d, In %s(), itemNumInReadPosArray=%ld != itemNumInReadPosTmpArray=%ld, error!\n", __LINE__, __func__, *itemNumInReadPosArray, itemNumInReadPosTmpArray);
		return FAILED;
	}
	// ###################### Debug information ########################

	// ###################### Debug information ########################
	//printf("itemNumInReadListArray=%ld, itemNumInReadPosArray=%ld\n", *itemNumInReadListArray, *itemNumInReadPosArray);
	// ###################### Debug information ########################

	return SUCCESSFUL;
}

/**
 * Get the valid insert size, and output their fragment size to binary file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int outputInsertSizeToFile(FILE *fpFragSize, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray)
{
	int64_t i;
	int32_t fragmentSize;
	readList_t *pLeftRead, *pRightRead;
	readPos_t *pLeftReadPos, *pRightReadPos;

	i = 0;
	while(i<itemNumInReadListArray-1)
	{
		if(readListArray[i].readID % 2 == 1 && readListArray[i].readID+1==readListArray[i+1].readID)
		{ // the odd readID
			//if(readListArray[i].matchNum==1 && readListArray[i+1].matchNum==1)
			if(readListArray[i].matchNum==1 && readListArray[i+1].matchNum==1
				&& readPosArray[readListArray[i].firstRow].matchBaseNum==readLen && readPosArray[readListArray[i+1].firstRow].matchBaseNum==readLen)
			{
				if(readPosArray[ readListArray[i].firstRow ].contigPos < readPosArray[ readListArray[i+1].firstRow ].contigPos)
				{
					pLeftRead = readListArray + i;
					pRightRead = readListArray + i + 1;
					pLeftReadPos = readPosArray + readListArray[i].firstRow;
					pRightReadPos = readPosArray + readListArray[i+1].firstRow;
				}else
				{
					pLeftRead = readListArray + i + 1;
					pRightRead = readListArray + i;
					pLeftReadPos = readPosArray + readListArray[i+1].firstRow;
					pRightReadPos = readPosArray + readListArray[i].firstRow;
				}

				if(pLeftReadPos->orientation==ORIENTATION_PLUS && pRightReadPos->orientation==ORIENTATION_MINUS)
				{
					fragmentSize = pRightReadPos->contigPos + pRightReadPos->matchBaseNum - 1 - pLeftReadPos->contigPos + 1;
					if(PEGivenType>=INSERT_PE_GIVEN_TYPE && fragmentSize<=MAX_INSERT_SIZE_FACTOR*meanSizeInsert)
					{
						if(fwrite(&fragmentSize, sizeof(int32_t), 1, fpFragSize)!=1)
						{
							printf("line=%d, In %s(), (%lu,%u,%u,%c), (%lu,%u,%u,%c), fwrite error!\n", __LINE__, __func__, pLeftRead->readID, pLeftReadPos->contigPos, pLeftReadPos->matchBaseNum, pLeftReadPos->orientation, pRightRead->readID, pRightReadPos->contigPos, pRightReadPos->matchBaseNum, pRightReadPos->orientation);
							return FAILED;
						}
					}else
					{
						if(fwrite(&fragmentSize, sizeof(int32_t), 1, fpFragSize)!=1)
						{
							printf("line=%d, In %s(), (%lu,%u,%u,%c), (%lu,%u,%u,%c), fwrite error!\n", __LINE__, __func__, pLeftRead->readID, pLeftReadPos->contigPos, pLeftReadPos->matchBaseNum, pLeftReadPos->orientation, pRightRead->readID, pRightReadPos->contigPos, pRightReadPos->matchBaseNum, pRightReadPos->orientation);
							return FAILED;
						}
					}
				}else
				{
					//printf("line=%d, In %s(), (%lu,%u,%u,%c), (%lu,%u,%u,%c), invalid paired end!\n", __LINE__, __func__, pLeftRead->readID, pLeftReadPos->contigPos, pLeftReadPos->matchBaseNum, pLeftReadPos->orientation, pRightRead->readID, pRightReadPos->contigPos, pRightReadPos->matchBaseNum, pRightReadPos->orientation);
				}
			}

			i += 2;
		}else
		{
			i ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * Reload the file and estimate the insert size and standard deviation.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int computeInsertSizeAndSdev(double *meanSizeInsert, double *standardDev, const char *fragmentSizeFile)
{
	int64_t fragNum, fragNum2;
	FILE *fpFragSize;
	int32_t fragmentSize;
	double tmpInsertSize, tmpSdev;

	fpFragSize = fopen(fragmentSizeFile, "rb");
	if(fpFragSize==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, fragmentSizeFile);
		return FAILED;
	}

	// estimate the insert size
	fragNum = 0;
	tmpInsertSize = 0;
	while(1)
	{
		if(fread(&fragmentSize, sizeof(int32_t), 1, fpFragSize)!=1)
		{
			if(feof(fpFragSize))
				break;
			else
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		tmpInsertSize += fragmentSize;
		fragNum ++;
	}
	*meanSizeInsert = tmpInsertSize / fragNum;

	if(fragNum<=1)
	{
		printf("line=%d, In %s(), fragNum=%ld, error!\n", __LINE__, __func__, fragNum);
		return FAILED;
	}

	// estimate the standard deviation
	fseek(fpFragSize, 0, SEEK_SET);   // seek to the beginning of the file
	tmpSdev = 0;
	fragNum2 = 0;
	while(1)
	{
		if(fread(&fragmentSize, sizeof(int32_t), 1, fpFragSize)!=1)
		{
			if(feof(fpFragSize))
				break;
			else
			{
				printf("line=%d, In %s(), fread error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		tmpSdev += (fragmentSize - *meanSizeInsert) * (fragmentSize - *meanSizeInsert);
		fragNum2 ++;
	}
	*standardDev = sqrt(tmpSdev / (fragNum-1));

	fclose(fpFragSize);
	fpFragSize = NULL;

	// ########################## Debug information ###########################
	if(fragNum!=fragNum2)
	{
		printf("line=%d, In %s(), fragNum=%ld != fragNum2=%ld, error!\n", __LINE__, __func__, fragNum, fragNum2);
		return FAILED;
	}
	// ########################## Debug information ###########################

	printf("The estimated insert size  : %.2f\n", *meanSizeInsert);
	printf("The estimated standard dev : %.2f\n", *standardDev);
	printf("Pairs of used paired ends  : %ld\n", fragNum);

	return SUCCESSFUL;
}
