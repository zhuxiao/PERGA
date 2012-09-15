/*
 * PEAssembly.c
 *
 *  Created on: Dec 8, 2011
 *      Author: xiao
 */


#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildEstContigs(char *contigFile)
{
	int i, turnContigIndex;
	int64_t localContigID;
	double averOccNumNaviOccQueue;

	contigNumEstContigArr = 0;
	basesNum = 0; //拼接的所有contigs的碱基的总和

	localContigID = 0;
	kmerIndex = 0;
	contigsNum = 0;
	while(kmerIndex < hashTableSize)
	{
		contighead = NULL;
		contigtail = NULL;
		contig36 = NULL;
		contigIndex = 1;
		lastseq36[0] = '\0';
		itemNumDecisionTable = 0;
		successContig = NULL;
		assemblyRound = FIRST_ROUND_ASSEMBLY;  //第1轮拼接
		lockedReadsNum = 0; //初始化清空锁定的reads数量
		this_successReadNum = 0; //本次拼接中的成功reads数量置为0
		readsNumInPEHashArr = 0;
		regLenPEHash = 0;
		turnContigIndex = 0;
		allowedUpdatePEHashArrFlag = YES;
		localContigID ++;

		//取得拼接的首个kmers及其正向的碱基序列
		if(getFirstKmers(&kmerIndex, &firstKmer)==FAILED)
		{
			printf("line=%d, In %s(), cannot get first kmers, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(kmerIndex>=hashTableSize)
		{ //已经到达哈希表的尾部时, 则拼接结束
			break;
		}

		// initialize the contig chain
		if(initContig(&contighead, &contigtail)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the contig nodes, error!\n", __LINE__, __func__);
			return FAILED;
		}
		contigIndex = kmerSize;

		// initialize the decision table
		if(addFirstKmerToDecisionTable(kmers)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}


		//将kmer碱基添加进碱基数组
		strcpy(lastseq36, getKmerBaseByInt(kmerSeqIntAssembly));

		if(setEmptyNaviOccQueue(naviOccQueue, &itemNumNaviOccQueue, &frontRowNaviOccQueue, &rearRowNaviOccQueue)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the empty navigation occurrence queue, error!\n", __LINE__, __func__);
			return FAILED;
		}

		while(kmers[0]||kmers[1])
		{
			// ############################ Debug information ##############################
			//if(contigsNum+1==162 && contigIndex>=3597 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			//if(localContigID==3 && contigIndex>=73040 && assemblyRound!=FIRST_ROUND_ASSEMBLY)
			//{
			//	printf("localContigID=%ld, contigID=%d, contigIndex=%d\n", localContigID, contigsNum+1, contigIndex);
			//}
			// ############################ Debug information ##############################

			// initialize or update the PE hash table
			if(PEGivenType>NONE_PE_GIVEN_TYPE && contigIndex>=minContigLenUsingPE)
			{
				if(updatePEHashTable(contigIndex, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot update the PE hash table, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}

				//取正反向kmer
				if(readsNumInPEHashArr>=minReadsNumPEHashThres && regLenPEHash>=minRegLenUsingPE)
				{
					if(getNextKmerByMix(contigIndex, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannot get the next kmer by mix, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}
				}else
				{
					navigationFlag = NAVI_SE_FLAG;
					if(getNextKmerBySE(contigIndex)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}
				}
			}else
			{
				//取正反向kmer
				navigationFlag = NAVI_SE_FLAG;
				if(getNextKmerBySE(contigIndex)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
			}


			// check the reads number in the sub read region
			if(kmers[0] || kmers[1])
			{
				if(contigIndex>=minContigLenCheckingReadsNum)
				{
					if(updateReadsNumReg(itemNumSuccessReadsArr, contigIndex, assemblyRound)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannnot check the reads number in reads number region, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}
				}
			}


			if(kmers[0]==NULL && kmers[1]==NULL)
			{
				if(successContig==NULL)
				{ //没有成功的reads, 该contig拼接失败
					//printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, the successContig==NULL!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
					break;
				}

//				printf("localContigID=%ld, assemblyRound=%d, contigIndex=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, contigIndex, itemNumDecisionTable);
//				printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
//				printf("\tdistance=%d, readsNumRatio=%.2f\n", contigIndex-successContig->index, readsNumRatio);

				//开始下一轮的拼接的预处理
				if(assemblyRound==FIRST_ROUND_ASSEMBLY)
				{ //该contig的第一轮拼接结束, 将进行第二轮拼接

					assemblyRound ++;
					turnContigIndex = contigIndex;

					int returnCode = initSecondAssembly();
					if(returnCode==FAILED)
					{
						break;
					}else if(returnCode==ERROR)
					{
						return FAILED;
					}

					continue;  //开始执行第二轮的拼接
					//break; //只进行第一轮单向拼接

				}else
				{ //该contig的第2轮拼接结束, 该contig的拼接结束

					// ############################ Debug information ##############################
					//if(successContig==NULL && contigIndex>=CONTIG_LEN_THRESHOLD)
					//{ //出错
					//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, successContig==NULL, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
					//	return FAILED;
					//}
					// ############################ Debug information ##############################

					break;
				}
			}

			if(navigationFlag==NAVI_PE_FLAG)
			{
				if(updateNaviOccQueue(naviOccQueue, &itemNumNaviOccQueue, &frontRowNaviOccQueue, &rearRowNaviOccQueue, maxOccPE)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot update the navigation occurrence queue, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
			}else
			{
				if(updateNaviOccQueue(naviOccQueue, &itemNumNaviOccQueue, &frontRowNaviOccQueue, &rearRowNaviOccQueue, maxOccSE)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot update the navigation occurrence queue, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
			}

			contigIndex ++;


			// Append a base to contig tail
			if(addContigBase(&contigtail, kmerSeqIntAssembly[entriesPerKmer-1] & 3, contigIndex)==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, cannot add a contig base, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
				return FAILED;
			}

			// update the decision table according to kmers
			if(updateDecisionTable(kmers)==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, cannot update decision table, error!\n", __LINE__, __func__, localContigID);
				return FAILED;
			}

			// Update the reads status in decision table
			if(updateAssemblingreadsStatus()==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, cannot update reads status in decision table, error!\n", __LINE__, __func__, localContigID);
				return FAILED;
			}

			// Update the locked reads and their total number
			if(updateLockedReads()==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, cannot update locked reads, error!\n", __LINE__, __func__, localContigID);
				return FAILED;
			}

			//if(localContigID==2)
			//	outputReadsInDecisionTableToFile(decisionTable, itemNumDecisionTable, (int)localContigID, contigIndex);
			//outputFailedReadsInDecisionTable(decisionTable, itemNumDecisionTable, (int)localContigID, contigIndex);

			// update the finished reads in decision table, and record the successful reads into successful reads array
			if(updateFinishedReadsInDecisionTable()==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, cannot update finished reads in decision table, error!\n", __LINE__, __func__, localContigID);
				return FAILED;
			}

			if(contigIndex<=readLen)
			{ //如果contig的节点数目<=36，则contig36指向contighead节点，并将当次拼接kmer的序列的最后一个碱基添加进lastseq36数组
				contig36 = contighead;
				switch(kmerSeqIntAssembly[entriesPerKmer-1] & 3)
				{
					case 0: lastseq36[contigIndex-1] = 'A'; break;
					case 1: lastseq36[contigIndex-1] = 'C'; break;
					case 2: lastseq36[contigIndex-1] = 'G'; break;
					case 3: lastseq36[contigIndex-1] = 'T'; break;
				}
				lastseq36[contigIndex] = '\0';
			}else
			{  //如果contig中节点数目>36，则contig36后移一个位置，并将lastseq36数组的碱基前移一位，并将当次拼接的kmer的序列添加进lastseq36数组
				contig36 = contig36->next;
				for(i=1; i<readLen; i++) lastseq36[i-1] = lastseq36[i];
				switch(kmerSeqIntAssembly[entriesPerKmer-1] & 3)
				{
					case 0: lastseq36[readLen-1] = 'A'; break;
					case 1: lastseq36[readLen-1] = 'C'; break;
					case 2: lastseq36[readLen-1] = 'G'; break;
					case 3: lastseq36[readLen-1] = 'T'; break;
				}
				lastseq36[readLen] = '\0';
			}

			contigtype *tmp_successContig = successContig;
			if(itemNumSuccessReadsArr>0)//如果有成功结束的reads
			{

				// delete reads from De Bruijn graph
				if(delReadsFromGraph(successReadsArr, itemNumSuccessReadsArr, lastseq36)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, contigIndex=%d, cannot delete the reads from graph, error!\n", __LINE__, __func__, localContigID, contigsNum+1, contigIndex);
					outputSuccessReads(successReadsArr, itemNumSuccessReadsArr);
					return FAILED;
				}

				// add the successful reads information to contig chain
				if(addRidposToContig(successReadsArr, &itemNumSuccessReadsArr, contig36, contigIndex)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, cannot add a contig base, error!\n", __LINE__, __func__, localContigID, contigsNum+1);
					return FAILED;
				}

				tmp_successContig = getSuccessContig(contig36, successContig, contigIndex);

				if(tmp_successContig==NULL)
				{
					//printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, the tmp_successContig==NULL!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
					//outputContig(contig36);
					//return FAILED;
				}else
				{
					successContig = tmp_successContig;

					// Update the successful reads number
					this_successReadNum += itemNumSuccessReadsArr;
				}
			}

			//处理死循环
			if(successContig==NULL)
			{ //还未有成功的reads, 并且拼接长度已经超过60, 则该contig拼接失败, 退出
				if(contigIndex>2*readLen-kmerSize)
				{
					//printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, successContig==NULL!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
					break;
				}
			}else if(tmp_successContig==NULL)
			{
				break;
			//}else if((contigIndex-successContig->index > readLen-MIN_OVERLAP_LEN) || (contigIndex-successContig->index > 8 && ((readsNumRatio>2 || readsNumRatio<0.3) || secondOccSE>=minKmerOccSE)))
			}else if(contigIndex-successContig->index > readLen-MIN_OVERLAP_LEN )
			{ //已经有成功的reads, 则根据拼接的情况, 确定是否需要继续拼接

				number_of_overlap_less_than_threshold ++;

//				printf("===localContigID=%ld, assemblyRound=%d, contigIndex=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, contigIndex, itemNumDecisionTable);
//				printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
//				printf("\tdistance=%d, readsNumRatio=%.2f\n", contigIndex-successContig->index, readsNumRatio);

				//第二轮拼接时, contig长度大于100时才进行衔接操作
				if(assemblyRound==SECOND_ROUND_ASSEMBLY && contigIndex<CONTIG_LEN_THRESHOLD)
				{
					break;
				}

				//第二轮拼接的预处理
				if(assemblyRound==FIRST_ROUND_ASSEMBLY)
				{ //现在处于第1论拼接, 需要进行第二轮拼接

					assemblyRound ++; //第二轮拼接标记
					turnContigIndex = contigIndex;

					int returnCode = initSecondAssembly();
					if(returnCode==FAILED)
					{
						break;
					}else if(returnCode==ERROR)
					{
						return FAILED;
					}

				}else
				{ //现在已经处于第2论拼接, 则该contig链表的拼接结束

					// ############################ Debug information ##############################
					//if(successContig==NULL && contigIndex>=CONTIG_LEN_THRESHOLD)
					//{ //出错
					//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, successContig==NULL, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
					//	return FAILED;
					//}
					// ############################ Debug information ##############################

					break;
				}
			}
/*
			else if(contigIndex-successContig->index >= 7)
			{ //已经有成功的reads, 则根据拼接的情况, 确定是否需要继续拼接

				if(calcAverOccNaviOccQueue(&averOccNumNaviOccQueue, naviOccQueue, itemNumNaviOccQueue)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, cannot compute the average occurrence in navigation occurrence queue, error!\n", __LINE__, __func__, localContigID, contigsNum+1);
					return FAILED;
				}

				if((averOccNumNaviOccQueue<2.5*minKmerOccSE) || (readsNumRatio>2 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
				{
					number_of_overlap_less_than_threshold ++;

					printf("===localContigID=%ld, assemblyRound=%d, contigIndex=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, contigIndex, itemNumDecisionTable);
					printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					printf("\tdistance=%d, readsNumRatio=%.2f, averOccNumNaviOccQueue=%.2f\n", contigIndex-successContig->index, readsNumRatio, averOccNumNaviOccQueue);

					//第二轮拼接时, contig长度大于100时才进行衔接操作
					if(assemblyRound==SECOND_ROUND_ASSEMBLY && contigIndex<CONTIG_LEN_THRESHOLD)
					{
						break;
					}

					//第二轮拼接的预处理
					if(assemblyRound==FIRST_ROUND_ASSEMBLY)
					{ //现在处于第1论拼接, 需要进行第二轮拼接

						assemblyRound ++; //第二轮拼接标记
						turnContigIndex = contigIndex;

						int returnCode = initSecondAssembly();
						if(returnCode==FAILED)
						{
							break;
						}else if(returnCode==ERROR)
						{
							return FAILED;
						}

					}else
					{ //现在已经处于第2论拼接, 则该contig链表的拼接结束

						// ############################ Debug information ##############################
						//if(successContig==NULL && contigIndex>=CONTIG_LEN_THRESHOLD)
						//{ //出错
						//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, successContig==NULL, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
						//	return FAILED;
						//}
						// ############################ Debug information ##############################

						break;
					}
				}
			}
*/
		}//end while(kmer)

		if(successContig)
		{
			// 将contig节点退回到最近成功的contigIndex的位置, 并将之后的contig节点删掉.
			if(updateContigtailnodes(contighead, successContig, &contigIndex)==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, cannot update Contigtail nodes, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
				return FAILED;
			}

			//====================================================
			// trim a read length of contig nodes at tail
			if(contigIndex>=3*readLen)
			{
				if(trimContigTailByReadLen(contighead, &contigtail, &successContig, &contigIndex, SECOND_ROUND_ASSEMBLY)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot trim contig nodes at contig tail by a read length, error!\n", __LINE__, __func__, localContigID);
					return ERROR;
				}
			}
			//====================================================

			if(PEGivenType>NONE_PE_GIVEN_TYPE && readsNumInPEHashArr>0)
			{
				if(cleanReadsFromPEHashtable()==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot clean PE hash table, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
			}

			//只考虑长度大于100的contig链表
			if(contigIndex>=minContigLenEst)
			{ // if the contig length is larger than minContigLenEst, then save it to the contig array
				contigsNum ++;

				// ############################ Debug information ##############################
				//printf("==== localContigID=%ld, contigID=%d, contigLen=%d, turnContigIndex=%d.\n", localContigID, contigsNum, contigIndex, turnContigIndex);
				// ############################ Debug information ##############################

				//save the contig to the contig array
				estContigArr[contigNumEstContigArr].contigID = contigsNum;
				estContigArr[contigNumEstContigArr].contighead = contighead;
				estContigArr[contigNumEstContigArr].contigLen = contigIndex;

				successReadNum += this_successReadNum;
				basesNum += contigIndex;
				contigNumEstContigArr ++;

				if(contigNumEstContigArr>=MAX_NUM_EST_CONTIG || basesNum>=TOTAL_CONTIG_LEN_EST_THRES)
					break;
			}else
			{ //contig长度小于100, 则恢复该contig链表中的reads

				successReadNum -= this_successReadNum; //更新成功的reads数量
/*
				//长度小于100, 则将该contig中的reads的删除标记还原
				if(recoverDeledReads(contighead)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, contigIndex=%d, cannot recover the deleted reads, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, contigIndex);
					break;
				}
*/

				//释放该contig拼接过程中占用的内存，并初始化新的contig
				releaseContig(contighead);
			}
		}else
		{
			// ############################ Debug information ##############################
			//if(contigIndex>=CONTIG_LEN_THRESHOLD)
			//{
			//	printf("line=%d, In %s(), localContigID=%ld, contigID=%d, assemblyRound=%d, contigIndex=%d, successContig==NULL, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
			//	return FAILED;
			//}
			// ############################ Debug information ##############################

			//释放该contig拼接过程中占用的内存，并初始化新的contig
			releaseContig(contighead);
		}

		// clean the PE hash table
		if(PEGivenType>NONE_PE_GIVEN_TYPE && readsNumInPEHashArr>0)
		{
			if(cleanReadsFromPEHashtable()==FAILED)
			{
				printf("line=%d, In %s(), localContigID=%ld, cannot clean reads from PE hash table, error!\n", __LINE__, __func__, localContigID);
				return FAILED;
			}
		}

	} //end while(kmerIndex < TABLE_SIZE_DE_BRUIJN)

	// ############################ Debug information ##############################
	if(contigNumEstContigArr!=contigsNum)
	{
		printf("line=%d, In %s(), contigNumEstContigArr=%d != contigsNum=%d, error!\n", __LINE__, __func__, contigNumEstContigArr, contigsNum);
		return FAILED;
	}
	// ############################ Debug information ##############################

	// ############################ Debug information ##############################
	FILE *fpContigBaseTmp, *fpContigHangingTmp;
	char hangingFile[256];
	strcpy(hangingFile, contigFile);
	strcat(hangingFile, ".hang");
	fpContigBaseTmp = fopen(contigFile, "w");
	if(fpContigBaseTmp==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFile);
		return FAILED;
	}
	if(hangingContigOutFlag==YES)
	{
		fpContigHangingTmp = fopen(hangingFile, "w");
		if(fpContigHangingTmp==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, hangingFile);
			return FAILED;
		}
	}
	for(i=0; i<contigNumEstContigArr; i++)
	{
		if(outputContigToFile(fpContigBaseTmp, BASE_TYPE_FASTA_CONTIG_FILE, estContigArr[i].contighead, estContigArr[i].contigID, estContigArr[i].contigLen)==FAILED)
		{
			printf("line=%d, In %s(), cannot output contig nodes to file, error!\n", __LINE__, __func__);
			return FAILED;
		}
		if(hangingContigOutFlag==YES)
		{
			if(outputContigToFile(fpContigHangingTmp, HANGING_READ_TYPE_CONTIG_FILE, estContigArr[i].contighead, estContigArr[i].contigID, estContigArr[i].contigLen)==FAILED)
			{
				printf("line=%d, In %s(), cannot output contig nodes to file, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
	}
	fclose(fpContigBaseTmp);
	fpContigBaseTmp = NULL;
	if(hangingContigOutFlag==YES)
	{
		fclose(fpContigHangingTmp);
		fpContigHangingTmp = NULL;
	}
	// ############################ Debug information ##############################

	if(contigNumEstContigArr>0)
		return SUCCESSFUL;
	else
		return FAILED;
}

/**
 * Get the next kmer by mixture of PE and SE.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getNextKmerByMix(int contigNodesNum, int assemblyRound)
{
	unsigned int tmp_kmerseqInt, tmp_kmerseqPE, tmp_kmerseqSE;
	kmertype *tmp_kmers[2]/*, *tmp_kmersPE[2], *tmp_kmersSE[2]*/;
	int i;
	double sumSecondOccSE;

	if(itemNumDecisionTable>MAX_DECISION_TABLE_SIZE_HTRES)
	{
		kmers[0] = kmers[1] = NULL;
		return SUCCESSFUL;
	}

	navigationFlag = NAVI_PE_FLAG;
	tmp_kmers[0] = kmers[0];
	tmp_kmers[1] = kmers[1];
	if(memcpy(tmpKmerSeqIntAssembly, kmerSeqIntAssembly, entriesPerKmer*sizeof(uint64_t))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// get the score and occsNum of PE
	if(getNextKmerByPE(contigNodesNum, assemblyRound)==FAILED)
	{
		printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(kmers[0] || kmers[1])
	{
		if((readsNumRatio<0.2 || readsNumRatio>3) && (successContig!=NULL && contigIndex-successContig->index >= 15))
		{
			kmers[0] = kmers[1] = NULL;
		}else if((maxOccPE>0 && maxOccPE<=minKmerOccPE) && (readsNumRatio<0.3 || readsNumRatio>3))
		{
			kmers[0] = kmers[1] = NULL;
		}else if(readsNumRatio>2.5 && (successContig!=NULL && contigIndex-successContig->index >= 20))
		{
			kmers[0] = kmers[1] = NULL;
		}else if(readsNumRatio>6 && (maxOccPE>maxOccNumFaiedPE && secondOccPE/maxOccPE>0.2))
		{
			kmers[0] = kmers[1] = NULL;
		}else if(readsNumRatio>5)
		{
			kmers[0] = kmers[1] = NULL;
		}else if(maxOccPE==1 && readsNumRatio < 0.7)
		{
			kmers[0] = kmers[1] = NULL;
		}
//		else if(readsNumRatio < 0.6 && (successContig!=NULL &&  contigIndex-successContig->index >= 13))
//		{
//			kmers[0] = kmers[1] = NULL;
//		}
	}


	if(kmers[0]==NULL && kmers[1]==NULL)
	{
//		if(maxOccPE<2)
//		{
			// get the score and occsNum of SE
			navigationFlag = NAVI_SE_FLAG;
			kmers[0] = tmp_kmers[0];
			kmers[1] = tmp_kmers[1];
			if(memcpy(kmerSeqIntAssembly, tmpKmerSeqIntAssembly, entriesPerKmer*sizeof(uint64_t))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(getNextKmerBySE(contigNodesNum)==FAILED)
			{
				printf("line=%d, In %s(), cannot get next kmer, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(kmers[0] || kmers[1])
			{
				sumSecondOccSE = 0;
				for(i=0; i<4; i++) if(i!=maxOccIndexSE) sumSecondOccSE += occsNumSE[i];

				// check the scores and occsNums of PE and SE
				//if(maxOccIndexPE==-1 && maxOccSE>OCCS_NUM_SE_FAILED_PE_FACTOR*averKmerOcc)  //==================================
				//if(maxOccIndexPE==-1 && maxOccSE>maxFirstOcc)
				//if(maxOccIndexPE==-1 && maxOccSE>maxOccNumFaiedPE)
				//if((maxOccIndexPE==-1 && (secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO && maxOccSE>maxOccNumFaiedPE)) && navigationID==0)
				//if((maxOccIndexPE==-1 && ((secondOccSE>=minKmerOccSE && maxOccSE>maxOccNumFaiedPE) || maxOccSE<minKmerOccSE)) && navigationID==0)
				//if((maxOccIndexPE==-1 && ((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)) || maxOccSE<minKmerOccSE)) && navigationID==0)
				//if((maxOccIndexPE==-1 && ((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>0.5)) || maxOccSE<minKmerOccSE)) && navigationID==0)
				//if(((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>0.5)) || maxOccSE<minKmerOccSE) && navigationID==0)
				if((secondOccSE>=minKmerOccSE && (maxOccSE>maxOccNumFaiedPE || secondOccSE/maxOccSE>0.5)) || maxOccSE<minKmerOccSE)
				{
					if(maxOccIndexPE==-1)
					{
						if(secondOccSE/maxOccSE>0.5)
							kmers[0] = kmers[1] = NULL;
						else if(maxOccSE==1 && (successContig!=NULL &&  contigIndex-successContig->index >= 3))
						{
							if(maxOccIndexPE!=maxOccIndexSE)
								kmers[0] = kmers[1] = NULL;
						}
					}else
					{
						if(maxOccPE>=minKmerOccPE && maxOccSE>=minKmerOccSE && maxOccIndexPE==maxOccIndexSE && secondOccSE/maxOccSE<=0.5)
						{
							if(readsNumRatio>6 && secondOccSE/maxOccSE>0.15)
								kmers[0] = kmers[1] = NULL;
						}else
						{
							if(secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>0.1)
								kmers[0] = kmers[1] = NULL;
						}
					}

				//}else if(secondOccPE>minKmerOccPE && secondOccSE>=2*minKmerOccSE && secondOccSE/maxOccSE>0.5)
				}
//				else if(secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>0.15)
//				{
//					kmers[0] = kmers[1] = NULL;
//				}
				//else if(maxOccSE>maxOccNumFaiedPE && (readsNumRatio<0.2 || readsNumRatio>2.5 || (successContig!=NULL && contigIndex-successContig->index >= 15)))
				//else if(maxOccSE>maxOccNumFaiedPE && (readsNumRatio<0.2 || readsNumRatio>2.5 || (successContig!=NULL && contigIndex-successContig->index >= 20)))
				else if(maxOccSE>maxOccNumFaiedPE && (readsNumRatio<0.2 || readsNumRatio>3.5 || (successContig!=NULL && contigIndex-successContig->index >= 15)))
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}else if(readsNumRatio>2 && (successContig!=NULL && contigIndex-successContig->index >= 20))
				//}else if(readsNumRatio>2 && (successContig!=NULL && contigIndex-successContig->index >= 30))
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}
				else if(secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>0.2 && (successContig!=NULL &&  contigIndex-successContig->index >= 7))
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}
				//===============
				else if(maxOccSE>maxOccNumFaiedPE && readsNumRatio<0.25 && (successContig!=NULL &&  contigIndex-successContig->index >= 7))
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}

				else if(secondOccSE>=2*minKmerOccSE && secondOccSE/maxOccSE>0.35 && readsNumRatio < 0.8)
				//else if(secondOccSE>=2*minKmerOccSE && secondOccSE/maxOccSE>0.35)
				//else if(secondOccSE>=2*minKmerOccSE && secondOccSE/maxOccSE>0.1 && readsNumRatio < 1)
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}
				else if(secondOccSE>=3*minKmerOccSE && sumSecondOccSE/maxOccSE>0.2 && (readsNumRatio >1 || readsNumRatio < 0.7))
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}else if(secondOccSE<minKmerOccSE && sumSecondOccSE/maxOccSE>0.4 && maxOccSE-sumSecondOccSE<2)
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}

				else if(secondOccSE>minKmerOccSE && sumSecondOccSE/maxOccSE>0.2 && readsNumRatio >1.6)
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}
				else if(maxOccSE<10*minKmerOccSE && secondOccSE>minKmerOccSE && sumSecondOccSE/maxOccSE>0.25)
				{
					if(maxOccIndexPE!=maxOccIndexSE)
						kmers[0] = kmers[1] = NULL;
				}

//				else if(maxOccIndexPE==-1 && maxOccSE > 19*minKmerOccSE && readsNumRatio < 0.615 && (successContig!=NULL &&  contigIndex-successContig->index >= 5))
//				{
//					//if(maxOccIndexPE!=maxOccIndexSE)
//						kmers[0] = kmers[1] = NULL;
//				}
//
//				else if(maxOccIndexPE==-1 && maxOccSE > 12*minKmerOccSE && readsNumRatio < 0.87 && (successContig!=NULL &&  contigIndex-successContig->index >= 10))
//				{
//					kmers[0] = kmers[1] = NULL;
//				}

//				else if(maxOccSE > 8*minKmerOccSE && readsNumRatio < 0.9 && (successContig!=NULL &&  contigIndex-successContig->index >= 13))
//				{
//					if(maxOccIndexPE!=maxOccIndexSE)
//						kmers[0] = kmers[1] = NULL;
//				}

//				else if(readsNumRatio < 0.6 && (successContig!=NULL &&  contigIndex-successContig->index >= 13))
//				{
//					kmers[0] = kmers[1] = NULL;
//				}
//				else if(readsNumRatio>6 && (maxOccSE>maxOccNumFaiedPE && secondOccSE/maxOccSE>0.2))
//					kmers[0] = kmers[1] = NULL;

				navigationID = 0;
				navigationNumSE ++;
			}
//		}
	}else
	{
		navigationID = 1;
		navigationNumSE = 0;
	}

	return SUCCESSFUL;
}

/**
*    取下次拼接的kmer.
*    策略, 两个过程:
*        (1) 首先, 以决策表为基础, 计算得分. 该过程只考虑决策表中存在的reads.
*        (2) 其次, 计算决策表中不存在的reads的得分. 也就是新的reads的得分.
*/
short getNextKmerByPE(int contigNodesNum, int assemblyRound)
{
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	kmertype *tmp_kmers[4][2] = {{0,0},{0,0},{0,0},{0,0}};  //tmp_kmers[i][0]为正向的kmer, tmp_kmers[i][1]为反向的kmer
	short validKmerNum = 0, base_index = 0; //有效的kmer数目
	int i = 0, j;
	//double maxOcc = 0, secondOcc = 0;
	//int maxOccIndex = -1, secondOccIndex = -1;
	kmer_len = 0;

	if(itemNumDecisionTable>MAX_DECISION_TABLE_SIZE_HTRES)
	{
		kmers[0] = kmers[1] = NULL;
		return SUCCESSFUL;
	}

	//将8个正反向kmer添加进临时数组tmp_kmers
	for(i=0; i<4; i++)
	{
		occsNumPE[i] = 0;

		//开始计算
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | i) & lastEntryMask;

		tmp_kmers[i][0] = getKmer(tmp_kmerseq, deBruijnGraph); //取得kmer的指针放入数组
		tmp_kmers[i][1] = getReverseKmer(tmp_kmerseqRev, tmp_kmerseq, deBruijnGraph); //取得反向互补的kmer的指针放入数组
		if(tmp_kmers[i][0] || tmp_kmers[i][1])
		{  //两个kmer中只要有一个存在, 就统计其数量
			validKmerNum ++;
			base_index = i;
		}

	}

	//检测kmer的数量
	if(validKmerNum==1)
	{ //只有一个kmer, 则将该kmer返回
		//计算该kmer得分,
		//score[base_index] = computeKmerScoreByPE(tmp_kmers[base_index], occsNum+base_index, assemblingreads, numassemblingreads);
		if(computeKmerOccNumUnlockedByPE(tmp_kmers[base_index], occsNumPE+base_index, contigNodesNum, assemblyRound)==FAILED)
		{
			printf("In %s(), cannot compute the occsNum of a kmer, error!\n", __func__);
			return FAILED;
		}

		//========================= condition 1 ==============================
		if(occsNumPE[base_index]>0)
		//if(occsNumPE[base_index]>=minKmerOccPE)
		{
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
				}
				kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | base_index) & lastEntryMask;

			kmers[0] = tmp_kmers[base_index][0];
			kmers[1] = tmp_kmers[base_index][1];

			maxOccIndexPE = base_index;
			maxOccPE = occsNumPE[maxOccIndexPE];
			secondOccIndexPE = -1;
			secondOccPE = 0;
		}else
		{
			kmers[0] = kmers[1] = NULL;
		}

		return SUCCESSFUL;
	}else if(validKmerNum==0)
	{ //没有可以扩展的kmer, 则返回NULL
		kmers[0] = kmers[1] = NULL;
		return SUCCESSFUL;
	}


	//****************************************************************************************************

	//开始计算每个kmer得分
	validKmerNum = 0;
	maxOccPE = 0, secondOccPE = 0;
	maxOccIndexPE = -1, secondOccIndexPE = -1;
	//tmp_max = 0, tmp_maxOcc = 0, tmp_maxIndex = -1, tmp_maxOccIndex = -1;
	for(i=0; i<4; i++)
	{
		occsNumPE[i] = 0;

		//kmer连接数目限制
		//######################### begin #############################//
		if(tmp_kmers[i][0]==NULL && tmp_kmers[i][1]==NULL)
		{  //如果该kmer不存在, 并且反向互补的kmer也不存在, 则检测下一个kmer
			continue;
		}

		//========================= condition 2 ==============================
		else
		{
			//修剪kmer连接数小于阈值的kmers
			if(tmp_kmers[i][0]!=NULL && tmp_kmers[i][1]!=NULL)
			{
				if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<minKmerOccPE)
				//if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<=0)
				{
					continue;
				}
			}else if(tmp_kmers[i][0]!=NULL)
			{
				if(tmp_kmers[i][0]->multiplicity<minKmerOccPE)
				//if(tmp_kmers[i][0]->multiplicity<=0)
				{
					continue;
				}
			}else if(tmp_kmers[i][1]!=NULL)
			{
				if(tmp_kmers[i][1]->multiplicity<minKmerOccPE)
				//if(tmp_kmers[i][1]->multiplicity<=0)
				{
					continue;
				}
			}
		}

		//########################## end ########################//


		//if(contigNodesNum>READ_LEN)
//			if(contigNodesNum>kmer_len)
//			//if(contigNodesNum>=kmer_len)  //--bad result
//			{
//				if(computeLongKmerScoreByPE(tmp_kmers[i], scorePE+i, occsNumPE+i, kmer_len, contigNodesNum, assemblyRound)==FAILED)
//				{
//					printf("line=%d, In %s(), cannot compute score and occsNum of long k-mers, error!\n", __LINE__, __func__);
//					return FAILED;
//				}
//			}else
//			{

			//if(computeKmerScoreByPE(tmp_kmers[i], score+i, occsNum+i, contigNodesNum, assemblyRound)==FAILED)		// 100_1
			if(computeKmerOccNumUnlockedByPE(tmp_kmers[i], occsNumPE+i, contigNodesNum, assemblyRound)==FAILED)	//100_2
			{
				printf("line=%d, In %s(), cannot compute occsNum of k-mers, error!\n", __LINE__, __func__);
				return FAILED;
			}
//			}

		//score[i] = computeKmerScore(tmp_kmers[i], occsNum+i, assemblingreads, numassemblingreads);

		//========================= condition 3 ==============================
		//if(score[i]>=MIN_SCORE_THRESHOLD/**MIN_SCORE_FACTOR*/)
		//if(score[i]>=MIN_SCORE_THRESHOLD && occsNum[i]>=MIN_CONNECT_KMER_NUM)
		//if(occsNum[i]>=MIN_CONNECT_KMER_NUM)
		if(occsNumPE[i]>0)
		{ //该kmer得分>0
			validKmerNum ++;  //有效的kmer数目+1

//				if(score[i]>tmp_max)
//				{
//					tmp_max = score[i];
//					tmp_maxIndex = i;
//				}
		}
//			else
//			{
//				score[i] = 0;
//				occsNum[i] = 0;
//			}
	}

	if(validKmerNum>0)
	{
		maxOccPE = 0, maxOccIndexPE = -1, secondOccPE = 0, secondOccIndexPE = -1;
		for(j=0; j<4; j++)
		{
			if(maxOccPE<occsNumPE[j])
			{
				secondOccPE = maxOccPE;
				secondOccIndexPE = maxOccIndexPE;
				maxOccPE = occsNumPE[j];
				maxOccIndexPE = j;
			}else if(secondOccPE<occsNumPE[j])
			{
				secondOccPE = occsNumPE[j];
				secondOccIndexPE = j;
			}
		}
	}

//	if(validKmerNum>1)
//	{
//		printf("validKmerNum=%d, contigNodexNum=%d\n", validKmerNum, contigNodesNum);
//	}

/*
	//========================= condition 4 ==============================
	//if(validKmerNum>1 && maxIndex1!=maxOccIndex)
	//if(validKmerNum>1 && occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
	//if(validKmerNum>0 && occsNumPE[maxIndex1]<MIN_CONNECT_KMER_NUM) //--best result
	if(validKmerNum>0 && occsNumPE[maxIndex1]<minKmerOccPE)
	//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || (validKmerNum>1 && maxIndex1!=maxOccIndex))
	//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || occsNum[secondOccIndex]>=MIN_CONNECT_KMER_NUM)
	{
		validKmerNum = 0;
	}


	//========================= condition 5 ==============================
	//if(validKmerNum>0 && (float)occsNum[maxIndex1]/numassemblingreads < VALID_OCC_RATIO)
	//if(validKmerNum>0 && occsNum[maxIndex1] > MAX_OCC_NUM)
	//if((float)occsNum[thiskmerseq & 3]/numassemblingreads < VALID_OCC_RATIO || occsNum[thiskmerseq & 3] > MAX_OCC_NUM)
	//if((float)occsNum[thiskmerseq & 3]/numassemblingreads < VALID_OCC_RATIO && occsNum[thiskmerseq & 3] > MAX_OCC_NUM)
	//if(validKmerNum>1 && (occsNum[secondIndex1] > 7 || occsNum[maxIndex1] > 80)) //--best result
	//if(validKmerNum>1 && (occsNum[secondIndex1] > MAX_SECOND_OCC_NUM || occsNum[maxIndex1] > MAX_FIRST_OCC_NUM))
	//if(validKmerNum>1 && kmer_len>=MIN_KMER_SIZE && (occsNum[secondIndex1] > 10 || occsNum[maxIndex1] > MAX_OCC_NUM))
	//if(validKmerNum>1 && (occsNum[secondIndex1] > MAX_SECOND_OCC_NUM || occsNum[maxIndex1] > MAX_FIRST_OCC_NUM
	//		|| ((float)occsNum[secondIndex1]/occsNum[maxIndex1]>0.5 && occsNum[secondOccIndex]>=MIN_CONNECT_KMER_NUM))) //-- good result
	if(validKmerNum>1 && (occsNumPE[secondIndex1] > maxSecondOcc || (occsNumPE[maxIndex1] > maxFirstOcc && occsNumPE[secondIndex1]>=maxSecondOcc)
			|| ((float)occsNumPE[secondIndex1]/occsNumPE[maxIndex1]>SECOND_FIRST_OCC_RATIO && occsNumPE[secondIndex1]>=maxSecondOcc))) //-- best result
	{
		validKmerNum = 0;
	}


	//========================= condition 6 ==============================
	if(validKmerNum>1 && ((float)occsNumPE[maxIndex1]/itemNumDecisionTable < VALID_OCC_RATIO && occsNumPE[secondIndex1] >= maxSecondOcc))
	{
		validKmerNum = 0;
	}


	//========================= condition 7 ==============================
	if(validKmerNum>1 && second/max > SECOND_FIRST_SECORE_RATIO)
	{
		validKmerNum = 0;
	}
*/


	//========================= condition 8 ==============================
	//=====these several lines have bad result, thus they are omitted. =======//
	//if(validKmerNum>1 && (secondOcc/maxOcc>SECOND_FIRST_OCC_RATIO || (kmer_len==25 && secondOcc>SECOND_OCC_THRESHOLD)))
	//if(validKmerNum>1 && kmer_len==25 && ((secondOcc/maxOcc>SECOND_FIRST_OCC_RATIO && secondOcc>MIN_CONNECT_KMER_NUM) || (secondOcc>SECOND_OCC_THRESHOLD)))
	//if(validKmerNum>1 && ((secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO/2 && secondOccPE>minKmerOccPE) || (secondOccPE>maxSecondOcc)))
	//if(maxOccPE<minKmerOccPE || (validKmerNum>1 && ((secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO/2 && secondOccPE>minKmerOccPE) || (secondOccPE>maxSecondOcc))))
	if(maxOccPE<minKmerOccPE || (validKmerNum>1 && (secondOccPE/maxOccPE>SECOND_FIRST_OCC_RATIO/2 || (secondOccPE>maxSecondOcc))))
	{
		if(secondOccPE/maxOccPE>SECOND_FIRST_OCC_FAILED_RATIO)
			validKmerNum = 0;
	}



	if(validKmerNum==1)
	{
		//if(occsNum[tmp_maxIndex]<MIN_CONNECT_KMER_NUM)
//			if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
//			{
//				kmers[0] = kmers[1] = NULL;
//			}else
//			{
			//更新kmerseq，并返回得分最大的kmer
			if(entriesPerKmer>=2)
			{
				for(j=0; j<entriesPerKmer-2; j++)
				{
					kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
				}
				kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
			}
			kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexPE) & lastEntryMask;

			kmers[0] = tmp_kmers[maxOccIndexPE][0];
			kmers[1] = tmp_kmers[maxOccIndexPE][1];
//			}

		return SUCCESSFUL;

	}
	//***********************************************************************************
	else if(validKmerNum==0)
	{
		kmers[0] = kmers[1] = NULL;
		return SUCCESSFUL;
	}
	else
	{
		//更新kmerseq，并返回得分最大的kmer
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				kmerSeqIntAssembly[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			}
			kmerSeqIntAssembly[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexPE) & lastEntryMask;

		kmers[0] = tmp_kmers[maxOccIndexPE][0];
		kmers[1] = tmp_kmers[maxOccIndexPE][1];

		return SUCCESSFUL;
	}
	//***********************************************************************************

	kmers[0] = kmers[1] = NULL;
	return SUCCESSFUL;
}

/**
 * Compute the kmer occNum.
 *  	If there some locked reads, then the occNum will be computed by considering the locked reads;
 *  	otherwise, consider all the reads in decision table.
 *
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
*/
int computeKmerOccNumByPE(kmertype *tmp_kmers[2], int *occNum, int contigNodesNum, int assemblyRound)
{
	//if(lockedReadsNum>=KOCKED_READS_NUM_THRESHOLD)
	if(lockedReadsNum>=lockedReadsNumThres)
	{ //如果决策表中含有锁定的reads, 则只考虑锁定的reads.
		if(computeKmerOccNumLockedByPE(tmp_kmers, occNum)==FAILED)
		{
			printf("In %s(), cannot compute kmer occNum in local assembly, error!\n", __func__);
			return FAILED;
		}
	}else
	{ //否则,考虑全部的reads.
		//计算决策表中的reads的得分
		if(computeKmerOccNumUnlockedByPE(tmp_kmers, occNum, contigNodesNum, assemblyRound)==FAILED)
		{
			printf("In %s(), cannot compute kmer occNum in local assembly, error!\n", __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * 按照一定的策略计算kmer得分.
 * 考虑全部的reads.
 *
 *   当前只考虑:
 *   	(1) 上次拼接出现的reads;
 *   	(2) 当次拼接出现, 上次拼接未出现并连续12次未出现的情况.
 *   	(3) ridpostable表中未考虑的reads的得分, 也即是新的reads的得分.
 */
int computeKmerOccNumUnlockedByPE(kmertype *tmp_kmers[2], int *occNum, int contigNodesNum, int assemblyRound)
{
	assemblingreadtype *this_assemblingRead = NULL;
	ridpostype *rid_pos = NULL;
	ridpostype *rid_pos_table = NULL;
	int posNum = 0;
	unsigned int this_rid = 0;
	unsigned short this_pos = 0;
	int i = 0, j = 0, properIndex = -1, startIndex = -1, exectIndex = -1, limitLastpos = -1, exitFlag = NO;

	*occNum = 0;

	for(i=0; i<itemNumDecisionTable; i++)
	{ //顺序访问决策表中的每行
		//找合适的lastpos的read
		if(decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0 && decisionTable[i].lastpos>0)
		{
			properIndex = getProperIndex(decisionTable+i, decisionTable, itemNumDecisionTable);

			// ############################ Debug information ##############################
			if(properIndex<0)
			{
				printf("line=%d, In %s(), properIndex=%d, i=%d, Error!\n", __LINE__, __func__, properIndex, i);
				//continue;
				return FAILED;
			}
			// ############################ Debug information ##############################

			this_assemblingRead = decisionTable + properIndex;
			limitLastpos = this_assemblingRead->lastpos;

			// ############################ Debug information ##############################
			//if(properIndex!=i)
			//{
			//	printf("line=%d, In %s(), properIndex=%d, i=%d\n", __LINE__, __func__, properIndex, i);
			//}
			// ############################ Debug information ##############################

			exitFlag = NO;
			while(exitFlag==NO)
			{
				for(j=0; j<2; j++) //j==0为正向kmer, j==1为反向互补kmer
				{
					if(tmp_kmers[j] && ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS)||(j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS)))
					{ //计算正向和反向互补kmer的得分
						//判断ridpostable表中该read是否存在
						rid_pos_table = tmp_kmers[j]->ppos;
						posNum = tmp_kmers[j]->arraysize;
						startIndex = findStartIndex(this_assemblingRead->rid, rid_pos_table, posNum);//二分查找rid_pos_table
						if(startIndex>=0)
						{  //存在, 继续查找精确位置
							if(this_assemblingRead->lastpos>0)
							{ //该read上次拼接出现
								if(j==0) //正向kmer
									exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->lastpos+1, startIndex, rid_pos_table, posNum);
								else //反向互补kmer
									exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->lastpos-1, startIndex, rid_pos_table, posNum);
							}else
							{  //该read上次拼接未出现
								if(j==0) //正向kmer
									exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->firstpos+this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
								else //反向互补kmer
									exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->firstpos-this_assemblingRead->kmerappeartimes-this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
							}

							if(exectIndex>=0)
							{  //该read位置合理
								rid_pos = &rid_pos_table[exectIndex];
								if(rid_pos->delsign==0) //该read未被删除
								{
									this_rid = rid_pos->rid;  //取得read的rid
									this_pos = rid_pos->pos;  //取得pos
									if(this_assemblingRead->lastpos==0)
									{ //该read上次拼接未出现
										if(this_assemblingRead->kmerunappeartimes==kmerSize)
										{
											(*occNum) ++;
										}else  //连续未出现的次数不等于kmer_size，说明该read拼接不合理
										{
											continue;
										}
									}else
									{
										(*occNum) ++;

									}
								}
							}
							exitFlag = YES;
						}else
						{
							exitFlag = YES;
						}
					}
				} //end for(j)


				if(exitFlag==NO)
				{
					//该read位置不合理, 寻找合理的位置
					properIndex = getProperIndexLimited(decisionTable+i, decisionTable, itemNumDecisionTable, limitLastpos);
					if(properIndex<0)
					{ //没有了合适的reads, 退出while循环
						exitFlag = YES; //退出标记置为YES
						break;
					}
					this_assemblingRead = decisionTable + properIndex;
					limitLastpos = this_assemblingRead->lastpos;
				}

			} //end while(exitFlag)
		} //end if(reserved)
	}// end for(i)

#if 1
	int returnCode;
	if(tmp_kmers[0])
	{ //正向kmer不为空
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				returnCode = validReadPair(rid_pos_table[i].rid);
				if(returnCode==YES)
				{
					(*occNum) ++;
				}else if(returnCode==ERROR)
				{
					printf("In %s(), cannot valid read pair, error!\n", __func__);
					return FAILED;
				}
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}

	if(tmp_kmers[1])
	{ //反向互补的kmer不为空
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>readLen-kmerSize+1-errorRegLenEnd3)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				returnCode = validReadPair(rid_pos_table[i].rid);
				if(returnCode==YES)
				{
					(*occNum) ++;
				}else if(returnCode==ERROR)
				{
					printf("In %s(), cannot valid read pair, error!\n", __func__);
					return FAILED;
				}
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0; //决策表中保留位清零

	return SUCCESSFUL;
}

/**
 * 按照一定的策略计算kmer得分.
 * 决策表中含有锁定的reads, 则只考虑锁定的reads.
 *
 *   当前只考虑:
 *   	(1) 上次拼接出现的reads;
 *   	(2) 当次拼接出现, 上次拼接未出现并连续12次未出现的情况.
 *   	(3) ridpostable表中未考虑的reads的得分, 也即是新的reads的得分.
 */

int computeKmerOccNumLockedByPE(kmertype *tmp_kmers[2], int *occNum)
{
	assemblingreadtype *this_assemblingRead = NULL;
	ridpostype *rid_pos = NULL;
	ridpostype *rid_pos_table = NULL;
	int posNum = 0;
	unsigned int this_rid = 0;
	unsigned short this_pos = 0;
	int i, j, properIndex, startIndex, exectIndex, limitLastpos = -1, exitFlag = NO;

	*occNum = 0;

	for(i=0; i<itemNumDecisionTable; i++)
	{ //顺序访问决策表中的每行
		//找合适的lastpos的read, 并计算得分
		if(decisionTable[i].locked)
		{
			//找合适的lastpos的read
			if(decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0 && decisionTable[i].lastpos>0)
			{
				properIndex = getProperIndex(decisionTable+i, decisionTable, itemNumDecisionTable);

				// ############################ Debug information ##############################
				if(properIndex<0)
				{
					printf("line=%d, In %s(), properIndex=%d, i=%d, Error!\n", __LINE__, __func__, properIndex, i);
					//continue;
					return FAILED;
				}
				// ############################ Debug information ##############################

				this_assemblingRead = decisionTable + properIndex;
				limitLastpos = this_assemblingRead->lastpos;

				// ############################ Debug information ##############################
				//if(properIndex!=i)
				//{
				//	printf("line=%d, In %s(), properIndex=%d, i=%d\n", __LINE__, __func__, properIndex, i);
				//}
				// ############################ Debug information ##############################

				exitFlag = NO;
				while(exitFlag==NO)
				{
					for(j=0; j<2; j++)
					{ //j==0为正向kmer, j==1为反向互补kmer
						if(tmp_kmers[j] && ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS) || (j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS)))
						{ //计算正向和反向互补kmer的得分
							//判断ridpostable表中该read是否存在
							rid_pos_table = tmp_kmers[j]->ppos;
							posNum = tmp_kmers[j]->arraysize;
							startIndex = findStartIndex(this_assemblingRead->rid, rid_pos_table, posNum);//二分查找rid_pos_table
							if(startIndex>=0)
							{ //存在, 继续查找精确位置
								if(this_assemblingRead->lastpos>0)
								{ //该read上次拼接出现
									if(j==0) //正向kmer
										exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->lastpos+1, startIndex, rid_pos_table, posNum);
									else //反向互补kmer
										exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->lastpos-1, startIndex, rid_pos_table, posNum);
								}else
								{  //该read上次拼接未出现
									if(j==0) //正向kmer
										exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->firstpos+this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
									else //反向互补kmer
										exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->firstpos-this_assemblingRead->kmerappeartimes-this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
								}

								if(exectIndex>=0)
								{ //该read位置合理
									rid_pos = &rid_pos_table[exectIndex];
									if(rid_pos->delsign==0) //该read未被删除
									{
										this_rid = rid_pos->rid;  //取得read的rid
										this_pos = rid_pos->pos;  //取得pos
										if(this_assemblingRead->lastpos==0)
										{ //该read上次拼接未出现
											if(this_assemblingRead->kmerunappeartimes==kmerSize)
											{
												(*occNum) ++;
											}else  //连续未出现的次数不等于kmer_size，说明该read拼接不合理
											{
												continue;
											}
										}else
										{
											(*occNum) ++;
										}
									}
									exitFlag = YES; //退出标记置为YES
								}
							}else
							{
								exitFlag = YES;
							}
						}

					} //end for(j)

					if(exitFlag==NO)
					{
						//该read位置不合理, 寻找合理的位置
						properIndex = getProperIndexLimited(decisionTable+i, decisionTable, itemNumDecisionTable, limitLastpos);
						if(properIndex<0)
						{ //没有了合适的reads, 退出while循环
							exitFlag = YES; //退出标记置为YES
							break;
						}
						this_assemblingRead = decisionTable + properIndex;
						limitLastpos = this_assemblingRead->lastpos;
					}
				} //end while(exitFlag)

			} //end if(reserved)
		} //end if (locked)
	}// end for(i)

#if 1
	int returnCode;
	if(tmp_kmers[0])
	{ //正向kmer不为空
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				returnCode = validReadPair(rid_pos_table[i].rid);
				if(returnCode==YES)
				{
					(*occNum) ++;
				}else if(returnCode==ERROR)
				{
					printf("In %s(), cannot valid read pair, error!\n", __func__);
					return FAILED;
				}
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}

	if(tmp_kmers[1])
	{ //反向互补的kmer不为空
		rid_pos_table = tmp_kmers[1]->ppos;
		posNum = tmp_kmers[1]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>readLen-kmerSize+1-errorRegLenEnd3)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				returnCode = validReadPair(rid_pos_table[i].rid);
				if(returnCode==YES)
				{
					(*occNum) ++;
				}else if(returnCode==ERROR)
				{
					printf("In %s(), cannot valid read pair, error!\n", __func__);
					return FAILED;
				}
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0; //决策表中保留位清零

	return SUCCESSFUL;
}

/**
 * Compute the occNum according to long K-mers.
 * 	@return:
 * 		If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int computeLongKmerOccNumByPE(kmertype *tmp_kmers[2], int *occNum, int length_k, int contigNodesNum, int assemblyRound)
{
	assemblingreadtype *this_assemblingRead = NULL;
	ridpostype *rid_pos = NULL;
	ridpostype *rid_pos_table = NULL;
	int posNum = 0;
	unsigned int this_rid = 0;
	unsigned short this_pos = 0;
	int i = 0, j = 0, properIndex = -1, startIndex = -1, exectIndex = -1, limitLastpos = -1, exitFlag = NO;

	*occNum = 0;

	for(i=0; i<itemNumDecisionTable; i++)
	{ //顺序访问决策表中的每行
		//找合适的lastpos的read
		if((decisionTable[i].matedFlag==YES && decisionTable[i].reserved==0 && decisionTable[i].lastpos>0)
				&& (decisionTable[i].kmerunappeartimes==0 && decisionTable[i].kmerunappearblocks==0)
				&& ((decisionTable[i].orientation==ORIENTATION_PLUS && decisionTable[i].lastpos>length_k-kmerSize)
						||(decisionTable[i].orientation==ORIENTATION_MINUS && decisionTable[i].lastpos<readLen-kmerSize-(length_k-kmerSize)+1)))
		{
			properIndex = getProperIndex(decisionTable+i, decisionTable, itemNumDecisionTable);

			// ############################ Debug information ##############################
			if(properIndex<0)
			{
				printf("line=%d, In %s(), properIndex=%d, i=%d, Error!\n", __LINE__, __func__, properIndex, i);
				//continue;
				return FAILED;
			}
			// ############################ Debug information ##############################

			this_assemblingRead = decisionTable + properIndex;
			limitLastpos = this_assemblingRead->lastpos;

			// ############################ Debug information ##############################
			//if(properIndex!=i)
			//{
			//	printf("line=%d, In %d(), properIndex=%d, i=%d\n", __LINE__, __func__, properIndex, i);
			//}
			// ############################ Debug information ##############################

			exitFlag = NO;
			while(exitFlag==NO)
			{
				for(j=0; j<2; j++) //j==0为正向kmer, j==1为反向互补kmer
				{
					if(tmp_kmers[j] && (this_assemblingRead->kmerunappeartimes==0 && this_assemblingRead->kmerunappearblocks==0)
							&& ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS && this_assemblingRead->lastpos>length_k-kmerSize)
									||(j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS && this_assemblingRead->lastpos<readLen-kmerSize-(length_k-kmerSize)+1)))
					{ //计算正向和反向互补kmer的得分
						//判断ridpostable表中该read是否存在
						rid_pos_table = tmp_kmers[j]->ppos;
						posNum = tmp_kmers[j]->arraysize;
						startIndex = findStartIndex(this_assemblingRead->rid, rid_pos_table, posNum);//二分查找rid_pos_table
						if(startIndex>=0)
						{  //存在, 继续查找精确位置
							if(this_assemblingRead->lastpos>0)
							{ //该read上次拼接出现
								if(j==0) //正向kmer
									exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->lastpos+1, startIndex, rid_pos_table, posNum);
								else //反向互补kmer
									exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->lastpos-1, startIndex, rid_pos_table, posNum);
							}else
							{  //该read上次拼接未出现
								if(j==0) //正向kmer
									exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->firstpos+this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
								else //反向互补kmer
									exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->firstpos-this_assemblingRead->kmerappeartimes-this_assemblingRead->kmerunappeartimes, startIndex, rid_pos_table, posNum);
							}

							if(exectIndex>=0)
							{  //该read位置合理
								rid_pos = rid_pos_table + exectIndex;
								this_rid = rid_pos->rid;  //取得read的rid
								this_pos = rid_pos->pos;  //取得pos
								if((this_assemblingRead->orientation==ORIENTATION_PLUS && this_pos>length_k-kmerSize)
									|| (this_assemblingRead->orientation==ORIENTATION_MINUS && this_pos<readLen-kmerSize-(length_k-kmerSize)+1))
								{
									if(rid_pos->delsign==0) //该read未被删除
									{
										if(this_assemblingRead->lastpos==0)
										{ //该read上次拼接未出现
											if(this_assemblingRead->kmerunappeartimes==kmerSize)
											{
												(*occNum) ++;
											}else  //连续未出现的次数不等于kmer_size，说明该read拼接不合理
											{
												continue;
											}
										}else if(this_assemblingRead->kmerunappeartimes==0)
										{
											(*occNum) ++;
										}
									}
								}
							}
							exitFlag = YES;
						}else
						{
							exitFlag = YES;
						}
					}
				} //end for(j)

				if(exitFlag==NO)
				{
					//该read位置不合理, 寻找合理的位置
					properIndex = getProperIndexLimited(decisionTable+i, decisionTable, itemNumDecisionTable, limitLastpos);
					if(properIndex<0)
					{ //没有了合适的reads, 退出while循环
						exitFlag = YES; //退出标记置为YES
						break;
					}
					this_assemblingRead = decisionTable + properIndex;
					limitLastpos = this_assemblingRead->lastpos;
				}

			} //end while(exitFlag)
		} //end if(reserved)
	}// end for(i)

#if 1
	int returnCode;
	//if(length_k<kmerSize+longKmerStepSize)
	if(length_k==kmerSize)
	{
		if(tmp_kmers[0])
		{ //正向kmer不为空
			rid_pos_table = tmp_kmers[0]->ppos;
			posNum = tmp_kmers[0]->arraysize;
			//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
			for(i=0; i<posNum; i++)//遍历ridpostable表
			{
				if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
				{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
					returnCode = validReadPair(rid_pos_table[i].rid);
					if(returnCode==YES)
					{
						(*occNum) ++;
					}else if(returnCode==ERROR)
					{
						printf("In %s(), cannot valid read pair, error!\n", __func__);
						return FAILED;
					}
				}
				//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
			}
		}

		if(tmp_kmers[1])
		{ //反向互补的kmer不为空
			rid_pos_table = tmp_kmers[1]->ppos;
			posNum = tmp_kmers[1]->arraysize;
			//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
			for(i=0; i<posNum; i++)//遍历ridpostable表
			{
				if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos>readLen-kmerSize+1-errorRegLenEnd3)
				{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
					returnCode = validReadPair(rid_pos_table[i].rid);
					if(returnCode==YES)
					{
						(*occNum) ++;
					}else if(returnCode==ERROR)
					{
						printf("In %s(), cannot valid read pair, error!\n", __func__);
						return FAILED;
					}
				}
				//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
			}
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0; //决策表中保留位清零

	return SUCCESSFUL;
}

/**
 * Check whether the read has valid read pair.
 *  @return:
 *  	(1) If the read has valid read pair, return YES; if invalid, return NO.
 *  	(2) If errors occurred, return ERROR.
 */
short validReadPair(uint64_t readID)
{
	uint64_t readID_paired;
	int readOrient_paired, contigPos_paired;
	int fragmentSize, validReadOrient, validReadOrient_paired;
	PERead_t *pPERead;

	// generate the paired readID
	if(readID%2==1)
	{ // odd --> even
		readID_paired = readID + 1;
	}else
	{ // even --> odd
		readID_paired = readID - 1;
	}

	if(getReadFromPEHashtable(&pPERead, readID_paired)==FAILED)
	{
		printf("In %s(), cannot get the paired read %lu, error!\n", __func__, readID_paired);
		return ERROR;
	}

	if(pPERead)
	{
		return YES;
/*
		contigPos_paired = pPERead->cpos;
		readOrient_paired = pPERead->orient;

		// get the contig sequence length and valid reads orientations
		if(assemblyRound==FRIST_ROUND_ASSEMBLY)
		{ // the first round
			validReadOrient = ORIENTATION_MINUS;
			validReadOrient_paired = ORIENTATION_PLUS;
		}else
		{ // the second round
			validReadOrient = ORIENTATION_MINUS;
			validReadOrient_paired = ORIENTATION_PLUS;
		}

		if(readOrient==validReadOrient && readOrient_paired==validReadOrient_paired)
		{
			fragmentSize = contigNodesNum - contigPos_paired + 1 + READ_LEN - KMER_SIZE + 1;
			if(fragmentSize<=meanSizeInsert+standardDevFactor*standardDev && fragmentSize>=meanSizeInsert-standardDevFactor*standardDev)
			{ // m - 5*sdev <= fragSize <= m + 5*sdev
				return YES;
			}
		}
*/
	}

	return NO;
}

/**
 * Trim a read length of contig nodes at tail.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int trimContigTailByReadLen(contigtype *contighead, contigtype **contigtail, contigtype **successContig, int *contigNodesNum, int assemblyRound)
{
	int i, j, divContigIndex, startCheckContigIndex, ridposnum, newSuccessContigIndex;
	contigtype *contig, *startCheckContig, *newSuccessContig;
	successRead_t *ridposorient;
	char seq[1000] = {0}, reversed_seq[1000] = {0}; //正向的碱基序列和反向互补的碱基序列
	int seq_len, seqIndex;

	// ######################### Debug information ##############################
	//printf("Before trim tail: contigNodesNum=%d, assemblyRound=%d\n", *contigNodesNum, assemblyRound);
	// ######################### Debug information ##############################

	divContigIndex = *contigNodesNum - readLen + (contighead->index - 1);
	startCheckContigIndex = *contigNodesNum - readLen * 3 + 1 + (contighead->index - 1);
	if(startCheckContigIndex<=contighead->index)
		startCheckContigIndex = contighead->index;

	contig = contighead;
	while(contig)
	{
		if(contig->index==startCheckContigIndex)
			break;
		contig = contig->next;
	}
	startCheckContig = contig;

	i = 0;
	contig = startCheckContig;
	while(contig)
	{
		switch(contig->base)
		{
			case 0: seq[i] = 'A'; break;
			case 1: seq[i] = 'C'; break;
			case 2: seq[i] = 'G'; break;
			case 3: seq[i] = 'T'; break;
			default:
				printf("line=%d, In %s(), error base int: %d\n", __LINE__, __func__, contig->base);
				return FAILED;
		}

		i ++;

		// ############################ Debug information ##############################
		if(i>=1000)
		{
			printf("line=%d, In %s(), such a long contig here, contig len>=%d, error!\n", __LINE__, __func__, i);
			return FAILED;
		}
		// ############################ Debug information ##############################

		contig = contig->next;
	}
	seq[i] = '\0';
	seq_len = i;

	if(getReversedSeq(reversed_seq, seq, seq_len)==FAILED)
	{ //取得该contig的反向互补的碱基序列
		printf("line=%d, In %s(), cannot get the reverse contig sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first assembly round
		seqIndex = 0;
		contig = startCheckContig;
		while(contig)
		{
			if(contig->index>divContigIndex && contig->ridposnum>0)
			{ // remove the reads from the contig nodes list
/*
//				deledReadNum = 0;
				ridposnum = contig->ridposnum;
				ridposorient = contig->pridposorientation;
				for(i=0; i<ridposnum; i++)
				{
					// recover the read
					if(ridposorient[i].orientation==ORIENTATION_PLUS)
					{ //该read为正向, 实际为反向
						if(recoverReadFromGraph(reversed_seq+seq_len-1-seqIndex, ridposorient[i].rid, deBruijnGraph)==FAILED)
						{ //恢复该read
							printf("line=%d, In %s(), cannot recover the read (%d,%d,%d,%c), seqIndex=%d, error!\n", __LINE__, __func__, ridposorient[i].rid, ridposorient[i].pos, ridposorient[i].matchnum, ridposorient[i].orientation, seqIndex);
							return FAILED;
						}
					}else
					{ //该read为反向, 实际为正向
						if(recoverReadFromGraph(seq+seqIndex+1-ridposorient[i].matchnum, ridposorient[i].rid, deBruijnGraph)==FAILED)
						{ //恢复该read
							printf("line=%d, In %s(), cannot recover the read (%d,%d,%d,%c), seqIndex=%d, error!\n", __LINE__, __func__, ridposorient[i].rid, ridposorient[i].pos, ridposorient[i].matchnum, ridposorient[i].orientation, seqIndex);
							return FAILED;
						}
					}

//					// remove the read from the contig node list
//					for(j=i+1; j<ridposnum; j++)
//					{
//						if(memcpy(ridposorient+j-1, ridposorient+j, sizeof(successRead_t))==NULL)
//						{
//							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
//							return FAILED;
//						}
//					}
//					deledReadNum ++;
				}
*/

				contig->ridposnum = 0;
				free(contig->pridposorientation);
				contig->pridposorientation = NULL;

//				contig->ridposnum -= deledReadNum;
//				if(contig->ridposnum==0)
//				{ // free the memory
//					free(contig->pridposorientation);
//					contig->pridposorientation = NULL;
//				}else if(contig->ridposnum<0)
//				{ // error
//					printf("line=%d, In %s(), ridposnum=%d, error!\n", __LINE__, __func__, contig->ridposnum);
//					return FAILED;
//				}
			} // end if(contig->index>divContigIndex && contig->ridposnum>0)

			seqIndex ++;
			contig = contig->next;
		}
	}else
	{ // the second assembly round
		seqIndex = 0;
		contig = startCheckContig;
		while(contig)
		{
			if(contig->ridposnum>0)
			{
				ridposnum = contig->ridposnum;
				ridposorient = contig->pridposorientation;
				i = 0;
				while(i<ridposnum)
				{
					if(contig->index+ridposorient[i].matchnum-1>divContigIndex)
					{
						// remove the read from the contig node list
						for(j=i+1; j<ridposnum; j++)
						{
							if(memcpy(ridposorient+j-1, ridposorient+j, sizeof(successRead_t))==NULL)
							{
								printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
								return FAILED;
							}
						}

						ridposnum --;
						continue;
					}

					i ++;
				}

				contig->ridposnum = ridposnum;
				if(contig->ridposnum==0)
				{ // free the memory
					free(contig->pridposorientation);
					contig->pridposorientation = NULL;
				}else if(contig->ridposnum<0)
				{ // error
					printf("line=%d, In %s(), ridposnum=%d, error!\n", __LINE__, __func__, contig->ridposnum);
					return FAILED;
				}
			}

			seqIndex ++;
			contig = contig->next;
		}
	}

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ // the first assembly round
		// get newSuccessContig
		newSuccessContig = NULL;
		newSuccessContigIndex = 0;
		contig = startCheckContig;
		while(contig)
		{
			if(contig->ridposnum>0 && contig->index>newSuccessContigIndex)
			{
				newSuccessContig = contig;
				newSuccessContigIndex = contig->index;
			}
			contig = contig->next;
		}

		// ##################### Debug information #######################
		if(newSuccessContigIndex>divContigIndex)
		{
			printf("line=%d, In %s(), newSuccessContigIndex=%d > divContigIndex=%d, error!\n", __LINE__, __func__, newSuccessContigIndex, divContigIndex);
			return FAILED;
		}
		// ##################### Debug information #######################
	}else
	{ // the second assembly round
		// get newSuccessContig
		newSuccessContigIndex = 0;
		contig = startCheckContig;
		while(contig)
		{
			if(contig->ridposnum>0)
			{
				ridposnum = contig->ridposnum;
				ridposorient = contig->pridposorientation;
				for(i=0; i<ridposnum; i++)
				{
					if(contig->index+ridposorient[i].matchnum-1>newSuccessContigIndex)
					{
						newSuccessContigIndex = contig->index + ridposorient[i].matchnum - 1;
					}
				}
			}

			contig = contig->next;
		}

		// ##################### Debug information #######################
		if(newSuccessContigIndex>divContigIndex)
		{
			printf("line=%d, In %s(), newSuccessContigIndex=%d > divContigIndex=%d, error!\n", __LINE__, __func__, newSuccessContigIndex, divContigIndex);
			return FAILED;
		}
		// ##################### Debug information #######################

		if(newSuccessContigIndex==0)
		{
			newSuccessContig = NULL;
		}else if(newSuccessContigIndex>0)
		{
			newSuccessContig = NULL;
			contig = startCheckContig;
			while(contig)
			{
				if(contig->index==newSuccessContigIndex)
					break;

				contig = contig->next;
			}
			newSuccessContig = contig;
		}else
		{
			printf("line=%d, In %s(), newSuccessContigIndex=%d, error!\n", __LINE__, __func__, newSuccessContigIndex);
			return FAILED;
		}
	}

	// free the memory of tail contig nodes
	if(newSuccessContig!=NULL)
	{
		// delete the contig nodes at tail
		contig = startCheckContig;
		while(contig)
		{
			if(contig==newSuccessContig)
			{
				releaseContig(contig->next);
				contig->next = NULL;
				break;
			}
			contig = contig->next;
		}
	}else
	{
		printf("line=%d, In %s(), newSuccessContig==NULL, error!\n", __LINE__, __func__);
		return FAILED;
	}

	*contigtail = newSuccessContig;
	*successContig = newSuccessContig;
	*contigNodesNum = (*contigtail)->index - (contighead->index - 1);

	// ######################### Debug information ##############################
	//printf("After trim tail: contigNodesNum=%d, assemblyRound=%d\n", *contigNodesNum, assemblyRound);
	// ######################### Debug information ##############################

	return SUCCESSFUL;
}

/**
 * Set the navigation occurrence queue to empty.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short setEmptyNaviOccQueue(double *naviOccQueuePara, int *itemNumNaviOccQueuePara, int *frontRowNaviOccQueuePara, int *rearRowNaviOccQueuePara)
{
	int i;

	for(i=0; i<*itemNumNaviOccQueuePara; i++) naviOccQueuePara[i] = 0;
	*itemNumNaviOccQueuePara = 0;
	*frontRowNaviOccQueuePara = 0;
	*rearRowNaviOccQueuePara = 0;

	return SUCCESSFUL;
}

/**
 * Update the navigation occurrence queue.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateNaviOccQueue(double *naviOccQueuePara, int *itemNumNaviOccQueuePara, int *frontRowNaviOccQueuePara, int *rearRowNaviOccQueuePara, double maxOccNum)
{
	if(*itemNumNaviOccQueuePara>=maxItemNumNaviOccQueue)
	{ // full queue, then add the first item
		(*rearRowNaviOccQueuePara) = ((*rearRowNaviOccQueuePara) + 1) % maxItemNumNaviOccQueue;
		(*frontRowNaviOccQueuePara) = ((*frontRowNaviOccQueuePara) + 1) % maxItemNumNaviOccQueue;
		naviOccQueuePara[*rearRowNaviOccQueuePara] = maxOccNum;
	}else if(*itemNumNaviOccQueuePara==0)
	{ // empty queue, then add the first item
		*frontRowNaviOccQueuePara = *rearRowNaviOccQueuePara = 0;
		naviOccQueuePara[*rearRowNaviOccQueuePara] = maxOccNum;
		*itemNumNaviOccQueuePara = 1;
	}else if(*itemNumNaviOccQueuePara<maxItemNumNaviOccQueue)
	{ // non-full queue, add new item
		(*rearRowNaviOccQueuePara) ++;
		naviOccQueuePara[*rearRowNaviOccQueuePara] = maxOccNum;
		(*itemNumNaviOccQueuePara) ++;
	}else
	{ // Exception errors
		printf("line=%d, In %s(), itemNumNaviOccQueue=%d, error!\n", __LINE__, __func__, *itemNumNaviOccQueuePara);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Compute the average number in navigation occurrence queue.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short calcAverOccNaviOccQueue(double *averOccNum, double *naviOccQueuePara, int itemNumNaviOccQueuePara)
{
	int i;

	*averOccNum = 0;
	for(i=0; i<itemNumNaviOccQueuePara; i++)
		*averOccNum += naviOccQueuePara[i];

	if(itemNumNaviOccQueuePara>0)
		*averOccNum /= itemNumNaviOccQueuePara;
	else
		*averOccNum = 0;

	return SUCCESSFUL;
}

/**
 * Get the length of low occurrence region in navigation occurrence queue.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getLowOccLenNaviOccQueue(int *lowLen, double *naviOccQueuePara, int itemNumNaviOccQueuePara, int frontRowNaviOccQueuePara)
{
	int i, j, maxLen, len;

	*lowLen = 0;
	for(i=0; i<itemNumNaviOccQueuePara; i++)
	{
		if(naviOccQueuePara[i]<=lowOccThresNaviOccQueue)
		{
			(*lowLen) ++;
		}
	}

/*
	len = 0;
	maxLen = 0;
	j = frontRowNaviOccQueuePara;
	for(i=0; i<itemNumNaviOccQueuePara; i++)
	{
		if(naviOccQueuePara[j]<=lowOccThresNaviOccQueue)
		{
			len ++;
		}else
		{
			if(len>maxLen)
				maxLen = len;
			len = 0;
		}
		j = (j+1) % maxItemNumNaviOccQueue;
	}

	*lowLen = maxLen;
*/

	return SUCCESSFUL;
}


