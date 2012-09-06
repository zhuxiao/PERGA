
#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Build contigs.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short buildContigs(char *contigFile, char *graphFile)
{
	printf("\n============= Begin building contigs, please wait ... =============\n");

	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);

	char hangingFile[256];
	int i, turnContigIndex;
	int64_t localContigID;
	double averOccNumNaviOccQueue;
	int lowOccNum;

	// load the graph to memory
	if(loadGraph(&deBruijnGraph, graphFile)==FAILED)
	{
		printf("line=%d, In %s(), cannot load graph to memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	// ############################ Debug information ##############################
	//if(checkGraph(deBruijnGraph)==FAILED)
	//{
	//	printf("line=%d, In %s(), checking graph error!\n", __LINE__, __func__);
	//	return FAILED;
	//}
	// ############################ Debug information ##############################

	initFirstKmerThreshold();
	if(initMemory()==FAILED)
	{
		printf("line=%d, In %s(), cannot init the memory of two tables, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//++++++++++++++++++++++++++++++++++++
	if(PEGivenType>=NONE_PE_GIVEN_TYPE)
	{
		// estimate the insert size and standard deviation of fragment library
		if(estimateInsertSizeAndSdev()==FAILED)
		{
			printf("line=%d, In %s(), cannot estimate the insert size and standard deviation of fragment library, error!\n", __LINE__, __func__);
			return FAILED;
		}
		// ############################ Debug information ##############################
		//return SUCCESSFUL;
		// ############################ Debug information ##############################

		if(initPEHashParas()==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the PEHash table parameters, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	fpContigsBase = fopen(contigFile, "w");
	if(fpContigsBase==NULL)
	{
		printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, contigFile);
		return FAILED;
	}
	if(hangingContigOutFlag==YES)
	{
		strcpy(hangingFile, contigFile);
		strcat(hangingFile, ".hang");
		fpContigsHanging = fopen(hangingFile, "w");
		if(fpContigsHanging==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, hangingFile);
			return FAILED;
		}
	}

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
//			if(localContigID==100 && contigIndex>=11860 && assemblyRound==FIRST_ROUND_ASSEMBLY)
//			{
//				printf("localContigID=%ld, contigID=%d, contigIndex=%d, assemblyRound=%d\n", localContigID, contigsNum+1, contigIndex, assemblyRound);
//			}
			// ############################ Debug information ##############################


			// initialize or update the PE hash table
			if(PEGivenType>NONE_PE_GIVEN_TYPE && contigIndex>=minContigLenUsingPE)
			{
				if(updatePEHashTable(contigIndex, assemblyRound)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, contigIndex=%d, cannot update the PE hash table, error!\n", __LINE__, __func__, localContigID, contigsNum+1, contigIndex);
					return FAILED;
				}

				// ########################### Debug information ########################
				//if(assemblyRound==SECOND_ROUND_ASSEMBLY && contigIndex-hashRegRightContig->index<minContigLenUsingPE-1)
				//{
				//	printf("line=%d, In %s(), assemblyRound=%d, contigIndex=%d, hashRegRightContig->index=%d, error!\n", __LINE__, __func__, assemblyRound, contigIndex, hashRegRightContig->index);
				//	return FAILED;
				//}
				// ########################### Debug information ########################

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
						printf("line=%d, In %s(), localContigID=%ld, cannot get next kmer, error!\n", __LINE__, __func__, localContigID);
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
						printf("line=%d, In %s(), localContigID=%ld, contigIndex=%d, cannnot check the reads number in reads number region, error!\n", __LINE__, __func__, localContigID, contigIndex);
						return FAILED;
					}
				}
			}

			if(kmers[0]==NULL && kmers[1]==NULL)
			{
				if(successContig==NULL)
				{ //没有成功的reads, 该contig拼接失败
					//printf("line=%d, In %s(), localContigID=%ld, contigsNum=%d, assemblyRound=%d, contigIndex=%d, the successContig==NULL!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
					break;
				}

//				printf("localContigID=%ld, assemblyRound=%d, contigIndex=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, contigIndex, itemNumDecisionTable);
//				if(navigationFlag==NAVI_PE_FLAG)
//					printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
//				else if(navigationFlag==NAVI_SE_FLAG)
//				{
//					printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
//					printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
//				}
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

			// ############################ Debug information ##############################
//			if(localContigID==100 && contigIndex>=11860 && assemblyRound==FIRST_ROUND_ASSEMBLY)
//			{
//				printf("localContigID=%ld, assemblyRound=%d, contigIndex=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, contigIndex, itemNumDecisionTable);
//				if(navigationFlag==NAVI_PE_FLAG)
//					printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
//				else if(navigationFlag==NAVI_SE_FLAG)
//				{
//					printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
//					printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
//				}
//				printf("\tdistance=%d, readsNumRatio=%.2f\n", contigIndex-successContig->index, readsNumRatio);
//			}
			// ############################ Debug information ##############################

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
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, contigIndex=%d, assemblyRound=%d, cannot delete the reads from graph, error!\n", __LINE__, __func__, localContigID, contigsNum+1, contigIndex, assemblyRound);
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
					//printf("line=%d, In %s(), localContigID=%ld, contigsNum=%d, assemblyRound=%d, contigIndex=%d, the tmp_successContig==NULL!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
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
			}else if(contigIndex-successContig->index > readLen-MIN_OVERLAP_LEN)
			{ //已经有成功的reads, 则根据拼接的情况, 确定是否需要继续拼接

				number_of_overlap_less_than_threshold ++;

//				printf("===localContigID=%ld, assemblyRound=%d, contigIndex=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, contigIndex, itemNumDecisionTable);
//				if(navigationFlag==NAVI_PE_FLAG)
//					printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
//				else if(navigationFlag==NAVI_SE_FLAG)
//				{
//					printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
//					printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
//				}
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

				// get the length of low occurrence region
				if(getLowOccLenNaviOccQueue(&lowOccNum, naviOccQueue, itemNumNaviOccQueue, frontRowNaviOccQueue)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, contigID=%d, cannot get the low occurrence number in navigation occurrence queue, error!\n", __LINE__, __func__, localContigID, contigsNum+1);
					return FAILED;
				}
				printf("#### contigIndex=%d, distance=%d, readsNumRatio=%.2f, averOccNumNaviOccQueue=%.2f, lowOccNum=%d\n", contigIndex, contigIndex-successContig->index, readsNumRatio, averOccNumNaviOccQueue, lowOccNum);

				//if((lowOccNum>2) || ((navigationFlag==NAVI_PE_FLAG && averOccNumNaviOccQueue<2*minKmerOccPE) || (navigationFlag==NAVI_SE_FLAG && averOccNumNaviOccQueue<2*minKmerOccSE)) || (readsNumRatio>2 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
				//if((lowOccNum>2) || (readsNumRatio>2 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
				//if((averOccNumNaviOccQueue<2.5*minKmerOccSE) || (readsNumRatio>2 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
				if((readsNumRatio>3 || readsNumRatio<0.3) || ((navigationFlag==NAVI_PE_FLAG && secondOccPE>=minKmerOccPE && secondOccPE/maxOccPE>=0.2) || (navigationFlag==NAVI_SE_FLAG && secondOccSE>=minKmerOccSE && secondOccSE/maxOccSE>=0.15)))
				{
					number_of_overlap_less_than_threshold ++;

					printf("===localContigID=%ld, assemblyRound=%d, contigIndex=%d, itemNumDecisionTable=%d\n", localContigID, assemblyRound, contigIndex, itemNumDecisionTable);
					if(navigationFlag==NAVI_PE_FLAG)
						printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
					else if(navigationFlag==NAVI_SE_FLAG)
					{
						printf("\toccsNumPE: (%d, %d, %d, %d)\n", occsNumPE[0], occsNumPE[1], occsNumPE[2], occsNumPE[3]);
						printf("\toccsNumSE: (%d, %d, %d, %d)\n", occsNumSE[0], occsNumSE[1], occsNumSE[2], occsNumSE[3]);
					}
					printf("\tdistance=%d, readsNumRatio=%.2f, averOccNumNaviOccQueue=%.2f, lowOccNum=%d\n", contigIndex-successContig->index, readsNumRatio, averOccNumNaviOccQueue, lowOccNum);

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
			// ############################ Debug information ##############################
//			if(localContigID==98)
//			{
//				outputContigToTmpFile(contighead, HANGING_READ_TYPE_CONTIG_FILE);
//			}
			// ############################ Debug information ##############################

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
			if(contigIndex>=minContigLen)
			{ //contig长度大于100, 则写进文件中
				contigsNum ++;

				// ############################ Debug information ##############################
				//printf("contigID=%d, contigLen=%d, turnContigIndex=%d.\n", contigsNum, contigIndex, turnContigIndex);
				// ############################ Debug information ##############################

				successReadNum += this_successReadNum;

				// output contig nodes to file
				if(outputContigToFile(fpContigsBase, BASE_TYPE_FASTA_CONTIG_FILE, contighead, contigsNum, contigIndex)==FAILED)
				{
					printf("line=%d, In %s(), localContigID=%ld, cannot output contig nodes to file, error!\n", __LINE__, __func__, localContigID);
					return FAILED;
				}
				if(hangingContigOutFlag==YES)
				{
					if(outputContigToFile(fpContigsHanging, HANGING_READ_TYPE_CONTIG_FILE, contighead, contigsNum, contigIndex)==FAILED)
					{
						printf("line=%d, In %s(), localContigID=%ld, cannot output contig nodes to file, error!\n", __LINE__, __func__, localContigID);
						return FAILED;
					}
				}

				basesNum += contigIndex;
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
			}
		}else
		{
			// ############################ Debug information ##############################
			//if(contigIndex>=CONTIG_LEN_THRESHOLD)
			//{
			//	printf("line=%d, In %s(), localContigID=%ld, contigsNum=%d, assemblyRound=%d, contigIndex=%d, successContig==NULL, Error!\n", __LINE__, __func__, localContigID, contigsNum+1, assemblyRound, contigIndex);
			//	return FAILED;
			//}
			// ############################ Debug information ##############################
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

		//释放该contig拼接过程中占用的内存，并初始化新的contig
		releaseContig(contighead);

	} //end while(kmerIndex < TABLE_SIZE_DE_BRUIJN)

	fclose(fpContigsBase);
	fpContigsBase = NULL;
	if(hangingContigOutFlag==YES)
	{
		fclose(fpContigsHanging);
		fpContigsHanging = NULL;
	}

	freeMemory();
	releaseGraph(deBruijnGraph);  //free graph

	printf("contigsNum=%d, basesNum=%ld\n", contigsNum, basesNum);


	gettimeofday(&tp_end,NULL);
	time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;
	printf("Assembly Used Time: %f Seconds.\n", time_used);

	printf("============= End build contigs. =============\n");

	return SUCCESSFUL;
}


/**
 * Initialize the memory for assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initMemory()
{
	//longKmerSize = ceil((readLen - 2*errorRegLenEnd3 - kmerSize) * LONG_KMER_SIZE_FACTOR) + kmerSize;
	longKmerSize = ceil((readLen - 1.5*errorRegLenEnd3 - kmerSize) * LONG_KMER_SIZE_FACTOR) + kmerSize;
	//longKmerSize = ceil((readLen - errorRegLenEnd3 - kmerSize) * LONG_KMER_SIZE_FACTOR) + kmerSize;
	//longKmerSize = ceil((readLen - errorRegLenEnd3) * 0.9);

	if((longKmerSize & 1) == 0)
		longKmerSize --;
	longKmerStepSize = floor((longKmerSize - kmerSize) / 4.0);
	if((longKmerStepSize & 1) == 1)
		longKmerStepSize ++;
	if(longKmerStepSize<1)
		longKmerStepSize = 2;

	printf("averKmerOcc=%.2f\n", averKmerOcc);
	printf("longKmerSize=%d, longKmerStepSize=%d\n", longKmerSize, longKmerStepSize);

	if(PEGivenType>=NONE_PE_GIVEN_TYPE)
	{
		minKmerOccPE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);
		minKmerOccSE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);

		if(minKmerOccPE<MIN_KMER_OCC_THRES)
			minKmerOccPE = MIN_KMER_OCC_THRES;
		if(minKmerOccSE<MIN_KMER_OCC_THRES)
			minKmerOccSE = MIN_KMER_OCC_THRES;

		//maxSecondOcc = minKmerOccSE * OCCS_NUM_FACTOR;
		//maxSecondOcc = ceil(minKmerOccSE * MAX_SECOND_OCC_FACTOR);
		maxSecondOcc = averKmerOcc;
		//maxFirstOcc = ceil(minKmerOccSE * MAX_FIRST_OCC_FACTOR);
		//minLongKmerOcc = floor(minKmerOccSE * LONG_KMER_OCC_FACTOR);
		minLongKmerOcc = floor(averKmerOcc * 0.5);
		minReadsNumPEHashThres = ceil(averKmerOcc * MIN_READ_NUM_PE_HASH_FACTOR);

		//if(maxSecondOcc>MAX_SECOND_OCC_THRES)
		//{
		//	maxSecondOcc = MAX_SECOND_OCC_THRES;
		//}

		if(minLongKmerOcc<minKmerOccSE)
			minLongKmerOcc = minKmerOccSE;
		if(minLongKmerOcc>MIN_LONG_KMER_OCC_THRES)
		{
			minLongKmerOcc = MIN_LONG_KMER_OCC_THRES;
		}


		maxOccNumFaiedPE = ceil(OCCS_NUM_SE_FAILED_PE_FACTOR * averKmerOcc);
		if(maxOccNumFaiedPE>MAX_OCC_NUM_FAILED_PE_THRES)
			maxOccNumFaiedPE = MAX_OCC_NUM_FAILED_PE_THRES;
//		else
//			maxOccNumFaiedPE *= 0.8;
		//maxNavigationNumSE = MAX_NAVI_NUM_SE_THRES;

		printf("minKmerOccSE=%.2f, minKmerOccPE=%.2f\n", minKmerOccSE, minKmerOccPE);
		printf("maxSecondOcc=%.2f\n", maxSecondOcc);
		//printf("maxFirstOcc=%.2f\n", maxFirstOcc);
		printf("minLongKmerOcc=%.2f\n", minLongKmerOcc);
		printf("minReadsNumPEHashThres=%.2f\n", minReadsNumPEHashThres);
		printf("maxOccNumFaiedPE=%.2f\n", maxOccNumFaiedPE);
		//printf("maxNavigationNumSE=%d\n", maxNavigationNumSE);
	}else
	{
		minKmerOccSE = ceil(averKmerOcc * MIN_KMER_OCC_FACTOR);

		if(minKmerOccSE<MIN_KMER_OCC_THRES)
			minKmerOccSE = MIN_KMER_OCC_THRES;

		//maxSecondOcc = minKmerOccSE * OCCS_NUM_FACTOR;
		//maxSecondOcc = ceil(minKmerOccSE * MAX_SECOND_OCC_FACTOR);
		maxSecondOcc = averKmerOcc;
		//maxFirstOcc = ceil(minKmerOccSE * MAX_FIRST_OCC_FACTOR);
		//minLongKmerOcc = floor(minKmerOccSE * LONG_KMER_OCC_FACTOR);
		minLongKmerOcc = floor(averKmerOcc * 0.5);


		//if(maxSecondOcc>MAX_SECOND_OCC_THRES)
		//{
		//	maxSecondOcc = MAX_SECOND_OCC_THRES;
		//}

		if(minLongKmerOcc<minKmerOccSE)
			minLongKmerOcc = minKmerOccSE;
		if(minLongKmerOcc>MIN_LONG_KMER_OCC_THRES)
		{
			minLongKmerOcc = MIN_LONG_KMER_OCC_THRES;
		}

		printf("minKmerOccSE=%.2f\n", minKmerOccSE);
		printf("maxSecondOcc=%.2f\n", maxSecondOcc);
		//printf("maxFirstOcc=%.2f\n", maxFirstOcc);
		printf("minLongKmerOcc=%.2f\n", minLongKmerOcc);
	}

	lockedReadsNumThres = averKmerOcc;
	printf("lockedReadsNumThres=%.2f\n", lockedReadsNumThres);

	// the global variables of reads number region
	maxRegLenReadsNumReg = ceil(readLen * REG_LEN_READS_NUM_REG_FACTOR);
	minContigLenCheckingReadsNum = readLen + maxRegLenReadsNumReg;
	maxReadsNumRatioThres = MAX_READS_NUM_RATIO_THRES;
	minReadsNumRatioThres = MIN_READS_NUM_RATIO_THRES;
	solvedRepeatsNum = 0;

	lowOccThresNaviOccQueue = 0.4 * averKmerOcc;

	printf("maxRegLenReadsNumReg=%d\n", maxRegLenReadsNumReg);
	printf("minContigLenCheckingReadsNum=%d\n", minContigLenCheckingReadsNum);
	printf("maxReadsNumRatioThres=%.2f\n", maxReadsNumRatioThres);
	printf("minReadsNumRatioThres=%.2f\n", minReadsNumRatioThres);
	printf("lowOccThresNaviOccQueue=%.2f\n", lowOccThresNaviOccQueue);


	hangingContigOutFlag = NO;
	maxItemNumDecisionTable = TABLE_SIZE_ASSEMBLINGREAD;
	//决策表相关变量初始化
	decisionTable = (assemblingreadtype*) malloc(maxItemNumDecisionTable*sizeof(assemblingreadtype));
	if(decisionTable==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	maxItemNumSuccessReadsArr = TABLE_SIZE_RIDPOSORIENTATION;
	successReadsArr = (successRead_t*) malloc(maxItemNumSuccessReadsArr*sizeof(successRead_t));
	if(successReadsArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//++++++++++++++++++++++++++++++++++++
	PEHashArr = (PERead_t **) calloc(TABLE_SIZE_HASH_PE, sizeof(PERead_t **));
	if(PEHashArr==NULL)
	{
		printf("In %s(), cannot allocate memory, error!\n", __func__);
		return FAILED;
	}

	lastseq36 = (char*) malloc((readLen+1)*sizeof(char));
	if(lastseq36==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	kmerSeqIntAssembly = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(kmerSeqIntAssembly==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	kmerSeqIntAssemblyRev = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(kmerSeqIntAssemblyRev==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	tmpKmerSeqIntAssembly = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(tmpKmerSeqIntAssembly==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//maxItemNumNaviOccQueue = MAX_ITEM_NUM_NAVI_OCC_QUEUE;
	maxItemNumNaviOccQueue = readLen * 0.3;
	naviOccQueue = (double*) malloc(maxItemNumNaviOccQueue*sizeof(double));
	if(naviOccQueue==NULL)
	{
		printf("line=%d, In %s(), cannot allocate the memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	int i;
	for(i=0; i<4; i++) tabooSeqInt[i] = 0;
	for(i=0; i<32; i++)
	{
		tabooSeqInt[1] = (tabooSeqInt[1] << 2) | 1;
		tabooSeqInt[2] = (tabooSeqInt[2] << 2) | 2;
		tabooSeqInt[3] = (tabooSeqInt[3] <<  2) | 3;
	}

	return SUCCESSFUL;
}

void freeMemory()
{
	maxItemNumDecisionTable = 0;
	maxItemNumSuccessReadsArr = 0;

	free(decisionTable);
	decisionTable = NULL;

	free(successReadsArr);
	successReadsArr = NULL;

	free(PEHashArr);
	PEHashArr = NULL;

	free(lastseq36);
	lastseq36 = NULL;

	free(kmerSeqIntAssembly);
	kmerSeqIntAssembly = NULL;

	free(kmerSeqIntAssemblyRev);
	kmerSeqIntAssemblyRev = NULL;

	free(tmpKmerSeqIntAssembly);
	tmpKmerSeqIntAssembly = NULL;

	free(naviOccQueue);
	naviOccQueue = NULL;
}


/**
 * 初始化contig.
 *   @ return:
 *     成功, 返回成功标记; 失败, 返回失败标记.
 */
short initContig(contigtype **contighead,contigtype **tailContig)
{
	int i;
	contigtype *contig, *pre_contig;

	pre_contig = NULL;
	i = 1;
	for(; i<=kmerSize; i++)
	{
		contig = (contigtype *)malloc(sizeof(contigtype));
		if(contig==NULL)
		{
			printf("out of memory in function initContig()!\n");
			return FAILED;
		}
		contig->index = i;
		if((i-1)/32==entriesPerKmer-1)
			contig->base = (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*(lastEntryBaseNum-i))) & 3;
		else
			contig->base = (kmerSeqIntAssembly[(i-1)/32] >> (64 - 2*i)) & 3;
		contig->next = NULL;
		contig->ridposnum = 0;
		contig->pridposorientation = NULL;

		if(i==1)
		{
			*contighead = contig;
			pre_contig = contig;
		}else
		{
			pre_contig->next = contig;
			pre_contig = contig;
		}
	}
	*tailContig = pre_contig;

	return SUCCESSFUL;
}

/**
 * set first kmer threshold.
 */
short initFirstKmerThreshold()
{
	int i, filledBucketNum;
	uint64_t tmp_kmerSum, tmp_kmerNum;
	kmertype *kmer;

	filledBucketNum = 0;
	tmp_kmerSum = tmp_kmerNum = 0;
	for(i=0; i<hashTableSize; i++)
	{
		kmer = deBruijnGraph->pkmers[i];
		if(kmer!=NULL)
			filledBucketNum ++;

		while(kmer)
		{
			tmp_kmerSum += kmer->arraysize;
			tmp_kmerNum ++;

			kmer = kmer->next;
		}
	}

	averKmerOcc = (double)tmp_kmerSum / tmp_kmerNum;
	if(averKmerOcc<1)
		averKmerOcc = 1;

	firstKmerThres = averKmerOcc;

	//printf("filledBucketNum=%d, filledRatio=%.2f, kmerNum=%lu\n", filledBucketNum, (double)filledBucketNum/hashTableSize, tmp_kmerNum);

	//printf("averKmerOcc=%.2f.\n", averKmerOcc);
	//printf("firstKmerThres=%.2f.\n", firstKmerThres);

	return SUCCESSFUL;
}


/**
 * 取得起始的kmers.
 */
short getFirstKmers(uint64_t *kmerIndex, kmertype **firstKmer)
{
	kmertype *kmer;
	uint64_t i;

	i = *kmerIndex;
	kmer = *firstKmer;
	kmers[0] = kmers[1] = NULL;
	while(i<hashTableSize)
	{
		if(kmer==NULL)
			kmer = deBruijnGraph->pkmers[i];
		else
			kmer = kmer->next;

		while(kmer)
		{
			//########################## Debug information #############################
//			if(kmer->multiplicity>=firstKmerThres && kmer->multiplicity<=FIRSTKMER_FACTOR*firstKmerThres)
//				printf("########### i=%lu, multiplicity=%u\n", i, kmer->multiplicity);
			//########################## Debug information #############################

			if(isValidFirstKmer(kmer)==YES)
			{
				//########################## Debug information #############################
//				printf("########### i=%lu\n",i);
				//########################## Debug information #############################

				//if(kmer->multiplicity>=firstKmerThres && kmer->multiplicity>=FIRSTKMER_SUBTRACT_THRESHOLD*kmer->arraysize)
				if(kmer->multiplicity>=firstKmerThres && kmer->multiplicity<=FIRSTKMER_FACTOR*firstKmerThres && kmer->multiplicity>=FIRSTKMER_SUBTRACT_THRESHOLD*kmer->arraysize)
				{
					if(containFirstPos(kmer)==YES)
					{
						kmers[0] = kmer;

						if(memcpy(kmerSeqIntAssembly, kmer->kmerseq, entriesPerKmer*sizeof(uint64_t))==NULL)
						{
							printf("line=%d, In %s(), can not memcpy the kmerseq. error!\n", __LINE__, __func__);
							return FAILED;
						}
						kmers[1] = getReverseKmer(kmerSeqIntAssemblyRev, kmerSeqIntAssembly, deBruijnGraph);

						*kmerIndex = i;
						*firstKmer = kmer;

						return SUCCESSFUL;
					}
					else if(containLastPos(kmer)==YES)
					{
						kmers[1] = kmer;

						if(memcpy(kmerSeqIntAssemblyRev, kmer->kmerseq, entriesPerKmer*sizeof(uint64_t))==NULL)
						{
							printf("line=%d, In %s(), can not memcpy the kmerseq. error!\n", __LINE__, __func__);
							return FAILED;
						}
						kmers[0] = getReverseKmer(kmerSeqIntAssembly, kmerSeqIntAssemblyRev, deBruijnGraph);

						*kmerIndex = i;
						*firstKmer = kmer;

						return SUCCESSFUL;
					}
				}
			}
			kmer = kmer->next;
		}
		i++;
	}
	if(i>=hashTableSize)
	{
		*kmerIndex = i;
		*firstKmer = NULL;
	}

	return SUCCESSFUL;
}

/**
 * whether a kmer is a valid first kmer.
 *  @return:
 *  	If valid, return YES; otherwise return NO.
 *  	If errors, return ERROR.
 */
short isValidFirstKmer(kmertype *kmer)
{
	int i, j, tabooFlag;
	uint64_t tmpTabooSeq;
	uint64_t *seqInt;

	if(kmer==NULL)
		return NO;

	seqInt = kmer->kmerseq;
	if(seqInt==NULL)
	{
		printf("line=%d, In %s(), can not get the seqInt. Error.\n", __LINE__, __func__);
		return ERROR;
	}

	tabooFlag = NO;
	if(entriesPerKmer >= 2)
	{
		for(j=0; j<4; j++)
		{
			if(seqInt[0] == tabooSeqInt[j])
			{
				tabooFlag = YES;
				tmpTabooSeq = seqInt[0];
				break;
			}
		}

		if(tabooFlag==YES)
		{
			for(i=1; i<entriesPerKmer-1; i++)
			{
				if(seqInt[i]!=tmpTabooSeq)
					return YES;
			}

			if(seqInt[entriesPerKmer-1] == (tmpTabooSeq & lastEntryMask))
				return NO;
		}else
			return YES;
	}else
	{
		for(j=0; j<4; j++)
		{
			if(seqInt[0] == (tabooSeqInt[j] & lastEntryMask))
				return NO;
		}
	}

	return YES;
}

/**
 * 检查kmer是否包含在read的第一个位置，如果是，返回YES，否则返回NO。
 */
short containFirstPos(kmertype *kmer)
{
	int i, posNum;
	ridpostype *ridpostable;

	if(kmer->multiplicity==0)
		return NO;

	posNum = kmer->arraysize;
	ridpostable = kmer->ppos; //取得ridpos表格
	for(i=0; i<posNum; i++)
	{
		if(ridpostable[i].delsign==0 && ridpostable[i].pos==1)
			return YES;
	}
	return NO;
}

/**
 * 检查kmer是否包含在read的第一个位置，如果是，返回YES，否则返回NO。
 */
short containLastPos(kmertype *kmer)
{
	int i, posNum;
	ridpostype *ridpostable;  //取得ridpos表格

	if(kmer->multiplicity==0)
		return NO;

	posNum = kmer->arraysize;
	ridpostable = kmer->ppos;  //取得ridpos表格
	for(i=0; i<posNum; i++)
	{
		if(ridpostable[i].delsign==0 && ridpostable[i].pos==readLen-kmerSize+1)
			return YES;
	}
	return NO;
}

/**
 * Add first kmer of a contig to decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int addFirstKmerToDecisionTable(kmertype **kmers)
{
	ridpostype *ridpostable;
	unsigned int i, rpos, posNum;

	//处理正向的kmer
	if(kmers[0])
	{ //正向kmer不为空, 则添加pos==1的reads到决策表
		ridpostable = kmers[0]->ppos;
		posNum = kmers[0]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;
			if(rpos==1 && ridpostable[i].delsign==0) //该位置为起始位置，并且未被删除
			{
				// ############################ Debug information ##############################
				//if(ridpostable[i].rid==664454)
				//{
				//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
				//}
				// ############################ Debug information ##############################

				if(addReadToDecisionTable(ridpostable[i].rid, rpos, ORIENTATION_PLUS, NO)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read %lu, to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
					return FAILED;
				}
			}
		}
	}

	//处理反向互补的kmer
	if(kmers[1])
	{ //反向互补kmer不为空, 则添加pos==25的reads到决策表
		ridpostable = kmers[1]->ppos;
		posNum = kmers[1]->arraysize;
		for(i=0; i<posNum; i++)
		{
			rpos = ridpostable[i].pos;
			if(rpos>readLen-kmerSize+1-errorRegLenEnd3 && ridpostable[i].delsign==0) //该位置为起始位置，并且未被删除
			//if(rpos==READ_LEN-KMER_SIZE+1 && ridpostable[i].delsign==0) //该位置为起始位置，并且未被删除
			{
				// ############################ Debug information ##############################
				//if(ridpostable[i].rid==664454)
				//{
				//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
				//}
				// ############################ Debug information ##############################

				if(addReadToDecisionTable(ridpostable[i].rid, rpos, ORIENTATION_MINUS, NO)==FAILED)
				{
					printf("line=%d, In %s(), cannot add read %lu, to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
					return FAILED;
				}
			}
		}
	}

	return SUCCESSFUL;
}

/**
 * Add a read into decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int addReadToDecisionTable(uint64_t rid, int rpos, int orientation, int matedFlag)
{
	assemblingreadtype *this_assemblingRead = decisionTable + itemNumDecisionTable;

	// ####################### Debug information ##########################
	//if(rid==8399857)
	//{
	//	printf("line=%d, In %s(), rid=%lu, rpos=%d, orient=%d, matedFlag=%d\n", __LINE__, __func__, rid, rpos, orientation, matedFlag);
	//}
	// ####################### Debug information ##########################

	this_assemblingRead->rid = rid;
	this_assemblingRead->firstpos = rpos;
	this_assemblingRead->orientation = orientation;
	this_assemblingRead->status = ASSEMBLING_STATUS;
	this_assemblingRead->kmerappeartimes = 1;
	this_assemblingRead->kmerunappeartimes = 0;
	this_assemblingRead->lastappearpos = rpos;
	this_assemblingRead->lastpos = rpos;
	this_assemblingRead->kmerunappearblocks = 0;
	this_assemblingRead->delsign = 0;
	this_assemblingRead->reserved = 0;
	this_assemblingRead->locked = 0;
	this_assemblingRead->matedFlag = matedFlag;

	itemNumDecisionTable ++;

	// check the allowed number of reads in the array
	if(itemNumDecisionTable==maxItemNumDecisionTable)
	{
		if(reallocateDecisionTable()==FAILED)
		{
			printf("line=%d, In %s(), cannot reallocate memory for the reads in decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

    return SUCCESSFUL;
}

/**
 * Reallocate the memory of decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int reallocateDecisionTable()
{
	assemblingreadtype *tmp_pDecisionTable;
	tmp_pDecisionTable = (assemblingreadtype*) malloc (2*maxItemNumDecisionTable * sizeof(assemblingreadtype));
	if(tmp_pDecisionTable==NULL)
	{
		printf("In %s(), cannot reallocate memory for the reads in decision table, error!\n", __func__);
		return FAILED;
	}

	// copy memory
	if(memcpy(tmp_pDecisionTable, decisionTable, maxItemNumDecisionTable * sizeof(assemblingreadtype))==NULL)
	{
		printf("In %s(), cannot copy memory for the reads in decision table, error!\n", __func__);
		return FAILED;
	}

	free(decisionTable);
	decisionTable = tmp_pDecisionTable;
	maxItemNumDecisionTable *= 2;

	return SUCCESSFUL;
}

/**
*    取下次拼接的kmer.
*    策略, 两个过程:
*        (1) 首先, 以决策表为基础, 计算得分. 该过程只考虑决策表中存在的reads.
*        (2) 其次, 计算决策表中不存在的reads的得分. 也就是新的reads的得分.
*/
short getNextKmerBySE(int contigNodesNum)
{
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	kmertype *tmp_kmers[4][2] = {{0,0},{0,0},{0,0},{0,0}};  //tmp_kmers[i][0]为正向的kmer, tmp_kmers[i][1]为反向的kmer
	short validKmerNum = 0, base_index = 0; //有效的kmer数目
	int i, j;
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
		occsNumSE[i] = 0;

		//开始计算
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (kmerSeqIntAssembly[j] << 2) | (kmerSeqIntAssembly[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (kmerSeqIntAssembly[entriesPerKmer-2] << 2) | (kmerSeqIntAssembly[entriesPerKmer-1] >> (2*(lastEntryBaseNum-1)));
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
		//score[base_index] = computeKmerScore(tmp_kmers[base_index], occsNum+base_index, assemblingreads, numassemblingreads);
		//score[base_index] = computeKmerScoreUnlocked(tmp_kmers[base_index], occsNum+base_index, assemblingreads, numassemblingreads);
		if(computeKmerOccNumUnlocked(tmp_kmers[base_index], occsNumSE+base_index)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum by all reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
		//========================= condition 1 ==============================
		if(occsNumSE[base_index]>0)
		//if(occsNumSE[base_index]>=minKmerOccSE)
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

			maxOccIndexSE = base_index;
			maxOccSE = occsNumSE[maxOccIndexSE];
			secondOccIndexSE = -1;
			secondOccSE = 0;
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

	//kmer_len = longKmerSize;

	if(contigNodesNum>=longKmerSize)
		kmer_len = longKmerSize;
	else
		kmer_len = contigNodesNum;

	//while(kmer_len>kmerSize)
	while(kmer_len>=kmerSize)
	//while(kmer_len>=MIN_KMER_SIZE)
	{
		//开始计算每个kmer得分
		validKmerNum = 0;
		maxOccSE = 0, secondOccSE = 0;
		maxOccIndexSE = -1, secondOccIndexSE = -1;
		for(i=0; i<4; i++)
		{
			occsNumSE[i] = 0;

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
					if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<minKmerOccSE)
					//if(tmp_kmers[i][0]->multiplicity+tmp_kmers[i][1]->multiplicity<=0)
					{
						continue;
					}
				}else if(tmp_kmers[i][0]!=NULL)
				{
					if(tmp_kmers[i][0]->multiplicity<minKmerOccSE)
					//if(tmp_kmers[i][0]->multiplicity<=0)
					{
						continue;
					}
				}else if(tmp_kmers[i][1]!=NULL)
				{
					if(tmp_kmers[i][1]->multiplicity<minKmerOccSE)
					//if(tmp_kmers[i][1]->multiplicity<=0)
					{
						continue;
					}
				}
			}

			//########################## end ########################//

			//if(contigNodesNum>READ_LEN)
			if(contigNodesNum>kmer_len)
			//if(contigNodesNum>=kmer_len)  //--bad result
			{
				if(computeLongKmerOccNum(tmp_kmers[i], occsNumSE+i, kmer_len)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute kmer occNum by long kmer, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}else
			{
				if(computeKmerOccNum(tmp_kmers[i], occsNumSE+i)==FAILED)
				{
					printf("line=%d, In %s(), cannot compute kmer occNum, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}

			//score[i] = computeKmerScore(tmp_kmers[i], occsNum+i, assemblingreads, numassemblingreads);

			//========================= condition 3 ==============================
			//if(score[i]>=MIN_SCORE_THRESHOLD/**MIN_SCORE_FACTOR*/)
			//if(score[i]>=MIN_SCORE_THRESHOLD && occsNum[i]>=MIN_CONNECT_KMER_NUM)
			//if(occsNumSE[i]>=MIN_CONNECT_KMER_NUM)
			if(occsNumSE[i]>0)
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
			maxOccSE = 0, maxOccIndexSE = -1, secondOccSE = 0, secondOccIndexSE = -1;
			for(j=0; j<4; j++)
			{
				if(maxOccSE<occsNumSE[j])
				{
					secondOccSE = maxOccSE;
					secondOccIndexSE = maxOccIndexSE;
					maxOccSE = occsNumSE[j];
					maxOccIndexSE = j;
				}else if(secondOccSE<occsNumSE[j])
				{
					secondOccSE = occsNumSE[j];
					secondOccIndexSE = j;
				}
			}
		}

/*
		//========================= condition 4 ==============================
		//if(validKmerNum>1 && maxIndex1!=maxOccIndex)
		//if(validKmerNum>1 && occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM)
		//if(validKmerNum>0 && occsNumSE[maxIndex1]<minKmerOccSE) //--best result
		//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || (validKmerNum>1 && maxIndex1!=maxOccIndex))
		//if(occsNum[maxIndex1]<MIN_CONNECT_KMER_NUM || occsNum[secondOccIndex]>=MIN_CONNECT_KMER_NUM)
		if(validKmerNum>0 && maxOcc<minKmerOccSE)  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		{
			validKmerNum = 0;
		}
*/

/*
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
		if(validKmerNum>1 && (occsNumSE[secondIndex1] > maxSecondOcc || (occsNumSE[maxIndex1] > maxFirstOcc && occsNumSE[secondIndex1]>=maxSecondOcc)
				|| ((float)occsNumSE[secondIndex1]/occsNumSE[maxIndex1]>SECOND_FIRST_OCC_RATIO && occsNumSE[secondIndex1]>=maxSecondOcc))) //-- best result
		{
			validKmerNum = 0;
		}
*/

/*
		//========================= condition 6 ==============================
		//if(validKmerNum>1 && ((float)occsNumSE[maxIndex1]/itemNumDecisionTable < VALID_OCC_RATIO && occsNumSE[secondIndex1] >= maxSecondOcc))
		if(validKmerNum>1 && ((float)maxOcc/itemNumDecisionTable < VALID_OCC_RATIO && secondOcc > maxSecondOcc))
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
		//if(validKmerNum>1 && ((secondOccSE/maxOccSE>SECOND_FIRST_OCC_RATIO && secondOccSE>minKmerOccSE) || (secondOccSE>maxSecondOcc)))
		if(maxOccSE<minKmerOccSE || (validKmerNum>1 && ((secondOccSE/maxOccSE>SECOND_FIRST_OCC_RATIO && secondOccSE>minKmerOccSE) || (secondOccSE>maxSecondOcc))))
		//if(maxOccSE<minKmerOccSE || (validKmerNum>1 && (secondOccSE/maxOccSE>SECOND_FIRST_OCC_RATIO || (secondOccSE>maxSecondOcc))))
		{
			if(secondOccSE/maxOccSE>SECOND_FIRST_OCC_FAILED_RATIO)
				validKmerNum = 0;
		}


		//========================= condition 9 ==============================
		//if(kmer_len > longKmerSize - longKmerStepSize && occsNumSE[maxIndex1] <= minLongKmerOcc)
		if(kmer_len > longKmerSize - longKmerStepSize && maxOccSE < minLongKmerOcc)
		{
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
				kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexSE) & lastEntryMask;

				kmers[0] = tmp_kmers[maxOccIndexSE][0];
				kmers[1] = tmp_kmers[maxOccIndexSE][1];
//			}

			return SUCCESSFUL;

		}
		//***********************************************************************************
		else if(validKmerNum==0)
		{
//			if(length_k<=KMER_SIZE+KMER_SIZE_STEP)
//				length_k -= KMER_SIZE_STEP / 2;
//			else
//				length_k -= KMER_SIZE_STEP;

			kmer_len -= longKmerStepSize;
			continue;
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
			kmerSeqIntAssembly[entriesPerKmer-1] = ((kmerSeqIntAssembly[entriesPerKmer-1] << 2) | maxOccIndexSE) & lastEntryMask;

			kmers[0] = tmp_kmers[maxOccIndexSE][0];
			kmers[1] = tmp_kmers[maxOccIndexSE][1];

			return SUCCESSFUL;
		}
		//***********************************************************************************
	}
	kmers[0] = kmers[1] = NULL;
	return SUCCESSFUL;
}


/**
 * Compute the occNum according to long K-mers.
 * 	@return:
 * 		If succeeds, return SUCCSEEFUL; otherwise, return FAILED.
 */
int computeLongKmerOccNum(kmertype *tmp_kmers[2], int *occNum, int length_k)
{
	assemblingreadtype *this_assemblingRead = NULL;
	ridpostype *rid_pos = NULL;
	ridpostype *rid_pos_table = NULL;
	int posNum = 0;
	unsigned int this_rid = 0;
	unsigned int this_pos = 0;
	int i = 0, j = 0, properIndex = -1, startIndex = -1, exectIndex = -1, limitLastpos = -1, exitFlag = NO;

	*occNum = 0;

	for(i=0; i<itemNumDecisionTable; i++)
	{ //顺序访问决策表中的每行
		//找合适的lastpos的read
		if((decisionTable[i].reserved==0 && decisionTable[i].lastpos>0)
				&& (decisionTable[i].kmerunappeartimes==0 && decisionTable[i].kmerunappearblocks==0)
				&& ((decisionTable[i].orientation==ORIENTATION_PLUS && decisionTable[i].lastpos>length_k-kmerSize)
						||(decisionTable[i].orientation==ORIENTATION_MINUS && decisionTable[i].lastpos<readLen-kmerSize-(length_k-kmerSize)+1)))
		{
			properIndex = getProperIndex(decisionTable+i, decisionTable, itemNumDecisionTable);

			// ############################ Debug information ##############################
			if(properIndex<0)
			{
				printf("line=%d, In %s(), properIndex=%d, i=%d, Error!\n", __LINE__, __func__, properIndex, i);
				continue;
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
					(*occNum) ++;
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
					(*occNum) ++;
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
 * Compute kmer occurrence number.
 *  If there are some locked reads, the occNum will be computed by the locked reads;
 *  otherwise, consider all reads.
 *
 *  @return:
 *   	If succeeds, return SUCCSEEFUL; otherwise, return FAILED.
*/
int computeKmerOccNum(kmertype *tmp_kmers[2], int *occNum)
{
	//if(lockedReadsNum>=KOCKED_READS_NUM_THRESHOLD)
	if(lockedReadsNum>=lockedReadsNumThres)
	{ //如果决策表中含有锁定的reads, 则只考虑锁定的reads.
		if(computeKmerOccNumLocked(tmp_kmers, occNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum by locked reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}else
	{ //否则,考虑全部的reads.
		//计算决策表中的reads的occNum
		if(computeKmerOccNumUnlocked(tmp_kmers, occNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot compute kmer occNum by all reads, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * Compute kmer occNum by all reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int computeKmerOccNumUnlocked(kmertype *tmp_kmers[2], int *occNum)
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
		if(decisionTable[i].reserved==0 && decisionTable[i].lastpos>0)
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
								rid_pos = rid_pos_table + exectIndex;
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
	if(tmp_kmers[0])
	{ //正向kmer不为空
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				(*occNum) ++;
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
				(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0; //决策表中保留位清零

	return SUCCESSFUL;
}


/**
 * Compute the kmer occurrence number by locked reads.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int computeKmerOccNumLocked(kmertype *tmp_kmers[2], int *occNum)
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
			if(decisionTable[i].reserved==0 && decisionTable[i].lastpos>0)
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
	if(tmp_kmers[0])
	{ //正向kmer不为空
		rid_pos_table = tmp_kmers[0]->ppos;
		posNum = tmp_kmers[0]->arraysize;
		//计算ridpostable表中未考虑的reads的得分, 也即是新的reads的得分: 只要可以拼接，则按照1分计算
		for(i=0; i<posNum; i++)//遍历ridpostable表
		{
			if(rid_pos_table[i].delsign==0 && rid_pos_table[i].pos==1)
			{  //未被删除, 未标记并且位置合理,则将该read可以拼接, 并计算该read的得分
				(*occNum) ++;
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
				(*occNum) ++;
			}
			//rid_pos_table[i].reserved = 0; //ridpostable表中保留位清零
		}
	}
#endif

	for(i=0; i<itemNumDecisionTable; i++) decisionTable[i].reserved = 0; //决策表中保留位清零

	return SUCCESSFUL;
}

/**
 * Append a base to contig tail.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
*/
int addContigBase(contigtype **contigtail, unsigned char base, int contigIndex)
{
	contigtype *contignode = (contigtype*) malloc(sizeof(contigtype));
	if(contignode==NULL)
	{  //分配内存失败，返回NULL
		printf("line=%d, In %s(), can not allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}else
	{
		contignode->index = contigIndex;
		contignode->base = base;
		contignode->ridposnum = 0;
		contignode->pridposorientation = NULL;
		contignode->next = NULL;
		contignode->reserved = 0;
	}

	(*contigtail)->next = contignode;
	*contigtail = contignode;

	return SUCCESSFUL;
}

/**
    在contig中，追加拼接成功结束的reads的信息.
		@return:
			成功, 返回成功标记; 失败, 返回失败标记.
*/
short addRidposToContig(successRead_t *successReadArray, int *successReadNum, contigtype *contig36, int contigNodesNum)
{
	contigtype *contig;
	successRead_t *ridposorient = NULL;
	int i, j, num, k, Num, successNum, sumNewReads, oldReadsNum, tmp_contigNodesNum, //成功的reads数量, 当次添加的新的reads的总数量, 该节点上的原有的reads数量
			checkContigNum; //限制检测contig的数量
	short indicator[*successReadNum];

	//如果contig节点个数大于36, 则contig36所指向的节点个数为36
	tmp_contigNodesNum = contigNodesNum;
	if(tmp_contigNodesNum>readLen)
		tmp_contigNodesNum = readLen;

	Num = readLen - kmerSize + 1;
	i = 0;
	successNum = 0;
	sumNewReads = 0;
	checkContigNum = 0;
	contig = contig36;
	while(contig)
	{
		i++;
		if(assemblyRound==FIRST_ROUND_ASSEMBLY && i<tmp_contigNodesNum-1)
		{ //在第一轮拼接时, 只检测最后的contig
			contig = contig->next; //指向下一个contig
			continue;
		}

		if(assemblyRound==SECOND_ROUND_ASSEMBLY)
		{ //在第二轮拼接时, 需要限定检测的数量
			checkContigNum ++;
			if(checkContigNum>errorRegLenEnd3+1)  //限制检测数量为5
				break;
		}

		if(memset(indicator, 0, (*successReadNum)*sizeof(short))==NULL)
		{
			printf("line=%d, In %s(), cannot memset indicator, error!\n", __LINE__, __func__);
			return FAILED;
		}

		num = 0;
		j = 0;
		oldReadsNum = 0; //设该节点上的原有的reads数量为0
		while(j<(*successReadNum))
		{ //查找该contig中的read个数

			// ############################ Debug information ##############################
			//if(successReadArray[j].rid==4882354)
			//{
			//	printf("line=%d, In %s(), rid=%d, pos=%d, matchnum=%d, orientation=%c, i=%d\n", __LINE__, __func__, successReadArray[j].rid, successReadArray[j].pos, successReadArray[j].matchnum, successReadArray[j].orientation, i);
			//}
			// ############################ Debug information ##############################

			if(assemblyRound==FIRST_ROUND_ASSEMBLY)
			{ //第一轮拼接的判断
				//if((ridposorientation[j].orientation==ORIENTATION_PLUS && ridposorientation[j].matchnum+KMER_SIZE+ridposorientation[j].pos-2==i)
						//||(ridposorientation[j].orientation==ORIENTATION_MINUS && contigNodesNum+ridposorientation[j].matchnum-ridposorientation[j].pos==i))
				if((i==tmp_contigNodesNum-1 && successReadArray[j].orientation==ORIENTATION_PLUS && successReadArray[j].matchnum<Num)
						|| (i==tmp_contigNodesNum && ((successReadArray[j].orientation==ORIENTATION_PLUS  && successReadArray[j].matchnum==Num) || successReadArray[j].orientation==ORIENTATION_MINUS)))
				{
					num ++;
					indicator[j] = 1;
				}
			}else
			{ //第二轮拼接的判断
				//if((ridposorientation[j].orientation==ORIENTATION_PLUS && ridposorientation[j].pos==i)
						//||(ridposorientation[j].orientation==ORIENTATION_MINUS && contigNodesNum-KMER_SIZE+2-ridposorientation[j].pos==i))
				//if((i==1 && ridposorientation[j].orientation==ORIENTATION_PLUS)
						//|| (i==Num-ridposorientation[j].matchnum+1 && ridposorientation[j].orientation==ORIENTATION_MINUS)) --2010年8月22日23点注释掉
				if((i==1 && successReadArray[j].matchnum==Num)
						|| (i==tmp_contigNodesNum-kmerSize+1-successReadArray[j].matchnum && successReadArray[j].orientation==ORIENTATION_PLUS)
						|| (i==tmp_contigNodesNum-kmerSize+1-successReadArray[j].matchnum+1 && successReadArray[j].orientation==ORIENTATION_MINUS))
				{
					num ++;
					indicator[j] = 1;
				}
			}
			j++;
		}

		successNum += oldReadsNum; //将原有的reads算作成功的数量

		// ############################ Debug information ##############################
		if(successNum>(*successReadNum))
		{
			printf("line=%d, In %s(), successNum(%d) > successReadNum(%d), Error!\n", __LINE__, __func__, successNum, *successReadNum);
			return FAILED;
		}
		// ############################ Debug information ##############################

		if(num>0)
		{ //如果read个数大于0，则添加这些reads的信息
			//如果该contig节点已经有成功的reads, 则需要将这些已经成功的reads合并到新的reads集合中
			ridposorient = (successRead_t*) malloc((contig->ridposnum + num)*sizeof(successRead_t));
			if(ridposorient==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}

			if(contig->ridposnum>0)
			{ //该contig已经有成功的reads, 则将其复制到新开辟的数组中
				//复制原来的成功的reads
				if(memcpy(ridposorient, contig->pridposorientation, contig->ridposnum*sizeof(successRead_t))==NULL)
				{
					printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				//释放掉原来的位置数组占用的内存
				free(contig->pridposorientation);
				contig->pridposorientation = NULL;
			}

			//将开辟出来的ridposorientation结构添加到contig中
			contig->pridposorientation = ridposorient;

			j = 0;
			k = contig->ridposnum;
			while(j<(*successReadNum)) //添加新的成功的reads
			{ //循环遍历成功的reads表，将以该contig开始的reads添加进该contig
				if(indicator[j]==1)
				{ //该read需要添加到该contig节点
					if(assemblyRound==FIRST_ROUND_ASSEMBLY)
					{ //第一轮拼接需要变换成功reads的pos和方向
						if(successReadArray[j].orientation==ORIENTATION_PLUS)
						{ //标示为正向, 需要转化为反向
							//按照read中的碱基位置编号, 由条件1得到
							ridposorient[k].pos = successReadArray[j].matchnum + kmerSize + successReadArray[j].pos - 2;
							//按照read中的kmer位置编号
							//ridposorient[k].pos = ridposorientation[j].matchnum+ridposorientation[j].pos-1;
							ridposorient[k].orientation = ORIENTATION_MINUS;
						}else
						{ //标示为反向, 需要转化为正向
							//按照read中的碱基位置编号, 由条件2得到: READ_LEN+1-(matchnum+READ_LEN-pos)
							ridposorient[k].pos = successReadArray[j].pos - successReadArray[j].matchnum + 1;
							//按照read中的kmer位置编号, 与按照碱基位置编号相同
							//ridposorient[k].pos = ridposorientation[j].pos - ridposorientation[j].matchnum + 1;
							ridposorient[k].orientation = ORIENTATION_PLUS;
						}
					}else
					{ //第二轮方向不需要改变
						if(successReadArray[j].orientation==ORIENTATION_PLUS)
						{ //正向的reads, 位置直接赋值
							ridposorient[k].pos = successReadArray[j].pos;
						}else
						{ //反向的reads, 位置需要变换
							ridposorient[k].pos = successReadArray[j].pos + kmerSize - 1;
						}
						ridposorient[k].orientation = successReadArray[j].orientation;
					}
					ridposorient[k].rid = successReadArray[j].rid;
					//ridposorient[k].matchnum = ridposorientation[j].matchnum; //按照read中的kmer位置计算
					ridposorient[k].matchnum = successReadArray[j].matchnum + kmerSize - 1; //按照read中的碱基位置计算

					//ridposorient[k].deledlink = successReadArray[j].deledlink;

					//统计新的reads数量
					//if(successReadArray[j].deledlink==NO)
					//{
						sumNewReads ++;
					//}

					successNum ++;

					k++;
				}
				j++;
			}

			contig->ridposnum = k; //更新成功的reads数量

			if(successNum==(*successReadNum))
			{ //所有的reads都已经添加完毕, 则退出循环
				break;
			}

		} //end if(num>0)

		contig = contig->next; //指向下一个contig
	} //end while(contig)

	*successReadNum = sumNewReads; //新有的reads数量

	return SUCCESSFUL;
}

/**
 * Release contig nodes.
 */
void releaseContig(contigtype *contighead)
{
	contigtype *contig = contighead;
	while(contig)
	{
		contig = contighead->next;

		if(contighead->ridposnum>0)
			free(contighead->pridposorientation);
		free(contighead);

		contighead = contig;
	}
	contighead = NULL;
}

/**
 * Output the contig nodes to file.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short outputContigToFile(FILE *fpContig, int outFileType, contigtype *contighead, int contigID, int contigNodeNum)
{
	contigtype *contig;
	successRead_t *ridposorientation;
	int i, num, index;
	char tmp_base;

	if(fpContig)
	{
		if(outFileType==BASE_TYPE_FASTA_CONTIG_FILE)
		{
			fprintf(fpContig, ">%d length: %d\n", contigID, contigNodeNum); //输出标题信息
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

				fprintf(fpContig, "%c", tmp_base); //输出碱基
				contig = contig->next;
				//contigID++;
			}
			fprintf(fpContig, "\n");
			//fflush(pFileResult);
		}else if(outFileType==HANGING_READ_TYPE_CONTIG_FILE)
		{ // output the contigs in hanging reads format
			fprintf(fpContig, ">%d length: %d\n", contigID, contigNodeNum); //输出标题信息

			index = 1;
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

				fprintf(fpContig, "%d\t%c\t%d", index++, tmp_base, contig->ridposnum); //输出碱基
				ridposorientation = contig->pridposorientation;
				num = contig->ridposnum;
				for(i=0; i<num; i++)
				{
					fprintf(fpContig, "\t(%lu,%d,%d,%c)",
							ridposorientation[i].rid,
							ridposorientation[i].pos,
							ridposorientation[i].matchnum,
							ridposorientation[i].orientation);
				}
				fprintf(fpContig, "\n");

				contig = contig->next;
				//contigID++;
			}
			//fflush(pFileResult);
		}else
		{
			printf("line=%d, In %s(), unknown output contig file type %d, error!\n", __LINE__, __func__, outFileType);
			return FAILED;
		}
	}else
	{
		printf("line=%d, In %s(), file pointer is NULL, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Delete the successful reads from de Bruijn graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delReadsFromGraph(successRead_t *successReadArray, int successReadNum, char *lastseq36)
{
	int i, j, rpos, len;
	char /*tmpseq[KMER_SIZE+1] = {0},*/ reversed_lastseq36[readLen+1];
	short Num, perfectMatchFlag;  //kmer数目, read完全与contig匹配的标记
	uint64_t hashcode, rid;
	uint64_t tmp_kmerseq[entriesPerKmer];
	int baseInt, first_kmer_index;//记录read在lastseq36中第一个kmer的下标

	len = strlen(lastseq36);

	//将lastseq36复制到临时数组
	strcpy(reversed_lastseq36, lastseq36);
	//reversed_lastseq36[len] = '\0';

	//将临时碱基反向互补
	if(reverseReadseq(reversed_lastseq36)==FAILED)
	{
		printf("line=%d, In %s(), cannot reverse the read sequence: %s", __LINE__, __func__, reversed_lastseq36);
		return FAILED;
	}

	//取得正向的序列的第一个kmer的哈希值
	//strncpy(tmpseq, lastseq36, KMER_SIZE); //取得该kmer的碱基序列
	//tmpseq[KMER_SIZE] = '\0';
	//hashcode = kmerhash(tmpseq);
	//first_hashcode = kmerhash(lastseq36);

	//取得反向互补的序列的第一个kmer的哈希值
	//strncpy(tmpseq, reversed_lastseq36, KMER_SIZE); //取得该kmer的碱基序列
	//tmpseq[KMER_SIZE] = '\0';


	i = 0;
	Num = len - kmerSize + 1;
	while(i<successReadNum)
	{
		perfectMatchFlag = YES;
		rid = successReadArray[i].rid; //取得read ID

		// ############################ Debug information ##############################
		//if(rid==8399857)
		//{
		//	printf("line=%d, In %s(), rid=%lu, pos=%d, matchnum=%d, orientation=%c\n", __LINE__, __func__, rid, successReadArray[i].pos, successReadArray[i].matchnum, successReadArray[i].orientation);
		//}
		// ############################ Debug information ##############################

		if(successReadArray[i].orientation==ORIENTATION_PLUS)
		{ //正向read
			//开始删除该read中的kmers
			//hashcode = first_hashcode;
			if(successReadArray[i].matchnum==readLen-kmerSize+1)
			{
				first_kmer_index = Num - successReadArray[i].matchnum;
				if(first_kmer_index<0)
				{
					printf("line=%d, In %s(), first_kmer_index=%d, error!\n", __LINE__, __func__, first_kmer_index);
					return FAILED;
				}

			}else
			{
				first_kmer_index = Num - successReadArray[i].matchnum - 1;
				if(first_kmer_index<0)
				{
					printf("line=%d, In %s(), first_kmer_index=%d, error!\n", __LINE__, __func__, first_kmer_index);
					return FAILED;
				}
			}

			// generate the kmer integer sequence
			if(generateKmerSeqInt(tmp_kmerseq, lastseq36+first_kmer_index)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			rpos = 1;
			while(rpos<=successReadArray[i].matchnum)
			{
				// 删除该read的一个kmer, 也即是该read的一个ridpos位置
				hashcode = kmerhashInt(tmp_kmerseq);
				if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
				{
					perfectMatchFlag = NO; //该read与contig并不是完全匹配
					if(delRemainedKmers(lastseq36+first_kmer_index+rpos-1, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), can not delete the read [%lu, %c] from position %d. Error!\n", __LINE__, __func__, rid, successReadArray[i].orientation, rpos);
						return FAILED;
					}

					break; //该read的删除操作终止
				}

				if(rpos==successReadArray[i].matchnum)
					break;

				//取下一个kmer
				switch(lastseq36[first_kmer_index+kmerSize+rpos-1])
				{
					case 'A':
					case 'a':
						baseInt = 0;
						break;
					case 'C':
					case 'c':
						baseInt = 1;
						break;
					case 'G':
					case 'g':
						baseInt = 2;
						break;
					case 'T':
					case 't':
						baseInt = 3;
						break;
					default:
						printf("line=%d, In %s(), unknown base, Error!\n", __LINE__, __func__);
						return FAILED;
				}

				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
						tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
					tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

				rpos ++;
			}
		}else
		{ //反向read
			// generate the kmer integer sequence
			if(generateKmerSeqInt(tmp_kmerseq, reversed_lastseq36)==FAILED)
			{
				printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
				return FAILED;
			}

			//开始删除该read中的kmers
			//hashcode = reversed_first_hashcode;
			rpos = 1;
			while(rpos<=Num)
			{
				//删除该read的一个kmer, 也即是该read的一个ridpos位置
				hashcode = kmerhashInt(tmp_kmerseq);
				if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
				{
					perfectMatchFlag = NO; //该read与contig并不是完全匹配
					if(delRemainedKmers(reversed_lastseq36+rpos-1, tmp_kmerseq, rid, rpos, deBruijnGraph)==FAILED)
					{
						printf("line=%d, In %s(), can not delete the read [%lu, %c] from position %d. Error!\n", __LINE__, __func__, rid, successReadArray[i].orientation, rpos);
						return FAILED;
					}

					break; //该read的删除操作终止
				}

				if(rpos==Num)
					break;

				//取下一个kmer
				switch(reversed_lastseq36[kmerSize+rpos-1])
				{
					case 'A':
					case 'a':
						baseInt = 0;
						break;
					case 'C':
					case 'c':
						baseInt = 1;
						break;
					case 'G':
					case 'g':
						baseInt = 2;
						break;
					case 'T':
					case 't':
						baseInt = 3;
						break;
					default:
						printf("line=%d, In %s(), unknown base, Error!\n", __LINE__, __func__);
						return FAILED;
				}

				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
						tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
					tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

				rpos ++;
			}
		}

		//在lastseq36长度小于36个碱基的情况下, 将read中的其他kmer删除
		if(rpos<readLen-kmerSize+1 && perfectMatchFlag==YES)
		{
			rpos ++;
			while(rpos<=readLen-kmerSize+1)
			{
				//下一个kmer的碱基的整数表示
				if(entriesPerKmer>=2)
				{
					for(j=0; j<entriesPerKmer-2; j++)
						tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
					tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
				}
				tmp_kmerseq[entriesPerKmer-1] = (tmp_kmerseq[entriesPerKmer-1] << 2) & lastEntryMask;

				for(j=0; j<4; j++)
				{ //循环探测该kmer, 直到找到并将其删除为止
					hashcode = kmerhashInt(tmp_kmerseq);
					if(delKmerByHash(hashcode, tmp_kmerseq, rid, rpos, deBruijnGraph)==SUCCESSFUL)
						break;
					tmp_kmerseq[entriesPerKmer-1] ++;
				}

				// ############################ Debug information ##############################
				if(j==4)
				{
					printf("line=%d, In %s(), cannot delete the kmer (%lu,%d), Error!\n", __LINE__, __func__, rid, rpos);
					return FAILED;
				}
				// ############################ Debug information ##############################

				rpos ++;
			}
		}

		i ++;
	}

	return SUCCESSFUL;
}

/**
 * Reverse the read sequence.
 *  @return:
 *  	If succeeds, return SUCCESFUL; otherwise, return FAILED.
 */
short reverseReadseq(char *str)
{
	if(!str) { printf("line=%d, In %s(): the tuple is NULL! Error!\n", __LINE__, __func__); return FAILED; }
	int i = 0, len = strlen(str);
	char tmp[len+1];
	strcpy(tmp, str);

	while(i<len)
	{
		switch(tmp[i])
		{
			case 'A': str[len-i-1] = 'T'; break;
			case 'C': str[len-i-1] = 'G'; break;
			case 'G': str[len-i-1] = 'C'; break;
			case 'T': str[len-i-1] = 'A'; break;
			case 'a': str[len-i-1] = 't'; break;
			case 'c': str[len-i-1] = 'g'; break;
			case 'g': str[len-i-1] = 'c'; break;
			case 't': str[len-i-1] = 'a'; break;
			default:  str[len-i-1] = tmp[i];  return FAILED;
		}
		i++;
	}
	return SUCCESSFUL;
}

/**
 * Delete the remained kmer from de Bruijn graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short delRemainedKmers(char *seq, uint64_t *tmp_kmerseq, uint64_t rid, uint16_t rpos, graphtype *graph)
{
	int len = strlen(seq);
	int i, j, baseInt;
	uint64_t hash;

	//删除该read的从rpos开始的第1个kmer
	for(j=0; j<3; j++)
	{
		if((tmp_kmerseq[entriesPerKmer-1] & 3)==3)
			tmp_kmerseq[entriesPerKmer-1] -= 3;
		else
			tmp_kmerseq[entriesPerKmer-1] ++;

		hash = kmerhashInt(tmp_kmerseq);
		if(delKmerByHash(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
			break;
	}

	// ############################ Debug information ##############################
	if(j==3)
	{ //该kmer不能删除, 打印出错
		printf("line=%d, In %s(), can not delete this kmer (%lu, %d) Error!\n", __LINE__, __func__, rid, rpos);
		return FAILED;
	}
	// ############################ Debug information ##############################

	//删除该read的从rpos开始的其他kmer
	i = 1;
	rpos++;
	while(i+kmerSize-1<len)
	{
		//取下一个kmer
		switch(seq[i+kmerSize-1])
		{
			case 'A':
			case 'a':
				baseInt = 0;
				break;
			case 'C':
			case 'c':
				baseInt = 1;
				break;
			case 'G':
			case 'g':
				baseInt = 2;
				break;
			case 'T':
			case 't':
				baseInt = 3;
				break;
			default:
				printf("line=%d, In %s(), unknown base, Error!\n", __LINE__, __func__);
				return FAILED;
		}

		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

		//删除该kmer
		hash = kmerhashInt(tmp_kmerseq);
		if(delKmerByHash(hash, tmp_kmerseq, rid, rpos, graph)==FAILED)
		{ //删除该kmer失败, 该kmer需要重新确定
			//重新确定该kmer
			for(j=0; j<3; j++)
			{
				if((tmp_kmerseq[entriesPerKmer-1] & 3)==3)
					tmp_kmerseq[entriesPerKmer-1] -= 3;
				else
					tmp_kmerseq[entriesPerKmer-1] ++;

				hash = kmerhashInt(tmp_kmerseq);
				if(delKmerByHash(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
					break;
			}

			// ############################ Debug information ##############################
			if(j==3)
			{ //该kmer不能删除, 打印出错
				printf("line=%d, In %s(), can not delete this kmer (%lu, %d) Error!\n", __LINE__, __func__, rid, rpos);
				return FAILED;
			}
			// ############################ Debug information ##############################
		}
		i++;
		rpos++;
	}

	//如果存在剩余的kmers, 则将其删除
	while(rpos <= readLen-kmerSize+1)
	{
		//下一个kmer的碱基的整数表示
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = (tmp_kmerseq[entriesPerKmer-1] << 2) & lastEntryMask;

		for(j=0; j<4; j++)
		{ //循环探测该kmer, 直到找到并将其删除为止
			hash = kmerhashInt(tmp_kmerseq);
			if(delKmerByHash(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
				break;
			tmp_kmerseq[entriesPerKmer-1] ++;
		}

		// ############################ Debug information ##############################
		if(j==4)
		{
			printf("line=%d, In %s(), cannot delete the kmer (%lu,%d), Error!\n", __LINE__, __func__, rid, rpos);
			return FAILED;
		}
		// ############################ Debug information ##############################

		rpos++;
	}

	return SUCCESSFUL;
}

/**
 * Update the locked reads and their total number.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int updateLockedReads()
{
	assemblingreadtype *this_assemblingRead;

	int i = 0, appearNum = 0;
	if(lockedReadsNum>0)
	{
		this_assemblingRead = decisionTable;
		//更新锁定的reads数量
		for(i=0; i<itemNumDecisionTable; i++)
		{
			if(this_assemblingRead->locked==1)
			{ //针对锁定的reads
				if(this_assemblingRead->status!=ASSEMBLING_STATUS)
				{ //状态不为拼接状态, 则锁定的数目-1
					lockedReadsNum --;
				}else if(this_assemblingRead->lastpos>0)
				{ //本次拼接出现, 出现次数+1
					appearNum++;
				}
			}
			this_assemblingRead ++;
		}
	}

	//更新锁定标记
	//if(appearNum<minKmerOccSE*2 || lockedReadsNum<KOCKED_READS_NUM_THRESHOLD)
	if(appearNum<minKmerOccSE*2 || lockedReadsNum<lockedReadsNumThres)
	{ //出现的reads数量为0, 或者锁定的reads数量小于阈值, 则重新更新reads的锁定标记

//		if(appearNum>=MIN_CONNECT_KMER_NUM*2)
//		{
//			printf("appearNum=%d, lockedReadsNum=%d\n", appearNum, lockedReadsNum);
//		}else
//		{
//			printf("++++++++++appearNum=%d, lockedReadsNum=%d\n", appearNum, lockedReadsNum);
//		}

		lockedReadsNum = 0; //重新设定初始锁定数量为0

		//更新决策表中的reads的锁定标记
		this_assemblingRead = decisionTable;
		for(i=0; i<itemNumDecisionTable; i++)
		{ //遍历决策表, 寻找可以锁定的reads

			// ############################ Debug information ##############################
			//if(this_assemblingRead->rid==1655993)
			//{
			//	printf("line=%d, In %s(), rid=%d\n", __LINE__, __func__, this_assemblingRead->rid);
			//}
			// ############################ Debug information ##############################

			if(this_assemblingRead->lastpos==0 && this_assemblingRead->locked==1)
			{ //当次拼接未出现, 并且处于锁定状态, 则将该read解锁
				this_assemblingRead->locked = 0;
			}else if(this_assemblingRead->lastpos>0	&& this_assemblingRead->status==ASSEMBLING_STATUS)
			{ //该read当次拼接出现, 并且处于拼接中, 则需要锁定该read, 并更新锁定数量
				this_assemblingRead->locked = 1; //锁定该read
				lockedReadsNum ++; //锁定的reads数量+1
			}
			this_assemblingRead ++;
		}
	}

	return SUCCESSFUL;
}

/**
 * 从决策表中相同rid的read中找lastpos合适的read，返回其下标，并标记rid的这些read.
 *	正向reads取lastpos最大的, 反向reads取lastpos最小的.
 *   @ return:
 *     成功, 返回lastpos合适的read在决策表中的下标; 失败, 返回-1.
 */
int getProperIndex(assemblingreadtype *assemblingread, assemblingreadtype *assemblingreads, int numassemblingreads)
{
	int proLastpos = 0, i, proIndex = -1, rid = assemblingread->rid;

	if(assemblingread->orientation==ORIENTATION_PLUS)
	{ //read正向
		for(i=0; i<numassemblingreads; i++)
		{
			if(assemblingreads[i].rid==rid && assemblingreads[i].orientation==ORIENTATION_PLUS)
			{
				assemblingreads[i].reserved = 1;
				if(assemblingreads[i].lastpos>proLastpos)
				{
					proLastpos = assemblingreads[i].lastpos;
					proIndex = i;
				}
			}
		}
	}else
	{ //read反向
		proLastpos = INT_MAX;
		for(i=0;i<numassemblingreads;i++)
		{
			if(assemblingreads[i].rid==rid && assemblingreads[i].orientation==ORIENTATION_MINUS)
			{
				assemblingreads[i].reserved = 1;
				if(assemblingreads[i].lastpos>0 && assemblingreads[i].lastpos<proLastpos)
				{
					proLastpos = assemblingreads[i].lastpos;
					proIndex = i;
				}
			}
		}
	}

	return proIndex;
}

/**
 * 取得具有限制的合适的位置.
 * 正向reads取lastpos最大的, 反向reads取lastpos最小的.
 *   @ return:
 *     成功, 返回lastpos合适的read的下标; 失败, 返回-1.
 */
int getProperIndexLimited(assemblingreadtype *assemblingread, assemblingreadtype *assemblingreads, int numassemblingreads, int limitLastpos)
{
	int pro, i, proIndex;
	uint64_t rid;

	rid = assemblingread->rid;
	proIndex = -1;
	pro = limitLastpos;

	if(assemblingread->orientation==ORIENTATION_PLUS)
	{ //read正向
		for(i=0; i<numassemblingreads; i++)
		{
			if(assemblingreads[i].rid==rid)
			{
				//assemblingreads[i].reserved = 1;
				if(assemblingreads[i].lastpos>pro)
				{
					pro = assemblingreads[i].lastpos;
					proIndex = i;
				}
			}
		}
	}else
	{ //read反向
		for(i=0;i<numassemblingreads;i++)
		{
			if(assemblingreads[i].rid==rid)
			{
				//assemblingreads[i].reserved = 1;
				if(assemblingreads[i].lastpos>0 && assemblingreads[i].lastpos<pro)
				{
					pro = assemblingreads[i].lastpos;
					proIndex = i;
				}
			}
		}
	}
	return proIndex;
}

/**
 * Trim contig tail nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateContigtailnodes(contigtype *contighead, contigtype *successContig, int *contigIndex)
{
	contigtype *tmp_contig;
	int startIndex;
	int i, j, divisionIndex;
	successRead_t *ridposorient, *tmp_ridposorient;
	int numridposorient, delNum;

	startIndex = successContig->index - readLen + 1;
	if(startIndex<=0)
		startIndex = 1;

	tmp_contig = contighead;
	while(tmp_contig)
	{
		if(tmp_contig->index==startIndex)
			break;

		tmp_contig = tmp_contig->next;
	}

	divisionIndex = successContig->index;
	while(tmp_contig)
	{
		if(tmp_contig->index>divisionIndex)
			break;

		if(tmp_contig->ridposnum>0)
		{
			delNum = 0;
			ridposorient = tmp_contig->pridposorientation;
			numridposorient = tmp_contig->ridposnum;
			short indicator[numridposorient];
			if(memset(indicator, 0, sizeof(short)*numridposorient)==NULL)
			{
				printf("line=%d, In %s(), cannot set indicator to zero, error!\n", __LINE__, __func__);
				return FAILED;
			}

			for(i=0; i<numridposorient; i++)
			{
				if(tmp_contig->index+ridposorient->matchnum-1 > divisionIndex)
				{
					indicator[i] = 1;
					delNum ++;
				}
			}

			if(delNum==numridposorient)
			{
				free(tmp_contig->pridposorientation);
				tmp_contig->pridposorientation = NULL;
				tmp_contig->ridposnum = 0;

			}else if(delNum>0)
			{
				tmp_ridposorient = (successRead_t*) malloc(sizeof(successRead_t)*(numridposorient-delNum));
				if(tmp_ridposorient==NULL)
				{
					printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
					return FAILED;
				}
				tmp_contig->ridposnum = numridposorient - delNum;
				tmp_contig->pridposorientation = tmp_ridposorient;

				j = 0;
				for(i=0; i<numridposorient; i++)
				{
					if(indicator[i]==0)
					{
						if(memcpy(tmp_ridposorient+j, ridposorient+i, sizeof(successRead_t))==NULL)
						{
							printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
							return FAILED;
						}
						j ++;
					}
				}

				free(ridposorient);
			}
		}

		tmp_contig = tmp_contig->next;
	}

	int count = 0;
	while(successContig->next!=NULL)
	{
		tmp_contig = successContig->next;
		successContig->next = tmp_contig->next;
		if(tmp_contig->ridposnum>0)
			free(tmp_contig->pridposorientation);
		free(tmp_contig);
		count ++;
	}
	(*contigIndex) -= count;

	return SUCCESSFUL;
}

/**
 * Trim contig nodes before second round assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short trimContigBeforeCycle2(contigtype **contighead, contigtype **contigtail, contigtype **successContig, int *contigNodeNum)
{
	contigtype *contig = NULL, *firstSuccessContig = NULL;
	firstSuccessContig = getFirstSuccessContig(*contighead);

	if(firstSuccessContig==NULL) //找不到首个成功的contig节点, 返回失败标记
	{
		printf("line=%d, In %s(), cannot get the first successful contig nodes, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//去除头部的contig节点
	if(*contighead != firstSuccessContig)
	{ //contig链表头节点不是成功的节点, 删除头部的不准确的碱基
		while(*contighead)
		{ //循环查找标记成功reads的contig节点, 将其前面的节点删除掉
			if(*contighead==firstSuccessContig)
			{
				break;
			}
			contig = (*contighead)->next;
			if((*contighead)->ridposnum>0)
				free((*contighead)->pridposorientation);
			free(*contighead);
			*contighead = contig;
			(*contigNodeNum) --; //contig节点数量自减
		}
	}

	//去除尾部的contig节点
	if(*contigtail != *successContig)
	{ //尾节点不是成功的节点, 需要删除不准确的碱基

		*contigtail = *successContig; //将contig链表中最后的contig节点赋值为成功的contig节点
		//去除末端的contig节点
		while((*successContig)->next!=NULL)
		{ //循环删除
			contig = (*successContig)->next;
			(*successContig)->next = contig->next;
			if(contig->ridposnum>0)
				free(contig->pridposorientation);
			free(contig);
			(*contigNodeNum) --; //contig节点数量自减
		}

		//将successContig指向contig链表尾节点
		*successContig = *contigtail;
	}

	return SUCCESSFUL;
}

/**
 * Get the first successful contig nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
contigtype *getFirstSuccessContig(contigtype *contighead)
{
	contigtype *contig = contighead;
	successRead_t *ridposorient = NULL;
	short maxMatchNum = 0, successFlag = NO, successDistance = 0; //成功contig节点出现的标记, 从当前节点到第1个成功节点的距离
	int i = 0, numridposorient = 0, contigIndex = 1, successContigIndex = INT_MAX, this_successContigIndex = INT_MAX;
	while(contig)
	{ //遍历contig链表
		if(contig->ridposnum > 0)
		{ //找到标记有成功reads的contig节点

			successFlag = YES;  //更新成功contig节点出现的标记为'YES'

			ridposorient = contig->pridposorientation;
			numridposorient = contig->ridposnum;
			for(i=0; i<numridposorient; i++)
			{
				if(maxMatchNum < ridposorient[i].matchnum)
				{
					maxMatchNum = ridposorient[i].matchnum;
				}
			}

			// ############################ Debug information ##############################
			if(maxMatchNum==0)
			{
				printf("line=%d, In %s(), maxMatchNum==0, error!\n", __LINE__, __func__);
				return NULL;
			}
			// ############################ Debug information ##############################

			this_successContigIndex = contigIndex - maxMatchNum + 1; //按照碱基位置计算
			//this_successContigIndex = contigIndex - maxMatchNum - KMER_SIZE + 2; //按照kmer位置计算

			// ############################ Debug information ##############################
			if(this_successContigIndex<=0)
			{
				printf("In getFirstSuccessContig(), this_successContigIndex<=0, error!\n");
				return NULL;
			}
			// ############################ Debug information ##############################

			if(this_successContigIndex < successContigIndex)
			{ //更新成功的contig的index
				successContigIndex = this_successContigIndex;
			}
		}

		if(successFlag==YES)
		{ //若果已经出现成功的contig, 则最多再检测5个节点
			//if(successDistance>=5)
			if(successDistance>errorRegLenEnd3+1)
			{ //超出距离, 则停止检测
				break;
			}
			successDistance ++;
		}

		contig = contig->next;
		contigIndex ++;
	}

	//查找第1个成功的contig节点, 并将其返回
	contig = contighead;
	i = 1;
	while(contig)
	{
		if(i==successContigIndex)
		{ //找到第一个成功的reads的起始匹配的contig节点
			return contig;
		}
		contig = contig->next;
		i++;
	}

	return NULL;
}

/**
 * 从contig链表中取得contig36.
 * 	@return:
 * 		成功, 返回contig36地址; 失败, 返回NULL.
 */
contigtype *getLastContig(contigtype *contighead, int contigNodesNum, int lastNodesNum)
{
	if(contigNodesNum<=lastNodesNum)
	{ //contig节点数量<=36, 则头节点即为contig36, 直接返回头节点
		return contighead;
	}else
	{
		int i = contigNodesNum;
		contigtype *tmp_contig = contighead;
		//取得第36个contig节点
		while(tmp_contig)
		{
			if(tmp_contig->next==NULL || i==lastNodesNum)
				break;
			tmp_contig = tmp_contig->next;
			i--;
		}
		return tmp_contig;  //将找到的节点即为contig36, 直接返回该节点
	}
}

/**
 * Get the last 36 bp sequence from contig tail.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getLastseq(char *lastseq, contigtype *startContig)
{
	//更新lastseq
	contigtype *tmp_contig;
	int i;

	i = 0;
	tmp_contig = startContig;
	while(tmp_contig)
	{
		switch(tmp_contig->base)
		{
			case 0: lastseq[i] = 'A'; break;
			case 1: lastseq[i] = 'C'; break;
			case 2: lastseq[i] = 'G'; break;
			case 3: lastseq[i] = 'T'; break;
			default:
				printf("line=%d, In %s(), error base %u !\n", __LINE__, __func__, tmp_contig->base);
				return FAILED;
		}
		tmp_contig = tmp_contig->next;
		i++;
	}
	lastseq[i] = '\0';

	// ############################ Debug information ##############################
	if(strlen(lastseq)<kmerSize)
	{
		printf("line=%d, In %s(), lastseq_len==%d\n", __LINE__, __func__, (int)strlen(lastseq));
	}
	// ############################ Debug information ##############################

	return SUCCESSFUL;
}

/**
 * 取得第二轮拼接的kmer.
 *   @ return :
 *   	成功, 返回成功标记; 失败, 返回失败标记.
 */
short getSecondAssemblyFirstKmers(contigtype *contighead, graphtype *graph)
{
	int i, j;
	contigtype *contig;

	for(i=0; i<entriesPerKmer; i++)	kmerSeqIntAssemblyRev[i] = 0;

	contig = contighead;
	i = 0;
	j = 0;
	while(i<kmerSize)
	{
		kmerSeqIntAssemblyRev[j] = (kmerSeqIntAssemblyRev[j] << 2) | contig->base;
		i ++;
		if(i%32==0)
			j++;

		contig = contig->next;
	}

	kmers[1] = getKmer(kmerSeqIntAssemblyRev, graph);
	kmers[0] = getReverseKmer(kmerSeqIntAssembly, kmerSeqIntAssemblyRev, graph);

	//判断是否可以进行第二轮的拼接
	if(kmers[0]==NULL && kmers[1]==NULL) //取得的kmer为空, 则不能进行第二轮的拼接
		return FAILED;

	return SUCCESSFUL;
}

/**
 * Reverse the contig nodes.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short reverseContig(contigtype **contighead, contigtype **contigtail)
{
	//交换头节点和尾节点
	contigtype *tmp, *head, *work;
	int index;

	tmp = *contighead;
	*contighead = *contigtail;
	*contigtail = tmp;

	//将contig的链表反转
	head = (*contigtail)->next;
	work = *contigtail;
	while(work!=(*contighead))
	{
		work->next = (*contighead)->next;
		(*contighead)->next = work;
		tmp = head;
		head = head->next;
		work = tmp;
	}

	//将碱基互补, 并将index反转
	index = 1;
	tmp = *contighead;
	while(tmp)
	{
		switch(tmp->base)
		{
			case 0: tmp->base = 3; break;
			case 1: tmp->base = 2; break;
			case 2: tmp->base = 1; break;
			case 3: tmp->base = 0; break;
			default: printf("line=%d, In %s(), cannot reverse the base, error!\n", __LINE__, __func__); return FAILED;
		}
		tmp->index = index++;
		tmp = tmp->next;
	}

	return SUCCESSFUL;
}

/**
 * 第二轮拼接时, 初始化决策表.
 *   将决策表构造成已经拼接到当前kmers时的状态.
 */
int initAssemblingTableSecondAssembly(char *lastseq36, graphtype *graph)
{
	uint64_t tmp_kmerseq[entriesPerKmer], tmp_kmerseqRev[entriesPerKmer];
	kmertype *tmp_kmers[2];		//tmp_kmers[0]存放正向kmer地址, tmp_kmers[1]存放反向互补kmer地址
	int i, j, lastseqLen, baseInt;

	lastseqLen = strlen(lastseq36);
	// ############################ Debug information ##############################
	if(strlen(lastseq36)<kmerSize)
	{
		printf("line=%d, In %s(), the lastseq36 length==%d, error!\n", __LINE__, __func__, (int)strlen(lastseq36));
		return FAILED;
	}
	// ############################ Debug information ##############################

	if(generateKmerSeqInt(tmp_kmerseq, lastseq36)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence,error!\n", __LINE__, __func__);
		return FAILED;
	}

	tmp_kmers[0] = getKmer(tmp_kmerseq, graph);
	tmp_kmers[1] = getReverseKmer(tmp_kmerseqRev, tmp_kmerseq, graph);

	// initialize the decision table
	itemNumDecisionTable = 0;
	if(addFirstKmerToDecisionTable(tmp_kmers)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	for(i=kmerSize; i<lastseqLen; i++)
	{
		if(PEGivenType>NONE_PE_GIVEN_TYPE && contigIndex-(readLen-i)+1>=minContigLenUsingPE)
		{
			if(updatePEHashTable(contigIndex-(readLen-i)+1, assemblyRound)==FAILED)
			{
				printf("In %s(), cannot update the PE hash table, error!\n", __func__);
				return FAILED;
			}
		}

		//取得正向kmer的碱基序列
		switch(lastseq36[i])
		{
			case 'A': baseInt = 0; break;
			case 'C': baseInt = 1; break;
			case 'G': baseInt = 2; break;
			case 'T': baseInt = 3; break;
			default: printf("line=%d, error base: %c\n", __LINE__, lastseq36[i]); return FAILED;
		}
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

		//取得临时的kmers
		tmp_kmers[0] = getKmer(tmp_kmerseq, graph);
		tmp_kmers[1] = getReverseKmer(tmp_kmerseqRev, tmp_kmerseq, graph);

		// update the decision table according to kmers
		if(updateDecisionTable(tmp_kmers)==FAILED)
		{
			printf("line=%d, In %s(), cannot update decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// Update the reads status in decision table
		if(updateAssemblingreadsStatus()==FAILED)
		{
			printf("line=%d, In %s(), cannot update reads status in decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}

		//将结束的reads从决策表中删除
		//itemNumDecisionTable = resortAssemblingtable(decisionTable, itemNumDecisionTable);
		if(removeFinishedReadsFromDecisionTable()==FAILED)
		{
			printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	return SUCCESSFUL;
}

/**
 * reset the contig index.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; other
 */
short resetContigIndex(contigtype *contighead)
{
	if(contighead==NULL)
		return FAILED;
	if(contighead->index==1)
		return SUCCESSFUL;

	int i = 1;
	contigtype *contig = contighead;
	while(contig)
	{
		contig->index = i;
		contig = contig->next;
		i ++;
	}
	return SUCCESSFUL;
}

/**
 * Get the successful contig.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
contigtype *getSuccessContig(contigtype *contig36, contigtype *successContig, int contigIndex)
{
	contigtype *contig, *new_successContig = NULL;
	successRead_t *ridposorient = NULL;
	int i, numridposorient, minContigIndex, tmp_contigIndex;

	if(assemblyRound==FIRST_ROUND_ASSEMBLY)
	{ //第一轮拼接过程中查找successContig
		minContigIndex = contigIndex;
		contig = contig36;
		while(contig)
		{ //遍历contig链表
			if(contig->ridposnum > 0)
			{ //找到标记有成功reads的contig节点
				ridposorient = contig->pridposorientation;
				numridposorient = contig->ridposnum;
				for(i=0; i<numridposorient; i++)
				{
					tmp_contigIndex = contig->index - ridposorient[i].matchnum + 1;
					if(minContigIndex>tmp_contigIndex)
					{
						minContigIndex = tmp_contigIndex;
					}
				}
				new_successContig = contig;
			}
			contig = contig->next;
		}

		if(successContig!=NULL && minContigIndex>successContig->index-MIN_OVERLAP_LEN+1)
		{ //最小交叠 MIN_OVERLAP_LEN 个碱基
			return NULL;
		}
	}else
	{ //第二轮拼接过程中查找successContig

		int checkNum = 0;
		tmp_contigIndex = 0;
		contig = contig36;
		while(contig)
		{ //遍历contig链表
			checkNum ++;
			//if(checkNum>5)  //限制检测5次
			if(checkNum>errorRegLenEnd3+1)
				break;

			if(contig->ridposnum > 0)
			{ //找到标记有成功reads的contig节点
				tmp_contigIndex = contig->index;
				break;
			}

			contig = contig->next;
		}

		if(successContig!=NULL && tmp_contigIndex>0 && tmp_contigIndex>successContig->index-MIN_OVERLAP_LEN+1)
		{
			return NULL;
		}


		int maxContigIndex = 0;
		int maxMatchNum = 0;
		checkNum = 0;
		tmp_contigIndex = 1;
		contig = contig36;
		while(contig)
		{ //遍历contig链表
			checkNum ++;
			//if(checkNum>5)  //限制检测5次
			if(checkNum>errorRegLenEnd3+1)
				break;

			if(contig->ridposnum > 0)
			{ //找到标记有成功reads的contig节点
				maxMatchNum = 0;
				ridposorient = contig->pridposorientation;
				numridposorient = contig->ridposnum;
				for(i=0; i<numridposorient; i++)
				{
					if(maxMatchNum < ridposorient[i].matchnum)
					{
						maxMatchNum = ridposorient[i].matchnum;
					}
				}

				//找最大的successContigIndex
				if(tmp_contigIndex+maxMatchNum-1 > maxContigIndex)
					maxContigIndex = tmp_contigIndex + maxMatchNum - 1;
			}
			contig = contig->next;
			tmp_contigIndex ++;
		}

		// ############################ Debug information ##############################
		if(maxContigIndex<=0)
		{
			printf("line=%d, In %s(), assemblyRound=%d, the maxContigIndex==%d\n", __LINE__, __func__, assemblyRound, maxContigIndex);
			return NULL;
		}
		// ############################ Debug information ##############################

		//找最大index的成功contig节点
		contig = contig36;
		tmp_contigIndex = 1;
		while(contig)
		{
			if(tmp_contigIndex==maxContigIndex)
			{ //找到最大index的成功的reads的contig节点
				new_successContig = contig;
				break;
			}
			contig = contig->next;
			tmp_contigIndex ++;
		}
	}

	//检测新的successContig是否合理
	if(successContig==NULL)
	{
		return new_successContig;
	}else if(successContig->index > new_successContig->index)
	{
		return successContig;
	}else
	{
		return new_successContig;
	}
}

/**
 * Recover the deleted reads in De Bruijn graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short recoverDeledReads(contigtype *startContig)
{
	char seq[1000] = {0}, reversed_seq[1000] = {0}; //正向的碱基序列和反向互补的碱基序列
	contigtype *contig = startContig;
	successRead_t *read = NULL;
	int i = 0, j = 0, seq_len = 0;
	//char tmp[READ_LEN+1] = {0}; //保存反向互补read

	//取得contig链表的碱基序列
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
		contig = contig->next;

		// ############################ Debug information ##############################
		if(i>=1000)
		{
			printf("line=%d, In %s(), such a long contig here, contig len>=%d, error!\n", __LINE__, __func__, i);
			return FAILED;
		}
		// ############################ Debug information ##############################
	}
	seq[i] = '\0';
	seq_len = i;
	if(getReversedSeq(reversed_seq, seq, seq_len)==FAILED)
	{ //取得该contig的反向互补的碱基序列
		printf("line=%d, In %s(), cannot get the reverse contig sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	i = 0;
	contig = startContig;
	while(contig)
	{ //循环恢复以该contig节点开始的reads
		read = contig->pridposorientation;
		for(j=0; j<contig->ridposnum; j++)
		{ //恢复每一条read

			// ############################ Debug information ##############################
			//if(read->rid==2556467)
			//{
			//	printf("line=%d, In %s(), rid=%d, pos=%d, matchnum=%d, orientation=%c, i=%d\n", __LINE__, __func__, read->rid, read->pos, read->matchnum, read->orientation, i);
			//}
			// ############################ Debug information ##############################

			if(read->orientation==ORIENTATION_PLUS)
			{ //该read为正向
				if(recoverReadFromGraph(seq+i, read->rid, deBruijnGraph)==FAILED)
				{ //恢复该read
					printf("line=%d, In %s(), cannot recover the read (%lu,%d,%d,%c), i=%d, error!\n", __LINE__, __func__, read->rid, read->pos, read->matchnum, read->orientation, i);
					return FAILED;
				}
			}else
			{ //该read为反向
				//getReverseSeq(tmp, seq+index, read[j].matchnum);
				if(recoverReadFromGraph(reversed_seq+seq_len-i-read->matchnum, read->rid, deBruijnGraph)==FAILED)
				{ //恢复该read
					printf("line=%d, In %s(), cannot recover the read (%lu,%d,%d,%c), i=%d, error!\n", __LINE__, __func__, read->rid, read->pos, read->matchnum, read->orientation, i);
					return FAILED;
				}
			}

			read ++;
		}

		contig = contig->next;
		i++;
	}

	return SUCCESSFUL;
}

/**
 * Recover a read from De Bruijn graph.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short recoverReadFromGraph(char *seq, uint64_t rid, graphtype *graph)
{
	//计算哈希值
	uint64_t hashcode;
	uint64_t tmp_kmerseq[entriesPerKmer];
	int Num, j, rpos, baseInt;

	Num = readLen - kmerSize + 1;

	// generate the kmer integer sequence
	if(generateKmerSeqInt(tmp_kmerseq, seq)==FAILED)
	{
		printf("line=%d, In %s(), cannot generate the kmer integer sequence, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//开始恢复该read中的kmers
	rpos = 1;
	while(rpos<=Num)
	{
		// 恢复该read的一个kmer, 也即是该read的一个ridpos位置
		hashcode = kmerhashInt(tmp_kmerseq);
		if(recoverKmerByHash(hashcode, tmp_kmerseq, rid, rpos, graph)==FAILED)
		{
			if(recoverRemainedKmers(seq+rpos-1, tmp_kmerseq, rid, rpos, graph)==FAILED)
			{
				printf("line=%d, In %s(), can not recover the read %lu from position %d. Error!\n", __LINE__, __func__, rid, rpos);
				return FAILED;
			}
			return SUCCESSFUL; //该read的恢复操作终止, 并直接返回成功标记
		}

		if(rpos==Num || seq[kmerSize+rpos-1]=='\0')
		{ //该read中的kmers已经全部被恢复, 或者seq已经到达尾部, 则退出循环
			break;
		}

		//取下一个kmer
		switch(seq[kmerSize+rpos-1])
		{
			case 'A':
			case 'a':
				baseInt = 0;
				break;
			case 'C':
			case 'c':
				baseInt = 1;
				break;
			case 'G':
			case 'g':
				baseInt = 2;
				break;
			case 'T':
			case 't':
				baseInt = 3;
				break;
			default:
				printf("line=%d, In %s(), unknown base, Error!\n", __LINE__, __func__);
				return FAILED;
		}

		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

		rpos++;
	}

	//在seq长度小于36个碱基的情况下, 将read中的其他kmers恢复
	rpos++;
	while(rpos<=Num)
	{
		//下一个kmer的碱基的整数表示
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = (tmp_kmerseq[entriesPerKmer-1] << 2) & lastEntryMask;

		for(j=0; j<4; j++)
		{ //循环探测该kmer, 直到找到并将其恢复为止
			hashcode = kmerhashInt(tmp_kmerseq);
			if(recoverKmerByHash(hashcode, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
				break;
			tmp_kmerseq[entriesPerKmer-1] ++;
		}

		// ############################ Debug information ##############################
		if(j==4)
		{ //不能恢复该kmer, 打印出错信息
			printf("line=%d, In %s(), cannot recover the kmer (%lu,%d), Error!\n", __LINE__, __func__, rid, rpos);
			return FAILED;
		}
		// ############################ Debug information ##############################

		rpos++;
	}

	return SUCCESSFUL;
}

/**
 * Recover the remained kmers.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short recoverRemainedKmers(char *seq, uint64_t *tmp_kmerseq, uint64_t rid, uint16_t rpos, graphtype *graph)
{
	//short len = strlen(seq);
	int i, j, baseInt;
	uint64_t hash;

	//恢复该read的从rpos开始的第1个kmer
	for(j=0; j<3; j++)
	{
		if((tmp_kmerseq[entriesPerKmer-1] & 3)==3)
			tmp_kmerseq[entriesPerKmer-1] -= 3;
		else
			tmp_kmerseq[entriesPerKmer-1] ++;

		hash = kmerhashInt(tmp_kmerseq);
		if(recoverKmerByHash(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
			break;
	}

	// ############################ Debug information ##############################
	if(j==3)
	{ //该kmer不能恢复, 打印出错信息
		printf("line=%d, In %s(), can not recover this kmer (%lu, %d) Error!\n", __LINE__, __func__, rid, rpos);
		return FAILED;
	}
	// ############################ Debug information ##############################

	//恢复该read的从rpos开始的其他kmer
	i = 1;
	rpos ++;
	while(seq[i+kmerSize-1] && rpos<=readLen-kmerSize+1)
	{
		//取下一个kmer
		switch(seq[i+kmerSize-1])
		{
			case 'A':
			case 'a':
				baseInt = 0;
				break;
			case 'C':
			case 'c':
				baseInt = 1;
				break;
			case 'G':
			case 'g':
				baseInt = 2;
				break;
			case 'T':
			case 't':
				baseInt = 3;
				break;
			default:
				printf("line=%d, In %s(), unknown base, Error!\n", __LINE__, __func__);
				return FAILED;
		}

		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = ((tmp_kmerseq[entriesPerKmer-1] << 2) | baseInt) & lastEntryMask;

		//恢复该kmer
		hash = kmerhashInt(tmp_kmerseq);
		if(recoverKmerByHash(hash, tmp_kmerseq, rid, rpos, graph)==FAILED)
		{ //恢复该kmer失败, 该kmer需要重新确定
			//重新确定该kmer
			for(j=0;j<3;j++)
			{
				if((tmp_kmerseq[entriesPerKmer-1] & 3)==3)
					tmp_kmerseq[entriesPerKmer-1] -= 3;
				else
					tmp_kmerseq[entriesPerKmer-1] ++;

				hash = kmerhashInt(tmp_kmerseq);
				if(recoverKmerByHash(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
					break;
			}

			// ############################ Debug information ##############################
			if(j==3)
			{ //该kmer不能恢复, 打印出错
				printf("line=%d, In %s(), can not recover this kmer (%lu, %d) Error!\n", __LINE__, __func__, rid, rpos);
				return FAILED;
			}
			// ############################ Debug information ##############################
		}
		i++;
		rpos++;
	}

	//如果存在剩余的kmers, 则将其恢复
	while(rpos <= readLen-kmerSize+1)
	{
		//下一个kmer的碱基的整数表示
		if(entriesPerKmer>=2)
		{
			for(j=0; j<entriesPerKmer-2; j++)
			{
				tmp_kmerseq[j] = (tmp_kmerseq[j] << 2) | (tmp_kmerseq[j+1] >> 62);
			}
			tmp_kmerseq[entriesPerKmer-2] = (tmp_kmerseq[entriesPerKmer-2] << 2) | (tmp_kmerseq[entriesPerKmer-1] >> (2*lastEntryBaseNum-2));
		}
		tmp_kmerseq[entriesPerKmer-1] = (tmp_kmerseq[entriesPerKmer-1] << 2) & lastEntryMask;

		for(j=0; j<4; j++)
		{ //循环探测该kmer, 直到找到并将其恢复为止
			hash = kmerhashInt(tmp_kmerseq);
			if(recoverKmerByHash(hash, tmp_kmerseq, rid, rpos, graph)==SUCCESSFUL)
				break;
			tmp_kmerseq[entriesPerKmer-1] ++;
		}

		// ############################ Debug information ##############################
		if(j==4)
		{
			printf("line=%d, In %s(), can not recover this kmer (%lu, %d) Error!\n", __LINE__, __func__, rid, rpos);
			return FAILED;
		}
		// ############################ Debug information ##############################

		rpos++;
	}

	return SUCCESSFUL;
}

/**
 * Get the reversed sequence.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short getReversedSeq(char *reversed_seq, char *seq, int seq_len)
{
	short i = 0;
	for(; i<seq_len; i++)
	{
		switch(seq[i])
		{
			case 'A':
			case 'a':
				reversed_seq[seq_len-i-1] = 'T';
				break;
			case 'C':
			case 'c':
				reversed_seq[seq_len-i-1] = 'G';
				break;
			case 'G':
			case 'g':
				reversed_seq[seq_len-i-1] = 'C';
				break;
			case 'T':
			case 't':
				reversed_seq[seq_len-i-1] = 'A';
				break;
			default:
				printf("line=%d, In %s(), unknown base: %c, Error!\n", __LINE__, __func__, seq[i]);
				return FAILED;
		}
	}

	return SUCCESSFUL;
}


/**
 * Get the contig tail.
 */
contigtype *getContigtail(contigtype *startContig)
{
	if(startContig==NULL)
	{
		printf("In getContigtail(), the startContig==NULL, Error!\n");
		return NULL;
	}

	contigtype *contig = startContig;
	while(contig->next)
	{
		contig = contig->next;
	}

	return contig;
}

/**
 * Initialize the second round assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; if it cannot do second round assembly, return FAILED;
 *  	otherwise, return ERROR.
 */
short initSecondAssembly()
{
	//更新contig链表两端的碱基, 并更新contigIndex到合适的值
	if(successContig)
	//if(successContig!=NULL && successContig!=contigtail)
	{ //第一轮拼接有成功的reads, 则进行第二轮的拼接
		//删除contig链表两端的不准确的碱基
		if(trimContigBeforeCycle2(&contighead, &contigtail, &successContig, &contigIndex)==FAILED)
		{
			printf("In buildContigds(), contigsNum==%d, contigIndex=%d, cannot trim the Contig before 2nd round assembly, Error!\n", contigsNum+1, contigIndex);
			return ERROR;
		}
	}else if(successContig==NULL)
	{ //第一轮拼接没有成功的reads, 则该contig的拼接结束, 并且该contig链表拼接失败

		// ############################ Debug information ##############################
		//if(contigIndex>=CONTIG_LEN_THRESHOLD)
		//{
		//	printf("line=%d, In %s(), contigsNum=%d, contigIndex=%d, the successContig==NULL, error!\n", __LINE__, __func__, contigsNum+1, contigIndex);
		//	return ERROR;
		//}
		// ############################ Debug information ##############################

		return FAILED;
	}

/*
	if(successContig==NULL)
	{ //第一轮拼接没有成功的reads, 则该contig的拼接结束, 并且该contig链表拼接失败
		break;
	}
*/

	//====================================================
	// trim a read length of contig nodes at tail
	if(contigIndex>=3*readLen)
	{
		if(trimContigTailByReadLen(contighead, &contigtail, &successContig, &contigIndex, FIRST_ROUND_ASSEMBLY)==FAILED)
		{
			printf("line=%d, In %s(), cannot trim contig nodes at contig tail by a read length, error!\n", __LINE__, __func__);
			return ERROR;
		}
	}
	//====================================================

	//更新successContig到contig头节点上
	successContig = contighead;

	//更新contig36
	int i = 1;
	contigtype *tmp_contig = contighead;
	//取得第36个contig节点
	while(tmp_contig)
	{
		if(tmp_contig->next==NULL || i==readLen)
			break;
		tmp_contig = tmp_contig->next;
		i++;
	}
	contig36 = tmp_contig;  //将找到的节点赋给contig36

	// ############################ Debug information ##############################
	if(contig36==NULL)
	{
		printf("line=%d, In %s(), contigsNum=%d, contigIndex=%d, contig36==NULL\n", __LINE__, __func__, contigsNum+1, contigIndex);
		return FAILED;
	}
	// ############################ Debug information ##############################

	//取得第二轮拼接的第一个kmer及其碱基序列, 如果不存在该kmer, 则该contig拼接结束
	if(getSecondAssemblyFirstKmers(contighead, deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), contigsNum=%d, cannot get the second round assembly first kmer, Error!\n", __LINE__, __func__, contigsNum+1);
		return FAILED;
	}

	//reverse contig nodes
	if(reverseContig(&contighead, &contigtail)==FAILED)
	{
		printf("line=%d, In %s(), contigsNum=%d, cannot reverse Contig, Error!\n", __LINE__, __func__, contigsNum+1);
		return ERROR;
	}

	// update lastseq36
	if(getLastseq(lastseq36, contig36)==FAILED)
	{
		printf("line=%d, In %s(), contigsNum=%d, cannot get lastseq36, Error!\n", __LINE__, __func__, contigsNum+1);
		return ERROR;
	}

	if(PEGivenType>NONE_PE_GIVEN_TYPE && contigIndex>=minContigLenUsingPE)
	{
		//initialize PEhashTable
		if(initPEHashtableSecondAssembly(contighead, contigIndex-readLen+kmerSize)==FAILED)
		{
			if(cleanReadsFromPEHashtable()==FAILED)
			{
				printf("In %s, cannot clean PE hash table, error!\n", __func__);
				return FAILED;
			}
			printf("In %s, cannot initialize the PE hash table before second round assembly, error!\n", __func__);
			return FAILED;
		}
	}

	//添加reads到决策表中
	if(initAssemblingTableSecondAssembly(lastseq36, deBruijnGraph)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize assembly table when the second assembly, error!\n", __LINE__, __func__);
		return ERROR;
	}
	if(itemNumDecisionTable<0)
	{
		printf("line=%d, In buildContigs(), cannot init assembly table when the second assembly, error!\n", __LINE__);
		return ERROR;
	}

	// initialize the reads number region
	if(contigIndex>=minContigLenCheckingReadsNum)
	{
		if(initReadsNumRegSecondAssembly(contigIndex)==FAILED)
		{
			printf("line=%d, In %s(), cannot initialize the reads number region, error!\n", __LINE__, __func__);
			return FAILED;
		}
	}

	if(setEmptyNaviOccQueue(naviOccQueue, &itemNumNaviOccQueue, &frontRowNaviOccQueue, &rearRowNaviOccQueue)==FAILED)
	{
		printf("line=%d, In %s(), cannot initialize the empty navigation occurrence queue, error!\n", __LINE__, __func__);
		return FAILED;
	}

	//将锁定的reads数量置为0
	lockedReadsNum = 0;

	if(PEGivenType>NONE_PE_GIVEN_TYPE && contigIndex>=minContigLenUsingPE)
		allowedUpdatePEHashArrFlag = NO;

	return SUCCESSFUL;
}

/**
 * Check the reads number in the reads number region.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int updateReadsNumReg(int itemNumSuccessReadsArr, int contigNodesNum, int assemblyRound)
{
	int contigIndexLeft, contigIndexRight;
	contigtype *contig;

	if(contigNodesNum>minContigLenCheckingReadsNum)
	{ // only slide the region

		// ############################ Debug information ########################
		if(regLenReadsNumReg<maxRegLenReadsNumReg)
		{
			printf("line=%d, In %s(), regLenReadsNumReg=%d < maxRegLenReadsNumReg=%d, error!\n", __LINE__, __func__, regLenReadsNumReg, maxRegLenReadsNumReg);
			return FAILED;
		}
		// ############################ Debug information ########################

		// get the number of reads leaving the region, and subtract the number
		readsNumReadsNumReg -= leftContigReadsNumReg->ridposnum;
		leftContigReadsNumReg = leftContigReadsNumReg->next;

		// get the number of reads entering the region, and add the number
		readsNumReadsNumReg += itemNumSuccessReadsArr;
		rightContigReadsNumReg = rightContigReadsNumReg->next;

		readsNumTotal += itemNumSuccessReadsArr;


	}else if(contigNodesNum==minContigLenCheckingReadsNum)
	{ // initialize the region with the maximal size
		readsNumTotal = readsNumReadsNumReg = 0;
		regLenReadsNumReg = 0;
		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			contigIndexLeft = contigNodesNum - maxRegLenReadsNumReg + 1;
			contigIndexRight = contigNodesNum;
		}else
		{ // the second round
			contigIndexLeft = contigNodesNum - maxRegLenReadsNumReg + 1 - readLen + 1;
			contigIndexRight = contigNodesNum - readLen + 1;
		}

		if(contigIndexLeft<=0)
			contigIndexLeft = 1;
		if(contigIndexRight<contigIndexLeft)
			contigIndexRight = contigIndexRight;

		regLenReadsNumReg = contigIndexRight - contigIndexLeft + 1;

		// ############################ Debug information ########################
		if(regLenReadsNumReg!=maxRegLenReadsNumReg)
		{
			printf("line=%d, In %s(), regLenReadsNumReg=%d != maxRegLenReadsNumReg=%d, error!\n", __LINE__, __func__, regLenReadsNumReg, maxRegLenReadsNumReg);
			return FAILED;
		}
		// ############################ Debug information ########################

		contig = contighead;
		while(contighead)
		{
			if(contig->index==contigIndexLeft)
				break;
			contig = contig->next;
		}
		leftContigReadsNumReg = contig;
		readsNumReadsNumReg = 0;
		while(contig)
		{
			readsNumReadsNumReg += contig->ridposnum;
			if(contig->index==contigIndexRight)
				break;
			contig = contig->next;
		}
		rightContigReadsNumReg = contig;

		readsNumTotal = 0;
		contig = contighead;
		while(contig)
		{
			readsNumTotal += contig->ridposnum;
			contig = contig->next;
		}
	}

	// compute the ratio of read numbers, and decide whether the assembly should be terminated
	readsNumRatio = (double)readsNumReadsNumReg * contigNodesNum / (readsNumTotal * regLenReadsNumReg);
	if(readsNumRatio>maxReadsNumRatioThres || readsNumRatio<minReadsNumRatioThres)
	//if(readsNumRatio>maxReadsNumRatioThres)
	{
		solvedRepeatsNum ++;
		readsNumTotal = readsNumReadsNumReg = 0;
		leftContigReadsNumReg = rightContigReadsNumReg = NULL;
		kmers[0] = kmers[1] = NULL;
	}

	return SUCCESSFUL;
}


/**
 * Initialize the reads number region when second assembly.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int initReadsNumRegSecondAssembly(int contigNodesNum)
{
	int contigIndexLeft, contigIndexRight;
	contigtype *contig;

	if(contigNodesNum>=minContigLenCheckingReadsNum)
	{ // initialize the region with the maximal size

		readsNumTotal = readsNumReadsNumReg = 0;
		regLenReadsNumReg = 0;
		if(assemblyRound==FIRST_ROUND_ASSEMBLY)
		{ // the first round
			contigIndexLeft = contigNodesNum - maxRegLenReadsNumReg + 1;
			contigIndexRight = contigNodesNum;
		}else
		{ // the second round
			contigIndexLeft = contigNodesNum - maxRegLenReadsNumReg + 1 - readLen + 1;
			contigIndexRight = contigNodesNum - readLen + 1;
		}

		if(contigIndexLeft<=0)
			contigIndexLeft = 1;
		if(contigIndexRight<contigIndexLeft)
			contigIndexRight = contigIndexRight;

		regLenReadsNumReg = contigIndexRight - contigIndexLeft + 1;

		// ############################ Debug information ########################
		if(regLenReadsNumReg!=maxRegLenReadsNumReg)
		{
			printf("line=%d, In %s(), regLenReadsNumReg=%d != maxRegLenReadsNumReg=%d, error!\n", __LINE__, __func__, regLenReadsNumReg, maxRegLenReadsNumReg);
			return FAILED;
		}
		// ############################ Debug information ########################

		contig = contighead;
		while(contighead)
		{
			if(contig->index==contigIndexLeft)
				break;
			contig = contig->next;
		}
		leftContigReadsNumReg = contig;
		readsNumReadsNumReg = 0;
		while(contig)
		{
			readsNumReadsNumReg += contig->ridposnum;
			if(contig->index==contigIndexRight)
				break;
			contig = contig->next;
		}
		rightContigReadsNumReg = contig;

		readsNumTotal = 0;
		contig = contighead;
		while(contig)
		{
			readsNumTotal += contig->ridposnum;
			contig = contig->next;
		}

	}

	// compute the ratio of read numbers, and decide whether the assembly should be terminated
	readsNumRatio = (double)readsNumReadsNumReg * contigNodesNum / (readsNumTotal * regLenReadsNumReg);
	if(readsNumRatio>maxReadsNumRatioThres || readsNumRatio<minReadsNumRatioThres)
	//if(readsNumRatio>maxReadsNumRatioThres)
	{
		readsNumTotal = readsNumReadsNumReg = 0;
		leftContigReadsNumReg = rightContigReadsNumReg = NULL;
		kmers[0] = kmers[1] = NULL;

		//====================================================
		// trim a read length of contig nodes at tail
		if(contigIndex>=3*readLen)
		{
			if(trimContigTailByReadLen(contighead, &contigtail, &successContig, &contigIndex, SECOND_ROUND_ASSEMBLY)==FAILED)
			{
				printf("line=%d, In %s(), cannot trim contig nodes at contig tail by a read length, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}
		//====================================================
	}

	return SUCCESSFUL;
}
