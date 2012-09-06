#include "inc/stdinc.h"
#include "inc/extvab.h"


/**
 * Update decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int updateDecisionTable(kmertype *tmp_kmers[2])
{
	int i, j, startIndex, exectIndex, posNum;
	ridpostype *ridpostable = NULL;  //ridpostable表
	ridpostype rid_pos = {0};            //ridpostable表中的一行
	posNum = 0;
	assemblingreadtype *this_assemblingRead = decisionTable;

	// 遍历决策表，更新决策表中存在的read信息
	for(i=0; i<itemNumDecisionTable; i++)//遍历决策表
	{

		// ######################### Debug information ###########################
		//if(this_assemblingRead->rid==212805)
		//{
		//	printf("line=%d, In %s(), rid=%lu\n", __LINE__, __func__, (uint64_t)this_assemblingRead->rid);
		//}
		// ######################### Debug information ###########################

		for(j=0; j<2; j++)
		{ //j==0为正向kmer, j==1为反向互补kmer
			if(tmp_kmers[j] && ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS) || (j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS)))
			{ //read为正向且正向kmer不为空; 或read为反向且反向互补kmer不为空
				//判断ridpostable表中该read是否存在
				ridpostable = tmp_kmers[j]->ppos;
				posNum = tmp_kmers[j]->arraysize;
				startIndex = findStartIndex(this_assemblingRead->rid, ridpostable, posNum);//二分查找rid_pos_table
				if(startIndex>=0)
				{  //存在, 继续查找精确位置
					if(this_assemblingRead->lastpos>0)
					{ //该read上次拼接出现
						if(j==0) //正向kmer
							exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->lastpos+1, startIndex, ridpostable, posNum);
						else //反向互补kmer
							exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->lastpos-1, startIndex, ridpostable, posNum);
					}else
					{  //该read上次拼接未出现
						if(j==0) //正向kmer
							exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->firstpos+this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes, startIndex, ridpostable, posNum);
						else //反向互补kmer
							exectIndex = getExectIndex(this_assemblingRead->rid, this_assemblingRead->firstpos-this_assemblingRead->kmerappeartimes-this_assemblingRead->kmerunappeartimes, startIndex, ridpostable, posNum);
					}

					if(exectIndex>=0)
					{  //该read位置合理
						rid_pos = ridpostable[exectIndex];
						if(rid_pos.delsign==0) //该read未被删除
						{
							if(this_assemblingRead->lastpos!=0)
							{  //该read上次拼接出现, 则按照该read当次能够拼接处理
								this_assemblingRead->kmerappeartimes++;
								this_assemblingRead->lastpos = rid_pos.pos;
								this_assemblingRead->lastappearpos = rid_pos.pos;
							}else
							{  //该read上次拼接未出现, 则按照该read当次能够拼接处理
								//if(assemblingreadtable[i].kmerunapperblocks>0)
								//{  //连续未出现块数>0, 即该连续块在此结束
									//assemblingreadtable[i].kmerappeartimes += assemblingreadtable[i].kmerunappeartimes + 1;
									//assemblingreadtable[i].kmerunappeartimes = 0;

									this_assemblingRead->kmerappeartimes++;
									this_assemblingRead->lastpos = rid_pos.pos;
									this_assemblingRead->lastappearpos = rid_pos.pos;

									//count the number_of_corrected_reads
									number_of_corrected_reads++;
								//}
							}
							ridpostable[exectIndex].reserved = 1;
						}else
						{ //该read已经被删除
							this_assemblingRead->delsign = 1;
						}
					}else
					{ //该read位置不合理, 按照该read当次不能够拼接处理
						if(this_assemblingRead->lastpos!=0)
						{ //该read上次拼接出现, 当次拼接未出现, 开始新的未出现连续块, 未出现次数重新计算
							this_assemblingRead->kmerappeartimes += this_assemblingRead->kmerunappeartimes;
							this_assemblingRead->kmerunappeartimes = 1;
							this_assemblingRead->lastpos = 0;
							this_assemblingRead->kmerunappearblocks++;
						}else
						{ //该read上次拼接两次未出现
							this_assemblingRead->kmerunappeartimes++;
						}

					} // end if(exectIndex>=0)

				}else //在ridpos表中该read不存在，按照该read当次不能拼接处理
				{

					if(this_assemblingRead->lastpos!=0)
					{ //该read上次拼接出现, 当次拼接未出现, 开始新的未出现连续块, 未出现次数重新计算
						this_assemblingRead->kmerappeartimes += this_assemblingRead->kmerunappeartimes;
						this_assemblingRead->kmerunappeartimes = 1;
						this_assemblingRead->lastpos = 0;
						this_assemblingRead->kmerunappearblocks++;
					}else
					{
						this_assemblingRead->kmerunappeartimes++;
					}

/*
					if(assemblingreadtable[i].lastpos>0)
					{  //该read上次拼接出现，当次拼接未出现, 开始新的未出现连续块
						assemblingreadtable[i].kmerunappeartimes++;
						assemblingreadtable[i].lastpos = 0;
						assemblingreadtable[i].kmerunapperblocks++;  //新块数+1
					}else
					{  //该read上次拼接两次未出现
						assemblingreadtable[i].kmerunappeartimes++;
					}
*/
				}
			}else if(tmp_kmers[j]==NULL && ((j==0 && this_assemblingRead->orientation==ORIENTATION_PLUS) || (j==1 && this_assemblingRead->orientation==ORIENTATION_MINUS)))
			{ //read为正向且正向kmer为空; 或read为反向且反向互补kmer为空
				//按照kmer不能拼接处理
				if(this_assemblingRead->lastpos!=0)
				{ //该read上次拼接出现, 当次拼接未出现, 开始新的未出现连续块, 未出现次数重新计算
					this_assemblingRead->kmerappeartimes += this_assemblingRead->kmerunappeartimes;
					this_assemblingRead->kmerunappeartimes = 1;
					this_assemblingRead->lastpos = 0;
					this_assemblingRead->kmerunappearblocks++;
				}else
				{
					this_assemblingRead->kmerunappeartimes++;
				}
			}
		}

		this_assemblingRead ++;
	}

	// 遍历kmer的ridpos位置表格，将合理的新出现的reads添加进决策表
	int /*newNum = 0,*/ pos = 0;
	int returnCode, matedFlag;
	//if(PEGivenType>NONE_PE_GIVEN_TYPE && readsNumInPEHashArr>=MIN_READ_NUM_PE_HASH)
	if(PEGivenType>=NONE_PE_GIVEN_TYPE)
	{
		if(tmp_kmers[0])
		{ //正向kmer不为空
			ridpostable = tmp_kmers[0]->ppos;
			posNum = tmp_kmers[0]->arraysize;
			for(i=0; i<posNum; i++)//再次遍历kmer
			{
				pos = ridpostable[i].pos;
				if(ridpostable[i].delsign==0 && ridpostable[i].reserved==0 && pos==1)
				{  //未被删除, 未标记并且位置为5'端起始位置, 则将该read添加进决策表

					// ############################ Debug information ##############################
					//if(ridpostable[i].rid==212805)
					//{
					//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_PLUS, NO)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
					}
					//newNum++;
				}//end if
				ridpostable[i].reserved = 0;//标志位清零
			}//end for
		}

		if(tmp_kmers[1])
		{ //反向互补kmer不为空
			ridpostable = tmp_kmers[1]->ppos;
			posNum = tmp_kmers[1]->arraysize;
			for(i=0; i<posNum; i++)//再次遍历kmer
			{
				pos = ridpostable[i].pos;
				//this_assemblingRead = &assemblingreadtable[assemblingtable_size];
				if(ridpostable[i].delsign==0 && ridpostable[i].reserved==0 && pos>=readLen-kmerSize+1-errorRegLenEnd3)
				{  //未被删除, 未标记并且位置为3'端起始位置, 则将该read添加进决策表

					// ############################ Debug information ##############################
					//if(ridpostable[i].rid==212805)
					//{
					//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					returnCode = validReadPair(ridpostable[i].rid);
					if(returnCode==YES)
					{
						matedFlag = YES;
					}else if(returnCode==NO)
					{
						matedFlag = NO;
					}else
					{
						printf("line=%d, In %s(), cannot valid read pair, error!\n", __LINE__, __func__);
						return FAILED;
					}
					if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_MINUS, matedFlag)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
					}

					//newNum++;

				}//end if
				ridpostable[i].reserved = 0;//标志位清零
			}//end for
		}

	}else
	{  // single end assembly
		if(tmp_kmers[0])
		{ //正向kmer不为空
			ridpostable = tmp_kmers[0]->ppos;
			posNum = tmp_kmers[0]->arraysize;
			for(i=0; i<posNum; i++)//再次遍历kmer
			{
				pos = ridpostable[i].pos;
				if(ridpostable[i].delsign==0 && ridpostable[i].reserved==0 && pos==1)
				{  //未被删除, 未标记并且位置为5'端起始位置, 则将该read添加进决策表

					// ############################ Debug information ##############################
					//if(ridpostable[i].rid==212805)
					//{
					//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_PLUS, NO)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
					}
					//newNum++;

				}//end if
				ridpostable[i].reserved = 0;//标志位清零
			}//end for
		}

		if(tmp_kmers[1])
		{ //反向互补kmer不为空
			ridpostable = tmp_kmers[1]->ppos;
			posNum = tmp_kmers[1]->arraysize;
			for(i=0; i<posNum; i++)//再次遍历kmer
			{
				pos = ridpostable[i].pos;
				//this_assemblingRead = &assemblingreadtable[assemblingtable_size];
				if(ridpostable[i].delsign==0 && ridpostable[i].reserved==0 && pos>=readLen-kmerSize+1-errorRegLenEnd3)
				{  //未被删除, 未标记并且位置为3'端起始位置, 则将该read添加进决策表

					// ############################ Debug information ##############################
					//if(ridpostable[i].rid==212805)
					//{
					//	printf("line=%d, In %s(), rid=%lu, pos=%d, delsign=%d, reserved=%d\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid, ridpostable[i].pos, ridpostable[i].delsign, ridpostable[i].reserved);
					//}
					// ############################ Debug information ##############################

					if(addReadToDecisionTable(ridpostable[i].rid, pos, ORIENTATION_MINUS, NO)==FAILED)
					{
						printf("line=%d, In %s(), cannot add read %lu to decision table, error!\n", __LINE__, __func__, (uint64_t)ridpostable[i].rid);
						return FAILED;
					}
					//newNum++;

				}//end if
				ridpostable[i].reserved = 0;//标志位清零
			}//end for
		}
	}

	// ############################ Debug information ##############################
	//printf("In %s(), new reads Num=%d, numassemblingreads=%d\n", __func__, newNum, itemNumDecisionTable);
	// ############################ Debug information ##############################

	return SUCCESSFUL;
}


/**
 * Update finished reads and record the succefful reads to successReadsArr.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int updateFinishedReadsInDecisionTable()
{
	int i;

	itemNumSuccessReadsArr = 0;
	for(i = 0; i<itemNumDecisionTable; i++)
	{
		if(decisionTable[ i ].status==SUCCESSFUL_STATUS)
		{  //该read成功结束
			successReadsArr[ itemNumSuccessReadsArr ].rid = decisionTable[i].rid;
			successReadsArr[ itemNumSuccessReadsArr ].pos = decisionTable[i].firstpos;
			if(decisionTable[i].orientation==ORIENTATION_PLUS)
			{ //正向reads的kmer匹配数量, 为出现的两端的kmer相减+1
				successReadsArr[ itemNumSuccessReadsArr ].matchnum = decisionTable[i].lastappearpos - decisionTable[i].firstpos + 1;
			}else
			{
				successReadsArr[ itemNumSuccessReadsArr ].matchnum = decisionTable[i].firstpos - decisionTable[i].lastappearpos + 1;
			}
			successReadsArr[ itemNumSuccessReadsArr ].orientation = decisionTable[i].orientation;

			itemNumSuccessReadsArr ++;

			if(itemNumSuccessReadsArr>=maxItemNumSuccessReadsArr)
			{
				if(reallocateSuccessReadsArr()==FAILED)
				{
					printf("line=%d, In %s(), cannot reallocate memory for success reads array, error!\n", __LINE__, __func__);
					return FAILED;
				}
			}
		}else if(decisionTable[ i ].status==FAILED_STATUS)
		{  //该read失败结束
		}

	}

	if(removeFinishedReadsFromDecisionTable()==FAILED)
	{
		printf("line=%d, In %s(), cannot resort decision table, error!\n", __LINE__, __func__);
		return FAILED;
	}

	return SUCCESSFUL;
}

/**
 * Reallocate memory for successful reads array.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int reallocateSuccessReadsArr()
{
	successRead_t *pSuccessReadsArr;

	pSuccessReadsArr = (successRead_t *) malloc(2 * maxItemNumSuccessReadsArr * sizeof(successRead_t));
	if(pSuccessReadsArr==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory for successful reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(memcpy(pSuccessReadsArr, successReadsArr, maxItemNumSuccessReadsArr * sizeof(successRead_t))==NULL)
	{
		printf("line=%d, In %s(), cannot copy memory for successful reads, error!\n", __LINE__, __func__);
		return FAILED;
	}

	free(successReadsArr);
	successReadsArr = pSuccessReadsArr;
	maxItemNumSuccessReadsArr *= 2;

	return SUCCESSFUL;
}

/**
 * Remove the finished element in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
int removeFinishedReadsFromDecisionTable()
{
	int i , j ;
	i = 0;
	j= itemNumDecisionTable - 1;
	while(i <= j)
	{
		if(decisionTable[i].status !=  ASSEMBLING_STATUS && decisionTable[j].status == ASSEMBLING_STATUS)
		{ //第i个为空，第j个不为空时，将用第j个覆盖第i个，然后数组元素个数减1, 尾指针上移一位, 头指针下移一位。
			if(memcpy(decisionTable+i, decisionTable+j, sizeof(assemblingreadtype))==NULL)
			{
				printf("line=%d, In %s(), cannot copy memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			decisionTable[i].status = ASSEMBLING_STATUS;
			itemNumDecisionTable --;  //数组元素个数减1
			j--;  //尾指针上移一位
			i++;  //头指针下移一位
		}

		if(decisionTable[i].status == ASSEMBLING_STATUS)
			i++; //头指针下移一位

		if(decisionTable[j].status != ASSEMBLING_STATUS)
		{
			j--; //尾指针上移一位
			itemNumDecisionTable --;  //数组元素个数减1
		}
	}

	return SUCCESSFUL;
}

/**
 * Update reads status in decision table.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short updateAssemblingreadsStatus()
{
	assemblingreadtype *this_assemblingRead = decisionTable;
	int i;

	//int successNum = 0, failedNum = 0;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		// ######################### Debug information ###########################
		//if(this_assemblingRead->rid==212805)
		//{
		//	printf("line=%d, In %s(), rid=%lu\n", __LINE__, __func__, (uint64_t)this_assemblingRead->rid);
		//}
		// ######################### Debug information ###########################

		if(this_assemblingRead->orientation==ORIENTATION_PLUS)
		{ //正向read
			if(this_assemblingRead->delsign==1)
			{//该read已经被删除，则结束，并且失败；
				this_assemblingRead->status = FAILED_STATUS;
				//failedNum++;
			}else if(this_assemblingRead->kmerunappearblocks>1)
			{  //未出现次数连续块数>1，则结束，并且失败；
				this_assemblingRead->status = FAILED_STATUS;
				//failedNum++;
			}else if(this_assemblingRead->kmerunappearblocks==1)
			{  //未出现次数连续块数==1
				if(this_assemblingRead->lastpos==0 && this_assemblingRead->lastappearpos>=readLen-kmerSize+1-errorRegLenEnd3)
				{//正向read 3'端2个碱基未出现,该read成功结束
					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailed(this_assemblingRead->rid, i);
					//successNum++;
				}else if(this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes==readLen-kmerSize+1)
				{  //拼接次数达到read最大拼接次数, 则结束;
					if(this_assemblingRead->lastpos==readLen-kmerSize+1)
					{  //当次拼接出现，并且该位置为该read的最后一个kmer，则成功;
						this_assemblingRead->status = SUCCESSFUL_STATUS;
						updateSameReadStatusToFailed(this_assemblingRead->rid, i);
						//successNum++;

					}
					else
					{  //否则, 失败;
						this_assemblingRead->status = FAILED_STATUS;
						//failedNum++;
					}
				}else
				{  //拼接次数未达到read最大拼接次数
					if(this_assemblingRead->lastpos>0 && this_assemblingRead->kmerunappeartimes!=kmerSize)
					{  //当次拼接出现且未出现次数!=KMER_SIZE, 则结束，并且失败;
						this_assemblingRead->status = FAILED_STATUS;
						//failedNum++;
					}else if(this_assemblingRead->lastpos==0 && this_assemblingRead->kmerunappeartimes>kmerSize)
					{ //当次拼接未出现，并且未出现次数>KMER_SIZE, 则结束，并且失败;
						this_assemblingRead->status = FAILED_STATUS;
						//failedNum++;
					}
				}
			}else
			{  //未出现次数连续块数==0
				if(this_assemblingRead->kmerappeartimes==readLen-kmerSize+1)
				{ //首位置不为1, 且出现次数==该read的最后一个位置，则结束, 并且成功；
					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailed(this_assemblingRead->rid, i);
					//successNum++;

				}
			}
		}else
		{ //反向read
			if(this_assemblingRead->delsign==1)
			{//该read已经被删除，则结束，并且失败；
				this_assemblingRead->status = FAILED_STATUS;
				//failedNum++;
			}else if(this_assemblingRead->kmerunappearblocks>1)
			{  //未出现次数连续块数>1，则结束，并且失败；
				this_assemblingRead->status = FAILED_STATUS;
				//failedNum++;
			}else if(this_assemblingRead->kmerunappearblocks==1)
			{  //未出现次数连续块数==1
				if(this_assemblingRead->firstpos!=readLen-kmerSize+1)
				{ //首位置不为25, 则结束，并且失败;
					this_assemblingRead->status = FAILED_STATUS;
					//failedNum++;
				}else
				{  //首位置为25
					if(this_assemblingRead->kmerappeartimes+this_assemblingRead->kmerunappeartimes==readLen-kmerSize+1)
					{  //拼接次数达到read最大拼接次数, 则结束
						if(this_assemblingRead->lastpos==1)
						{  //当次拼接出现，并且该位置为该read的最后一个kmer，则成功;
							this_assemblingRead->status = SUCCESSFUL_STATUS;
							updateSameReadStatusToFailed(this_assemblingRead->rid, i);
							//successNum++;

						}
						else
						{  //否则, 失败;
							this_assemblingRead->status = FAILED_STATUS;
							//failedNum++;
						}
					}else
					{  //拼接次数未达到read最大拼接次数
						if(this_assemblingRead->lastpos>0 && this_assemblingRead->kmerunappeartimes!=kmerSize)
						{  //当次拼接出现且未出现次数!=KMER_SIZE, 则结束，并且失败;
							this_assemblingRead->status = FAILED_STATUS;
							//failedNum++;
						}else if(this_assemblingRead->lastpos==0 && this_assemblingRead->kmerunappeartimes>kmerSize)
						{ //当次拼接未出现，并且未出现次数>KMER_SIZE, 则结束，并且失败;
							this_assemblingRead->status = FAILED_STATUS;
							//failedNum++;
						}
					}
				}
			}else
			{  //未出现次数连续块数==0
				if(this_assemblingRead->firstpos!=readLen-kmerSize+1 && this_assemblingRead->lastpos==1)
				{ //首位置不为25, 且最后的pos==该read的最后一个位置，则结束, 并且成功；
					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailed(this_assemblingRead->rid, i);
					//successNum++;


				}else if(this_assemblingRead->firstpos==readLen-kmerSize+1 && this_assemblingRead->kmerappeartimes==readLen-kmerSize+1)
				{  //首位置为25, 且该read所有的kmer都拼接上, 则结束, 并且成功;
					this_assemblingRead->status = SUCCESSFUL_STATUS;
					updateSameReadStatusToFailed(this_assemblingRead->rid, i);
					//successNum++;


				}
			}
		}

		this_assemblingRead ++;
	}

	return SUCCESSFUL;
}

/**
 * 更新决策表中相同rid的reads的状态标记为失败标记.
 *  当决策表中有成功的reads的时候，要执行该操作.
 *
 *   @ return:
 *     成功, 返回成功标记; 失败, 返回失败标记.
 */
short updateSameReadStatusToFailed(uint64_t rid/*assemblingreadtype *assemblingRead*/, int assemblingreadIndex)
{
	assemblingreadtype *this_assemblingRead = decisionTable;
	int i;
	for(i=0; i<itemNumDecisionTable; i++)
	{
		if(this_assemblingRead->rid==rid && i!=assemblingreadIndex)
		{ //rid相同的read, 数组下标不同
			this_assemblingRead->status = FAILED_STATUS;
		}
		this_assemblingRead ++;
	}
	return SUCCESSFUL;
}
