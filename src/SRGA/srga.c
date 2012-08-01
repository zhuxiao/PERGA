
#include "inc/stdinc.h"
#include "inc/extvab.h"

/**
 *  @return:
 *  	If succeeds, return SUCCESSFUL;
 *  	if the job is executed and failed, return ERROR;
 *  	otherwise, return FAILED.
 **/
short startSRGA(int operationModePara, int kmerSizePara, int readLenCutOffPara, int pairedModePara, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, char *graphFilePara, double meanSizeInsertPara, double standardDevPara, char *outputPathPara, char *outputPrefixPara, int minContigLenPara)
{
	struct timeval tp_start,tp_end;
	double time_used;
	gettimeofday(&tp_start,NULL);

	// initialize the global parameters
	//if(initGlobalParas("../output/", "S.pombe_ERR018885", "../reads/SRR023412_1.fastq", "../reads/SRR023412_2.fastq", 21, 50)==FAILED)
	if(initGlobalParas(operationModePara, outputPathPara, outputPrefixPara, readFilesPara, readFileNumPara, singleBaseQualThresPara, pairedModePara, kmerSizePara, readLenCutOffPara, graphFilePara, meanSizeInsertPara, standardDevPara, minContigLenPara)==FAILED)
	{
		//printf("line=%d, In %s(), cannot initialize the global parameters, error!\n", __LINE__, __func__);
		return FAILED;
	}


	if(operationMode!=2)
	{
		//if(constructGraphByPEFastq(graphFile, readFilesFastqPE[0], readFilesFastqPE[1])==FAILED)
		if(constructGraph(graphFile, readFilesInput, readFileNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot construct the graph, Error!\n", __LINE__, __func__);
			return ERROR;
		}
	}

	if(operationMode!=1)
	{
		if(buildContigs(contigsFileFasta, graphFile)==FAILED)
		{ //构建contigs
			printf("line=%d, In %s(), cannot build contigs, error!\n", __LINE__, __func__);
			return ERROR;
		}
	}


	freeGlobalParas();

	if(operationMode==0 || operationMode==2)
	{
//		printf("totalReadNum=%ld, validReadNum=%ld, deledReadNum=%ld, validRatio=%.4f, "
//				"successReadNum=%ld, failedReadNum=%ld, failedRatio=%.4f\n",
//				totalReadNum, validReadNum, totalReadNum-validReadNum, (float)validReadNum/totalReadNum,
//				successReadNum, validReadNum-successReadNum, (float)(validReadNum-successReadNum)/validReadNum);

		//printf("Corrected read number: %d\n", number_of_corrected_reads);
		//printf("Overlaps less than %d: %d\n", MIN_OVERLAP_LEN, number_of_overlap_less_than_threshold);
		//printf("searchNum=%ld, totalSearchLen=%llu, fold=%lld\n", searchNum, totalSearchLen, totalSearchLen/searchNum);
	}

	gettimeofday(&tp_end,NULL);
	time_used = tp_end.tv_sec-tp_start.tv_sec + (double)(tp_end.tv_usec-tp_start.tv_usec)/1000000;

	//printf("Graph Used Mem : %f MB.\n", (double)countMemory/1024/1024);

	printf("\nTotal Used Time: %f Seconds.\n", time_used);


    return SUCCESSFUL;
}


/**
 * Initialize the global parameters.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short initGlobalParas(int operationModePara, char *outputPathName, char *prefix, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, int pairedModePara, int kmerLen, int readLenCut, char *graphFilePara, double meanSizeInsertPara, double standardDevPara, int minContigLenPara)
{
	int i, prefixLen;
	char kmerSizeStr[20], readLenStr[20];
	struct stat st;

	printf("SRGA version : %s\n", VERSION_STR);
	printf("Release date : %s\n", RELEASE_DATE_STR);

	//printf("\n============= Begin setting global parameters, please wait ... =============\n");

	operationMode = operationModePara;

	if(setGlobalPath(outputPathName)==FAILED)
	{
		printf("line=%d, In %s(), cannot set global paths, error!\n", __LINE__, __func__);
		return FAILED;
	}

	if(operationMode==0 || operationMode==1)
	{
		readFileNum = readFileNumPara;
		for(i=0; i<readFileNum; i++)
		{
			readFilesInput[i] = (char *) calloc (256, sizeof(char));
			if(readFilesInput[i]==NULL)
			{
				printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
				return FAILED;
			}
			strcpy(readFilesInput[i], readFilesPara[i]);
		}

		// get reads file format type
		if(getReadsFileFormat(&readsFileFormatType, readFilesInput, readFileNum)==FAILED)
		{
			printf("line=%d, In %s(), cannot get reads file formats, error!\n", __LINE__, __func__);
			return FAILED;
		}

		// get the readlenInFile from fastq file
		if(readsFileFormatType==FILE_FORMAT_FASTA)
		{
			if(getMinReadLenFromFastaFiles(&readLenInFile, readFilesInput, readFileNum)==FAILED)
			{
				//printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}else if(readsFileFormatType==FILE_FORMAT_FASTQ)
		{
			if(getMinReadLenFromFastqFiles(&readLenInFile, readFilesInput, readFileNum)==FAILED)
			{
				//printf("line=%d, In %s(), cannot get readlenInFile, error!\n", __LINE__, __func__);
				return FAILED;
			}
		}

		readLenCutOff = readLenCut;
		if(readLenCutOff==0)
			readLen = readLenInFile;
		else if(readLenCutOff>readLenInFile)
			readLen = readLenInFile;
		else
			readLen = readLenCutOff;

		if(kmerLen==0)
		{
			kmerSize = DEFAULT_KMER_SIZE;
		}else
		{
			kmerSize = kmerLen;
		}

		if(kmerSize%2==0)
			kmerSize --;

		if(kmerSize>readLen)
		{
			printf("Exception: invalid k-mer size, it is larger than the read length.\n");
			return FAILED;
		}

/*
		if(kmerSize>readLen-errorRegLenEnd3)
		{
			kmerSize = readLen - errorRegLenEnd3;
			if(kmerSize%2==0)
				kmerSize --;
		}
*/

		hashTableSize = HASH_TABLE_SIZE;
		//hashTableSize = 1LLU << (kmerSize << 1);

	}else if(operationMode==2)
	{
		strcpy(graphFile, graphFilePara);

		if(stat(graphFile, &st)==-1)
		{
			printf("Exception: please specify correct graph file.\n");
			return FAILED;
		}

		// load the graph to memory
		if(GlobalParasFromGraph(&readLen, &kmerSize, &hashTableSize, &pairedMode, graphFile)==FAILED)
		{
			//printf("line=%d, In %s(), cannot load graph to memory, error!\n", __LINE__, __func__);
			return FAILED;
		}

	}else
	{
		printf("Exception: please specify correct command, error!\n");
		return FAILED;
	}

	qualityBaseNumEnd3 = ceil(readLen * QUAL_BASE_NUM_3End_FACTOR);
	qualityBaseNumEnd5 = readLen - qualityBaseNumEnd3;
	errorRegLenEnd3 = ceil(readLen * ERROR_REGION_LEN_3End_FACTOR);
	if(singleBaseQualThresPara>0)
		singleBaseQualThres = singleBaseQualThresPara;
	else
		singleBaseQualThres = SINGLE_QUAL_THRESHOLD;

	prefixLen = strlen(prefix);

	sprintf(kmerSizeStr, "%d", kmerSize);
	sprintf(readLenStr, "%d", readLen);

	if(operationMode==0 || operationMode==2)
	{
		strcpy(contigsFileFasta, outputPathStr);
		if(prefixLen>0)
		{
			strcat(contigsFileFasta, prefix);
			strcat(contigsFileFasta, "_");
		}
		strcat(contigsFileFasta, "contigs");
		strcat(contigsFileFasta, "_K");
		strcat(contigsFileFasta, kmerSizeStr);
		strcat(contigsFileFasta, "_R");
		strcat(contigsFileFasta, readLenStr);
		strcat(contigsFileFasta, ".fa");
		strcpy(contigsFileHanging, outputPathStr);
		if(prefixLen>0)
		{
			strcat(contigsFileHanging, prefix);
			strcat(contigsFileHanging, "_");
		}
		strcat(contigsFileHanging, "contigs");
		strcat(contigsFileHanging, "_K");
		strcat(contigsFileHanging, kmerSizeStr);
		strcat(contigsFileHanging, "_R");
		strcat(contigsFileHanging, readLenStr);
		strcat(contigsFileHanging, "_hanging.fa");

		strcpy(fragmentSizeFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(fragmentSizeFile, prefix);
			strcat(fragmentSizeFile, "_");
		}
		strcat(fragmentSizeFile, "fragmentSize.bin");

		strcpy(sampleContigsFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(sampleContigsFile, prefix);
			strcat(sampleContigsFile, "_");
		}
		strcat(sampleContigsFile, "sampleContigs.fa");
	}

	if(operationMode==0 || operationMode==1)
	{
		strcpy(graphFile, outputPathStr);
		if(prefixLen>0)
		{
			strcat(graphFile, prefix);
			strcat(graphFile, "_");
		}
		strcat(graphFile, "hashtable");
		strcat(graphFile, "_K");
		strcat(graphFile, kmerSizeStr);
		strcat(graphFile, "_R");
		strcat(graphFile, readLenStr);
		strcat(graphFile, ".bin");
	}

//	printf("outputPathStr=%s\n", outputPathStr);
//	printf("contigsFileFasta=%s\n", contigsFileFasta);
//	printf("contigsFileHanging=%s\n", contigsFileHanging);
//	printf("fragmentSizeFile=%s\n", fragmentSizeFile);
//	printf("graphFile=%s\n", graphFile);
//	printf("sampleContigsFile=%s\n", sampleContigsFile);
//
//	printf("qualityBaseNumEnd3=%d\n", qualityBaseNumEnd3);
//	printf("qualityBaseNumEnd5=%d\n", qualityBaseNumEnd5);
//	printf("errorRegLenEnd3=%d\n", errorRegLenEnd3);
//	printf("AVERAGE_QUAL_THRESHOLD_3End=%.4f\n", AVERAGE_QUAL_THRESHOLD_3End);
//	printf("AVERAGE_QUAL_THRESHOLD_5End=%.4f\n", AVERAGE_QUAL_THRESHOLD_5End);
//	printf("SINGLE_QUAL_THRESHOLD=%d\n", SINGLE_QUAL_THRESHOLD);
//	printf("ARTIFACTS_BASE_A_THRESHOLD=%.4f\n", ARTIFACTS_BASE_A_THRESHOLD);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//PEGivenType = SE_GIVEN_TYPE;			// 0

	//PEGivenType = NONE_PE_GIVEN_TYPE;		// 1

	//PEGivenType = INSERT_PE_GIVEN_TYPE;		// 2
	//meanSizeInsert = 200;

	//PEGivenType = BOTH_PE_GIVEN_TYPE;		// 3
	//meanSizeInsert = 213.59;
	//standardDev = 14.53;

	//kmerSize = 12;
	//kmerSize = 13;
	//kmerSize = 15;
	//kmerSize = 17;
	//kmerSize = 19;
	//kmerSize = 21;
	//kmerSize = 23;
	//kmerSize = 25;
	//kmerSize = 27;
	//kmerSize = 29;
	//kmerSize = 31;

	if(operationMode==0 || operationMode==1)
		pairedMode = pairedModePara;

	if(operationMode==0 || operationMode==2)
	{
		if(pairedMode==0)
		{
			PEGivenType = SE_GIVEN_TYPE;
			meanSizeInsert = 0;
			standardDev = 0;
		}else if(meanSizeInsertPara==0)
		{
			PEGivenType = NONE_PE_GIVEN_TYPE;
			meanSizeInsert = 0;
			standardDev = 0;
		}else if(meanSizeInsertPara!=0 && standardDevPara==0)
		{
			PEGivenType = INSERT_PE_GIVEN_TYPE;
			meanSizeInsert = meanSizeInsertPara;
			standardDev = meanSizeInsert * DRAFT_SDEV_FACTOR;
		}else if(meanSizeInsertPara!=0 && standardDevPara!=0)
		{
			PEGivenType = BOTH_PE_GIVEN_TYPE;
			meanSizeInsert = meanSizeInsertPara;
			standardDev = standardDevPara;
		}

		minContigLen = minContigLenPara;
		if(minContigLen==0)
			minContigLen = CONTIG_LEN_THRESHOLD;
	}

	entriesPerKmer = ((kmerSize-1) / 32) + 1;
	if(kmerSize%32==0)
	{
		lastEntryMask = (uint64_t) -1;
		lastEntryBaseNum = 32;
	}else
	{
		lastEntryMask = (1LLU << ((kmerSize%32)<<1)) - 1;
		lastEntryBaseNum = kmerSize % 32;
	}



	//printf("PEGivenType=%d\n", PEGivenType);
	//if(PEGivenType==INSERT_PE_GIVEN_TYPE)
	//{
		//printf("meanSizeInsert=%.2f\n", meanSizeInsert);
	//}else if(PEGivenType==BOTH_PE_GIVEN_TYPE)
	//{
		//printf("meanSizeInsert=%.2f\n", meanSizeInsert);
		//printf("standardDev=%.2f\n", standardDev);
	//}

	//printf("readLen=%d\n", readLen);
	//printf("kmerSize=%d\n", kmerSize);
	//printf("entriesPerKmer=%d\n", entriesPerKmer);
	//printf("lastEntryBaseNum=%d\n", lastEntryBaseNum);
	//printf("lastEntryMask=0x%lX\n", lastEntryMask);
	//printf("hashTableSize=%lu\n", hashTableSize);


	// allocate memory
	kmerSeqInt = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(kmerSeqInt==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}

	kmerSeqIntRev = (uint64_t*) malloc(entriesPerKmer*sizeof(uint64_t));
	if(kmerSeqIntRev==NULL)
	{
		printf("line=%d, In %s(), cannot allocate memory, error!\n", __LINE__, __func__);
		return FAILED;
	}


	// output the variables
	if(operationMode==0)
	{
		printf("\nkmer size          : %d\n", kmerSize);
		//printf("read length cutoff : %d\n", readLenCutOff);
		printf("read length        : %d\n", readLen);
		printf("paired mode        : %d\n", pairedMode);
		for(i=0; i<readFileNumPara; i++)
		printf("read files[%d]      : %s\n", i, readFilesInput[i]);
		printf("single qual thres  : %d\n", singleBaseQualThres);
		if(meanSizeInsert>0)
		{
			printf("insert size:       : %.2f\n", meanSizeInsert);
			printf("insert size sdev.  : %.2f\n", standardDev);
		}
		printf("output directory   : %s\n", outputPathStr);
		printf("graph file         : %s\n", graphFile);
		printf("contig file        : %s\n", contigsFileFasta);
		printf("minimal contig size: %d\n\n", minContigLen);
	}else if(operationMode==1)
	{
		printf("\nkmer size          : %d\n", kmerSize);
		//printf("read length cutoff : %d\n", readLenCutOff);
		printf("read length        : %d\n", readLen);
		printf("paired mode        : %d\n", pairedMode);
		for(i=0; i<readFileNumPara; i++)
		printf("read files[%d]      : %s\n", i, readFilesInput[i]);
		printf("single qual thres  : %d\n", singleBaseQualThres);
		printf("output directory   : %s\n", outputPathStr);
		printf("graph file         : %s\n", graphFile);
	}else if(operationMode==2)
	{
		printf("\ngraph file         : %s\n", graphFile);
		printf("kmer size          : %d\n", kmerSize);
		//printf("read length cutoff : %d\n", readLenCutOff);
		printf("read length        : %d\n", readLen);
		printf("paired mode        : %d\n", pairedMode);
		if(meanSizeInsert>0)
		{
			printf("insert size:       : %.2f\n", meanSizeInsert);
			printf("insert size sdev.  : %.2f\n", standardDev);
		}
		printf("output directory   : %s\n", outputPathStr);
		printf("contig file        : %s\n", contigsFileFasta);
		printf("minimal contig size: %d\n\n", minContigLen);
	}else
	{
		printf("Exception: invalid command %d\n", operationMode);
		return FAILED;
	}


	return SUCCESSFUL;
}

/**
 * Set the global input and output path directory.
 *  @return:
 *  	If succeeds, return SUCCESSFUL; otherwise, return FAILED.
 */
short setGlobalPath(const char *outPathStr)
{
	int inPathLen, outPathLen;
	struct stat st;

	strcpy(outputPathStr, outPathStr);

	inPathLen = strlen(inputPathStr);
	if(inPathLen>=255)
	{
		printf("line=%d, In %s(), input path length=%d, error!\n", __LINE__, __func__, inPathLen);
		return FAILED;
	}
	outPathLen = strlen(outputPathStr);
	if(outPathLen>=255)
	{
		printf("line=%d, In %s(), output path length=%d, error!\n", __LINE__, __func__, outPathLen);
		return FAILED;
	}

	if(inputPathStr[inPathLen-1]!='/')
	{
		inputPathStr[inPathLen] = '/';
		inputPathStr[inPathLen+1] = '\0';
	}

	if(outPathLen>0)
	{
		if(outputPathStr[outPathLen-1]!='/')
		{
			outputPathStr[outPathLen] = '/';
			outputPathStr[outPathLen+1] = '\0';
		}

		if(stat(outPathStr, &st)==-1)
		{
			if(mkdir(outPathStr, 0755)==-1)
			{
				printf("line=%d, In %s(), cannot create directory [ %s ], error!\n", __LINE__, __func__, outPathStr);
				return FAILED;
			}
		}
	}else
	{
		strcpy(outputPathStr, "./");
	}

	return SUCCESSFUL;
}

/**
 * Free the global parameters.
 */
void freeGlobalParas()
{
	int i;
	if(operationMode==0 || operationMode==1)
		for(i=0; i<readFileNum; i++) { free(readFilesInput[i]); readFilesInput[i] = NULL; }
	free(kmerSeqInt);
	kmerSeqInt = NULL;
	free(kmerSeqIntRev);
	kmerSeqIntRev = NULL;
}

/**
 *
 */
short getReadsFileFormat(int *readsFileFormatType, char **readFilesInput, int readFileNum)
{
	int i, readFormat;
	FILE *fpRead;
	char ch;

	readFormat = -1;
	for(i=0; i<readFileNum; i++)
	{
		fpRead = fopen(readFilesInput[i], "r");
		if(fpRead==NULL)
		{
			printf("line=%d, In %s(), cannot open file [ %s ], error!\n", __LINE__, __func__, readFilesInput[i]);
			return FAILED;
		}

		ch = fgetc(fpRead);
		if(ch=='>')
		{
			if(readFormat!=-1 && readFormat!=FILE_FORMAT_FASTA)
			{
				printf("Reads files cannot be in multiple format comblined together, error!\n");
				return FAILED;
			}

			readFormat = FILE_FORMAT_FASTA;
		}else if(ch=='@')
		{
			if(readFormat!=-1 && readFormat!=FILE_FORMAT_FASTQ)
			{
				printf("Reads files cannot be in multiple format combined together, error!\n");
				return FAILED;
			}

			readFormat = FILE_FORMAT_FASTQ;
		}else
		{
			printf("Reads files must be in fasta or fastq format, error!\n");
			return FAILED;
		}

		fclose(fpRead);
	}

	if(readFormat!=-1)
		*readsFileFormatType = readFormat;
	else
	{
		*readsFileFormatType = -1;
		printf("Reads files must be in fasta or fastq format, error!\n");
		return FAILED;
	}

	return SUCCESSFUL;
}
