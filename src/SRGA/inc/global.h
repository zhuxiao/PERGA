/*
 * inc/global.h
 *
 *  Created on: May 9, 2011
 *      Author: zhuxiao
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_ 1




//=================== Global files or paths ==================
char inputPathStr[256];
char outputPathStr[256];
char outputPrefixStr[256];

char *readFilesInput[256];
//char readFilesFastqPE[256][256];				// [0]-- xxx_1.fastq, [1]-- xxx_2.fastq
int readFileNum;
int readsFileFormatType;
char contigsFileFasta[256];
char contigsFileHanging[256];

char fragmentSizeFile[256];
char graphFile[256];
char sampleContigsFile[256];

int operationMode;
int kmerSize;
int readLenCutOff;					// the read length after cutting at the 3' end of reads
int pairedMode;
double meanSizeInsert, standardDev;
int minContigLen;
int trimReadLenFlag;


//*********************** graph.c **************************
graphtype *deBruijnGraph;

uint64_t totalKmerNum;
uint64_t totalRidposNum;

uint64_t totalReadNum, validReadNum;  //reads总数, 有效reads总数
double timeuse_deBruijn;


//*********************** contig.c **************************
int occsNumSE[4], occsNumPE[4];  //下一个kmer的得分，对应ACGT的顺序

int maxOccIndexSE, maxOccIndexPE, secondOccIndexSE, secondOccIndexPE;
double maxOccSE, maxOccPE, secondOccSE, secondOccPE;

//int readLenInFile;					// the read length in the data file
int maxReadLenInFile;				// the maximal length of read length in file
int minReadLenInFile;				// the minimal length of read length in file
int averReadLenInFile;				// the average length of read length in file
//int readLenCutOff;					// the read length after cutting at the 3' end of reads

int qualityBaseNumEnd3;				// checked base number at 3' end of reads
int qualityBaseNumEnd5;				// checked base number at 5' end of reads
int errorRegLenEnd3;				// erroneous region length of 3' end of a read
int singleBaseQualThres;			// single base quality threshold

int longKmerSize;
int longKmerStepSize;
double averKmerOcc;
double firstKmerThres;
double minKmerOccSE;
double minKmerOccPE;
double maxSecondOcc;
double maxFirstOcc;
double minLongKmerOcc;
double lockedReadsNumThres;
double minReadsNumPEHashThres;
double maxOccNumFaiedPE;
int navigationFlag;

double *naviOccQueue;
int itemNumNaviOccQueue;
int maxItemNumNaviOccQueue;
int frontRowNaviOccQueue;
int rearRowNaviOccQueue;
double lowOccThresNaviOccQueue;

//unsigned int thiskmerseq = 0;

FILE *fpContigsBase, *fpContigsHanging;
int hangingContigOutFlag;		// whether output the hanging contig file: YES / NO (default)

short assemblyCycle;	// values: 1 for 0.5x <= firstKmerThres <= 15x, 2 for firstKmerThres > 15x, 3 for 2 <= firstKmerThres < 0.5x
short assemblyRound; //FRIST_ROUND_ASSEMBLY  or SECOND_ROUND_ASSEMBLY
kmertype *kmers[2]; //kmers[0]存放正向的kmer, kmers[1]存放反向互补的kmer
int kmer_len;
double lowerBoundFirstKmer, upperBoundFirstKmer;

assemblingreadtype *decisionTable;	// the decision table
int itemNumDecisionTable;			// item number in decision table
int maxItemNumDecisionTable;		// the maximal number of reads in decision table

contigtype *contighead;				//contig头指针
contigtype *contigtail;				//contig链表中的最后一个contig的指针
contigtype *contig36;				//指向当前contig的前36个节点的指针
int contigIndex;
int contigsNum;						//contig的数量
int64_t basesNum = 0;				//拼接的所有contigs的碱基的总和
int64_t this_successReadNum;
uint64_t kmerIndex;					//标识拼接初始kmer位于de Bruijn图中的哈希表中的位置
contigtype *successContig;			//指向具有成功reads的contig节点
char *lastseq36;					//暂存contig的最后36个碱基的数组


int itemNumSuccessReadsArr;
int maxItemNumSuccessReadsArr;
successRead_t *successReadsArr;


int lockedReadsNum = 0;

int64_t successReadNum = 0; //成功拼接的reads总数
int number_of_corrected_reads = 0;
int number_of_overlap_less_than_threshold = 0;


//++++++++++++++++++++++++++++++++++++
int PEGivenType, oldPEGivenType;
int estimateSuccessFlag;
PERead_t **PEHashArr;
int readsNumInPEHashArr;
int allowedUpdatePEHashArrFlag;		// YES or NO
int regLenPEHash, maxRegLenPEHash, minRegLenUsingPE, minMarginLenPEHash, maxMarginLenPEHash;// int leftMarginLenPEHash;
int minContigLenUsingPE, shiftLenRound1, validReadOrientPEHash;
contigtype *hashRegLeftContig, *hashRegRightContig, *shiftedRegLeftContig, *shiftedRegRightContig;

double oldMeanSizeInsert, oldStandardDev, standardDevFactor, meanSizeInsertEst, standardDevEst;
estContig_t *estContigArr;
int contigNumEstContigArr;
int minContigLenEst;			// the minimal contig length that can be used to estimate the insert size and standard deviation

readPosTemp_t *readPosTmpArr, *readPosTmpArrBuf;
int64_t readsNumSingleContig, totalReadsNumContigs;

readList_t *readListArr;
readPos_t *readPosArr;
int64_t itemNumInReadListArr, itemNumInReadPosArr;



//++++++++++++++++++++++++++++++++++++
int readLen;
//int kmerSize;
int entriesPerKmer;
uint64_t lastEntryMask;
int lastEntryBaseNum;		// the base number of the last entry of its sequence array

uint64_t hashTableSize;

char kmerBaseSeq[1000];

uint64_t *kmerSeqInt;
uint64_t *kmerSeqIntRev;

uint64_t *kmerSeqIntAssembly;
uint64_t *kmerSeqIntAssemblyRev;
uint64_t *tmpKmerSeqIntAssembly;

kmertype *firstKmer;
uint64_t tabooSeqInt[4];

int navigationID;			// 0 -- Single End, 1 -- Paired End, -1 -- unknown
int navigationNumSE;
int maxNavigationNumSE;


//++++++++++++++++++++++++++++++++++++
int regLenReadsNumReg, maxRegLenReadsNumReg;
int minContigLenCheckingReadsNum;
contigtype *leftContigReadsNumReg, *rightContigReadsNumReg;
int readsNumTotal, readsNumReadsNumReg;
double readsNumRatio;
double maxReadsNumRatioThres, minReadsNumRatioThres;
int solvedRepeatsNum;


#endif /* GLOBAL_H_ */
