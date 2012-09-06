/*
 * inc/extvab.h
 *
 *  Created on: May 16, 2011
 *      Author: zhuxiao
 */

#ifndef EXTVAB_H_
#define EXTVAB_H_

#include "constants.h"


//=================== Global files or paths ==================
extern char inputPathStr[256];
extern char outputPathStr[256];
extern char outputPrefixStr[256];

extern char *readFilesInput[256];
//extern char readFilesFastqPE[2][256];				// [0]-- xxx_1.fastq, [1]-- xxx_2.fastq
extern int readFileNum;
extern int readsFileFormatType;
extern char contigsFileFasta[256];
extern char contigsFileHanging[256];

extern char fragmentSizeFile[256];
extern char graphFile[256];
extern char sampleContigsFile[256];

extern int operationMode;
extern int kmerSize;
int readLenCutOff;					// the read length after cutting at the 3' end of reads
extern int pairedMode;
extern double meanSizeInsert, standardDev;
extern int minContigLen;


//**************** for graph.c *******************
extern graphtype *deBruijnGraph;

extern double timeuse_deBruijn;
//extern long int countMemory;
extern long newKmerNum;
extern long newRidposNum;
extern long totalReadNum, validReadNum;  //reads总数, 有效reads总数, 未拼接的reads总数

//**************** for contig.c *******************
extern int readLenInFile;
extern int readLenCutOff;

extern int qualityBaseNumEnd3;
extern int qualityBaseNumEnd5;
extern int errorRegLenEnd3;
extern int singleBaseQualThres;			// single base quality threshold

extern int longKmerSize;
extern int longKmerStepSize;
extern double averKmerOcc;
extern double firstKmerThres;
extern double minKmerOccSE;
extern double minKmerOccPE;
extern double maxSecondOcc;
extern double maxFirstOcc;
extern double minLongKmerOcc;
extern double lockedReadsNumThres;
extern double minReadsNumPEHashThres;
extern double maxOccNumFaiedPE;
extern int navigationFlag;

extern double *naviOccQueue;
extern int itemNumNaviOccQueue;
extern int maxItemNumNaviOccQueue;
extern int frontRowNaviOccQueue;
extern int rearRowNaviOccQueue;
extern double lowOccThresNaviOccQueue;

extern int kmer_len;
extern int occsNumSE[4], occsNumPE[4];  //下一个kmer的得分，对应ACGT的顺序

extern int maxOccIndexSE, maxOccIndexPE, secondOccIndexSE, secondOccIndexPE;
extern double maxOccSE, maxOccPE, secondOccSE, secondOccPE;

extern short assemblyRound; //FRIST_ROUND_ASSEMBLY  or SECOND_ROUND_ASSEMBLY
extern int lockedReadsNum;
extern kmertype *kmers[2]; //kmers[0]存放正向的kmer, kmers[1]存放反向互补的kmer

extern assemblingreadtype *decisionTable;		// the decision table
extern int itemNumDecisionTable;				// item number in decision table
extern int maxItemNumDecisionTable;				// the maximal number of reads in decision table

extern contigtype *contighead;  //contig头指针
extern contigtype *contigtail;  //contig链表中的最后一个contig的指针
extern contigtype *contig36;  //指向当前contig的前36个节点的指针
extern int contigIndex;
extern int contigsNum; //contig的数量
extern int64_t basesNum; //拼接的所有contigs的碱基的总和
extern int64_t this_successReadNum;
extern uint64_t kmerIndex;  //标识拼接初始kmer位于de Bruijn图中的哈希表中的位置
extern contigtype *successContig; //指向具有成功reads的contig节点
extern char *lastseq36;  //暂存contig的最后36个碱基的数组

extern int itemNumSuccessReadsArr;
extern int maxItemNumSuccessReadsArr;
extern successRead_t *successReadsArr;

extern FILE *fpContigsBase, *fpContigsHanging;
extern int hangingContigOutFlag;		// whether output the hanging contig file: YES / NO (default)

extern int number_of_overlap_less_than_threshold;
extern int64_t successReadNum;


//**************** for update.c *******************
int number_of_corrected_reads;

//++++++++++++++++++++++++++++++++++++
extern int PEGivenType, oldPEGivenType;
extern int estimateSuccessFlag;;
extern PERead_t **PEHashArr;
extern int readsNumInPEHashArr;
extern int allowedUpdatePEHashArrFlag;
extern int regLenPEHash, maxRegLenPEHash, minRegLenUsingPE, minMarginLenPEHash, maxMarginLenPEHash;// int leftMarginLenPEHash;
extern int minContigLenUsingPE, shiftLenRound1, validReadOrientPEHash;
extern contigtype *hashRegLeftContig, *hashRegRightContig, *shiftedRegLeftContig, *shiftedRegRightContig;

extern double oldMeanSizeInsert, oldStandardDev, standardDevFactor, meanSizeInsertEst, standardDevEst;
extern estContig_t *estContigArr;
extern int contigNumEstContigArr;
extern int minContigLenEst;			// the minimal contig length that can be used to estimate the insert size and standard deviation

extern readPosTemp_t *readPosTmpArr, *readPosTmpArrBuf;
extern int64_t readsNumSingleContig, totalReadsNumContigs;

extern readList_t *readListArr;
extern readPos_t *readPosArr;
extern int64_t itemNumInReadListArr, itemNumInReadPosArr;



//++++++++++++++++++++++++++++++++++++
extern int readLen;
//extern int kmerSize;
extern int entriesPerKmer;
extern uint64_t lastEntryMask;
extern int lastEntryBaseNum;		// the base number of the last entry of its sequence array

extern uint64_t hashTableSize;

extern char kmerBaseSeq[1000];

extern uint64_t *kmerSeqInt;
extern uint64_t *kmerSeqIntRev;

extern uint64_t *kmerSeqIntAssembly;
extern uint64_t *kmerSeqIntAssemblyRev;
extern uint64_t *tmpKmerSeqIntAssembly;

extern kmertype *firstKmer;
extern uint64_t tabooSeqInt[4];
extern int navigationID;		// 0 -- Single End, 1 -- Paired End, -1 -- unknown
extern int navigationNumSE;
extern int maxNavigationNumSE;


//++++++++++++++++++++++++++++++++++++
extern int regLenReadsNumReg, maxRegLenReadsNumReg;
extern int minContigLenCheckingReadsNum;
extern contigtype *leftContigReadsNumReg, *rightContigReadsNumReg;
extern int readsNumTotal, readsNumReadsNumReg;
extern double readsNumRatio;
extern double maxReadsNumRatioThres, minReadsNumRatioThres;
extern int solvedRepeatsNum;


#endif /* EXTVAB_H_ */
