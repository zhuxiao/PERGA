#ifndef METHODS_H_INCLUDED
#define METHODS_H_INCLUDED

#include "structure.h"

//================= srga.c 函数声明 开始 ================/
short parseCommandParasAndExe(int argc, char **argv);
short showUsageInfo();
//================= srga.c 函数声明 结束 ================/

//================= srga.c 函数声明 开始 ================/
short startSRGA(int operationModePara, int kmerSizePara, int readLenCutOffPara, int pairedModePara, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, char *graphFilePara, double meanSizeInsertPara, double standardDevPara, char *outputPathPara, char *outputPrefixPara, int minContigLenPara);
short initGlobalParas(int operationModepara, char *outputPathName, char *prefix, char **readFilesPara, int readFileNumPara, int singleBaseQualThresPara, int pairedModePara, int kmerLen, int readLenCut, char *graphFilePara, double meanSizeInsertPara, double standardDevPara, int minContigLenPara);
short setGlobalPath(const char *outPathStr);
void freeGlobalParas();
short getReadsFileFormat(int *readsFileFormatType, char **readFilesInput, int readFileNum);
short getReadLenAndCountFromFilesFasta(uint64_t *tmpTotalReadNum, int *tmpAverReadLenInFile, int *tmpMaxReadLenInFile, int *tmpMinReadLenInFile, char **readsFileNames, int readsFileNum);
short getReadLenAndCountFromFilesFastq(uint64_t *tmpTotalReadNum, int *tmpAverReadLenInFile, int *tmpMaxReadLenInFile, int *tmpMinReadLenInFile, char **readsFileNames, int readsFileNum);
short getReadLenAndCountFromSingleFasta(int64_t *tmpReadCount, int64_t *tmpSumReadLen, int64_t *tmpMaxReadLen, int64_t *tmpMinReadLen, char *fastaFile);
short getReadLenAndCountFromSingleFastq(int64_t *tmpReadCount, int64_t *tmpSumReadLen, int64_t *tmpMaxReadLen, int64_t *tmpMinReadLen, char *fastqFile);
//================= srga.c 函数声明 结束 ================/


//================= util.c 函数声明 开始 ================/
short outputContig(contigtype *contighead);
short outputContigEnd3(contigtype *contighead, int nodeNum);
short outputContigEnd5(contigtype *contighead, int nodeNum);
short outputUndelKmerpos(graphtype *graph);
short outputRemainedKmers(graphtype *graph);
short releaseGraph(graphtype *graph);
short outputReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputReadsInDecisionTableToFile(assemblingreadtype *decisionTable, int readsNum, int contigID, int contigIndex);
short outputFailedReadsInDecisionTable(assemblingreadtype *decisionTable, int itemNumDecisionTable, int contigID, int contigIndex);
short outputLockedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputMatedReadsInDecisionTable(assemblingreadtype *decisionTable, int readsNum);
short outputKmer(graphtype *graph, int hashcode, uint64_t *kmerSeqInt); //输出kmer中的内容
short outputRidpos(ridpostype *ridpos, int posNum);
void outputSuccessReads(successRead_t *successReadArray, int successReadNum);
short checkGraph(graphtype *graph);
void print_info();
short outputContigToTmpFile(contigtype *contighead, int outFileType);
short outputPEHashArray(PERead_t **PEHashArray);
short checkReadListArr(readList_t *readListArray, int64_t itemNumInReadListArray);
short outputReadListToFile(char *readListFile, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray);
short outputMatedReadsInReadListToFile(char *readListFile, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray);
short convertFragmentSizeFileInText(const char *fragmentSizeFile);
short outputReadPosInGraph(int64_t readID, graphtype *graph);
short outputNaviOccQueue(double *naviOccQueuePara, int itemNumNaviOccQueuePara, int frontRowNaviOccQueuePara);
//================= util.c 函数声明 结束 ================/

//================= graph.c 函数声明 开始 ================/
short constructGraph(char *graphFileName, char **readsFileNames, int readsFileNum);
short constructGraphBySEFasta(char *graphFile, char **readsFileNames, int readsFileNum);
short constructGraphBySEFastq(char *graphFile, char **readsFileNames, int readsFileNum);
short constructGraphByPEFastaSeparate(char *graphFile, char **readsFileNames, int readsFileNum);
short constructGraphByPEFastqSeparate(char *graphFile, char **readsFileNames, int readsFileNum);
short constructGraphByPEFastaInterleaved(char *graphFile, char **readsFileNames, int readsFileNum);
short constructGraphByPEFastqInterleaved(char *graphFile, char **readsFileNames, int readsFileNum);
short getMinReadLenFromFastaFiles(int *readLenInFile, char **readFilesInput, int readFileNum);
short getMinReadLenFromFastqFiles(int *readLenInFile, char **readFilesInput, int readFileNum);
short getReadLenFromFasta(int *tmpReadLen, char *fastqFile);
short getReadLenFromFastq(int *tmpReadLen, char *fastqFile);
graphtype *ReadDataFilePreBySEFasta(char **readsFileNames, int readsFileNum);
graphtype *ReadDataFilePreBySEFastq(char **readsFileNames, int readsFileNum);
graphtype *ReadDataFilePreByPEFastaSeparate(char **readsFileNames, int readsFileNum);
graphtype *ReadDataFilePreByPEFastqSeparate(char **readsFileNames, int readsFileNum);
graphtype *ReadDataFilePreByPEFastaInterleaved(char **readsFileNames, int readsFileNum);
graphtype *ReadDataFilePreByPEFastqInterleaved(char **readsFileNames, int readsFileNum);
short fillReadsToBufFasta(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum);
short fillReadsToBuf(FILE *fpReads, readBuf_t *pBuf, uint64_t *readsNum);
short getSingleReadFasta(FILE *fpPE, readBuf_t *pReadBuf);
short getSingleReadFastq(FILE *fpPE, readBuf_t *pReadBuf);
short containUnknownBase(char *seq);
float calcAverQual3End(char *qual_data);
float calcAverQual5End(char *qual_data);
short qualSatisfied(char *qual_data);
float getRatioBaseA(char *seq);
short addReadPre(char *seq, graphtype *graph);
short generateKmerSeqInt(uint64_t *seqInt, char *seq);
graphtype *ReadDataFileBySEFasta(char **readsFileNames, int readsFileNum, graphtype *graph);
graphtype *ReadDataFileBySEFastq(char **readsFileNames, int readsFileNum, graphtype *graph);
graphtype *ReadDataFileByPEFastaSeparate(char **readsFileNames, int readsFileNum, graphtype *graph);
graphtype *ReadDataFileByPEFastqSeparate(char **readsFileNames, int readsFileNum, graphtype *graph);
graphtype *ReadDataFileByPEFastaInterleaved(char **readsFileNames, int readsFileNum, graphtype *graph);
graphtype *ReadDataFileByPEFastqInterleaved(char **readsFileNames, int readsFileNum, graphtype *graph);
short addRead(char *seq, uint64_t rid, graphtype *graph);
uint64_t kmerhashInt(uint64_t *seqInt);
graphtype *initgraph();
short countKmer(uint64_t hashcode, uint64_t *kmerSeqInt, graphtype *graph);
short addKmer(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph);
short delKmer(char *str, uint64_t rid, unsigned short rpos, graphtype *graph);
short delKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, unsigned short rpos, graphtype *graph);

int findStartIndex(uint64_t rid, ridpostype *rid_pos_table, int posNum);
int getExectIndex(uint64_t rid, uint16_t rpos, int startIndex, ridpostype *ridpostable, int posNum);
int findDirectIndex(uint64_t rid, uint16_t rpos, ridpostype *ridpostable, int posNum);
inline kmertype *getKmerByBase(char *tmp_kmerseq, graphtype *graph);
inline kmertype *getKmer(uint64_t *kmerSeqInt, graphtype *graph);
inline kmertype *getKmerByHash(uint64_t hashvalue, uint64_t *kmerSeqInt, graphtype *graph);
//char *getKmerseqInBase(kmertype *kmer);
short identicalKmerSeq(uint64_t *kmerSeqInt1, uint64_t *kmerSeqInt2);
char *getKmerBaseByInt(uint64_t *kmerSeqInt);
//char *getKmerseqByHash(unsigned int hashvalue);
kmertype *getReverseKmer(uint64_t *kmerSeqIntRev, uint64_t *kmerSeqInt, graphtype *graph);
//kmertype *getReverseKmerByHash(unsigned int hashcode, graphtype *graph);
unsigned int getReverseKmerseqIntByHash(unsigned int hashcode);
void reverse(char *str);
short outputGraphToFile(char *graphFile, graphtype *graph);
short loadGraph(graphtype **graph, char *graphFile);
short GlobalParasFromGraph(int *readLenPara, int *kmerSizePara, uint64_t *hashTableSizePara, int *pairedModePara, char *graphFileName);
short resetGraph(graphtype *graph);
//================= graph.c 函数声明 结束 ================/

//================= contig.c 函数声明 开始 ================/
short buildContigs(char *contigFile);
short initMemory();
void freeMemory();
short initFirstKmerBounder(double *lowerBoundFirstKmer, double *upperBoundFirstKmer, short assemblyCycle, double averKmerOccNum);
short initContig(contigtype **contighead, contigtype ** tailContig);
struct kmertype *getFirstKmer(graphtype *graph, unsigned int *kmerIndex);
short initFirstKmerThreshold();
short getFirstKmers(uint64_t *kmerIndex, kmertype **firstKmer);
short isValidFirstKmer(kmertype *kmer);
short containFirstPos(kmertype *kmer);
short containLastPos(kmertype *kmer);
short initAssemblingReadTable(assemblingreadtype *assemblingreads);
short getNextKmerBySE(int contigNodesNum);
int addContigBase(contigtype **contigtail, unsigned char base, int contigIndex);
short addRidposToContig(successRead_t *successReadArray, int *successReadNum, contigtype *contig36, int contigNodesNum);
void releaseContig(contigtype *contighead);
short outputContigToFile(FILE *fpContig, int outFileType, contigtype *contighead, int contigID, int contigNodeNum);
short addFirstKmerToDecisionTable(kmertype **kmers);
short addReadToDecisionTable(uint64_t rid, int rpos, int orientation, int matedFlag);
short reallocateDecisionTable();
//int initAssemblingTable(kmertype *kmers[2], assemblingreadtype *assemblingreads, int numassemblingreads);
short computeLongKmerOccNum(kmertype *tmp_kmers[2], int *occNum, int length_k);
short computeKmerOccNum(kmertype *tmp_kmers[2], int *occNum);
short computeKmerOccNumUnlocked(kmertype *tmp_kmers[2], int *occNum);
short computeKmerOccNumLocked(kmertype *tmp_kmers[2], int *occNum);
short delReadsFromGraph(successRead_t *successReadArray, int successReadNum, char *sequence);
short reverseReadseq(char *seq);
short delRemainedKmers(char *seq, uint64_t *tmp_kmerseq, uint64_t rid, uint16_t rpos, graphtype *graph);
short updateLockedReads();
short updateContigtailnodes(contigtype *contig36, contigtype *successContig, int *contigIndex);
void updateContigheadnodes(contigtype **contighead, int *contigNodeNum);
int getProperIndex(assemblingreadtype *assemblingread,assemblingreadtype *assemblingreads,int numassemblingreads);
int getProperIndexLimited(assemblingreadtype *assemblingread,assemblingreadtype *assemblingreads,int numassemblingreads, int limitLastpos);

contigtype *getLastContig(contigtype *contighead, int contigNodesNum, int lastNodesNum);
short getLastseq(char *lastseq, contigtype *startContig);
short getSecondAssemblyFirstKmers(contigtype *contighead, graphtype *graph);
contigtype *getFirstSuccessContig(contigtype *contighead);
short trimContigBeforeCycle2(contigtype **contighead, contigtype **contigtail, contigtype **successContig, int *contigNodeNum);
short reverseContig(contigtype **contighead, contigtype **contigtail);
int initAssemblingTableSecondAssembly(char *lastseq36, graphtype *graph);
short resetContigIndex(contigtype *contighead);
contigtype *getSuccessContig(contigtype *contighead, contigtype *successContig, int contigIndex);

short getReversedSeq(char *reversed_seq, char *seq, int seq_len);
short recoverRemainedKmers(char *seq, uint64_t *tmp_kmerseq, uint64_t rid, uint16_t rpos, graphtype *graph);
short recoverKmerByHash(uint64_t hashcode, uint64_t *kmerSeqInt, uint64_t rid, uint16_t rpos, graphtype *graph);
short recoverReadFromGraph(char *seq, uint64_t rid, graphtype *graph);
short recoverDeledReads( contigtype *startContig);

contigtype *getContigtail(contigtype *startContig);
short initSecondAssembly();
short updateReadsNumReg(int itemNumSuccessReadsArr, int contigNodesNum, int assemblyRound);
short initReadsNumRegSecondAssembly(int contigNodesNum);
//================= contig.c 函数声明 结束 ================/

//================= update.c 函数声明 开始 ================/
short updateDecisionTable(kmertype *tmp_kmers[2]);
ridpostype *getRidpos(assemblingreadtype assemblingRead, ridpostype *rid_pos_table, int posNum);
short reallocateSuccessReadsArr();
short removeFinishedReadsFromDecisionTable();
short updateFinishedReadsInDecisionTable();
short updateAssemblingreadsStatus(); //更新决策表中reads状态标记
short updateSameReadStatusToFailed(uint64_t rid, int assemblingreadIndex);//当决策表中有成功的reads的时候,更新决策表中相同rid的reads的状态标记为失败标记
//================= update.c 函数声明 结束 ================/

//================= hashPE.c 函数声明 开始 ================/
short estimateInsertSizeAndSdev();
short initPEHashParas();
short updatePEHashTable(int contigNodesNum, int assemblyRound);
short getReadFromPEHashtable(PERead_t **pRead, uint64_t readID);
short addReadToPEHashtable(successRead_t *ridposOrient, int contigPos, int assemblyRound);
short delReadfromPEHashtable(uint64_t readID);
short cleanReadsFromPEHashtable();
short initPEHashtableSecondAssembly(contigtype *contighead, int contigNodesNum);

short meanSizeInsertAndSdevEstimation(const char *fragmentSizeFile, estContig_t *estContigArray, int contigNumEstContigArray);
short getPairedEndsFromSingleContig(FILE *fpFragSize, estContig_t *estContig);
short initMemGetPESingleContig(estContig_t *estContig);
void freeMemGetPESingleContig();
short getTotalReadsNumOfSingleContig(int64_t *totalReadsNum, contigtype *contighead);
short fillDataReadPosTmpArr(readPosTemp_t *readPosTmpArray, contigtype *contighead);
short radixSortReadPosTmpArr(readPosTemp_t *readPosTmpArray, readPosTemp_t *readPosTmpArrBuf, int64_t readsNumInArray);
short fillDataReadList(readList_t *readListArray, readPos_t *readPosArray, int64_t *itemNumInReadListArray, int64_t *itemNumInReadPosArray, readPosTemp_t *readPosTmpArray, int64_t itemNumInReadPosTmpArray);
short outputInsertSizeToFile(FILE *fpFragSize, readList_t *readListArray, readPos_t *readPosArray, int64_t itemNumInReadListArray);
short computeInsertSizeAndSdev(double *meanSizeInsert, double *standardDev, const char *fragmentSizeFile);
//================= hashPE.c 函数声明 结束 ================/

//================= PEAssembly.c 函数声明 开始 ================/
short buildEstContigs(char *contigFile);
short getNextKmerByMix(int contigNodesNum, int assemblyRound);
short getNextKmerByPE(int contigNodesNum, int assemblyRound);
short computeKmerOccNumByPE(kmertype *tmp_kmers[2], int *occNum, int contigNodesNum, int assemblyRound);
short computeKmerOccNumUnlockedByPE(kmertype *tmp_kmers[2], int *occNum, int contigNodesNum, int assemblyRound);
short computeKmerOccNumLockedByPE(kmertype *tmp_kmers[2], int *occNum);
short computeLongKmerOccNumByPE(kmertype *tmp_kmers[2], int *occNum, int length_k, int contigNodesNum, int assemblyRound);
short validReadPair(uint64_t readID);
short trimContigTailByReadLen(contigtype *contighead, contigtype **contigtail, contigtype **successContig, int *contigNodesNum, int assemblyRound);
short setEmptyNaviOccQueue(double *naviOccQueuePara, int *itemNumNaviOccQueuePara, int *frontRowNaviOccQueuePara, int *rearRowNaviOccQueuePara);
short updateNaviOccQueue(double *naviOccQueuePara, int *itemNumNaviOccQueuePara, int *frontRowNaviOccQueuePara, int *rearRowNaviOccQueuePara, double maxOccNum);
short calcAverOccNaviOccQueue(double *averOccNum, double *naviOccQueuePara, int itemNumNaviOccQueuePara);
short getLowOccLenNaviOccQueue(int *lowLen, double *naviOccQueuePara, int itemNumNaviOccQueuePara, int frontRowNaviOccQueuePara);
short computeGapSizeInContig(int *gapSize, contigtype *contighead, contigtype *contig36, int contigNodesNum, int assemblyRound);
//================= PEAssembly.c 函数声明 结束 ================/

#endif // METHODS_H_INCLUDED
