#ifndef CONSTANTS_H_INCLUDED
#define CONSTANTS_H_INCLUDED 1

#define DEBUG_OUTPUT 0
#define DEBUG_PRINTF 0

#define __USE_FILE_OFFSET64
#define __USE_LARGEFILE64
//#define _FILE_OFFSET_BITS 64

#define VERSION_STR						("v0.3.04.05")
#define RELEASE_DATE_STR				("Jul 29, 2012")

#ifndef NULL
#define NULL ((void *)0)
#endif

#define SUCCESSFUL 0
#define FAILED	1
#define ERROR	-1

#define YES 1
#define NO 0

/* the Strand flag */
#define ORIENTATION_PLUS				'+'
#define ORIENTATION_MINUS				'-'

#define DEFAULT_KMER_SIZE				21
#define LONG_KMER_SIZE_FACTOR			0.7f

#define HASH_TABLE_SIZE					(157286611LLU)

#define TABLE_SIZE_ASSEMBLINGREAD		50000
#define TABLE_SIZE_RIDPOSORIENTATION 	50000
#define MAX_DECISION_TABLE_SIZE_HTRES	50000

#define FIRST_ROUND_ASSEMBLY			1		// the first assembly round
#define SECOND_ROUND_ASSEMBLY			2		// the second assembly round

#define FIRSTKMER_FACTOR 				5

#define FIRSTKMER_SUBTRACT_THRESHOLD 	0.9f

#define AVERAGE_QUAL_THRESHOLD_3End 	2.0f
#define AVERAGE_QUAL_THRESHOLD_5End 	20.0f
#define SINGLE_QUAL_THRESHOLD			2

#define ARTIFACTS_BASE_A_THRESHOLD 		0.9f

#define QUAL_BASE_NUM_3End_FACTOR		0.2f
#define ERROR_REGION_LEN_3End_FACTOR	0.1f

//#define VALID_OCC_RATIO					0.05   //==================================
//#define MAX_SECOND_OCC_FACTOR			2  //==================================
#define MAX_SECOND_OCC_FACTOR			4.0f
//#define MAX_SECOND_OCC_THRES			8  //==================================
#define MAX_SECOND_OCC_THRES			6
//#define MAX_FIRST_OCC_FACTOR			30  //==================================
#define MAX_FIRST_OCC_FACTOR			20

//#define SECOND_FIRST_SECORE_RATIO		0.5f   //==================================
//#define SECOND_FIRST_SECORE_RATIO		0.7f
//#define SECOND_FIRST_OCC_RATIO			0.5f   //==================================
//#define SECOND_FIRST_OCC_RATIO			0.7f // best
#define SECOND_FIRST_OCC_RATIO			0.8f

#define CONTIG_LEN_THRESHOLD			100


//#define MIN_OVERLAP_LEN					23   //==================================
//#define MIN_OVERLAP_LEN					11
#define MIN_OVERLAP_LEN					7

//#define MIN_KMER_OCC_FACTOR			0.2f  // best : 0.17f
//#define MIN_KMER_OCC_FACTOR			0.17f    //==================================
#define MIN_KMER_OCC_FACTOR				0.15f
//#define MIN_KMER_OCC_THRES				3    //==================================
#define MIN_KMER_OCC_THRES				2

//#define MIN_LONG_KMER_OCC_THRES			8  //==================================
#define MIN_LONG_KMER_OCC_THRES			6
//#define LONG_KMER_OCC_FACTOR			2   //==================================
#define LONG_KMER_OCC_FACTOR			2.5f

//决策表中reads的3种状态
#define ASSEMBLING_STATUS				1		//拼接中状态
#define SUCCESSFUL_STATUS				2		//拼接成功状态
#define FAILED_STATUS					3		//拼接失败状态

#define BASE_TYPE_FASTA_CONTIG_FILE				0
#define HANGING_READ_TYPE_CONTIG_FILE			1

//++++++++++++++++++++++++++++++++++++
#define MAX_READ_BUF_SIZE				1000000


#define FILE_FORMAT_FASTA				1
#define FILE_FORMAT_FASTQ				2

// PE given types: 0--singeEnd, 1--none, 2--only insert, 3-both insert and sdev.
#define SE_GIVEN_TYPE					0
#define NONE_PE_GIVEN_TYPE				1
#define INSERT_PE_GIVEN_TYPE			2
#define BOTH_PE_GIVEN_TYPE				3

#define RID_LOW_BITS_NUM				20
#define TABLE_SIZE_HASH_PE				(1 << RID_LOW_BITS_NUM)
#define RID_LOW_BITS_MASK				((1 << RID_LOW_BITS_NUM) - 1)

#define MIN_READ_NUM_PE_HASH_FACTOR		0.5f

#define MAX_NUM_EST_CONTIG				100
#define TOTAL_CONTIG_LEN_EST_THRES		500000
#define MIN_CONTIG_LEN_EST				1000
#define MIN_CONTIG_LEN_EST_FACTOR		5.0f

#define MAX_INSERT_SIZE_FACTOR			5.0f
#define SDEV_FACTOR						3.0f  //===================================
//#define SDEV_FACTOR						5.0f
#define DRAFT_SDEV_FACTOR				0.15f
#define REG_LEN_PE_HASH_FACTOR			0.6f


//++++++++++++++++++++++++++++++++++++
//#define REG_LEN_READS_NUM_REG_FACTOR			1.3f
#define REG_LEN_READS_NUM_REG_FACTOR			2.5f
//#define MAX_READS_NUM_RATIO_THRES				2.0f   //===================================
//#define MAX_READS_NUM_RATIO_THRES				5.0f
#define MAX_READS_NUM_RATIO_THRES				10.0f
#define MIN_READS_NUM_RATIO_THRES				0.1f   //===================================
//#define MIN_READS_NUM_RATIO_THRES				0.2f
//#define MIN_READS_NUM_RATIO_THRES				0.05f

//#define OCCS_NUM_SE_FAILED_PE_FACTOR			1.0f  //===================================
//#define OCCS_NUM_SE_FAILED_PE_FACTOR			5.0f
#define OCCS_NUM_SE_FAILED_PE_FACTOR			6.0f
//#define OCCS_NUM_SE_FAILED_PE_FACTOR			4.5f
//#define MAX_OCC_NUM_FAILED_PE_THRES				51.0f
#define MAX_OCC_NUM_FAILED_PE_THRES				60.0f

#define MAX_NAVI_NUM_SE_THRES					2  // === not used ===

#endif // CONSTANTS_H_INCLUDED
