#ifndef STRUCTURE_H_INCLUDED
#define STRUCTURE_H_INCLUDED 1

//kmer结构
typedef struct kmerNode
{
	struct ridposnode *ppos;			//8-bytes; ָ指向kmer位置结构数组到指针
	uint32_t multiplicity;			//4-bytes; kmer重复度
	uint32_t arraysize;				//4-bytes; 位置数组的大小
	uint64_t *kmerseq;
	struct kmerNode *next;
}kmertype;

//de Bruijn图的结构
typedef struct graphNode
{
	struct kmerNode **pkmers;	//指向哈希表数组的指针
}graphtype;


//ridpos结构
typedef struct ridposnode
{
	uint64_t delsign:1;		//删除标记
	uint64_t pos:16;			//位置
	uint64_t reserved:1;	//保留位
	uint64_t rid:46;		//read ID
}ridpostype;

//contig结构, 共24字节
typedef struct contignode
{
	int index;											//4-bytes; contig ID
	unsigned char base;									//1-byte;  contig碱基, 整数
	unsigned short ridposnum;							//2-bytes; rid_pos_orientation数组中元素的个数
	unsigned char reserved;								//1-byte;  保留1字节
	struct successReadNode *pridposorientation;			//8-bytes; 指向ridposorientation数组到指针
	struct contignode *next;							//8-bytes; ָ指向下一个contig节点到指针
}contigtype;

//successReadNode结构, 共16字节
typedef struct successReadNode
{
	uint64_t rid;				//4-bytes; read ID
	uint16_t pos;				//1-byte;  read中的首个pos
	uint16_t matchnum;			//1-byte;  匹配的碱基个数
	uint16_t orientation;		//1-byte;  read方向
}successRead_t;

//assemblingread结构
typedef struct assemblingreadnode
{
	uint64_t rid:48;						// 4-bytes; read ID
	uint64_t firstpos:16;					// 1-bytes; 首个pos位置
	uint16_t kmerappeartimes;			// 1-bytes; kmer出现次数
	uint16_t kmerunappeartimes;		// 1-bytes; kmer未出现次数
	uint16_t lastappearpos;			// 1-bytes; 最近出现的pos位置
	uint16_t lastpos;					// 1-bytes; 上次出现kmer的pos位置
	uint16_t orientation:8;			// 1-bytes; readƴ方向
	uint16_t status: 2;				// 2-bits; 拼接状态:拼接中，成功，失败
	uint16_t kmerunappearblocks: 2;	// 2-bits; kmer连续未出现块的个数
	uint16_t delsign: 1;				// 1-bit; 删除标记
	uint16_t locked: 1;				// 1-bit; 1位锁定标记
	uint16_t reserved: 1;				// 1-bit; 1位保留
	uint16_t matedFlag: 1;				// the mated flag
}assemblingreadtype;


//++++++++++++++++++++++++++++++++++++
typedef struct PEReadNode
{
	int64_t rid;
	int32_t cpos;
	int32_t orient;
	struct PEReadNode *next;
}PERead_t;

typedef struct readBufNode{
	char *seq;
	char *qual;
	int len;
}readBuf_t;

typedef struct estContigNode
{
	uint32_t contigID;
	uint32_t contigLen;
	struct contignode *contighead;
}estContig_t;

// temporary match information of read
typedef struct readPosTempNode
{
	uint64_t readID;
	uint32_t contigPos;
	uint16_t orientation;
	uint16_t matchBaseNum;
}readPosTemp_t;

//======================== structures for read list (RL) ===================
//read list Index
typedef struct readListNode
{
	uint64_t readID;
	uint64_t matchNum: 16;
	uint64_t firstRow: 48;
}readList_t;

//read position array
typedef struct readPosNode
{
	uint32_t contigPos;			// the contig position
	uint16_t orientation;		// the read orientation
	uint16_t matchBaseNum;		// the matched base number
}readPos_t;





#endif // STRUCTURE_H_INCLUDED
