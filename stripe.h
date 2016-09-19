#ifndef Parallele_Filling_H
#define Parallele_Filling_H
#include <iostream>
#include <string>
#include <fstream>
#include <queue>
#include <algorithm>
#include "dem.h"
#include "node.h"
#include "utils.h"
#include <time.h>
#include <list>
#include <stack>
#include <unordered_map>
using namespace std;
typedef std::vector<Node> NodeVector;
typedef std::priority_queue<Node, NodeVector, Node::Greater> PriorityQueue;
class BorderNode
{
public:
	float spill; //current spill elevation value of the border cell
	int isDepressionCell;
	int isPushedToPQAlready;
public:
	BorderNode(){
 	   isDepressionCell=false;
	   isPushedToPQAlready=false;
	}
};
class StripeBorder
{
public:
	std::vector<int> arrSpillChangedNodeIndex; //index in  arrBorderNodes, marking those cells whose spill has changed
	std::vector<BorderNode> arrBorderNodes; // border nodes array for stroring intermediate spill elevations, 
public:
	void Init(int size){
		arrBorderNodes.resize(size);
		arrSpillChangedNodeIndex.resize(size);
		arrSpillChangedNodeIndex.clear(); //clear the array
	}
	void CopyFrom(StripeBorder& stripBorder){
		arrBorderNodes=stripBorder.arrBorderNodes;
		arrSpillChangedNodeIndex=stripBorder.arrSpillChangedNodeIndex;
	}
};
class Stripe
{
public:
	int totalWidth,totalHeight; //the width and height of the original entire DEM
	int isRowMajor;
	int rowOffset, colOffset;
	int width,height; //the width and height of this stripe
	int stripeIndex,stripeSize,stripeNumber; //stripeSize is the standard strip size
	PriorityQueue priorityQueue;
	queue<Node> traceQueue;
	queue<Node> depressionQue;
	Flag mainFlag;//flag for slope cells and depression cells with known spill elevations
	Flag depressionFlag; //flag array for depression cells with uncertain spill elevations
	CDEM* pDemStripe;
	StripeBorder border1, border2, borderInNeighborStrip1, borderInNeighborStrip2;
public:


	void ComputeOffsets();
	float asFloat(int row, int col);//get elevation using relative index in the strip
	void SetData(int row, int col,float val);
	int is_NoData(int row, int col);
	int is_InDEM(int row, int col);//decide whether a cell in this strip is in DEM or not, not in the strip 
	int is_InStrip(int row, int col);//decide whether a cell in this strip 
	int isCellOnInterStripBorder(int row, int col,int& cellIndex,int& borderIndex,float& spill);
	
	void Initialize();
	void PriorityFlood(int isFirstTime);
	void AddCellToBorderDCList(int cellIndex, int borderIndex,float spill);

	void FillDepression(Node node);
	void PushBorderCellsIntoPQ();
	void ProcessSlopeCells(int isFirstTime);
	void ProcessDepressionCells();
	int PushCellIntoPQ_FirstTime(StripeBorder& border, StripeBorder& borderInNeighboringStrip, int borderIndex);
	int UpdatePotentialSpillandPushToPQ_FirstTime(Stripe* pStrip);
	
	int UpdatePotentialSpillandPushToPQ_OtherTimes(Stripe* pStrip);
	void AddSlopeCellToListWithLowerSpillEle(int cellIndex, int borderIndex,float spill);
	int PushCellIntoPQ_OtherTimes(StripeBorder& border, StripeBorder& borderInNeighboringStrip, int borderIndex);
	void FillDepressionFromStripBorder();

	void Share(int isFirstTime);//
	void PackBorder(int isFirstTime,StripeBorder& border, unsigned char*& pBuff, int& byteCount);
	void UnpackBorder(int isFirstTime,StripeBorder& border, unsigned char* pBuff);
};

#endif

