#ifndef CDEM_HEADER_H
#define CDEM_HEADER_H

#include <string>
#include <algorithm>
#include <fstream>
#include <queue>
#include <functional>

#define NO_DATA_VALUE -9999.0f


/*
*	reverse of flow directions
*	2	4	8
*	1	0	16
*	128	64	32
*/
static unsigned char inverse[8] = {16, 32, 64, 128, 1, 2, 4, 8};
/*
*	flow direction		
*	32	64	128		
*	16	0	1		
*	8	4	2		
*/
static unsigned char	dir[8] = {1, 2, 4, 8, 16, 32, 64, 128};
class CDEM
{
protected:
	float* pDem;
	float* borderInNeighbor1, *borderInNeighbor2;
	unsigned int width, height;
public:
	CDEM()
	{
		pDem=NULL;
		borderInNeighbor1=NULL;
		borderInNeighbor2=NULL;

	}
	~CDEM()
	{
		delete[] pDem;
		delete[] borderInNeighbor1;
		delete[] borderInNeighbor2;
	}
	bool Allocate();

	void freeMem();

	void initialElementsNodata();
	float asFloat(unsigned int row,unsigned int col) const;
	void Set_Value(unsigned int row,unsigned int col, float z);
	bool is_NoData(unsigned int row, unsigned int col) const;
	void Assign_NoData();
	int Get_NY() const;
	int Get_NX() const;
	float* getDEMdata() const;
	float* getBorderInNeighbor1();
	float* getBorderInNeighbor2(); 
	void SetHeight(unsigned int height);
	void SetWidth(unsigned int width);
	void readDEM(const std::string& filePath);
	bool is_InGrid(unsigned int row,unsigned  int col) const;
	float getLength(unsigned int dir);
	unsigned char getDirction(unsigned int row, unsigned int col, float spill);
};
#endif
