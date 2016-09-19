#ifndef UTILS_HEAD_H
#define UTILS_HEAD_H

#include "gdal_priv.h"
#include <queue>
#include <algorithm>
#include "dem.h"

class Stripe;
#ifdef WIN32
#pragma comment(lib,"gdal_i.lib")
#endif


void calculateStatistics(const CDEM& dem, double* min, double* max, double* mean, double* stdDev);

extern int	ix[8];
extern int	iy[8];
inline int Get_rowTo(int dir, int row)
{
	return( row + ix[dir] );
}
inline int Get_colTo(int dir, int col){
	return( col + iy[dir] );
}
void setNoData(unsigned char* data, int length, unsigned char noDataValue);
void setNoData(float* data, int length, float noDataValue);
void setFlag(int index, unsigned char* flagArray);
bool isProcessed(int index, const unsigned char* flagArray);

bool  CreateGeoTIFF(char* path,int height, int width,void* pData, GDALDataType type, double* geoTransformArray6Eles,
					double* min, double* max, double* mean, double* stdDev, double nodatavalue);

bool  WriteTIFFStripe_MPI(char* path,int height, int width,GDALDataType type, Stripe* pStripe,double nodatavalue);

bool readTIFF(const char* path, GDALDataType type, CDEM& dem, double* geoTransformArray6Eles);
CDEM* diff(CDEM& demA, CDEM& demB);

bool readTIFFStripe(const char* path, GDALDataType type, CDEM& dem, int rowStart, int rowSize, int totalHeight);
bool readTIFFSize(const char* path,int& width, int& height);
int CreateDiffImage(const char *demA, const char *demB, char *resultPath, GDALDataType type, double nodatavalue);
extern const unsigned char value[8];
class Flag
{
public:
	unsigned int width, height;
	unsigned char* flagArray;
public:
	~Flag()
	{
		Free();
	}
	bool Init(int width,int height)
	{
		flagArray=NULL;
		this->width=width;
		this->height=height;
		flagArray = new unsigned char[(this->width * this->height + 7) / 8]();
		return flagArray!=NULL;
	}
	void Free()
	{
		delete[] flagArray;
		flagArray=NULL;
	}
	void ClearAllFlags()
	{
		unsigned int i=0; 
		unsigned int total=(width * height + 7) / 8;
		for (i=0; i<total; i++) flagArray[i]=0;
	}
	void RemoveFlag(int row, int col)
	{
		unsigned int index=(unsigned int)row*width+(unsigned int)col;
		flagArray[index / 8] &= ~value[index % 8];
	}

	void SetFlag(int row, int col)
	{
		unsigned int index=(unsigned int)row*width+(unsigned int)col;
		flagArray[index / 8] |= value[index % 8];
	}
	void SetFlags(int row, int col, Flag& flag)
	{
		unsigned int index=(unsigned int)row*width+(unsigned int)col;
		int bIndex=index / 8;
		int bShift=index % 8;
		flagArray[bIndex] |= value[bShift];
		flag.flagArray[bIndex] |= value[bShift];
	}
	int IsProcessed(int row, int col)
	{
		//if the cell is outside the DEM, is is regared as processed
		if (row<0 ||row>=height || col<0|| col>=width) return true;
		unsigned int index=(unsigned int)row*width+(unsigned int)col;
		return flagArray[index / 8]&value[index % 8];
	}
	int IsProcessedDirect(int row, int col)
	{
		unsigned int index=(unsigned int)row*width+(unsigned int)col;
		return flagArray[index / 8]&value[index % 8];
	}
};


#endif
