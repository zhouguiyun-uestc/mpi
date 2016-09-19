#include "dem.h"
#include "utils.h"
bool CDEM::Allocate()
{
	delete[] borderInNeighbor1;
	delete[] borderInNeighbor2;
	borderInNeighbor1=new float[width];
	borderInNeighbor2=new float[width];

	delete[] pDem;
	pDem=new float[width*height];
	if (pDem == NULL)
	{
		return false;
	}
	else
	{
		setNoData(pDem, width * height, NO_DATA_VALUE);
		return true;
	}

}
void CDEM::freeMem()
{
	delete[] pDem;
	pDem = NULL;
}
void CDEM::initialElementsNodata()
{
	setNoData(pDem, width * height, NO_DATA_VALUE);
}
float CDEM::asFloat(unsigned int row,unsigned int col) const
{
	return pDem[row*width+col];
}
void CDEM::Set_Value(unsigned int row,unsigned int col, float z)
{
	pDem[row*width+col]=z;
}

bool CDEM::is_NoData(unsigned int row, unsigned int col) const
{
	//for row-major only, to be updated
	if (row==-1)
	{
		if (fabs(borderInNeighbor1[col]-NO_DATA_VALUE)<0.00001) return true;
		return false;
	}
	else if (row==height){
		if (fabs(borderInNeighbor2[col]-NO_DATA_VALUE)<0.00001) return true;
		return false;
	}
	else {
		if (fabs(pDem[row*width+col]-NO_DATA_VALUE)<0.00001) return true;
		return false;
	}
}
void CDEM::Assign_NoData()
{
	for (int i=0; i<width*height; i++)
		pDem[i]=NO_DATA_VALUE;
}
int CDEM::Get_NY() const
{
	return height;
}
int CDEM::Get_NX() const
{
	return width;
}
float* CDEM::getDEMdata() const
{
	return pDem;
}

float* CDEM::getBorderInNeighbor1()
{
	return this->borderInNeighbor1;
}

float* CDEM::getBorderInNeighbor2()
{
	return this->borderInNeighbor2;
}

void CDEM::SetHeight(unsigned int height) 
{
	this->height = height;
}
void CDEM::SetWidth(unsigned int width)
{
	this->width = width;
}

void CDEM::readDEM(const std::string& filePath)
{
	//std::ifstream is;
	//is.open (filePath, std::ios::binary);
	//is.read((char*)pDem,sizeof(float)*width*height);
	//is.close();
}
bool CDEM::is_InGrid(unsigned int row, unsigned int col) const
{
	if ((row >=0 && row < height) && (col >= 0 && col < width))
		return true;

	return false;
}

float CDEM::getLength(unsigned int dir)
{
	if ((dir & 0x1) == 1)
	{
		return 1.41421f;
	}
	else return 1.0f;
}

//compute flow directoin 
unsigned char CDEM::getDirction(unsigned int row, unsigned int col, float spill)
{

	int iRow, iCol;
	float iSpill, max, gradient;
	unsigned char steepestSpill;
	max = 0.0f;
	steepestSpill = 255;
	unsigned char lastIndexINGridNoData = 0;
	for (int i = 0; i < 8; i++)
	{
		iRow = Get_rowTo(i, row);
		iCol = Get_colTo(i, col);

		if (is_InGrid(iRow, iCol) && !is_NoData(iRow, iCol) && (iSpill = asFloat(iRow, iCol)) < spill)
		{
			gradient = (spill - iSpill) / getLength(i);
			if (max < gradient)
			{
				max = gradient;
				steepestSpill = i;
			}
		}
		if (!is_InGrid(iRow, iCol) || is_NoData(iRow, iCol))
		{
			lastIndexINGridNoData = i;
		}

	}
	return steepestSpill != 255 ? dir[steepestSpill] : dir[lastIndexINGridNoData];
}
