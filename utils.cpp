#include "utils.h"
#include "gdal_priv.h"
#include "dem.h"
#include <string>
#include "stripe.h"
#include "mpi.h"
//create a new GeoTIFF file
bool  CreateGeoTIFF(char* path,int height, int width,void* pData, GDALDataType type, double* geoTransformArray6Eles,
					double* min, double* max, double* mean, double* stdDev, double nodatavalue)
{
    GDALDataset *poDataset;   
    GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
 
	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    char **papszOptions = NULL;
    poDataset = poDriver->Create(path,width, height, 1, type, 
                                papszOptions );
	

	if (geoTransformArray6Eles != NULL)
		poDataset->SetGeoTransform(geoTransformArray6Eles);


	GDALRasterBand* poBand;
	poBand= poDataset->GetRasterBand(1);
	
	poBand->SetNoDataValue(nodatavalue);

	if (min != NULL && max != NULL && mean != NULL && stdDev != NULL)
	{
		poBand->SetStatistics(*min, *max, *mean, *stdDev);
	}
	poBand->RasterIO( GF_Write, 0, 0, width, height, 
                      pData, width, height, type, 0, 0 );    

	GDALClose( (GDALDatasetH) poDataset );

	return true;
}

//create a new GeoTIFF file
bool  WriteTIFFStripe_MPI(char* path,int height, int width,GDALDataType type, Stripe* pStripe,double nodatavalue)
{
    GDALDataset *poDataset;   
    GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
 
	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
    char **papszOptions = NULL;
	if (pStripe->stripeIndex==0) //rank 0
	{
		poDataset = poDriver->Create(path,width, height, 1, type, 
									papszOptions );
		if (poDataset == NULL)
		{
			printf("Failed to create the GeoTIFF file\n");
			return false;
		}

		GDALRasterBand* poBand;
		poBand= poDataset->GetRasterBand(1);	
		poBand->SetNoDataValue(nodatavalue);
	
		poBand->RasterIO(GF_Write, 0, pStripe->stripeIndex*pStripe->stripeSize, pStripe->width, pStripe->height, 
						(void*)pStripe->pDemStripe->getDEMdata(), pStripe->width, pStripe->height, type, 0, 0 );    
		GDALClose( (GDALDatasetH) poDataset );
		if (pStripe->stripeIndex+1<pStripe->stripeNumber)
		{
			//send message to rank 1
			int d=0;
			MPI_Send(&d, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
		}
	}
	else
	{
		MPI_Status status;
		int d;
		MPI_Recv(&d, 1, MPI_INT, pStripe->stripeIndex - 1, 1, MPI_COMM_WORLD, &status);  //DGT check status to see that receive was correct.  Print status and rank
		
		//update this stripe 
		poDataset = (GDALDataset* )GDALOpen(path, GA_Update);
		if (poDataset == NULL)
		{
			printf("Failed to write the GeoTIFF file\n");
			return false;
		}
		GDALRasterBand* poBand;
		poBand = poDataset->GetRasterBand(1);
		poBand->RasterIO(GF_Write, 0, pStripe->stripeIndex*pStripe->stripeSize, pStripe->width, pStripe->height, 
						(void*)pStripe->pDemStripe->getDEMdata(), pStripe->width, pStripe->height, type, 0, 0 ); 

		GDALClose((GDALDatasetH)poDataset);
		if (pStripe->stripeIndex+1<pStripe->stripeNumber)
		{
			//send message to next rank 
			int d=0;
			MPI_Send(&d, 1, MPI_INT, pStripe->stripeIndex+1, 1, MPI_COMM_WORLD);
		}
	}
	return true;
}
bool readTIFFSize(const char* path,int& width, int& height)
{
	GDALDataset *poDataset;   
    GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	poDataset = (GDALDataset* )GDALOpen(path, GA_ReadOnly);
	if (poDataset == NULL)
	{
		printf("Failed to read the GeoTIFF file\n");
		return false;
	}

	GDALRasterBand* poBand;
	poBand = poDataset->GetRasterBand(1);

	width=poBand->GetXSize();
	height=poBand->GetYSize();

	GDALClose((GDALDatasetH)poDataset);

	return true;
}
//read a DEM stripe in GeoTIFF file 
//also read neighboring rows in adjacent rows
bool readTIFFStripe(const char* path, GDALDataType type, CDEM& dem, int rowStart, int rowSize, int totalHeight)
{
	GDALDataset *poDataset;   
    GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	poDataset = (GDALDataset* )GDALOpen(path, GA_ReadOnly);
	if (poDataset == NULL)
	{
		printf("Failed to read the GeoTIFF file\n");
		return false;
	}

	GDALRasterBand* poBand;
	poBand = poDataset->GetRasterBand(1);
	GDALDataType dataType = poBand->GetRasterDataType();
	if (dataType != type)
	{
		return false;
	}

	dem.SetWidth(poBand->GetXSize());
	dem.SetHeight(rowSize);
	
	if (!dem.Allocate()) return false;
	if (rowStart>0)
	{
		poBand->RasterIO(GF_Read, 0, rowStart-1, dem.Get_NX(), 1, 
			(void *)dem.getBorderInNeighbor1(), dem.Get_NX(), 1, dataType, 0, 0);
	}
	poBand->RasterIO(GF_Read, 0, rowStart, dem.Get_NX(), dem.Get_NY(), 
		(void *)dem.getDEMdata(), dem.Get_NX(), dem.Get_NY(), dataType, 0, 0);

	if (rowStart+rowSize<totalHeight)
	{
		poBand->RasterIO(GF_Read, 0, rowStart+rowSize, dem.Get_NX(), 1, 
			(void *)dem.getBorderInNeighbor2(), dem.Get_NX(), 1, dataType, 0, 0);
	}
	GDALClose((GDALDatasetH)poDataset);
	return true;
}
//read a DEM GeoTIFF file 
bool readTIFF(const char* path, GDALDataType type, CDEM& dem, double* geoTransformArray6Eles)
{
	GDALDataset *poDataset;   
    GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	poDataset = (GDALDataset* )GDALOpen(path, GA_ReadOnly);
	if (poDataset == NULL)
	{
		printf("Failed to read the GeoTIFF file\n");
		return false;
	}

	GDALRasterBand* poBand;
	poBand = poDataset->GetRasterBand(1);
	GDALDataType dataType = poBand->GetRasterDataType();
	if (dataType != type)
	{
		return false;
	}
	if (geoTransformArray6Eles == NULL)
	{
		printf("Transformation parameters can not be NULL\n");
		return false;
	}

	memset(geoTransformArray6Eles, 0, 6);
	poDataset->GetGeoTransform(geoTransformArray6Eles);

	dem.SetWidth(poBand->GetXSize());
	dem.SetHeight(poBand->GetYSize());
	
	if (!dem.Allocate()) return false;

	poBand->RasterIO(GF_Read, 0, 0, dem.Get_NX(), dem.Get_NY(), 
		(void *)dem.getDEMdata(), dem.Get_NX(), dem.Get_NY(), dataType, 0, 0);

	GDALClose((GDALDatasetH)poDataset);
	return true;
}
/*
*	neighbor index
*	5  6  7
*	4     0
*	3  2  1
*/
int	ix[8]	= { 0, 1, 1, 1, 0,-1,-1,-1 };
int	iy[8]	= { 1, 1, 0,-1,-1,-1, 0, 1 };

void setNoData(unsigned char* data, int length, unsigned char noDataValue)
{
	if (data == NULL || length == 0)
	{
		return;
	}

	for (int i = 0; i < length; i++)
	{
		data[i] = noDataValue;
	}
}
void setNoData(float* data, int length, float noDataValue)
{
	for (int i = 0; i < length; i++)
	{
		data[i] = noDataValue;
	}
}


const unsigned char value[8] = {128, 64, 32, 16, 8, 4, 2, 1};
