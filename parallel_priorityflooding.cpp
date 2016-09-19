#include "stripe.h"
#include <omp.h>
#include <mpi.h>
#include <unistd.h>

Stripe* CreateStrips(int height, int width, int stripNumber,CDEM* pDEM)
{
	int i;
	int isRowMajor;
	int linePerStrip;
	Stripe* pStrips=new Stripe[stripNumber];
	int total;
	if (height>=width) {
		isRowMajor=true;
		total=height;
	}
	else { 
		isRowMajor=false;
		total=width;
	}

	linePerStrip=total/stripNumber;	
	for (i=0; i<stripNumber; i++){
		pStrips[i].pDemStripe=pDEM;
		pStrips[i].stripeIndex=i;
		pStrips[i].stripeNumber=stripNumber;
		pStrips[i].isRowMajor=isRowMajor;
		pStrips[i].stripeSize=linePerStrip;
		if (isRowMajor)
		{
			pStrips[i].width=width;
			pStrips[i].height=linePerStrip;
			if (i==stripNumber-1)
				pStrips[i].height=total-i*linePerStrip;
		}
		else
		{
			pStrips[i].height=height;
			pStrips[i].width=linePerStrip;
			if (i==stripNumber-1)
				pStrips[i].width=total-i*linePerStrip;
		}

	}

	return pStrips;
}

Stripe* CreateStripe(int totalWidth, int totalHeight, int rank, int size)
{
	int i;
	int isRowMajor=true;
	int total=totalHeight;
	int linePerStripe=total/size;	

	Stripe* pStripe=new Stripe();
	pStripe->totalHeight=totalHeight;
	pStripe->totalWidth=totalWidth;
	pStripe->stripeIndex=rank;
	pStripe->stripeNumber=size;
	pStripe->isRowMajor=isRowMajor;
	pStripe->stripeSize=linePerStripe;
	pStripe->width=totalWidth;
	pStripe->height=linePerStripe;
	if (rank==size-1) pStripe->height=total-pStripe->stripeIndex*linePerStripe;
	return pStripe;
}

//parallel filling with MPI
void FillDEM_Parallel_MPI(char* inputFile, char* outputFilledPath)
{
	int rank,size;
	MPI_Init(NULL,NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);


	int totalWidth, totalHeight;

	double timeBeforeRead=MPI_Wtime();
	readTIFFSize(inputFile,totalWidth,totalHeight);
	
	if (rank==0)
		printf("Total size of DEM, width:%d, height:%d", totalWidth,totalHeight);
	
	double timeAfterRead = MPI_Wtime();
	//Create Stripe on this processor
	Stripe* pStripe=CreateStripe(totalWidth,totalHeight,rank,size);

	//read stripe into DEM (should be called DEMStripe)
	CDEM dem;
	readTIFFStripe(inputFile,GDALDataType::GDT_Float32, dem,pStripe->stripeSize*pStripe->stripeIndex,pStripe->height,totalHeight);
	if (rank==0) printf("\nstripe read, rank %d\n",rank);
	pStripe->pDemStripe=&dem;

	int iterationNumber=1;
	////fill each stripe for the first time
	pStripe->Initialize();
	pStripe->PushBorderCellsIntoPQ();
	pStripe->PriorityFlood(true);

	int total=0;
	int totalGrand;

	double communicationTime=0;
	double time1, time2;

	time1=MPI_Wtime();
	//update border and reduce by summation
	total=pStripe->UpdatePotentialSpillandPushToPQ_FirstTime(NULL);
	MPI_Allreduce (&total, &totalGrand, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	
	time2=MPI_Wtime();
	communicationTime+=time2-time1;
	
	//process stripes at other times
	while (true){
		if (totalGrand==0) break;
		pStripe->PriorityFlood(false);	
		
		time1=MPI_Wtime();
		total=pStripe->UpdatePotentialSpillandPushToPQ_OtherTimes(pStripe);
		MPI_Allreduce (&total, &totalGrand, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		time2=MPI_Wtime();
		communicationTime+=time2-time1;
		iterationNumber++;
	}
	pStripe->FillDepressionFromStripBorder();

	double timeBeforeWrite = MPI_Wtime();	
	//write filled stripe to output DEM in turn
	WriteTIFFStripe_MPI(outputFilledPath,totalHeight,totalWidth,GDALDataType::GDT_Float32, pStripe,-9999);	
	double timeAfterWrite = MPI_Wtime();	

	double computeTime=timeBeforeWrite-timeAfterRead;
	double readTime=timeAfterRead-timeBeforeRead;
	double writeTime=timeAfterWrite-timeBeforeWrite;
	double totalTime=timeAfterWrite-timeBeforeRead;

	double temp;
	MPI_Allreduce (&computeTime, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	computeTime = temp/size;

	MPI_Allreduce (&readTime, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	readTime = temp/size;

	MPI_Allreduce (&writeTime, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	writeTime = temp/size;

	MPI_Allreduce (&totalTime, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	totalTime = temp/size;

	MPI_Allreduce (&communicationTime, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	communicationTime = temp/size;

	delete pStripe;
	MPI_Finalize();
	if (rank==0)
	{
		printf("\nTime used: \nCompute Time: %f\nRead Time: %f\nWrite Time: %f\nCommunication Time: %f\nTotal Time: %f\n", computeTime,readTime, writeTime, communicationTime,totalTime);
		printf("\nNumber of iterations: %d\n",iterationNumber);
	}
}