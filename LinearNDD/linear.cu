#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>

#define Mu0 1.256637e-6
#define pi 3.1415927
#define MuFeCore 0.0020 //linear region of BH curve
#define MaxNode 8000
#define MaxElem 16000
#define Init0 0.0
#define NRCount 1

#define CudaThrdNum 128
#define CudaBlckNum 128

#ifndef _TIMER_H_
#define _TIMER_H_
typedef struct
{
	double Y;
	double Y1;
	double Y2;
}interp_t;
#define GET_TIME(now){ \
	double t; \
	t=clock(); \
	now =t/CLOCKS_PER_SEC ;\
}
#endif

typedef struct
{
	int Id;
	double x;
	double y;
	int NumEle;
	int EleID[10];
	int EleOrd[10];
	int NeiborNode[20];
	int NumNeiborNodes;
	double TotalArea;
	int Type;
	double A0;
	double A1;
	double K;
	double SumRHSContri; // Sum from all elems involved with this node
	double JsSum;
	double SumNeiborJsSum;//from Js
}FEMNode;
typedef struct
{
	int Id;
	int Nodes[3];
	double Area;
	int Type;
	double Ve;
	double Js;
	double sigma;
	double Me[3][3];
	double ElmRowSum[3][3];// weighted row sum by gamma1 for 3 nodes
	double RHSContri[3]; // RHS or b calculated by elmK*A on 3 nodes
}FEMElem;


//Host parameters for FEM
int NumNodes, NumElem;
FEMNode MyNode[MaxNode];
FEMElem MyElem[MaxElem];
double gamma1 = 100.0;
double CurrentDensity = 1e6;
//Device parameters for FEM
int *d_NumNodes, *d_NumElem;
double *d_gamma1;
double *d_CurrentDensity;
FEMNode *d_MyNode;
FEMElem *d_MyElem;

void LoadMeshInfo()
{
	FILE* ip;
	int i, flag = 0;
	char filename[50];

	char line[50];


	sprintf(filename, "3733.mphtxt");

	if ((ip = fopen(filename, "r")) == NULL) {
		printf("error opening the input data.\n");
	}
	while (fgets(line, sizeof(line), ip) != NULL) {

		if (strstr(line, "# number of mesh points") != NULL)
		{
			sscanf(line, "%d", &(NumNodes));
		}
	}
	fclose(ip);


	if ((ip = fopen(filename, "r")) == NULL) {
		printf("error opening the input data.\n");
	}

	while (fgets(line, sizeof(line), ip) != NULL) {

		if (strstr(line, "# Mesh point coordinates") != NULL)
		{
			for (i = 0; i < NumNodes; i++)
			{
				fgets(line, sizeof(line), ip);
				sscanf(line, "%lf %lf\n", &(MyNode[i].x), &(MyNode[i].y));
			}
		}
	}
	fclose(ip);



	if ((ip = fopen(filename, "r")) == NULL) {
		printf("error opening the input data.\n");
	}

	while (fgets(line, sizeof(line), ip) != NULL) {

		if (strstr(line, "# number of elements") != NULL)
		{
			flag = flag + 1;
			if (flag == 3)
			{
				sscanf(line, "%d", &(NumElem));
				fgets(line, sizeof(line), ip);

				for (i = 0; i < NumElem; i++)
				{
					fgets(line, sizeof(line), ip);
					sscanf(line, "%d %d %d\n", &(MyElem[i].Nodes[0]), &(MyElem[i].Nodes[1]), &(MyElem[i].Nodes[2]));

				}
			}
		}
	}
	fclose(ip);

	flag = 0;
	if ((ip = fopen(filename, "r")) == NULL) {
		printf("error opening the input data.\n");
	}
	while (fgets(line, sizeof(line), ip) != NULL) {

		if (strstr(line, "# Geometric entity indices") != NULL)
		{
			flag = flag + 1;
			if (flag == 3)
			{
				for (i = 0; i < NumElem; i++)
				{
					fgets(line, sizeof(line), ip);
					sscanf(line, "%d\n", &(MyElem[i].Type));
					//	printf("\nType %d", MyElem[i].Type);
				}
			}
		}
	}
	fclose(ip);



}

void FEM_Host_Data_Prepare()
{
	int i, j, k = 0;
	double x1, x2, x3, y1, y2, y3;
	double b1, b2, b3, c1, c2, c3;
	static int dummy[MaxNode];
	LoadMeshInfo();

	//Boundary node ID and values
	for (i = 0; i < NumNodes; i++)
	{
		MyNode[i].A0 = Init0;
		MyNode[i].Type = 0;
		//MyNode[i].A0 = (double)i + 1.0;
		if (fabs(MyNode[i].x - 3.5) < 1e-5 || fabs(MyNode[i].x + 3.5) < 1e-5 || fabs(MyNode[i].y - 2.5) < 1e-5 || fabs(MyNode[i].y + 2.5) < 1e-5)

		{
			MyNode[i].Type = 1;


		}

	}

	//set initial Ve config
	for (i = 0; i < NumElem; i++)

	{
		MyElem[i].Js = 0;
		if (1)//MyElem[i].Type != 2)
			MyElem[i].Ve = 1.0 / Mu0;
		else
			MyElem[i].Ve = 1.0 / MuFeCore;

	}


	//Element K matrix 
	for (i = 0; i < NumElem; i++)
	{
		x1 = MyNode[MyElem[i].Nodes[0]].x; x2 = MyNode[MyElem[i].Nodes[1]].x; x3 = MyNode[MyElem[i].Nodes[2]].x;
		y1 = MyNode[MyElem[i].Nodes[0]].y; y2 = MyNode[MyElem[i].Nodes[1]].y; y3 = MyNode[MyElem[i].Nodes[2]].y;
		MyElem[i].Area = 0.5*(x1*(y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
		b1 = y2 - y3; c1 = x3 - x2;
		b2 = y3 - y1; c2 = x1 - x3;
		b3 = y1 - y2; c3 = x2 - x1;

		MyElem[i].Me[0][0] = 1.0 / 4 / MyElem[i].Area*(b1*b1 + c1 * c1);
		MyElem[i].Me[0][1] = 1.0 / 4 / MyElem[i].Area*(b1*b2 + c1 * c2);
		MyElem[i].Me[0][2] = 1.0 / 4 / MyElem[i].Area*(b1*b3 + c1 * c3);

		MyElem[i].Me[1][0] = 1.0 / 4 / MyElem[i].Area*(b1*b2 + c1 * c2);
		MyElem[i].Me[1][1] = 1.0 / 4 / MyElem[i].Area*(b2*b2 + c2 * c2);
		MyElem[i].Me[1][2] = 1.0 / 4 / MyElem[i].Area*(b2*b3 + c2 * c3);

		MyElem[i].Me[2][0] = 1.0 / 4 / MyElem[i].Area*(b1*b3 + c1 * c3);
		MyElem[i].Me[2][1] = 1.0 / 4 / MyElem[i].Area*(b3*b2 + c3 * c2);
		MyElem[i].Me[2][2] = 1.0 / 4 / MyElem[i].Area*(b3*b3 + c3 * c3);

		//get ElmRowSum[i][...] is gamma1 weighted row sum for ith node. 
				//set 0
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				MyElem[i].ElmRowSum[j][k] = 0;

		double temp;
		for (j = 0; j < 3; j++)
			for (k = 0; k < 3; k++)
				if (MyNode[MyElem[i].Nodes[k]].Type != 1) //first type bdry node row is set to 0
				{
					if (k == j)
						temp = 1.0;
					else
						temp = 1.0 / gamma1;
					MyElem[i].ElmRowSum[j][0] += temp * MyElem[i].Me[k][0];
					MyElem[i].ElmRowSum[j][1] += temp * MyElem[i].Me[k][1];
					MyElem[i].ElmRowSum[j][2] += temp * MyElem[i].Me[k][2];
				}

	}


	// topology link info of Nodes
	for (i = 0; i < NumNodes; i++)
	{
		MyNode[i].NumEle = 0;
		MyNode[i].TotalArea = 0;
		MyNode[i].NumNeiborNodes = 0;
	}
	for (i = 0; i < NumElem; i++)
		for (j = 0; j < 3; j++)
		{
			MyNode[MyElem[i].Nodes[j]].EleID[MyNode[MyElem[i].Nodes[j]].NumEle] = i;
			MyNode[MyElem[i].Nodes[j]].EleOrd[MyNode[MyElem[i].Nodes[j]].NumEle] = j;
			MyNode[MyElem[i].Nodes[j]].NumEle++;
			MyNode[MyElem[i].Nodes[j]].TotalArea = MyNode[MyElem[i].Nodes[j]].TotalArea + MyElem[i].Area;
		}
	for (i = 0; i < NumNodes; i++)
	{
		for (j = 0; j < NumNodes; j++)
			dummy[j] = 0;
		for (j = 0; j < MyNode[i].NumEle; j++)
			for (k = 0; k < 3; k++)
				dummy[MyElem[MyNode[i].EleID[j]].Nodes[k]] = 1;
		for (j = 0; j < NumNodes; j++)
			if (dummy[j] == 1 && j != i)
			{
				MyNode[i].NeiborNode[MyNode[i].NumNeiborNodes] = j;
				(MyNode[i].NumNeiborNodes)++;
			}

	}

	for (i = 0; i < NumNodes; i++)
	{
		if (MyNode[i].Type == 0)
		{
			MyNode[i].K = 0;
			for (j = 0; j < MyNode[i].NumEle; j++)
				MyNode[i].K = MyNode[i].K + MyElem[MyNode[i].EleID[j]].Me[MyNode[i].EleOrd[j]][MyNode[i].EleOrd[j]];
		}
	}
	printf("Number of Nodes: %d   Number of Elements: %d\n", NumNodes, NumElem);

}
__global__ void ApplyCurrent(double *d_CurrentDensity, int *d_NumElem, FEMElem* d_MyElem)
{
	int e;
	e = threadIdx.x + blockIdx.x*blockDim.x;
	if (e >= *d_NumElem)
		return;
	if (d_MyElem[e].Type == 5)  d_MyElem[e].Js = (*d_CurrentDensity);
	if (d_MyElem[e].Type == 7)  d_MyElem[e].Js = -(*d_CurrentDensity);
}
__global__ void JsSumCalculate(int *d_NumNodes, FEMElem* d_MyElem, FEMNode* d_MyNode)
{
	int node = threadIdx.x + blockIdx.x*blockDim.x;
	if (node >= *d_NumNodes)
		return;
	int e;
	d_MyNode[node].JsSum = 0;
	for (e = 0; e < d_MyNode[node].NumEle; e++)
		d_MyNode[node].JsSum += (d_MyElem[d_MyNode[node].EleID[e]].Js)*(d_MyElem[d_MyNode[node].EleID[e]].Area) / 3;
}
__global__ void SumNeiborJsSumCalculate(int *d_NumNodes, FEMElem* d_MyElem, FEMNode* d_MyNode)
{
	int node = threadIdx.x + blockIdx.x*blockDim.x;
	if (node >= *d_NumNodes)
		return;
	int n;
	d_MyNode[node].SumNeiborJsSum = 0;
	for (n = 0; n < d_MyNode[node].NumNeiborNodes; n++)
		d_MyNode[node].SumNeiborJsSum += (d_MyNode[d_MyNode[node].NeiborNode[n]].JsSum);
}
__global__ void SumNodeRHSContriCalculate(int *d_NumNodes, FEMElem* d_MyElem, FEMNode* d_MyNode)
{
	int node = threadIdx.x + blockIdx.x*blockDim.x;
	if (node >= *d_NumNodes)
		return;
	int e;
	d_MyNode[node].SumRHSContri = 0;
	for (e = 0; e < d_MyNode[node].NumEle; e++)
		d_MyNode[node].SumRHSContri += d_MyElem[d_MyNode[node].EleID[e]].RHSContri[d_MyNode[node].EleOrd[e]];

}
__global__ void ElmRHSContriCalculate(int *d_NumElem, FEMElem* d_MyElem, FEMNode* d_MyNode)
{

	int e;
	e = threadIdx.x + blockIdx.x*blockDim.x;
	if (e >= *d_NumElem)
		return;

	int row, col;
	for (row = 0; row < 3; row++)
		d_MyElem[e].RHSContri[row] = 0;
	for (row = 0; row < 3; row++)
	{
		for (col = 0; col < 3; col++)
			d_MyElem[e].RHSContri[row] += d_MyElem[e].Me[row][col] * d_MyNode[d_MyElem[e].Nodes[col]].A0;
		d_MyElem[e].RHSContri[row] = d_MyElem[e].RHSContri[row] * d_MyElem[e].Ve;
	}

}
__global__ void UpdateSolutiontoA1(int *d_NumNodes, FEMElem* d_MyElem, FEMNode* d_MyNode, double *d_gamma1)
{
	int node = threadIdx.x + blockIdx.x*blockDim.x;
	if (node >= *d_NumNodes)
		return;

	double RHS, LHS, dF_dA, NRsolition, temp;
	int i, j, NeiborID, e, LocalPos;
	double ALocal[3], B2, B, V, VB2, B2A;
	int NRct;
	if (d_MyNode[node].Type != 1)
	{
		//----------Get RHS
		RHS = 0;
		// subdomain RHS sum  = SumofAllneibornodes-Sumof subdomain elems's contribution on neibornodes 
		for (i = 0; i < d_MyNode[node].NumNeiborNodes; i++)
		{
			NeiborID = d_MyNode[node].NeiborNode[i];
			if (d_MyNode[NeiborID].Type != 1)
				RHS -= d_MyNode[NeiborID].SumRHSContri;
		}

		for (i = 0; i < d_MyNode[node].NumEle; i++)
		{
			e = d_MyNode[node].EleID[i];
			for (j = 0; j < 3; j++)
				if (node != d_MyElem[e].Nodes[j] && d_MyNode[d_MyElem[e].Nodes[j]].Type != 1)
					RHS += d_MyElem[e].RHSContri[j];
		}

		RHS = (RHS + d_MyNode[node].SumNeiborJsSum) / (*d_gamma1) + d_MyNode[node].JsSum;
		//----------Get LHS and dF_dA

		NRsolition = d_MyNode[node].A0;
		//----------NRiteration
		for (NRct = 0; NRct < NRCount; NRct++)
		{
			LHS = 0;
			dF_dA = 0;
			for (i = 0; i < d_MyNode[node].NumEle; i++)
			{
				e = d_MyNode[node].EleID[i];
				LocalPos = d_MyNode[node].EleOrd[i];
				ALocal[0] = d_MyNode[d_MyElem[e].Nodes[0]].A0;
				ALocal[1] = d_MyNode[d_MyElem[e].Nodes[1]].A0;
				ALocal[2] = d_MyNode[d_MyElem[e].Nodes[2]].A0;
				ALocal[LocalPos] = NRsolition;
				if (0)
				{
					B2 = -1 / d_MyElem[e].Area*(d_MyElem[e].Me[0][1] * pow(ALocal[0] - ALocal[1], 2)
						+ d_MyElem[e].Me[1][2] * pow(ALocal[1] - ALocal[2], 2)
						+ d_MyElem[e].Me[2][0] * pow(ALocal[2] - ALocal[0], 2));
					B = sqrt(B2);

					if (B <= 0.6)
					{
						for (j = 0; j < 3; j++)
							LHS += ALocal[j] * d_MyElem[e].ElmRowSum[LocalPos][j] / MuFeCore;
						dF_dA += d_MyElem[e].ElmRowSum[LocalPos][LocalPos] / MuFeCore;
						//		if (NRct == (NRCount - 1))
								//	d_MyElem[e].Ve = 1 / MuFeCore;
					}
					else
					{
						V = 1 / MuFeCore + 3000.0*pow(B - 0.6, 3) / B;
						VB2 = (B*9000.0*pow(B - 0.6, 2) - 3000.0*pow(B - 0.6, 3)) / B / B / 2 / B;
						B2A = 0;

						for (j = 0; j < 3; j++)
							if (j != LocalPos)
								B2A = B2A + d_MyElem[e].Me[LocalPos][j] * (ALocal[LocalPos] - ALocal[j]);
						B2A = -B2A * 2 / d_MyElem[e].Area;
						// LHS+=V* Alocal dot ElmRowSum(Localpos,:) 
						temp = 0;
						for (j = 0; j < 3; j++)
						{
							dF_dA += VB2 * B2A * ALocal[j] * d_MyElem[e].ElmRowSum[LocalPos][j];
							temp += ALocal[j] * d_MyElem[e].ElmRowSum[LocalPos][j];
						}
						dF_dA += V * d_MyElem[e].ElmRowSum[LocalPos][LocalPos];
						LHS += temp * V;
						//		if (NRct == (NRCount - 1))
							//		d_MyElem[e].Ve = V;
					}
				}
				else
				{
					for (j = 0; j < 3; j++)
						LHS += ALocal[j] * d_MyElem[e].ElmRowSum[LocalPos][j] / Mu0;
					dF_dA += d_MyElem[e].ElmRowSum[LocalPos][LocalPos] / Mu0;
				}
			}
			NRsolition += (RHS - LHS) / dF_dA;
			//	printf("\nB%e NR%e RHS%e LHS%e dFda %e\n",B, NRsolition, RHS,LHS, dF_dA);
		}
		d_MyNode[node].A1 = NRsolition;

	}
}
__global__ void CopyA1ToA0(int *d_NumNodes, FEMElem* d_MyElem, FEMNode* d_MyNode)
{
	int i = threadIdx.x + blockIdx.x*blockDim.x;
	if (i >= *d_NumNodes)
		return;
	d_MyNode[i].A0 = d_MyNode[i].A1;
}
__global__ void UpdateVe(int *d_NumElem, FEMElem* d_MyElem, FEMNode* d_MyNode)
{
	int e;
	e = threadIdx.x + blockIdx.x*blockDim.x;
	if (e >= *d_NumElem || 1)
		return;
	double ALocal[3], B2, B, V;
	ALocal[0] = d_MyNode[d_MyElem[e].Nodes[0]].A0;
	ALocal[1] = d_MyNode[d_MyElem[e].Nodes[1]].A0;
	ALocal[2] = d_MyNode[d_MyElem[e].Nodes[2]].A0;
	B2 = -1 / d_MyElem[e].Area*(d_MyElem[e].Me[0][1] * pow(ALocal[0] - ALocal[1], 2)
		+ d_MyElem[e].Me[1][2] * pow(ALocal[1] - ALocal[2], 2)
		+ d_MyElem[e].Me[2][0] * pow(ALocal[2] - ALocal[0], 2));
	B = sqrt(B2);
	if (B <= 0.6)
		V = 1 / MuFeCore;
	else
		V = 1 / MuFeCore + 3000.0*pow(B - 0.6, 3) / B;
	d_MyElem[e].Ve = V;
}
__global__ void Updategamma1OnDevice(int *d_NumElem, FEMElem* d_MyElem, FEMNode* d_MyNode, double *d_gamma1)
{
	int e;
	e = threadIdx.x + blockIdx.x*blockDim.x;
	if (e >= *d_NumElem)
		return;
	int j, k;
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
			d_MyElem[e].ElmRowSum[j][k] = 0;

	double temp;
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
			if (d_MyNode[d_MyElem[e].Nodes[k]].Type != 1) //first type bdry node row is set to 0
			{
				if (k == j)
					temp = 1.0;
				else
					temp = 1.0 / (*d_gamma1);
				d_MyElem[e].ElmRowSum[j][0] += temp * d_MyElem[e].Me[k][0];
				d_MyElem[e].ElmRowSum[j][1] += temp * d_MyElem[e].Me[k][1];
				d_MyElem[e].ElmRowSum[j][2] += temp * d_MyElem[e].Me[k][2];
			}
}

void GPUInitialMallocCopy()
{
	//cuda initiallze
	int num_devices, device;
	cudaGetDeviceCount(&num_devices);
	if (num_devices > 1) {
		int max_multiprocessors = 0;
		for (device = 0; device < num_devices; device++) {
			cudaDeviceProp properties;
			cudaGetDeviceProperties(&properties, device);
			if (max_multiprocessors < properties.multiProcessorCount) {
				max_multiprocessors = properties.multiProcessorCount;
				
			}
		}

	}
	cudaSetDevice(0);
	cudaMalloc((void**)&d_MyNode, NumNodes * sizeof(FEMNode));
	cudaMalloc((void**)&d_MyElem, NumElem * sizeof(FEMElem));
	cudaMalloc((void**)&d_NumNodes, sizeof(int));
	cudaMalloc((void**)&d_NumElem, sizeof(int));
	cudaMalloc((void**)&d_gamma1, sizeof(double));
	cudaMalloc((void**)&d_CurrentDensity, sizeof(double));
	cudaMemcpy(d_MyNode, MyNode, NumNodes * sizeof(FEMNode), cudaMemcpyHostToDevice);
	cudaMemcpy(d_MyElem, MyElem, NumElem * sizeof(FEMElem), cudaMemcpyHostToDevice);
	cudaMemcpy(d_NumNodes, &NumNodes, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_NumElem, &NumElem, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_gamma1, &gamma1, sizeof(double), cudaMemcpyHostToDevice);
}
void GPUFree()
{
	cudaFree(d_NumElem);
	cudaFree(d_NumNodes);
	cudaFree(d_MyElem);
	cudaFree(d_MyNode);
	cudaFree(d_gamma1);
	cudaFree(d_CurrentDensity);
}
//for the test
double My_abs(double x)
{
	if (x >= 0)
		return x;
	else
		return -x;
}
double PorbeValueVsItration[4000];
void MyUpdategamma1(double danteng)
{
	cudaMemcpy(d_gamma1, &danteng, sizeof(double), cudaMemcpyHostToDevice);
	Updategamma1OnDevice << <CudaBlckNum, CudaThrdNum >> > (d_NumElem, d_MyElem, d_MyNode, d_gamma1);
}
int main() {

	FEM_Host_Data_Prepare();
	GPUInitialMallocCopy();
	int probeID = -1;
	//find the beloved probe node
	for (int i = 1; i < NumNodes; i++)
		if (My_abs(0.25 - MyNode[i].x) < 1e-13)
			if (My_abs(1.30 - MyNode[i].y) < 1e-13)
				probeID = i;
	printf("\nprobeID=%d\n", probeID);


	//Set currentdensity on device
	CurrentDensity = 1e6;
	cudaMemcpy(d_CurrentDensity, &CurrentDensity, sizeof(double), cudaMemcpyHostToDevice);
	ApplyCurrent << <CudaBlckNum, CudaThrdNum >> > (d_CurrentDensity, d_NumElem, d_MyElem);

	//Js contribution prepare
	JsSumCalculate << <CudaBlckNum, CudaThrdNum >> > (d_NumNodes, d_MyElem, d_MyNode);
	SumNeiborJsSumCalculate << <CudaBlckNum, CudaThrdNum >> > (d_NumNodes, d_MyElem, d_MyNode);
	//gamma1 prepare
	gamma1 = 5.0;
	cudaMemcpy(d_gamma1, &gamma1, sizeof(double), cudaMemcpyHostToDevice);
	Updategamma1OnDevice << <CudaBlckNum, CudaThrdNum >> > (d_NumElem, d_MyElem, d_MyNode, d_gamma1);

	cudaEvent_t start, stop;//unit: ms

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	for (int RelaxCount = 0; RelaxCount < 1000; RelaxCount++)
	{
		//if (RelaxCount == 2000) MyUpdategamma1(4.9);
		//if (RelaxCount ==40) MyUpdategamma1(3.3);
		//if (RelaxCount == 30) MyUpdategamma1(5.0);
		//if (RelaxCount == 2000) MyUpdategamma1(6.2);
		ElmRHSContriCalculate << <CudaBlckNum, CudaThrdNum >> > (d_NumElem, d_MyElem, d_MyNode);

		SumNodeRHSContriCalculate << <CudaBlckNum, CudaThrdNum >> > (d_NumNodes, d_MyElem, d_MyNode);

		UpdateSolutiontoA1 << <CudaBlckNum, CudaThrdNum >> > (d_NumNodes, d_MyElem, d_MyNode, d_gamma1);

		CopyA1ToA0 << <CudaBlckNum, CudaThrdNum >> > (d_NumNodes, d_MyElem, d_MyNode);

		UpdateVe << <CudaBlckNum, CudaThrdNum >> > (d_NumElem, d_MyElem, d_MyNode);

		cudaMemcpy(&PorbeValueVsItration[RelaxCount], &(d_MyNode[probeID].A0), sizeof(double), cudaMemcpyDeviceToHost);
	}
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	printf("\nElapsed time for NDDR Iteration is %.4f ms\n", elapsedTime);
	cudaMemcpy(&MyElem, d_MyElem, NumElem * sizeof(FEMElem), cudaMemcpyDeviceToHost);
	cudaMemcpy(&MyNode, d_MyNode, NumNodes * sizeof(FEMNode), cudaMemcpyDeviceToHost);
	GPUFree();

	printf("\n%e\n", MyNode[probeID].A0);
	//for the test
	FILE *write_ptr;
	write_ptr = fopen("PorbeValueVsItration.bin", "wb");  // w for write, b for binary
	fwrite(PorbeValueVsItration, sizeof(PorbeValueVsItration), 1, write_ptr);
	//for the test

	return 0;
}


