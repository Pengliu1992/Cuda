#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//define GET_TIME
#ifndef _TIMER_H_
#define _TIMER_H_
typedef struct
{
	double Y;
	double Y1;
	double Y2;
}interp_t;
#include <time.h>

#define GET_TIME(now){ \
	double t; \
	t=clock(); \
	now =t/CLOCKS_PER_SEC ;\
}
#endif

#define CudaBlckNum 32
#define CudaThrdNum 32
#define NumOfNodes 320
#define NumOfElems 578
#define Delt (0.01)

typedef struct
{
	int NumEle;
	int EleID[10];
	int EleOrd[10];
	int Type;
	float Jext;
	float Jeddy = 0.0;
	float A;// not in use-> relpaced by float *d_x outside of struct
	float A1;// not in use-> relpaced by float *d_x outside of struct
	float TLM_Eq_Jsource = 0.0;

}FEMNode;
typedef struct
{
	int I, J, K;
	int Type;
	float Ve;
	float Kc[3][3];//Kc ��G01 G12����ϵġ�ע��ֻ��������һ��
	float Dp[3][3];
	float MatrixMultiVecBuff[3];
	float Dp1;
	float Dp2;
	//TLM vars
	float Area;
	float Vr[3] = { 0.0,0.0,0.0 };// Vr[0]=Vr01,[1]=12,[2]=20
	float Vi[3] = { 0.0,0.0,0.0 };
	// Vr , Vi is with respect to nonlinear element resistors
	float Ieq[3];// equvlt currents by Vr from nonlinear resistors
}FEMElem;

//FEM vars
//------------FEM vars host
FEMNode MyNode[NumOfNodes];
FEMElem MyElem[NumOfElems];
int NumNodes, NumElems;
//------------FEM vars device
FEMNode *d_MyNode;
FEMElem *d_MyElem;
int *d_NumNodes, *d_NumElems;
int *d_WannaUpdateVe;
float *d_yinuonuo, *d_anuonuo;//dummys watchers for debug


//PCG deviece vars solving Ax=b
typedef struct
{
	float r1, r0, z1, z0, b, Ap, PreCondi;
}PcgVec;
PcgVec *d_PcgVec;// PCG solving vector var package floating clouds 
float *d_x;// x= A,the mag potential on each nodes. cmmunication vects
float *d_p;// cmmunication vects
float *d_buffer1, *d_buffer2; //reduction buffers
float *d_alpha, *d_beta, *d_pap, *d_rz1, *d_rz0; //floating clouds PCG steping vars pap=p'*A*p; rz1=r1'*z1 ...

// prepare data and memory control
void LoadMeshInfoAndPrepareFEM()
{
	FILE* ip;
	int i, flag = 0;
	char filename[50];
	
	char line[50];

	sprintf(filename, "Untitled.mphtxt");


	//load mesh connections and coord
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
	printf("Num_nodes = %d \n", NumNodes);

	double *X, *Y;
	X = (double *)malloc(NumNodes * sizeof(double));
	Y = (double *)malloc(NumNodes * sizeof(double));
	while (fgets(line, sizeof(line), ip) != NULL) {

		if (strstr(line, "# Mesh point coordinates") != NULL)
		{
			for (i = 0; i < NumNodes; i++)
			{
				fgets(line, sizeof(line), ip);
				sscanf(line, "%lf %lf\n", &X[i], &Y[i]);
			}
		}
	}
	fclose(ip);

	if ((ip = fopen(filename, "r")) == NULL) {
		printf("error opening the input data.\n");
	}
	flag = 0;
	while (fgets(line, sizeof(line), ip) != NULL) {

		if (strstr(line, "# number of elements") != NULL)
		{
			flag = flag + 1;
			if (flag == 3)
			{
				sscanf(line, "%d", &(NumElems));
				fgets(line, sizeof(line), ip);
				printf("Num_elems = %d \n", NumElems);
				for (i = 0; i < NumElems; i++)
				{
					fgets(line, sizeof(line), ip);
					sscanf(line, "%d %d %d\n", &(MyElem[i].I), &(MyElem[i].J), &(MyElem[i].K));
				}
			}
		}
	}
	fclose(ip);

	if ((ip = fopen(filename, "r")) == NULL) {
		printf("error opening the input data.\n");
	}
	flag = 0;
	while (fgets(line, sizeof(line), ip) != NULL) {

		if (strstr(line, "# number of geometric entity indices") != NULL)
		{
			flag = flag + 1;
			if (flag == 3)
			{
				sscanf(line, "%d", &(NumElems));
				fgets(line, sizeof(line), ip);
				//printf("Num_elems = %d \n", NumElems);
				for (i = 0; i < NumElems; i++)
				{
					fgets(line, sizeof(line), ip);
					sscanf(line, "%d\n", &(MyElem[i].Type));
				}
			}
		}
	}
	fclose(ip);
	//load mesh connections and coord finished

	//link NDD topology according to mesh connecitons
	for (i = 0; i < NumNodes; i++)
		MyNode[i].NumEle = 0;
	int nownode;
	for (i = 0; i < NumElems; i++)
	{
		nownode = MyElem[i].I;
		MyNode[nownode].EleID[MyNode[nownode].NumEle] = i;
		MyNode[nownode].EleOrd[MyNode[nownode].NumEle] = 0;
		MyNode[nownode].NumEle += 1;
		nownode = MyElem[i].J;
		MyNode[nownode].EleID[MyNode[nownode].NumEle] = i;
		MyNode[nownode].EleOrd[MyNode[nownode].NumEle] = 1;
		MyNode[nownode].NumEle += 1;
		nownode = MyElem[i].K;
		MyNode[nownode].EleID[MyNode[nownode].NumEle] = i;
		MyNode[nownode].EleOrd[MyNode[nownode].NumEle] = 2;
		MyNode[nownode].NumEle += 1;
	}
	//link NDD topology according to mesh connecitons finshed

	//mark boundary nodes and set Jext to 0
	for (i = 0; i < NumNodes; i++)
	{
		if (abs((long)X[i]) >= 0.99 || abs((long)Y[i]) >= 0.99)
			MyNode[i].Type = 1;
		else
			MyNode[i].Type = 0;
		MyNode[i].Jext = 0.0;
		MyNode[i].A = 0.0;
		MyNode[i].A1 = 0.0;
	}
	//mark boundary nodes and set Jext to 0 finished

	//Get Stiff and Dampinmg matrix and add to Node.Jext
	double x1, y1, x2, y2, x3, y3, Area;
	double b1, c1, b2, c2, b3, c3;
	double sigma, Jext;
	int I, J, K;
	double G01, G20, G12;
	for (i = 0; i < NumElems; i++)
	{
		sigma = 1.0;
		MyElem[i].Ve = 1.0;
		Jext = 0.0;
		if (MyElem[i].Type != -1000)
		{
			sigma = 1.0;
			MyElem[i].Ve = 1.0;
			Jext = 3.21;
		}
		I = MyElem[i].I; J = MyElem[i].J; K = MyElem[i].K;
		x1 = X[I]; x2 = X[J]; x3 = X[K];
		y1 = Y[I]; y2 = Y[J]; y3 = Y[K];
		Area = 0.5*(x1*(y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
		MyElem[i].Area = Area;
		b1 = y2 - y3; c1 = x3 - x2;
		b2 = y3 - y1; c2 = x1 - x3;
		b3 = y1 - y2; c3 = x2 - x1;

		MyElem[i].Dp1 = (float)(sigma *Area / 12.0 / Delt);
		MyElem[i].Dp2 = (float)(sigma * Area / 6.0 / Delt);

		for (int ii = 0; ii < 3; ii++)
			for (int jj = 0; jj < 3; jj++)
				MyElem[i].Dp[ii][jj] = MyElem[i].Dp1;
		for (int ii = 0; ii < 3; ii++)
			MyElem[i].Dp[ii][ii] = MyElem[i].Dp2;

		G01 = 1.0 / 4  / Area * (b1 * b2 + c1 * c2);
		G20 = 1.0 / 4 / Area * (b1 * b3 + c1 * c3);
		G12 = 1.0 / 4 / Area * (b3 * b2 + c3 * c2);

		MyElem[i].Kc[0][0] = (float)(-G01 - G20);
		MyElem[i].Kc[0][1] = (float)G01;
		MyElem[i].Kc[0][2] = (float)G20;

		MyElem[i].Kc[1][0] = (float)G01;
		MyElem[i].Kc[1][1] = (float)(-G01 - G12);
		MyElem[i].Kc[1][2] = (float)G12;

		MyElem[i].Kc[2][0] = (float)G20;
		MyElem[i].Kc[2][1] = (float)G12;
		MyElem[i].Kc[2][2] = (float)(-G20 - G12);

		MyNode[I].Jext += (float)Area / 3.0*Jext;
		MyNode[J].Jext += (float)Area / 3.0*Jext;
		MyNode[K].Jext += (float)Area / 3.0*Jext;

	}

	for (i = 0; i < NumNodes; i++)
		if (MyNode[i].Type == 1)
			MyNode[i].Jext = 0.0;
	//Get Stiff and Dampinmg matrix and add to Node.Jext done

	free(X); free(Y);
}
void GPUMallocCopy()
{
	cudaMalloc((void**)&d_WannaUpdateVe, sizeof(int));
	cudaMalloc((void**)&d_NumElems, sizeof(int));
	cudaMalloc((void**)&d_NumNodes, sizeof(int));
	cudaMalloc((void**)&d_MyElem, NumElems * sizeof(FEMElem));
	cudaMalloc((void**)&d_MyNode, NumNodes * sizeof(FEMNode));
	cudaMemcpy(d_MyElem, MyElem, NumElems * sizeof(FEMElem), cudaMemcpyHostToDevice);
	cudaMemcpy(d_MyNode, MyNode, NumNodes * sizeof(FEMNode), cudaMemcpyHostToDevice);
	cudaMemcpy(d_NumNodes, &NumNodes, sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_NumElems, &NumElems, sizeof(int), cudaMemcpyHostToDevice);


	//message vecs between FEM and PCG
	cudaMalloc((void**)&d_p, NumNodes * sizeof(float));

	cudaMalloc((void**)&d_x, NumNodes * sizeof(float));
	float *hhh = (float *)malloc(NumNodes * sizeof(float));
	for (int ii = 0; ii < NumNodes; ii++)
		hhh[ii] = 0.0;
	cudaMemcpy(d_x, hhh, NumNodes * sizeof(float), cudaMemcpyHostToDevice);
	free(hhh);

	//message vecs between FEM and PCG

	//PCG vars
	//-------single  vars
	cudaMalloc((void**)&d_alpha, sizeof(float));
	cudaMalloc((void**)&d_beta, sizeof(float));
	cudaMalloc((void**)&d_pap, sizeof(float));
	cudaMalloc((void**)&d_rz1, sizeof(float));
	cudaMalloc((void**)&d_rz0, sizeof(float));
	//-------reduction buffers
	cudaMalloc((void**)&d_buffer1, NumNodes * sizeof(float));
	cudaMalloc((void**)&d_buffer2, NumNodes * sizeof(float));
	//-------vector var package
	cudaMalloc((void**)&d_PcgVec, NumNodes * sizeof(PcgVec));
	PcgVec *aaa = (PcgVec *)malloc(NumNodes * sizeof(PcgVec));
	for (int ii = 0; ii < NumNodes; ii++)
		aaa[ii].PreCondi = 1.0;
	cudaMemcpy(d_PcgVec, aaa, NumNodes * sizeof(PcgVec), cudaMemcpyHostToDevice);
	free(aaa);

	cudaMalloc((void**)&d_yinuonuo, 6000 * sizeof(float));
	cudaMalloc((void**)&d_anuonuo, NumNodes * sizeof(float));
}
void GPUFree()
{
	cudaFree(d_MyElem);
	cudaFree(d_MyNode);
	cudaFree(d_NumNodes);
	cudaFree(d_NumElems);
	cudaFree(d_WannaUpdateVe);
	cudaFree(d_x);
	cudaFree(d_p);

	//PCG vars
	cudaFree(d_x);
	cudaFree(d_p);
	cudaFree(d_buffer1);
	cudaFree(d_buffer2);
	cudaFree(d_alpha);
	cudaFree(d_beta);
	cudaFree(d_pap);
	cudaFree(d_rz1);
	cudaFree(d_rz0);
	cudaFree(d_PcgVec);
	cudaFree(d_yinuonuo);
	cudaFree(d_anuonuo);

}



//NDD matrix multi Funcs Vec A*x or A*p
__global__ void ElmentLocal_MtxA_Dot_Vector(FEMElem* d_MyElem, FEMNode* d_MyNode, int *d_NumElems, float *Vector)
{

	int e = threadIdx.x + blockIdx.x*blockDim.x;
	if (e >= *d_NumElems)
		return;

	int i, j;
	float ALocal[3];
	ALocal[0] = Vector[d_MyElem[e].I];
	ALocal[1] = Vector[d_MyElem[e].J];
	ALocal[2] = Vector[d_MyElem[e].K];
	for (i = 0; i < 3; i++)
	{
		d_MyElem[e].MatrixMultiVecBuff[i] = 0;
		for (j = 0; j < 3; j++)
			d_MyElem[e].MatrixMultiVecBuff[i] += (d_MyElem[e].Ve*d_MyElem[e].Kc[i][j] + d_MyElem[e].Dp[i][j]) * ALocal[j];
	}
}
__global__ void EachNodeExtractFromElemBuffInto_Ap(FEMElem* d_MyElem, FEMNode* d_MyNode, PcgVec *d_PcgVec, float *Vector, int *d_NumNodes)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	int e, LocalPos, i;
	float temp = 0;
	if (d_MyNode[nd].Type == 0)
	{
		temp = 0.0;
		for (i = 0; i < d_MyNode[nd].NumEle; i++)
		{
			e = d_MyNode[nd].EleID[i];
			LocalPos = d_MyNode[nd].EleOrd[i];
			temp += d_MyElem[e].MatrixMultiVecBuff[LocalPos];
		}
	}
	else
		temp = Vector[nd];
	d_PcgVec[nd].Ap = temp;

}

//math opertations to Vector and Number
__global__ void Num3isNum1dividedbyNum2(float *Num3, float *Num1, float *Num2)
{
	*Num3 = (*Num1) / (*Num2);
}
//-----------aim=sum_Vector
__global__ void reduction(float * Vector, float* aim, int *d_NumNodes) {
	__shared__ float sdata[5120];
	// each thread loads one element from global to shared mem
	unsigned int tid = threadIdx.x;

	for (int ii = 0; ii < 40; ii++)
		sdata[tid + 128 * ii] = 0.0;

	for (int ii = 0; ii < 40; ii++)
		if ((tid + 128 * ii) < (*d_NumNodes))
			sdata[tid + 128 * ii] = Vector[tid + 128 * ii];
	__syncthreads();

	for (int ii = 1; ii < 40; ii++)
		sdata[tid] += sdata[tid + 128 * ii];
	__syncthreads();
	if (tid < 32)
	{
		for (int ii = 1; ii < 4; ii++)
			sdata[tid] += sdata[tid + 32 * ii];
	}

	__syncthreads();
	if (tid == 0)
	{
		for (int ii = 1; ii < 32; ii++)
			sdata[0] += sdata[ii];
		*aim = sdata[0];
	}
}
__global__ void reduction_2X(float * Vector1, float* aim1, float * Vector2, float* aim2, int *d_NumNodes) {
	 __shared__ float sdata1[5120];
	 __shared__ float sdata2[5120];
	// each thread loads one element from global to shared mem
	unsigned int tid = threadIdx.x;

	for (int ii = 0; ii < 40; ii++)
	{
		sdata1[tid + 128 * ii] = 0.0;
		sdata2[tid + 128 * ii] = 0.0;
	}

	for (int ii = 0; ii < 40; ii++)
		if ((tid + 128 * ii) < (*d_NumNodes))
		{
			sdata1[tid + 128 * ii] = Vector1[tid + 128 * ii];
			sdata2[tid + 128 * ii] = Vector2[tid + 128 * ii];
		}
	__syncthreads();
	
	for (int ii = 1; ii < 40; ii++)
	{
		sdata1[tid] += sdata1[tid + 128 * ii];
		sdata2[tid] += sdata2[tid + 128 * ii];
	}
	__syncthreads();
	if (tid < 32)
	{
		for (int ii = 1; ii < 4; ii++)
		{
			sdata1[tid] += sdata1[tid + 32 * ii];
			sdata2[tid] += sdata2[tid + 32 * ii];
		}
	}

	__syncthreads();
	if (tid == 0)
	{
		for (int ii = 1; ii < 32; ii++)
		{
			sdata1[0] += sdata1[ii];
			sdata2[0] += sdata2[ii];
		}
		*aim1 = sdata1[0];
		*aim2 = sdata2[0];
	}
}


// FEM JExcitation->PcgVec.b
__global__ void UpdateExcitationSumTo_PCGSolver(FEMNode* d_MyNode, PcgVec *d_PcgVec, int *d_NumNodes)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;

	d_PcgVec[nd].b = d_MyNode[nd].Jext + d_MyNode[nd].TLM_Eq_Jsource + d_MyNode[nd].Jeddy;
}





float yinuonuo1[6000];
//PGC core Funcs
//----------make M=diag(A)^-1
__global__ void GPUmakePreconditoner(FEMElem* d_MyElem, FEMNode* d_MyNode, int *d_NumNodes, PcgVec *d_PcgVec)
{

	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	float temp;
	int e, LocalPos, i;
	if (d_MyNode[nd].Type == 0)
	{
		temp = 0.0;
		for (i = 0; i < d_MyNode[nd].NumEle; i++)
		{
			e = d_MyNode[nd].EleID[i];
			LocalPos = d_MyNode[nd].EleOrd[i];
			temp += (d_MyElem[e].Kc[LocalPos][LocalPos] * d_MyElem[e].Ve + d_MyElem[e].Dp[LocalPos][LocalPos]);
		}
		d_PcgVec[nd].PreCondi = 1.0 / temp;
	}

}
//----------initalize Pcg vectors based on x
__global__ void GPU_PCG_Prepare(int *d_NumNodes, float *d_p, PcgVec *d_PcgVec)
{

	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	float r0, z0;

	r0 = d_PcgVec[nd].b - d_PcgVec[nd].Ap;
	//r_old = (b - Ax);
	z0 = d_PcgVec[nd].PreCondi * r0;
	//z_old = M * r_old;
	d_PcgVec[nd].r0 = r0;
	//r_new = r_old;
	d_PcgVec[nd].r1 = r0;
	d_PcgVec[nd].z0 = z0;
	//z_new = z_old;
	d_PcgVec[nd].z1 = z0;
	//p = z_old;
	d_p[nd] = z0;

}
__global__ void Set0the_d_x(float *d_x, int *d_NumNodes)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	d_x[nd] = 0;
}

void CPU_PCG_Prepare()
{
	Set0the_d_x << <CudaBlckNum, CudaThrdNum >> > (d_x, d_NumNodes);
	GPUmakePreconditoner << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_MyNode, d_NumNodes, d_PcgVec);
	//PcgVec.Ap=Ax;
	ElmentLocal_MtxA_Dot_Vector << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_MyNode, d_NumElems, d_x);
	EachNodeExtractFromElemBuffInto_Ap << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_MyNode, d_PcgVec, d_x, d_NumNodes);
	GPU_PCG_Prepare << <CudaBlckNum, CudaThrdNum >> > (d_NumNodes, d_p, d_PcgVec);
}

//----------PCG middle processes. seprated into differnet functions because need sync
__global__ void pAp_into_Buffer1_And_r0z0_into_Buffer2(PcgVec *d_PcgVec, int *d_NumNodes, float *d_p, float *d_buffer1, float *d_buffer2)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	d_buffer1[nd] = (d_PcgVec[nd].Ap)*d_p[nd];
	d_buffer2[nd] = (d_PcgVec[nd].r0)*(d_PcgVec[nd].z0);
}
__global__ void r1z1_into_Buffer1(PcgVec *d_PcgVec, int *d_NumNodes, float *d_buffer1)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	d_buffer1[nd] = (d_PcgVec[nd].r1)*(d_PcgVec[nd].z1);
}
__global__ void GPU_PCG_Process1(PcgVec *d_PcgVec, float *d_p, float* d_x, float *d_alpha, int *d_NumNodes)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	float r1;
	//x = x + alpha * p;
	d_x[nd] += (*d_alpha)*d_p[nd];
	//r_new = r_old - alpha * A * p;
	r1 = d_PcgVec[nd].r0 - (*d_alpha)* d_PcgVec[nd].Ap;
	d_PcgVec[nd].r1 = r1;
	//z_new = M * r_new;
	d_PcgVec[nd].z1 = d_PcgVec[nd].PreCondi * r1;
}
__global__ void GPU_PCG_Process2(PcgVec *d_PcgVec, float *d_p, float* d_x, float *d_beta, int *d_NumNodes)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	//p = z_new + beta *p;
	d_p[nd] = d_PcgVec[nd].z1 + (*d_beta)*d_p[nd];
	//r_old = r_new;
	d_PcgVec[nd].r0 = d_PcgVec[nd].r1;
	//z_old= z_new;
	d_PcgVec[nd].z0 = d_PcgVec[nd].z1;
}
//---------PCGsolver call
void PCGSolve_Ax_b()
{
	//alpha = r_old.' * z_old/(p.' * A * p);
	//-----------PcgVec.Ap=A*p;
	ElmentLocal_MtxA_Dot_Vector << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_MyNode, d_NumElems, d_p);
	EachNodeExtractFromElemBuffInto_Ap << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_MyNode, d_PcgVec, d_p, d_NumNodes);
	//-----------buffer1=pAp buffer2=r0dotz0
	pAp_into_Buffer1_And_r0z0_into_Buffer2 << <CudaBlckNum, CudaThrdNum >> > (d_PcgVec, d_NumNodes, d_p, d_buffer1, d_buffer2);
	//-----------pap=sum buffer1 rz0==sum buffer2
	reduction_2X << <1, 128 >> > (d_buffer1, d_pap, d_buffer2, d_rz0, d_NumNodes);
	//-----------alpha=rz0/pap
	Num3isNum1dividedbyNum2 << <1, 1 >> > (d_alpha, d_rz0, d_pap);
	//cudaMemcpy(&yinuonuo1, d_alpha, sizeof(float), cudaMemcpyDeviceToHost);

	//x = x + alpha * p;
	//r_new = r_old - alpha * A * p;
	//z_new = M * r_new;
	GPU_PCG_Process1 << <CudaBlckNum, CudaThrdNum >> > (d_PcgVec, d_p, d_x, d_alpha, d_NumNodes);
	//cudaMemcpy(&yinuonuo1, d_x, NumNodes * sizeof(float), cudaMemcpyDeviceToHost);

	//beta = (z_new.'  * r_new)/(z_old.' * r_old);
	//-----------buffer1=z1*r1;
	r1z1_into_Buffer1 << <CudaBlckNum, CudaThrdNum >> > (d_PcgVec, d_NumNodes, d_buffer1);
	//-----------rz1=sum_buffer1
	reduction << <1, 128 >> > (d_buffer1, d_rz1, d_NumNodes);
	//-----------beta=rz1/rz0
	Num3isNum1dividedbyNum2 << <1, 1 >> > (d_beta, d_rz1, d_rz0);
	//cudaMemcpy(&yinuonuo1, d_beta, sizeof(float), cudaMemcpyDeviceToHost);

	//p = z_new + beta *p;
	//r_old = r_new;
	//z_old= z_new;
	GPU_PCG_Process2 << <CudaBlckNum, CudaThrdNum >> > (d_PcgVec, d_p, d_x, d_beta, d_NumNodes);
}

//TLM Funcs
__device__ void GaussElimSolve(float *X, float Mtrx[3][3], float *RHS)//Mtrx changed after use X=Mtrx\RHS;
{

	int i, j, k; float Aik; float S;
	for (k = 0; k < 2; k++)
	{
		if (!Mtrx[k][k])
			printf("Matrix is not good\n");
		for (i = k + 1; i < 3; i++)
		{
			Aik = Mtrx[i][k] / Mtrx[k][k];
			for (j = k; j < 3; j++)
			{
				Mtrx[i][j] = Mtrx[i][j] - Aik * Mtrx[k][j];
			}
			RHS[i] = RHS[i] - Aik * RHS[k];
		}
	}

	X[2] = RHS[2] / Mtrx[2][2];
	for (k = 1; k >= 0; k--)
	{
		S = RHS[k];
		for (j = k + 1; j < 3; j++)
		{
			S = S - Mtrx[k][j] * X[j];
		}
		X[k] = S / Mtrx[k][k];
	}



}
//-----------  Vi[3]---inject--->ElementReistors---reflect---->Vr[3] and return value: Ve at given Vi
__device__ float Newton_Raphson_Element(float *Vr, float *Vi, float Ve, float *Kc, float *B2terms)
//Vi is x0,y,z0 Vr is x y z // Kc+Ve restore G and Y //B2terms+X+Y restore B^2 Thus mag flux density B
{


	int count = 0, i, j;
	float dVr[3];
	float Y[3];
	float B = 0, B2;
	float dB2_dVr[3];
	float VeAtB, dVeAtB_dB2;

	float Jacobian[3][3];
	float Residual[3];

	if (Vi[0] == 0 && Vi[1] == 0 && Vi[2] == 0)
	{
		Vr[0] = 0; Vr[1] = 0; Vr[2] = 0;
		return Ve;
	}
	else
	{
		for (i = 0; i < 3; i++)
			Y[i] = Kc[i] * Ve;//ʵ���и����š�ֻ��Ϊ�˺�ԭ���ĳ��򱣳�һ��

		while (count < 10)
			//while(count<10)
		{
			//for (i = 0; i < 3; i++)
				//Vrtemp[i] = Vr[i];
			B2 = 0;
			for (i = 0; i < 3; i++)
				B2 += B2terms[i] * (Vi[i] + Vr[i])*(Vi[i] + Vr[i]);
			B = sqrt(B2);
			//---------------------------------------modify nonliear curve here below
			if (B < 0.6)
			{
				VeAtB = 2.0;
				dVeAtB_dB2 = 0.0;
			}
			else
			{
				VeAtB = 2.0 + 1e5 * (B - 0.6) * (B - 0.6)* (B - 0.6)/ B;
				dVeAtB_dB2 = (B * 3e5 * (B - 0.6) * (B - 0.6)- 1e5 * (B - 0.6)* (B - 0.6)* (B - 0.6)) / B / B / 2 / B;
			}
			//---------------------------------------modify nonliear curve here above

			//--------make residual and Jacobian matrix
			for (i = 0; i < 3; i++)
				Residual[i] = (VeAtB*Kc[i] * (Vi[i] + Vr[i]) - Y[i] * (Vi[i] - Vr[i]));

			for (i = 0; i < 3; i++)
				dB2_dVr[i] = (Vi[i] + Vr[i]) * 2 * B2terms[i];

			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					Jacobian[i][j] = -Kc[i] * (Vi[i] + Vr[i])*dVeAtB_dB2*dB2_dVr[j];
			for (i = 0; i < 3; i++)
				Jacobian[i][i] -= Kc[i] * VeAtB + Y[i];
			//--------solve for Vr step dVr=Jcb\Rsd
			GaussElimSolve(dVr, Jacobian, Residual);
			for (i = 0; i < 3; i++)
				Vr[i] += dVr[i];

			count++;
		}
	}

	//printf("%f %f\n", Vi[0], Vr[0]);
	return VeAtB;
}

__global__ void NonlinearElementsUpdate_Ve_And_EquvJSource(FEMElem* d_MyElem, int *d_NumElems, float *d_x, int *d_WannaUpdateVe)
{

	int e = threadIdx.x + blockIdx.x*blockDim.x;
	if (e >= *d_NumElems)
		return;
	//	if (d_MyElem[e].Type != 1)//nonliaer elements geometry group
		//	return;
		//Local Node ID
	int I = d_MyElem[e].I; int J = d_MyElem[e].J; int K = d_MyElem[e].K;

	float Vr[3];//copy to Local
	Vr[0] = d_MyElem[e].Vr[0];
	Vr[1] = d_MyElem[e].Vr[1];
	Vr[2] = d_MyElem[e].Vr[2];

	//Collect Inject wave amplitude from the global network (volatge)
	float Vi[3];
	Vi[0] = (d_x[J] - d_x[I] - Vr[0]);
	Vi[1] = (d_x[K] - d_x[J] - Vr[1]);
	Vi[2] = (d_x[I] - d_x[K] - Vr[2]);
	//NR solve reflection waves
	float Kc[3];// 3 Kc terms to calculate G func and Y value in the 3*3 nonliear matrix group
	Kc[0] = d_MyElem[e].Kc[0][1];// ----------note minus sign!!!!!!!!!!!
	Kc[1] = d_MyElem[e].Kc[1][2];
	Kc[2] = d_MyElem[e].Kc[2][0];
	float B2terms[3];// ---------restore B^2 from Vr Vi
	B2terms[0] = -Kc[0] / d_MyElem[e].Area;
	B2terms[1] = -Kc[1] / d_MyElem[e].Area;
	B2terms[2] = -Kc[2] / d_MyElem[e].Area;
	float Ve = Newton_Raphson_Element(Vr, Vi, d_MyElem[e].Ve, Kc, B2terms);// ----------elemental solve
				//for (int ii = 0; ii < 3; ii++)//linear test replcement for above func
				//Vr[ii] = -Vi[ii] / 3.0;
	//Update Ve here if in need
	//chagne reflected volage waves into eqiuvalent inject current sources on 3 element nodes to ground
	float temp1, temp2, temp3;
	temp1 = Kc[0] * d_MyElem[e].Ve * Vr[0] * 2;
	temp2 = Kc[1] * d_MyElem[e].Ve * Vr[1] * 2;
	temp3 = Kc[2] * d_MyElem[e].Ve * Vr[2] * 2;
	d_MyElem[e].Ieq[0] = temp1 - temp3;
	d_MyElem[e].Ieq[1] = temp2 - temp1;
	d_MyElem[e].Ieq[2] = temp3 - temp2;
	d_MyElem[e].Vr[0] = Vr[0];
	d_MyElem[e].Vr[1] = Vr[1];
	d_MyElem[e].Vr[2] = Vr[2];
	if (*d_WannaUpdateVe == 1)
		d_MyElem[e].Ve = Ve;
}

__global__ void EachNodeExtractFromElemIeqInto_NodeIeq(FEMElem* d_MyElem, FEMNode* d_MyNode, int *d_NumNodes)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	int e, LocalPos, i;
	float temp = 0;
	if (d_MyNode[nd].Type == 0)
	{
		temp = 0.0;
		for (i = 0; i < d_MyNode[nd].NumEle; i++)
		{
			e = d_MyNode[nd].EleID[i];
			LocalPos = d_MyNode[nd].EleOrd[i];
			temp += d_MyElem[e].Ieq[LocalPos];
		}
	}
	else
		temp = 0.0;
	d_MyNode[nd].TLM_Eq_Jsource = temp;

}





PcgVec PcgV[320];

void TLM_MatrixFree_Solve_for_this_timepoint()
{
	int wannaupdateve = 0;
	cudaMemcpy(d_WannaUpdateVe, &wannaupdateve, sizeof(int), cudaMemcpyHostToDevice);
	for (int TLMct = 0; TLMct < 6; TLMct++)
	{
		if (TLMct == 5)
		{
			wannaupdateve = 0;
			cudaMemcpy(d_WannaUpdateVe, &wannaupdateve, sizeof(int), cudaMemcpyHostToDevice);
		}
		UpdateExcitationSumTo_PCGSolver << <CudaBlckNum, CudaThrdNum >> > (d_MyNode, d_PcgVec, d_NumNodes);
		cudaMemcpy(&PcgV, d_PcgVec, 320 * sizeof(PcgVec), cudaMemcpyDeviceToHost);

		CPU_PCG_Prepare();
		for (int pcgit = 0; pcgit < 40; pcgit++)
			PCGSolve_Ax_b();
		NonlinearElementsUpdate_Ve_And_EquvJSource << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_NumElems, d_x, d_WannaUpdateVe);
		EachNodeExtractFromElemIeqInto_NodeIeq << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_MyNode, d_NumNodes);
		cudaMemcpy(&MyNode, d_MyNode, NumNodes * sizeof(FEMNode), cudaMemcpyDeviceToHost);
		cudaMemcpy(&yinuonuo1, d_x, NumNodes * sizeof(float), cudaMemcpyDeviceToHost);
		//printf("%f\n", yinuonuo1[2]);
		//		for (int i = 0; i < 320; i++)
		//		printf("  %f", MyNode[i].TLM_Eq_Jsource);
	}
}


//Jeddy=Dp*A;
__global__ void ElmentLocal_MtxDp_Dot_d_x(FEMElem* d_MyElem, FEMNode* d_MyNode, int *d_NumElems, float *d_x)
{

	int e = threadIdx.x + blockIdx.x*blockDim.x;
	if (e >= *d_NumElems)
		return;

	int i, j;
	float ALocal[3];
	ALocal[0] = d_x[d_MyElem[e].I];
	ALocal[1] = d_x[d_MyElem[e].J];
	ALocal[2] = d_x[d_MyElem[e].K];
	for (i = 0; i < 3; i++)
	{
		d_MyElem[e].MatrixMultiVecBuff[i] = 0;
		for (j = 0; j < 3; j++)
			d_MyElem[e].MatrixMultiVecBuff[i] += (d_MyElem[e].Dp[i][j]) * ALocal[j];
	}
}
__global__ void EachNodeExtractFromElemBuffInto_Jeddy(FEMElem* d_MyElem, FEMNode* d_MyNode, float *d_x, int *d_NumNodes)
{
	int nd = threadIdx.x + blockIdx.x*blockDim.x;
	if (nd >= *d_NumNodes)
		return;
	int e, LocalPos, i;
	float temp = 0;
	if (d_MyNode[nd].Type == 0)
	{
		temp = 0.0;
		for (i = 0; i < d_MyNode[nd].NumEle; i++)
		{
			e = d_MyNode[nd].EleID[i];
			LocalPos = d_MyNode[nd].EleOrd[i];
			temp += d_MyElem[e].MatrixMultiVecBuff[LocalPos];
		}
	}
	else
		temp = d_x[nd];
	d_MyNode[nd].Jeddy = temp;

}



int main()
{

	LoadMeshInfoAndPrepareFEM();
	GPUMallocCopy();

	cudaEvent_t start, stop;//unit: ms

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	for (int t = 0; t < 10; t++)
	{
		ElmentLocal_MtxDp_Dot_d_x << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_MyNode, d_NumElems, d_x);
		EachNodeExtractFromElemBuffInto_Jeddy << <CudaBlckNum, CudaThrdNum >> > (d_MyElem, d_MyNode, d_x, d_NumNodes);

		TLM_MatrixFree_Solve_for_this_timepoint();
	}


	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	


	cudaMemcpy(&yinuonuo1, d_x, NumNodes * sizeof(float), cudaMemcpyDeviceToHost);
	// for (int i = 0; i < 320; i++)
	// printf("%f \n", yinuonuo1[i]);


	printf("\nElapsed time for NDDR Iteration is %.4f ms\n", elapsedTime);



	GPUFree();



}