
#include <stdio.h>
#include <cuda_runtime.h>


__global__ void MatrixMut(int *A, int *B, int *C, int N)
{

int i=threadIdx.x;
int j=blockIdx.x;
int temp=0;


for(int k=0;k<N;k++)
temp+= A[i*N+k]*B[k*N+j];



C[i*N+j]=temp;

}


int main()
{

const int dim=1<<10;
const int size=dim*dim;

int *A= (int *)malloc(size*sizeof(int));
int *B= (int *)malloc(size*sizeof(int));
int *C= (int *)malloc(size*sizeof(int));


for (int i = 0; i < size; ++i)
{
	
	A[i]=1;
	B[i]=2;

}


int *d_A=NULL;
int *d_B=NULL;
int *d_C=NULL;

cudaMalloc((void**) &d_A,size*sizeof(int));
cudaMalloc((void**) &d_B,size*sizeof(int));
cudaMalloc((void**) &d_C,size*sizeof(int));



cudaMemcpy(d_A,A,size*sizeof(int),cudaMemcpyHostToDevice);
cudaMemcpy(d_B,B,size*sizeof(int),cudaMemcpyHostToDevice);




dim3 grid_size(dim);
dim3 block_size(dim);

cudaEvent_t start, stop;


cudaEventCreate(&start);
cudaEventCreate(&stop);
cudaEventRecord(start, 0);


MatrixMut<<<grid_size,block_size>>>(d_A,d_B,d_C,dim);


cudaEventRecord(stop, 0);
cudaEventSynchronize(stop);

float msecTotal = 0.0f;
cudaEventElapsedTime(&msecTotal, start, stop);

printf("Eclapsed time is %f ms \n", msecTotal);

cudaMemcpy(C,d_C,size*sizeof(int),cudaMemcpyDeviceToHost);




for (int i = 0; i < size; ++i)
{
	
	if(C[i]!=dim*2)
	{
		printf("A[i]= %d, B[i]=%d, C[i]=%d \n", A[i], B[i], C[i]);
		exit(-1);
	}

}

printf("Test Passed \n");

free(A);
free(B);
free(C);

cudaFree(d_A);
cudaFree(d_B);
cudaFree(d_C);

return 0;

}



