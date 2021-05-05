
#include <stdio.h>
#include <cuda_runtime.h>


__global__ void MatrixMut(int *A, int *B, int *C, int N)
{

	__shared__ int As[256];
	__shared__ int Bs[256];


int tx=threadIdx.x, ty=threadIdx.y;

int row=blockDim.x*blockIdx.x+threadIdx.x;

int col=blockDim.y*blockIdx.y+threadIdx.y;


int temp=0;
// 16*16 subblock every time, 64 times, for each time
// each thread get 1 element, a block get 256 elements, then sysnchronize
// loop 16 multiply and summation, complete 1/64. Sync for next subblock
for(int k=0;k<64;k++)
{
	As[tx*16+ty]=A[row*N+k*16+ty];
	Bs[tx*16+ty]=B[col+(k*16+tx)*N];
	__syncthreads();

    for(int x=0;x<16;x++)
    temp+=As[tx*16+x]*Bs[x*16+ty];
	__syncthreads();
}


C[row*N+col]=temp;

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




dim3 grid_size(64,64);
dim3 block_size(16,16);

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
		printf("Test Failed!\n");
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



