#include <stdio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define GET_TIME(now){ \
	double t; \
	t=clock(); \
	now =t/CLOCKS_PER_SEC ;\
}
#define N 40000
__global__ void TreeSum(int*c, int*d);
__global__ void Multiply(int*a, int* b, int* c);
__global__ void LoopSum(int*c, int*d);
void printDevProp(cudaDeviceProp devProp);
int main()
{
	// Number of CUDA devices
	int devCount;
	double T_start, T_end;;
	cudaGetDeviceCount(&devCount);
	printf("CUDA Device Query...\n");
	printf("There are %d CUDA devices.\n", devCount);

	// Iterate through devices
	for (int i = 0; i < devCount; ++i)

	{
		// Get device properties
		printf("\nCUDA Device #%d\n", i);
		cudaDeviceProp devProp;
		cudaGetDeviceProperties(&devProp, i);
		printDevProp(devProp);
	}
	cudaSetDevice(0);
	int i, a[N], b[N], c[N], d;
	int *dev_a, *dev_b, *dev_c, *dev_d;

	cudaMalloc((void**)&dev_a, N * sizeof(int));
	cudaMalloc((void**)&dev_b, N * sizeof(int));
	cudaMalloc((void**)&dev_c, N * sizeof(int));
	cudaMalloc((void**)&dev_d, 1 * sizeof(int));
	for (i = 0; i<N; i++)
	{
		a[i] = i;
		b[i] = 1;
	}
	cudaMemcpy(dev_a, a, N * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_b, b, N * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_c, c, N * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_d, &d,  sizeof(int), cudaMemcpyHostToDevice);

	Multiply << <200, 200 >> >(dev_a, dev_b, dev_c);


	GET_TIME(T_start);

	TreeSum << <1, 1 >> >(dev_c, dev_d);
	cudaDeviceSynchronize();

	GET_TIME(T_end);

	cudaMemcpy(&d, dev_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);
	printf("\nResult of TreeAdder is %d\n", d);
	printf("Elapsed time for TreeAdder is %.4f s\n", (T_end - T_start));



	GET_TIME(T_start);

	LoopSum << <1, 1>> >(dev_c, dev_d);
	cudaDeviceSynchronize();

	GET_TIME(T_end);
	cudaMemcpy(&d, dev_d, 1 * sizeof(int), cudaMemcpyDeviceToHost);
	printf("\nResult of LoopAdder is %d\n", d);
	printf("Elapsed time for LoopAdder is %.4f s\n", (T_end - T_start));
	

	cudaFree(dev_a);
	cudaFree(dev_b);
	cudaFree(dev_c);
	cudaFree(dev_d);

	return 0;
}

__global__ void Multiply(int*a, int* b, int* c)
{
	int bid = threadIdx.x + blockIdx.x*blockDim.x;
	if (bid<N)
		c[bid] = a[bid] * b[bid];
}

__global__ void LoopSum(int*c, int*d)
{
	int i = 0;
	int temp = 0;
	for (i = 0; i < N; i++)
		temp += c[i];

	*d = temp;
}

__global__ void TreeSum(int*c, int*d)
{
	int i = 0;
	int temp = 0;
	for (i = 0; i < N; i++)
		temp += c[i];

	*d = temp;
}

// Print device properties
void printDevProp(cudaDeviceProp devProp)
{
	printf("Major revision number:         %d\n", devProp.major);
	printf("Minor revision number:         %d\n", devProp.minor);
	printf("Name:                          %s\n", devProp.name);
	printf("Total global memory:           %u\n", devProp.totalGlobalMem);
	printf("Total shared memory per block: %u\n", devProp.sharedMemPerBlock);
	printf("Total registers per block:     %d\n", devProp.regsPerBlock);
	printf("Warp size:                     %d\n", devProp.warpSize);
	printf("Maximum memory pitch:          %u\n", devProp.memPitch);
	printf("Maximum threads per block:     %d\n", devProp.maxThreadsPerBlock);
	for (int i = 0; i < 3; ++i)
		printf("Maximum dimension %d of block:  %d\n", i, devProp.maxThreadsDim[i]);
	for (int i = 0; i < 3; ++i)
		printf("Maximum dimension %d of grid:   %d\n", i, devProp.maxGridSize[i]);
	printf("Clock rate:                    %d\n", devProp.clockRate);
	printf("Total constant memory:         %u\n", devProp.totalConstMem);
	printf("Texture alignment:             %u\n", devProp.textureAlignment);
	printf("Concurrent copy and execution: %s\n", (devProp.deviceOverlap ? "Yes" : "No"));
	printf("Number of multiprocessors:     %d\n", devProp.multiProcessorCount);
	printf("Kernel execution timeout:      %s\n", (devProp.kernelExecTimeoutEnabled ? "Yes" : "No"));
	return;
}