#include <bits/stdc++.h>
#include <ctime>

void CPUcompute(float* grid, float* retGrid, int N){
  for (int i = 0; i < N; i++){
    retGrid[i] = grid[i] / 0.02;
  }
}

__global__ void cudaKernel(float* grid, float* retGrid, int N){
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  retGrid[i] = grid[i] / 0.02;
}

int main(){
  int N = 13 * 13 * 13;

  float* A = new float[N];
  float* B = new float[N];
  for(int i = 0; i < N; i++){
    A[i] = i;
    B[i] = i;
  }

  float* dA = NULL;
  float* dB = NULL;
  cudaMalloc(&dA, sizeof(float)*N);
  cudaMalloc(&dB, sizeof(float)*N);
  cudaMemcpy(dA, A, sizeof(float)*N, cudaMemcpyHostToDevice);
  cudaMemcpy(dB, B, sizeof(float)*N, cudaMemcpyHostToDevice);


  clock_t cpuStart = clock();
  CPUcompute(A, B, N);
  printf("CPU time: %f\n", (float)(clock()-cpuStart)/CLOCKS_PER_SEC);

  clock_t kernelStart = clock();
  cudaKernel<<<N/1000, 1000>>>(dA, dB, N);
  printf("Cuda time: %f\n", (float)(clock()-kernelStart)/CLOCKS_PER_SEC);
  return 0;
}
