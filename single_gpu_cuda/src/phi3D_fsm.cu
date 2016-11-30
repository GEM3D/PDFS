/**
 * @file phi3D_fsm.c
 * @brief      Source file for 3D Phi Function that implements the parallel fast
 *             sweeping method for solving the Eikonal equation in CUDA The
 *             algorithm implemented for parallel fast sweeping method is from a
 *             paper in the Journal of Computational Physics titled "A parallel
 *             fast sweeping method for the Eikonal Equation" by Miles Detrixhe,
 *             Federic Gibou, and Chohong Min.
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
 *
 *
 *
 * @see        http://www.sciencedirect.com/science/article/pii/S002199911200722X
 * 
 * Copyright (c) 2016
 * Mechanical and Bio-medical Engineering Department
 * Boise State University
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "phi3D.h"
#include "phi3D_fsm.h"

#define max(a, b) ((a > b) ? a : b)
#define min(a, b) ((a < b) ? a : b)

// Private method definitions
static void fast_sweep(Phi *p, int itr, cudaPitchedPtr dPitchPtr);
static void _cudaMemcpy3D(cudaPitchedPtr src, cudaPitchedPtr dst,
                          cudaExtent dExt, cudaMemcpyKind kind);
static int iDivUp(int a, int b);

// CUDA functions
__global__ void fast_sweep_kernel(cudaPitchedPtr dPitchPtr, SweepInfo s);
__device__ double solve_eikonal(double cur_dist, double minX, double minY,
                                double minZ, double dx, double dy, double dz);

/**
 * @brief         Calls the fast sweeping method a number of times specified by
 *                the iterations argument.
 *
 * @param[in,out] pf          Pointer to phi function.
 * @param[in]     iterations  Max iterations.
 */
void run_fsm(Phi *pf, int iterations) {

  int max_x = pf->x + 2;
  int max_y = pf->y + 2;
  int max_z = pf->z + 2;

  /*********************** CUDA ***********************/

  // time cuda code
  cudaEvent_t start, stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  cudaPitchedPtr hostPtr =
      make_cudaPitchedPtr(pf->distance, max_x * sizeof(double), max_x, max_y);

  cudaPitchedPtr devicePtr;
  cudaExtent dExt = make_cudaExtent(max_x * sizeof(double), max_y, max_z);
  cudaCheckError();

  cudaMalloc3D(&devicePtr, dExt);
  cudaCheckError();

  _cudaMemcpy3D(hostPtr, devicePtr, dExt, cudaMemcpyHostToDevice);
  cudaCheckError();

  fast_sweep(pf, iterations, devicePtr);

  _cudaMemcpy3D(devicePtr, hostPtr, dExt, cudaMemcpyDeviceToHost);
  cudaCheckError();

  cudaFree(devicePtr.ptr);
  cudaFree(hostPtr.ptr);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);

  printf("Parallel FSM time: %f s.\n", elapsedTime / 1000.0);
  /****************************************************/
}

/**
 * @brief         Calculates the distance field for a 3D grid by solving the
 *                Eikonal equation at each grid point using the parallel Fast
 *                Sweeping Method.
 *
 *                Sweeping Directions:
 *                (1) i = 1:I, j = 1:J, k = 1:K
 *                (2) i = I:1, j = 1:J, k = K:1
 *                (3) i = I:1, j = 1:J, k = 1:K
 *                (4) i = 1:I, j = 1:J, k = K:1
 *                (5) i = I:1, j = J:1, k = K:1
 *                (6) i = 1:I, j = J:1, k = 1:K
 *                (7) i = 1:I, j = J:1, k = K:1
 *                (8) i = I:1, j = J:1, k = 1:K
 *
 * @param[in.out] p          Pointer to phi function.
 * @param[in]     itr        Max iterations.
 * @param[in,out] dPitchPtr  Pointer to distance array in device memory.
 */
static void fast_sweep(Phi *p, int itr, cudaPitchedPtr dPitchPtr) {

  // Information regarding sweeping and linear indexing
  int meshDim = 3;

  SweepInfo sw;
  sw.xDim = p->x;
  sw.dx = p->dx;
  sw.yDim = p->y;
  sw.dy = p->dy;
  sw.zDim = p->z;
  sw.dz = p->dz;

  int totalLevels = sw.xDim + sw.yDim + sw.zDim;

  // loop till the number of times to sweep
  int loop = 1;
  while (loop <= itr) {

    printf("Please wait. Sweeping...[%d/%d]\n", loop, itr);

    for (int swCount = 1; swCount <= 8; ++swCount) {
      int start = (swCount == 2 || swCount == 5 || swCount == 7 || swCount == 8)
                      ? totalLevels
                      : meshDim;
      int end = (start == meshDim) ? totalLevels + 1 : meshDim - 1;
      int incr = (start == meshDim) ? true : false;

      // sweep offset is used for translating the 3D coordinates
      // to perform sweeps from different directions
      sw.xSweepOff = (swCount == 4 || swCount == 8) ? sw.xDim + 1 : 0;
      sw.ySweepOff = (swCount == 2 || swCount == 6) ? sw.yDim + 1 : 0;
      sw.zSweepOff = (swCount == 3 || swCount == 7) ? sw.zDim + 1 : 0;

      for (int level = start; level != end;
           level = (incr) ? level + 1 : level - 1) {
        int xs = max(1, level - (sw.yDim + sw.zDim)),
            ys = max(1, level - (sw.xDim + sw.zDim));
        int xe = min(sw.xDim, level - (meshDim - 1)),
            ye = min(sw.yDim, level - (meshDim - 1));

        int xr = xe - xs + 1, yr = ye - ys + 1;
        int tth = xr * yr; // Total number of threads needed

        dim3 bs(16, 16, 1);
        if (tth < 256) {
          bs.x = xr;
          bs.y = yr;
        }
        dim3 gs(iDivUp(xr, bs.x), iDivUp(yr, bs.y), 1);

        sw.level = level;
        sw.xOffSet = xs;
        sw.yOffset = ys;

        fast_sweep_kernel << <gs, bs>>> (dPitchPtr, sw);
        cudaThreadSynchronize();
        cudaCheckError();
      }
    }
    printf("Sweeping finished!......[%d/%d]\n", loop, itr);
    ++loop;
  }
}


/**
 * @brief         CUDA kernel for the fast sweeping method.
 *
 * @param[in,out] dPitchPtr  Pointer to distance array in device memory.
 * @param[in]     s          Sweep information.
 */
__global__ void fast_sweep_kernel(cudaPitchedPtr dPitchPtr, SweepInfo s) {
  int x = (blockIdx.x * blockDim.x + threadIdx.x) + s.xOffSet;
  int y = (blockIdx.y * blockDim.y + threadIdx.y) + s.yOffset;
  if (x <= s.xDim && y <= s.yDim) {
    int z = s.level - (x + y);
    if (z > 0 && z <= s.zDim) {
      int i = abs(z - s.zSweepOff);
      int j = abs(y - s.ySweepOff);
      int k = abs(x - s.xSweepOff);

      char *devPtr = (char *)dPitchPtr.ptr;
      size_t pitch = dPitchPtr.pitch;
      size_t slicePitch = pitch * (s.yDim + 2);

      double *c_row =
          (double *)((devPtr + i * slicePitch) + j * pitch); // center row
      double center = c_row[k];                              // center distance
      double left = c_row[k - 1];                            // left distance
      double right = c_row[k + 1];                           // right distance
      double up = ((double *)((devPtr + i * slicePitch) +
                              (j - 1) * pitch))[k]; // upper distance
      double down = ((double *)((devPtr + i * slicePitch) +
                                (j + 1) * pitch))[k]; // lower distance
      double front = ((double *)((devPtr + (i - 1) * slicePitch) +
                                 j * pitch))[k]; // front distance
      double back = ((double *)((devPtr + (i + 1) * slicePitch) +
                                j * pitch))[k]; // back distance

      double minX = min(left, right);
      double minY = min(up, down);
      double minZ = min(front, back);
      c_row[k] = solve_eikonal(center, minX, minY, minZ, s.dx, s.dy, s.dz);
    }
  }
}

/**
 * @brief      Solves Eikonal equation at linearized 3D index. Returns the
 *             minimum of calculated and old distance values.
 *
 * @param[in]  cur_dist  Current distance value.
 * @param[in]  minX      Minimum distance in the x-direction.
 * @param[in]  minY      Minimum distance in the y-direction.
 * @param[in]  minZ      Minimum distance in the z-direction.
 * @param[in]  dx        Spacing in the x-direction.
 * @param[in]  dy        Spacing in the y-direction.
 * @param[in]  dz        Spacing in the z-direction.
 *
 * @return     Minimum value of the solution at given index.
 */
__device__ double solve_eikonal(double cur_dist, double minX, double minY,
                                double minZ, double dx, double dy, double dz) {
  double dist_new = 0;
  double m[] = { minX, minY, minZ };
  double d[] = { dx, dy, dz };

  // sort the mins
  for (int i = 1; i < 3; i++) {
    for (int j = 0; j < 3 - i; j++) {
      if (m[j] > m[j + 1]) {
        double tmp_m = m[j];
        double tmp_d = d[j];
        m[j] = m[j + 1];
        d[j] = d[j + 1];
        m[j + 1] = tmp_m;
        d[j + 1] = tmp_d;
      }
    }
  }

  // simplifying the variables
  double m_0 = m[0], m_1 = m[1], m_2 = m[2];
  double d_0 = d[0], d_1 = d[1], d_2 = d[2];
  double m2_0 = m_0 * m_0, m2_1 = m_1 * m_1, m2_2 = m_2 * m_2;
  double d2_0 = d_0 * d_0, d2_1 = d_1 * d_1, d2_2 = d_2 * d_2;

  dist_new = m_0 + d_0;
  if (dist_new > m_1) {

    double s = sqrt(-m2_0 + 2 * m_0 * m_1 - m2_1 + d2_0 + d2_1);
    dist_new = (m_1 * d2_0 + m_0 * d2_1 + d_0 * d_1 * s) / (d2_0 + d2_1);

    if (dist_new > m_2) {

      double a =
          sqrt(-m2_0 * d2_1 - m2_0 * d2_2 + 2 * m_0 * m_1 * d2_2 - m2_1 * d2_0 -
               m2_1 * d2_2 + 2 * m_0 * m_2 * d2_1 - m2_2 * d2_0 - m2_2 * d2_1 +
               2 * m_1 * m_2 * d2_0 + d2_0 * d2_1 + d2_0 * d2_2 + d2_1 * d2_2);

      dist_new = (m_2 * d2_0 * d2_1 + m_1 * d2_0 * d2_2 + m_0 * d2_1 * d2_2 +
                  d_0 * d_1 * d_2 * a) /
                 (d2_0 * d2_1 + d2_0 * d2_2 + d2_1 * d2_2);
    }
  }

  return min(cur_dist, dist_new);
}

/*
 * Copies 3D memory from host to device and device to host.
 *
 * Arguments:
 *   cudaPitchedPtr [in]  - pointer to distance array
 *   cudaPitchedPtr [out] - pointer to distance array
 *   cudaMemcpyKind [in]  - specifies the direction of copy
 * Returns:
 *
 */


/**
 * @brief      Copies 3D memory from host to device and device to host.
 *
 * @param[in]  src   Pointer to source distance array.
 * @param[out] dst   Pointer to destination disance array
 * @param[in]  dExt  Cuda extent.
 * @param[in]  kind  Specifies the direction of copy.
 */
static void _cudaMemcpy3D(cudaPitchedPtr src, cudaPitchedPtr dst,
                          cudaExtent dExt, cudaMemcpyKind kind) {
  cudaMemcpy3DParms mcp = { 0 };

  mcp.kind = kind;
  mcp.extent = dExt;

  mcp.srcPtr = src;
  mcp.dstPtr = dst;

  cudaMemcpy3D(&mcp);
  cudaCheckError();
}

/**
 * @brief      Calculates number of threads in each dimension of a thread block
 *
 * @param[in]  a
 * @param[in]  b
 *
 * @return     { description_of_the_return_value }
 */
static int iDivUp(int a, int b) {
  return ((a % b) != 0) ? (a / b + 1) : (a / b);
}
