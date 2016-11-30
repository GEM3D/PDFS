/**
 * @file phi3D_fsm.h
 * @brief      Header file for 3D Phi Function that implements the parallel fast
 *             sweeping method for solving the Eikonal equation in CUDA
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
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

#ifndef PHI3D_FSM_H
#define PHI3D_FSM_H

/**
 * @brief      Macro for checking CUDA errors following a CUDA launch or API
 *             call
 */
#define cudaCheckError()                                                                           \
  {                                                                                              \
    cudaError_t e = cudaGetLastError();                                                        \
    if (e != cudaSuccess) {                                                                    \
      printf("\nCuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e));     \
      exit(EXIT_FAILURE);                                                                    \
    }                                                                                          \
  }

/**
 * @brief      Structure for storing information used during sweeping to manage
 *             internal grid dimensions, sweep directions, position of the node
 *             on the array and its offset in each kernel block
 */
typedef struct {
  int    level;
  int    xDim;
  int    yDim;
  int    zDim;
  int    xOffSet;
  int    yOffset;
  int    xSweepOff;
  int    ySweepOff;
  int    zSweepOff;
  double dx;
  double dy;
  double dz;
} SweepInfo;

// External Linkage Method Definition
extern "C" {
void run_fsm(Phi *, int);
}

#endif /* PHI3D_FSM_H */
