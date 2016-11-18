/**
 * @file phi3D_fsm.h
 * @brief      Header file for 3D Phi Function that implements the parallel fast
 *             sweeping method for solving the Eikonal equation in CUDA
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
 */

#ifndef PHI3D_FSM_H
#define PHI3D_FSM_H

/**
 * @brief      Macro for checking CUDA errors following a CUDA launch or API
 *             call
 */
#define cudaCheckError()                                                                     \
  {                                                                                          \
    cudaError_t e = cudaGetLastError();                                                      \
    if (e != cudaSuccess) {                                                                  \
      printf("\nCuda failure %s:%d: '%s'\n", __FILE__, __LINE__, cudaGetErrorString(e));     \
      exit(EXIT_FAILURE);                                                                    \
    }                                                                                        \
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
