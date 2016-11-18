/**
 * @file block_fsm.h
 * @brief      Header file for 3D Block Fast Sweeping Method (FSM).
 *
 *             The basic parallel algorithm for FSM is taken from the paper
 *             published in the Journal of Computational Physics titled
 *             "A parallel fast sweeping method for the Eikonal equation" by
 *             Miles Detrixhe, Deferic Gibou, and Chohong Min.
 *
 *             The functions within this file is responsible for dividing the
 *             grid into manageable blocks forming a coarser grid. The coarser
 *             grid is traversed and executed parallely similar to the Detrixhe
 *             et al. algorithm for FSM using MPI processes. For each sweep
 *             iteration blocks are retrieved from the NetCDF file, the sweep
 *             update is performed and the updated values are written back to
 *             the NetCDF file.
 *
 * @author     Shrestha, Anup
 * @date       18 JUL 2016
 *
 * @see        http://www.sciencedirect.com/science/article/pii/S002199911200722X
 */
#ifndef BLOCK_FSM_H
#define BLOCK_FSM_H

#include <stddef.h>

/**
 * @brief      Starts the parallel block FSM routine.
 *
 * @param[in]  ncid      NetCDF ID.
 * @param[in]  blockDim  The dimensions of the block grid.
 */
void startBlockFSM(int ncid, size_t *blockDim);

#endif /* BLOCK_FSM_H */
