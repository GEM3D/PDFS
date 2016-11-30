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
