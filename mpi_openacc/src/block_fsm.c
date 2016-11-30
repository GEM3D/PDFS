/**
 * @file block_fsm.c
 * @brief      Source file for 3D Block Fast Sweeping Method (FSM).
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

#include "block_fsm.h"
#include "macrologger.h"
#include "ncwrap.h"
#include "phi3D_fsm.h"
#include "utilities.h"

#include <math.h>
#include <mpi.h>
#include <openacc.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/* MPI routine variables */
static int mpiSize;
static int mpiRank;
/*-----------------------------------------*/

/* Dimension of decomposed blocks */
static int BOX_X;
static int BOX_Y;
static int BOX_Z;

static size_t *blockDim;
static size_t  count[NDIMS];
static int     x, y, z;
/*-----------------------------------------*/

static void runFSMAndPropagate(Phi *pf);
static void receive(Phi *pf, int rx, int ry, int rz, int pCount, int tag, int xs, int xe, int ys,
                    int ye, int zs, int ze);
static void send(Phi *pf, int sx, int sy, int sz, int pCount, int tag, int xs, int xe, int ys,
                 int ye, int zs, int ze);
static void receiveI(Phi *pf, int rx, int ry, int rz, int pCount, int xs, int xe, int ys, int ye,
                     int zs, int ze, MPI_Request req);
static void sendI(Phi *pf, int sx, int sy, int sz, int pCount, int xs, int xe, int ys, int ye,
                  int zs, int ze, double *plane, MPI_Request req);

/**
 * @brief      Convert a 1-Dimensional 0-indexed coordinate to 3-Dimensional
 *             1-indexed coordinate.
 *
 * @param[in]  i       1-Dimensional index.
 * @param[in]  width   Width of the grid.
 * @param[in]  height  Height of the grid.
 * @param[out] x       x-coordinate.
 * @param[out] y       y-coordinate.
 * @param[out] z       z-coordinate
 */
static void linearTo3D(int i, int width, int height, int *x, int *y, int *z)
{
  *x = (i % width) + 1;
  *y = ((i / width) % height) + 1;
  *z = (i / (width * height)) + 1;
}

/**
 * @brief      Convert a 3-Dimensional 1-indexed coordinate to 1-Dimensional
 *             0-indexed coordinate.
 *
 * @param[in]  width   Width of the grid.
 * @param[in]  height  Height of the grid.
 * @param[in]  x       x-coordinate.
 * @param[in]  y       y-coordinate.
 * @param[in]  z       z-coordinate
 *
 * @return     1-Dimensional index starting at position 0.
 */
static int linearFrom3D(int width, int height, int x, int y, int z)
{
  return (x - 1) + width * (y - 1) + width * height * (z - 1);
}

/**
 * @brief      Starts a block fsm.
 *
 * @param[in]  ncid    The ncid
 * @param      blkDim  The block dim
 */
void startBlockFSM(int ncid, size_t *blkDim)
{
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  blockDim = blkDim;
  linearTo3D(mpiRank, blockDim[0], blockDim[1], &x, &y, &z);
  LOG_DEBUG("[%d/%d] Block start coordinate: (%d, %d, %d)", mpiRank, mpiSize, x, y, z);

  size_t dim[NDIMS];
  double dlt[NDIMS];
  ncGetDim(ncid, dim, dlt);

  BOX_X = (blockDim[0] == 1) ? dim[0] : (dim[0] / blockDim[0]) + 1;
  BOX_Y = (blockDim[1] == 1) ? dim[1] : (dim[1] / blockDim[1]) + 1;
  BOX_Z = (blockDim[2] == 1) ? dim[2] : (dim[2] / blockDim[2]) + 1;

  size_t start[NDIMS];
  start[2] = (x - 1) * BOX_X;
  start[1] = (y - 1) * BOX_Y;
  start[0] = (z - 1) * BOX_Z;

  count[2] = min(dim[0] - 1, start[2] + BOX_X) - start[2];
  count[1] = min(dim[1] - 1, start[1] + BOX_Y) - start[1];
  count[0] = min(dim[2] - 1, start[0] + BOX_Z) - start[0];

  start[2] = (start[2] == 0) ? 0 : start[2] - 1;
  start[1] = (start[1] == 0) ? 0 : start[1] - 1;
  start[0] = (start[0] == 0) ? 0 : start[0] - 1;

  count[2] = (start[2] == 0) ? count[2] + 1 : count[2] + 2;
  count[1] = (start[1] == 0) ? count[1] + 1 : count[1] + 2;
  count[0] = (start[0] == 0) ? count[0] + 1 : count[0] + 2;

  /* -- Phi function routine and FSM -- */
  Phi *pf = (Phi *) malloc(sizeof(Phi));
  if (!pf) {
    printf("Error allocating memory for the Phi function.\n");
    exit(1);
  }

  makePhi3D(count, dlt, pf);

  /* Check for available device memory */
  int    totalCount = count[0] * count[1] * count[2];
  size_t reqMem     = sizeof(Phi) + (sizeof(double) + sizeof(int)) * totalCount;
  size_t devMem     = acc_get_free_memory();

  if (devMem < reqMem) {
    printf("[%d/%d] Not enough device memory ...\n", mpiRank, mpiSize);
    printf("[%d/%d] Required Memory(~): %zu\n Available Memory: %zu\n", mpiRank, mpiSize,
           reqMem, devMem);
    exit(1);
  }

  int dstVarId; /* Distance Variable ID */
  int flgVarId; /* Flags Variable ID */
  ncGetVarid(ncid, "Distance", &dstVarId);
  ncGetVarid(ncid, "Flags", &flgVarId);

  /* Grab slab of data from the NetCDF file */
  LOG_DEBUG("[%d/%d] Get Double ...............", mpiRank, mpiSize);
  ncGetDouble(ncid, dstVarId, start, count, pf->distance);

  double startTime, finishTime;
  MPI_Barrier(MPI_COMM_WORLD);

  GET_TIME(startTime);

  runFSMAndPropagate(pf);

  MPI_Barrier(MPI_COMM_WORLD);

  GET_TIME(finishTime);
  if (mpiRank == 0) {
    double elapTime = finishTime - startTime;
    printf("Parallel MPI/OpenACC Block FSM time: %f s.\n", elapTime);
  }

  LOG_DEBUG("[%d/%d] Put Double ...............", mpiRank, mpiSize);
  ncPutDouble(ncid, dstVarId, start, count, pf->distance);

  destroyPhi3D(pf);
  free(pf);
}

static void receive(Phi *pf, int rx, int ry, int rz, int pCount, int tag, int xs, int xe, int ys,
                    int ye, int zs, int ze)
{
  int     recRank = linearFrom3D(blockDim[0], blockDim[1], rx, ry, rz);
  double *plane   = (double *) calloc(pCount, sizeof(double));

  MPI_Request req;
  MPI_Status  status;

  MPI_Irecv(plane, pCount, MPI_DOUBLE, recRank, tag, MPI_COMM_WORLD, &req);

  MPI_Wait(&req, &status);

  // int recCount;
  // MPI_Get_count(&status, MPI_DOUBLE, &recCount);

  /* Update */
  int itr = 0;
  for (int i = zs; i < ze; i++) {
    for (int j = ys; j < ye; j++) {
      for (int k = xs; k < xe; k++) {
        int idx           = k + count[2] * j + count[2] * count[1] * i;
        pf->distance[idx] = plane[itr++];
      }
    }
  }
  free(plane);
}

static void send(Phi *pf, int sx, int sy, int sz, int pCount, int tag, int xs, int xe, int ys,
                 int ye, int zs, int ze)
{
  int     sendRank = linearFrom3D(blockDim[0], blockDim[1], sx, sy, sz);
  double *plane    = (double *) calloc(pCount, sizeof(double));

  /* Fetch */
  int itr = 0;
  for (int i = zs; i < ze; i++) {
    for (int j = ys; j < ye; j++) {
      for (int k = xs; k < xe; k++) {
        int idx      = k + count[2] * j + count[2] * count[1] * i;
        plane[itr++] = pf->distance[idx];
      }
    }
  }

  MPI_Request req;
  MPI_Status  status;

  int ierr = MPI_Isend(plane, pCount, MPI_DOUBLE, sendRank, tag, MPI_COMM_WORLD, &req);
  MPI_Wait(&req, &status);

  if (ierr != MPI_SUCCESS) {
    int  errClass, resultLen;
    char err_buffer[MPI_MAX_ERROR_STRING];
    MPI_Error_class(ierr, &errClass);
    MPI_Error_string(ierr, err_buffer, &resultLen);
    printf("[%d] Error: %s\n", mpiRank, err_buffer);
  }

  // int sendCount;
  // MPI_Get_count(&status, MPI_DOUBLE, &sendCount);

  free(plane);
}

static void runFSMAndPropagate(Phi *pf)
{
  int nx = blockDim[0];
  int ny = blockDim[1];
  int nz = blockDim[2];

  int xtag = 2;
  int ytag = 1;
  int ztag = 0;

  /* Each row defines the start coordinate for each sweep iteration */
  int sCoord[8][3] = {{1, 1, 1},    {nx, 1, nz}, {1, 1, nz},  {nx, 1, 1},
                      {nx, ny, nz}, {1, ny, 1},  {nx, ny, 1}, {1, ny, nz}};

  int yzCount = count[1] * count[0];
  int xzCount = count[2] * count[0];
  int xyCount = count[2] * count[1];

  int rcount = 0;
  for (int sn = 1; sn <= 8; sn++) {
    if (x == sCoord[sn - 1][0] && y == sCoord[sn - 1][1] && z == sCoord[sn - 1][2]) {
      runFSM(pf, sn);
    } else {
      /* Receive */
      int rx = (sn == 2 || sn == 4 || sn == 5 || sn == 7) ? (x + 1) : (x - 1);
      int ry = (sn == 5 || sn == 6 || sn == 7 || sn == 8) ? (y + 1) : (y - 1);
      int rz = (sn == 2 || sn == 3 || sn == 5 || sn == 8) ? (z + 1) : (z - 1);

      if (rx > 0 && rx <= blockDim[0]) {
        if (rx == (x - 1)) {
          receive(pf, rx, y, z, yzCount, xtag, 0, 1, 0, count[1], 0, count[0]);
        } else {
          receive(pf, rx, y, z, yzCount, xtag, count[2] - 1, count[2], 0, count[1], 0,
                  count[0]);
        }
      }

      if (ry > 0 && ry <= blockDim[1]) {
        if (ry == (y - 1)) {
          receive(pf, x, ry, z, xzCount, ytag, 0, count[2], 0, 1, 0, count[0]);
        } else {
          receive(pf, x, ry, z, xzCount, ytag, 0, count[2], count[1] - 1, count[1], 0,
                  count[0]);
        }
      }

      if (rz > 0 && rz <= blockDim[2]) {
        if (rz == (z - 1)) {
          receive(pf, x, y, rz, xyCount, ztag, 0, count[2], 0, count[1], 0, 1);
        } else {
          receive(pf, x, y, rz, xyCount, ztag, 0, count[2], 0, count[1], count[0] - 1,
                  count[0]);
        }
      }

      runFSM(pf, sn);
    }

    /* Send */
    int sx = (sn == 2 || sn == 4 || sn == 5 || sn == 7) ? (x - 1) : (x + 1);
    int sy = (sn == 5 || sn == 6 || sn == 7 || sn == 8) ? (y - 1) : (y + 1);
    int sz = (sn == 2 || sn == 3 || sn == 5 || sn == 8) ? (z - 1) : (z + 1);

    if (sx > 0 && sx <= blockDim[0]) {
      if (sx == (x + 1)) {
        send(pf, sx, y, z, yzCount, xtag, count[2] - 2, count[2] - 1, 0, count[1], 0,
             count[0]);
      } else {
        send(pf, sx, y, z, yzCount, xtag, 1, 2, 0, count[1], 0, count[0]);
      }
    }

    if (sy > 0 && sy <= blockDim[1]) {
      if (sy == (y + 1)) {
        send(pf, x, sy, z, xzCount, ytag, 0, count[2], count[1] - 2, count[1] - 1, 0,
             count[0]);
      } else {
        send(pf, x, sy, z, xzCount, ytag, 0, count[2], 1, 2, 0, count[0]);
      }
    }

    if (sz > 0 && sz <= blockDim[2]) {
      if (sz == (z + 1)) {
        send(pf, x, y, sz, xyCount, ztag, 0, count[2], 0, count[1], count[0] - 2,
             count[0] - 1);
      } else {
        send(pf, x, y, sz, xyCount, ztag, 0, count[2], 0, count[1], 1, 2);
      }
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }
}
