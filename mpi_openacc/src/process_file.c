/**
 * @file process_file.c
 * @brief      Source file for pre/post processing of the NetCDF file for FSM.
 *
 * @author     Shrestha, Anup
 * @date       12 JUL 2016
 */

#include "process_file.h"
#include "macrologger.h"
#include "ncwrap.h"
#include "phi3D.h"
#include "utilities.h"

#include <math.h>
#include <mpi.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#define min(X, Y) ((X) < (Y) ? (X) : (Y))
#define max(X, Y) ((X) > (Y) ? (X) : (Y))

/* MPI routine variables */
static int mpiSize;
static int mpiRank;
/*-----------------------------------------*/

/* Dimension of decomposed blocks */
static int BOX_X;
static int BOX_Y;
static int BOX_Z;

static size_t *blockDim;
/*-----------------------------------------*/

static inline void copyAndUpdateProcess(int ncidIn, int fVidIn, int dVidIn, int ncidOut,
                                        int fVidOut, int dVidOut, size_t *dim, size_t *blockDim);
static inline void updateDistance(double *dst, int totalNodes);
static inline void fillIntArray(int *arr, int n, const int val);
static inline void fillDblArray(double *arr, int n, const double val);
static void setBorderFace(size_t *start, size_t *count, int ncid, int fId, int dId);
static inline void setBorderDistance(int ncid, int flgVarId, int dstVarId, size_t *dim);
static inline void setInsideDistance(int *flg, double *dst, int totalNodes);
static inline void adjustBoundary(int ncid, int dstVarId, size_t *dim);
static void fixBorder(size_t *start_r, size_t *start_w, size_t *count, int ncid, int dstVarId);
void linearTo3D(int i, int width, int height, int *x, int *y, int *z);

/**
 * @brief      Pre-processing for generating the input for FSM.
 *
 * @param[in]  ncidIn   NetCDF ID of input file.
 * @param[in]  ncidOut  NetCDF ID of output file.
 * @param      blkDim   Length of block dimensions.
 */
void preProcessFile(int ncidIn, int ncidOut, size_t *blkDim)
{
  MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiRank);

  blockDim = blkDim;

  int dstVarIdIn; /* Distance Variable ID for input file */
  int flgVarIdIn; /* Flags Variable ID for input file */
  ncGetVarid(ncidIn, "Distance", &dstVarIdIn);
  ncGetVarid(ncidIn, "Flags", &flgVarIdIn);

  /* Retrieve the dimension lengths and delta values */
  size_t dimIn[NDIMS];
  double dlt[NDIMS];
  ncGetDim(ncidIn, dimIn, dlt);

  size_t dimOut[NDIMS] = {dimIn[0] + 2, dimIn[1] + 2, dimIn[2] + 2};

  int dimids[NDIMS];
  ncPutDim(ncidOut, dimids, dimOut, dlt);

  int dstVarIdOut;
  int flgVarIdOut;
  ncGetVarid(ncidOut, "Distance", &dstVarIdOut);
  ncGetVarid(ncidOut, "Flags", &flgVarIdOut);

  /* INDEPENDENT parallel access
  ncParAccess(ncidOut, dstVarIdIn);
  ncParAccess(ncidOut, flgVarIdIn);
  ncParAccess(ncidOut, dstVarIdOut);
  ncParAccess(ncidOut, flgVarIdOut);
  */

  BOX_X = (blockDim[0] == 1) ? dimOut[0] : (dimOut[0] / blockDim[0]) + 1;
  BOX_Y = (blockDim[1] == 1) ? dimOut[1] : (dimOut[1] / blockDim[1]) + 1;
  BOX_Z = (blockDim[2] == 1) ? dimOut[2] : (dimOut[2] / blockDim[2]) + 1;

  if (mpiRank == 0) printf("BOX_X: %d, BOX_Y: %d, BOX_Z: %d\n", BOX_X, BOX_Y, BOX_Z);

  LOG_DEBUG("[%d/%d] Copy and Update ..", mpiRank, mpiSize);
  copyAndUpdateProcess(ncidIn, flgVarIdIn, dstVarIdIn, ncidOut, flgVarIdOut, dstVarIdOut, dimIn,
                       blockDim);

  LOG_DEBUG("[%d/%d] Set Border .......", mpiRank, mpiSize);
  setBorderDistance(ncidOut, flgVarIdOut, dstVarIdOut, dimOut);
}

/**
 * @brief      Copies the values from the input file and updates the border
 *             values.
 *
 * @param[in]  ncidIn    NetCDF ID of input file.
 * @param[in]  fVidIn    Flags variable id of input file.
 * @param[in]  dVidIn    Distance variable id of input file.
 * @param[in]  ncidOut   NetCDF ID of output file.
 * @param[in]  fVidOut   Flags variable id of input file.
 * @param[in]  dVidOut   Distance variable id of output file.
 * @param      dim       Dimensions of the whole domain.
 * @param      blockDim  Dimensions of the decomposed domain.
 */
static inline void copyAndUpdateProcess(int ncidIn, int fVidIn, int dVidIn, int ncidOut,
                                        int fVidOut, int dVidOut, size_t *dim, size_t *blockDim)
{
  int zcount = ceil((double) (dim[2]) / mpiSize);

  size_t startIdxIn[NDIMS] = {mpiRank * zcount, 0, 0};
  size_t count[NDIMS]      = {0, dim[1], dim[0]};

  count[0] = min(dim[2], startIdxIn[0] + zcount) - startIdxIn[0];

  size_t startIdxOut[NDIMS] = {startIdxIn[0] + 1, startIdxIn[1] + 1, startIdxIn[2] + 1};

  int     totalNodes = count[2] * count[1] * count[0];
  int *   flg        = (int *) malloc(totalNodes * sizeof(int));
  double *dst        = (double *) malloc(totalNodes * sizeof(double));

  LOG_DEBUG("IN: [%d/%d] (%zu-%zu, %zu-%zu, %zu-%zu)", mpiRank, mpiSize, startIdxIn[2],
            startIdxIn[2] + count[2] - 1, startIdxIn[1], startIdxIn[1] + count[1] - 1,
            startIdxIn[0], startIdxIn[0] + count[0] - 1);
  LOG_DEBUG("OU: [%d/%d] (%zu-%zu, %zu-%zu, %zu-%zu)", mpiRank, mpiSize, startIdxOut[2],
            startIdxOut[2] + count[2] - 1, startIdxOut[1], startIdxOut[1] + count[1] - 1,
            startIdxOut[0], startIdxOut[0] + count[0] - 1);

  ncGetInt(ncidIn, fVidIn, startIdxIn, count, flg);

  ncGetDouble(ncidIn, dVidIn, startIdxIn, count, dst);

  ncPutInt(ncidOut, fVidOut, startIdxOut, count, flg);

  updateDistance(dst, totalNodes);

  ncPutDouble(ncidOut, dVidOut, startIdxOut, count, dst);

  free(flg);
  free(dst);
}

/**
 * @brief         Update unknown distance values to a default high number.
 *
 * @param[in,out] dst         The distance value array.
 * @param[in]     totalNodes  Total number of values in the distance array.
 */
static inline void updateDistance(double *dst, int totalNodes)
{
  double *dval = dst;
  for (int t = 0; t < totalNodes; t++) {
    if (!(*dval > 0.0 || *dval < 0.0)) *dval = DEFAULT_INTERIOR_DISTANCE;
    dval++;
  }
}

/**
 * @brief      Fills an array of int.
 *
 * @param      arr   Pointer to the int array.
 * @param[in]  n     Length of the array.
 * @param[in]  val   The value to fill in.
 */
static inline void fillIntArray(int *arr, int n, const int val)
{
  for (int i = 0; i < n; i++) {
    *(arr++) = val;
  }
}

/**
 * @brief      Fills an array of double.
 *
 * @param      arr   Pointer to the double array.
 * @param[in]  n     Length of the array.
 * @param[in]  val   The value to fill in.
 */
static inline void fillDblArray(double *arr, int n, const double val)
{
  for (int i = 0; i < n; i++) {
    *(arr++) = val;
  }
}

/**
 * @brief      Sets the border face.
 *
 * @param      start  Start index vector.
 * @param      count  Count vector.
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  fId    Flags variable id.
 * @param[in]  dId    Distance variable id.
 */
static void setBorderFace(size_t *start, size_t *count, int ncid, int fId, int dId)
{
  int     n  = count[2] * count[1] * count[0];
  int *   bl = (int *) malloc(n * sizeof(int));
  double *bd = (double *) malloc(n * sizeof(double));
  fillIntArray(bl, n, DEFAULT_BORDER_LOCATION);
  fillDblArray(bd, n, DEFAULT_BORDER_DISTANCE);
  ncPutInt(ncid, fId, start, count, bl);
  ncPutDouble(ncid, dId, start, count, bd);
  free(bl);
  free(bd);
}

/**
 * @brief      Sets the border distance.
 *
 * @param[in]  ncid      NetCDF ID.
 * @param[in]  flgVarId  Flags variable id.
 * @param[in]  dstVarId  Distance variable id.
 * @param      dim       Dimensions
 */
static inline void setBorderDistance(int ncid, int flgVarId, int dstVarId, size_t *dim)
{
  size_t start[NDIMS];
  size_t count[NDIMS];

  if (mpiSize == 1 || (mpiSize >= 2 && mpiRank == 0)) {
    // left face
    start[2] = 0;
    start[1] = 0;
    start[0] = 0;
    count[2] = 1;
    count[1] = dim[1];
    count[0] = dim[2];
    LOG_DEBUG("[%d/%d] Set left border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start[2],
              start[2] + count[2] - 1, start[1], start[1] + count[1] - 1, start[0],
              start[0] + count[0] - 1);
    setBorderFace(start, count, ncid, flgVarId, dstVarId);
  }

  if (mpiSize == 1 || (mpiSize >= 2 && mpiRank == 1)) {
    // right face
    start[2] = dim[0] - 1;
    start[1] = 0;
    start[0] = 0;
    count[2] = 1;
    count[1] = dim[1];
    count[0] = dim[2];
    LOG_DEBUG("[%d/%d] Set right border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start[2],
              start[2] + count[2] - 1, start[1], start[1] + count[1] - 1, start[0],
              start[0] + count[0] - 1);
    setBorderFace(start, count, ncid, flgVarId, dstVarId);
  }

  if (mpiSize == 1 || (mpiSize == 2 && mpiRank == 0) || (mpiSize > 2 && mpiRank > 1)) {
    // bottom face
    start[2] = 0;
    start[1] = 0;
    start[0] = 0;
    count[2] = dim[0];
    count[1] = 1;
    count[0] = dim[2];
    LOG_DEBUG("[%d/%d] Set bottom border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start[2],
              start[2] + count[2] - 1, start[1], start[1] + count[1] - 1, start[0],
              start[0] + count[0] - 1);
    setBorderFace(start, count, ncid, flgVarId, dstVarId);

    // top face
    start[2] = 0;
    start[1] = dim[1] - 1;
    start[0] = 0;
    count[2] = dim[0];
    count[1] = 1;
    count[0] = dim[2];
    LOG_DEBUG("[%d/%d] Set top border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start[2],
              start[2] + count[2] - 1, start[1], start[1] + count[1] - 1, start[0],
              start[0] + count[0] - 1);
    setBorderFace(start, count, ncid, flgVarId, dstVarId);
  }

  if (mpiSize == 1 || (mpiSize == 2 && mpiRank == 1) || (mpiSize > 2 && mpiRank > 1)) {
    // front face
    start[2] = 0;
    start[1] = 0;
    start[0] = 0;
    count[2] = dim[0];
    count[1] = dim[1];
    count[0] = 1;
    LOG_DEBUG("[%d/%d] Set front border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start[2],
              start[2] + count[2] - 1, start[1], start[1] + count[1] - 1, start[0],
              start[0] + count[0] - 1);
    setBorderFace(start, count, ncid, flgVarId, dstVarId);

    // back face
    start[2] = 0;
    start[1] = 0;
    start[0] = dim[2] - 1;
    count[2] = dim[0];
    count[1] = dim[1];
    count[0] = 1;
    LOG_DEBUG("[%d/%d] Set back border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start[2],
              start[2] + count[2] - 1, start[1], start[1] + count[1] - 1, start[0],
              start[0] + count[0] - 1);
    setBorderFace(start, count, ncid, flgVarId, dstVarId);
  }
}

/**
 * @brief      Post-processing of the output file.
 *
 * @param[in]  ncid  NetCDF ID.
 */
void postProcessFile(int ncid)
{
  size_t dim[NDIMS];
  double dlt[NDIMS];
  ncGetDim(ncid, dim, dlt);

  int dstVarId; /* Distance Variable ID */
  int flgVarId; /* Flags Variable ID */
  ncGetVarid(ncid, "Distance", &dstVarId);
  ncGetVarid(ncid, "Flags", &flgVarId);

  /* INDEPENDENT parallel access */
  ncParAccess(ncid, dstVarId);
  ncParAccess(ncid, flgVarId);

  int zcount = ceil((double) (dim[2]) / mpiSize);

  size_t start[NDIMS] = {mpiRank * zcount, 0, 0};
  size_t count[NDIMS] = {0, dim[1], dim[0]};

  count[0] = min(dim[2], start[0] + zcount) - start[0];

  int     totalNodes = count[0] * count[1] * count[2];
  int *   flg        = (int *) malloc(totalNodes * sizeof(int));
  double *dst        = (double *) malloc(totalNodes * sizeof(double));

  ncGetInt(ncid, flgVarId, start, count, flg);

  ncGetDouble(ncid, dstVarId, start, count, dst);

  LOG_DEBUG("[%d/%d] Set distance negative ...", mpiRank, mpiSize);
  LOG_DEBUG("[%d/%d] (%zu-%zu, %zu-%zu, %zu-%zu)", mpiRank, mpiSize, start[2],
            start[2] + count[2] - 1, start[1], start[1] + count[1] - 1, start[0],
            start[0] + count[0] - 1);

  setInsideDistance(flg, dst, totalNodes);

  ncPutDouble(ncid, dstVarId, start, count, dst);

  free(flg);
  free(dst);

  MPI_Barrier(MPI_COMM_WORLD);

  adjustBoundary(ncid, dstVarId, dim);
}

/**
 * @brief         Sets the distance values for the positions that are inside an
 *                object is to -1
 *
 * @param[in]     flg         The Flags array
 * @param[in,out] dst         The Distance array
 * @param[in]     totalNodes  Total number of nodes
 */
static inline void setInsideDistance(int *flg, double *dst, int totalNodes)
{
  int *   fval = flg;
  double *dval = dst;
  for (int t = 0; t < totalNodes; t++) {
    if (*fval == 1) {
      *dval = -1;
    }
    fval++;
    dval++;
  }
}

/**
 * @brief      Copies the values of the inner border to outer border.
 *
 * @param      start_r   Start index vector to copy from.
 * @param      start_w   Start index vector to write to.
 * @param      count     Count vector.
 * @param[in]  ncid      NetCDF ID.
 * @param[in]  dstVarId  Distance variable id.
 */
static void fixBorder(size_t *start_r, size_t *start_w, size_t *count, int ncid, int dstVarId)
{
  int     n = count[2] * count[1] * count[0];
  double *d = (double *) malloc(n * sizeof(double));
  ncGetDouble(ncid, dstVarId, start_r, count, d);
  ncPutDouble(ncid, dstVarId, start_w, count, d);
  free(d);
}

/**
 * @brief      Adjusts the boundary values.
 *
 * @param[in]  ncid      NetCDF ID.
 * @param[in]  dstVarId  Distance variable id.
 * @param      dim       Dimensions.
 */
static inline void adjustBoundary(int ncid, int dstVarId, size_t *dim)
{
  size_t start_r[NDIMS];
  size_t start_w[NDIMS];
  size_t count[NDIMS];

  // left face
  start_r[2] = 1;
  start_r[1] = 0;
  start_r[0] = 0;
  start_w[2] = 0;
  start_w[1] = 0;
  start_w[0] = 0;
  count[2]   = 1;
  count[1]   = dim[1];
  count[0]   = dim[2];
  LOG_DEBUG("[%d/%d] Adjust left border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start_w[2],
            start_w[2] + count[2] - 1, start_w[1], start_w[1] + count[1] - 1, start_w[0],
            start_w[0] + count[0] - 1);
  fixBorder(start_r, start_w, count, ncid, dstVarId);

  // right face
  start_r[2] = dim[0] - 2;
  start_r[1] = 0;
  start_r[0] = 0;
  start_w[2] = dim[0] - 1;
  start_w[1] = 0;
  start_w[0] = 0;
  count[2]   = 1;
  count[1]   = dim[1];
  count[0]   = dim[2];
  LOG_DEBUG("[%d/%d] Adjust right border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start_w[2],
            start_w[2] + count[2] - 1, start_w[1], start_w[1] + count[1] - 1, start_w[0],
            start_w[0] + count[0] - 1);
  fixBorder(start_r, start_w, count, ncid, dstVarId);

  // bottom face
  start_r[2] = 0;
  start_r[1] = 1;
  start_r[0] = 0;
  start_w[2] = 0;
  start_w[1] = 0;
  start_w[0] = 0;
  count[2]   = dim[0];
  count[1]   = 1;
  count[0]   = dim[2];
  LOG_DEBUG("[%d/%d] Adjust bottom border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize,
            start_w[2], start_w[2] + count[2] - 1, start_w[1], start_w[1] + count[1] - 1,
            start_w[0], start_w[0] + count[0] - 1);
  fixBorder(start_r, start_w, count, ncid, dstVarId);

  // top face
  start_r[2] = 0;
  start_r[1] = dim[1] - 2;
  start_r[0] = 0;
  start_w[2] = 0;
  start_w[1] = dim[1] - 1;
  start_w[0] = 0;
  count[2]   = dim[0];
  count[1]   = 1;
  count[0]   = dim[2];
  LOG_DEBUG("[%d/%d] Adjust top border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start_w[2],
            start_w[2] + count[2] - 1, start_w[1], start_w[1] + count[1] - 1, start_w[0],
            start_w[0] + count[0] - 1);
  fixBorder(start_r, start_w, count, ncid, dstVarId);

  // front face
  start_r[2] = 0;
  start_r[1] = 0;
  start_r[0] = 1;
  start_w[2] = 0;
  start_w[1] = 0;
  start_w[0] = 0;
  count[2]   = dim[0];
  count[1]   = dim[1];
  count[0]   = 1;
  LOG_DEBUG("[%d/%d] Adjust front border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start_w[2],
            start_w[2] + count[2] - 1, start_w[1], start_w[1] + count[1] - 1, start_w[0],
            start_w[0] + count[0] - 1);
  fixBorder(start_r, start_w, count, ncid, dstVarId);

  // back face
  start_r[2] = 0;
  start_r[1] = 0;
  start_r[0] = dim[2] - 2;
  start_w[2] = 0;
  start_w[1] = 0;
  start_w[0] = dim[2] - 1;
  count[2]   = dim[0];
  count[1]   = dim[1];
  count[0]   = 1;
  LOG_DEBUG("[%d/%d] Adjust back border (%zu-%zu,%zu-%zu,%zu-%zu)", mpiRank, mpiSize, start_w[2],
            start_w[2] + count[2] - 1, start_w[1], start_w[1] + count[1] - 1, start_w[0],
            start_w[0] + count[0] - 1);
  fixBorder(start_r, start_w, count, ncid, dstVarId);

  MPI_Barrier(MPI_COMM_WORLD);
}

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
void linearTo3D(int i, int width, int height, int *x, int *y, int *z)
{
  *x = (i % width) + 1;
  *y = ((i / width) % height) + 1;
  *z = (i / (width * height)) + 1;
}
