/**
 * @file ncwrap.c
 * @brief      Source file for the wrapper functions to NetCDF-4 API calls.
 *
 *             A custom wrapper for some of the NetCDF-4 API functions,
 *             especially designed to work with the parallel-io distance field
 *             calculation code. Some of the arguments are hard-coded to shorten
 *             the parameter list for the wrapper functions. The wrapper
 *             functions also check if the NetCDF-4 API call was successful. If
 *             the API call failed then an error message is printed and the
 *             program is quit.
 *
 * @author     Shrestha, Anup
 * @date       18 JUL 2016
 *
 * @note       The wrapper functions are especially designed to work with the
 *             parallel-io distance field calculation code. Do not use this
 *             wrapper for general use.
 *
 * @see
 * https://www.unidata.ucar.edu/software/NetCDF/NetCDF-4/newdocs/NetCDF-c.pdf
 */
#include "ncwrap.h"
#include "phi3D.h"

#include <mpi.h>
#include <netcdf.h>
#include <netcdf_par.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief      Open a parallel NetCDF-4 file.
 *
 *             Wrapper for nc_open_par() which opens a NetCDF-4 file for
 *             parallel access using MPI/IO.
 *
 * @param[in]  filename  The file name of the NetCDF dataset to be opened.
 * @param[out] ncid      Pointer to location where returned NetCDF id is to be
 */
void ncOpen(char *filename, int *ncid)
{
  int res; /* NetCDF return status code */

  if ((res = nc_open_par(filename, NC_NETCDF4 | NC_MPIIO | NC_WRITE, MPI_COMM_WORLD,
                         MPI_INFO_NULL, ncid)))
    ERR(res);
}

/**
 * @brief      Create a parallel NetCDF-4 file
 *
 *             Wrapper for nc_create_par() which creates a NetCDF-4 file on a
 *             MPI/IO parallel file system.
 *
 * @param[in]  filename  The file name of the new NetCDF dataset.
 * @param[out] ncid      Pointer to location where returned NetCDF id is to be
 *                       stored.
 */
void ncCreate(char *filename, int *ncid)
{
  int res; /* NetCDF return status code */

  if ((res = nc_create_par(filename, NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, ncid)))
    ERR(res);
}

/**
 * @brief Get the id of a NetCDF variable, given its name.
 *
 * Wrapper for nc_inq_varid().
 *
 * @param[in]  ncid    NetCDF ID.
 * @param[in]  varName The name of the variable to look for.
 * @param[out] varid   Pointer to location where returned NetCDF variable id is
 * to be stored.
 */
void ncGetVarid(int ncid, char *varName, int *varid)
{
  int res; /* NetCDF return status code */

  if ((res = nc_inq_varid(ncid, varName, varid))) ERR(res);
}

/**
 * @brief      Get dimensions and delta values.
 *
 *             Stores the length of the dimensions in the dim variable and the
 *             values of the delta attribute for each dimension in the dlt
 *             variable.
 *
 * @param[in]  ncid  NetCDF ID.
 * @param[out] dim   Pointer to location where length of the dimensions are to
 *                   be stored.
 * @param[out] dlt   Pointer to location where values of the delta attribute are
 *                   to be stored.
 */
void ncGetDim(int ncid, size_t *dim, double *dlt)
{
  int res; /* NetCDF return status code */
  int ndims;
  int dimids[NDIMS];

  if ((res = nc_inq_dimids(ncid, &ndims, dimids, 0))) ERR(res);

  /* Dimensions */
  if ((res = nc_inq_dimlen(ncid, dimids[0], &dim[0]))) ERR(res);

  if ((res = nc_inq_dimlen(ncid, dimids[1], &dim[1]))) ERR(res);

  if ((res = nc_inq_dimlen(ncid, dimids[2], &dim[2]))) ERR(res);

  /* Delta */
  if ((res = nc_get_att_double(ncid, dimids[0], "delta", &dlt[0]))) ERR(res);

  if ((res = nc_get_att_double(ncid, dimids[1], "delta", &dlt[1]))) ERR(res);

  if ((res = nc_get_att_double(ncid, dimids[2], "delta", &dlt[2]))) ERR(res);
}

/**
 * @brief      Define and populate dimension and delta values.
 *
 *             Writes the length and the data values for the dimension and the
 *             delta values for each dimension to the NetCDF file.
 *
 * @param[in]  ncid    NetCDF ID.
 * @param[out] dimids  Pointer to location where returned NetCDF dimension ids
 *                     are to be stored.
 * @param[out] dim     Length of the dimensions.
 * @param[out] dlt     Delta values for each dimension.
 */
void ncPutDim(int ncid, int *dimids, size_t *dim, double *dlt)
{
  int res; /* NetCDF return status code */
  int dimvarids[NDIMS];
  int flgVarId, dstVarId;

  /* Dimensions */
  if ((res = nc_def_dim(ncid, "x", dim[0], &dimids[2]))) ERR(res);

  if ((res = nc_def_dim(ncid, "y", dim[1], &dimids[1]))) ERR(res);

  if ((res = nc_def_dim(ncid, "z", dim[2], &dimids[0]))) ERR(res);

  if ((res = nc_def_var(ncid, "x", NC_FLOAT, 1, &dimids[2], &dimvarids[2]))) ERR(res);

  if ((res = nc_def_var(ncid, "y", NC_FLOAT, 1, &dimids[1], &dimvarids[1]))) ERR(res);

  if ((res = nc_def_var(ncid, "z", NC_FLOAT, 1, &dimids[0], &dimvarids[0]))) ERR(res);

  if ((res = nc_def_var(ncid, "Flags", NC_INT, NDIMS, dimids, &flgVarId))) ERR(res);
  if ((res = nc_def_var(ncid, "Distance", NC_DOUBLE, NDIMS, dimids, &dstVarId))) ERR(res);

  // int dbl = DEFAULT_BORDER_LOCATION;
  // int dbd = DEFAULT_BORDER_DISTANCE;
  // size_t chunkSize[NDIMS] = {129, 97, 129};

  // if((res = nc_def_var_fill(ncid, flgVarId, 0, &dbl))) ERR(res);
  // if((res = nc_def_var_fill(ncid, dstVarId, 0, &dbd))) ERR(res);

  // if((res = nc_def_var_chunking(ncid, flgVarId, NC_CHUNKED, chunkSize))) ERR(res);
  // if((res = nc_def_var_chunking(ncid, dstVarId, NC_CHUNKED, chunkSize))) ERR(res);

  /* End define mode */
  if ((res = nc_enddef(ncid))) ERR(res);

  /* Data for x, y, z variables */
  float lons[dim[0]], lats[dim[1]], deps[dim[2]];
  for (int k  = 0; k < dim[0]; k++)
    lons[k] = dlt[0] * k;
  for (int j  = 0; j < dim[1]; j++)
    lats[j] = dlt[1] * j;
  for (int i  = 0; i < dim[2]; i++)
    deps[i] = dlt[2] * i;

  /* Put the data values for the dimensions */
  if ((res = nc_put_var_float(ncid, dimvarids[2], &lons[0]))) ERR(res);

  if ((res = nc_put_var_float(ncid, dimvarids[1], &lats[0]))) ERR(res);

  if ((res = nc_put_var_float(ncid, dimvarids[0], &deps[0]))) ERR(res);

  if ((res = nc_redef(ncid))) ERR(res);

  /* Put the delta values */
  if ((res = nc_put_att_double(ncid, dimvarids[2], "delta", NC_DOUBLE, 1, &dlt[0]))) ERR(res);

  if ((res = nc_put_att_double(ncid, dimvarids[1], "delta", NC_DOUBLE, 1, &dlt[1]))) ERR(res);

  if ((res = nc_put_att_double(ncid, dimvarids[0], "delta", NC_DOUBLE, 1, &dlt[2]))) ERR(res);

  if ((res = nc_enddef(ncid))) ERR(res);
}

/**
 * @brief      NetCDF parallel independent access.
 *
 *             Wrapper for nc_var_par_access().
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 */
void ncParAccess(int ncid, int varid)
{
  int res; /* NetCDF return status code */

  if ((res = nc_var_par_access(ncid, varid, NC_INDEPENDENT))) ERR(res);
}

/**
 * @brief      Define variable.
 *
 *             Wrapper for nc_def_var() which adds a new variable to an open
 *             NetCDF dataset in define mode.
 *
 * @param[in]  ncid     NetCDF ID.
 * @param[in]  varName  Variable name.
 * @param[in]  xtype    Predefined NetCDF external data type.
 * @param[in]  dimids   Dimension IDs.
 * @param[out] varid    Pointer to location for the returned variable ID to be
 * stored.
 */
void ncDefVar(int ncid, char *varName, nc_type xtype, const int dimids[], int *varid)
{
  int res; /* NetCDF return status code */
  if ((res = nc_def_var(ncid, varName, xtype, NDIMS, dimids, varid))) ERR(res);
}

/**
 * @brief      Read a single integer value.
 *
 *             Wrapper for nc_get_var1_int() which reads a single integer value
 *             from an open NetCDF dataset that is in data mode.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  index  Index of the value to get.
 * @param[out] ip     Pointer to location for the returned value to be stored.
 */
void ncGet1Int(int ncid, int varid, const size_t index[], int *ip)
{
  int res; /* NetCDF return status code */

  if ((res = nc_get_var1_int(ncid, varid, index, ip))) ERR(res);
}

/**
 * @brief      Read an array of integer values.
 *
 *             Wrapper for nc_get_vara_int() which reads an array of integer
 *             values from an open NetCDF dataset that is in data mode.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  start  Index position to start at.
 * @param[in]  count  Vector of edge lengths.
 * @param[out] ip     Pointer to location for the returned values to be stored.
 */
void ncGetInt(int ncid, int varid, const size_t start[], const size_t count[], int *ip)
{
  int res; /* NetCDF return status code */

  if ((res = nc_get_vara_int(ncid, varid, start, count, ip))) ERR(res);
}

/**
 * @brief      Write a single integer value.
 *
 *             Wrapper for nc_put_var1_int() which writes a single integer value
 *             to an open NetCDF dataset that is in data mode.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  index  Index where the value to write.
 * @param[out] ip     Pointer to location of the value to write.
 */
void ncPut1Int(int ncid, int varid, const size_t index[], int *ip)
{
  int res; /* NetCDF return status code */

  if ((res = nc_put_var1_int(ncid, varid, index, ip))) ERR(res);
}

/**
 * @brief      Write an array of integer values.
 *
 *             Wrapper for nc_put_vara_int() which writes an array of integer
 *             values from an open NetCDF dataset that is in data mode.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  start  Index position to start at.
 * @param[in]  count  Vector of edge lengths.
 * @param[out] ip     Pointer to location of the values to write.
 */
void ncPutInt(int ncid, int varid, const size_t start[], const size_t count[], int *ip)
{
  int res; /* NetCDF return status code */

  if ((res = nc_put_vara_int(ncid, varid, start, count, ip))) ERR(res);
}

/**
 * @brief      Read a single double value.
 *
 *             Wrapper for nc_get_var1_double() which reads a single double
 *             value from an open NetCDF dataset that is in data mode.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  index  Index of the value to get.
 * @param[out] ip     Pointer to location for the returned value to be stored.
 */
void ncGet1Double(int ncid, int varid, const size_t index[], double *ip)
{
  int res; /* NetCDF return status code */

  if ((res = nc_get_var1_double(ncid, varid, index, ip))) ERR(res);
}

/**
 * @brief      Read an array of double values.
 *
 *             Wrapper for nc_get_vara_double() which reads an array of double
 *             values from an open NetCDF dataset that is in data mode.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  start  Index position to start at.
 * @param[in]  count  Vector of edge lengths.
 * @param[out] ip     Pointer to location for the returned values to be stored.
 */
void ncGetDouble(int ncid, int varid, const size_t start[], const size_t count[], double *ip)
{
  int res; /* NetCDF return status code */

  if ((res = nc_get_vara_double(ncid, varid, start, count, ip))) ERR(res);
}

/**
 * @brief      Write a single double value.
 *
 *             Wrapper for nc_put_var1_double() which writes a single double
 *             value to an open NetCDF dataset that is in data mode.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  index  Index where the value to write.
 * @param[out] ip     Pointer to location of the value to write.
 */
void ncPut1Double(int ncid, int varid, const size_t index[], double *ip)
{
  int res; /* NetCDF return status code */

  if ((res = nc_put_var1_double(ncid, varid, index, ip))) ERR(res);
}

/**
 * @brief      Write an array of double values.
 *
 *             Wrapper for nc_put_vara_double() which writes an array of double
 *             values from an open NetCDF dataset that is in data mode.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  start  Index position to start at.
 * @param[in]  count  Vector of edge lengths.
 * @param[out] ip     Pointer to location of the values to write.
 */
void ncPutDouble(int ncid, int varid, const size_t start[], const size_t count[], double *ip)
{
  int res; /* NetCDF return status code */

  if ((res = nc_put_vara_double(ncid, varid, start, count, ip))) ERR(res);
}

/**
 * @brief      Close the NetCDF file.
 *
 * @param[in]  ncid  NetCDF ID.
 */
void ncClose(int ncid)
{
  int res; /* NetCDF return status code */

  if ((res = nc_close(ncid))) ERR(res);
}
