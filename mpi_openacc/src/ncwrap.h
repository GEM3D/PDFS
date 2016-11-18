/**
 * @file ncwrap.h
 * @brief      Header file for the wrapper functions to NetCDF-4 API calls.
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
 * https://www.unidata.ucar.edu/software/netcdf/netcdf-4/newdocs/netcdf-c.pdf
 */

#ifndef NCWRAP_H
#define NCWRAP_H

#include <netcdf.h>
#include <stddef.h>
#include <stdio.h>

/**
 * @brief      NetCDF API call error handler.Calculates the number of seconds
 *             with microsecond accuracy.
 *
 *             Handle NetCDF API call errors by printing an error message and
 *             exiting with a non-zero status.
 *
 * @param[in]  e     Error-code
 */
#define ERR(e)                                                                                     \
	{                                                                                              \
		printf("\nError: %s: %d %s\n", __FILE__, __LINE__, nc_strerror(e));                        \
		exit(2);                                                                                   \
	}

#define NDIMS 3 /* Number of spatial dimensions */

/**
 * @brief      Open a parallel NetCDF-4 file.
 *
 * @param[in]  filename  The file name of the NetCDF dataset to be opened.
 * @param[out] ncid      Pointer to location where returned NetCDF id is to be
 */
void ncOpen(char *filename, int *ncid);

/**
 * @brief      Create a parallel NetCDF-4 file
 *
 * @param[in]  filename  The file name of the new NetCDF dataset.
 * @param[out] ncid      Pointer to location where returned NetCDF id is to be
 *                       stored.
 */
void ncCreate(char *filename, int *ncid);

/**
 * @brief      Get the id of a NetCDF variable, given its name.
 *
 * @param[in]  ncid     NetCDF ID.
 * @param[in]  varName  The name of the variable to look for.
 * @param[out] varid    Pointer to location where returned NetCDF variable id is
 *                      to be stored.
 */
void ncGetVarid(int ncid, char *varName, int *varid);

/**
 * @brief      Get dimensions and delta values.
 *
 * @param[in]  ncid  NetCDF ID.
 * @param[out] dim   Pointer to location where length of the dimensions are to
 *                   be stored.
 * @param[out] dlt   Pointer to location where values of the delta attribute are
 *                   to be stored.
 */
void ncGetDim(int ncid, size_t *dim, double *dlt);

/**
 * @brief      Define and populate dimension and delta values.
 *
 * @param[in]  ncid    NetCDF ID.
 * @param[out] dimids  Pointer to location where returned NetCDF dimension ids
 *                     are to be stored.
 * @param[out] dim     Length of the dimensions.
 * @param[out] dlt     Delta values for each dimension.
 */
void ncPutDim(int ncid, int *dimids, size_t *dim, double *dlt);

/**
 * @brief      NetCDF parallel independent access.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 */
void ncParAccess(int ncid, int varid);

/**
 * @brief      Define variable.
 *
 * @param[in]  ncid     NetCDF ID.
 * @param[in]  varName  Variable name.
 * @param[in]  xtype    Predefined NetCDF external data type.
 * @param[in]  dimids   Dimension IDs.
 * @param[out] varid    Pointer to location for the returned variable ID to be
 *                      stored.
 */
void ncDefVar(int ncid, char *varName, nc_type xtype, const int dimids[], int *varid);

/**
 * @brief      Read a single integer value.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  index  Index of the value to get.
 * @param[out] ip     Pointer to location for the returned value to be stored.
 */
void ncGet1Int(int ncid, int varid, const size_t index[], int *ip);

/**
 * @brief      Read an array of integer values.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  start  Index position to start at.
 * @param[in]  count  Vector of edge lengths.
 * @param[out] ip     Pointer to location for the returned values to be stored.
 */
void ncGetInt(int ncid, int varid, const size_t start[], const size_t count[], int *ip);

/**
 * @brief      Write a single integer value.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  index  Index where the value to write.
 * @param[out] ip     Pointer to location of the value to write.
 */
void ncPut1Int(int ncid, int varid, const size_t index[], int *ip);

/**
 * @brief      Write an array of integer values.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  start  Index position to start at.
 * @param[in]  count  Vector of edge lengths.
 * @param[out] ip     Pointer to location of the values to write.
 */
void ncPutInt(int ncid, int varid, const size_t start[], const size_t count[], int *ip);

/**
 * @brief      Read a single double value.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  index  Index of the value to get.
 * @param[out] ip     Pointer to location for the returned value to be stored.
 */
void ncGet1Double(int ncid, int varid, const size_t index[], double *ip);

/**
 * @brief      Read an array of double values.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  start  Index position to start at.
 * @param[in]  count  Vector of edge lengths.
 * @param[out] ip     Pointer to location for the returned values to be stored.
 */
void ncGetDouble(int ncid, int varid, const size_t start[], const size_t count[], double *ip);

/**
 * @brief      Write a single double value.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  index  Index where the value to write.
 * @param[out] ip     Pointer to location of the value to write.
 */
void ncPut1Double(int ncid, int varid, const size_t index[], double *ip);

/**
 * @brief      Write an array of double values.
 *
 * @param[in]  ncid   NetCDF ID.
 * @param[in]  varid  Variable ID.
 * @param[in]  start  Index position to start at.
 * @param[in]  count  Vector of edge lengths.
 * @param[out] ip     Pointer to location of the values to write.
 */
void ncPutDouble(int ncid, int varid, const size_t start[], const size_t count[], double *ip);

/**
 * @brief      Close the NetCDF file.
 *
 * @param[in]  ncid  NetCDF ID.
 */
void ncClose(int ncid);

#endif /* NCWRAP_H */
