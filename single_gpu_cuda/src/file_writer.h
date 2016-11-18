/**
 * @file file_writer.h
 * @brief      Header file for generating output files.
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
 */
#ifndef FILEWRITER_H
#define FILEWRITER_H

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief      Enumerator for types of output files
 *             1) VTI - VTK Image Data
 *             2) DAT - Data file
 */
typedef enum { VTI = 0x01, DAT = 0x02 } FileOut;

/**
 * @brief      Structure for storing the dimenstion and spacing of the grid.
 */
typedef struct {
  int    x;
  int    y;
  int    z;
  double dx;
  double dy;
  double dz;
} FileGrid;

/**
 * @brief      Structure for information required to generate an output file.
 */
typedef struct {
  char *   name; /**< name of the output file */
  double * data; /**< data to be put into file */
  FileOut  out;  /**< type of file to output */
  FileGrid grid; /**< grid dimensions/spacing */
} FileType;

// NetCDF Macros
// Handle errors by printing an error message and exiting with a
// non-zero status.
#define ERRCODE 2
#define ERR(e)                                                                                     \
  {                                                                                              \
    printf("Error: %s\n", nc_strerror(e));                                                     \
    exit(ERRCODE);                                                                             \
  }

/**
 * @brief      Creates a FileGrid data type.
 *
 * @param[in]  x     Length of x-dimension of grid.
 * @param[in]  y     Length of y-dimension of grid.
 * @param[in]  z     Length of z-dimension of grid.
 * @param[in]  dx    Grid spacing in x-dimension.
 * @param[in]  dy    Grid spacing in y-dimension.
 * @param[in]  dz    Grid spacing in z-dimension.
 *
 * @return     Struct with dimension and spacing information.
 */
FileGrid make_fileGrid(int x, int y, int z, double dx, double dy, double dz);

/**
 * @brief      Create a FileType struct.
 *
 * @param[in]  n     Name of the file.
 * @param[in]  d     Data to be written.
 * @param[in]  o     Type of file to generate.
 * @param[in]  g     Dimensions/spacing of the grid.
 *
 * @return     Struct containing all the information to generate a file.
 */
FileType make_fileType(char *n, double *d, FileOut o, FileGrid g);

/**
 * @brief      Determines what types of file to output generates filename for
 *             each type and calls the appropriate methods to generate the file.
 *
 * @param[in]  f     Pointer to FileType.
 */
void file_generate(FileType *f);

#endif /* FILEWRITER_H */
