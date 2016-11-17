/**
 * @file vti_parser.h
 * @brief      Header file for parsing VTI input file.
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
 */

#ifndef VTI_Parser_H
#define VTI_Parser_H

#include <stdio.h>
#include <string.h>

/**
 * @brief      Structure for storing the dimension and spacing of the grid.
 */
typedef struct {
	int x;
	int y;
	int z;
} Grid3D;

/**
 * @brief      Creates a Grid3D struct.
 *
 * @param[in]  x     x-dimensions.
 * @param[in]  y     y-dimensions.
 * @param[in]  z     z-dimensions.
 *
 * @return     Struct representing a 3D grid.
 */
Grid3D make_grid3D(int x, int y, int z);

/**
 * @brief      Searches for x, y, z, dx, dy, dz information in the VTI file and
 *             puts it into the array taken as the second argument.
 *
 * @param[in]  vti   Pointer to a vti input file.
 * @param[out] d     Array to store dimension and spacing data.
 */
void vti_get_dimensions(FILE *vti, double *d);

/**
 * @brief      Calls the appropriate methods to parse the file for location and
 *             distance values and stores them in an int and double array
 *             respectively.
 *
 * @param[in]  vti   Pointer to a vti input file.
 * @param[out] l     Pointer to an array to store location values.
 * @param[in]  b_l   Default border value for location.
 * @param[out] d     Pointer to an array to store distance values.
 * @param[in]  b_d   Default border value for distance.
 * @param[in]  g     Grid dimensions.
 */
void vti_get_data(FILE *vti, int *l, int b_l, double *d, double b_d, Grid3D g);

#endif /* VTI_Parser_H */
