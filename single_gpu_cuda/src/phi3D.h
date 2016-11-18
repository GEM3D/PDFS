/**
 * @file phi3D.h
 * @brief      Header file for 3D Phi function.
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
 */
#ifndef PHI3D_H_
#define PHI3D_H_

#include <math.h>

#include "file_writer.h"

#define SPEED 1 /* Speed of propagation [Equation (1) in the report] */
#define DEFAULT_BORDER_LOCATION -1
#define DEFAULT_BORDER_DISTANCE INFINITY
#define DEFAULT_INTERIOR_DISTANCE 90000

/**
 * @brief      Structure of the Phi3D datatype.
 */
typedef struct {
	int     x;        /**< Length of x-dimension*/
	int     y;        /**< Length of y-dimension*/
	int     z;        /**< Length of z-dimension*/
	double  dx;       /**< Node spacing in the x-direction */
	double  dy;       /**< Node spacing in the y-direction */
	double  dz;       /**< Node spacing in the z-direction */
	int *   location; /**< Location values array */
	double *distance; /**< Distance values array */
	double  F;        /**< Speed of propagation*/
} Phi;

/**
 * @brief      Creates a phi function.
 *
 * @param[in]  f  Pointer to the VTI file.
 * @param[out] pf Pointer to phi function.
 */
void create_phi_function(FILE *f, Phi *pf);

/**
 * @brief      Deallocates memory allocated for the phi function
 *
 * @param[out] pf    Pointer to phi function.
 */
void destroy_phi_function(Phi *pf);

/**
 * @brief      Calculates the distance field.
 *
 * @param[in,out] pf    Pointer to phi function.
 */
void calc_dist_field(Phi *pf);

/**
 * @brief      Generates the output file for the distance field.
 *
 * @param[in]  pf     Pointer to phi function.
 * @param[in]  out    The type of output file.
 * @param[in]  fname  Filename of the output file.
 */
void phi_write_file(Phi *pf, FileOut out, char *fname);

#endif /* PHI3D_H_ */
