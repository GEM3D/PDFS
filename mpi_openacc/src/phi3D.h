/**
 * @file phi3D.h
 * @brief      Header file for 3D Phi function.
 *
 * @author     Shrestha, Anup
 * @date       12 JUL 2016
 */

#ifndef PHI3D_H
#define PHI3D_H

#include <math.h>
#include <stddef.h>

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
	double *distance; /**< Distance values array */
	double  F;        /**< Speed of propagation*/
} Phi;

/**
 * @brief      Create the Phi3D data structure.
 *
 * @param[in]  dim   Array of dimension length.
 * @param[in]  dlt   Array of delta values.
 * @param[out] pf    Pointer to location of Phi function.
 *
 * @return     1 - Success; 0 - Failure
 */
int makePhi3D(size_t *dim, double *dlt, Phi *pf);

/**
 * @brief      Destroy the Phi3D data structure.
 *
 * @param[out] pf    Pointer to location of Phi function.
 *
 * @return     1 - Success; 0 - Failure
 */
int destroyPhi3D(Phi *pf);

#endif /* PHI3D_H */
