/**
 * @file phi3D.c
 * @brief      Source file for 3D Phi function.
 *
 * @author     Shrestha, Anup
 * @date       12 JUL 2016
 */

#include "phi3D.h"

#include <stdio.h>
#include <stdlib.h>
/**
 * @brief      Create the Phi3D data structure.
 *
 *             Assigns proper attributes to the Phi3D struct and allocates
 *             memory for the location and distance attribute. Returns 0 if the
 *             memory allocation failed.
 *
 * @note       This function makes the assumption that the dimension length
 *             includes the exterior nodes and therefore subtracts 2 to assign
 *             the dimension length of the interior nodes. The index 0,1,2 for
 *             dim param is the z,y,x-dimension.
 *
 * @param[in]  dim   Array of dimension length.
 * @param[in]  dlt   Array of delta values.
 * @param[out] pf    Pointer to location of Phi function.
 *
 * @return     1 - Success; 0 - Failure
 */
int makePhi3D(size_t *dim, double *dlt, Phi *pf)
{
	pf->x = dim[2] - 2;
	pf->y = dim[1] - 2;
	pf->z = dim[0] - 2;

	pf->dx = dlt[0];
	pf->dy = dlt[1];
	pf->dz = dlt[2];

	pf->F = SPEED;

	int totalNodes = dim[0] * dim[1] * dim[2];
	pf->distance   = (double *) malloc(sizeof(double) * totalNodes);
	if (pf->distance == NULL) {
		printf("Error allocating memory: %s:%d\n", __FILE__, __LINE__);
		return 0;
	}

	return 1;
}

/**
 * @brief      Destroy the Phi3D data structure.
 *
 *             Deallocates memory allocated by the makePhi3D function.
 *
 * @param[out] pf    Pointer to location of Phi function.
 *
 * @return     1 - Success; 0 - Failure
 */
int destroyPhi3D(Phi *pf)
{
	free(pf->distance);

	return 1;
}
