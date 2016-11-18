/**
 * @file phi3D_fsm.h
 * @brief      Header file for 3D Fast Sweeping Method (FSM).
 *
 *             The basic parallel algorithm for FSM is taken from the paper
 *             published in the Journal of Computational Physics titled
 *             "A parallel fast sweeping method for the Eikonal equation" by
 *             Miles Detrixhe, Deferic Gibou, and Chohong Min.
 *
 * @author     Shrestha, Anup
 * @date       12 JUL 2016
 *
 * @see        http://www.sciencedirect.com/science/article/pii/S002199911200722X
 */

#ifndef PHI3D_FSM_H
#define PHI3D_FSM_H

#include "phi3D.h"

/**
 * @brief      Structure created for ease of passing variables for FSM
 *             calculations.
 *
 *             Variables in this struct are used during sweeping part of FSM to
 *             manage internal grid dimensions, sweep directions and the
 *             position of the node on the array.
 */
typedef struct {
	int firstLevel; /**< Starting level, is always equal to NDIMS */
	int lastLevel;  /**< Final level, is equals to sum of dimensions */
	int level;      /**< Current level FSM is running on */
	int start;      /**< Level number to start at */
	int end;        /**< Level number to end at */
	int incr;       /**< True: Increment, False: Decrement loop counter */
	int xDim;       /**< Length of x-dimension w/o exterior nodes */
	int yDim;       /**< Length of y-dimension w/o exterior nodes */
	int zDim;       /**< Length of z-dimension w/o exterior nodes */
	int xSweepOff;  /**< Sweeping offset x-direction used for rotation */
	int ySweepOff;  /**< Sweeping offset y-direction used for rotation */
	int zSweepOff;  /**< Sweeping offset z-direction used for rotation */
} SweepInfo;

/**
 * @brief         Fast Sweeping Method
 *
 * @param[in,out] pf    Pointer to Phi function
 * @param[in]     sn    Sweep iteration number (represents the direction for
 *                      sweeping)
 */
void runFSM(Phi *pf, int sn);

/**
 * @brief      Sweep information.
 *
 * @param[in]  x     Dimension length in x-direction.
 * @param[in]  y     Dimension length in y-direction.
 * @param[in]  z     Dimension length in z-direction.
 * @param[in]  sn    Sweep iteration number (represents the direction for
 *                   sweeping)
 *
 * @return     SweepInfo
 */
SweepInfo makeSweepInfo(int x, int y, int z, int sn);

#endif /* PHI3D_FSM_H */
