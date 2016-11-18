/**
 * @file phi3D_fsm.c
 * @brief      Source file for 3D Fast Sweeping Method (FSM).
 *
 *             The basic parallel algorithm for FSM is taken from the paper
 *             published in the Journal of Computational Physics titled
 *             "A parallel fast sweeping method for the Eikonal equation" by
 *             Miles Detrixhe, Deferic Gibou, and Chohong Min.
 *
 * @author     Shrestha, Anup
 * @date       12 JUL 2016
 *
 * @see http://www.sciencedirect.com/science/article/pii/S002199911200722X
 */
#include "phi3D_fsm.h"
#include "utilities.h"

#include <stdio.h>
#include <stdlib.h>

/* Private method definitions */
static void fastSweep(Phi *pf, int sn);
static double solveEikonal(Phi *pf, int index, int extX, int extXY);
/*-----------------------------------------*/

/**
 * @brief         Fast Sweeping Method
 *
 *                Executes the Fast Sweeping Method routine to calculate the
 *                distance field for a 3D grid by solving the Eikonal equation.
 *                Sweeping Directions:
 *                  (1) i = 1:I, j = 1:J, k = 1:K
 *                  (2) i = I:1, j = 1:J, k = K:1
 *                  (3) i = I:1, j = 1:J, k = 1:K
 *                  (4) i = 1:I, j = 1:J, k = K:1
 *                  (5) i = I:1, j = J:1, k = K:1
 *                  (6) i = 1:I, j = J:1, k = 1:K
 *                  (7) i = 1:I, j = J:1, k = K:1
 *                  (8) i = I:1, j = J:1, k = 1:K
 *
 * @param[in,out] pf    Pointer to Phi function
 * @param[in]     sn    Sweep iteration number (represents the direction for
 *                      sweeping)
 */
void runFSM(Phi *pf, int sn) { fastSweep(pf, sn); }
/**
 * @brief         Fast Sweeping Method (FSM)
 *
 *                Calculates the distance field for a 3D grid by solving the
 *                Eikonal equation at each grid point using the parallel Fast
 *                Sweeping Method. For 3D there are 8 different sweeping
 *                directions each represented by a sweep order counter from 1-8,
 *                the sweeping is performed by alternating the start and end
 *                levels and translating the coordinates based on the value of
 *                the sweep order counter
 *
 * @param[in,out] pf    Pointer to Phi Function
 * @param[in]     sn    Sweep iteration number (represents the direction for
 *                      sweeping)
 */
static void fastSweep(Phi *pf, int sn)
{
  /* Length of x-dimension including exterior nodes */
  int extX = pf->x + 2;
  /* Number of nodes in the XY-plane including exterior nodes */
  int extXY = extX * (pf->y + 2);

  SweepInfo s = makeSweepInfo(pf->x, pf->y, pf->z, sn);

  for (int level = s.start; level != s.end; level = (s.incr) ? level + 1 : level - 1) {
    int xs = max(1, level - (s.yDim + s.zDim));
    int ys = max(1, level - (s.xDim + s.zDim));
    int xe = min(s.xDim, level - (s.firstLevel - 1));
    int ye = min(s.yDim, level - (s.firstLevel - 1));

    int xSO = s.xSweepOff;
    int ySO = s.ySweepOff;
    int zSO = s.zSweepOff;
#pragma acc kernels
    {
#pragma acc loop independent
      for (int x = xs; x <= xe; x++) {
        int k = abs(x - xSO);
#pragma acc loop independent
        for (int y = ys; y <= ye; y++) {
          int z = level - (x + y);
          if (z > 0 && z <= pf->z) {
            int j             = abs(y - ySO);
            int i             = abs(z - zSO);
            int idx           = i * extXY + j * extX + k; /* linearized 3D index */
            pf->distance[idx] = solveEikonal(pf, idx, extX, extXY);
          }
        }
      } /* End of acc kernels */
    }
  }
}

/**
 * @brief      Eikonal solver.
 *
 *             Solves Eikonal equation at linearized 3D index. Returns the
 *             minimum of calculated and previous distance values.
 *
 * @param      pf     Pointer to Phi function.
 * @param[in]  index  Linear index of 3D grid.
 *
 * @return     Minimum of calculated and previos solutio at the given index.
 */
#pragma acc routine seq
static double solveEikonal(Phi *pf, int index, int extX, int extXY)
{
  double dist_new = 0;
  double dist_old = pf->distance[index];

  double minX = min(pf->distance[index - 1], pf->distance[index + 1]);
  double minY = min(pf->distance[abs(index - extX)], pf->distance[abs(index + extX)]);
  double minZ = min(pf->distance[abs(index - extXY)], pf->distance[abs(index + extXY)]);

  double m[] = {minX, minY, minZ};
  double d[] = {pf->dx, pf->dy, pf->dz};

  // sort the mins
  int    i, j;
  double tmp_m, tmp_d;
  for (i = 1; i < 3; i++) {
    for (j = 0; j < 3 - i; j++) {
      if (m[j] > m[j + 1]) {
        tmp_m    = m[j];
        tmp_d    = d[j];
        m[j]     = m[j + 1];
        d[j]     = d[j + 1];
        m[j + 1] = tmp_m;
        d[j + 1] = tmp_d;
      }
    }
  }

  // simplifying the variables
  double m2_0 = m[0] * m[0], m2_1 = m[1] * m[1], m2_2 = m[2] * m[2];
  double d2_0 = d[0] * d[0], d2_1 = d[1] * d[1], d2_2 = d[2] * d[2];

  dist_new = m[0] + d[0];
  if (dist_new > m[1]) {
    double s = sqrt(-m2_0 + 2 * m[0] * m[1] - m2_1 + d2_0 + d2_1);
    dist_new = (m[1] * d2_0 + m[0] * d2_1 + d[0] * d[1] * s) / (d2_0 + d2_1);

    if (dist_new > m[2]) {
      double a = sqrt(-m2_0 * d2_1 - m2_0 * d2_2 + 2 * m[0] * m[1] * d2_2 - m2_1 * d2_0
                      - m2_1 * d2_2 + 2 * m[0] * m[2] * d2_1 - m2_2 * d2_0 - m2_2 * d2_1
                      + 2 * m[1] * m[2] * d2_0 + d2_0 * d2_1 + d2_0 * d2_2 + d2_1 * d2_2);

      dist_new = (m[2] * d2_0 * d2_1 + m[1] * d2_0 * d2_2 + m[0] * d2_1 * d2_2
                  + d[0] * d[1] * d[2] * a)
                 / (d2_0 * d2_1 + d2_0 * d2_2 + d2_1 * d2_2);
    }
  }

  return min(dist_old, dist_new);
}

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
SweepInfo makeSweepInfo(int x, int y, int z, int sn)
{
  SweepInfo s;

  s.firstLevel = 3;
  s.lastLevel  = x + y + z;

  s.xDim = x;
  s.yDim = y;
  s.zDim = z;

  s.xSweepOff = (sn == 4 || sn == 8) ? x + 1 : 0;
  s.ySweepOff = (sn == 2 || sn == 6) ? y + 1 : 0;
  s.zSweepOff = (sn == 3 || sn == 7) ? z + 1 : 0;

  if (sn == 2 || sn == 5 || sn == 7 || sn == 8) {
    s.start = s.lastLevel;
    s.end   = s.firstLevel - 1;
    s.incr  = 0;
  } else {
    s.start = s.firstLevel;
    s.end   = s.lastLevel + 1;
    s.incr  = 1;
  }

  return s;
}
