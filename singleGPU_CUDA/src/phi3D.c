/**
 * @file phi3D.c
 * @brief      Source file for 3D Phi function.
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
 */

#include "phi3D.h"
#include "vti_parser.h"

extern void run_fsm(Phi *, int);

// Private function definitions
static void update_distance(Phi *pf, int totalNodes);
static void set_distance_negative_inside(Phi *pf, int totalNodes);
static void adjust_boundary(Phi *pf);

/**
 * @brief      Creates phi function, calls the file parser methods to get the
 *             dimensions, location and distance data from the file. Also,
 *             allocates memory and initializes variables required by phi
 *             function.
 *
  * @param[in]  f  Pointer to the VTI file.
 * @param[out] pf Pointer to phi function.
 */
void create_phi_function(FILE *f, Phi *pf)
{
  // parse the file for dimension data
  double d[6];
  vti_get_dimensions(f, d);

  // initialize fields for the phi function
  pf->x  = (int) d[0] + 1;
  pf->dx = d[3];
  pf->y  = (int) d[1] + 1;
  pf->dy = d[4];
  pf->z  = (int) d[2] + 1;
  pf->dz = d[5];
  pf->F  = SPEED;

  // allocate memory for location and distance arrays
  int totalNodes = (pf->x + 2) * (pf->y + 2) * (pf->z + 2);
  pf->location   = (int *) malloc(sizeof(int) * totalNodes);
  pf->distance   = (double *) malloc(sizeof(double) * totalNodes);

  if (pf->location == NULL || pf->distance == NULL) {
    printf("Error allocating memory: %s:%d\n", __FILE__, __LINE__);
    exit(1);
  }

  // parse the file and fill the location
  // and distance arrays
  vti_get_data(f, pf->location, DEFAULT_BORDER_LOCATION, pf->distance, DEFAULT_BORDER_DISTANCE,
               make_grid3D(pf->x + 2, pf->y + 2, pf->z + 2));
}

/*
 * Deallocates memory allocated for the phi function.
 *
 * Arguments:
 *    Phi* [in/out] - pointer to phi function
 * Returns:
 *
 */

/**
 * @brief      Deallocates memory allocated for the phi function.
 *
 * @param[out] pf    Pointer to phi function.
 */
void destroy_phi_function(Phi *pf)
{
  free(pf->location);
  free(pf->distance);
}

/**
 * @brief      Calls appropriate functions required to calculate the distance
 *             field.
 *
 * @param[in,out] pf    Pointer to phi function.
 */
void calc_dist_field(Phi *pf)
{
  // get the total number of nodes
  // in the grid
  int totalNodes = (pf->x + 2) * (pf->y + 2) * (pf->z + 2);

  // update the distance values
  update_distance(pf, totalNodes);

  // use the fast sweeping method
  // to get the solution for the Eikonal Equation
  int itr = 1;
  run_fsm(pf, itr);

  // set the distance values to negative
  // for inside region
  set_distance_negative_inside(pf, totalNodes);

  adjust_boundary(pf);
}

/*
 * Calls the method that writes distance field
 * to different files based on the FileOut parameter
 *
 * Arguments:
 *      Phi* [in] - pointer to phi function
 *   FileOut [in] - Output file type enumerator
 *     char* [in] - name of the output file
 * Returns:
 *
 */

/**
 * @brief      Calls the method that writes distance field to different files
 *             based on the FileOut parameter
 *
 * @param[in]  pf     Pointer to phi function.
 * @param[in]  out    The type of output file.
 * @param[in]  fname  Filename of the output file.
 */
void phi_write_file(Phi *pf, FileOut out, char *fname)
{
  Grid3D   g3d = make_grid3D(pf->x + 2, pf->y + 2, pf->z + 2);
  FileGrid fg  = make_fileGrid(g3d.x, g3d.y, g3d.z, pf->dx, pf->dy, pf->dz);
  FileType ft  = make_fileType(fname, pf->distance, out, fg);

  // generate file(s)
  file_generate(&ft);
}

/**
 * @brief         Set statuses of the grid rim to converged. And sets unknown
 *                values of the grid points to default values.
 *
 * @param[in.out] pf          Pointer to phi function.
 * @param[in]     totalNodes  The total grid points.
 */
static void update_distance(Phi *pf, int totalNodes)
{
  int *   l = &pf->location[0];
  double *d = &pf->distance[0];

  int i;
  for (i = 0; i < totalNodes; i++) {
    if (*l != DEFAULT_BORDER_LOCATION && *d != DEFAULT_BORDER_DISTANCE) {
      *d = (*l == 1 && *d == DEFAULT_BORDER_DISTANCE) ? -1 : (*d > 0.0 || *d < 0.0)
                                                             ? *d
                                                             : DEFAULT_INTERIOR_DISTANCE;
    }
    l++;
    d++;
  }
}

/**
 * @brief         Sets distances for points inside the interface (location = 1)
 *                to -1. This function needs to be called once the sweeping is
 *                finished.
 *
 * @param[in,out] pf          Pointer to phi function.
 * @param[in]     totalNodes  The total grid points.
 */
static void set_distance_negative_inside(Phi *pf, int totalNodes)
{
  int *   l = &pf->location[0];
  double *d = &pf->distance[0];

  int i;
  for (i = 0; i < totalNodes; i++) {
    if (*l != DEFAULT_BORDER_LOCATION && *d != DEFAULT_BORDER_DISTANCE) {
      if (*l == 1) {
        *d = -1;
      }
    }
    l++;
    d++;
  }
}

/*
 * Adjusts the boundary values
 *
 * Arguments:
 *     Phi * [in/out] - pointer to phi function
 * Returns:
 *
 */

/**
 * @brief         Adjusts the boundary values
 *
 * @param[in,out] pf    Pointer to phi function.
 */
static void adjust_boundary(Phi *pf)
{
  int x, y, z, i, j, k, xy;
  x  = pf->x + 2;
  y  = pf->y + 2;
  z  = pf->z + 2;
  xy = x * y;

  for (i = 0; i < z; i++) {
    for (j = 0; j < y; j++) {
      for (k = 0; k < x; k++) {
        int I = i, J = j, K = k;
        I = (i == z - 1) ? I - 1 : (!i) ? I + 1 : I;
        J = (j == y - 1) ? J - 1 : (!j) ? J + 1 : J;
        K = (k == x - 1) ? K - 1 : (!k) ? K + 1 : K;
        if (i != I || j != J || k != K) {
          pf->distance[i * xy + j * x + k] = pf->distance[I * xy + J * x + K];
        }
      }
    }
  }
}
