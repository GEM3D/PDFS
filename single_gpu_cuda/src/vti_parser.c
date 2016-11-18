
/**
 * @file vti_parser.c
 * @brief      Source file for parsing VTI input file. Parses information from
 *             the VTI file and stores them in the appropriate output
 *             parameters. Also responsible for setting the border values for
 *             the arrays.
 *
 * @author     Shrestha, Anup
 * @date       09 OCT 2015
 */

#include "vti_parser.h"

// private method declarations
static void move_file_pointer(FILE *file_ptr, int lineNumber, int r);
static void get_location(FILE *vti, int *l, int b_l, Grid3D g);
static void get_distance(FILE *vti, double *d, double b_d, Grid3D g);

/**
 * @brief      Creates a Grid3D struct.
 *
 * @param[in]  x     x-dimensions.
 * @param[in]  y     y-dimensions.
 * @param[in]  z     z-dimensions.
 *
 * @return     Struct representing a 3D grid.
 */
Grid3D make_grid3D(int x, int y, int z)
{
  Grid3D g;
  g.x = x;
  g.y = y;
  g.z = z;

  return g;
}

/**
 * @brief          Moves the pointer of a FILE to specified line.
 *
 * @param[in, out] file_ptr    The file pointer to be moved.
 * @param[in]      lineNumber  The line number to move the pointer to.
 * @param[in]      r           0 - rewind; 1 - no rewind.
 */
static void move_file_pointer(FILE *file_ptr, int lineNumber, int r)
{
  char tmpStr[512];
  if (r) rewind(file_ptr);
  while (lineNumber > 0) {
    fgets(tmpStr, 511, file_ptr);
    lineNumber--;
  }
}

/**
 * @brief      Searches for x, y, z, dx, dy, dz information in the VTI file and
 *             puts it into the array taken as the second argument.
 *
 * @param[in]  vti   Pointer to a vti input file.
 * @param[out] d     Array to store dimension and spacing data.
 */
void vti_get_dimensions(FILE *vti, double *d)
{
  char tmpStr[512];
  rewind(vti);
  while (1) {
    fgets(tmpStr, 511, vti);
    if (strstr(tmpStr, "ImageData WholeExtent")) {
      sscanf(tmpStr, "    <ImageData WholeExtent=\"0 %lf 0 %lf 0 %lf\" "
                     "Spacing=\"%lf %lf %lf\">",
             &d[0], &d[1], &d[2], &d[3], &d[4], &d[5]);
      break;
    }
  }
}

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
void vti_get_data(FILE *vti, int *l, int b_l, double *d, double b_d, Grid3D g)
{
  // move the file pointer to
  // line 6 from beginning
  move_file_pointer(vti, 6, 1);

  get_location(vti, l, b_l, g);

  // move the file pointer 2 lines
  // forward from its last position
  move_file_pointer(vti, 2, 0);

  get_distance(vti, d, b_d, g);
}

/**
 * @brief      Parses the file for location values and stores them in the int
 *             array. Also adds the default border values for location in the
 *             array.
 *
 * @param[in]  vti   Pointer to a vti input file.
 * @param[out] l     Pointer to an array to store location values.
 * @param[in]  b_l   Default border value for location.
 * @param[in]  g     Grid dimensions.
 */
static void get_location(FILE *vti, int *l, int b_l, Grid3D g)
{
  int i, j, k, *t = &l[0];
  for (i = 0; i < g.z; i++) {
    for (j = 0; j < g.y; j++) {
      for (k = 0; k < g.x; k++) {
        // Border
        if (k == 0 || k == g.x - 1 || j == 0 || j == g.y - 1 || i == 0 || i == g.z - 1) {
          *(t++) = b_l;
        } else { // Interior
          fscanf(vti, "%d ", t++);
        }
      }
    }
  }
}

/**
 * @brief      Parses the file for distance values and stores them in the double
 *             array. Also adds the default border values for distance in the
 *             array.
 *
 * @param[in]  vti   Pointer to a vti input file.
 * @param[out] d     Pointer to an array to store distance values.
 * @param[in]  b_d   Default border value for distance.
 * @param[in]  g     Grid dimensions.
 */
static void get_distance(FILE *vti, double *d, double b_d, Grid3D g)
{
  int     i, j, k;
  double *t = &d[0];
  for (i = 0; i < g.z; i++) {
    for (j = 0; j < g.y; j++) {
      for (k = 0; k < g.x; k++) {
        // Border distance
        if (k == 0 || k == g.x - 1 || j == 0 || j == g.y - 1 || i == 0 || i == g.z - 1) {
          *(t++) = b_d;
        } else { // Interior distance
          fscanf(vti, "%lf ", t++);
        }
      }
    }
  }
}
