// rattle.cc

#include <ShapeConvexPolyhedron.h> 
#include <iostream>
#include "QuaternionMultiply.h"
#include <vector>

using namespace std;
using namespace hpmc;
using namespace hpmc::detail;

poly3d_verts tetrahedron_verts;
vector<void *> v_scp;
vector<vec3<Scalar>* > v_offsets;
vector<void *> v_verlet;

int readPOSFile();

double box_x=0, box_y=0, box_z=0;

int main(){

  // create the poly3d_verts object to represent tetrahedra
  // poly3d_verts tetrahedron_verts;
  tetrahedron_verts.N = 4;
  tetrahedron_verts.v[0] = vec3<Scalar>(1,1,1);
  tetrahedron_verts.v[1] = vec3<Scalar>(-1,-1,1);
  tetrahedron_verts.v[2] = vec3<Scalar>(-1,1,-1);
  tetrahedron_verts.v[3] = vec3<Scalar>(1,-1,-1);

  // Create a reference shape to probe with...
  quat<Scalar> q(1, vec3<Scalar>(1, 1, 1)); // not a real representation, but just to test
  ShapeConvexPolyhedron scp_probe(getQuaternion(0, vec3<Scalar>(1,0,0)), tetrahedron_verts);
  //ShapeConvexPolyhedron scp_probe(q, tetrahedron_verts);

  // read the .pos file
  int n_shapes = readPOSFile();

  cout << n_shapes << " shapes read.\n";

  rattleShape(unsigned int shape_number); 

/*
  // check overlap with dummy probe
  for (int i=0; i<n_shapes; i++) {
    ShapeConvexPolyhedron *p_scp = (ShapeConvexPolyhedron *)(v_scp[i]);
    ShapeConvexPolyhedron scp = *p_scp;
    vec3<Scalar>r_ab = *v_offsets[i];
    unsigned int err;
    if (test_overlap<ShapeConvexPolyhedron,ShapeConvexPolyhedron>(r_ab, scp_probe, scp, err)) cout << "bump! " << i << endl;
  }   
*/

  // replicate to form supercell
  // loop through each shape in main cell
  // translate shape in recursive search of space
  // give output for each contiguous point found
  // if (!test_overlap<ShapeConvexPolyhedron,ShapeConvexPolyhedron>(r_ab, shape1, shape2, err)) break;

} // end main

void rattleShape(unsigned int shape_number){

  Scalar x = v_offsets[shape_number].x;
  Scalar y = v_offsets[shape_number].y;
  Scalar z = v_offsets[shape_number].z;

  for (int i=0; i<n_shapes; i++) { // build a verlet list around shape of interest

    if (i==shape_number) continue; // cos we don't count ourself
    Scalar dx = v_offsets[i].x - x; Scalar dy = v_offsets[i].y - y; Scalar dz = v_offsets[i].z - z;
    if ((dx*dx + dy*dy + dz*dz) < verlet_sq) v_verlet.push_back(v_scp[i]);

  } // done building Verlet list

  // loop over set of rattlers
  //for (rattler = 0; rattler < number_of_rattlers; rattler++)
  //{
  
  int i, j, k; // indices for recursion matrix
  
  // clear recursion_matrix
  for (i=0; i<256; i++) for (j=0; j<256; j++) for (k=0; k<256; k++) recursion_matrix[i][j][k] = 0;
  // start recursion with rattler 0
  //                         test_x = x[rattler];
  //                             test_y = y[rattler];
  //                                 test_z = z[rattler];
  //                                     makeVerletList();
  //                                         visit(128, 128, 128);
  //                                           }
  //
  //                                             return 0;
  //                                             }
  //
  // visiting implies that site has already been determined to be included
  void visit(int _i, int _j, int _k)
  {
  recursion_matrix[_i][_j][_k] = 1; // mark the visited point
  test_x = verlet_center_x + (_i - 128) * resolution;
  test_y = verlet_center_y + (_j - 128) * resolution;
  test_z = verlet_center_z + (_k - 128) * resolution;
  printf("%06d\t%f\t%f\t%f\n", rattler, test_x, test_y, test_z);
  
  // check each neighbor for inclusion, and if included, then visit it...
  test_x = verlet_center_x + (_i - 129) * resolution;
  if (!recursion_matrix[_i - 1][_j][_k] && checkInclusion(test_x, test_y, test_z)) visit(_i - 1, _j, _k);
  test_x = verlet_center_x + (_i - 127) * resolution;
  if (!recursion_matrix[_i + 1][_j][_k] && checkInclusion(test_x, test_y, test_z)) visit(_i + 1, _j, _k);
  
  test_x = verlet_center_x + (_i - 128) * resolution;
  test_y = verlet_center_y + (_j - 129) * resolution;
  if (!recursion_matrix[_i][_j - 1][_k] && checkInclusion(test_x, test_y, test_z)) visit(_i, _j - 1, _k);
  test_y = verlet_center_y + (_j - 127) * resolution;
  if (!recursion_matrix[_i][_j + 1][_k] && checkInclusion(test_x, test_y, test_z)) visit(_i, _j + 1, _k);
  
  test_y = verlet_center_y + (_j - 128) * resolution;
  test_z = verlet_center_z + (_k - 129) * resolution;
  if (!recursion_matrix[_i][_j][_k - 1] && checkInclusion(test_x, test_y, test_z)) visit(_i, _j, _k - 1);
  test_z = verlet_center_z + (_k - 127) * resolution;
  if (!recursion_matrix[_i][_j][_k + 1] && checkInclusion(test_x, test_y, test_z)) visit(_i, _j, _k + 1);
  }
  
  int checkInclusion(double test_x, double test_y, double test_z)
  {
  printf("Checking inclusion of %lf, %lf, %lf\n", test_x, test_y, test_z);
  int i;
  double dx, dy, dz, dd;
  for (i=0; i < close_rattlers; i++)
  {
  dx = test_x - close_x[i];
  dy = test_y - close_y[i];
  dz = test_z - close_z[i];
  dd = dx*dx + dy*dy + dz*dz;
  if (dd < diameter_sq) return 0; // inside the diameter of another, so we don't count it.
  }
  
  return 1; // it's not inside of any of the surrounding rattlers, so it's 'free' 
  }
  
 

} // end rattleShape()

int readPOSFile(){

  // read the pos file
  char line[256];
  char *first, *shape_x, *shape_y, *shape_z, *quat_s, *quat_x, *quat_y, *quat_z;
  double x,y,z,qs,qx,qy,qz;
  int shapes_read=0;

  while (1) { // loop over all lines

    // grab a line; if end of file, then move on.
    fgets(line, 256, stdin);
    if (feof(stdin)) break;

    // get first token
    first = strtok(line, "\t");
    
    // did we get the box size spec line?
    if (!strcmp(first, "box")) {
      box_x = strtod(strtok(NULL, "\t"), NULL);
      box_y = strtod(strtok(NULL, "\t"), NULL);
      box_z = strtod(strtok(NULL, "\n"), NULL);
      cout << "got box = " << box_x << " x " << box_y << " x " << box_z << endl;
      continue;
    }
    // if it doesn't look like a shapeID then toss it. 
    if (first[0] != 'f' || first[1] != 'f') continue;

    shapes_read++;

    // grab the position and quaternion components
    shape_x = strtok(NULL, "\t"); shape_y = strtok(NULL, "\t"); shape_z = strtok(NULL, "\t");
    quat_s  = strtok(NULL, "\t"); quat_x  = strtok(NULL, "\t"); quat_y  = strtok(NULL, "\t"); quat_z  = strtok(NULL, "\n");

    // convert to floating point
    x = strtod(shape_x, NULL); y = strtod(shape_y, NULL); z = strtod(shape_z, NULL);
    qs = strtod(quat_s, NULL); qx = strtod(quat_x, NULL); qy = strtod(quat_y, NULL); qz = strtod(quat_z, NULL);

    // construct the quaternion, shape, and offset
    quat<Scalar> q(qs, vec3<Scalar>(qx, qy, qz));
    ShapeConvexPolyhedron *scp = new ShapeConvexPolyhedron(q, tetrahedron_verts);
    v_scp.push_back((void *)scp);

    vec3<Scalar> *offset = new vec3<Scalar>(x, y, z);
    v_offsets.push_back(offset);
    // Making a BAD design decision here, to allocate dynamically and store as generic pointers
    // because something in HOOMD doesn't like for the ShapeConvexPolyhedron object to be copied into an array
    // Talk to Josh and see what he recommends, but for now just proof-of-concept
    
  } // end while !EOF

  return shapes_read;

} // end readPOSFile()

