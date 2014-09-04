// rattle.cc

#include <ShapeConvexPolyhedron.h> 
#include <iostream>
#include "QuaternionMultiply.h"
#include <vector>

using namespace std;
using namespace hpmc;
using namespace hpmc::detail;

poly3d_verts verts;
vector<void *> v_scp;
vector<vec3<Scalar>* > v_offsets;
vector<void *> v_verlet_scp;
vector<vec3<Scalar>* > v_verlet_offsets;
unsigned int recursion_matrix[256][256][256];
Scalar resolution = 0.01;
Scalar verlet_sq = 9.0;

int readPOSFile();
void rattleShape(unsigned int shape_of_interest);
void visit(unsigned int shape_of_interest, int _i, int _j, int _k);

double box_x=0, box_y=0, box_z=0;


int main(){

  int n_shapes = readPOSFile();
  cerr << n_shapes << " shapes read.\n";

  for (int i=0; i<n_shapes; i++) rattleShape(i); // rattle each shape

} // end main


void rattleShape(unsigned int shape_of_interest){

  // clear old Verlet lists
  v_verlet_scp.clear();
  v_verlet_offsets.clear();

  for (int i=0; i<v_scp.size(); i++) { // build a verlet list around shape of interest

    if (i==shape_of_interest) continue; // cos we don't count ourself
    Scalar dx = v_offsets[i]->x - v_offsets[shape_of_interest]->x; Scalar dy = v_offsets[i]->y - v_offsets[shape_of_interest]->y; Scalar dz = v_offsets[i]->z - v_offsets[shape_of_interest]->z;
    if ((dx*dx + dy*dy + dz*dz) < verlet_sq) {
      v_verlet_scp.push_back(v_scp[i]);
      v_verlet_offsets.push_back(v_offsets[i]);
    }
  } // done building Verlet list

  // clear recursion_matrix
  for (int i=0; i<256; i++) for (int j=0; j<256; j++) for (int k=0; k<256; k++) recursion_matrix[i][j][k] = 0;

  // initiate recursive search of space
  visit(shape_of_interest, 128, 128, 128);

} // end rattleShape()
  

void visit(unsigned int shape_of_interest, int _i, int _j, int _k) {

  if (recursion_matrix[_i][_j][_k]) return; // already  visited
  recursion_matrix[_i][_j][_k] = 1; // mark the visited point

  ShapeConvexPolyhedron scp_of_interest = *((ShapeConvexPolyhedron *)v_scp[shape_of_interest]); // retrieve the shape object of interest

  for (int i_shape=0; i_shape<v_verlet_scp.size(); i_shape++) {

    ShapeConvexPolyhedron scp = *((ShapeConvexPolyhedron *)v_verlet_scp[i_shape]); // retrieve the shape object i_shape
    vec3<Scalar> r_ab; 
    r_ab.x = v_offsets[shape_of_interest]->x - v_verlet_offsets[i_shape]->x + (_i - 128) * resolution;
    r_ab.y = v_offsets[shape_of_interest]->y - v_verlet_offsets[i_shape]->y + (_j - 128) * resolution;
    r_ab.z = v_offsets[shape_of_interest]->z - v_verlet_offsets[i_shape]->z + (_k - 128) * resolution;
    unsigned int err;
    
    if (test_overlap<ShapeConvexPolyhedron,ShapeConvexPolyhedron>(r_ab, scp, scp_of_interest, err)) return;

  } // end for i_shape loop over verlet

  // no overlap was found, so the point is good. Record it and visit the neighbors.
  printf("%06d\t%f\t%f\t%f\n", shape_of_interest, v_offsets[shape_of_interest]->x + (_i - 128) * resolution, 
           v_offsets[shape_of_interest]->y + (_j - 128) * resolution, v_offsets[shape_of_interest]->z + (_k - 128) * resolution);

  // visit each neighbor
  visit(shape_of_interest, _i - 1, _j, _k); visit(shape_of_interest, _i + 1, _j, _k);
  visit(shape_of_interest, _i, _j - 1, _k); visit(shape_of_interest, _i, _j + 1, _k);
  visit(shape_of_interest, _i, _j, _k - 1); visit(shape_of_interest, _i, _j, _k + 1);

}
  

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
      cerr << "got box = " << box_x << " x " << box_y << " x " << box_z << endl;
      continue;
    }

    // did we get the shape declarator spec line?
    if (!strcmp(first, "shape")) {

      // this will grab rest of line, chop out the quotes
      char *tokens = strtok(NULL, "\"");
      char *shape_type = strtok(tokens, " ");
      if (!strcmp(shape_type, "poly3d")) {
	unsigned int n_vertices = strtol(strtok(NULL, " "), NULL, 10);
	verts.N = n_vertices;
	for(int i=0; i<n_vertices; i++){ // read triples
	  Scalar s1 = strtod(strtok(NULL, "\t"), NULL); 
	  Scalar s2 = strtod(strtok(NULL, "\t"), NULL); 
	  Scalar s3 = strtod(strtok(NULL, "\t"), NULL); 
	  verts.v[i] = vec3<Scalar>(s1, s2, s3);
	}

      } // end if poly3d

      continue; // keep reading the POS file... 

    } // end if shape

    // From here on, we assume we're reading shapes; if it doesn't look like a shapeID then toss it. 
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
    ShapeConvexPolyhedron *scp = new ShapeConvexPolyhedron(q, verts);
    v_scp.push_back((void *)scp);

    vec3<Scalar> *offset = new vec3<Scalar>(x, y, z);
    v_offsets.push_back(offset);
    // Making a BAD design decision here, to allocate dynamically and store as generic pointers
    // because something in HOOMD doesn't like for the ShapeConvexPolyhedron object to be copied into an array
    // Talk to Josh and see what he recommends, but for now just proof-of-concept
    
  } // end while !EOF

  return shapes_read;

} // end readPOSFile()

