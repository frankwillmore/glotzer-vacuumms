// rattle.cc

#include <ShapeConvexPolyhedron.h> 
#include <iostream>
#include "QuaternionMultiply.h"

#include <vector>

using namespace std;
using namespace hpmc;
using namespace hpmc::detail;

//quat<Scalar> getQuaternion(Scalar angle, vec3<Scalar> axis){
//  Scalar st2 = sin(angle/2);
//  return quat<Scalar>(cos(angle/2), vec3<Scalar>(axis.x*st2, axis.y*st2, axis.z*st2));
//}

//poly3d_verts p3dv;
//list<ShapeConvexPolyhedron*> list_of_shapes;
//list<vec3<Scalar>* > list_of_positions;
vector<void *> v_scp;
vector<void *> v_offsets;

int readPOSFile();

double box_x=0, box_y=0, box_z=0;

int main(){

  // create the poly3d_verts object to represent tetrahedra
  poly3d_verts p3dv;
  p3dv.N = 4;
  p3dv.v[0] = vec3<Scalar>(1,1,1);
  p3dv.v[1] = vec3<Scalar>(-1,-1,1);
  p3dv.v[2] = vec3<Scalar>(-1,1,-1);
  p3dv.v[3] = vec3<Scalar>(1,-1,-1);

  // read the .pos file
  int n_shapes = readPOSFile();

  cout << n_shapes << " shapes read.\n";

  // replicate to form supercell
  
  // loop through each shape in main cell
  
  // translate shape in recursive search of space
  
  // give output for each contiguous point found
  // if (!test_overlap<ShapeConvexPolyhedron,ShapeConvexPolyhedron>(r_ab, shape1, shape2, err)) break;
  

  // when I want to get them back
  list<ShapeConvexPolyhedron*>::iterator i;
  //for(i=list_of_shapes.begin(); i != list_of_shapes.end(); ++i) cout << *i << " ";
  for(i=list_of_shapes.begin(); i != list_of_shapes.end(); ++i) cout << *i->quat.s << " ";



  vector<void *> v_pvoid;
  v_pvoid.push_back((void*)&scp);
  v_pvoid.push_back((void*)&scp);
  v_pvoid.push_back((void*)&scp);

  for( vector<void *>::iterator it = v_pvoid.begin(); it != v_pvoid.end(); ++it) cout << "x" << endl;



} // end main

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
cout << box_x << " x " << box_y << " x " << box_z << endl;
      continue;
    }
    // if it doesn't look like a shapeID then toss it. 
    if (first[0] != 'f' || first[1] != 'f') continue;

    shapes_read++;

    // grab the position and quaternion components
    shape_x = strtok(NULL, "\t");
    shape_y = strtok(NULL, "\t");
    shape_z = strtok(NULL, "\t");
    quat_s  = strtok(NULL, "\t");
    quat_x  = strtok(NULL, "\t");
    quat_y  = strtok(NULL, "\t");
    quat_z  = strtok(NULL, "\n");

    // convert to floating point
    x = strtod(shape_x, NULL);
    y = strtod(shape_y, NULL);
    z = strtod(shape_z, NULL);
    qs = strtod(quat_s, NULL);
    qx = strtod(quat_x, NULL);
    qy = strtod(quat_y, NULL);
    qz = strtod(quat_z, NULL);

    // construct the quaternion
    quat<Scalar> q(qs, vec3<Scalar>(qx, qy, qz));

    // Making a BAD design decision here, to allocate dynamically and store as generic pointers
    // because something in HOOMD doesn't like for the ShapeConvexPolyhedron object to be copied into an array
    // Talk to Josh and see what he recommends, but for now just proof-of-concept
    
// I'm getting a list of positions/offsets and quaternions for each object read.

    // add the offsets...  need to make a list to store these values, too
//    vec3<Scalar> offsets(x,y,z);
//    v0.x += x; v0.y += y; v0.z += z;
//    v1.x += x; v1.y += y; v1.z += z;
//    v2.x += x; v2.y += y; v2.z += z; 
//    v3.x += x; v3.y += y; v3.z += z;
    
    // create object for each record, then push it onto the list
//    ShapeConvexPolyhedron shape(q, p3dv); 
//    list_of_shapes.push_back(shape);              // Insert a new element at the end
    list_of_shapes.push_back(new ShapeConvexPolyhedron(q, p3dv)); 
//    list_of_positions.push_back(new vec3d<Scalar>(x,y,z));

    
  } // end while !EOF

  return shapes_read;

} // end readPOSFile()

