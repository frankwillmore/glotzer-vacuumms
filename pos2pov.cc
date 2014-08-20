/* pos2pov.cc */
#include <ShapeConvexPolyhedron.h> 
#include <iostream>

#include "QuaternionMultiply.h"

using namespace hpmc;
using namespace hpmc::detail;
using namespace std;

char *write_tuple(vec3<Scalar> vvv, char *s){
  sprintf(s, "<%g, %g, %g>", vvv.x, vvv.y, vvv.z);
  return s;
}

int main(){

  // reference vectors for tetrahedra
  vec3<Scalar> canonical0(1,1,1);
  vec3<Scalar> canonical1(-1,-1,1);       
  vec3<Scalar> canonical2(-1,1,-1);      
  vec3<Scalar> canonical3(1,-1,-1);

  // read the pos file
  char line[256];
  char *first, *shape_x, *shape_y, *shape_z, *quat_s, *quat_x, *quat_y, *quat_z;
  double x,y,z,qs,qx,qy,qz;

  char *finish = " finish {phong 1} ";
  char *transmit = " transmit 0.5 ";
  char *color = " color rgb<1,0.5,0.25> ";

  // loop over all lines
  while (1) {
    // grab a line; if end of file, then move on.
    fgets(line, 256, stdin);
    if (feof(stdin)) break;

    // get first token
    first = strtok(line, "\t");
    
    // if it doesn't look like a shapeID then toss it. 
    if (first[0] != 'f' || first[1] != 'f') continue;

    // grab the position and quaternion components
    shape_x = strtok(NULL, "\t");
    shape_y = strtok(NULL, "\t");
    shape_z = strtok(NULL, "\t");
    quat_s  = strtok(NULL, "\t");
    quat_x  = strtok(NULL, "\t");
    quat_y  = strtok(NULL, "\t");
    quat_z  = strtok(NULL, "\n");

//printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", shape_x, shape_y, shape_z, quat_s, quat_x, quat_y, quat_z); 

    // convert to floating point
    x = strtod(shape_x, NULL);
    y = strtod(shape_y, NULL);
    z = strtod(shape_z, NULL);
    qs = strtod(quat_s, NULL);
    qx = strtod(quat_x, NULL);
    qy = strtod(quat_y, NULL);
    qz = strtod(quat_z, NULL);

    // get the quaternion representing the rotation
    quat<Scalar> q_rot(qs, vec3<Scalar>(qx, qy, qz));

    // generate the vertices by applying the rotation
    vec3<Scalar> v0 = applyRotation(canonical0, q_rot); 
    vec3<Scalar> v1 = applyRotation(canonical1, q_rot); 
    vec3<Scalar> v2 = applyRotation(canonical2, q_rot); 
    vec3<Scalar> v3 = applyRotation(canonical3, q_rot); 

// add the offsets... 
    v0.x += x; v0.y += y; v0.z += z;
    v1.x += x; v1.y += y; v1.z += z;
    v2.x += x; v2.y += y; v2.z += z; 
    v3.x += x; v3.y += y; v3.z += z;
    
    // four vertices, taken three at a time, will give me triangles which can then be used to generate pov info.
    char s1[256],s2[256],s3[256];

    cout << "triangle {" << write_tuple(v0, s1) << ", " << write_tuple(v1, s2) << ", " << write_tuple(v2, s3) << " texture{ pigment{ color " << color << " " << transmit << " }" << finish << "}}" << endl;
    cout << "triangle {" << write_tuple(v0, s1) << ", " << write_tuple(v1, s2) << ", " << write_tuple(v3, s3) << " texture{ pigment{ color rgb<1,0,0> }}}" << endl;
    cout << "triangle {" << write_tuple(v0, s1) << ", " << write_tuple(v2, s2) << ", " << write_tuple(v3, s3) << " texture{ pigment{ color rgb<1,0,0> }}}" << endl;
    cout << "triangle {" << write_tuple(v1, s1) << ", " << write_tuple(v2, s2) << ", " << write_tuple(v3, s3) << " texture{ pigment{ color rgb<1,0,0> }}}" << endl;

  }
}

triangle {<-12.628, 12.7058, 3.24325>, <-15.4378, 12.3828, 3.22268>, <-13.9821, 12.2576, 0.800844> texture{ pigment{ color rgb<1,0.5,0.25> transmit 0.5 }finish { phong 1}}}
