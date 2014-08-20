// rattle.cc

#include <ShapeConvexPolyhedron.h> 
#include <iostream>
#include "QuaternionMultiply.h"

using namespace std;
using namespace hpmc;
using namespace hpmc::detail;

quat<Scalar> getQuaternion(Scalar angle, vec3<Scalar> axis){
  Scalar st2 = sin(angle/2);
  return quat<Scalar>(cos(angle/2), vec3<Scalar>(axis.x*st2, axis.y*st2, axis.z*st2));
}

int main(){

  // read the .pos file
  readPOSFile();

  // replicate to form supercell
  
  // loop through each shape in main cell
  
  // translate shape in recursive search of space
  
  // give output for each contiguous point found
    int overlaptest(){
      quat<Scalar> q2(1, vec3<Scalar>(0,0,0)); 

      poly3d_verts p3dv;

      // create two identical tri-rectangular tetrahedra
      p3dv.N = 4;
      p3dv.v[0] = vec3<Scalar>(-0.25,-0.25,-0.25);
      p3dv.v[1] = vec3<Scalar>(0.75,-0.25,-0.25);
      p3dv.v[2] = vec3<Scalar>(-0.25,0.75,-0.25);
      p3dv.v[3] = vec3<Scalar>(-0.25,-0.25,0.75);

      ShapeConvexPolyhedron shape2(q2, p3dv); 

      Scalar sqrt3 = sqrt(3.0)/3;

      for (Scalar angle = 0; angle < 2; angle +=0.025) {
	quat<Scalar> q = getQuaternion(3.14159*angle, vec3<Scalar>(sqrt3,sqrt3,sqrt3)); 
	ShapeConvexPolyhedron shape1(q, p3dv); 
	Scalar x;
	for (x=0.0; x<1.1; x+=.01){ // slide to the right until no overlap
	  vec3<Scalar>r_ab(x, 0, 0);
	  unsigned int err=0;
	  if (!test_overlap<ShapeConvexPolyhedron,ShapeConvexPolyhedron>(r_ab, shape1, shape2, err)) break;
	} // end for x
	cout << "for angle = " << angle << ", overlap stops at x = " << x << endl;
      } // end for angle


} // end main

//------------------------------------

int readPOSFile(){

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
  char *color = " rgb<1,0.5,0.25> ";

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
    cout << "triangle {" << write_tuple(v0, s1) << ", " << write_tuple(v1, s2) << ", " << write_tuple(v3, s3) << " texture{ pigment{ color " << color << " " << transmit << " }" << finish << "}}" << endl;
    cout << "triangle {" << write_tuple(v0, s1) << ", " << write_tuple(v2, s2) << ", " << write_tuple(v3, s3) << " texture{ pigment{ color " << color << " " << transmit << " }" << finish << "}}" << endl;
    cout << "triangle {" << write_tuple(v1, s1) << ", " << write_tuple(v2, s2) << ", " << write_tuple(v3, s3) << " texture{ pigment{ color " << color << " " << transmit << " }" << finish << "}}" << endl;

    // create object for each record

    list<int> L;
    L.push_back(0);              // Insert a new element at the end
    L.push_front(0);             // Insert a new element at the beginning
    L.insert(++L.begin(),2);     // Insert "2" before position of first argument
    L.push_back(5);
    L.push_back(6);
    list<int>::iterator i;
    for(i=L.begin(); i != L.end(); ++i) cout << *i << " ";
    cout << endl;

  }
}

