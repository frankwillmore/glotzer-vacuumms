#include <ShapeConvexPolyhedron.h> 
#include <iostream>

#include "QuaternionMultiply.h"

using namespace hpmc;
using namespace hpmc::detail;
using namespace std;

int main(){

  Scalar angle = 3.1415926535/2;

  quat<Scalar> rotate_around_z = getQuaternion(angle, vec3<Scalar>(0,0,1));
  vec3<Scalar> v_in(1, 0, 0);
  vec3<Scalar> v_out = applyRotation(v_in, rotate_around_z); 

  cout << v_out.x << "\t" << v_out.x << "\t" << v_out.y << "\t" << v_out.z << "\n"; 
}

