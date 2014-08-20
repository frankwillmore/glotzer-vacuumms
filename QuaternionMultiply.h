/* QuaternionMultiply.h */

#include <ShapeConvexPolyhedron.h> 
#include <iostream>

using namespace hpmc;
using namespace hpmc::detail;
using namespace std;

inline quat<Scalar> getQuaternion(Scalar angle, vec3<Scalar> axis){
  Scalar st2 = sin(angle/2);
  return quat<Scalar>(cos(angle/2), vec3<Scalar>(axis.x*st2, axis.y*st2, axis.z*st2));
}

inline quat<Scalar> QuaternionMultiply(quat<Scalar> q1, quat<Scalar> q2){
  Scalar s = q1.s*q2.s - q1.v.x*q2.v.x - q1.v.y*q2.v.y - q1.v.z*q2.v.z;
  Scalar x = q1.s*q2.v.x + q1.v.x*q2.s + q1.v.y*q2.v.z - q1.v.z*q2.v.y;
  Scalar y = q1.s*q2.v.y - q1.v.x*q2.v.z + q1.v.y*q2.s + q1.v.z*q2.v.x;
  Scalar z = q1.s*q2.v.z + q1.v.x*q2.v.y - q1.v.y*q2.v.x + q1.v.z*q2.s;
  return quat<Scalar>(s, vec3<Scalar>(x,y,z));
}

inline vec3<Scalar> applyRotation(vec3<Scalar> v, quat<Scalar> q){
  quat<Scalar> q_star(q.s, vec3<Scalar>(-q.v.x, -q.v.y, -q.v.z));
  quat<Scalar> q_out = QuaternionMultiply(q, QuaternionMultiply(quat<Scalar>(0,v), q_star));
  return vec3<Scalar>(q_out.v);
}

