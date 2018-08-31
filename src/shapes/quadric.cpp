#include "shapes/quadric.h"
#include "core/transform.h"
#include "stats.h"
#include "efloat.h"

namespace pbrt {

Quadric::Quadric(const Transform *ObjectToWorld, const Transform *WorldToObject,
	             bool reverseOrientation, const Matrix4x4 &coeff): Shape(ObjectToWorld, ObjectToWorld, reverseOrientation) {

	// World space coefficient
    Matrix4x4 w1 = Transpose(*ObjectToWorld).GetInverseMatrix();
    Matrix4x4 w2 = Matrix4x4::Mul(w1, coeff);
    Matrix4x4 wcoeff = Matrix4x4::Mul(w2, WorldToObject->GetMatrix());
	 
    A = wcoeff.m[0][0];
    B = wcoeff.m[1][1];
    C = wcoeff.m[2][2];
    D = (wcoeff.m[0][1] + wcoeff.m[1][0]) * 0.5;
    E = (wcoeff.m[1][2] + wcoeff.m[2][1]) * 0.5;
    F = (wcoeff.m[2][0] + wcoeff.m[0][2]) * 0.5;
    G = (wcoeff.m[0][3] + wcoeff.m[3][0]) * 0.5;
    H = (wcoeff.m[1][3] + wcoeff.m[3][1]) * 0.5;
    I = (wcoeff.m[2][3] + wcoeff.m[3][2]) * 0.5;
    J = wcoeff.m[3][3];
}

bool Quadric::PreIntersect(const Ray &r, EFloat *t0, EFloat *t1) const {
	EFloat dx(r.d.x), dy(r.d.y), dz(r.d.z);
    EFloat ox(r.o.x), oy(r.o.y), oz(r.o.z);

	EFloat a = A * dx * dx + B * dy * dy + C * dz * dz + 2 * (D * dx * dy +
                    E * dy * dz + F * dx * dz);
    EFloat b =
        2 * (A * ox * dx + B * oy * dy + C * oz * dz + D * ox * dy +
             D * oy * dx + E * oy * dz + E * oz * dy + F * ox * dz +
             F * oz * dx + G * dx + H * dy + I * dz);

    EFloat c = A * ox * ox + B * oy * oy + C * oz * oz +
               2 * (D * ox * oy + E * oy * oz + F * ox * oz + G * ox + H * oy +
                    I * oz) +
               J;

	return Quadratic(a, b, c, t0, t1);
}

}