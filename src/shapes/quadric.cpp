#include "shapes/quadric.h"
#include "core/transform.h"
#include "stats.h"
#include "efloat.h"

namespace pbrt {

Quadric::Quadric(const Transform *ObjectToWorld, const Transform *WorldToObject,
	             bool reverseOrientation, const Matrix4x4 &coeff): Shape(ObjectToWorld, ObjectToWorld, reverseOrientation) {

	// World space coefficient
    Matrix4x4 wcoeff = Matrix4x4::Mul(
        Matrix4x4::Mul(Transpose(*ObjectToWorld).GetInverseMatrix(), coeff),
        WorldToObject->GetMatrix());

    A = wcoeff.m[0][0];
    B = wcoeff.m[1][1];
    C = wcoeff.m[2][2];
    D = wcoeff.m[0][1];
    E = wcoeff.m[1][2];
    F = wcoeff.m[2][0];
    G = wcoeff.m[0][3];
    H = wcoeff.m[1][3];
    I = wcoeff.m[2][3];
    J = wcoeff.m[3][3];
}

bool Quadric::PreIntersect(const Ray &r, EFloat *t0, EFloat *t1) const {
	EFloat dx(r.d.x), dy(r.d.y), dz(r.d.z);
    EFloat ox(r.o.x), oy(r.o.y), oz(r.o.z);

	EFloat a = 2 * (A * dx * dx + B * dy * dy + C * dz * dz + D * dx * dy +
                    E * dy * dz + F * dx * dz);
    EFloat b =
        2 * (A * ox * dx + B * oy * dy + C * oz * dz + D * ox * dy +
             D * oy * dx + E * oy * dz + E * oz * dy + F * ox * dz +
             F * oz * dx + G * dx + H * dy + I * dz);

    EFloat c = A * ox * ox + B * oy * oy + C * oz * oz +
               2 * (D * ox * oy + E * oy * oz + F * ox * oz + G * ox + H * oy +
                    I * oz) +
               J;

	return !Quadratic(a, b, c, t0, t1);
}

}