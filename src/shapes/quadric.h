#ifndef PBRT_CORE_QUADRIC_H
#define PBRT_CORE_QUADRIC_H

#include "shape.h"
#include "efloat.h"
namespace pbrt {

class Quadric: public Shape {
  public:
    Quadric(const Transform *ObjectToWorld, const Transform *WorldToObject,
            bool reverseOrientation);
    virtual ~Quadric();
    bool PreIntersect(const Ray &ray, EFloat *t0, EFloat *t1) const;
    virtual void InitCoefficient(Matrix4x4 &m) = 0;
  private:
    EFloat A, B, C, D, E, F, G, H, I, J;
};

}
#endif	// PBRT_CORE_QUADRIC_H