#ifndef BOUNDINGELLIPSOID_H
#define BOUNDINGELLIPSOID_H
#include "glm/glm.hpp"
#include "BoundingShape.h"

struct Ellipsoid{
    glm::vec3 R = glm::vec3(0.f); // holds Rx, Ry, Rz radii
    glm::vec3 center_pos;
};

struct UnitSphere{
    glm::vec3 center_pos;
};

class BoundingEllipsoid : public BoundingShape
{
public:
    BoundingEllipsoid();
};

#endif // BOUNDINGELLIPSOID_H
