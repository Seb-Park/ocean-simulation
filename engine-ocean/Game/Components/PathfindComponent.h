#ifndef PATHFINDCOMPONENT_H
#define PATHFINDCOMPONENT_H


#include "Game/Systems/Pathfinding/pathfinder.h"
#include "glm/fwd.hpp"
#include <vector>
#include "Component.h"

class PathfindComponent : public Component
{
public:
    PathfindComponent(std::vector<glm::vec3> vertices, std::vector<glm::ivec3> triangles);
    std::vector<glm::vec3> getPath(const glm::vec3 &A, const glm::vec3 &B);

private:
    std::vector<glm::vec3> m_vertices;
    std::vector<glm::ivec3> m_triangles;
    std::unique_ptr<Pathfinder> m_pathfinder;
};

#endif // PATHFINDCOMPONENT_H
