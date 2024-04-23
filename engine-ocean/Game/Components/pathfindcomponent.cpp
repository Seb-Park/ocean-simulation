#include "pathfindcomponent.h"
#include "glm/glm.hpp"
#include <vector>

PathfindComponent::PathfindComponent(std::vector<glm::vec3> vertices, std::vector<glm::ivec3> triangles):
    m_vertices(vertices),
    m_triangles(triangles),
    m_pathfinder(std::make_unique<Pathfinder>(vertices, triangles))
{

}


std::vector<glm::vec3> PathfindComponent::getPath(const glm::vec3 &A, const glm::vec3 &B){
    return m_pathfinder->findPath(A,B);
}
