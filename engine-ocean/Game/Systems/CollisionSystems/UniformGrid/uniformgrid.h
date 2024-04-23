#ifndef UNIFORMGRID_H
#define UNIFORMGRID_H


#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Game/GameObjects/GameObject.h"
#include <set>
#include <vector>

struct GridSquare {
    Bounds3f bounds;
    glm::vec3 centroid = (bounds.max - bounds.min)/ 2.f;


//    bool operator<( const GridSquare & s ) const {
//       return glm::all(glm::lessThan(this->bounds.min,s.bounds.min));   // for example
//     }


};

//bool operator<(const GridSquare& l, const GridSquare& r)
//{
//    return glm::all(glm::lessThan(l.bounds.min, r.bounds.min));
//}

//bool operator==(const GridSquare& l, const GridSquare& r)
//{
//    return l.bounds.min == r.bounds.min;
//}

//bool operator<(const glm::vec3& l, const glm::vec3& r)
//{
//    return glm::all(glm::lessThan(l, r));;
//}

//bool operator==(const glm::vec3& l, const glm::vec3& r)
//{
//    return glm::all(glm::equal(l, r));;
//}


class UniformGrid
{
public:
    UniformGrid(std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects);
    std::set<std::set<std::string>> detectPossibleCollisions();
    void updateAllGrid();



private:
    void populateContainingCell(glm::vec3 bound_corner, std::string entity_id);
    void updateUniformGrid(const std::pair<std::string, std::shared_ptr<GameObject>> &go);
    void moveEntityGridPosition(const std::pair<std::string, std::shared_ptr<GameObject>> &go);
    void initializeGrid();
    std::string vecToString(glm::vec3 v);


    std::map<std::string, std::shared_ptr<GameObject>>& m_dynamic_gameobjects;
    // maps square corner to set of strings
    std::map<std::string, std::set<std::string>> m_grid_map;
    glm::vec3 m_square_dimensions;

};

#endif // UNIFORMGRID_H
