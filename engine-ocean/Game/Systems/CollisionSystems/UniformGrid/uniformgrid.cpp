#include "uniformgrid.h"
#include "Game/Components/CollisionComponents/BoundingDynamicMesh.h"
#include "Game/Components/CollisionComponents/BoundingTriangle.h"
#include "Game/Components/CollisionComponents/CollisionComponent.h"
#include "Game/Components/CollisionComponents/CylinderCollider.h"
#include <algorithm>

UniformGrid::UniformGrid(std::map<std::string, std::shared_ptr<GameObject>>& dynamic_gameobjects):
    m_dynamic_gameobjects(dynamic_gameobjects)
{
    initializeGrid();
}

void UniformGrid::initializeGrid(){
 // grid square: size it based on the largest entity

    // set default square size to be same dimensions as player
    m_square_dimensions = ceil(m_dynamic_gameobjects.at("player")->getComponent<CollisionComponent>()->getCollisionShape<BoundingDynamicMesh>()->getCylinder().aabbDimensions);


    // find largest shape dimension in order to set grid square dimensions
    for (const auto &go : m_dynamic_gameobjects){
        if (go.second->hasComponent<CollisionComponent>()){
            glm::vec3 shape_dimensions = go.second->getComponent<CollisionComponent>()->getCollisionShape<BoundingDynamicMesh>()->getCylinder().aabbDimensions;
            // find the largest Bounds3f for a shape
            if (glm::all(glm::greaterThan(shape_dimensions, m_square_dimensions))){
                m_square_dimensions = ceil(shape_dimensions);
            }
        }
    }

    // initialize grid based on where entities are initially placed
    for (const auto &go : m_dynamic_gameobjects){
        updateUniformGrid(go);
    }
}

std::string UniformGrid::vecToString(glm::vec3 v){
    std::string square_bound = std::to_string(v.x) + std::to_string(v.y) + std::to_string(v.z);
    return square_bound;
}


void UniformGrid::populateContainingCell(glm::vec3 bound_corner, std::string entity_id){
    glm::vec3 factor = floor(bound_corner/m_square_dimensions);
    glm::vec3 container_min = factor*m_square_dimensions;

    for (auto &square : m_grid_map){
        //if the containing square is already in map, then append entity id to its associate list
        if (square.first == vecToString(container_min)){
            m_grid_map.at(square.first).insert(entity_id);
            return;
        }
    }
    // if not yet in map, add to map
    std::set<std::string> ids;
    ids.insert(entity_id);
    m_grid_map.insert(std::pair(vecToString(container_min), ids));
}

void UniformGrid::updateUniformGrid(const std::pair<std::string, std::shared_ptr<GameObject>> &go){
    if (go.second->hasComponent<CollisionComponent>()){
        Bounds3f shapeBounds = go.second->getComponent<CollisionComponent>()->getCollisionShape<BoundingDynamicMesh>()->getCylinder().bounds;
        // if bounds min is within an interval of square dimensions, then occupy map with that entity
        glm::vec3 shape_dimensions = go.second->getComponent<CollisionComponent>()->getCollisionShape<BoundingDynamicMesh>()->getCylinder().aabbDimensions;


        // check all 6 aabb corners
        float x_val = shapeBounds.min.x;
        for (int x=0; x<2; x++){
            float y_val = shapeBounds.min.y;
            for (int y=0; y<2; y++){
                float z_val = shapeBounds.min.z;
                for (int z=0; z<2; z++){
                    glm::vec3 corner = glm::vec3(x_val,y_val,z_val);
                    populateContainingCell(corner, go.first);
                    z_val += shape_dimensions.z;
                }
                y_val += shape_dimensions.y;
            }
            x_val += shape_dimensions.x;
        }
    }
}

// only for entities that moved:
void UniformGrid::moveEntityGridPosition(const std::pair<std::string, std::shared_ptr<GameObject>> &go){
    // remove previous cells that entity occupied
    for (auto &cell : m_grid_map){
         cell.second.erase(go.first);
    }

    // if there are no more ids associated with that cell, remove it from map
    for (auto cell = m_grid_map.cbegin(), next_i = cell; cell != m_grid_map.cend(); cell = next_i){
      ++next_i;
      if (cell->second.empty()){
        m_grid_map.erase(cell);
      }
    }

    // then update based on new pos
    updateUniformGrid(go);
}

void UniformGrid::updateAllGrid(){
    // initialize grid based on where entities are initially placed
    for (const auto &go : m_dynamic_gameobjects){
            moveEntityGridPosition(go);
    }

    //std::cout << "new ugrid size: " << m_grid_map.size() << std::endl;
}

std::set<std::set<std::string>> UniformGrid::detectPossibleCollisions(){
    // store list of names of objects that are near e/o
    // uses set there are no duplicate inner sets of the same objects
    std::set<std::set<std::string>> candidate_collision_sets;

    for (const auto &cell : m_grid_map){
        if (cell.second.size() > 1){
            candidate_collision_sets.insert(cell.second);
        }
    }

    // then for each inner set, detect collisions between the contained objects
    return candidate_collision_sets;
}
