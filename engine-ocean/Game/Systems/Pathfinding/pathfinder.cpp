#include "pathfinder.h"
#include "glm/glm.hpp"
#include "glm/gtx/hash.hpp"

#include <iostream>
#include <set>
#include <vector>

// there can be multiple pathfinders for multiple environemnt meshes
// pathfinding/environment component for an environment mesh?
Pathfinder::Pathfinder(std::vector<glm::vec3> vertices, std::vector<glm::ivec3> triangles):
    m_vertices(vertices),
    m_triangles(triangles)
{
    initializeEdges();
}

// for entire scene, generate data structure that holds nav mesh
void Pathfinder::initializeEdges(){

    std::map<std::pair<int,int>, TriangleEdge> temp_data;

    // indices of vertex pairs
    std::vector<std::pair<int,int>> range = {std::make_pair(0, 1),
                                            std::make_pair(1, 2),
                                            std::make_pair(2, 0)};

    // make edges for all triangles and populate temp_data map
    for (glm::vec3 triangle : m_triangles){
        std::vector<std::pair<int,int>> keys;
        for (auto index_pair : range ){
            // make ordered pair of the index of vertex held by triangle , for all edges
            std::pair<int, int> key = makeOrderedPair(triangle[index_pair.first], triangle[index_pair.second]);
            keys.push_back(key);

            if (temp_data.count(key) != 0){
                temp_data.at(key).count ++;
            } else {
                // make a new entry that populates the midpoint
                TriangleEdge edge = {m_vertices[key.first], m_vertices[key.second]};
                temp_data.insert(std::make_pair(key, edge));
            }
        }

        for (auto key : keys){
            // add the two neighbors for each key
            for (auto neighbor_id : keys){
                if (key != neighbor_id){
                    temp_data.at(key).adjacentEdges.insert(neighbor_id);
                }
            }
        }
    }

    std::set<std::pair<int,int>> interiorEdges;

    // remove exterior edges from adjacent lists
    for (auto &entry : temp_data){
        // make new adjacency list
        std::set<std::pair<int,int>> edited_adj_list;

        for (const auto adjNode : entry.second.adjacentEdges){
            // if adjacent node is an interior node, add to new adj list
            if (temp_data.at(adjNode).count > 1){
                edited_adj_list.insert(adjNode);
            }
        }

        // reassign adj list
        entry.second.adjacentEdges.swap(edited_adj_list);
    }

    // add interior edges
    for (const auto &entry : temp_data){
        if (entry.second.count > 1){
            ANode node;
            node.adjacentNodes = entry.second.adjacentEdges;
            node.pos = entry.second.midpoint;
            m_navdata.insert(std::make_pair(entry.first, node));
        }
    }


}


std::pair<int, int> Pathfinder::makeOrderedPair(int i, int j){
    std::pair<int, int> pair = std::pair(std::min(i,j), std::max(i,j));
    return pair;
}

// traverse navmeshgraph
// A --> start ; B --> destination
std::vector<glm::vec3> Pathfinder::findPath(const glm::vec3 A, const glm::vec3 B){
    // start and end indices are double pair because no existing edge has two of the same vertices
    m_startNodePos = A;
    m_endNodePos = B;

    // initialize empty path
    std::vector<glm::vec3> empty;
    empty.push_back(glm::vec3(0.f));

    // initialize start and ending node
    ANode startNode, endNode;
    startNode.pos = m_startNodePos;
    endNode.pos = m_endNodePos;
    startNode.Gcost = 0;
    m_navdata.insert(std::pair(startNodeID, startNode));
    m_navdata.insert(std::pair(endNodeID, endNode));


    // calculate distance to end point for each item in map
    for (auto &node : m_navdata){
        node.second.Hcost = getDistance(node.second.pos, B);
        updateFCost(node.first);
    }

    // determine if one or both start/end points are in navmesh or not
    std::pair<bool, std::set<std::pair<int, int>>> validStartPoint = findEnclosingTriangle(A);
    std::pair<bool, std::set<std::pair<int, int>>> validEndPoint = findEnclosingTriangle(B);

    // if neither or both points can be found in navmesh, end early
    if (!validStartPoint.first || !validEndPoint.first){
        std::cout << "-- a point is not in navmesh" << std::endl;
        return empty;
    }

    // if both points are on same triangle, then return only the destination
    if (validStartPoint.second == validEndPoint.second){
        std::vector<glm::vec3> destination;
        destination.push_back(B);
        std::cout << "-- SAME TRIANGLE" << std::endl;
        return destination;
    }

    std::set<std::pair<int,int>> open;
    std::set<std::pair<int,int>> closed;

    // initialize open list
    // check if id is a valid interior node
    for (const auto &id : validStartPoint.second){
        if (m_navdata.contains(id)){
            m_navdata[id].previousV = startNodeID;
            m_navdata[id].Gcost = getDistance(m_navdata[id].pos, m_startNodePos);
            updateFCost(id);
            open.insert(id);
        }
    }

    // calculate paths with Astar
    bool reachable = traverseAStar(startNodeID, validEndPoint.second, open, closed);

    if (reachable){
        std::cout << "reachable!" << std::endl;
        return getNavigablePath(startNodeID, endNodeID);
    }

    // otherwise, return glm::vec3(0) so that entity doesnt move
    std::cout << "not reachable" << std::endl;
    return empty;
}

void Pathfinder::updateFCost(const std::pair<int,int> &node){
    m_navdata[node].Fcost = m_navdata[node].Gcost + m_navdata[node].Hcost;

}

bool Pathfinder::traverseAStar(const std::pair<int,int> &currNodeID, const std::set<std::pair<int,int>> &endpointNodeIDs,
                               std::set<std::pair<int,int>> &open, std::set<std::pair<int,int>> &closed){

    glm::vec3 currNodePos = m_navdata.at(currNodeID).pos;

    // base case: if looking at nodes that are way too far from entity start pos, then just end
    if (getDistance(currNodePos, m_startNodePos) > m_maxDistance){
        std::cout << "-- DISTANCE TOO FAR" << std::endl;
        return false;
    }

    float lowestF = INFINITY;
    float lowestH = INFINITY;
    std::pair<int,int> C; // where C is the next node to visit

    for (const std::pair<int,int> &nodeID : open){
        // calculate G, which also populates F
        // --------------- do i also need to check if the previous f cost was less than the updated f cost?
        float new_Gcost = getDistance(m_navdata[nodeID].pos, currNodePos) + m_navdata[currNodeID].Gcost;
        m_navdata[nodeID].Gcost = new_Gcost;
        updateFCost(nodeID);

        if (m_navdata[nodeID].Fcost == lowestF){
            // go by lowest H if F costs are equal
            if (m_navdata[nodeID].Hcost < lowestH){
                C = nodeID;
                lowestF = m_navdata[nodeID].Fcost;
                lowestH = m_navdata[nodeID].Hcost;
            }
        }

        if (m_navdata[nodeID].Fcost < lowestF){
            C = nodeID;
            lowestF = m_navdata[nodeID].Fcost;
            lowestH = m_navdata[nodeID].Hcost;
        }
    }

    // move C from open to closed list
    open.erase(C);
    closed.insert(C);

    //if C is equal to any of the endpoint indices, then path is found
    for (const std::pair<int,int> &id : endpointNodeIDs){
        if (C == id){
            m_navdata[endNodeID].previousV = C;
            return true;
        }
    }

    // for neighbor N of C
    for (const std::pair<int,int> &N: m_navdata[C].adjacentNodes){
        if (closed.contains(N)) continue;

        // calculate potential F cost of C-->N
        float new_Gcost = getDistance(m_navdata[N].pos, m_navdata[C].pos) + m_navdata[C].Gcost;
        float new_Fcost = new_Gcost + m_navdata[N].Hcost;

        if (new_Fcost < m_navdata[N].Fcost || !open.contains(N)){
            // update N's G and F cost
            m_navdata[N].Gcost = new_Gcost;
            updateFCost(N);
            // set previous node of N to be C
            m_navdata[N].previousV = C;

            // add N to open if its not already in it
            open.insert(N);
        }
    }

    // if there are still open nodes, visit the next one (recurse)
    if (!open.empty()){
        traverseAStar(C, endpointNodeIDs, open, closed);
    }

    // otherwise no more open nodes, thus a* has finished
    return true;
}


// distance is the un-squarerooted distance between two points
float Pathfinder::getDistance(glm::vec3 a, glm::vec3 b){
    return pow(a.x-b.x, 2) + pow(a.y-b.y, 2) + pow(a.z-b.z, 2);
}


// referenced from https://blackpawn.com/texts/pointinpoly/default.html
bool Pathfinder::sameSide(glm::vec3 p1, glm::vec3 p2, glm::vec3 a, glm::vec3 b){
    glm::vec3 cp1 = glm::cross(b-a, p1-a);
    glm::vec3 cp2 = glm::cross(b-a, p2-a);

    glm::vec3 a1 = glm::vec3(cp1.x, cp1.y + .001, cp1.z);
    glm::vec3 a2 = glm::vec3(cp1.x, cp1.y - .001, cp1.z);

    glm::vec3 b1 = glm::vec3(cp2.x, cp2.y + .001, cp2.z);
    glm::vec3 b2 = glm::vec3(cp2.x, cp2.y - .001, cp2.z);


    if (glm::dot(a1, b1) >= 0 || glm::dot(a2, b2) >= 0
        || glm::dot(a1, b2) >= 0 || glm::dot(a2, b1) >= 0) return true;
    return false;
}

bool Pathfinder::pointInTriangle(glm::vec3 p, glm::vec3 v1, glm::vec3 v2, glm::vec3 v3){
    if (sameSide(p,v1, v2,v3) && sameSide(p,v2, v1,v3) && sameSide(p,v3, v1,v2)) return true;
    return false;
}

std::pair<bool, std::set<std::pair<int, int>>> Pathfinder::findEnclosingTriangle(glm::vec3 point){

    // for each triangle in navmesh, get its vertices
    for (glm::vec3 vertex_triple : m_triangles){
        // if point is in any navmesh triangle, then return the 3 nodes associated with that tri
        if (pointInTriangle(point, m_vertices[vertex_triple[0]], m_vertices[vertex_triple[1]], m_vertices[vertex_triple[2]])){
            return std::make_pair(true, getEnclosingTriangleEdges(vertex_triple));
        }
    }

    // if no triangles contained the point, return false
    return std::make_pair(false, getEnclosingTriangleEdges(glm::vec3(0.f)));
}

// gets node indices
std::set<std::pair<int, int>> Pathfinder::getEnclosingTriangleEdges(glm::vec3 indices){
    std::set<std::pair<int, int>> ids;
    ids.insert(makeOrderedPair(indices[0], indices[1]));
    ids.insert(makeOrderedPair(indices[1], indices[2]));
    ids.insert(makeOrderedPair(indices[2], indices[0]));

    return ids;
}

std::vector<glm::vec3> Pathfinder::getNavigablePath(const std::pair<int,int> &start, const std::pair<int,int> &end){
    // start at B, and use previousV pointers to get back to A
    std::vector<glm::vec3> path;
    path.push_back(m_navdata[end].pos);

    std::pair<int,int> currV = end;

    while (currV != start){
        // add current point to path
        glm::vec3 pos = m_navdata[m_navdata[currV].previousV].pos;
        path.push_back(pos);

        // reassign variable to where currV points, and enter loop again
        currV = m_navdata[currV].previousV;
    }

    // path where the very last entry is the start point, and the very first entry is the end point
    return path;

}


