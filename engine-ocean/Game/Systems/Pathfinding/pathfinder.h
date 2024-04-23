#ifndef PATHFINDER_H
#define PATHFINDER_H

#include "glm/fwd.hpp"
#include "glm/gtx/hash.hpp"


#include <map>
#include <set>
#include <vector>

struct ANode{
    glm::vec3 pos;
    // A* related:
    float Hcost = 0; // distance to destination
    float Gcost = 0; // distance from start to current node (accumulative G)
    float Fcost = Gcost+Hcost;
    std::pair<int,int> previousV;
    std::set<std::pair<int,int>> adjacentNodes;
};

struct TriangleEdge{
    glm::vec3 midpoint;
    std::set<std::pair<int,int>> adjacentEdges;

    TriangleEdge(glm::vec3 v1, glm::vec3 v2):
        midpoint(.5f*(v2+v1))
    {}
    // num of times it appears in the navgraph
    int count = 1;
};

class Pathfinder
{
public:
    Pathfinder(std::vector<glm::vec3> vertices, std::vector<glm::ivec3> triangles);
    std::vector<glm::vec3> findPath(const glm::vec3 A, const glm::vec3 B);



private:

    void updateFCost(const std::pair<int,int> &node);




    bool traverseAStar(const std::pair<int,int> &currNodeID, const std::set<std::pair<int,int>> &endpointNodeIDs,
                                   std::set<std::pair<int,int>> &open, std::set<std::pair<int,int>> &closed);

    std::set<std::pair<int, int>> getEnclosingTriangleEdges(glm::vec3 indices);

    void initializeEdges();
    float getDistance(glm::vec3 a, glm::vec3 b);
    std::vector<glm::vec3> getNavigablePath(const std::pair<int,int> &start, const std::pair<int,int> &end);
    std::pair<bool, std::set<std::pair<int, int>>> findEnclosingTriangle(glm::vec3 point);
    bool pointInTriangle(glm::vec3 p, glm::vec3 v1, glm::vec3 v2, glm::vec3 v3);
    bool sameSide(glm::vec3 p1, glm::vec3 p2, glm::vec3 a, glm::vec3 b);

    std::set<std::vector<float>> getAdjacentNodes(const std::set<std::pair<int,int>> &id_pairs, std::map<std::pair<int,int>, TriangleEdge> &temp_data);


    std::pair<int, int> makeOrderedPair(int i, int j);

    // maps pos (midpoint of a edge) to node
    std::map<std::pair<int,int>, ANode> m_navdata;
    std::vector<glm::vec3> m_vertices;
    std::vector<glm::ivec3> m_triangles;
    glm::vec3 m_startNodePos;
    glm::vec3 m_endNodePos;

    std::pair<int,int> startNodeID = std::make_pair(0,0);
    std::pair<int,int> endNodeID = std::make_pair(1,1);


    float m_maxDistance = 500; // 20 meters max distance that an ai can pathfind




};

#endif // PATHFINDER_H
