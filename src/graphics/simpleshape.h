#pragma once

#include <GL/glew.h>
#include <set>
#include <vector>
#include <unordered_set>

#define EIGEN_DONT_VECTORIZE
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
//#include "Eigen/StdVector"
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix2f)
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3f)
//EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::Matrix3i)
#include "Eigen/Dense"

enum SelectMode
{
    None     = 0,
    Anchor   = 1,
    Unanchor = 2
};




struct Edge{
    std::pair<int, int> ij;
    float weight = 0.f;
    std::set<int> opposite_two;
};

struct Cell{
    int center_index = 0; // this is i
    std::set<int> neighbor_js; // this is neighboring j's
    std::set<int> edge_ids; // this is ids if ij edges in edge vector
};


class Shader;

class SimpleShape
{
public:
    SimpleShape();

    void init(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3i> &triangles);
    void setVertices(const std::vector<Eigen::Vector3f> &vertices);

    void setModelMatrix(const Eigen::Affine3f &model);

    void draw(Shader *shader, GLenum mode);
    SelectMode select(Shader *shader, int vertex);
    bool selectWithSpecifiedMode(Shader *shader, int vertex, SelectMode mode);
    int  getClosestVertex(Eigen::Vector3f start, Eigen::Vector3f ray, float threshold);
    bool getAnchorPos(int lastSelected, Eigen::Vector3f& pos, Eigen::Vector3f ray, Eigen::Vector3f start);

    const std::vector<Eigen::Vector3f>& getVertices();
    const std::vector<Eigen::Vector3i>& getFaces();
    const std::unordered_set<int>& getAnchors();

private:
    GLuint m_surfaceVao;
    GLuint m_surfaceVbo;
    GLuint m_surfaceIbo;

    unsigned int m_numSurfaceVertices;
    unsigned int m_verticesSize;
    float m_red;
    float m_blue;
    float m_green;
    float m_alpha;

    std::vector<Eigen::Vector3i> m_faces;
    std::vector<Eigen::Vector3f> m_vertices;
    std::unordered_set<int>      m_anchors;

    Eigen::Matrix4f m_modelMatrix;
    int lastSelected = -1;

    // Helpers

    void selectHelper();
    Eigen::Vector3f getNormal(const Eigen::Vector3i& face);
    void updateMesh(const std::vector<Eigen::Vector3i> &triangles,
                    const std::vector<Eigen::Vector3f> &vertices,
                           std::vector<Eigen::Vector3f>& verts,
                           std::vector<Eigen::Vector3f>& normals,
                           std::vector<Eigen::Vector3f>& colors);
};
