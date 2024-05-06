#pragma once

#include "graphics/shape.h"
//#include "graphics/oceanshape.h"
#include "Eigen/StdList"
#include "Eigen/StdVector"
#include "ocean/ocean.h"
#include "ocean/ocean_alt.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <QList>
#include <QtConcurrent>

class Shader;

class ARAP
{
private:
    Shape m_shape;
    Shape m_foam_shape;


public:
    ARAP();

    void init(Eigen::Vector3f &min, Eigen::Vector3f &max);
    void move(int vertex, Eigen::Vector3f pos);
    void update(double seconds);


    // ================== Students, If You Choose To Modify The Code Below, It's On You

    int getClosestVertex(Eigen::Vector3f start, Eigen::Vector3f ray, float threshold)
    {
        return m_shape.getClosestVertex(start, ray, threshold);
    }

    void draw(Shader *shader, GLenum mode)
    {


        m_shape.draw(shader, mode);
    }

    void drawFoam(Shader *shader, GLenum mode)
    {


        m_foam_shape.draw(shader, mode);
    }

    double getTime() {
        return m_time;
    }

    void initGroundPlane(std::string texturePath, float depth, Shader* shader) {
        m_shape.initGroundPlane(texturePath, depth, shader);
    }

    void initSkyPlane(std::string texturePath, float height, Shader* shader) {
        m_shape.initSkyPlane(texturePath, height, shader);
    }

    SelectMode select(Shader *shader, int vertex)
    {
        return m_shape.select(shader, vertex);
    }

    bool selectWithSpecifiedMode(Shader *shader, int vertex, SelectMode mode)
    {
        return m_shape.selectWithSpecifiedMode(shader, vertex, mode);
    }

    bool getAnchorPos(int lastSelected, Eigen::Vector3f& pos, Eigen::Vector3f ray, Eigen::Vector3f start)
    {
        return m_shape.getAnchorPos(lastSelected, pos, ray, start);
    }

	// for determinig when to recompute
	int num_anchors;

	typedef Eigen::Matrix<float, 3, Eigen::Dynamic> PM; // position matrix
	typedef Eigen::Matrix<float, 3, 3> RM; // rotation matrix
	typedef Eigen::SparseMatrix<float, Eigen::RowMajor> Sparse;

	std::vector<RM> m_rotations;
	Sparse m_edge_weights;
	PM m_b, m_b_fixed;
	std::vector<Eigen::Index> m_vtx_to_free_vtx;
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<float>> m_solver;

	int m_num_iterations;
	const char * m_mesh_path;

    ocean_alt m_ocean;
        double m_time = 0.00;
        double m_timestep = 0.03;

    Eigen::Vector3f minCorner, maxCorner;
};

