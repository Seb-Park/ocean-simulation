#include "arap.h"
#include "graphics/meshloader.h"

#include <iostream>
#include <set>
#include <map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <QImage>

#include "ocean/ocean.h"

using namespace std;
using namespace Eigen;

ARAP::ARAP() {}

void ARAP::init
	(
	Eigen::Vector3f &coeffMin,
	Eigen::Vector3f &coeffMax
	)
{
	m_num_iterations = 1000;
	m_mesh_path = "meshes/bunny.obj";

	vector<Vector3f> vertices;
	vector<Vector3i> triangles;

	// If this doesn't work for you, remember to change your working directory
//	if (MeshLoader::loadTriMesh(m_mesh_path, vertices, triangles)) {
//		m_shape.init(vertices, triangles);
//	}
	m_ocean.update_ocean();
    vertices = m_ocean.get_vertices();
    triangles = m_ocean.get_faces();
	m_shape.init(vertices, triangles);
    m_foam_shape.init(vertices, triangles);

    m_shape.setColor(0.27f, .803f, .96f);

	// Students, please don't touch this code: get min and max for viewport stuff
	MatrixX3f all_vertices = MatrixX3f(vertices.size(), 3);
	int i = 0;
	for (unsigned long i = 0; i < vertices.size(); ++i) {
		all_vertices.row(i) = vertices[i];
	}
	coeffMin = all_vertices.colwise().minCoeff();
	coeffMax = all_vertices.colwise().maxCoeff();

    minCorner = coeffMin;
    maxCorner = coeffMax;
//
	std::cout << "minCorner: " << minCorner << std::endl;
	std::cout << "maxCorner: " << maxCorner << std::endl;

	GLint texSize; glGetIntegerv(GL_MAX_TEXTURE_BUFFER_SIZE, &texSize);
	std::cout << texSize << std::endl;

//	std::cout << GL_MAX_TEXTURE_BUFFER_SIZE << std::endl;
//	std::cout << GL_PROXY_TEXTURE_2D << std::endl;
//
	minCorner = Vector3f(0.0f, 0.0f, 0.0f);
	maxCorner = Vector3f(1.0f, 0.0f, 1.0f);
//
//	std::cout << "minCorner: " << minCorner << std::endl;
//	std::cout << "maxCorner: " << maxCorner << std::endl;
//
//	minCorner = Vector3f(-1.0f, -1.0f, -1.0f);
//	maxCorner = Vector3f(1.0f, 1.0f, 1.0f);

    initCausticsShape(1000);
}

void ARAP::update(double seconds)
{
    // STUDENTS: This method should contain all the time-stepping logic for your simulation.
    //   Specifically, the code you write here should compute new, updated vertex positions for your
    //   simulation mesh, and it should then call m_shape.setVertices to update the display with those
    //   newly-updated vertices.

    // STUDENTS: As currently written, the program will just continually compute simulation timesteps as long
    //    as the program is running (see View::tick in view.cpp) . You might want to e.g. add a hotkey for pausing
    //    the simulation, and perhaps start the simulation out in a paused state.

    // Note that the "seconds" parameter represents the amount of time that has passed since
    // the last update
	m_ocean.fft_prime(m_time);
    m_ocean.update_ocean();
    m_shape.setVertices_and_Normals(m_ocean.get_vertices(), m_ocean.getNormals());
//      m_shape.setVertices(m_ocean.get_vertices());

//	auto tmp = m_ocean.get_vertices();
//	// print the min and max of the vertices
	// print the min and max of the vertices
//	Vector3f min = Vector3f::Ones() * 1000000;
//	Vector3f max = Vector3f::Ones() * -1000000;
//	for (int i = 0; i < tmp.size(); i++) {
//		min = min.cwiseMin(tmp[i]);
//		max = max.cwiseMax(tmp[i]);
//	}
//	}w
//	std::cout << "min: " << min << std::endl;
//std::cout << "max: " << max << std::endl;

	FoamConstants foam = m_ocean.getFoamConstants();
    m_foam_shape.setFoamInputs(m_ocean.m_vertices, foam.wavelengths, foam.k_vectors, foam.texCoords);


     m_time += m_timestep;
//   std::cout << m_time << std::endl;
}

// Move an anchored vertex, defined by its index, to targetPosition
void ARAP::move
	(
	int vertex,
	Vector3f targetPosition
	)
{
    std::cout << "moving vertex: " << vertex << std::endl;
	// Here are some helpful controls for the application
	//
	// - You start in first-person camera mode
	//   - WASD to move, left-click and drag to rotate
	//   - R and F to move vertically up and down
	//
	// - C to change to orbit camera mode
	//
	// - Right-click (and, optionally, drag) to anchor/unanchor points
	//   - Left-click an anchored point to move it around
	//
	// - Minus and equal keys (click repeatedly) to change the size of the vertices
}

void ARAP::initCausticsShape(int res) {
//    std::vector<Eigen::Vector3f> gridPoints;
//    float step = 2.f / ((float) res);

//    for (int i = 0; i < res; ++i) {
//        for (int j = 0; j < res; ++j) {
//            float x = -1.f + i * step; // calculate x coordinate
//            float y = -1.f + j * step; // calculate y coordinate
//            gridPoints.push_back(Eigen::Vector3f(x, y, 0.f)); // add point to grid
//        }
//    }
//    std::vector<Eigen::Vector3f> verts;
//    float step = 2.f / ((float) res);

//    for (int i = 0; i < res; ++i) {
//        for (int j = 0; j < res; ++j) {
//            float x = -1.f + i * step; // calculate x coordinate
//            float y = -1.f + j * step; // calculate y coordinate
//            Eigen::Vector3f bottomLeft = Eigen::Vector3f(x, y, 0.f);
//            Eigen::Vector3f bottomRight = Eigen::Vector3f(x + step, y, 0.f);
//            Eigen::Vector3f topRight = Eigen::Vector3f(x + step, y + step, 0.f);
//            Eigen::Vector3f topLeft = Eigen::Vector3f(x, y + step, 0.f);
//            verts.push_back(topLeft);
//            verts.push_back(bottomLeft);
//            verts.push_back(bottomRight);
//            verts.push_back(topLeft);
//            verts.push_back(bottomRight);
//            verts.push_back(topRight);
//        }
//    }
//    m_causticsShape.setVertices(verts);
    std::vector<Eigen::Vector3f> verts;
    std::vector<Eigen::Vector3i> faces;
    float size = 2.f;
    float step = size / ((float) res);

    for (int i = 0; i <= res; ++i) {
        for (int j = 0; j <= res; ++j) {
            float x = -(size / 2.f) + i * step; // calculate x coordinate
            float y = -(size / 2.f) + j * step; // calculate y coordinate
            Eigen::Vector3f bottomLeft = Eigen::Vector3f(x, y, 0.f);
            verts.push_back(bottomLeft);
        }
    }

    for (int i = 0; i < res; ++i) {
        for (int j = 0; j < res; ++j) {
            int bottomLeft = i * (res + 1) + j;
            int bottomRight = i * (res + 1) + j + 1;
            int topRight = (i + 1) * (res + 1) + j + 1;
            int topLeft = (i + 1) * (res + 1) + j;
            faces.push_back(Eigen::Vector3i(topLeft, bottomLeft, bottomRight));
            faces.push_back(Eigen::Vector3i(topLeft, bottomRight, topRight));
            faces.push_back(Eigen::Vector3i(bottomRight, bottomLeft, topLeft));
            faces.push_back(Eigen::Vector3i(topRight, bottomRight, topLeft));
        }
    }
    m_causticsShape.init(verts, faces);
    m_causticsShape.setColor(0.27f, .803f, .96f);
}

