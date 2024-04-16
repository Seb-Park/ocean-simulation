#include "arap.h"
#include "graphics/meshloader.h"

#include <iostream>
#include <set>
#include <map>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SVD>

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

    vertices = m_ocean.get_vertices();
    triangles = m_ocean.get_faces();
	m_shape.init(vertices, triangles);
    m_shape.setColor(0.27f, .803f, .96f);

	// Students, please don't touch this code: get min and max for viewport stuff
	MatrixX3f all_vertices = MatrixX3f(vertices.size(), 3);
	int i = 0;
	for (unsigned long i = 0; i < vertices.size(); ++i) {
		all_vertices.row(i) = vertices[i];
	}
	coeffMin = all_vertices.colwise().minCoeff();
	coeffMax = all_vertices.colwise().maxCoeff();
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

   //m_ocean.updateVertexAmplitudes(m_time);
    m_ocean.fft_prime(m_time);
   m_shape.setVertices(m_ocean.get_vertices());

  m_time += m_timestep;
  // std::cout << m_time << std::endl;
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

