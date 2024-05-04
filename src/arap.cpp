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

    minCorner = coeffMin;
    maxCorner = coeffMax;


//    m_shape.initGroundPlane("cornell_box_full_lighting.png")
//    QImage ocean_floor_image;
//    GLuint ocean_floor_texture;
//    // Prepare filepath
//    QString ocean_floor_filepath = QString(":/resources/images/kitten.png");

//    // TASK 1: Obtain image from filepath
//    ocean_floor_image = QImage(ocean_floor_filepath);

//    // TASK 2: Format image to fit OpenGL
//    ocean_floor_image = ocean_floor_image.convertToFormat(QImage::Format_RGBA8888).mirrored();

//    // TASK 3: Generate kitten texture
//    glGenTextures(1, &ocean_floor_texture);

//    // TASK 9: Set the active texture slot to texture slot 0
//    glActiveTexture(GL_TEXTURE0);

//    // TASK 4: Bind kitten texture
//    glBindTexture(GL_TEXTURE_2D, ocean_floor_texture);

//    // TASK 5: Load image into kitten texture
//    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ocean_floor_image.width(), ocean_floor_image.height(), 0, GL_RGBA, GL_UNSIGNED_BYTE, ocean_floor_image.bits());

//    // TASK 6: Set min and mag filters' interpolation mode to linear
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

//    // TASK 7: Unbind kitten texture
//    glBindTexture(GL_TEXTURE_2D, 0);

//    // TASK 10: set the texture.frag uniform for our texture
//    glUseProgram(m_texture_shader);
//    glUniform1i(glGetUniformLocation(m_texture_shader, "sampler"), 0);
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
	m_shape.setVertices_and_Normals(m_ocean.get_vertices(), m_ocean.getNormals());
    // m_shape.setVertices(m_ocean.get_vertices());


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

