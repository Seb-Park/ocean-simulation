#include "simulation.h"
#include "graphics/meshloader.h"

#include <unordered_map>
#include <unordered_set>
#include <iostream>

using namespace Eigen;

Simulation::Simulation() {}

int createFaceHash(int a, int b, int c, int n_vertices) {
    int &low = a;
    int &middle = b;
    int &high = c;
    if (low > middle)
    {
        std::swap(low, middle);
    }
    if (middle > high)
    {
        std::swap(middle, high);
    }
    if (low > middle)
    {
        std::swap(low, middle);
    }

    return (n_vertices * n_vertices * low) + (n_vertices * middle) + high;
}

Eigen::Vector3d Simulation::calculateFaceNormal(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c) {
    return (b - a).cross(c - a).normalized();
}

bool Simulation::calculatePointBehindNormal(Eigen::Vector3d a, Eigen::Vector3d b, Eigen::Vector3d c, Eigen::Vector3d p) {
    Eigen::Vector3d a_to_p = p - a; // Calculates a direction from point on plane to point in question
    return a_to_p.dot(calculateFaceNormal(a, b, c)) < 0;
}

void Simulation::init()
{
    // STUDENTS: This code loads up the tetrahedral mesh in 'example-meshes/single-tet.mesh'
    //    (note: your working directory must be set to the root directory of the starter code
    //    repo for this file to load correctly). You'll probably want to instead have this code
    //    load up a tet mesh based on e.g. a file path specified with a command line argument.
    std::vector<Vector3d> vertices;
    std::vector<Vector4i> tets;
    this->mainSystem = System();
    if (MeshLoader::loadTetMesh("./example-meshes/sphere.mesh", vertices, tets)) {
        // STUDENTS: This code computes the surface mesh of the loaded tet mesh, i.e. the faces
        //    of tetrahedra which are on the exterior surface of the object. Right now, this is
        //    hard-coded for the single-tet mesh. You'll need to implement surface mesh extraction
        //    for arbitrary tet meshes. Think about how you can identify which tetrahedron faces
        //    are surface faces...
        std::vector<Vector3i> faces;
//        std::unordered_map<int, std::pair<Eigen::Vector3i*, int>> includedFaces;
//        std::unordered_map<int, int> includedFaces;
        std::unordered_map<int, Eigen::Vector3i> includedFaces;
        // includes the currently parsed faces and their index
        std::vector<Vector4i> faceOrders;

        faceOrders.emplace_back(1, 0, 2, 3); // first three are the included points, last is the excluded point
        faceOrders.emplace_back(2, 0, 3, 1);
        faceOrders.emplace_back(3, 1, 2, 0);
        faceOrders.emplace_back(3, 0, 1, 2);

        for(Vector4i t : tets) {
            for(Vector4i fo : faceOrders) {
                int hash = createFaceHash(t[fo[0]], t[fo[1]], t[fo[2]], vertices.size()); // Finding which vertex indexes the face has and ordering them from least to smallest
                if(includedFaces.contains(hash)) {
                    includedFaces.erase(hash);
                } else {
                    Eigen::Vector3i orderedFace;
                    Vector3d a = vertices.at(t[fo[0]]);
                    Vector3d b = vertices.at(t[fo[1]]);
                    Vector3d c = vertices.at(t[fo[2]]);
                    Vector3d d = vertices.at(t[fo[3]]);
                    if(calculatePointBehindNormal(a, b, c, d)) {
                        orderedFace = Eigen::Vector3i(t[fo[0]], t[fo[1]], t[fo[2]]); // Wind it backwards if the excluded point is behind the face
                    } else {
                        orderedFace = Eigen::Vector3i(t[fo[2]], t[fo[1]], t[fo[0]]); // Wind it backwards if the excluded point is in front of the face
                    }

                    includedFaces.emplace(hash, orderedFace);
                }
            }
        }
        for (const auto & [ key, value ] : includedFaces) {
            faces.push_back(value);
        }
        m_shape.init(vertices, faces, tets);
        mainSystem.initFromVecs(vertices, tets);
    }
    m_shape.setModelMatrix(Affine3f(Eigen::Translation3f(0, 2, 0)));

    initGround();
    initExtra();
}

void Simulation::update(double seconds)
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
    if(this->estimationMode == EULER) {
        mainSystem.updateCalculations();
        mainSystem.updatePositions(seconds);
        mainSystem.resolveCollisions();
    }
    else if (this->estimationMode == MIDPOINT) {
        mainSystem.updateCalculations();
        std::vector<Eigen::Vector3d> originalPositions = mainSystem.getPositions();
        std::vector<Eigen::Vector3d> originalVelocities = mainSystem.getVelocities();
        mainSystem.updatePositions(seconds / 2);
        mainSystem.updateCalculations();
        mainSystem.setVelocities(originalVelocities);
        mainSystem.setPositions(originalPositions);
        mainSystem.updatePositions(seconds);
        mainSystem.resolveCollisions();
    //    m_shape.setVertices(mainSystem.getNodePos());
        m_shape.setVerticesF(mainSystem.getNodePos(), mainSystem.getNodeForces());
    } else if (this->estimationMode == ADAPTIVE) {
        double errorTolerance = 1e-4;

        std::vector<Eigen::Vector3d> originalPositions = mainSystem.getPositions();
        std::vector<Eigen::Vector3d> originalVelocities = mainSystem.getVelocities();

        mainSystem.updateCalculations();
        mainSystem.updatePositions(seconds);
        Eigen::VectorXd state1 = mainSystem.getState();

        mainSystem.setVelocities(originalVelocities);
        mainSystem.setPositions(originalPositions);
        mainSystem.updatePositions(seconds / 2);
        mainSystem.updateCalculations();
        mainSystem.updatePositions(seconds / 2);
        Eigen::VectorXd state2 = mainSystem.getState();

        double newStepsize = seconds * std::sqrt(errorTolerance / (state1 - state2).norm());

        mainSystem.setVelocities(originalVelocities);
        mainSystem.setPositions(originalPositions);
        mainSystem.updateCalculations();
        mainSystem.updatePositions(newStepsize);
        mainSystem.resolveCollisions();
        m_shape.setVerticesF(mainSystem.getNodePos(), mainSystem.getNodeForces());
    }
}

void Simulation::draw(Shader *shader)
{
    m_shape.draw(shader);
    m_ground.draw(shader);
    m_extra.draw(shader);
}

void Simulation::toggleWire()
{
    m_shape.toggleWireframe();
}

void Simulation::toggleForceRender() {
    m_shape.toggleForce();
}

void Simulation::initGround()
{
    std::vector<Vector3d> groundVerts;
    std::vector<Vector3i> groundFaces;
    groundVerts.emplace_back(-5, 0, -5);
    groundVerts.emplace_back(-5, 0, 5);
    groundVerts.emplace_back(5, 0, 5);
    groundVerts.emplace_back(5, 0, -5);
    groundFaces.emplace_back(0, 1, 2);
    groundFaces.emplace_back(0, 2, 3);
    m_ground.init(groundVerts, groundFaces);
}

void Simulation::initExtra()
{
    std::vector<Vector3d> extraVerts;
    std::vector<Vector4i> extraTets;
    Eigen::Vector3d pos = Eigen::Vector3d(0, 0, 0);
    if (MeshLoader::loadTetMesh("./example-meshes/sphere.mesh", extraVerts, extraTets)) {

        std::vector<Vector3i> faces;
        std::unordered_map<int, Eigen::Vector3i> includedFaces;
        std::vector<Vector4i> faceOrders;

        faceOrders.emplace_back(1, 0, 2, 3); // first three are the included points, last is the excluded point
        faceOrders.emplace_back(2, 0, 3, 1);
        faceOrders.emplace_back(3, 1, 2, 0);
        faceOrders.emplace_back(3, 0, 1, 2);

        for(Eigen::Vector3d& ev : extraVerts) {
            ev += pos;
        }

        for(Vector4i t : extraTets) {
            for(Vector4i fo : faceOrders) {
                int hash = createFaceHash(t[fo[0]], t[fo[1]], t[fo[2]], extraVerts.size()); // Finding which vertex indexes the face has and ordering them from least to smallest
                if(includedFaces.contains(hash)) {
                    includedFaces.erase(hash);
                } else {
                    Eigen::Vector3i orderedFace;
                    Vector3d a = extraVerts.at(t[fo[0]]);
                    Vector3d b = extraVerts.at(t[fo[1]]);
                    Vector3d c = extraVerts.at(t[fo[2]]);
                    Vector3d d = extraVerts.at(t[fo[3]]);
                    if(calculatePointBehindNormal(a, b, c, d)) {
                        orderedFace = Eigen::Vector3i(t[fo[0]], t[fo[1]], t[fo[2]]); // Wind it backwards if the excluded point is behind the face
                    } else {
                        orderedFace = Eigen::Vector3i(t[fo[2]], t[fo[1]], t[fo[0]]); // Wind it backwards if the excluded point is in front of the face
                    }

                    includedFaces.emplace(hash, orderedFace);
                }
            }
        }
        for (const auto & [ key, value ] : includedFaces) {
//            std::cout << key << ": " << value << std::endl;
            faces.push_back(value);
        }
        m_extra.init(extraVerts, faces);
        m_extra.setColor(0.9, 0.8, 0.1);
    }
}
