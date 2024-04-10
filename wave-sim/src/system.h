#pragma once

#include <Eigen/Dense>
#include <vector>

class System
{
public:
    System();

    struct Node;
    struct Tet;

    Eigen::Vector3d _gravity = Eigen::Vector3d(0., -1.0, 0.);
    double groundAbsorbance = 4e4;
    double groundTraction = 4e4;

    std::vector<Node*> nodes;
    std::vector<Tet*> tets;

    std::vector<Eigen::Vector3d> getPositions();
    void setPositions(std::vector<Eigen::Vector3d> positions);

    std::vector<Eigen::Vector3d> getVelocities();
    void setVelocities(std::vector<Eigen::Vector3d> velocities);

    void initFromVecs(std::vector<Eigen::Vector3d> v, std::vector<Eigen::Vector4i> t);

    struct Params {
        const Eigen::Vector3d gravity = Eigen::Vector3d(0, -1, 0);
        const double _lambda = 1e4; //incompressibility for the whole material
        const double _mu = 1e4; //rigidity for the whole material
        const double _rho = 1200.; //density
        const double _phi = 1e1; //coefficients of viscosity
        const double _psi = 1e1;
        const double _traction = 4e4;
        const double _absorbance = 1e1;
    };

    void setParams(Params params);

    struct Node {
        Eigen::Vector3d pos;
        double mass = 0;

        Eigen::Vector3d vel = Eigen::Vector3d(0., 0., 0.);
        Eigen::Vector3d force = Eigen::Vector3d(0., 0., 0.);
        Eigen::Vector3d acc = Eigen::Vector3d(0., 0., 0.);
    };

    struct Face {
        double area;
        Eigen::Vector3d normal;

        std::vector<Node*> nodes;
        Node* opp;

        // P is unused point we use to see if normal is facing right way.
        Face(Node* a, Node* b, Node* c, Node* p);

        void calculateAreaNorms();
    };

    struct Tet {
        double _lambda = 1e4; //incompressibility for the whole material
        double _mu = 1e4; //rigidity for the whole material
        double _rho = 1200.; //density
        double _phi = 1e1; //coefficients of viscosity
        double _psi = 1e1;

        std::vector<Node*> nodes;
//        Eigen::Vector4d faceAreas;
//        std::vector<Eigen::Vector3d> normals;
        std::vector<Face*> faces;
        Eigen::Matrix3d betaMat;
        Eigen::Matrix3d posMat;
        Eigen::Matrix3d velMat;

        Eigen::Matrix3d FMat;
        Eigen::Matrix3d strain;
        Eigen::Matrix3d strainRate;
        Eigen::Matrix3d stress;

        double mass;

        Tet(std::vector<Node*> in_nodes);

        void update();
    };

    void updateCalculations();
    void resolveCollisions();
    void updatePositions(double deltaTime);

    std::vector<Eigen::Vector3d> getNodePos();
    std::vector<Eigen::Vector3d> getNodeForces();

    Eigen::VectorXd getState();
};
