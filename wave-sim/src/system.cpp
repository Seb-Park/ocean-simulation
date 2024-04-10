#include "system.h"
#include <iostream>

System::System()
{
    nodes = std::vector<Node*>();
    tets = std::vector<Tet*>();
}

void System::initFromVecs(std::vector<Eigen::Vector3d> v, std::vector<Eigen::Vector4i> t) {
    nodes = std::vector<Node*>();
    tets = std::vector<Tet*>();

    for(Eigen::Vector3d vertPos : v) {
        nodes.push_back( new Node { vertPos } );
    }
    for(Eigen::Vector4i tetVerts : t) {
        std::vector<Node*> definingNodes;
        definingNodes.push_back(nodes.at(tetVerts[0]));
        definingNodes.push_back(nodes.at(tetVerts[1]));
        definingNodes.push_back(nodes.at(tetVerts[2]));
        definingNodes.push_back(nodes.at(tetVerts[3]));
        tets.push_back( new Tet(definingNodes) );
    }
    this->setParams( Params { } );
}

std::vector<Eigen::Vector3d> System::getPositions() {
    std::vector<Eigen::Vector3d> res;
    for(int i = 0 ; i < nodes.size(); ++i) {
        res.push_back(nodes.at(i)->pos);
    }

    return res;
}

void System::setPositions(std::vector<Eigen::Vector3d> positions) {
    assert(positions.size() == nodes.size());
    for(int i = 0 ; i < positions.size(); ++i) {
        nodes.at(i)->pos = positions.at(i);
    }
}

std::vector<Eigen::Vector3d> System::getVelocities() {
    std::vector<Eigen::Vector3d> res;
    for(int i = 0 ; i < nodes.size(); ++i) {
        res.push_back(nodes.at(i)->vel);
    }

    return res;
}

void System::setVelocities(std::vector<Eigen::Vector3d> velocities) {
    assert(velocities.size() == nodes.size());
    for(int i = 0 ; i < velocities.size(); ++i) {
        nodes.at(i)->vel = velocities.at(i);
    }
}

System::Face::Face(Node* a, Node* b, Node* c, Node* p){
    area = (b->pos - a->pos).cross(c->pos - a->pos).norm() / 2;
    normal = (b->pos - a->pos).cross(c->pos - a->pos).normalized();
    if(normal.dot(p->pos - a->pos) >= 0) {
        // If the normal faces towards the back point, flip it and wind nodes other way.
        normal = -normal;
        nodes.push_back(c);
        nodes.push_back(b);
        nodes.push_back(a);
    }
    else {
        nodes.push_back(a);
        nodes.push_back(b);
        nodes.push_back(c);
    }
    opp = p;
}

void System::Face::calculateAreaNorms() {
    Node* a = nodes.at(0);
    Node* b = nodes.at(1);
    Node* c = nodes.at(2);
    this->area = (b->pos - a->pos).cross(c->pos - a->pos).norm() / 2;
    this->normal = (b->pos - a->pos).cross(c->pos - a->pos).normalized();
}

System::Tet::Tet(std::vector<Node*> in_nodes)  {
    nodes = in_nodes;

    Eigen::Vector3d p0 = in_nodes.at(0)->pos;
    Eigen::Vector3d p1 = in_nodes.at(1)->pos;
    Eigen::Vector3d p2 = in_nodes.at(2)->pos;
    Eigen::Vector3d p3 = in_nodes.at(3)->pos;

    Node* a = nodes.at(0);
    Node* b = nodes.at(1);
    Node* c = nodes.at(2);
    Node* d = nodes.at(3); // Origin at d

    // Face 0 defined by 1, 0, 2
    // Face 1 defined by 2, 0, 3
    // Face 2 defined by 3, 1, 2
    // Face 3 defined by 3, 0, 1
    faces.push_back(new Face(b, a, c, d));
    faces.push_back(new Face(c, a, d, b));
    faces.push_back(new Face(d, b, c, a));
    faces.push_back(new Face(d, a, b, c));

    betaMat << (a->pos - d->pos), (b->pos - d->pos), (c->pos - d->pos);

    posMat = betaMat;
    betaMat = betaMat.inverse().eval();
    this->velMat << (a->vel - d->vel), (b->vel - d->vel), (c->vel - d->vel);

    Eigen::Matrix4d massDeterminer;

    massDeterminer << a->pos[0], a->pos[1], a->pos[2], 1,
                      b->pos[0], b->pos[1], b->pos[2], 1,
                      c->pos[0], c->pos[1], c->pos[2], 1,
                      d->pos[0], d->pos[1], d->pos[2], 1;

    mass = this->_rho * (massDeterminer.determinant()) / 6.;

    for(Node* n : this->nodes) {
        n->mass += mass / 4;
    }

    for(Face* f : this->faces) {
        f->calculateAreaNorms();
//            totalFaceArea += f->area;
    }
}

void System::updateCalculations() {
    for(Node* n : nodes) {
        n->force = Eigen::Vector3d::Zero();
        n->force = this->_gravity * n->mass;
    }

    for(Tet* t : tets) {
        t->update();
    }
}

void System::resolveCollisions() {
    double groundLevel = -2;
    double radius = 1;

    for(Node* n : nodes) {
        if(n->pos[1] < groundLevel && std::abs(n->pos[0]) < 5 && std::abs(n->pos[2]) < 5) {
//            n->vel[1] *= -1/(4) * (n->pos[1] + 2);
            n->pos[1] = groundLevel + (groundLevel - n->pos[1]);
            n->vel[0] /= this->groundTraction;
            n->vel[1] = std::abs(n->vel[1]) / this->groundAbsorbance;
        }
        if(n->pos[0] * n->pos[0] + n->pos[2] * n->pos[2] + std::pow(n->pos[1] + 2, 2) < radius * radius) {
            Eigen::Vector3d normal = (n->pos - Eigen::Vector3d(0, -2, 0)).normalized();
//            n->pos = 2 * normal * radius + Eigen::Vector3d(0, -2, 0) - (n->pos - Eigen::Vector3d(0, -2, 0));
            n->pos = normal * radius + Eigen::Vector3d(0, -2, 0);
            Eigen::Vector3d normalComponent = n->vel.dot(normal) * normal;
            Eigen::Vector3d tangentComponent = n->vel - normalComponent;
            n->vel = -(normalComponent / this->groundAbsorbance) - tangentComponent / this->groundTraction;
        }
    }
}

void System::setParams(Params params) {
    for(Tet* t : tets) {
        t->_lambda = params._lambda;
        t->_mu = params._mu;
        t->_rho = params._rho;
        t->_phi = params._phi;
        t->_psi = params._psi;
    }
    this->_gravity = params.gravity;
    this->groundAbsorbance = params._absorbance;
    this->groundTraction = params._traction;
}

void System::updatePositions(double deltaTime) {
//    double totalFaceArea = 0;
    for(Node* n : nodes) {
        n->acc = n->force * (1 / n->mass);
        n->vel += n->acc * deltaTime;
        n->pos += n->vel * deltaTime;
    }
//    this->resolveCollisions();
//    std::cout << totalFaceArea << std::endl;
}

void System::Tet::update() {
//    for(Face* f : this->faces) {
//        f->calculateAreaNorms();
////            totalFaceArea += f->area;
//    }

    Node* a = nodes.at(0);
    Node* b = nodes.at(1);
    Node* c = nodes.at(2);
    Node* d = nodes.at(3); // Origin at d

    this->posMat << (a->pos - d->pos), (b->pos - d->pos), (c->pos - d->pos);
    this->velMat << (a->vel - d->vel), (b->vel - d->vel), (c->vel - d->vel);

    this->FMat = posMat * betaMat;

//    Eigen::Matrix3d test;
//    test << Eigen::Vector3d(1, 2, 3), Eigen::Vector3d(4, 5, 6), Eigen::Vector3d(7, 8, 9);

    if(this->FMat.isApprox(Eigen::Matrix3d::Identity())) {
        this->FMat = Eigen::Matrix3d::Identity();
    }

    this->strain = (FMat.transpose() * FMat) - Eigen::Matrix3d::Identity();
    this->strainRate = ((posMat * betaMat).transpose() * (velMat * betaMat)) + ((velMat * betaMat).transpose() * (posMat * betaMat));

    Eigen::Matrix3d stressFromStrain = (this->_lambda * Eigen::Matrix3d::Identity() * strain.trace()) + (2 * this->_mu * strain);
    Eigen::Matrix3d stressFromStrainRate = (this->_phi * Eigen::Matrix3d::Identity() * strainRate.trace()) + (2 * this->_psi * strainRate);
    this->stress = stressFromStrain + stressFromStrainRate;

//    for(Node* n : this->nodes) {
////        n->force = Eigen::Vector3d(0., -.1, 0.) * n->mass;
////        n->force = Eigen::Vector3d(0, 0, 0);
////        n->acc = Ei
//    }

    for(Face* f : this->faces) {
        Eigen::Vector3d faceForce = this->FMat * this->stress * f->area * f->normal;
        for(Node* n : f->nodes) {
            n->force += faceForce / 3;
        }
//        f->opp->force -= faceForce;
    }
}

std::vector<Eigen::Vector3d> System::getNodePos() {
    std::vector<Eigen::Vector3d> res;
    for(Node* n : nodes) {
        res.push_back(n->pos);
    }
    return res;
}

std::vector<Eigen::Vector3d> System::getNodeForces() {
    std::vector<Eigen::Vector3d> res;
    for(Node* n : nodes) {
        res.push_back(n->force);
    }
    return res;
}

Eigen::VectorXd System::getState() {
    std::vector<double> state;
    for(Node* n : this->nodes) {
        state.push_back(n->pos[0]);
        state.push_back(n->pos[1]);
        state.push_back(n->pos[2]);
        state.push_back(n->vel[0]);
        state.push_back(n->vel[1]);
        state.push_back(n->vel[2]);
    }
    return Eigen::Map<Eigen::VectorXd>(state.data(), state.size());
}
