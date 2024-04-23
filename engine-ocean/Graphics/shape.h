#pragma once

#include "GLWrappers/vao.h"
#include "Graphics/material.h"
#include <Eigen/Dense>

class Shape{
public:
    Shape(std::shared_ptr<VAO> vao);
    Shape(std::shared_ptr<VAO> vao, std::shared_ptr<Material> shape_material);

    ~Shape();

    // Draw function
    void draw();
    std::shared_ptr<Material> getShapeMaterial();
    bool hasMaterial();
    void updateVAO(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3i> &faces);




private:
    std::shared_ptr<VAO> m_vao;
    std::shared_ptr<Material> m_shape_material;
    bool hasShapeSpecificMaterial = false;
};
