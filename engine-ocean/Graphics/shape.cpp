#include "shape.h"
#include <Eigen/Dense>
#include "Graphics/global.h"
#include "Graphics/material.h"
#include <iostream>

Shape::Shape(std::shared_ptr<VAO> vao):
    m_vao(vao)
{

}

Shape::Shape(std::shared_ptr<VAO> vao, std::shared_ptr<Material> shape_material):
    m_vao(vao),
    m_shape_material(shape_material)
{
    hasShapeSpecificMaterial = true;

}

Shape::~Shape(){

}

std::shared_ptr<Material> Shape::getShapeMaterial(){
    if (hasShapeSpecificMaterial){
       return m_shape_material;
    }
    std::cout << "this shape does not have material." << std::endl;
}

bool Shape::hasMaterial(){
    return hasShapeSpecificMaterial;
}

void Shape::draw(){
    m_vao->draw();
}

void Shape::updateVAO(const std::vector<Eigen::Vector3f> &vertices, const std::vector<Eigen::Vector3i> &faces){


    m_vao = Global::graphics.makeVAOFromData(vertices, faces, true, false).second;


}
