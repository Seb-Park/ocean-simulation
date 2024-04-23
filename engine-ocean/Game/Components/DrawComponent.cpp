#include "drawcomponent.h"
#include <string>

DrawComponent::DrawComponent(std::shared_ptr<Shape> shape, std::string shape_name)
{
    m_shape = shape;
    m_shape_name = shape_name;
}

DrawComponent::DrawComponent(std::shared_ptr<Shape> shape)
{
    m_shape = shape;
}

DrawComponent::DrawComponent(std::vector<std::shared_ptr<Shape>> shapes)
{
    m_shapes = shapes;
    hasMultipleShapes = true;
}


void DrawComponent::addMaterial(std::string material_name, std::string material_filepath){
    Global::graphics.addMaterial(material_name, material_filepath);
    m_material_name = material_name;
    hasMaterial = true;
}

std::shared_ptr<Shape> DrawComponent::getShape(){
    return m_shape;
}

std::vector<std::shared_ptr<Shape>> DrawComponent::getShapesWithMaterials(){
    if (hasMultipleShapes){
        return m_shapes;
    }
}

std::string DrawComponent::getShapeName(){
    return m_shape_name;
}

std::shared_ptr<Material> DrawComponent::getMaterial(){
    return Global::graphics.getMaterial(m_material_name);
}

bool DrawComponent::objHasMaterial(){
    return hasMaterial;
}

bool DrawComponent::objHasMultipleShapes(){
    return hasMultipleShapes;
}


