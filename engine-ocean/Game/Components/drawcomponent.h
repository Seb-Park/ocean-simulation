#ifndef DRAWCOMPONENT_H
#define DRAWCOMPONENT_H


#include "Graphics/global.h"
#include <string>
#include "Component.h"

class DrawComponent : public Component
{
public:
    DrawComponent(std::shared_ptr<Shape> shape);
    DrawComponent(std::shared_ptr<Shape> shape, std::string shape_name);

    // for materials, with multiple shape parts
    DrawComponent(std::vector<std::shared_ptr<Shape>> shapes);


    void draw(const std::shared_ptr<ModelTransform> &entity_mt);
    void addMaterial(std::string material_name, std::string material_filepath);
    std::shared_ptr<Material> getMaterial();
    std::shared_ptr<Shape> getShape();
    std::vector<std::shared_ptr<Shape>> getShapesWithMaterials();
    bool objHasMaterial();
    bool objHasMultipleShapes();

    std::string getShapeName();

private:
    std::shared_ptr<Shape> m_shape;
    std::vector<std::shared_ptr<Shape>> m_shapes;

    std::string m_material_name;
    bool hasMaterial = false;
    std::string m_shape_name = "empty";

    bool hasMultipleShapes = false;

};

#endif // DRAWCOMPONENT_H
