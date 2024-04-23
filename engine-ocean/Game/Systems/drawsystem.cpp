#include "drawsystem.h"
#include "Game/Components/DrawComponent.h"
#include "Game/Components/TransformComponent.h"

DrawSystem::DrawSystem(std::map<std::string, std::shared_ptr<GameObject>>& gameobjects):
    m_gameobjects(gameobjects)
{


}

void DrawSystem::update(double deltaTime){}
// store all game objects and get draw components, draw

void DrawSystem::draw(){
    Global::graphics.bindShader("phong");
    //std::cout << "IN DRAW SYSTEM" << std::endl;

    for (auto &pair: m_gameobjects){
        if (pair.second->hasComponent<DrawComponent>()){

            DrawComponent *draw_comp = pair.second->getComponent<DrawComponent>();
            TransformComponent *transform_comp = pair.second->getComponent<TransformComponent>();
            if (transform_comp->hasMultipleMT()){
                if (draw_comp->objHasMaterial()){
                    for (const std::shared_ptr<ModelTransform> &mt : transform_comp->getAllMT()){
                        Global::graphics.drawShape(draw_comp->getShape(), mt, draw_comp->getMaterial());
                    }
                } else {
                    for (const std::shared_ptr<ModelTransform> &mt : transform_comp->getAllMT()){
                        Global::graphics.drawShape(draw_comp->getShape(), mt);
                    }
                }
            } else {
                // case: object is seperated by material
                if (draw_comp->objHasMultipleShapes()){
                    for (auto &shape : draw_comp->getShapesWithMaterials()){
                        if (shape->hasMaterial()){
                            Global::graphics.drawShape(shape, transform_comp->getMT(), shape->getShapeMaterial());
                        } else {
                            Global::graphics.drawShape(shape, transform_comp->getMT());
                        }
                    }

                } else {
                    if (draw_comp->objHasMaterial()){
                        Global::graphics.drawShape(draw_comp->getShape(), transform_comp->getMT(), draw_comp->getMaterial());
                    } else {
                        //std::cout << "draw shape: " << pair.first << std::endl;
                        Global::graphics.drawShape(draw_comp->getShape(), transform_comp->getMT());
                    }
                }
            }
        }
    }

}

void DrawSystem::scrollEvent(double distance){}
void DrawSystem::mousePosEvent(double xpos, double ypos){}
