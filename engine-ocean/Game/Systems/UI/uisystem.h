#ifndef UISYSTEM_H
#define UISYSTEM_H
#include "Game/GameWorld.h"
#include "Game/Systems/UI/UITextures/UIButton.h"
#include "Game/Systems/UI/UITextures/uitexture.h"
#include "Game/Systems/UI/uielement.h"
#include "Graphics/global.h"
#include <GLFW/glfw3.h>
#include "Game/Systems/system.h"

struct UIScreen{
    bool isActive = false;
    GLuint screenVAOID;
    std::map<std::string, std::shared_ptr<UIElement>> screenElements;
    std::vector<std::string> elementsDepthOrder;
};


class UISystem : public System
{
public:
    UISystem(std::shared_ptr<Camera> camera,
             std::map<int, Input>& input_map,
             std::set<std::string>& shownScreens);
    ~UISystem();

    void draw() override;
    void update(double deltaTime) override;
    void scrollEvent(double distance) override;
    void mousePosEvent(double xpos, double ypos) override;
    void onWindowResize(int width, int height);
    GLuint makeVAO(std::vector<float> positions);

private:

    int m_quad_numVertices = 4;

    void initializeScreenMap();
    void initializeScreen();
    void initializeStartScreen();
    void initializeInventory();


    void renderScreen();
    glm::mat4 makeTransformationMat(glm::vec2 translation, glm::vec2 scale);
    glm::mat4 getScaleMatrix(glm::vec2 scale);
    bool checkInsideGUI(glm::vec2 mouseClickPos, Bounds2f guiBounds);
    void initializeProfileScreen();


    void makeButtonElement(std::shared_ptr<UIScreen> &screen, std::string elementName, const char* filename, const glm::vec2 pos, const glm::vec2 scale, AspectRatio aspectRatio, bool isCloseButton = false, CornerPosition corner = NONE);
    void makeDisplayElement(std::shared_ptr<UIScreen> &screen, std::string elementName, const char* filename, const glm::vec2 pos, const glm::vec2 scale, AspectRatio aspectRatio = FIT_SCREEN);
    void makeDisplayElement(std::shared_ptr<UIScreen> &screen, std::string elementName, TextureData &tex, const glm::vec2 pos, const glm::vec2 scale, AspectRatio aspectRatio = FIT_SCREEN);



    std::shared_ptr<Camera> m_camera;
    float m_windowWidth = 640;
    float m_windowHeight = 480;
    float m_aspectRatio = 480.f/640.f;

    glm::vec2 m_mousepos = glm::vec2(0.f);
    std::map<int, Input>& m_input_map;

    std::map<std::string, std::shared_ptr<UIScreen>> m_all_screens;

    GLFWcursor* m_pointerCursor;

    std::shared_ptr<ButtonAction> showProfileAction;

    std::set<std::string>& m_shownScreens;
    float m_showStartScreen = true;


    float dim = 1.0f;
    std::vector<float> m_quadPos = {
        -1.0f, 1.0f,
        -1.0f, -1.0f,
        1.0f, 1.0f,
        1.0f, -1.0f
    };

    glm::vec2 drawAspect(AspectRatio aspectType, float textureAspect);



};

#endif // UISYSTEM_H
