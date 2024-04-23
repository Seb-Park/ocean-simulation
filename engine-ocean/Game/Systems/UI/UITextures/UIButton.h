#ifndef UIBUTTON_H
#define UIBUTTON_H
#include "Game/Systems/UI/ButtonAction/buttonaction.h"
#include "Game/Systems/UI/UITextures/uidisplay.h"
#include "Graphics/global.h"
#include <GLFW/glfw3.h>
#include <set>
#include "uitexture.h"

enum CornerPosition {
    TOPLEFT,
    TOPRIGHT,
    BOTTOMLEFT,
    BOTTOMRIGHT,
    NONE
};

class UIButton : public UITexture
{
public:
    UIButton(TextureData tex, glm::vec2 pos, glm::vec2 scale, std::set<std::string>& shownScreens,
             bool isCloseButton = false, CornerPosition corner = NONE, AspectRatio ratio = LAND_FIT);
    ~UIButton();
    void draw() override;
    GLuint getTexID() override;
    glm::vec2 getPos() override;
    glm::vec2 getScale() override;
    float getTextureRatio() override;
    AspectRatio getAspectRatio();



    void setWindowPos(int width, int height) override;
    glm::vec2 getWindowPos();

    int getHeight() override;
    int getWidth() override;
    Bounds2f getBounds() override;
    float getTextureScaleAspect() override;

    void addButtonAction(std::shared_ptr<ButtonAction> &action);
    bool onButtonPress();
    void setWindowToClose(std::string windowID);

    void setParentDisplay(std::string parent);
    std::string getParentDisplay();

    CornerPosition getCornerPos();
    bool hasCornerPos = false;
    void setToCorner(const CornerPosition corner, int width, int height);

    void setTransformationMat(glm::vec2 translation, glm::vec2 scale);
    glm::mat4 getTransformationMat();



private:
    glm::mat4 getScaleMatrix(glm::vec2 scale);
    int getScreenHeight();
    int getScreenWidth();
    void setTexID(GLuint &newTexID);


    TextureData m_tex;
    glm::vec2 m_pos;
    glm::vec2 m_scale;
    Bounds2f m_bounds;
    glm::vec2 m_windowPos; // width, height
    float m_windowHeight = 480.f;
    float m_windowWidth = 640.f;
    float m_toScreenScale = 1.f;

    int m_screenImageHeight;
    int m_screenImageWidth;
    float m_tex_aspectRatio = 1.f;
    bool m_isCloseButton = false;
    std::string m_attachedWindow;

    std::vector<std::shared_ptr<ButtonAction>> m_actions;
    std::set<std::string>& m_shownScreens;

    CornerPosition m_cornerPos;
    std::string m_parentDisplay;

    glm::mat4 m_transformationMat;
    float m_textureAspect = 1.f;

    AspectRatio m_aspectRatio;



};

#endif // UIButton_H
