#ifndef UIDISPLAY_H
#define UIDISPLAY_H
#include "uitexture.h"
#include <set>



class UIDisplay : public UITexture
{
public:
    UIDisplay(TextureData tex, glm::vec2 pos, glm::vec2 scale, std::set<std::string>& shownScreens,
              AspectRatio ratio = FIT_SCREEN);
    ~UIDisplay();
    void draw() override;
    GLuint getTexID() override;
    glm::vec2 getPos() override;
    glm::vec2 getScale() override;

    void setWindowPos(int width, int height) override;
    glm::vec2 getWindowPos();

    int getHeight() override;
    int getWidth() override;
    Bounds2f getBounds() override;
    float getTextureRatio() override;
    float getTextureScaleAspect() override;
    void setPos(glm::vec2 pos);
    void setScale(glm::vec2 scale);

    AspectRatio getAspectRatio();


   // glm::vec2 getCornerPos(CornerPosition corner, glm::vec2 elementDimensions);
    void setTransformationMat(glm::vec2 translation, glm::vec2 scale);
    glm::mat4 getTransformationMat();

private:
    int getScreenHeight();
    int getScreenWidth();
    void setTexID(GLuint &newTexID);
    glm::mat4 getScaleMatrix(glm::vec2 scale);



    TextureData m_tex;
    glm::vec2 m_pos;
    glm::vec2 m_scale;
    Bounds2f m_bounds;
    glm::vec2 m_windowPos; // width, height
    int m_windowHeight = 480.f;
    int m_windowWidth = 640.f;
    float m_toScreenScale = 1.f;

    int m_screenImageHeight;
    int m_screenImageWidth;
    float m_tex_aspectRatio = 1.f;
    bool m_isCloseButton = false;
    std::string m_attachedWindow;

    std::set<std::string>& m_shownScreens;
    glm::mat4 m_transformationMat;
    float m_textureAspect = 1.f;

    AspectRatio m_aspectRatio;
};

#endif // UIDISPLAY_H
