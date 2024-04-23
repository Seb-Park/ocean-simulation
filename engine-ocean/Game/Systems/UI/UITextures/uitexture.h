#ifndef UITEXTURE_H
#define UITEXTURE_H
#include "Graphics/global.h"
#include <GLFW/glfw3.h>

struct Bounds2f{
    glm::vec2 min;
    glm::vec2 max;
};

enum AspectRatio {
    LAND_FIT,
    LAND_FILL,
    PORTRAIT_FIT,
    PORTRAIT_FILL,
    FIT_SCREEN
};

class UITexture
{
public:
    UITexture();
    virtual void draw() = 0;
    virtual GLuint getTexID() = 0;
    virtual glm::vec2 getPos() = 0;
    virtual glm::vec2 getScale() = 0;
    virtual int getHeight() = 0;
    virtual int getWidth() = 0;
    virtual Bounds2f getBounds() = 0;
    virtual void setWindowPos(int width, int height) = 0;
    virtual float getTextureRatio() = 0;
    virtual float getTextureScaleAspect() = 0;




};

#endif // UITEXTURE_H
