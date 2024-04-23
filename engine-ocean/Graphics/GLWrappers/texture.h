#pragma once

#include "GL/glew.h"
#include <string>

#include "glm/glm.hpp"

class Texture
{
public:
    Texture(int width, int height, GLenum texUnit = GL_TEXTURE0, GLint internalFormat = GL_RGBA, GLenum texTarget = GL_TEXTURE_2D);
    Texture(std::string filePath, GLenum texUnit = GL_TEXTURE0, GLint internalFormat = GL_RGBA, GLenum texTarget = GL_TEXTURE_2D);
    ~Texture();

    void bind();
    void bind(GLenum texUnit);
    void unbind();
    void unbind(GLenum texUnit);
    GLuint getHandle();
    GLuint getTexUnitUint();
    GLenum getTexUnitEnum();
    int getWidth();
    int getHeight();

private:
    GLuint m_handle;
    GLenum m_texUnit;
    GLenum m_texTarget;
    int m_width;
    int m_height;
};
