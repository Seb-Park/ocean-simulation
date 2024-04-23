#include "texture.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>

#include "../debug.h"

Texture::Texture(int width, int height, GLenum texUnit, GLint internalFormat, GLenum texTarget):
    m_texTarget(texTarget),
    m_texUnit(texUnit)
{
    glGenTextures(1, &m_handle);
    bind();
    glTexParameteri(m_texTarget, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(m_texTarget, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(m_texTarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(m_texTarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexImage2D(m_texTarget, 0, internalFormat, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    m_height = height;
    m_width = width;
    unbind();
}

Texture::Texture(std::string filepath, GLenum texUnit, GLint internalFormat, GLenum texTarget):
    m_texTarget(texTarget),
    m_texUnit(texUnit)
{
    glGenTextures(1, &m_handle);
    bind();
    glTexParameteri(m_texTarget, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(m_texTarget, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(m_texTarget, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(m_texTarget, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    stbi_set_flip_vertically_on_load(1);
    int width, height, numChannels;
    unsigned char* data = stbi_load(filepath.c_str(), &width, &height, &numChannels, 4);
    if (stbi_failure_reason()){
        std::cout << stbi_failure_reason() << std::endl;
    }
    glTexImage2D(m_texTarget, 0, internalFormat, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    m_height = height;
    m_width = width;
    stbi_image_free(data);
    unbind();
}

Texture::~Texture(){
    glDeleteTextures(1, &m_handle);
}

void Texture::bind(){
    glActiveTexture(m_texUnit);
    glBindTexture(m_texTarget, m_handle);
}

void Texture::bind(GLenum texUnit){
    glActiveTexture(texUnit);
    glBindTexture(m_texTarget, m_handle);
}

void Texture::unbind(){
    glActiveTexture(m_texUnit);
    glBindTexture(m_texTarget, 0);
}

void Texture::unbind(GLenum texUnit){
    glActiveTexture(texUnit);
    glBindTexture(m_texTarget, 0);
}

GLuint Texture::getHandle(){
    return m_handle;
}

GLuint Texture::getTexUnitUint(){
    return GLuint(m_texUnit) - GLuint(GL_TEXTURE0);
}

GLenum Texture::getTexUnitEnum(){
    return m_texUnit;
}

int Texture::getWidth(){
    return m_width;
}

int Texture::getHeight(){
    return m_width;
}
