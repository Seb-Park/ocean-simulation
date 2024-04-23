#pragma once

#include <GL/glew.h>
#include <vector>

class VBO
{
public:
    VBO(std::vector<float> data);
    ~VBO();

    void bind();
    void unbind();
    GLsizei getLength();

private:
    GLuint m_handle;
    GLsizei m_length;
};