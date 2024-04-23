#pragma once

// Credit to https://learnopengl.com/In-Practice/Text-Rendering

#include "ft2build.h"
#include "freetype/freetype.h"
#include <GL/glew.h>
#include "glm/glm.hpp"
#include <map>
#include <string>

struct Character {
    GLuint textureID;
    glm::ivec2   size;
    glm::ivec2   bearing;
    unsigned int advance;
};

class Font {
public:
    Font(const std::string filepath);
    ~Font();

    Character& getCharacter(char c);

private:
    void createGlyphTexture(unsigned char c, FT_Face &face);
    std::map<char, Character> m_characters;
};