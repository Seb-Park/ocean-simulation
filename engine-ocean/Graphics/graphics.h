#pragma once

#include "GLWrappers/vao.h"
#include "camera.h"
#include "GLWrappers/shader.h"
#include "material.h"
#include <map>
#include <GL/glew.h>
#include <memory>
#include <vector>
#include "shape.h"
#include "modeltransform.h"
#include "light.h"
#include "textrenderer.h"

#include <Eigen/Dense>

struct TextureData{
    GLuint textureID;
    int width;
    int height;
};

class Graphics
{
public:
    Graphics();
    ~Graphics();

    void initialize();
    void clearScreen(GLbitfield mask);
    void setClearColor(glm::vec3 clearColor);
    TextureData loadTextureFromFile(const char *path);
    TextureData loadTextureFromFile_Repeat(const char *path);

    GLuint loadCubeMap(std::vector<const char*> textureFiles);

    GLuint makeVAO(std::vector<float> positions);
    std::pair<std::vector<glm::vec3>, std::shared_ptr<VAO>> makeVAOFromData(const std::vector<Eigen::Vector3f> &verts,
                                                   const std::vector<Eigen::Vector3i> &faces,
                                                   bool hasUV, bool needsCollisionData = true);
    float generateRandomNumbers(float lb, float ub);

    void setClipPlane(glm::vec4 plane);








    void setCameraData(std::shared_ptr<Camera> camera);

    void addShader(std::string shaderName, std::vector<GLenum> shaderTypes, std::vector<const char*> filepaths);
    void removeShader(std::string shaderName);
    void bindShader(std::string shaderName = "phong");
    GLuint getShaderID(std::string shaderName);


    // Shape-related methods
    std::shared_ptr<Shape> addShape(std::string shapeName, std::vector<float> data, VAOAttrib attribs);
    std::vector<glm::vec3> addShape(std::string shapeName, std::string filepath, bool hasUV = true);
    std::vector<glm::vec3> addShape_manual(std::string shapeName,
                                                     const std::vector<Eigen::Vector3f> &verts,
                                                     const std::vector<Eigen::Vector3i> &faces,
                                                     bool hasUV = true);

    std::vector<glm::vec3> addShape_withMaterial(std::string shapeName, std::string filepath, std::string mtl_file, bool hasUV);

    // Returns pair containing positions followed by triangles indexed into the positions vector
    std::pair<std::vector<glm::vec3>, std::vector<glm::ivec3>> getNavmeshData(std::string filepath);
    void removeShape(std::string shapeName);
    std::shared_ptr<Shape> getShape(std::string shapeName);
    std::vector<std::shared_ptr<Shape>> getShapeGroup(std::string shapeName);
    void drawShape(std::shared_ptr<Shape> myShape, std::shared_ptr<ModelTransform> modelTransform, std::shared_ptr<Material> material = nullptr);
    void drawShape(std::shared_ptr<Shape> myShape, glm::mat4 modelMatrix, std::shared_ptr<Material> material = nullptr);

    // Obj-loading method
    std::vector<float> getObjData(std::string filepath);

    // Material-related methods
    std::shared_ptr<Material> addMaterial(std::string materialName, glm::vec3 color, float shininess = 1.f);
    std::shared_ptr<Material> addMaterial(std::string materialName, float shininess = 1.f);
    std::shared_ptr<Material> addMaterial(std::string materialName, std::string filePath, float shininess = 1.f);
    void removeMaterial(std::string materialName);
    std::shared_ptr<Material> getMaterial(std::string materialName);

    // Text and Font-related methods
    std::shared_ptr<Font> addFont(std::string fontName, std::string filepath);
    void removeFont(std::string fontName);
    std::shared_ptr<Font> getFont(std::string fontName);
    void drawUIText(std::shared_ptr<Font> font, std::string text, glm::vec2 anchorPosition, AnchorPoint anchorPoint, float textBoxWidth, float fontSize, float lineSpacing, glm::vec3 textColor);

    void setGlobalData(glm::vec3 globalCoeffs);

    void setLights(std::vector<std::shared_ptr<Light>> lights);
    void clearLights();

    void initializeGLEW();

    void setWindowSize(glm::ivec2 windowSize);
    glm::ivec2 getWindowSize();

    void setFramebufferSize(glm::ivec2 framebufferSize);
    glm::ivec2 getFramebufferSize();


    // WATER STUFF
    GLuint reflectionFBO, reflectionTexture, reflectionDepthBuffer;
    GLuint refractionFBO, refractionTexture, refractionDepthTexture;

    void bindRefraction();
    void bindReflection();
    GLuint getRefractionTexture();
    GLuint getReflectionTexture();
    void unbindFBO();
    GLuint getReflectionRBO();
    GLuint getRefractionDepth();




    const int REFLECTION_W = 512;
    const int REFLECTION_H = 256;
    const int REFRACTION_W = 512;
    const int REFRACTION_H = 256;

    void setWaterHeight(float height);
    float getWaterHeight();




private:
    void initializeReflectionFBO();
    void initializeRefractionFBO();
    void deleteFBOS();
    void bindFBO(GLuint fbo, int width, int height);
    GLuint createFBO();
    GLuint createTextureAttachment(int width, int height);
    GLuint createDepthTextureAttachment(int width, int height);
    GLuint createDepthBuffer(int width, int height);




    glm::ivec2 m_windowSize;
    glm::ivec2 m_framebufferSize;

    std::shared_ptr<TextRenderer> m_textRenderer;

    std::map<std::string, std::shared_ptr<Shader>> m_shaders;
    std::map<std::string, std::shared_ptr<Shape>> m_shapes;
    std::map<std::string, std::vector<std::string>> m_shape_groups; // stores the shape names kept in m_shapes associated with object
    std::map<std::string, std::shared_ptr<Material>> m_materials;
    std::map<std::string, std::shared_ptr<Font>> m_fonts;

    std::shared_ptr<Shader> m_active_shader;

    float m_water_height = 0.f;


    Eigen::Vector3f getFaceNormal(const std::vector<Eigen::Vector3f> &verts, const Eigen::Vector3i& face);

};
