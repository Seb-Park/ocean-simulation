#include "graphics.h"
#include "shapedata.h"
#include <Eigen/Dense>


#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"


#include "stb_image.h"

Graphics::Graphics():
    m_textRenderer(std::make_shared<TextRenderer>())
{


}

Graphics::~Graphics(){
    deleteFBOS();

}

void Graphics::initializeGLEW(){
    glewExperimental = GL_TRUE;
    std::cout<<"GLEWInit status: "<<glewInit()<<std::endl;
}

void Graphics::initialize(){
    glEnable(GL_DEPTH_TEST);
    //glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    std::cout<<"initialize text renderer"<<std::endl;
    m_textRenderer->initialize();

    std::cout<<"add shapes"<<std::endl;
    addShape("quad", quadVertexBufferData, VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV);
    addShape("cube", cubeVertexBufferData, VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV);
    addShape("sphere", sphereVertexBufferData, VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV);
    addShape("cylinder", cylinderVertexBufferData, VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV);
    addShape("cone", coneVertexBufferData, VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV);

    addShader("phong", {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER}, {"Resources/Shaders/phong.vert", "Resources/Shaders/phong.frag"});
    addShader("text", {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER}, {"Resources/Shaders/text.vert", "Resources/Shaders/text.frag"});
    addShader("grass", {GL_VERTEX_SHADER, GL_GEOMETRY_SHADER, GL_FRAGMENT_SHADER}, {"Resources/Shaders/grass.vert", "Resources/Shaders/grassgeom.geom", "Resources/Shaders/grass.frag"});
    addShader("skybox", {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER}, {"Resources/Shaders/skybox.vert", "Resources/Shaders/skybox.frag"});
    addShader("ui", {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER}, {"Resources/Shaders/UIShader.vert", "Resources/Shaders/UIShader.frag"});
    addShader("inventory", {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER}, {"Resources/Shaders/UIShader.vert", "Resources/Shaders/InventoryShader.frag"});
    addShader("water", {GL_VERTEX_SHADER, GL_FRAGMENT_SHADER}, {"Resources/Shaders/water.vert", "Resources/Shaders/water.frag"});





    bindShader("phong");

    addMaterial("default", glm::vec3(1));

    std::cout<<"add font"<<std::endl;
    addFont("opensans", "Resources/Fonts/OpenSans.ttf");
    getFont("opensans");

    // specify texture for grass shader
    bindShader("grass");
    glUniform1i(glGetUniformLocation(getShaderID("grass"), "grass_texture"),0);
    glUniform1i(glGetUniformLocation(getShaderID("grass"), "wind_texture"),1);
    glUniform1i(glGetUniformLocation(getShaderID("grass"), "meadow_texture"),7);

    //glUniform1f(glGetUniformLocation(getShaderID("grass"), "numRows"), 2);


    // specify texture for skybox
    bindShader("skybox");
    glUniform1i(glGetUniformLocation(getShaderID("skybox"), "cubeMap"), 3);

    bindShader("ui");
    glUniform1i(glGetUniformLocation(getShaderID("ui"), "guiTexture0"), 5);


    bindShader("inventory");
    glUniform1i(glGetUniformLocation(getShaderID("inventory"), "inventoryTexture0"), 6);

    initializeReflectionFBO();
    initializeRefractionFBO();

    bindShader("water");
    glUniform1i(glGetUniformLocation(getShaderID("water"), "reflect_texture"), 8);
    glUniform1i(glGetUniformLocation(getShaderID("water"), "refract_texture"), 9);
    glUniform1i(glGetUniformLocation(getShaderID("water"), "du_dv_map"), 10);


    glm::vec3 sc = glm::vec3(.77f, .85f, .99f);
    bindShader("skybox");
    glUniform3f(glGetUniformLocation(getShaderID("skybox"), "skyColor"), sc.x, sc.y, sc.z);
    bindShader("water");
    glUniform3f(glGetUniformLocation(getShaderID("water"), "skyColor"), sc.x, sc.y, sc.z);
    bindShader("grass");
    glUniform3f(glGetUniformLocation(getShaderID("grass"), "skyColor"), sc.x, sc.y, sc.z);
    bindShader("phong");
    glUniform3f(glGetUniformLocation(getShaderID("phong"), "skyColor"), sc.x, sc.y, sc.z);





}

void Graphics::setClipPlane(glm::vec4 plane){
    bindShader("phong");
    glUniform4f(glGetUniformLocation(getShaderID("phong"), "plane"), plane.x, plane.y, plane.z, plane.w);
    bindShader("grass");
    glUniform4f(glGetUniformLocation(getShaderID("grass"), "plane"), plane.x, plane.y, plane.z, plane.w);
    bindShader("water");
    glUniform4f(glGetUniformLocation(getShaderID("water"), "plane"), plane.x, plane.y, plane.z, plane.w);
}

float Graphics::generateRandomNumbers(float lb, float ub){
    float random = lb + static_cast <float>(rand()) / (static_cast <float>(RAND_MAX/(ub-lb)));

    return random;
}


GLuint Graphics::loadCubeMap(std::vector<const char*> textureFiles){
    // create empty texture
    GLuint textureID;
    glGenTextures(1, &textureID);

    //glActiveTexture(GL_TEXTURE3);
    glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    GLuint target = GL_TEXTURE_CUBE_MAP_POSITIVE_X;
    for (int i=0; i<6; i++){
        std::string filename = std::string(textureFiles[i]);//directory + '/' + filename;
        int width, height, nrChannels;
        unsigned char* data = stbi_load(filename.c_str(), &width, &height, &nrChannels, 0);

        if (data){
            stbi_set_flip_vertically_on_load(false);
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
            stbi_image_free(data);

        }    else {
            std::cout << "Texture failed to load at path: " << textureFiles[i] << std::endl;
            stbi_image_free(data);
        }
    }

    return textureID;
}


TextureData Graphics::loadTextureFromFile(const char *path)
{
    std::string filename = std::string(path);

    GLuint textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    stbi_set_flip_vertically_on_load(true);
    unsigned char *data = stbi_load(filename.c_str(), &width, &height, &nrComponents, 0);
    stbi_set_flip_vertically_on_load(false);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    TextureData newtex;
    newtex.textureID = textureID;
    newtex.height = height;
    newtex.width = width;
    return newtex;
}

TextureData Graphics::loadTextureFromFile_Repeat(const char *path)
{
    std::string filename = std::string(path);

    GLuint textureID;
    glGenTextures(1, &textureID);

    int width, height, nrComponents;
    stbi_set_flip_vertically_on_load(true);
    unsigned char *data = stbi_load(filename.c_str(), &width, &height, &nrComponents, 0);
    stbi_set_flip_vertically_on_load(false);
    if (data)
    {
        GLenum format;
        if (nrComponents == 1)
            format = GL_RED;
        else if (nrComponents == 3)
            format = GL_RGB;
        else if (nrComponents == 4)
            format = GL_RGBA;

        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, GL_UNSIGNED_BYTE, data);
        glGenerateMipmap(GL_TEXTURE_2D);

        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        stbi_image_free(data);
    }
    else
    {
        std::cout << "Texture failed to load at path: " << path << std::endl;
        stbi_image_free(data);
    }

    TextureData newtex;
    newtex.textureID = textureID;
    newtex.height = height;
    newtex.width = width;
    return newtex;
}



GLuint Graphics::makeVAO(std::vector<float> positions){
    int numVertices = positions.size()/2;
    GLuint vao, vbo;

    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);

    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, positions.size()*sizeof(float), positions.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(float), (void*)0);

    // unbind
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    return vao;
}

void Graphics::clearScreen(GLbitfield mask){
    glClear(mask);
}

void Graphics::setClearColor(glm::vec3 clearColor){
    glClearColor(clearColor.r, clearColor.g, clearColor.b, 1);
}

void Graphics::setCameraData(std::shared_ptr<Camera> camera){
    m_active_shader->setCamera(camera);
}

void Graphics::addShader(std::string shaderName, std::vector<GLenum> shaderTypes, std::vector<const char*> filepaths){
    m_shaders.insert({shaderName, std::make_shared<Shader>(shaderTypes, filepaths)});
}

void Graphics::removeShader(std::string shaderName){
    if(shaderName != "phong" && shaderName != "text"){
        m_shaders.erase(shaderName);
    }
}

void Graphics::bindShader(std::string shaderName){
    m_shaders.at(shaderName)->bind();
    m_active_shader = m_shaders.at(shaderName);
}

GLuint Graphics::getShaderID(std::string shaderName){
    return m_shaders.at(shaderName)->getHandle();

}

std::shared_ptr<Shape> Graphics::addShape(std::string shapeName, std::vector<float> data, VAOAttrib attribs){
    m_shapes.insert({shapeName, std::make_shared<Shape>(std::make_shared<VAO>(std::make_shared<VBO>(data), attribs))});
    return m_shapes.at(shapeName);
}
std::vector<glm::vec3> Graphics::addShape(std::string shapeName, std::string filepath, bool hasUV){
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str())) {
        throw std::runtime_error(warn + err);
    }

    int numTriangles = 0;
    for(size_t s = 0; s < shapes.size(); s++){
        numTriangles += shapes[s].mesh.num_face_vertices.size();
    }

    std::vector<float> drawData;
    drawData.resize(numTriangles * 3 * 8);
    std::vector<glm::vec3> collisionData;
    collisionData.resize(numTriangles * 3);

    int i = 0;
    int j = 0;
    for(size_t s = 0; s < shapes.size(); s++) {
        size_t index_offset = 0;
        for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            unsigned int fv = shapes[s].mesh.num_face_vertices[f];

            for(size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

                // Add position
                drawData[i] = attrib.vertices[3*idx.vertex_index];
                drawData[i + 1] = attrib.vertices[3*idx.vertex_index + 1];
                drawData[i + 2] = attrib.vertices[3*idx.vertex_index + 2];
                // Add normal
                drawData[i + 3] = attrib.normals[3*idx.normal_index];
                drawData[i + 4] = attrib.normals[3*idx.normal_index + 1];
                drawData[i + 5] = attrib.normals[3*idx.normal_index + 2];
                // Add uv
                drawData[i + 6] = 0;
                drawData[i + 7] = 0;
                if(hasUV){
                    drawData[i + 6] = attrib.texcoords[2*idx.texcoord_index];
                    drawData[i + 7] = attrib.texcoords[2*idx.texcoord_index + 1];
                }

                // Add collision position data
                collisionData[j] = glm::vec3(attrib.vertices[3*idx.vertex_index], attrib.vertices[3*idx.vertex_index + 1], attrib.vertices[3*idx.vertex_index + 2]);
                i += 8;
                j += 1;
            }

            index_offset += fv;
        }
    }

    m_shapes.insert({shapeName, std::make_shared<Shape>(std::make_shared<VAO>(std::make_shared<VBO>(drawData), VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV))});

    std::cout << "coliision data size normal: " << collisionData.size() << std::endl;
    return collisionData;
}

std::vector<glm::vec3> Graphics::addShape_manual(std::string shapeName, 
                                                 const std::vector<Eigen::Vector3f> &verts, 
                                                 const std::vector<Eigen::Vector3i> &faces,
                                                 bool hasUV){

    std::pair<std::vector<glm::vec3>, std::shared_ptr<VAO>> info = makeVAOFromData(verts, faces, hasUV);

    m_shapes.insert({shapeName, std::make_shared<Shape>(info.second)});

    std::cout << "ADD SHAPE MANUAL" << std::endl;

    return info.first;
}

std::pair<std::vector<glm::vec3>, std::shared_ptr<VAO>> Graphics::makeVAOFromData(const std::vector<Eigen::Vector3f> &verts,
                                               const std::vector<Eigen::Vector3i> &faces,
                                               bool hasUV, bool needsCollisionData){



  std::vector<float> drawData;
  drawData.resize(faces.size() * 3 * 8);
  std::vector<glm::vec3> collisionData;
  collisionData.resize(faces.size() * 3);

  int i = 0;
  int j = 0;

  float uv = 0;

  size_t index_offset = 0;
  for(size_t f = 0; f < faces.size(); f++) {
      Eigen::Vector3i face = faces[f];
      Eigen::Vector3f faceNormal = getFaceNormal(verts, face);


      for(int v = 0; v < 3; v++) {
          int idx = face[v];
          Eigen::Vector3f vertexPos = verts[idx];

          // Add position
          drawData[i] = vertexPos[0];
          drawData[i + 1] = vertexPos[1];
          drawData[i + 2] = vertexPos[2];
          // Add normal
          drawData[i + 3] = 0.f;
          drawData[i + 4] = 1.f;
          drawData[i + 5] = 0.f;
          // Add uv
          drawData[i + 6] = 0;
          drawData[i + 7] = 0;

          if(hasUV){
              drawData[i + 6] = uv;
              drawData[i + 7] = uv;

              uv += .01;
          }

          if (needsCollisionData){
              // Add collision position data
              collisionData[j] = glm::vec3(vertexPos[0], vertexPos[1], vertexPos[2]);
          }
          i += 8;
          j += 1;
      }
  }


  std::shared_ptr<VAO> vao = std::make_shared<VAO>(std::make_shared<VBO>(drawData), VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV);
  return std::make_pair(collisionData, vao);


}

Eigen::Vector3f Graphics::getFaceNormal(const std::vector<Eigen::Vector3f> &verts, const Eigen::Vector3i& face)
{
    Eigen::Vector3f v1 = verts[face[0]];
    Eigen::Vector3f v2 = verts[face[1]];
    Eigen::Vector3f v3 = verts[face[2]];
    Eigen::Vector3f e1 = v2 - v1;
    Eigen::Vector3f e2 = v3 - v1;
    Eigen::Vector3f n = e1.cross(e2);
    return n.normalized();
}



// MATERIALSS
std::vector<glm::vec3> Graphics::addShape_withMaterial(std::string shapeName, std::string filepath, std::string mtl_file, bool hasUV){
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str(), mtl_file.c_str())) {
        throw std::runtime_error(warn + err);
    }

    int numTriangles = 0;
    for(size_t s = 0; s < shapes.size(); s++){
        numTriangles += shapes[s].mesh.num_face_vertices.size();
    }

    std::cout << "materials size: " << materials.size() << std::endl;
    for (auto m : materials){
        std::cout << "material: " << m.name << std::endl;
    }


    std::vector<std::string> shape_names;
    std::vector<glm::vec3> collisionData;
    std::cout << "---- NUM Triangles: " << numTriangles * 3 << std::endl;
    collisionData.resize(numTriangles * 3);



    int j = 0;
    for(size_t s = 0; s < shapes.size(); s++) {
        int i = 0;
        bool hasMaterial = true;
        std::shared_ptr<Material> shape_material;
        std::vector<float> drawData;
        drawData.resize(shapes[s].mesh.num_face_vertices.size() * 3 * 8);

        size_t index_offset = 0;

        int shape_material_id = shapes[s].mesh.material_ids[0];

        if (shape_material_id < 0) hasMaterial = false;
        if (hasMaterial){
            std::string shape_material_name = materials[shape_material_id].name;
            //std::cout << "load: " << shape_material_name << std::endl;
            std::string shape_material_path = "/Users/jesswan/Desktop/cs1950u/cs1950u-jjesswan/Resources/Images/" + shape_material_name + ".png";
            shape_material = addMaterial(shape_material_name, shape_material_path, 1.f);
        }

//        std::cout << "shape name:" << shapes[s].name << std::endl;
//        std::cout << "mesh material id sample:" << shapes[s].mesh.material_ids[0] << std::endl;
//        std::cout << "mesh material id sample2:" << shapes[s].mesh.material_ids[5] << std::endl;

//        std::cout << "mesh material id size:" << shapes[s].mesh.material_ids.size() << std::endl;
//        std::cout << "num faces:" << shapes[s].mesh.num_face_vertices.size() << std::endl;

        for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            unsigned int fv = shapes[s].mesh.num_face_vertices[f];

            for(size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

                // Add position
                drawData[i] = attrib.vertices[3*idx.vertex_index];
                drawData[i + 1] = attrib.vertices[3*idx.vertex_index + 1];
                drawData[i + 2] = attrib.vertices[3*idx.vertex_index + 2];
                // Add normal
                drawData[i + 3] = attrib.normals[3*idx.normal_index];
                drawData[i + 4] = attrib.normals[3*idx.normal_index + 1];
                drawData[i + 5] = attrib.normals[3*idx.normal_index + 2];
                // Add uv
                drawData[i + 6] = 0;
                drawData[i + 7] = 0;
                if(hasUV){
                    drawData[i + 6] = attrib.texcoords[2*idx.texcoord_index];
                    drawData[i + 7] = attrib.texcoords[2*idx.texcoord_index + 1];
                }

                // Add collision position data
                collisionData[j] = glm::vec3(attrib.vertices[3*idx.vertex_index], attrib.vertices[3*idx.vertex_index + 1], attrib.vertices[3*idx.vertex_index + 2]);
                i += 8;
                j += 1;
            }

            index_offset += fv;
        }

        if (hasMaterial){
            m_shapes.insert({shapeName + std::to_string(s), std::make_shared<Shape>(std::make_shared<VAO>(std::make_shared<VBO>(drawData), VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV), shape_material)});
        } else {
            m_shapes.insert({shapeName + std::to_string(s), std::make_shared<Shape>(std::make_shared<VAO>(std::make_shared<VBO>(drawData), VAOAttrib::POS | VAOAttrib::NORM | VAOAttrib::UV))});
        }

        shape_names.push_back(shapeName + std::to_string(s));
    }

    std::cout << "DONE LOADING MTL" << std::endl;
    m_shape_groups.insert({shapeName, shape_names});
    std::cout << "collision data size: " << collisionData.size() << std::endl;
    return collisionData;
}

std::pair<std::vector<glm::vec3>, std::vector<glm::ivec3>> Graphics::getNavmeshData(std::string filepath){
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str())) {
        throw std::runtime_error(warn + err);
    }

    int numTriangles = 0;
    for(size_t s = 0; s < shapes.size(); s++){
        numTriangles += shapes[s].mesh.num_face_vertices.size();
    }

    std::vector<glm::vec3> positions;
    std::vector<glm::ivec3> faces;

    for(int i = 0; i<attrib.vertices.size(); i+=3){
        positions.push_back(glm::vec3(attrib.vertices[i], attrib.vertices[i+1], attrib.vertices[i+2]));
    }


    int i = 0;
    int j = 0;
    for(size_t s = 0; s < shapes.size(); s++) {
        for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            tinyobj::index_t idx1 = shapes[s].mesh.indices[3* f];
            tinyobj::index_t idx2 = shapes[s].mesh.indices[3* f + 1];
            tinyobj::index_t idx3 = shapes[s].mesh.indices[3* f + 2];
            faces.push_back(glm::ivec3(idx1.vertex_index, idx2.vertex_index, idx3.vertex_index));
        }
    }

    return std::make_pair(positions, faces);
}

void Graphics::removeShape(std::string shapeName){
    m_shapes.erase(shapeName);
}

std::shared_ptr<Shape> Graphics::getShape(std::string shapeName){
    return m_shapes.at(shapeName);
}

std::vector<std::shared_ptr<Shape>> Graphics::getShapeGroup(std::string shapeName){
    std::vector<std::string> shape_ids = m_shape_groups.at(shapeName);
    std::vector<std::shared_ptr<Shape>> all_shapes;

    for (std::string id : shape_ids){
        all_shapes.push_back(m_shapes.at(id));
    }

    return all_shapes;
}

void Graphics::drawShape(std::shared_ptr<Shape> myShape, std::shared_ptr<ModelTransform> modelTransform, std::shared_ptr<Material> material){
    if(material == nullptr){
        m_active_shader->setMaterial(getMaterial("default"));
    }
    else{
        m_active_shader->setMaterial(material);
    }
    m_active_shader->setModelTransform(modelTransform);
    myShape->draw();
}

void Graphics::drawShape(std::shared_ptr<Shape> myShape, glm::mat4 modelMatrix, std::shared_ptr<Material> material){
    if(material == nullptr){
        m_active_shader->setMaterial(getMaterial("default"));
    }
    else{
        m_active_shader->setMaterial(material);
    }
    m_active_shader->setModelTransform(modelMatrix);
    myShape->draw();
}

std::vector<float> Graphics::getObjData(std::string filepath){
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;
    std::string warn, err;

    if (!tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filepath.c_str())) {
        throw std::runtime_error(warn + err);
    }

    std::vector<glm::vec3> faces;
    std::vector<glm::vec3> vertices;

    for(size_t s = 0; s < shapes.size(); s++) {
        size_t index_offset = 0;
        for(size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
            unsigned int fv = shapes[s].mesh.num_face_vertices[f];

            glm::vec3 face;
            for(size_t v = 0; v < fv; v++) {
                tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];

                face[v] = idx.vertex_index;

            }
            faces.push_back(face);

            index_offset += fv;
        }
    }
    for (size_t i = 0; i < attrib.vertices.size(); i += 3) {
        vertices.emplace_back(attrib.vertices[i], attrib.vertices[i + 1], attrib.vertices[i + 2]);
    }

    std::vector<float> data;
    data.resize(18 * faces.size());
    for(int i = 0; i<faces.size(); i++){
        glm::vec3 v0 = vertices[static_cast<int>(faces[i].x)];
        glm::vec3 v1 = vertices[static_cast<int>(faces[i].y)];
        glm::vec3 v2 = vertices[static_cast<int>(faces[i].z)];
        glm::vec3 normal = glm::normalize(glm::cross(v1 - v0, v2 - v0));

        data[18*i] = v0.x;
        data[18*i+1] = v0.y;
        data[18*i+2] = v0.z;
        data[18*i+3] = normal.x;
        data[18*i+4] = normal.y;
        data[18*i+5] = normal.z;
        data[18*i+6] = v1.x;
        data[18*i+7] = v1.y;
        data[18*i+8] = v1.z;
        data[18*i+9] = normal.x;
        data[18*i+10] = normal.y;
        data[18*i+11] = normal.z;
        data[18*i+12] = v2.x;
        data[18*i+13] = v2.y;
        data[18*i+14] = v2.z;
        data[18*i+15] = normal.x;
        data[18*i+16] = normal.y;
        data[18*i+17] = normal.z;
    }

    return data;
}

std::shared_ptr<Material> Graphics::addMaterial(std::string materialName, glm::vec3 color, float shininess){
    std::shared_ptr<Material> newMaterial = std::make_shared<Material>(color, shininess);
    m_materials.insert({materialName, newMaterial});
    return m_materials.at(materialName);
}

std::shared_ptr<Material> Graphics::addMaterial(std::string materialName, float shininess){
    std::shared_ptr<Material> newMaterial = std::make_shared<Material>(shininess);
    m_materials.insert({materialName, newMaterial});
    return m_materials.at(materialName);
}

std::shared_ptr<Material> Graphics::addMaterial(std::string materialName, std::string filePath, float shininess){
    std::shared_ptr<Texture> newTexture = std::make_shared<Texture>(filePath);
    std::shared_ptr<Material> newMaterial = std::make_shared<Material>(newTexture, shininess);
    m_materials.insert({materialName, newMaterial});
    return m_materials.at(materialName);
}

void Graphics::removeMaterial(std::string materialName){
    if(materialName != "default"){
        m_materials.erase(materialName);
    }
}

std::shared_ptr<Material> Graphics::getMaterial(std::string materialName){
    return m_materials.at(materialName);
}

std::shared_ptr<Font> Graphics::addFont(std::string fontName, std::string filepath){
    std::shared_ptr<Font> newFont = std::make_shared<Font>(filepath);
    m_fonts.insert({fontName, newFont});
    return m_fonts.at(fontName);
}

void Graphics::removeFont(std::string fontName){
    m_fonts.erase(fontName);
}

std::shared_ptr<Font> Graphics::getFont(std::string fontName){
    return m_fonts.at(fontName);
}

void Graphics::drawUIText(std::shared_ptr<Font> font, std::string text, glm::vec2 anchorPosition, AnchorPoint anchorPoint, float textBoxWidth, float fontSize, float lineSpacing, glm::vec3 textColor){
    m_active_shader->setTextUniforms(m_windowSize.x, m_windowSize.y, textColor);
    m_textRenderer->renderUIText(font, text, anchorPosition, anchorPoint, textBoxWidth, fontSize, lineSpacing, textColor);
}

void Graphics::setGlobalData(glm::vec3 globalCoeffs){
    m_active_shader->setGlobalCoeffs(globalCoeffs);
}

void Graphics::setLights(std::vector<std::shared_ptr<Light>> lights){
    m_active_shader->setLights(lights);
}

void Graphics::clearLights(){
    m_active_shader->clearLights();
}

void Graphics::setWindowSize(glm::ivec2 windowSize){
    m_windowSize = windowSize;
}

glm::ivec2 Graphics::getWindowSize(){
    return m_windowSize;
}

void Graphics::setFramebufferSize(glm::ivec2 framebufferSize){
    m_framebufferSize = framebufferSize;
    glViewport(0, 0, m_framebufferSize.x, m_framebufferSize.y);
}

glm::ivec2 Graphics::getFramebufferSize(){
    return m_framebufferSize;
}


// WATERRRR //////

void Graphics::setWaterHeight(float height){
    m_water_height = height;

}

float Graphics::getWaterHeight(){
    return m_water_height;
}


void Graphics::bindReflection(){
    Debug::checkGLError();
    bindFBO(reflectionFBO, REFLECTION_W, REFLECTION_H);
    Debug::checkGLError();

}

void Graphics::bindRefraction(){
    Debug::checkGLError();

    bindFBO(refractionFBO, REFRACTION_W, REFRACTION_H);
    Debug::checkGLError();

}

GLuint Graphics::getReflectionTexture(){
    return reflectionTexture;
}

GLuint Graphics::getRefractionTexture(){
    return refractionTexture;
}

GLuint Graphics::getReflectionRBO(){
    return reflectionDepthBuffer;
}

GLuint Graphics::getRefractionDepth(){
    return refractionDepthTexture;
}


void Graphics::initializeReflectionFBO(){
    Debug::checkGLError();

    reflectionFBO = createFBO();
    Debug::checkGLError();

    reflectionTexture = createTextureAttachment(REFLECTION_W, REFLECTION_H);
    Debug::checkGLError();

    reflectionDepthBuffer = createDepthBuffer(REFLECTION_W, REFLECTION_H);
    Debug::checkGLError();

    unbindFBO();
    Debug::checkGLError();

}

void Graphics::initializeRefractionFBO(){
    Debug::checkGLError();

    refractionFBO = createFBO();
    Debug::checkGLError();

    refractionTexture = createTextureAttachment(REFRACTION_W, REFRACTION_H);
    Debug::checkGLError();

   // glDrawBuffer(GL_NONE);
    refractionDepthTexture = createDepthTextureAttachment(REFRACTION_W, REFRACTION_H);
    Debug::checkGLError();
    unbindFBO();
    Debug::checkGLError();


}
void Graphics::bindFBO(GLuint fbo, int width, int height){
    // make sure texture isnt already bound
    Debug::checkGLError();
    glBindTexture(GL_TEXTURE_2D, 0);
    Debug::checkGLError();

    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    Debug::checkGLError();

    glViewport(0, 0, width, height);
    Debug::checkGLError();

    auto fboStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
    if (fboStatus != GL_FRAMEBUFFER_COMPLETE){
        std::cout << "Framebuffer not complete: " << fboStatus << std::endl;
    }
    Debug::checkGLError();


}

void Graphics::unbindFBO(){
    Debug::checkGLError();

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glViewport(0, 0, m_framebufferSize.x, m_framebufferSize.y);
}

GLuint Graphics::createFBO(){
    GLuint fbo;
    glGenFramebuffers(1, &fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo);
    glDrawBuffer(GL_COLOR_ATTACHMENT0);
    return fbo;
}

GLuint Graphics::createTextureAttachment(int width, int height){
    GLuint texture;
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // adds texture attachement to currently bound fbo
    std::cout << "Color buffer: " << texture << std::endl;

    //glFramebufferTexture(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, texture, 0);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texture, 0);

    return texture;
}

GLuint Graphics::createDepthTextureAttachment(int width, int height){
    GLuint texture;
    glBindTexture(GL_TEXTURE_2D, 0);
    glGenTextures(1, &texture);
    glBindTexture(GL_TEXTURE_2D, texture);
    Debug::checkGLError();

    glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH24_STENCIL8, width, height, 0, GL_DEPTH_STENCIL, GL_UNSIGNED_INT_24_8, NULL);
    Debug::checkGLError();
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    std::cout << "Depth texture: " << texture << std::endl;

    Debug::checkGLError();

    glBindTexture(GL_TEXTURE_2D, 0);

    Debug::checkGLError();

    // adds texture attachement to currently bound fbo
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, texture, 0);
    Debug::checkGLError();

    return texture;
}

GLuint Graphics::createDepthBuffer(int width, int height){
    GLuint depthBuffer;
    glGenRenderbuffers(1, &depthBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    std::cout << "Depth buffer: " << depthBuffer << std::endl;
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);
    return depthBuffer;
}

void Graphics::deleteFBOS(){
    glDeleteFramebuffers(1, &reflectionFBO);
    glDeleteFramebuffers(1, &refractionFBO);

    glDeleteTextures(1, &reflectionTexture);
    glDeleteTextures(1, &refractionTexture);
    glDeleteTextures(1, &refractionDepthTexture);

    glDeleteRenderbuffers(1, &reflectionDepthBuffer);
}
