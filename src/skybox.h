#ifndef SKYBOX_H
#define SKYBOX_H
#include "graphics/shape.h"
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
#define EIGEN_DONT_VECTORIZE
#include "graphics/camera.h"
#include "graphics/shader.h"
#include <GL/glew.h>

//#ifdef __APPLE__
//#define GL_SILENCE_DEPRECATION
//#endif

#include <QOpenGLWidget>
#include <QElapsedTimer>
#include <QTimer>
#include <memory>


class skybox
{
public:
    skybox();
    void draw(Shader *skybox_shader, Camera  m_camera);
       void update(double deltaTime);

       GLuint getVAO(){
           return VAO;
       }

       GLuint getSkyboxTex(){
           return skybox_tex;
       }

        std::vector<GLfloat>  getVertices(){
           return m_vertices;
       }

        Eigen::Matrix4f getRotMat(){
            return m_rotation_mat;
        }

        std::vector<const char*> getTextureFiles(){
            return m_skyboxTextureFiles;
        }

        void initializeVAO();


   private:
       GLuint loadCubeMap(std::vector<const char*> textureFiles);


       GLuint VAO;
       GLuint VBO;
       GLuint EBO;
       GLuint skybox_tex;

       float SIZE = 100000000.f; //100000000.f;

       std::vector<int> m_vertex_indices = {
           // Right
               1, 2, 6,
               6, 5, 1,
               // Left
               0, 4, 7,
               7, 3, 0,
               // Top
               4, 5, 6,
               6, 7, 4,
               // Bottom
               0, 3, 2,
               2, 1, 0,
               // Back
               0, 1, 5,
               5, 4, 0,
               // Front
               3, 7, 6,
               6, 2, 3
       };

       std::vector<Eigen::Vector3i> m_faces = {
           // Right
              Eigen::Vector3i( 1, 2, 6),
             Eigen::Vector3i(  6, 5, 1),
               // Left
              Eigen::Vector3i( 0, 4, 7),
              Eigen::Vector3i( 7, 3, 0),
               // Top
             Eigen::Vector3i(  4, 5, 6),
             Eigen::Vector3i(  6, 7, 4),
               // Bottom
             Eigen::Vector3i(  0, 3, 2),
             Eigen::Vector3i(  2, 1, 0),
               // Back
             Eigen::Vector3i(  0, 1, 5),
             Eigen::Vector3i(  5, 4, 0),
               // Front
             Eigen::Vector3i(  3, 7, 6),
             Eigen::Vector3i(  6, 2, 3)
       };
       std::vector<GLfloat> m_vertices = {
           -SIZE,  SIZE, -SIZE,
           -SIZE, -SIZE, -SIZE,
            SIZE, -SIZE, -SIZE,
            SIZE, -SIZE, -SIZE,
            SIZE,  SIZE, -SIZE,
           -SIZE,  SIZE, -SIZE,

           -SIZE, -SIZE,  SIZE,
           -SIZE, -SIZE, -SIZE,
           -SIZE,  SIZE, -SIZE,
           -SIZE,  SIZE, -SIZE,
           -SIZE,  SIZE,  SIZE,
           -SIZE, -SIZE,  SIZE,

            SIZE, -SIZE, -SIZE,
            SIZE, -SIZE,  SIZE,
            SIZE,  SIZE,  SIZE,
            SIZE,  SIZE,  SIZE,
            SIZE,  SIZE, -SIZE,
            SIZE, -SIZE, -SIZE,

           -SIZE, -SIZE,  SIZE,
           -SIZE,  SIZE,  SIZE,
            SIZE,  SIZE,  SIZE,
            SIZE,  SIZE,  SIZE,
            SIZE, -SIZE,  SIZE,
           -SIZE, -SIZE,  SIZE,

           -SIZE,  SIZE, -SIZE,
            SIZE,  SIZE, -SIZE,
            SIZE,  SIZE,  SIZE,
            SIZE,  SIZE,  SIZE,
           -SIZE,  SIZE,  SIZE,
           -SIZE,  SIZE, -SIZE,

           -SIZE, -SIZE, -SIZE,
           -SIZE, -SIZE,  SIZE,
            SIZE, -SIZE, -SIZE,
            SIZE, -SIZE, -SIZE,
           -SIZE, -SIZE,  SIZE,
            SIZE, -SIZE,  SIZE
       };

       std::vector<Eigen::Vector3f> m_vertices_eigen = {
          Eigen::Vector3f( -SIZE,  SIZE, -SIZE),
          Eigen::Vector3f( -SIZE, -SIZE, -SIZE),
          Eigen::Vector3f(  SIZE, -SIZE, -SIZE),
          Eigen::Vector3f(  SIZE, -SIZE, -SIZE),
          Eigen::Vector3f(  SIZE,  SIZE, -SIZE),
          Eigen::Vector3f( -SIZE,  SIZE, -SIZE),

          Eigen::Vector3f( -SIZE, -SIZE,  SIZE),
         Eigen::Vector3f(  -SIZE, -SIZE, -SIZE),
          Eigen::Vector3f( -SIZE,  SIZE, -SIZE),
          Eigen::Vector3f( -SIZE,  SIZE, -SIZE),
         Eigen::Vector3f(  -SIZE,  SIZE,  SIZE),
         Eigen::Vector3f(  -SIZE, -SIZE,  SIZE),

          Eigen::Vector3f(  SIZE, -SIZE, -SIZE),
          Eigen::Vector3f(  SIZE, -SIZE,  SIZE),
          Eigen::Vector3f(  SIZE,  SIZE,  SIZE),
           Eigen::Vector3f( SIZE,  SIZE,  SIZE),
          Eigen::Vector3f(  SIZE,  SIZE, -SIZE),
          Eigen::Vector3f(  SIZE, -SIZE, -SIZE),

          Eigen::Vector3f( -SIZE, -SIZE,  SIZE),
          Eigen::Vector3f( -SIZE,  SIZE,  SIZE),
          Eigen::Vector3f(  SIZE,  SIZE,  SIZE),
          Eigen::Vector3f(  SIZE,  SIZE,  SIZE),
          Eigen::Vector3f(  SIZE, -SIZE,  SIZE),
          Eigen::Vector3f( -SIZE, -SIZE,  SIZE),

          Eigen::Vector3f( -SIZE,  SIZE, -SIZE),
           Eigen::Vector3f( SIZE,  SIZE, -SIZE),
          Eigen::Vector3f(  SIZE,  SIZE,  SIZE),
          Eigen::Vector3f(  SIZE,  SIZE,  SIZE),
          Eigen::Vector3f( -SIZE,  SIZE,  SIZE),
          Eigen::Vector3f( -SIZE,  SIZE, -SIZE),

          Eigen::Vector3f( -SIZE, -SIZE, -SIZE),
          Eigen::Vector3f( -SIZE, -SIZE,  SIZE),
          Eigen::Vector3f(  SIZE, -SIZE, -SIZE),
          Eigen::Vector3f(  SIZE, -SIZE, -SIZE),
          Eigen::Vector3f( -SIZE, -SIZE,  SIZE),
          Eigen::Vector3f(  SIZE, -SIZE,  SIZE)
       };

       std::vector<const char*> m_skyboxTextureFiles =

//       {
//           "/Users/jesswan/Desktop/cs2240/ocean-simulation/resources/images/foam3.png",
//           "/Users/jesswan/Desktop/cs2240/ocean-simulation/resources/images/foam3.png",
//           "/Users/jesswan/Desktop/cs2240/ocean-simulation/resources/images/foam3.png",
//           "/Users/jesswan/Desktop/cs2240/ocean-simulation/resources/images/foam3.png",
//           "/Users/jesswan/Desktop/cs2240/ocean-simulation/resources/images/foam3.png",
//           "/Users/jesswan/Desktop/cs2240/ocean-simulation/resources/images/foam3.png"

//       };




        {":resources/images/px.png",
     ":resources/images/nx.png",
               ":resources/images/py.png",
               ":resources/images/ny.png",
               ":resources/images/pz.png",
               ":resources/images/nz.png",
   };

       float ROTATE_SPEED = .01f; // 1 degree per sec
       float m_rotation = 0.f;
       Eigen::Matrix4f m_rotation_mat = Eigen::Matrix4f::Identity();

       Shape sky_shape;
};

#endif // SKYBOX_H
