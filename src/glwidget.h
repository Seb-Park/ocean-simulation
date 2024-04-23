#pragma once

#ifdef __APPLE__
#define GL_SILENCE_DEPRECATION
#endif

#include "arap.h"
#include "graphics/camera.h"
#include "graphics/shader.h"

#include <QOpenGLWidget>
#include <QElapsedTimer>
#include <QTimer>
#include <memory>

struct TextureData{
    GLuint textureID;
    int width;
    int height;
};

class GLWidget : public QOpenGLWidget
{
    Q_OBJECT

public:
    GLWidget(QWidget *parent = nullptr);
    ~GLWidget();

private:
    static const int FRAMES_TO_AVERAGE = 30;

    Eigen::Vector3f transformToWorldRay(int x, int y);

    // Basic OpenGL Overrides
    void initializeGL()         override;
    void paintGL()              override;
    void paintTexture(GLuint texture, bool filtered);
    void resizeGL(int w, int h) override;

    // Event Listeners
    void mousePressEvent  (QMouseEvent *event) override;
    void mouseMoveEvent   (QMouseEvent *event) override;
    void mouseReleaseEvent(QMouseEvent *event) override;
    void wheelEvent       (QWheelEvent *event) override;
    void keyPressEvent    (QKeyEvent   *event) override;
    void keyReleaseEvent  (QKeyEvent   *event) override;

    TextureData loadTextureFromFile(const char *path);
    void makeFBO();
    void initFullScreenQuad();
    void paintTexture(GLuint texture, bool postProcessOn);





private slots:
    // Physics Tick
    void tick();

private:
    ARAP    m_arap;
    Camera  m_camera;
    Shader *m_defaultShader;
    Shader *m_pointShader;
    Shader *m_texture_shader;

    GLuint m_fullscreen_vbo;
    GLuint m_fullscreen_vao;
    QImage m_ground_image;
    GLuint m_ground_texture;

    float m_movementScaling;
    float m_vertexSelectionThreshold;
    float m_vSize;

    // Timing
    QElapsedTimer m_deltaTimeProvider; // For measuring elapsed time
    QTimer        m_intervalTimer;     // For triggering timed events

    // Movement
    int m_forward;
    int m_sideways;
    int m_vertical;

    // Mouse handler stuff
    int m_lastX;
    int m_lastY;
    bool m_leftCapture;
    bool m_rightCapture;
    SelectMode m_rightClickSelectMode;
    int m_lastSelectedVertex = -1;


    int m_devicePixelRatio;
      GLuint m_defaultFBO;
      int m_fbo_width;
      int m_fbo_height;
      int m_screen_width;
      int m_screen_height;

        GLuint m_fullscreen_vbo;
        GLuint m_fullscreen_vao;
        GLuint m_fbo;
        GLuint m_fbo_texture;
        GLuint m_fbo_renderbuffer;
    GLuint m_texture0;

    std::vector<GLfloat> fullscreen_quad_data =
     { //     POSITIONS    //
         -1.f,  1.f, 0.0f,
          0.f, 1.f, //uv
         -1.f, -1.f, 0.0f,
          0.f, 0.f, //uv
          1.f, -1.f, 0.0f,
          1.f, 0.f, //uv
          1.f,  1.f, 0.0f,
          1.f, 1.f, //uv
         -1.f,  1.f, 0.0f,
          0.f, 1.f, //uv
          1.f, -1.f, 0.0f,
          1.f, 0.f //uv
     };

};
