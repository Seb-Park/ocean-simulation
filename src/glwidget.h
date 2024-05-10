#pragma once

#include "skybox.h"
#include <GL/glew.h>


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

    void makeFBO();
    void makeFBO1();

    void initCaustics();
    void paintCaustics();

    TextureData loadTextureFromFile(const char *path);
    GLuint loadCubeMap(std::vector<const char*> textureFiles);


private slots:
    // Physics Tick
    void tick();

private:
    ARAP    m_arap;
    Camera  m_camera;
    Shader *m_defaultShader;
    Shader *m_pointShader;
    Shader *m_texture_shader;
    Shader *m_causticsShader;

    Shader *m_colorShader;
    Shader *m_foamShader;
    Shader *m_skyboxShader;



    GLuint m_fullscreen_vbo;
    GLuint m_fullscreen_vao;
    QImage m_ground_image;
    GLuint m_ground_texture;

    int m_fbo_width;
    int m_fbo_height;
    float m_devicePixelRatio;

    GLuint m_fbo;
    GLuint m_fbo_texture;
    GLuint m_fbo_renderbuffer;

    GLuint m_fbo1;
    GLuint m_fbo_texture1;
    GLuint m_fbo_renderbuffer1;

    GLuint m_defaultFBO;

    GLuint m_floor_vbo;
    GLuint m_floor_vao;

    float m_movementScaling;
    float m_vertexSelectionThreshold;
    float m_vSize;

    // FOAM
    GLuint m_halftone_tex;
    GLuint m_foam_tex;


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

    skybox m_skybox;
};
