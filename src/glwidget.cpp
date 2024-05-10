#include "glwidget.h"

#include <QApplication>
#include <QKeyEvent>
#include <iostream>

#define SPEED 1.5
#define ROTATE_SPEED 0.0025

using namespace std;
using namespace Eigen;

GLWidget::GLWidget(QWidget *parent) :
    QOpenGLWidget(parent),
    m_arap(),
    m_camera(),
    m_defaultShader(),
    m_pointShader(),
    m_vSize(),
    m_movementScaling(),
    m_vertexSelectionThreshold(),
    // Movement
    m_deltaTimeProvider(),
    m_intervalTimer(),
    // Timing
    m_forward(),
    m_sideways(),
    m_vertical(),
    // Mouse handler stuff
    m_lastX(),
    m_lastY(),
    m_leftCapture(false),
    m_rightCapture(false),
    m_rightClickSelectMode(SelectMode::None),
    m_lastSelectedVertex(-1)
{
    // GLWidget needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    QApplication::setOverrideCursor(Qt::ArrowCursor);

    // GLWidget needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // Function tick() will be called once per interva
    connect(&m_intervalTimer, SIGNAL(timeout()), this, SLOT(tick()));
}

GLWidget::~GLWidget()
{
    if (m_defaultShader != nullptr) delete m_defaultShader;
    if (m_pointShader   != nullptr) delete m_pointShader;
}

// ================== Basic OpenGL Overrides

void GLWidget::initializeGL()
{
    // Initialize GL extension wrangler
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK) fprintf(stderr, "Error while initializing GLEW: %s\n", glewGetErrorString(err));
    fprintf(stdout, "Successfully initialized GLEW %s\n", glewGetString(GLEW_VERSION));

    // Set clear color to white
    glClearColor(1, 0.98f, 0.85f, 1);

    // Enable depth-testing and backface culling
    glEnable(GL_DEPTH_TEST);
//    glEnable(GL_CULL_FACE);
//    glCullFace(GL_BACK);
//    glShadeModel(GL_SMOOTH);


    // Initialize shaders
    m_defaultShader = new Shader(":resources/shaders/shader.vert",      ":resources/shaders/shader.frag");
    m_pointShader   = new Shader(":resources/shaders/anchorPoint.vert", ":resources/shaders/anchorPoint.geom", ":resources/shaders/anchorPoint.frag");
//    m_texture_shader = new Shader(":/resources/shaders/texture.vert", ":/resources/shaders/texture.frag");
    m_colorShader   = new Shader(":resources/shaders/color.vert",      ":resources/shaders/color.frag");
    m_causticsShader = new Shader(":resources/shaders/caustics.vert",      ":resources/shaders/caustics.frag");

    initCaustics();
    // INITIALIZE TEXTURE STUFF

    // Prepare filepath


//    // TASK 11: Fix this "fullscreen" quad's vertex data
//    // TASK 12: Play around with different values!
//    // TASK 13: Add UV coordinates
//    std::vector<GLfloat> fullscreen_quad_data =
//    { //     POSITIONS    //    UVs    //
//        -1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
//        -1.0f, -1.0f, 0.0f, 0.0f, 0.0f,
//         1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
//         1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
//        -1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
//         1.0f, -1.0f, 0.0f, 1.0f, 0.0f
//    };

    // START FBO STUFF

    m_devicePixelRatio = this->devicePixelRatio();

    m_defaultFBO = 3;
    m_fbo_width = size().width() * m_devicePixelRatio;
    m_fbo_height = size().height() * m_devicePixelRatio;

    makeFBO();
    makeFBO1();

    // FBO STUFF END

    m_arap.initSkyPlane(":/resources/images/sky_clouds.png", 2, m_defaultShader);
    m_arap.initGroundPlane(":/resources/images/daniel.jpg", 2, m_defaultShader);

//    // Generate and bind a VBO and a VAO for a fullscreen quad
//    glGenBuffers(1, &m_fullscreen_vbo);
//    glBindBuffer(GL_ARRAY_BUFFER, m_fullscreen_vbo);
//    glBufferData(GL_ARRAY_BUFFER, fullscreen_quad_data.size()*sizeof(GLfloat), fullscreen_quad_data.data(), GL_STATIC_DRAW);
//    glGenVertexArrays(1, &m_fullscreen_vao);
//    glBindVertexArray(m_fullscreen_vao);

//    // TASK 14: modify the code below to add a second attribute to the vertex attribute array
//    glEnableVertexAttribArray(0);
//    glEnableVertexAttribArray(1);
//    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), nullptr);
//    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 5 * sizeof(GLfloat), reinterpret_cast<void *>(3 * sizeof(GLfloat)));

//    // Unbind the fullscreen quad's VBO and VAO
//    glBindBuffer(GL_ARRAY_BUFFER, 0);
//    glBindVertexArray(0);

    // END INITIALIZE TEXTURE STUFF


    // Initialize ARAP, and get parameters needed to decide the camera position, etc
    Vector3f coeffMin, coeffMax;
    m_arap.init(coeffMin, coeffMax);

    Vector3f center = (coeffMax + coeffMin) / 2.0f;
    float extentLength  = (coeffMax - coeffMin).norm();

    // Screen-space size of vertex points
    m_vSize = 0.005 * extentLength;

    // Scale all movement by this amount
    m_movementScaling = extentLength * 0.5;

    // When raycasting, select closest vertex within this distance
    m_vertexSelectionThreshold = extentLength * 0.025;

    // Note for maintainers: Z-up
    float fovY = 120;
    float nearPlane = 0.1f;
    float farPlane  = 3 * extentLength;

    // Initialize camera with a reasonable transform
    Eigen::Vector3f eye    = center - Eigen::Vector3f::UnitZ() * extentLength;
    Eigen::Vector3f target = center;
    m_camera.lookAt(eye, target);
    m_camera.setOrbitPoint(target);
    m_camera.setPerspective(120, width() / static_cast<float>(height()), nearPlane, farPlane);

    m_deltaTimeProvider.start();
    m_intervalTimer.start(1000 / 60);

}

void GLWidget::paintCaustics() {
    glClearColor(0.68f, 0.58f, 0.38f, 1);
//    glClearColor(0., 0., 0., 1);

    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo1);
//    glBindFramebuffer(GL_FRAMEBUFFER, m_defaultFBO);

    // Clear Screen
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );

//    // Bind the shader
//    m_colorShader->bind();

////    m_colorShader->setUniform("asdf", m_arap.minCorner[0], m_arap.maxCorner[0]);
////    // Bind VAO
//    glBindVertexArray(m_floor_vao);

////    // Draw the VAO
//    glDrawArrays(GL_TRIANGLES, 0, 3);


    m_colorShader->bind();
    //
    m_colorShader->setUniform("proj", m_camera.getProjection());
    m_colorShader->setUniform("view", m_camera.getView());
    Eigen::Matrix4f inverseView = m_camera.getView().inverse();
    m_colorShader->setUniform("inverseView", inverseView);
    //
    m_colorShader->setUniform("widthBounds", m_arap.minCorner[0], m_arap.maxCorner[0]);
    m_colorShader->setUniform("lengthBounds", m_arap.minCorner[2], m_arap.maxCorner[2]);
    m_arap.draw(m_colorShader, GL_TRIANGLES);

    // Unbind the VAO
    glBindVertexArray(0);

    // Unbind the shader
    glUseProgram(0);

    m_causticsShader->bind();

//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_texture);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );
    glDisable( GL_DEPTH_TEST );

    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, m_fbo_texture1);
    glUniform1i(glGetUniformLocation(m_causticsShader ->id(), "normSamp"), 2);

    m_arap.m_causticsShape.draw(m_causticsShader, GL_TRIANGLES);
    glBindVertexArray(0);
    glUseProgram(0);
}

void GLWidget::initCaustics() {
    glGenBuffers(1, &m_floor_vbo);

    // Task 6: Bind the VBO you created here

    // TA SOLUTION
    glBindBuffer(GL_ARRAY_BUFFER, m_floor_vbo);

    // Task 8: Construct your std::vector of triangle data here | Task 10: Add colors to your triangle here

    // TA SOLUTION
    std::vector<GLfloat> triangle =
    {   //    POSITIONS    //    COLORS    //
        -0.5f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f,
         0.5f, -0.5f, 0.0f, 0.0f, 1.0f, 0.0f,
         0.0f,  0.5f, 0.0f, 0.0f, 0.0f, 1.0f
    };

    // Task 9: Pass the triangle vector into your VBO here

    // TA SOLUTION
    glBufferData(GL_ARRAY_BUFFER, triangle.size() * sizeof(GLfloat), triangle.data(), GL_STATIC_DRAW);

    // ================== Vertex Array Objects

    // Task 11: Generate a VAO here and store it in m_vao

    // TA SOLUTION
    glGenVertexArrays(1, &m_floor_vao);

    // Task 12: Bind the VAO you created here

    // TA SOLUTION
    glBindVertexArray(m_floor_vao);

    // Task 13: Add position and color attributes to your VAO here

    // TA SOLUTION
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    GLsizei stride = 6 * sizeof(GLfloat);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, nullptr);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, reinterpret_cast<void *>(3 * sizeof(GLfloat)));

    // ================== Returning to Default State

    // Task 14: Unbind your VBO and VAO here

    // TA SOLUTION
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

//void GLWidget::paintTexture(GLuint texture, bool filtered){
////    glUseProgram(m_texture_shader->id());
//    m_texture_shader->bind();
//    // TASK 32: Set your bool uniform on whether or not to filter the texture drawn
////    glUniform1i(glGetUniformLocation(m_texture_shader->id(), "filtered"), filtered);
//    m_texture_shader->setUniform("filtered", filtered);
//    // TASK 10: Bind "texture" to slot 0
//    glActiveTexture(GL_TEXTURE0);
//    glBindTexture(GL_TEXTURE_2D, texture);
//    glBindVertexArray(m_fullscreen_vao);
////    std::cout << texture << std::endl;
//    glDrawArrays(GL_TRIANGLES, 0, 6);
//    glBindTexture(GL_TEXTURE_2D, 0);
//    glBindVertexArray(0);
//    glUseProgram(0);
//    m_texture_shader->unbind();
//}

void GLWidget::makeFBO() {
    // Task 19
    glActiveTexture(GL_TEXTURE2);
    glGenTextures(1, &m_fbo_texture);
    glBindTexture(GL_TEXTURE_2D, m_fbo_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_fbo_width, m_fbo_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    // Task 20
    glGenRenderbuffers(1, &m_fbo_renderbuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, m_fbo_renderbuffer);
    glRenderbufferStorage(GL_RENDERBUFFER,
                          GL_DEPTH24_STENCIL8,
                          m_fbo_width,
                          m_fbo_height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    // Task 18
    glGenFramebuffers(1, &m_fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);

    // Task 21
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_fbo_texture, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, m_fbo_renderbuffer);

    // Task 22
    glBindFramebuffer(GL_FRAMEBUFFER, m_defaultFBO);
}

void GLWidget::makeFBO1() {
    // Task 19
    glActiveTexture(GL_TEXTURE2);
    glGenTextures(1, &m_fbo_texture1);
    glBindTexture(GL_TEXTURE_2D, m_fbo_texture1);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_fbo_width, m_fbo_height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);

    // Task 20
    glGenRenderbuffers(1, &m_fbo_renderbuffer1);
    glBindRenderbuffer(GL_RENDERBUFFER, m_fbo_renderbuffer1);
    glRenderbufferStorage(GL_RENDERBUFFER,
                          GL_DEPTH24_STENCIL8,
                          m_fbo_width,
                          m_fbo_height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);

    // Task 18
    glGenFramebuffers(1, &m_fbo1);
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo1);

    // Task 21
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_fbo_texture1, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, m_fbo_renderbuffer1);

    // Task 22
    glBindFramebuffer(GL_FRAMEBUFFER, m_defaultFBO);
}

void GLWidget::paintGL()
{
    paintCaustics();
    glClearColor(0.56f, 0.69f, 0.74f, 1);
    glBindFramebuffer(GL_FRAMEBUFFER, m_defaultFBO);
//    return;
//    paintTexture(m_ground_texture, false);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable( GL_BLEND );
    glEnable(GL_DEPTH_TEST);

    m_defaultShader->bind();
    m_defaultShader->setUniform("proj", m_camera.getProjection());
    m_defaultShader->setUniform("view", m_camera.getView());
    Eigen::Matrix4f inverseView = m_camera.getView().inverse();
    m_defaultShader->setUniform("inverseView", inverseView);
    m_defaultShader->setUniform("widthBounds", m_arap.minCorner[0], m_arap.maxCorner[0]);
    m_defaultShader->setUniform("lengthBounds", m_arap.minCorner[2], m_arap.maxCorner[2]);

    // Set the ocean floor to the painted caustics
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, m_fbo_texture);
    glUniform1i(glGetUniformLocation(m_defaultShader->id(), "groundSampler"), 2);

    m_arap.draw(m_defaultShader, GL_TRIANGLES);
    m_defaultShader->unbind();

    glClear(GL_DEPTH_BUFFER_BIT);

//    m_pointShader->bind();
//    m_pointShader->setUniform("proj",   m_camera.getProjection());
//    m_pointShader->setUniform("view",   m_camera.getView());
//    m_pointShader->setUniform("vSize",  m_vSize);
//    m_pointShader->setUniform("width",  width());
//    m_pointShader->setUniform("height", height());
//    m_arap.draw(m_pointShader, GL_POINTS);
//    m_pointShader->unbind();
}

void GLWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    m_camera.setAspect(static_cast<float>(w) / h);
    glDeleteTextures(1, &m_fbo_texture);
    glDeleteRenderbuffers(1, &m_fbo_renderbuffer);
    glDeleteFramebuffers(1, &m_fbo);
    glDeleteTextures(1, &m_fbo_texture1);
    glDeleteRenderbuffers(1, &m_fbo_renderbuffer1);
    glDeleteFramebuffers(1, &m_fbo1);
    makeFBO();
    makeFBO1();
}

// ================== Event Listeners

Eigen::Vector3f GLWidget::transformToWorldRay(int x, int y)
{
    Eigen::Vector4f clipCoords = Eigen::Vector4f(
                (float(x) / width()) * 2.f - 1.f,
                1.f - (float(y) / height()) * 2.f,
                -1.f,
                1.f);

    Eigen::Vector4f transformed_coords = m_camera.getProjection().inverse() * clipCoords;
    transformed_coords = Eigen::Vector4f(transformed_coords.x(), transformed_coords.y(), -1.f, 0.f);
    transformed_coords = m_camera.getView().inverse() * transformed_coords;

    return Eigen::Vector3f(transformed_coords.x(), transformed_coords.y(), transformed_coords.z()).normalized();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    // Get current mouse coordinates
    const int currX = event->position().x();
    const int currY = event->position().y();

    // Get closest vertex to ray
    const Vector3f ray = transformToWorldRay(currX, currY);
    const int closest_vertex = m_arap.getClosestVertex(m_camera.getPosition(), ray, m_vertexSelectionThreshold);

    // Switch on button
    switch (event->button()) {
    case Qt::MouseButton::RightButton: {
        // Capture
        m_rightCapture = true;
        // Anchor/un-anchor the vertex
        m_rightClickSelectMode = m_arap.select(m_pointShader, closest_vertex);
        break;
    }
    case Qt::MouseButton::LeftButton: {
        // Capture
        m_leftCapture = true;
        // Select this vertex
        m_lastSelectedVertex = closest_vertex;
        break;
    }
    default: break;
    }

    // Set last mouse coordinates
    m_lastX = currX;
    m_lastY = currY;
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    // Return if neither mouse button is currently held down
    if (!(m_leftCapture || m_rightCapture)) {
        return;
    }

    // Get current mouse coordinates
    const int currX = event->position().x();
    const int currY = event->position().y();

    // Find ray
    const Vector3f ray = transformToWorldRay(event->position().x(), event->position().y());
    Vector3f pos;

    // If right is held down
    if (m_rightCapture) {
        // Get closest vertex to ray
        const int closest_vertex = m_arap.getClosestVertex(m_camera.getPosition(), ray, m_vertexSelectionThreshold);

        // Anchor/un-anchor the vertex
        if (m_rightClickSelectMode == SelectMode::None) {
            m_rightClickSelectMode = m_arap.select(m_pointShader, closest_vertex);
        } else {
            m_arap.selectWithSpecifiedMode(m_pointShader, closest_vertex, m_rightClickSelectMode);
        }

        return;
    }

    // If the selected point is an anchor point
    if (m_lastSelectedVertex != -1 && m_arap.getAnchorPos(m_lastSelectedVertex, pos, ray, m_camera.getPosition())) {
        // Move it
        m_arap.move(m_lastSelectedVertex, pos);
    } else {
        // Rotate the camera
        const int deltaX = currX - m_lastX;
        const int deltaY = currY - m_lastY;
        if (deltaX != 0 || deltaY != 0) {
            m_camera.rotate(deltaY * ROTATE_SPEED, -deltaX * ROTATE_SPEED);
        }
    }

    // Set last mouse coordinates
    m_lastX = currX;
    m_lastY = currY;
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_leftCapture = false;
    m_lastSelectedVertex = -1;

    m_rightCapture = false;
    m_rightClickSelectMode = SelectMode::None;
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
    float zoom = 1 - event->pixelDelta().y() * 0.1f / 120.f;
    m_camera.zoom(zoom);
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  += SPEED; break;
    case Qt::Key_S: m_forward  -= SPEED; break;
    case Qt::Key_A: m_sideways -= SPEED; break;
    case Qt::Key_D: m_sideways += SPEED; break;
    case Qt::Key_F: m_vertical -= SPEED; break;
    case Qt::Key_R: m_vertical += SPEED; break;
    case Qt::Key_C: m_camera.toggleIsOrbiting(); break;
    case Qt::Key_Equal: m_vSize *= 11.0f / 10.0f; break;
    case Qt::Key_Minus: m_vSize *= 10.0f / 11.0f; break;
    case Qt::Key_Escape: QApplication::quit();
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  -= SPEED; break;
    case Qt::Key_S: m_forward  += SPEED; break;
    case Qt::Key_A: m_sideways += SPEED; break;
    case Qt::Key_D: m_sideways -= SPEED; break;
    case Qt::Key_F: m_vertical += SPEED; break;
    case Qt::Key_R: m_vertical -= SPEED; break;
    }
}

// ================== Physics Tick

void GLWidget::tick()
{
    float deltaSeconds = m_deltaTimeProvider.restart() / 1000.f;
    m_arap.update(deltaSeconds);

    // Move camera
    auto look = m_camera.getLook();
    look.y() = 0;
    look.normalize();
    Eigen::Vector3f perp(-look.z(), 0, look.x());
    Eigen::Vector3f moveVec = m_forward * look.normalized() + m_sideways * perp.normalized() + m_vertical * Eigen::Vector3f::UnitY();
    moveVec *= m_movementScaling;
    moveVec *= deltaSeconds;
    m_camera.move(moveVec);

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();
}
