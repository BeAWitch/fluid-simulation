#include <iostream>
#include "glut.h"
#include "Fluid.h"

using namespace std;

int scale = 5;
int N = 128;
Fluid fluid(N, 0.05, 0, 0);

// 显示函数
void display()
{
    glClear(GL_COLOR_BUFFER_BIT);

    fluid.fluidCubeStep();
    fluid.fadeDensity();

    const vector<float>& density = fluid.getDensity();

    glBegin(GL_QUADS);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            float d = density[fluid.IX(i, j)];
            d = min(d * 0.02f, 1.0f);

            //float r = fmod((i * 13 + j * 7 + d * 1000), 255) / 255.0f;
            //float g = fmod((i * 3 + j * 17 + d * 500), 255) / 255.0f;
            //float b = fmod((i * 23 + j * 5 + d * 800), 255) / 255.0f;
            //glColor3f(r * d, g * d, b * d);

            glColor3f(d, d, d);

            float x = i * scale;
            float y = j * scale;

            glVertex2f(x, y);
            glVertex2f(x + scale, y);
            glVertex2f(x + scale, y + scale);
            glVertex2f(x, y + scale);
        }
    }
    glEnd();

    glutSwapBuffers();
}


// 窗口调整
void reshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, N * scale, 0, N * scale);
    glMatrixMode(GL_MODELVIEW);
}

int lastMouseX, lastMouseY;
bool isDragging = false;
void mouseFunc(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON) {
        if (state == GLUT_DOWN) {
            isDragging = true;
            lastMouseX = x;
            lastMouseY = y;
        }
        else if (state == GLUT_UP) {
            isDragging = false;
        }
    }
}

void motionCallback(int x, int y) {
    if (isDragging) {
        int dx = x - lastMouseX;
        int dy = lastMouseY - y;

        lastMouseX = x;
        lastMouseY = y;

        int windowHeight = glutGet(GLUT_WINDOW_HEIGHT);
        int gridX = x / scale;
        int gridY = (windowHeight - y) / scale;

        if (gridX >= 0 && gridX < N && gridY >= 0 && gridY < N)
        {
            fluid.addDensity(gridX, gridY, 100.0f);
            fluid.addVelocity(gridX, gridY, dx, dy);
        }

        glutPostRedisplay();
    }
}

void idleFunc()
{
    glutPostRedisplay();
}

int main(int argc, char** argv) 
{
    // 初始化GLUT
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(N * scale, N * scale);
    glutCreateWindow("fluid");

    // 注册回调函数
    glutDisplayFunc(display);
    glutIdleFunc(idleFunc);
    //glutKeyboardFunc(keyboard);
    glutMouseFunc(mouseFunc);
    glutMotionFunc(motionCallback);
    glutReshapeFunc(reshape);

    // 设置背景色
    glClearColor(0.1, 0.1, 0.1, 1.0);

    glutMainLoop();
    return 0;
}