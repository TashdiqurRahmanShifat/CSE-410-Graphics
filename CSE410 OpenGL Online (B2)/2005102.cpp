#include <GL/glut.h>
#include <cmath>  // For cosf(), sinf(), etc.


GLfloat eyex = 10, eyey = 10, eyez = 10;
GLfloat centerx = 0, centery = 0, centerz = 0;
GLfloat upx = 0, upy = 1, upz = 0;
float rotang = 0,bladespeed=4,bladeang=0;

#define PI 3.14159265358979323846 

void initGL() {
    glClearColor(0, 0, 0, 1);
    glEnable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, 1.0, 1.0, 100.0);
}

void reshapeListener(int w, int h) {
    if (h == 0) h = 1;
    float aspect = (float)w / h;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60.0, aspect, 1.0, 100.0);
}

void drawCube() {
    float width = 1.5f, height = 2.0f, depth = 3.0f; // Cube dimensions
    glBegin(GL_QUADS);

    // Front face (z = +depth)
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f( width,  height,  depth);
    glVertex3f(-width,  height,  depth);
    glVertex3f(-width, -height,  depth);
    glVertex3f( width, -height,  depth);

    // // Back face (z = -depth)
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f( width, -height, -depth);
    glVertex3f(-width, -height, -depth);
    glVertex3f(-width,  height, -depth);
    glVertex3f( width,  height, -depth);

    // // Left face (x = -width)
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f(-width*2,  height,  depth/2);
    glVertex3f(-width*2,  height, -depth/2);
    glVertex3f(-width*2, -height, -depth/2);
    glVertex3f(-width*2, -height,  depth/2);

    // Right face (x = +width)
    // glColor3f(0.0f,1.0f,0.0f);
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f( width*2,  height, -depth/2);
    glVertex3f( width*2,  height,  depth/2);
    glVertex3f( width*2, -height,  depth/2);
    glVertex3f( width*2, -height, -depth/2);

    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(-width,  height,  depth);
    glVertex3f(-width, -height,  depth);
    glVertex3f( -width*2, -height,  depth/2);
    glVertex3f( -width*2, height, depth/2);

    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(width,  height,  depth);
    glVertex3f(width, -height,  depth);
    glVertex3f( width*2, -height,  depth/2);
    glVertex3f( width*2, height, depth/2);

    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(width,  height,  -depth);
    glVertex3f(width, -height,  -depth);
    glVertex3f( width*2, -height,  -depth/2);
    glVertex3f( width*2, height, -depth/2);

    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(-width,  -height,  -depth);
    glVertex3f(-width, height,  -depth);
    glVertex3f( -width*2, height,  -depth/2);
    glVertex3f( -width*2, -height, -depth/2);


    // // Top face (y = +height)
    // glColor3f(0.675f, 0.455f, 0.910f); // Purple
    // glVertex3f( width,  height, -depth);
    // glVertex3f(-width,  height, -depth);
    // glVertex3f(-width,  height,  depth);
    // glVertex3f( width,  height,  depth);

    // // Bottom face (y = -height)
    // glColor3f(0.1f, 0.1f, 0.1f); // Dark gray / black
    // glVertex3f( width, -height,  depth);
    // glVertex3f(-width, -height,  depth);
    // glVertex3f(-width, -height, -depth);
    // glVertex3f( width, -height, -depth);

    glEnd();
}





void drawAxes()
{
    glLineWidth(3); // Set line thickness

    glBegin(GL_LINES);

    // X axis (red)
    glColor3f(1, 1, 1);
    glVertex3f(-100, 0, 0);
    glVertex3f(100, 0, 0);

    // Y axis (green)
    glColor3f(1, 1, 1);
    glVertex3f(0, -100, 0);
    glVertex3f(0, 100, 0);

    // Z axis (blue)
    glColor3f(1, 1, 1);
    glVertex3f(0, 0, -100);
    glVertex3f(0, 0, 100);

    glEnd();
}

void drawCylinder(float baseRadius, float topRadius, float height)
{
    GLUquadric *quad = gluNewQuadric();
    glColor3f(0.867f,0.922f,0.122f);
    gluCylinder(quad, baseRadius, topRadius, height, 20, 20);
    gluDeleteQuadric(quad);
}
// void drawBlade(){
//     glBegin(GL_TRIANGLES);
//     glColor3f(1.0f, 0.0f, 0.0f);
//     glVertex3f(0,0,2.05);
//     glVertex3f(0,2,2.05);
//     glVertex3f(1,2,2.05);
//     glEnd();

// }
void drawBlade(){
    // glBegin(GL_QUADS);
    // glColor3f(1.0f, 0.0f, 0.0f);
    // glVertex3f(-0,5,2.05);
    // glVertex3f(0,2,2.05);
    // glVertex3f(1,2,2.05);
    // glVertex3f(1,2,2.05);
    // glEnd();
    float width = 0.2f, height = 6.0f, depth = 2.05f; // Cube dimensions
    glBegin(GL_QUADS);

    // Front face (z = +depth)
    glColor3f(1.0f, 0.0f, 0.0f);
    glVertex3f( width,  height,  depth);
    glVertex3f(-width,  height,  depth);
    glVertex3f(-width, -height,  depth);
    glVertex3f( width, -height,  depth);
    glEnd();

}
void drawRotatingBlade(){
   for(int i=0;i<6;i++){
    glPushMatrix();
    glRotatef(i*60,0,0,1);
    drawBlade();
    glPopMatrix();
   }
}


void drawCube2() {
    float width = .5f, height = 2.0f, depth = .9f; // Cube dimensions
    glBegin(GL_QUADS);

    // Front face (z = +depth)
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f( width,  height,  depth);
    glVertex3f(-width,  height,  depth);
    glVertex3f(-width, -height,  depth);
    glVertex3f( width, -height,  depth);

    // // Back face (z = -depth)
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f( width, -height, -depth);
    glVertex3f(-width, -height, -depth);
    glVertex3f(-width,  height, -depth);
    glVertex3f( width,  height, -depth);

    // // Left face (x = -width)
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f(-width*2,  height,  depth/2);
    glVertex3f(-width*2,  height, -depth/2);
    glVertex3f(-width*2, -height, -depth/2);
    glVertex3f(-width*2, -height,  depth/2);

    // Right face (x = +width)
    // glColor3f(0.0f,1.0f,0.0f);
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f( width*2,  height, -depth/2);
    glVertex3f( width*2,  height,  depth/2);
    glVertex3f( width*2, -height,  depth/2);
    glVertex3f( width*2, -height, -depth/2);

    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(-width,  height,  depth);
    glVertex3f(-width, -height,  depth);
    glVertex3f( -width*2, -height,  depth/2);
    glVertex3f( -width*2, height, depth/2);

    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(width,  height,  depth);
    glVertex3f(width, -height,  depth);
    glVertex3f( width*2, -height,  depth/2);
    glVertex3f( width*2, height, depth/2);

    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(width,  height,  -depth);
    glVertex3f(width, -height,  -depth);
    glVertex3f( width*2, -height,  -depth/2);
    glVertex3f( width*2, height, -depth/2);

    glColor3f(0.0f,1.0f,0.0f);
    glVertex3f(-width,  -height,  -depth);
    glVertex3f(-width, height,  -depth);
    glVertex3f( -width*2, height,  -depth/2);
    glVertex3f( -width*2, -height, -depth/2);


    // // Top face (y = +height)
    // glColor3f(0.675f, 0.455f, 0.910f); // Purple
    // glVertex3f( width,  height, -depth);
    // glVertex3f(-width,  height, -depth);
    // glVertex3f(-width,  height,  depth);
    // glVertex3f( width,  height,  depth);

    // // Bottom face (y = -height)
    // glColor3f(0.1f, 0.1f, 0.1f); // Dark gray / black
    // glVertex3f( width, -height,  depth);
    // glVertex3f(-width, -height,  depth);
    // glVertex3f(-width, -height, -depth);
    // glVertex3f( width, -height, -depth);

    glEnd();
}


void display() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eyex, eyey, eyez, centerx, centery, centerz, upx, upy, upz);
    drawAxes();
    glPushMatrix();
    glRotatef(rotang, -1, 0, 0); // Yaw for better visibility
    glColor3f(0.2f, 0.7f, 0.9f);
    drawCube();
    drawCube2();

    // glTranslatef(0,5,0);
        glPushMatrix();
        glRotatef(270,1.0f,0.0f,0.0f);
        // drawCylinder(0.2f, 0.2f, 2.05f); 
        
        
        glRotatef(bladeang,0,0,1);
        drawRotatingBlade();

        glPopMatrix();
    glPopMatrix();

    glutSwapBuffers();
}

void keyboardListener(unsigned char key, int x, int y) {
    if (key == 'a') rotang -= 2;
    if (key == 'd') rotang += 2;
    if (key == 'w') bladespeed += 3;
    if (key == 's') 
    {
        if(bladespeed>4)
            bladespeed -= 3;
    }
    glutPostRedisplay();
}
void timerFunc(int value)
{
    bladeang+=bladespeed;
    glutPostRedisplay();
    glutTimerFunc(16,timerFunc,0);
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(640, 640);
    glutCreateWindow("OpenGL Debug Cube");

    initGL();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboardListener);
    glutReshapeFunc(reshapeListener);
    glutTimerFunc(0,timerFunc,0);
    glutMainLoop();
    return 0;
}