#include <GL/glut.h>
#include <ctime>
#include <cmath>
#include <string>
#include <iostream>
#include <windows.h>


//In radian unit
double PI=acos(-1);

int screenWidth=500,screenHeight=500;

time_t startTime;

//Create square pointer 
void drawSquarePointer(float size,float red,float green,float blue) 
{
    glColor3f(red,green,blue);
    glBegin(GL_QUADS);
    glVertex3f(-size,-size,0.0f);
    glVertex3f(size,-size,0.0f);
    glVertex3f(size,size,0.0f);
    glVertex3f(-size,size,0.0f);
    glEnd();
}

void drawHand(float angleDeg,float length,float width,float red,float green,float blue) 
{
    glPushMatrix();

    glColor3f(red,green,blue);
    glLineWidth(width);

    //Rotate about -z axis
    glRotatef(angleDeg,0.0f,0.0f,1.0f);

    //Draw line from center to length
    glBegin(GL_LINES);
    glVertex3f(0.0f,0.0f,0.0f);
    glVertex3f(length,0.0f,0.0f);
    glEnd();

    //Move the tracker
    glTranslatef(0.8f,0.0f,0.0f); 
    drawSquarePointer(0.02f,red,green,blue);

    glPopMatrix();
}

void createClock() 
{
    //Creating white border for the clock
    glColor3f(1.0f,1.0f,1.0f);

    //Thickness for the border
    glLineWidth(2);

    //Connecting points
    // glBegin(GL_POINTS);
    glBegin(GL_LINE_LOOP);

    //Placing points in the window(0-360 deg)
    for(int i=0;i<=100;i++) 
    {
        float angle=2*PI*i/100;
        //Circle points will have a radius of 90% of window size
        glVertex2f(0.8f*cos(angle),0.8f*sin(angle));
    }
    glEnd();

    //Creating Hour Markers
    for(int i=0;i<12;i++) 
    {
        float angle=i*PI/6;
        float x1=0.65f*cos(angle);
        float y1=0.65f*sin(angle);
        float x2=0.75f*cos(angle);
        float y2=0.75f*sin(angle);
        glLineWidth(3);
        glBegin(GL_LINES);
        glVertex3f(x1,y1,0.0f);
        glVertex3f(x2,y2,0.0f);
        glEnd();
    }

    //Creating Minute Markers
    for(int i=0;i<60;i++)
    {
        //To avoid conflict with hour markers
        if(i%5==0) continue;

        float angle=i*PI/30;
        float x1=0.70f*cos(angle);
        float y1=0.70f*sin(angle);
        float x2=0.75f*cos(angle);
        float y2=0.75f*sin(angle);
        glLineWidth(1);
        glBegin(GL_LINES);
        glVertex3f(x1,y1,0.0f);
        glVertex3f(x2,y2,0.0f);
        glEnd();
    }

    //Draw numbers 1–12
    for(int i=1;i<=12;i++) 
    {
        float angle=(PI/6)*(3-i);
        float x=0.50f*cos(angle);
        float y =0.50f*sin(angle);
        std::string num=std::to_string(i);
        // std::cout << num << std::endl;

        //Placing the number below the hour markers
        glRasterPos2f(x-0.04f,y-0.04f);
        for(char c:num) 
        {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,c);
        }
    }
}


void display() 
{
    //Clear the color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //Reset transformations
    glLoadIdentity();
    
    //Get elapsed milliseconds since program start
    int elapsed_ms=glutGet(GLUT_ELAPSED_TIME);
    float elapsed_seconds=elapsed_ms/1000.0f;

    //Calculate current simulated time based on startTime
    time_t currentTime=startTime+(time_t)elapsed_seconds;
    struct tm* local=localtime(&currentTime);

    int hours=local->tm_hour%12;
    int minutes=local->tm_min;
    int seconds=local->tm_sec;

    float fractional_seconds=elapsed_seconds-(int)elapsed_seconds;

    float secAngle=90.0f-(seconds+fractional_seconds)*6.0f;
    float minAngle=90.0f-(minutes+(seconds+fractional_seconds)/60.0f)*6.0f;
    float hourAngle=90.0f-(hours+(minutes+seconds/60.0f)/60.0f)*30.0f;


    createClock();

    //Hour Hand-White
    drawHand(hourAngle,0.35f,4.0f,1.0f,1.0f,1.0f);

    //Minute Hand-White
    drawHand(minAngle,0.55f,3.0f,0.0f,0.0f,1.0f);

    //Second Hand-Red
    drawHand(secAngle,0.65f,2.0f,1.0f,0.0f,0.0f);
    glutSwapBuffers();
}


void timer(int value) 
{
    glutPostRedisplay();
    glutTimerFunc(100, timer, 0); 
}


void init() 
{
    //Black Background
    glClearColor(0.0f,0.0f,0.0f,1.0f);
}


int main(int argc,char** argv) 
{
    startTime=time(0);
    //Initialize GLUT
    glutInit(&argc, argv);

    //Set initial window size and position
    glutInitWindowSize(screenWidth,screenHeight);
    glutInitWindowPosition(100,100);


    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    
    //Create a window with the title
    glutCreateWindow("Analog Clock");

    glutDisplayFunc(display);

    glutTimerFunc(0,timer,0);

    init();

    glutMainLoop();
    return 0;
}