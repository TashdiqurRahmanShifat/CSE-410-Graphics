#include <GL/glut.h>
#include <math.h>
#include <iostream>
#include <random>
#define ROTATION_ANGLE 0.01f

//Camera position coordinates
GLfloat eyex=7.8,eyey=5.5,eyez=7.8; 
//Look-at point coordinates
GLfloat centerx=0,centery=-1,centerz=0;
//Up vector coordinates
GLfloat upx=0,upy=1,upz=0;             
//Rotation angles
GLfloat yaw=0,pitch=0,roll=0;

//Movement speed
GLfloat moveSpeed=0.05f;


bool simulate=false;
bool showVelocity=false;
float restitution=0.75f;
float initialSpeed=1.0f;


float spinAxis[3]={0.0f, 0.0f, 0.0f};
float spinAngle=0.0f;


//Random Number Generator
std::random_device rd;
std::mt19937 gen(rd());
float randomBetween(float low,float high) 
{
    std::uniform_real_distribution<> dis(low,high);
    return dis(gen);
}



struct Ball 
{
    float radius = 0.1f;
    float position[3]={0.0f,-1.4f,0.0f}; 
    float velocity[3];// = {1.0f, .40f, 2.0f};
    float rotation[3]={0.0f,0.0f,0.0f};
} ball;


void showMessage()
{
    std::cout<<"--- 3D Bouncing Ball Controls ---\n"
            <<"Space: Toggle simulation on/off\n"
            <<"v: Toggle velocity arrow display\n"
            <<"+/-: Increase/decrease ball speed (when paused)\n"
            <<"r: Reset ball to random position (when paused)\n"
            <<"\n"
            <<"--- Camera Controls ---\n"
            <<"Arrow keys: Move camera left/right/forward/backward\n"
            <<"Page Up/Down: Move camera up/down\n"
            <<"w/s: Move camera up/down (keeping reference point)\n"
            <<"1/2: Look left/right\n"
            <<"3/4: Look up/down\n"
            <<"5/6: Tilt camera clockwise/counterclockwise\n"
            <<"Reset ball with speed: "<<ball.velocity[1]<<".00f"<< std::endl;
}



void initVelocity() 
{
    ball.velocity[0]=randomBetween(-2.0f,2.0f);
    ball.velocity[1]=initialSpeed;
    ball.velocity[2]=randomBetween(1.0f,2.0f);
    ball.position[0]=randomBetween(-2.9f, 2.9f); 
    ball.position[2]=randomBetween(-2.9f, 2.9f); 
}



void initGL()
{
    glClearColor(0.2,0.2,0.2,1.0);
    glEnable(GL_DEPTH_TEST); 
}

void drawChessboardFace() 
{
    glPushMatrix();
    glLoadIdentity();

    int rows=20;
    int cols=20;
    float startX=-3.0f;
    float startZ=-3.0f;
    float tileWidth=(6.0f)/cols;//(3-(-3))=6
    float tileDepth=(6.0f)/rows;//(3-(-3))=6
    float y=-1.5f;

    for(int i=0;i<rows;++i) 
    {
        for(int j=0;j<cols;++j) 
        {
            if((i+j)%2==0)
                glColor3f(1.0f,1.0f,1.0f);//White tile
            else
                glColor3f(0.0f,0.0f,0.0f);//Black tile

            float x1=startX+j*tileWidth;
            float x2=x1+tileWidth;
            float z1=startZ+i*tileDepth;
            float z2=z1+tileDepth;

            glBegin(GL_QUADS);
                glVertex3f(x1,y,z1);
                glVertex3f(x2,y,z1);
                glVertex3f(x2,y,z2);
                glVertex3f(x1,y,z2);
            glEnd();
        }
    }
    glPopMatrix();
}

void drawCube()
{
    glBegin(GL_QUADS);

    //Left face
    glColor3f(0.392f,0.949f,0.949f);
    glVertex3f(-3.0f,-1.5f,3.0f);  
    glVertex3f(3.0f,-1.5f,3.0f);   
    glVertex3f(3.0f,1.5f,3.0f);    
    glVertex3f(-3.0f,1.5f,3.0f);

    //Right-Right face
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f(-3.0f,-1.5f,-3.0f);
    glVertex3f(-3.0f,1.5f,-3.0f);
    glVertex3f(3.0f,1.5f,-3.0f);
    glVertex3f(3.0f,-1.5f,-3.0f);

    //Top face
    glColor3f(0.675f,0.455f,0.910f);
    glVertex3f(-3.0f,1.5f,-3.0f);
    glVertex3f(-3.0f,1.5f,3.0f);
    glVertex3f(3.0f,1.5f,3.0f);
    glVertex3f(3.0f,1.5f,-3.0f);

    //Right face
    glColor3f(0.408f,0.859f,0.298f);
    glVertex3f(3.0f,-1.5f,-3.0f);
    glVertex3f(3.0f,1.5f,-3.0f);
    glVertex3f(3.0f,1.5f,3.0f);
    glVertex3f(3.0f,-1.5f,3.0f);

    //Left-Left face
    glColor3f(0.941f,0.322f,0.322f);
    glVertex3f(-3.0f,-1.5f,-3.0f);
    glVertex3f(-3.0f,-1.5f,3.0f);
    glVertex3f(-3.0f,1.5f,3.0f);
    glVertex3f(-3.0f,1.5f,-3.0f);

    drawChessboardFace();

    glEnd();
}



void crossProduct(float a[3],float b[3],float result[3])
{
    result[0]=a[1]*b[2]-a[2]*b[1];
    result[1]=a[2]*b[0]-a[0]*b[2];
    result[2]=a[0]*b[1]-a[1]*b[0];
}


//For storing previous position of the ball
float previousPosition[3]={0.0f,0.0f,0.0f};


void updateBall(float deltaTime) 
{
    if(!simulate) return;

    //Store previous position
    previousPosition[0]=ball.position[0];
    previousPosition[1]=ball.position[1];
    previousPosition[2]=ball.position[2];

    float* p=ball.position;
    float* v=ball.velocity;
    float gravity=-9.8f;

    //Apply gravity
    v[1]+=gravity*deltaTime*0.1f;

    //Update position
    for(int i=0;i<3;++i)
        p[i]+=v[i]*deltaTime*0.2f;

    //Bounce Handling
    for (int i=0;i<3;++i) 
    {
        float bound=(i==1)?1.5f:3.0f;//Y-axis bound is 1.5f,others are 3.0

        if(p[i]+ball.radius>bound) 
        {
            p[i]=bound-ball.radius;
            v[i]=-v[i]*restitution;
        } 
        else if(p[i]-ball.radius<-bound) 
        {
            p[i]=-bound+ball.radius;
            v[i]=-v[i]*restitution;
        }
    }


    if (fabs(v[1])<0.01f && p[1]-ball.radius<=-1.5f) 
    {
        v[1]=0;
    }

    //Ball Movement Vector
    float movement[3]={
        ball.position[0]-previousPosition[0],
        ball.position[1]-previousPosition[1],
        ball.position[2]-previousPosition[2]
    };

    //Cross product with up vector to get rotation axis
    float up[3]= {0.0f,1.0f,0.0f};
    crossProduct(up,movement,spinAxis);

    //Normalize axis
    float len=sqrt(spinAxis[0]*spinAxis[0]+spinAxis[1]*spinAxis[1]+spinAxis[2]*spinAxis[2]);
    if(len>0.001f) 
    {
        spinAxis[0]/=len;
        spinAxis[1]/=len;
        spinAxis[2]/=len;

        //Calculate distance ball has covered
        float dist=sqrt(movement[0]*movement[0]+movement[1]*movement[1]+movement[2]*movement[2]);
        //Set spin rate.Controls spin speed
        spinAngle+=dist*500.0f; 
    }
}



void drawBall() 
{
    glPushMatrix();
    glTranslatef(ball.position[0],ball.position[1],ball.position[2]);


    if(spinAngle>0.01f) 
    {
        glRotatef(spinAngle,spinAxis[0],spinAxis[1],spinAxis[2]);
    }
    

    const int slices=24;
    const int stacks=20;
    const float radius=ball.radius;

    
    for(int i=0;i<slices;++i) 
    {
        float theta1=i*2*M_PI/slices;
        float theta2=(i+1)*2*M_PI/slices;

        //Upper half pattern.Alternating red and green
        if(i%2==0)
            glColor3f(1.0f,0.0f,0.0f); 
        else
            glColor3f(0.0f,1.0f,0.0f);

        glBegin(GL_QUAD_STRIP);

            for(int j=0;j<=stacks/2;++j) 
            { 
                float phi=j*M_PI/stacks;
                float x1=radius*sin(phi)*cos(theta1);
                float y1=radius*cos(phi);
                float z1=radius*sin(phi)*sin(theta1);

                float x2=radius*sin(phi)*cos(theta2);
                float y2=radius*cos(phi);
                float z2=radius*sin(phi)*sin(theta2);

                glVertex3f(x1,y1,z1);
                glVertex3f(x2,y2,z2);
            }

        glEnd();

        //Lower half pattern.Alternating green and red
        if(i%2==0)
            glColor3f(0.0f,1.0f,0.0f); 
        else
            glColor3f(1.0f,0.0f,0.0f); 

        glBegin(GL_QUAD_STRIP);

            for(int j=stacks/2;j<=stacks;++j) 
            {
                float phi=j* M_PI/stacks;
                float x1=radius*sin(phi)*cos(theta1);
                float y1=radius*cos(phi);
                float z1=radius*sin(phi)*sin(theta1);

                float x2=radius*sin(phi)*cos(theta2);
                float y2=radius*cos(phi);
                float z2=radius*sin(phi)*sin(theta2);

                glVertex3f(x1, y1, z1);
                glVertex3f(x2, y2, z2);
            }

        glEnd();
    }
    
    glPopMatrix();


    if(showVelocity) 
    {
        float vx=ball.velocity[0];
        float vy=ball.velocity[1];
        float vz=ball.velocity[2];
        float vmag=sqrt(vx*vx+vy*vy+vz*vz);
    
        //Arrow length scaling
        float arrowScale=0.3f;  
        if(vmag>0.001f) 
        {
            //Normalize direction
            float dx=vx/vmag;
            float dy=vy/vmag;
            float dz=vz/vmag;
    

            dx*=arrowScale;
            dy*=arrowScale;
            dz*=arrowScale;
    

            float px=ball.position[0];
            float py=ball.position[1];
            float pz=ball.position[2];
    
            //Arrow length ends point
            float tipX=px+dx;
            float tipY=py+dy;
            float tipZ=pz+dz;
    

            glLineWidth(2);
            glBegin(GL_LINES);
                glColor3f(1.0f,1.0f,0.0f);
                glVertex3f(px,py,pz);
                glVertex3f(tipX,tipY,tipZ);
            glEnd();
    
            //Triangle head pointing in direction
            float headLength=0.05f;
            float baseOffset=0.02f;
    
            //Calculate vector which is perpendicular of velocity vector(dx,dy,dz) an up vector
            float ux=-dz;//dy*0-dz*1
            float uy=0;//dz*0-dx*0
            float uz=dx;//dx*1-dy*0

            float len=sqrt(ux*ux+uz*uz);
            if(len>0.001f) 
            {
                ux/=len;
                uz/=len;
            }
    
            //Triangle base corners
            float bx1=tipX-dx*headLength+ux*baseOffset;
            float by1=tipY-dy*headLength;
            float bz1=tipZ-dz*headLength+uz*baseOffset;
    
            float bx2=tipX-dx*headLength-ux*baseOffset;
            float by2=tipY-dy*headLength;
            float bz2=tipZ-dz*headLength-uz*baseOffset;
    

            glBegin(GL_TRIANGLES);
            
            glColor3f(1.0f, 0.647f, 0.0f);
                glVertex3f(tipX,tipY,tipZ);
                glVertex3f(bx1,by1,bz1);
                glVertex3f(bx2,by2,bz2);
            glEnd();
        }
    }
    
    
}


void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW); 
    glLoadIdentity();

    // glRotatef(pitch,1.0f,0.0f,0.0f);  
    // glRotatef(yaw,0.0f,1.0f,0.0f);    
    // glRotatef(roll,0.0f,0.0f,1.0f);   

    gluLookAt(eyex,eyey,eyez,centerx,centery,centerz,upx,upy,upz);


    drawCube();
    drawBall(); 

    glutSwapBuffers();
}


void timer(int value) 
{
    updateBall(0.016f); 
    glutPostRedisplay();
    glutTimerFunc(16,timer,0);
}



void reshapeListener(GLsizei width, GLsizei height)
{
    if(height==0)
        height=1;

    GLfloat aspect=(GLfloat)width/(GLfloat)height;
    glViewport(0,0,width,height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f,aspect,0.1f,100.0f);
}

 
void rotateView(float angle,float x,float y,float z) 
{
    float viewX=centerx-eyex;
    float viewY=centery-eyey;
    float viewZ=centerz-eyez;

    //Rodrigues' Rotation Formula
    //v'=v*cos(θ)+(k×v)*sin(θ)+k*(k⋅v)*(1-cos(θ))
 
    float c=cos(angle);
    float s=sin(angle);

    float dot=x*viewX+y*viewY+z*viewZ;
    float newViewX=viewX*c+(y*viewZ-z*viewY)*s+x*dot*(1-c);
    float newViewY=viewY*c+(z*viewX-x*viewZ)*s+y*dot*(1-c);
    float newViewZ=viewZ*c+(x*viewY-y*viewX)*s+z*dot*(1-c);

    centerx=eyex+newViewX;
    centery=eyey+newViewY;
    centerz=eyez+newViewZ;
}

void rotateUp(float angle,float x,float y,float z) 
{
    float c=cos(angle);
    float s=sin(angle);

    float dot=x*upx+y*upy+z*upz;
    float newUpX=upx*c+(y*upz-z*upy)*s+x*dot*(1-c);
    float newUpY=upy*c+(z*upx-x*upz)*s+y*dot*(1-c);
    float newUpZ=upz*c+(x*upy-y*upx)*s+z*dot*(1-c);

    upx=newUpX;
    upy=newUpY;
    upz=newUpZ;
}

bool track=false;

void keyboardListener(unsigned char key,int x,int y) 
{
    switch (key) 
    {
        //Rotate left
        case '1':
            rotateView(ROTATION_ANGLE,0.0f,1.0f,0.0f);
            break;

        //Rotate right
        case '2': 
            rotateView(-ROTATION_ANGLE,0.0f,1.0f,0.0f);
            break;



        //Look up
        case '3': 
        {
            float viewX=centerx-eyex;
            float viewY=centery-eyey;
            float viewZ=centerz-eyez;

            //Cross product of view and up
            float rightX=viewY*upz-viewZ*upy;
            float rightY=viewZ*upx-viewX*upz;
            float rightZ=viewX*upy-viewY*upx;

            //Normalize the right vector
            float length=sqrt(rightX*rightX+rightY*rightY+rightZ*rightZ);
            rightX/=length;
            rightY/=length;
            rightZ/=length;

            //Rotate view and up vectors around the right vector
            rotateView(ROTATION_ANGLE,rightX,rightY,rightZ); 
            rotateUp(ROTATION_ANGLE,rightX,rightY,rightZ);   
            break;
        }

        //Look down
        case '4': 
        {
            float viewX=centerx-eyex;
            float viewY=centery-eyey;
            float viewZ=centerz-eyez;

            //Cross product of view and up
            float rightX=viewY*upz-viewZ*upy;
            float rightY=viewZ*upx-viewX*upz;
            float rightZ=viewX*upy-viewY*upx;

            //Normalize the right vector
            float length=sqrt(rightX*rightX+rightY*rightY+rightZ*rightZ);
            rightX/=length;
            rightY/=length;
            rightZ/=length;

            //Rotate view and up vectors around the right vector
            rotateView(-ROTATION_ANGLE,rightX,rightY,rightZ); 
            rotateUp(-ROTATION_ANGLE,rightX,rightY,rightZ); 
            break;
        }

        //Roll right
        case '5': 
            rotateUp(ROTATION_ANGLE,(centerx-eyex),(centery-eyey),(centerz-eyez));
            break;

        //Roll left
        case '6':
            rotateUp(-ROTATION_ANGLE,(centerx-eyex),(centery-eyey),(centerz-eyez));
            break;
        
        case 'w':
            eyey+=moveSpeed;
            break;
        case 's':
            eyey-=moveSpeed;
            break;

        case 27:
            exit(0);
            break;


        case ' ':
            simulate=!simulate;
            
            if (simulate && !track) 
            {
                ball.position[1]=-1.2f;//Set to a new height when simulation starts
            }
            track=true;

            break;

        case 'r':
            if(!simulate)
            {
                ball.position[0]=randomBetween(-2.9f, 2.9f); 
                ball.position[2]=randomBetween(-2.9f, 2.9f);  
                ball.position[1]=-1.4f;//The ball is placed just above the floor
            
                initVelocity();
                spinAngle=0;
                spinAxis[0]=spinAxis[1]=spinAxis[2]=0;
                track=false;
                std::cout<<"Reset ball with speed: "<<ball.velocity[1]<<".00f"<<std::endl;
            }
            break;

        case '+':
            if(!track)
            {
                initialSpeed++;
                ball.velocity[1]=initialSpeed;
                std::cout<<"Unknown special key pressed"<<std::endl;
                std::cout<<"Initial speed increased to: "<<ball.velocity[1]<<".00f"<<std::endl;
               
            }
            break;

        case '-':
            if(!track && ball.velocity[1]>1.0f)
            {
                initialSpeed--;
                ball.velocity[1]=initialSpeed;
                std::cout<<"Initial speed decreased to: "<<ball.velocity[1]<<".00f"<<std::endl;
            }
            break;

        case 'v':
            showVelocity=!showVelocity; 
            break;
        
        default:
            std::cout<<"Unknown key pressed"<<std::endl;
            break;
    }

    glutPostRedisplay();
}



void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
        case GLUT_KEY_UP:
        {
            float dirX=centerx-eyex;
            float dirY=centery-eyey;
            float dirZ=centerz-eyez;
        
            float length=sqrt(dirX*dirX+dirY*dirY+dirZ*dirZ);
        
            if(length>0.0001f) 
            {
                dirX/=length;
                dirY/=length;
                dirZ/=length;
        
                eyex+=dirX*moveSpeed;
                eyey+=dirY*moveSpeed;
                eyez+=dirZ*moveSpeed;
        
                centerx+=dirX*moveSpeed;
                centery+=dirY*moveSpeed;
                centerz+=dirZ*moveSpeed;
            }
            break;
        }
        
        case GLUT_KEY_DOWN:
        {
            float dirX=centerx-eyex;
            float dirY=centery-eyey;
            float dirZ=centerz-eyez;
        
            float length=sqrt(dirX*dirX+dirY*dirY+dirZ*dirZ);
        
            if(length>0.0001f) 
            {
                dirX/=length;
                dirY/=length;
                dirZ/=length;
        
                eyex-=dirX*moveSpeed;
                eyey-=dirY*moveSpeed;
                eyez-=dirZ*moveSpeed;
        
                centerx-=dirX*moveSpeed;
                centery-=dirY*moveSpeed;
                centerz-=dirZ*moveSpeed;
            }
            break;
        }
        


        case GLUT_KEY_LEFT:
        {
            //Calculate view direction
            float viewX=centerx-eyex;
            float viewY=centery-eyey;
            float viewZ=centerz-eyez;
        
            //Calculate cross product(view×up)
            float rightX=viewY*upz-viewZ*upy;
            float rightY=viewZ*upx-viewX*upz;
            float rightZ=viewX*upy-viewY*upx;
        
            //Normalize right vector
            float length=sqrt(rightX*rightX+rightY*rightY+rightZ*rightZ);
            rightX/=length;
            rightY/=length;
            rightZ/=length;
        
            //Move camera to the left
            eyex-=rightX*moveSpeed;
            eyey-=rightY*moveSpeed;
            eyez-=rightZ*moveSpeed;
        
            //Move center
            centerx-=rightX*moveSpeed;
            centery-=rightY*moveSpeed;
            centerz-=rightZ*moveSpeed;
            break;
        }
        
        case GLUT_KEY_RIGHT: 
        {
            //Calculate view direction
            float viewX=centerx-eyex;
            float viewY=centery-eyey;
            float viewZ=centerz-eyez;
        
            //Calculate cross product(view×up)
            float rightX=viewY*upz-viewZ*upy;
            float rightY=viewZ*upx-viewX*upz;
            float rightZ=viewX*upy-viewY*upx;
        
            //Normalize right vector
            float length=sqrt(rightX*rightX+rightY*rightY+rightZ*rightZ);
            rightX/=length;
            rightY/=length;
            rightZ/=length;
        
            //Move camera to the right
            eyex+=rightX*moveSpeed;
            eyey+=rightY*moveSpeed;
            eyez+=rightZ*moveSpeed;

            //Move center
            centerx+=rightX*moveSpeed;
            centery+=rightY*moveSpeed;
            centerz+=rightZ*moveSpeed;
            break;
            
        }  
            

        //Move upward
        case GLUT_KEY_PAGE_UP:
            eyey+=moveSpeed;
            centery+=moveSpeed;
            break;

        //Move downward
        case GLUT_KEY_PAGE_DOWN: 
            eyey-=moveSpeed;
            centery-=moveSpeed;
            break;
  
        default:
            std::cout<<"Unknown key pressed"<<std::endl;
            break;
    }

    glutPostRedisplay(); 
}


int main(int argc,char **argv)
{
    srand(time(0));
    //Initialize GLUT
    glutInit(&argc,argv);
    
    initVelocity();
    showMessage();

    //Configure display mode and window
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(640,640);
    glutInitWindowPosition(50,50);
    glutCreateWindow("Bouncing Ball");

    //Register callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    //Initialize OpenGL settings
    initGL();
    glutTimerFunc(0, timer, 0);

    glutMainLoop();

    return 0;
}