#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#include "2005102.h"

const double fovY=80;
const double windowHeight=500;
const double windowWidth=500;

Point lightPosition(-20,-200,70);
Point u(0,0,1);
Point r(1,0,0);
Point l(0,1,0);

vector<Object*>objects;
vector<PointLight>pointLights;
vector<SpotLight>spotLights;

int recursionLevel;
int pixelCount;
int totalObjects;
int pointLightCount;
int spotLightCount;
int imageCount;

double floorWidth=1000; 
double tileWidth=20;

unsigned char* textureData=nullptr;
int textureWidth=0,textureHeight=0,textureChannels=0;
bool useTexture=false;


Color sampleTexture(double u,double v)
{
    if(!textureData||textureWidth<=0||textureHeight<=0)
    {
        return Color(0.5,0.5,0.5); 
    }

    u=std::max(0.0,std::min(1.0,u));
    v=std::max(0.0,std::min(1.0,v));

    int pixelX=(int)(u*(textureWidth-1));
    //Flip Y
    int pixelY=(int)((1.0-v)*(textureHeight-1));
    pixelX=std::max(0,std::min(textureWidth-1,pixelX));
    pixelY=std::max(0,std::min(textureHeight-1,pixelY));

    int index=(pixelY*textureWidth+pixelX)*textureChannels;
    int maxIndex=textureWidth*textureHeight*textureChannels;

    if(index<0||index+2>=maxIndex)
    {
        return Color(1.0,0.0,1.0);
    }
    Color color;
    color.r=textureData[index]/255.0;

    if(textureChannels>=2)
    {
        color.g=textureData[index+1]/255.0;
    } 
    else
    {
        color.g=color.r;
    }

    if(textureChannels>=3)
    {
        color.b=textureData[index+2]/255.0;
    }
    else
    {
        color.b=color.r;
    }
    return color;
}

void loadTexture(const char* filename)
{
    if(textureData)freeTexture();
    textureData=stbi_load(filename,&textureWidth,&textureHeight,&textureChannels,0);
    if(!textureData)
    {
        printf("Failed to load texture: %s\n",filename);
        textureWidth=textureHeight=textureChannels=0;
    }
}

void freeTexture() 
{
    if(textureData)
    {
        stbi_image_free(textureData);
        textureData=nullptr;
        textureWidth=textureHeight=textureChannels=0;
    }
}

void loadData()
{
    ifstream inputFile;
    inputFile.open("scene.txt");
    inputFile>>recursionLevel;
    inputFile>>pixelCount;
    inputFile>>totalObjects;

    for(int i=0;i<totalObjects;i++)
    {
        string objectName;
        inputFile>>objectName;
        if(objectName=="sphere")
        {
            Point center;
            inputFile>>center.x>>center.y>>center.z;
            double radius;
            inputFile>>radius;
            double color[3];
            inputFile>>color[0]>>color[1]>>color[2];
            double coEfficients[4];
            inputFile>>coEfficients[0]>>coEfficients[1]>>coEfficients[2]>>coEfficients[3];
            int shine;
            inputFile>>shine;

            Sphere *sphere=new Sphere(center,radius);

            sphere->setColor(color[0],color[1],color[2]);
            sphere->setCoEfficients(coEfficients[0],coEfficients[1],coEfficients[2],coEfficients[3]);
            sphere->setShine(shine);

            objects.push_back(sphere);
        }
        else if(objectName=="triangle")
        {
            Point A,B,C;
            inputFile>>A.x>>A.y>>A.z;
            inputFile>>B.x>>B.y>>B.z;
            inputFile>>C.x>>C.y>>C.z;

            double color[3];
            inputFile>>color[0]>>color[1]>>color[2];
            double coEfficients[4];
            inputFile>>coEfficients[0]>>coEfficients[1]>>coEfficients[2]>>coEfficients[3];
            int shine;
            inputFile>>shine;

            Triangle *triangle=new Triangle(A,B,C);

            triangle->setColor(color[0],color[1],color[2]);
            triangle->setCoEfficients(coEfficients[0],coEfficients[1],coEfficients[2],coEfficients[3]);
            triangle->setShine(shine);

            objects.push_back(triangle);
        }
        else if(objectName=="general")
        {
            double A,B,C,D,E,F,G,H,I,J;
            inputFile>>A>>B>>C>>D>>E>>F>>G>>H>>I>>J;
            Point referencePoint;
            inputFile>>referencePoint.x>>referencePoint.y>>referencePoint.z;
            double length,width,height;
            inputFile>>length>>width>>height;

            double color[3];
            inputFile>>color[0]>>color[1]>>color[2];
            double coEfficients[4];
            inputFile>>coEfficients[0]>>coEfficients[1]>>coEfficients[2]>>coEfficients[3];
            int shine;
            inputFile>>shine;

            General *general=new General(A,B,C,D,E,F,G,H,I,J,referencePoint,height,width,length);

            general->setColor(color[0],color[1],color[2]);
            general->setCoEfficients(coEfficients[0],coEfficients[1],coEfficients[2],coEfficients[3]);
            general->setShine(shine);

            objects.push_back(general);
        }
    }
    inputFile>>pointLightCount;

    for(int i=0;i<pointLightCount;i++)
    {
        Point lightPosition;
        inputFile>>lightPosition.x>>lightPosition.y>>lightPosition.z;
        double color[3];
        inputFile>>color[0]>>color[1]>>color[2];

        pointLights.push_back(PointLight(lightPosition,color[0],color[1],color[2]));
    }

    inputFile>>spotLightCount;

    for(int i=0;i<spotLightCount;i++)
    {
        Point lightPosition;
        inputFile>>lightPosition.x>>lightPosition.y>>lightPosition.z;
        double color[3];
        inputFile>>color[0]>>color[1]>>color[2];
        Point direction;
        inputFile>>direction.x>>direction.y>>direction.z;
        double cutOffAngle;
        inputFile>>cutOffAngle;

        PointLight pl(lightPosition,color[0],color[1],color[2]);
        spotLights.push_back(SpotLight(pl,direction,cutOffAngle));
    }
    inputFile.close();

    Floor *floor=new Floor(floorWidth,tileWidth);
    objects.push_back(floor);
    loadTexture("Texture3.jpg");
}

void capture()
{
    imageCount++;
    cout<<"Image "<<imageCount<<" Creating..."<<endl;

    //Create bitmap image
    bitmap_image image(pixelCount,pixelCount);

    for(int i=0;i<pixelCount;i++)
    {
        for(int j=0;j<pixelCount;j++)
        {
            image.set_pixel(i,j,0,0,0);
        }
    }

    
    double planeDistance=windowHeight/2.0/(tan(fovY/2.0*pi/180.0));
    Point topLeft=lightPosition+l*planeDistance-r*windowWidth/2+u*windowHeight/2;

    //Height and Width of each pixel
    double du=windowWidth/pixelCount;
    double dv=windowHeight/pixelCount;

    topLeft=topLeft+r*du*0.5-u*dv*0.5;

    for(int i=0;i<pixelCount;i++)
    {
        for(int j=0;j<pixelCount;j++)
        {
            Point curPixel=topLeft+r*du*i-u*dv*j;
            double tMin=INF;

            Point Rd=curPixel-lightPosition;
            Rd=Rd/vecLength(Rd);

            Ray ray(lightPosition,Rd);
            Object *which=NULL;

            for(Object *obj:objects)
            {
                double color[3]={0,0,0};
                double t=obj->intersect(ray,color,0);

                if(t>0&&t<tMin)
                {
                    tMin=t;
                    which=obj;
                }
            }

            if(which!=NULL)
            {
                double color[3]={0,0,0};
                which->intersect(ray,color,1);
                image.set_pixel(i,j,round(color[0]*255),round(color[1]*255),round(color[2]*255));
            }
        }
    }
    cout<<"Image Saved"<<endl;

    string path="Output_1";
    path+=to_string(imageCount);
    path+=".bmp";

    image.save_image(path);
}


void keyboardListener(unsigned char key,int x,int y)
{
    switch(key)
    {
        case '0':
            capture();
            break;
        case 't':
            useTexture=!useTexture;
            printf("Texture %s\n",useTexture?"Enabled":"Disbaled");
            break;
        case '1':
            l=rodriguesRotate(l,u,5);
            r=rodriguesRotate(r,u,5);
            break;
        case '2':
            l=rodriguesRotate(l,u,-5);
            r=rodriguesRotate(r,u,-5);
            break;
        case '3':
            l=rodriguesRotate(l,r,5);
            u=rodriguesRotate(u,r,5);
            break;
        case '4':
            l=rodriguesRotate(l,r,-5);
            u=rodriguesRotate(u,r,-5);
            break;
        case '5':
            u=rodriguesRotate(u,l,5);
            r=rodriguesRotate(r,l,5);
            break;
        case '6':
            u=rodriguesRotate(u,l,-5);
            r=rodriguesRotate(r,l,-5);
            break;

        default:
            break;
    }
}


void specialKeyListener(int key,int x,int y)
{
    switch(key)
    {
        case GLUT_KEY_DOWN:		
            lightPosition=lightPosition-l*5;
            break;
        case GLUT_KEY_UP:	
            lightPosition=lightPosition+l*5;
            break;

        case GLUT_KEY_RIGHT:
            lightPosition=lightPosition+r*5;
            break;
        case GLUT_KEY_LEFT:
            lightPosition=lightPosition-r*5;
            break;

        case GLUT_KEY_PAGE_UP:
            lightPosition=lightPosition+u*5;
            break;
        case GLUT_KEY_PAGE_DOWN:
            lightPosition=lightPosition-u*5;
            break;
        
        default:
            break;
    }
}

void display()
{
    assert(dotProduct(u,r)<1e-9);
    assert(dotProduct(r,l)<1e-9);
    assert(dotProduct(l,u)<1e-9);

    assert(crossProduct(r,l)==u);
    assert(crossProduct(u,r)==l);
    assert(crossProduct(l,u)==r);

    assert(abs(vecLength(u)-1)<1e-9);
    assert(abs(vecLength(r)-1)<1e-9);
    assert(abs(vecLength(l)-1)<1e-9);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);

    glLoadIdentity();

    Point dir=lightPosition+l;

    gluLookAt(lightPosition.x,lightPosition.y,lightPosition.z,dir.x,dir.y,dir.z,u.x,u.y,u.z);

    glMatrixMode(GL_MODELVIEW);

    for(auto u:objects)
        u->draw();
    // for(auto u:pointLights)
    //     u.draw();
    // for(auto u:spotLights)
    //     u.draw();
    
    glutSwapBuffers();
}


void animate()
{
    glutPostRedisplay();
}

void init()
{
    glClearColor(0,0,0,0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(fovY,1,1,1000.0);
}

int main(int argc, char **argv)
{
    glutInit(&argc,argv);
    glutInitWindowSize(windowWidth,windowHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH|GLUT_DOUBLE|GLUT_RGB);

    glutCreateWindow("Ray Tracing");

    init();
    glEnable(GL_DEPTH_TEST);

    glutDisplayFunc(display);
    glutIdleFunc(animate);	

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    loadData();
    glutMainLoop();
    return 0;
}