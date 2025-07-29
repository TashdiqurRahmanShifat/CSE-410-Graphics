#ifndef RAYTRACER_2005102_H
#define RAYTRACER_2005102_H

#include<bits/stdc++.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"
#define EPSILON 1e-6
#define INF 1e9

using namespace std;


const double pi=acos(-1);

struct Color
{
    double r,g,b;
    Color()
    {
        this->r=0;
        this->g=0;
        this->b=0;
    }

    Color(double r,double g,double b)
    {
        this->r=r;
        this->g=g;
        this->b=b;
    }
};

extern unsigned char* textureData;
extern int textureWidth,textureHeight,textureChannels;
extern bool useTexture;

Color sampleTexture(double u,double v);
void loadTexture(const char* filename);
void freeTexture();


extern vector<class PointLight>pointLights;
extern vector<class SpotLight>spotLights;
extern vector<class Object*>objects;
extern int recursionLevel;


class Point
{
public:
    double x,y,z;

    Point(){}
    Point(double x,double y,double z)
    {
        this->x=x;
        this->y=y;
        this->z=z;
    }

    Point operator+(const Point& secondPoint)
    {
        return Point(x+secondPoint.x,y+secondPoint.y,z+secondPoint.z);
    }

    Point operator-(const Point& secondPoint)
    {
        return Point(x-secondPoint.x,y-secondPoint.y,z-secondPoint.z);
    }

    //Check for speed
    bool operator==(const Point& secondPoint)
    {
        return abs(x-secondPoint.x)<1e-9 && abs(y-secondPoint.y)<1e-9 && abs(z-secondPoint.z)<1e-9;
    }

    Point operator*(double secondPoint)
    {
        return Point(x*secondPoint,y*secondPoint,z*secondPoint);
    }   

    Point operator/(double secondPoint)
    {
        return Point(x/secondPoint,y/secondPoint,z/secondPoint);
    }
};

bool isSamePoint(Point A,Point B)
{
    return fabs(A.x-B.x)<EPSILON && fabs(A.y-B.y)<EPSILON && fabs(A.z-B.z)<EPSILON;
}

double dotProduct(Point a,Point b)
{
    return a.x*b.x+a.y*b.y+a.z*b.z;
}

Point crossProduct(Point a,Point b)
{
    return Point(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);
}

double vecLength(Point a)
{
    return sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
}

//k-axis,v-vector to rotate
Point rodriguesRotate(Point v,Point k,double angle)
{
    k=k/vecLength(k);
    angle=angle*pi/180.0;
    return v*cos(angle)+crossProduct(k,v)*sin(angle)+k*dotProduct(k,v)*(1-cos(angle));
}

bool isSame(double a,double b)
{
    return fabs(a-b)<EPSILON;
}


double angleInDegree(Point A,Point B)
{
    double angle=acos(dotProduct(A,B)/(vecLength(A)*vecLength(B)))*180.0/pi;
    return angle;
}

class Matrix
{
public:
    double matrix[3][3];

    //Converting from 1D array to 2D matrix
    Matrix(double *mat) 
    {
        for(int i=0;i<3;i++)
        {
            for(int j=0;j<3;j++)
            {
                matrix[i][j]=mat[i*3+j];
            }
        }
    }

    double determinant()
    {
        return matrix[0][0]*(matrix[1][1]*matrix[2][2]-matrix[1][2]*matrix[2][1])
        -matrix[0][1]*(matrix[1][0]*matrix[2][2]-matrix[1][2]*matrix[2][0])
        +matrix[0][2]*(matrix[1][0]*matrix[2][1]-matrix[1][1]*matrix[2][0]);
    }
};

class Ray
{
public:
    Point start;
    Point direction;
    Ray(Point start,Point direction)
    {
        this->start=start;
        this->direction=direction;

        //Normalize the direction vector to ensure that,in start+direction*t equation,
        //t directly represents the distance from the origin
        this->direction=this->direction/vecLength(this->direction);
    }
};

//Object class
class Object
{
public:
    Point referencePoint;
    double height,width,length;
    double color[3];
    double coEfficients[4];

    int shine;

    Object(){}

    virtual void draw()=0;
    virtual double intersect(Ray ray,double* color,int level)
    {
        return -1;
    }

    void setColor(double r,double g,double b) 
    {
        color[0]=r;
        color[1]=g;
        color[2]=b;
    }
    
    void setShine(int shine) 
    {
        this->shine=shine;
    }

    void setCoEfficients(double ambient,double diffuse,double specular,double recursiveReflection)
    {
        coEfficients[0]=ambient;
        coEfficients[1]=diffuse;
        coEfficients[2]=specular;
        coEfficients[3]=recursiveReflection;
    }
};


class PointLight
{
public:
    Point lightPosition;
    double color[3];

    PointLight(){}

    PointLight(Point lightPosition,double r,double g,double b)
    {
        this->lightPosition=lightPosition;
        color[0]=r;
        color[1]=g;
        color[2]=b;
    }
    
    void draw()
    {
        const int stacks=45;
        const int slices=45;

        //radius of the light sphere
        const double radius=1.0;

        //To store the 3D coordinates of the sphere's surface points
        Point points[100][100];
        int i,j;
        double h,r;

        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }

        glPushMatrix();
        glTranslatef(lightPosition.x,lightPosition.y,lightPosition.z);

        //Draw using points
        for(i=0;i<stacks;i++)
        {
            glColor3f(color[0],color[1],color[2]);
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);
                {
                    //Upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                    //Lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }
                glEnd();
            }
        }
        glPopMatrix();
    }
};


class SpotLight
{
public:
    PointLight pointLight;
    Point lightDirection;
    double cutoffAngle;

    SpotLight(){}

    SpotLight(PointLight pointLight,Point lightDirection,double cutoffAngle)
    {
        this->pointLight=pointLight;
        this->lightDirection=lightDirection;
        this->cutoffAngle=cutoffAngle;
    }

    void draw()
    {
        this->pointLight.draw();
    }
};



class Sphere:public Object
{
public:
    Point center;
    double radius;

    Sphere(Point center, double radius)
    {
        this->center=center;
        this->radius=radius;
    }
    
    void draw()
    {
        const int stacks=45;
        const int slices=45;

        Point points[100][100];
        int i,j;
        double h,r;

        //generate points
        for(i=0;i<=stacks;i++)
        {
            h=radius*sin(((double)i/(double)stacks)*(pi/2));
            r=radius*cos(((double)i/(double)stacks)*(pi/2));
            for(j=0;j<=slices;j++)
            {
                points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
                points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
                points[i][j].z=h;
            }
        }

        glPushMatrix();
        glTranslatef(center.x,center.y,center.z);

        //Draw using points
        for(i=0;i<stacks;i++)
        {
            glColor3f(color[0],color[1],color[2]);
            for(j=0;j<slices;j++)
            {
                glBegin(GL_QUADS);
                {
                    //Upper hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);

                    //Lower hemisphere
                    glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
                    glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
                    glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
                    glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
                }
                glEnd();
            }
        }
        glPopMatrix();
    }

    double intersect(Ray ray,double *color,int level)
    {
        //origin implicitly moved to center, has no impact on value of parameter t
        ray.start=ray.start-center;

        double a=1;//d is normalized,so d^2=1
        double b=dotProduct(ray.start,ray.direction)*2;
        double c=dotProduct(ray.start,ray.start)-radius*radius;

        double d=b*b-4*a*c;
        double t;

        if(d<0)
        {
            t=-1;
        }
        else if(isSame(d,0))
        {
            t=-b/(a*2);
        }
        else
        {
            d=sqrt(d);
            double t_minus=(-b-d)/(a*2);
            double t_plus=(-b+d)/(a*2);

            if(t_minus>0)
            {
                t=t_minus;
            }
            else if(t_plus>0)
            {
                t=t_plus;
            }
            else
            {
                t=-1;
            }
        }
        if(level==0||t<0)
            return t;
        ray.start=ray.start+center;
        Point intersectionPoint=ray.start+ray.direction*t;

        double intersectionPointColor[3]= {0,0,0};

        intersectionPointColor[0]=this->color[0];
        intersectionPointColor[1]=this->color[1];
        intersectionPointColor[2]=this->color[2];

        //Ambient color
        color[0]=intersectionPointColor[0]*coEfficients[0];
        color[1]=intersectionPointColor[1]*coEfficients[0];
        color[2]=intersectionPointColor[2]*coEfficients[0];

        //Lighting calculations
        Point normal=intersectionPoint-center;
        normal=normal/vecLength(normal);

        if(dotProduct(normal,ray.direction*-1)<0)
        {
            normal=normal*-1;
        }

        for(PointLight pl:pointLights)
        {
            Ray L(pl.lightPosition,intersectionPoint-pl.lightPosition);

            double tBlock=INF;
            for(auto obj:objects)
            {
                double tempColor[3]={0,0,0};

                //Holds the nearest intersection point with the object
                double tempIntersect=obj->intersect(L,tempColor,0);
                if(tempIntersect>0&&tempIntersect<tBlock)
                {
                    //Store the nearest intersection point
                    tBlock=tempIntersect;
                }
            }

            if(tBlock>0)
            {
                Point blocked_at=L.start+L.direction*tBlock;
                if(!isSamePoint(blocked_at,intersectionPoint))
                {
                    continue;
                }
            }

            Point directionToLight=L.direction*-1;

            double lambert=dotProduct(directionToLight,normal);
            Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
            R=R/vecLength(R);
            lambert=max(0.0,lambert);

            double phong=dotProduct(ray.direction*-1,R);
            phong=max(0.0,phong);

            //Diffuse
            color[0]+=pl.color[0]*coEfficients[1]*lambert*intersectionPointColor[0];
            color[1]+=pl.color[1]*coEfficients[1]*lambert*intersectionPointColor[1];
            color[2]+=pl.color[2]*coEfficients[1]*lambert*intersectionPointColor[2];


            //Specular
            color[0]+=pl.color[0]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[0];
            color[1]+=pl.color[1]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[1];
            color[2]+=pl.color[2]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[2];
        }

        for(SpotLight sl:spotLights)
        {
            Ray L(sl.pointLight.lightPosition,intersectionPoint-sl.pointLight.lightPosition);
            
            //Check cutoff angle
            double angle=angleInDegree(sl.lightDirection,L.direction);
            if(angle>sl.cutoffAngle)
            {
                continue;
            }
            double tBlock=INF;
            for(auto obj:objects)
            {
                double tempColor[3]={0,0,0};
                double tempIntersect=obj->intersect(L,tempColor,0);
                if(tempIntersect>0&&tempIntersect<tBlock)
                {
                    tBlock=tempIntersect;
                }
            }

            if(tBlock>0)
            {
                Point blocked_at=L.start+L.direction*tBlock;
                if(!isSamePoint(blocked_at,intersectionPoint))
                {
                    continue;
                }
            }

            Point directionToLight=L.direction*-1;
            
            double lambert=dotProduct(directionToLight,normal);
            Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
            R=R/vecLength(R);
            lambert=max(0.0,lambert);

            double phong=dotProduct(ray.direction*-1,R);
            phong=max(0.0,phong);


            color[0]+=sl.pointLight.color[0]*coEfficients[1]*lambert*intersectionPointColor[0];
            color[1]+=sl.pointLight.color[1]*coEfficients[1]*lambert*intersectionPointColor[1];
            color[2]+=sl.pointLight.color[2]*coEfficients[1]*lambert*intersectionPointColor[2];

            color[0]+=sl.pointLight.color[0]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[0];
            color[1]+=sl.pointLight.color[1]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[1];
            color[2]+=sl.pointLight.color[2]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[2];
        }

        if(level==recursionLevel)
            return t;

        Point directionToLight=ray.direction*-1;
        Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;

        Ray reflectRay(intersectionPoint+R*EPSILON*3,R);
        Object *closestObj=NULL;

        //Minimum distance to the closest object
        double tMinimumDistance=INF;

        for(Object *obj:objects)
        {
            double tempColor[3]={0,0,0};
            double tempIntersect=obj->intersect(reflectRay,tempColor,0);

            if(tempIntersect>0&&tempIntersect<tMinimumDistance)
            {
                tMinimumDistance=tempIntersect;
                closestObj=obj;
            }
        }

        if(closestObj!=NULL)
        {
            double reflectedColor[3]={0,0,0};
            closestObj->intersect(reflectRay,reflectedColor,level+1);

            //Recursive reflection color
            color[0]+=reflectedColor[0]*coEfficients[3];
            color[1]+=reflectedColor[1]*coEfficients[3];
            color[2]+=reflectedColor[2]*coEfficients[3];
        }

        color[0]=min(1.0,max(0.0,color[0]));
        color[1]=min(1.0,max(0.0,color[1]));
        color[2]=min(1.0,max(0.0,color[2]));

        return t;
    }
};

class Triangle: public Object
{
public:
    Point vertex[3];

    Triangle(Point A,Point B,Point C)
    {
        vertex[0]=A;
        vertex[1]=B;
        vertex[2]=C;
    }
    
    void draw()
    {
        glPushMatrix();
        glColor3f(color[0],color[1],color[2]);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(vertex[0].x,vertex[0].y,vertex[0].z);
            glVertex3f(vertex[1].x,vertex[1].y,vertex[1].z);
            glVertex3f(vertex[2].x,vertex[2].y,vertex[2].z);
        }
        glEnd();
        glPopMatrix();
    }
    
    double intersect(Ray ray,double *color,int level)
    {
        double v1[9]=
        {
            -ray.direction.x,vertex[1].x-vertex[0].x,vertex[2].x-vertex[0].x,
            -ray.direction.y,vertex[1].y-vertex[0].y,vertex[2].y-vertex[0].y,
            -ray.direction.z,vertex[1].z-vertex[0].z,vertex[2].z-vertex[0].z
        };

        Matrix matrixA=Matrix(v1);

        double v2[9]=
        {
            -ray.direction.x,ray.start.x-vertex[0].x,vertex[2].x-vertex[0].x,
            -ray.direction.y,ray.start.y-vertex[0].y,vertex[2].y-vertex[0].y,
            -ray.direction.z,ray.start.z-vertex[0].z,vertex[2].z-vertex[0].z
        };

        Matrix betaMatrix=Matrix(v2);

        double v3[9]=
        {
            -ray.direction.x,vertex[1].x-vertex[0].x,ray.start.x-vertex[0].x,
            -ray.direction.y,vertex[1].y-vertex[0].y,ray.start.y-vertex[0].y,
            -ray.direction.z,vertex[1].z-vertex[0].z,ray.start.z-vertex[0].z
        };

        Matrix gammaMatrix=Matrix(v3);

        double v4[9]=
        {
            ray.start.x-vertex[0].x,vertex[1].x-vertex[0].x,vertex[2].x-vertex[0].x,
            ray.start.y-vertex[0].y,vertex[1].y-vertex[0].y,vertex[2].y-vertex[0].y,
            ray.start.z-vertex[0].z,vertex[1].z-vertex[0].z,vertex[2].z-vertex[0].z
        };

        Matrix tMatrix=Matrix(v4);

        double A=matrixA.determinant();

        if(isSame(A,0))
            return -1;

        double beta=betaMatrix.determinant()/A;
        double gamma=gammaMatrix.determinant()/A;
        double t=tMatrix.determinant()/A;
        // cout<<beta<<" "<<gamma<<" "<<t<<endl;

        if(!(beta+gamma<1&&beta>0&&gamma>0&&t>0))
            t=-1;
        if(level==0||t<0)
            return t;

        Point intersectionPoint=ray.start+ray.direction*t;

        double intersectionPointColor[3]= {0,0,0};

        intersectionPointColor[0]=this->color[0];
        intersectionPointColor[1]=this->color[1];
        intersectionPointColor[2]=this->color[2];

        color[0]=intersectionPointColor[0]*coEfficients[0];
        color[1]=intersectionPointColor[1]*coEfficients[0];
        color[2]=intersectionPointColor[2]*coEfficients[0];

        //Lighting calculations
        Point normal=crossProduct(vertex[1]-vertex[0],vertex[2]-vertex[0]);
        normal=normal/vecLength(normal);

        if(dotProduct(normal,ray.direction*-1)<0)
        {
            normal=normal*-1;
        }

        for(PointLight pl:pointLights)
        {
            Ray L(pl.lightPosition,intersectionPoint-pl.lightPosition);

            double tBlock=INF;
            for(auto obj:objects)
            {
                double tempColor[3]={0,0,0};
                double tempIntersect=obj->intersect(L,tempColor,0);
                if(tempIntersect>0&&tempIntersect<tBlock)
                {
                    tBlock=tempIntersect;
                }
            }

            if(tBlock>0)
            {
                Point blocked_at=L.start+L.direction*tBlock;
                if(!isSamePoint(blocked_at,intersectionPoint))
                {
                    continue;
                }
            }

            Point directionToLight=L.direction*-1;

            double lambert=dotProduct(directionToLight,normal);
            Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
            R=R/vecLength(R);

            double phong=dotProduct(ray.direction*-1,R);

            lambert=max(0.0,lambert);
            phong=max(0.0,phong);

            color[0]+=pl.color[0]*coEfficients[1]*lambert*intersectionPointColor[0];
            color[1]+=pl.color[1]*coEfficients[1]*lambert*intersectionPointColor[1];
            color[2]+=pl.color[2]*coEfficients[1]*lambert*intersectionPointColor[2];

            color[0]+=pl.color[0]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[0];
            color[1]+=pl.color[1]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[1];
            color[2]+=pl.color[2]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[2];
        }

        for(SpotLight sl:spotLights)
        {
            Ray L(sl.pointLight.lightPosition,intersectionPoint-sl.pointLight.lightPosition);
            
            //Check cutoff angle
            double angle=angleInDegree(sl.lightDirection,L.direction);
            if(angle>sl.cutoffAngle)
            {
                continue;
            }

            double tBlock=INF;
            for(auto obj:objects)
            {
                double tempColor[3]={0,0,0};
                double tempIntersect=obj->intersect(L,tempColor,0);
                if(tempIntersect>0&&tempIntersect<tBlock)
                {
                    tBlock=tempIntersect;
                }
            }

            if(tBlock>0)
            {
                Point blocked_at=L.start+L.direction*tBlock;
                if(!isSamePoint(blocked_at,intersectionPoint))
                {
                    continue;
                }
            }

            Point directionToLight=L.direction*-1;
            
            double lambert=dotProduct(directionToLight,normal);
            Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
            R=R/vecLength(R);
            lambert=max(0.0,lambert);

            double phong=dotProduct(ray.direction*-1,R);
            phong=max(0.0,phong);

            //Diffuse
            color[0]+=sl.pointLight.color[0]*coEfficients[1]*lambert*intersectionPointColor[0];
            color[1]+=sl.pointLight.color[1]*coEfficients[1]*lambert*intersectionPointColor[1];
            color[2]+=sl.pointLight.color[2]*coEfficients[1]*lambert*intersectionPointColor[2];

            //Specular
            color[0]+=sl.pointLight.color[0]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[0];
            color[1]+=sl.pointLight.color[1]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[1];
            color[2]+=sl.pointLight.color[2]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[2];
        }

        if(level==recursionLevel)
            return t;

        Point directionToLight=ray.direction*-1;
        Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
        Ray reflectRay(intersectionPoint+R*EPSILON*3,R);

        Object *closestObj=NULL;
        double tMinimumDistance=INF;

        for(Object *obj:objects)
        {
            double tempColor[3]={0,0,0};
            double tempIntersect=obj->intersect(reflectRay,tempColor,0);

            if(tempIntersect>0&&tempIntersect<tMinimumDistance)
            {
                tMinimumDistance=tempIntersect;
                closestObj=obj;
            }
        }

        if(closestObj!=NULL)
        {
            double reflectedColor[3]={0,0,0};
            closestObj->intersect(reflectRay,reflectedColor,level+1);

            color[0]+=reflectedColor[0]*coEfficients[3];
            color[1]+=reflectedColor[1]*coEfficients[3];
            color[2]+=reflectedColor[2]*coEfficients[3];
        }

        color[0]=min(1.0,max(0.0,color[0]));
        color[1]=min(1.0,max(0.0,color[1]));
        color[2]=min(1.0,max(0.0,color[2]));

        return t;
    }
};

class General: public Object
{
public:
    double A,B,C,D,E,F,G,H,I,J;
    Point referencePoint;
    double height,width,length;

    General(double A,double B,double C,double D,double E,double F,double G,double H,double I,double J,Point referencePoint,double height,double width,double length)
    {
        this->A=A;
        this->B=B;
        this->C=C;
        this->D=D;
        this->E=E;
        this->F=F;
        this->G=G;
        this->H=H;
        this->I=I;
        this->J=J;

        this->referencePoint=referencePoint;
        this->height=height;
        this->width=width;
        this->length=length;
    }
    
    void draw(){}
    
    double intersect(Ray ray,double *color,int level)
    {
        Point direction=ray.direction;
        Point start=ray.start;

        double a=A*direction.x*direction.x+B*direction.y*direction.y+C*direction.z*direction.z+D*direction.x*direction.y+E*direction.x*direction.z+F*direction.y*direction.z;
        double b=2.0*A*start.x*direction.x+2*B*start.y*direction.y+2*C*start.z*direction.z+D*(start.x*direction.y+start.y*direction.x)+E*(start.x*direction.z+start.z*direction.x)+F*(start.y*direction.z+direction.y*start.z)+G*direction.x+H*direction.y+I*direction.z;
        double c=A*start.x*start.x+B*start.y*start.y+C*start.z*start.z+D*start.x*start.y+E*start.x*start.z+F*start.y*start.z+G*start.x+H*start.y+I*start.z+J;

        //Initializing t0 and t1 to -1,closestObj indicates no intersection
        double t0=-1,t1=-1;

        if(isSame(a,0))
        {
            if(b==0)
                t0=t1=-1;
            else
                t0=t1=-c/b;
        }
        else
        {
            double d=b*b-4*a*c;
            if(d<0)
            {
                t0=t1=-1;
            }
            else
            {
                t0=(-b-sqrt(d))/(a*2);
                t1=(-b+sqrt(d))/(a*2);
            }
        }

        //Calculate the intersection points
        Point p0=start+direction*t0;
        Point p1=start+direction*t1;

        //Check if the intersection points are within the bounds of the object
        if(!isSame(length,0)&&(p0.x<referencePoint.x||p0.x>referencePoint.x+length))
            t0=-1;
        if(!isSame(width,0)&&(p0.y<referencePoint.y||p0.y>referencePoint.y+width))
            t0=-1;
        if(!isSame(height,0)&&(p0.z<referencePoint.z||p0.z>referencePoint.z+height))
            t0=-1;

        if(!isSame(length,0)&&(p1.x<referencePoint.x||p1.x>referencePoint.x+length))
            t1=-1;
        if(!isSame(width,0)&&(p1.y<referencePoint.y||p1.y>referencePoint.y+width))
            t1=-1;
        if(!isSame(height,0)&&(p1.z<referencePoint.z||p1.z>referencePoint.z+height))
            t1=-1;

        color[0]=this->color[0];
        color[1]=this->color[1];
        color[2]=this->color[2];

        double t=-1;

        if(t0>0)
        {
            t=t0;
        }
        else if(t1>0)
        {
            t=t1;
        }
        if(level==0||t<0)
        {
            return t;
        }

        Point intersectionPoint=ray.start+ray.direction*t;

        double intersectionPointColor[3]= {0,0,0};

        intersectionPointColor[0]=this->color[0];
        intersectionPointColor[1]=this->color[1];
        intersectionPointColor[2]=this->color[2];

        color[0]=intersectionPointColor[0]*coEfficients[0];
        color[1]=intersectionPointColor[1]*coEfficients[0];
        color[2]=intersectionPointColor[2]*coEfficients[0];

        //Lighting calculations
        Point normal;

        normal.x=2*A*intersectionPoint.x+D*intersectionPoint.y+E*intersectionPoint.z+G;
        normal.y=2*B*intersectionPoint.y+D*intersectionPoint.x+F*intersectionPoint.z+H;
        normal.z=2*C*intersectionPoint.z+E*intersectionPoint.x+F*intersectionPoint.y+I;

        normal=normal/vecLength(normal);

        if(dotProduct(normal,ray.direction*-1)<0)
        {
            normal=normal*-1;
        }

        for(PointLight pl:pointLights)
        {
            Ray L(pl.lightPosition,intersectionPoint-pl.lightPosition);

            double tBlock=INF;
            for(auto obj:objects)
            {
                double tempColor[3]={0,0,0};
                double tempIntersect=obj->intersect(L,tempColor,0);
                if(tempIntersect>0&&tempIntersect<tBlock)
                {
                    tBlock=tempIntersect;
                }
            }

            if(tBlock>0)
            {
                Point blocked_at=L.start+L.direction*tBlock;
                if(!isSamePoint(blocked_at,intersectionPoint))
                {
                    continue;
                }
            }

            Point directionToLight=L.direction*-1;

            double lambert=dotProduct(directionToLight,normal);
            Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
            R=R/vecLength(R);

            double phong=dotProduct(ray.direction*-1,R);

            lambert=max(0.0,lambert);
            phong=max(0.0,phong);

            color[0]+=pl.color[0]*coEfficients[1]*lambert*intersectionPointColor[0];
            color[1]+=pl.color[1]*coEfficients[1]*lambert*intersectionPointColor[1];
            color[2]+=pl.color[2]*coEfficients[1]*lambert*intersectionPointColor[2];

            color[0]+=pl.color[0]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[0];
            color[1]+=pl.color[1]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[1];
            color[2]+=pl.color[2]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[2];
        }

        for(SpotLight sl:spotLights)
        {
            Ray L(sl.pointLight.lightPosition,intersectionPoint-sl.pointLight.lightPosition);
            
            //Check cutoff angle
            double angle=angleInDegree(sl.lightDirection,L.direction);
            if(angle>sl.cutoffAngle)
            {
                continue;
            }

            double tBlock=INF;
            for(auto obj:objects)
            {
                double tempColor[3]={0,0,0};
                double tempIntersect=obj->intersect(L,tempColor,0);
                if(tempIntersect>0&&tempIntersect<tBlock)
                {
                    tBlock=tempIntersect;
                }
            }

            if(tBlock>0)
            {
                Point blocked_at=L.start+L.direction*tBlock;
                if(!isSamePoint(blocked_at,intersectionPoint))
                {
                    continue;
                }
            }

            Point directionToLight=L.direction*-1;
            
            double lambert=dotProduct(directionToLight,normal);
            Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
            R=R/vecLength(R);
            lambert=max(0.0,lambert);

            double phong=dotProduct(ray.direction*-1,R);
            phong=max(0.0,phong);

            //Diffuse
            color[0]+=sl.pointLight.color[0]*coEfficients[1]*lambert*intersectionPointColor[0];
            color[1]+=sl.pointLight.color[1]*coEfficients[1]*lambert*intersectionPointColor[1];
            color[2]+=sl.pointLight.color[2]*coEfficients[1]*lambert*intersectionPointColor[2];

            //Specular
            color[0]+=sl.pointLight.color[0]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[0];
            color[1]+=sl.pointLight.color[1]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[1];
            color[2]+=sl.pointLight.color[2]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[2];
        }

        if(level==recursionLevel)
            return t;

        Point directionToLight=ray.direction*-1;
        Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
        Ray reflectRay(intersectionPoint+R*EPSILON*3,R);

        Object *closestObj=NULL;
        double tMinimumDistance=INF;

        for(Object *obj:objects)
        {
            double tempColor[3]={0,0,0};
            double tempIntersect=obj->intersect(reflectRay,tempColor,0);

            if(tempIntersect>0&&tempIntersect<tMinimumDistance)
            {
                tMinimumDistance=tempIntersect;
                closestObj=obj;
            }
        }

        if(closestObj!=NULL)
        {
            double reflectedColor[3]={0,0,0};
            closestObj->intersect(reflectRay,reflectedColor,level+1);

            color[0]+=reflectedColor[0]*coEfficients[3];
            color[1]+=reflectedColor[1]*coEfficients[3];
            color[2]+=reflectedColor[2]*coEfficients[3];
        }

        color[0]=min(1.0,max(0.0,color[0]));
        color[1]=min(1.0,max(0.0,color[1]));
        color[2]=min(1.0,max(0.0,color[2]));

        return t;
    }
};

extern Point lightPosition;
extern Point u;
extern Point r;
extern Point l;

class Floor: public Object
{
public:
    double floorWidth;
    double tileWidth;

    Floor(double floorWidth,double tileWidth)
    {
        this->floorWidth=floorWidth;
        this->tileWidth=tileWidth;
        setCoEfficients(0.4,0.2,0.2,0.2);
    }
    
    void draw()
    {
        glPushMatrix();
        Point bottomLeft(-floorWidth/2,-floorWidth/2,0);

        int countCell=ceil(floorWidth/tileWidth);

        for(int i=0;i<countCell;i++)
        {
            for(int j=0;j<countCell;j++)
            {
                Point A=bottomLeft+Point(1,0,0)*i*tileWidth+Point(0,1,0)*j*tileWidth;
                Point B=bottomLeft+Point(1,0,0)*(i+1)*tileWidth+Point(0,1,0)*j*tileWidth;
                Point C=bottomLeft+Point(1,0,0)*(i+1)*tileWidth+Point(0,1,0)*(j+1)*tileWidth;
                Point D=bottomLeft+Point(1,0,0)*i*tileWidth+Point(0,1,0)*(j+1)*tileWidth;

                if((i+j)%2==0)
                {
                    glColor3f(1,1,1);
                }
                else
                {
                    glColor3f(0,0,0);
                }

                glBegin(GL_QUADS);
                {
                    glVertex3f(A.x,A.y,A.z);
                    glVertex3f(B.x,B.y,B.z);
                    glVertex3f(C.x,C.y,C.z);
                    glVertex3f(D.x,D.y,D.z);
                }
                glEnd();
            }
        }

        glPopMatrix();
    }
    
    double intersect(Ray ray,double *color,int level)
    {
        Point normal=Point(0,0,1);
        if(dotProduct(normal,ray.direction*-1)<0)
        {
            normal=normal*-1;
        }

        double down=dotProduct(normal,ray.direction);
        if(isSame(down,0))
            return -1;
        double up=-dotProduct(normal,ray.start);

        double t=up/down;

        Point intersectionPoint=ray.start+ray.direction*t;

        if(!(intersectionPoint.x>=-floorWidth/2&&intersectionPoint.x<=floorWidth/2&&intersectionPoint.y>=-floorWidth/2&&intersectionPoint.y<=floorWidth/2))
            return -1;

        if(t<0||level==0)
            return t;

        //Calculate the tile index
        int i=floor((intersectionPoint.x+floorWidth/2)/tileWidth);
        int j=floor((intersectionPoint.y+floorWidth/2)/tileWidth);

        double intersectionPointColor[3];

        if(useTexture) 
        {
            //Calculate local coordinates within the tile
            double localX=(intersectionPoint.x+floorWidth/2)-i*tileWidth;
            double localY=(intersectionPoint.y+floorWidth/2)-j*tileWidth;
            
            //Map local coordinates to [0,1] for the texture
            double u=localX/tileWidth;
            double v=localY/tileWidth;
            Color texColor=sampleTexture(u,v);
            intersectionPointColor[0]=texColor.r;
            intersectionPointColor[1]=texColor.g;
            intersectionPointColor[2]=texColor.b;
        }
        else 
        {
            //CheckerBoard Pattern
            if((i+j)%2==0)
            {
                intersectionPointColor[0]=1;
                intersectionPointColor[1]=1;
                intersectionPointColor[2]=1;
            }
            else
            {
                intersectionPointColor[0]=0;
                intersectionPointColor[1]=0;
                intersectionPointColor[2]=0;
            }
        }

        color[0]=intersectionPointColor[0]*coEfficients[0];
        color[1]=intersectionPointColor[1]*coEfficients[0];
        color[2]=intersectionPointColor[2]*coEfficients[0];

        for(PointLight pl:pointLights)
        {
            Ray L(pl.lightPosition,intersectionPoint-pl.lightPosition);

            double tBlock=INF;
            for(auto obj:objects)
            {
                double tempColor[3]={0,0,0};

                //Holds the nearest intersection point with the object
                double tempIntersect=obj->intersect(L,tempColor,0);
                if(tempIntersect>0&&tempIntersect<tBlock)
                {
                    //Store the nearest intersection point
                    tBlock=tempIntersect;
                }
            }

            if(tBlock>0)
            {
                Point blocked_at=L.start+L.direction*tBlock;
                if(!isSamePoint(blocked_at,intersectionPoint))
                {
                    continue;
                }
            }

            Point directionToLight=L.direction*-1;

            double lambert=dotProduct(directionToLight,normal);
            Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
            R=R/vecLength(R);
            lambert=max(0.0,lambert);

            double phong=dotProduct(ray.direction*-1,R);
            phong=max(0.0,phong);

            //Diffuse
            color[0]+=pl.color[0]*coEfficients[1]*lambert*intersectionPointColor[0];
            color[1]+=pl.color[1]*coEfficients[1]*lambert*intersectionPointColor[1];
            color[2]+=pl.color[2]*coEfficients[1]*lambert*intersectionPointColor[2];


            //Specular
            color[0]+=pl.color[0]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[0];
            color[1]+=pl.color[1]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[1];
            color[2]+=pl.color[2]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[2];
        }

        for(SpotLight sl:spotLights)
        {
            Ray L(sl.pointLight.lightPosition,intersectionPoint-sl.pointLight.lightPosition);
            
            //Check cutoff angle
            double angle=angleInDegree(sl.lightDirection,L.direction);
            if(angle>sl.cutoffAngle)
            {
                continue;
            }
            double tBlock=INF;
            for(auto obj:objects)
            {
                double tempColor[3]={0,0,0};
                double tempIntersect=obj->intersect(L,tempColor,0);
                if(tempIntersect>0&&tempIntersect<tBlock)
                {
                    tBlock=tempIntersect;
                }
            }

            if(tBlock>0)
            {
                Point blocked_at=L.start+L.direction*tBlock;
                if(!isSamePoint(blocked_at,intersectionPoint))
                {
                    continue;
                }
            }

            Point directionToLight=L.direction*-1;
            
            double lambert=dotProduct(directionToLight,normal);
            Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
            R=R/vecLength(R);
            lambert=max(0.0,lambert);

            double phong=dotProduct(ray.direction*-1,R);
            phong=max(0.0,phong);

            //diffuse
            color[0]+=sl.pointLight.color[0]*coEfficients[1]*lambert*intersectionPointColor[0];
            color[1]+=sl.pointLight.color[1]*coEfficients[1]*lambert*intersectionPointColor[1];
            color[2]+=sl.pointLight.color[2]*coEfficients[1]*lambert*intersectionPointColor[2];

            //Specular
            color[0]+=sl.pointLight.color[0]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[0];
            color[1]+=sl.pointLight.color[1]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[1];
            color[2]+=sl.pointLight.color[2]*coEfficients[2]*pow(phong,shine)*intersectionPointColor[2];
        }

        if(level==recursionLevel)
            return t;

        Point directionToLight=ray.direction*-1;
        Point R=normal*dotProduct(directionToLight,normal)*2-directionToLight;
        Ray reflectRay(intersectionPoint+R*EPSILON*3,R);

        Object *closestObj=NULL;
        double tMinimumDistance=INF;

        for(Object *obj:objects)
        {
            double tempColor[3]={0,0,0};
            double tempIntersect=obj->intersect(reflectRay,tempColor,0);

            if(tempIntersect>0&&tempIntersect<tMinimumDistance)
            {
                tMinimumDistance=tempIntersect;
                closestObj=obj;
            }
        }

        if(closestObj!=NULL)
        {
            double reflectedColor[3]={0,0,0};
            closestObj->intersect(reflectRay,reflectedColor,level+1);

            color[0]+=reflectedColor[0]*coEfficients[3];
            color[1]+=reflectedColor[1]*coEfficients[3];
            color[2]+=reflectedColor[2]*coEfficients[3];
        }

        color[0]=min(1.0,max(0.0,color[0]));
        color[1]=min(1.0,max(0.0,color[1]));
        color[2]=min(1.0,max(0.0,color[2]));

        return t;
    }
};

#endif