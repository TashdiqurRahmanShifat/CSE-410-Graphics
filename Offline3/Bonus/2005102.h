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
    
    //Additional operators needed for PBR
    Point operator+(double scalar)
    {
        return Point(x+scalar,y+scalar,z+scalar);
    }
    
    Point operator-(double scalar)
    {
        return Point(x-scalar,y-scalar,z-scalar);
    }
    
    //Component-wise multiplication operator
    Point operator*(const Point& other) const
    {
        return Point(x *other.x,y *other.y,z *other.z);
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

//PBR Helper Functions
double distributionGGX(Point N,Point H,double roughness)
{
    double a=roughness*roughness;
    double a2=a*a;
    double NdotH=max(dotProduct(N,H),0.0);
    double NdotH2=NdotH*NdotH;

    double num=a2;
    double denom=(NdotH2*(a2-1.0)+1.0);
    denom=pi*denom*denom;

    return num/denom;
}

double geometrySchlickGGX(double NdotV,double roughness)
{
    double r=(roughness+1.0);
    double k=(r*r)/8.0;

    double num=NdotV;
    double denom=NdotV*(1.0-k)+k;

    return num/denom;
}

double geometrySmith(Point N,Point V,Point L,double roughness)
{
    double NdotV=max(dotProduct(N,V),0.0);
    double NdotL=max(dotProduct(N,L),0.0);
    double ggx2=geometrySchlickGGX(NdotV,roughness);
    double ggx1=geometrySchlickGGX(NdotL,roughness);

    return ggx1*ggx2;
}

Point fresnelSchlick(double cosTheta,Point F0)
{
    double oneMinusCos=1.0-cosTheta;
    double pow5=oneMinusCos*oneMinusCos*oneMinusCos*oneMinusCos*oneMinusCos;
    return F0+(Point(1.0,1.0,1.0)-F0)*pow5;
}

//Enhanced Fresnel for better dielectric handling
Point fresnelSchlickRoughness(double cosTheta,Point F0,double roughness)
{
    double oneMinusCos=1.0 -cosTheta;
    double pow5=oneMinusCos*oneMinusCos*oneMinusCos*oneMinusCos*oneMinusCos;
    double smoothness=1.0-roughness;
    return F0+(Point(max(smoothness,F0.x),max(smoothness,F0.y),max(smoothness,F0.z))-F0)*pow5;
}

//Gamma correction functions
double toLinear(double value)
{
    return (value<=0.04045)?value/12.92:pow((value+0.055)/1.055,2.4);
}

double toGamma(double value)
{
    return (value<=0.0031308)?12.92*value:1.055*pow(value,1.0/2.4)-0.055;
}

Point toLinearSpace(Point color)
{
    return Point(toLinear(color.x),toLinear(color.y),toLinear(color.z));
}

Point toGammaSpace(Point color)
{
    return Point(toGamma(color.x),toGamma(color.y),toGamma(color.z));
}

// PBR BRDF calculation function - IMPROVED for better shininess
Point calculatePBR(Point normal, Point viewDir, Point lightDir, Point lightColor, Point albedo, double metallic, double roughness, double distance)
{
    Point H = lightDir + viewDir;
    H = H / vecLength(H);

    // Improved attenuation - less aggressive falloff for better highlights
    double attenuation = 1.0 / (1.0 + 0.045 * distance + 0.0075 * distance * distance);
    Point radiance = lightColor * attenuation;

    // Calculate F0 (base reflectance) - Enhanced for better reflectance
    Point F0 = Point(0.04, 0.04, 0.04);
    F0 = Point(F0.x * (1.0 - metallic), F0.y * (1.0 - metallic), F0.z * (1.0 - metallic)) + Point(albedo.x * metallic, albedo.y * metallic, albedo.z * metallic);

    // Cook-Torrance BRDF
    Point F = fresnelSchlick(max(dotProduct(H, viewDir), 0.0), F0);
    double NDF = distributionGGX(normal, H, roughness);
    double G = geometrySmith(normal, viewDir, lightDir, roughness);

    Point numerator = F * NDF * G;
    double denominator = 4.0 * max(dotProduct(normal, viewDir), 0.0) * max(dotProduct(normal, lightDir), 0.0) + 0.0001;
    Point specular = numerator / denominator;

    // Enhanced specular contribution for more prominent highlights
    specular = specular * (2.0 + metallic); // Boost specular, especially for metals

    // Energy conservation
    Point kS = F;
    Point kD = Point(1.0, 1.0, 1.0) - kS;
    kD = kD * (1.0 - metallic);

    double NdotL = max(dotProduct(normal, lightDir), 0.0);
    return (kD * albedo / pi + specular) * radiance * NdotL;
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
    
    // PBR Material Properties
    double albedo[3];      // Base color/albedo
    double metallic;       // Metallic factor (0.0 = dielectric, 1.0 = metallic)
    double roughness;      // Surface roughness (0.0 = mirror, 1.0 = completely rough)
    double ao;             // Ambient occlusion

    Object()
    {
        // Initialize PBR properties with default values
        metallic = 0.0;      // Default to dielectric
        roughness = 0.5;     // Medium roughness
        ao = 1.0;            // No ambient occlusion
        albedo[0] = albedo[1] = albedo[2] = 0.5;  // Medium gray
    }

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
        // Set albedo same as color for backward compatibility
        albedo[0]=r;
        albedo[1]=g;
        albedo[2]=b;
    }
    
    void setShine(int shine) 
    {
        this->shine=shine;
        // IMPROVED: More perceptually linear roughness conversion
        this->roughness = sqrt(1.0 - (shine / 255.0));
        this->roughness = max(0.01, min(1.0, this->roughness));
    }

    void setCoEfficients(double ambient,double diffuse,double specular,double recursiveReflection)
    {
        coEfficients[0]=ambient;
        coEfficients[1]=diffuse;
        coEfficients[2]=specular;
        coEfficients[3]=recursiveReflection;
    }
    
    // PBR specific setters
    void setPBRMaterial(double albedoR, double albedoG, double albedoB, double metallic, double roughness, double ao = 1.0)
    {
        this->albedo[0] = albedoR;
        this->albedo[1] = albedoG;
        this->albedo[2] = albedoB;
        this->metallic = max(0.0, min(1.0, metallic));
        this->roughness = max(0.01, min(1.0, roughness));
        this->ao = max(0.0, min(1.0, ao));
        
        // Set legacy color for compatibility
        this->color[0] = albedoR;
        this->color[1] = albedoG;
        this->color[2] = albedoB;
    }

    // Universal PBR lighting calculation (defined after PointLight and SpotLight classes)
    void calculatePBRLighting(Point intersectionPoint, Point normal, Point viewDir, double* color);
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

        //Ambient color (simplified ambient occlusion)
        color[0]=this->albedo[0]*coEfficients[0]*this->ao;
        color[1]=this->albedo[1]*coEfficients[0]*this->ao;
        color[2]=this->albedo[2]*coEfficients[0]*this->ao;

        //Lighting calculations
        Point normal=intersectionPoint-center;
        normal=normal/vecLength(normal);

        if(dotProduct(normal,ray.direction*-1)<0)
        {
            normal=normal*-1;
        }

        Point V = (ray.direction * -1);  // View direction
        V = V / vecLength(V);

        // Use unified PBR lighting calculation
        this->calculatePBRLighting(intersectionPoint, normal, V, color);

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

        color[0]=this->albedo[0]*coEfficients[0]*this->ao;
        color[1]=this->albedo[1]*coEfficients[0]*this->ao;
        color[2]=this->albedo[2]*coEfficients[0]*this->ao;

        //Lighting calculations
        Point normal=crossProduct(vertex[1]-vertex[0],vertex[2]-vertex[0]);
        normal=normal/vecLength(normal);

        if(dotProduct(normal,ray.direction*-1)<0)
        {
            normal=normal*-1;
        }

        Point V = (ray.direction * -1);  // View direction
        V = V / vecLength(V);

        // Use unified PBR lighting calculation
        this->calculatePBRLighting(intersectionPoint, normal, V, color);

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

        color[0]=this->albedo[0]*coEfficients[0]*this->ao;
        color[1]=this->albedo[1]*coEfficients[0]*this->ao;
        color[2]=this->albedo[2]*coEfficients[0]*this->ao;

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

        Point V=(ray.direction*-1);  //View direction
        V=V/vecLength(V);

        //Use unified PBR lighting calculation
        this->calculatePBRLighting(intersectionPoint,normal,V,color);

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

        //Set albedo based on texture or checkerboard pattern
        this->albedo[0]=intersectionPointColor[0];
        this->albedo[1]=intersectionPointColor[1];
        this->albedo[2]=intersectionPointColor[2];

        color[0]=this->albedo[0]*coEfficients[0]*this->ao;
        color[1]=this->albedo[1]*coEfficients[0]*this->ao;
        color[2]=this->albedo[2]*coEfficients[0]*this->ao;

        Point V=(ray.direction*-1);  //View direction
        V=V/vecLength(V);

        //Use unified PBR lighting calculation
        this->calculatePBRLighting(intersectionPoint,normal,V,color);

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

//HDR Tone Mapping (Simple Reinhard)
Point toneMap(Point color) 
{
    return Point(color.x/(color.x+1.0),color.y/(color.y+1.0),color.z/(color.z+1.0));
}

//Implementation of PBR lighting calculation
void Object::calculatePBRLighting(Point intersectionPoint,Point normal,Point viewDir,double* color)
{
    //Convert albedo to linear space for proper PBR calculations
    Point albedoLinear=toLinearSpace(Point(this->albedo[0],this->albedo[1],this->albedo[2]));
    
    //Point lights
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

        Point lightDir=L.direction*-1;
        lightDir=lightDir/vecLength(lightDir);
        
        double distance=vecLength(intersectionPoint-pl.lightPosition);
        Point lightColorLinear=toLinearSpace(Point(pl.color[0],pl.color[1],pl.color[2]));
        
        Point Lo=calculatePBR(normal,viewDir,lightDir,lightColorLinear,albedoLinear,this->metallic,this->roughness,distance);
        
        //Enhanced contribution for better visual impact
        double lightingMultiplier=1.0+this->metallic*0.5; //Boost metallic objects
        color[0]+=Lo.x*coEfficients[1]*lightingMultiplier;
        color[1]+=Lo.y*coEfficients[1]*lightingMultiplier;
        color[2]+=Lo.z*coEfficients[1]*lightingMultiplier;
    }

    //Spot lights
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

        Point lightDir=L.direction*-1;
        lightDir=lightDir/vecLength(lightDir);
        
        double distance=vecLength(intersectionPoint-sl.pointLight.lightPosition);
        Point lightColorLinear=toLinearSpace(Point(sl.pointLight.color[0],sl.pointLight.color[1],sl.pointLight.color[2]));
        
        Point Lo=calculatePBR(normal,viewDir,lightDir,lightColorLinear,albedoLinear,this->metallic,this->roughness,distance);
        
        //Enhanced contribution for better visual impact
        double lightingMultiplier=1.0+this->metallic*0.5; 
        color[0]+=Lo.x*coEfficients[1]*lightingMultiplier;
        color[1]+=Lo.y*coEfficients[1]*lightingMultiplier;
        color[2]+=Lo.z*coEfficients[1]*lightingMultiplier;
    }
}

#endif