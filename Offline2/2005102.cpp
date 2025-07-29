#include <bits/stdc++.h>
#include "bitmap_image.hpp"
using namespace std;


const double PI=acos(-1);
const double EPSILON=1e-10;

class Point 
{
public:
    double x,y,z,w;
    Point()
    {
        this->x=0;
        this->y=0;
        this->z=0;
        this->w=1;
    }

    Point(double x,double y,double z,double w=1.0)
    {
        this->x=x;
        this->y=y;
        this->z=z;
        this->w=w;
    }

    //Normalize vector
    void normalize() 
    {
        double length=sqrt(x*x+y*y+z*z);
        if(length>EPSILON) 
        {
            x/=length;
            y/=length;
            z/=length;
        }
    }

    //Homogenize(w=1)
    void homogenize()
    {
        if(fabs(w)>EPSILON) 
        {
            x/=w;
            y/=w;
            z/=w;
            w=1.0;
        }
    }


    //Subtract
    Point operator-(const Point& p) const
    {
        return Point(x-p.x,y-p.y,z-p.z,0.0);
    }
    

    friend ostream& operator<<(ostream& os,const Point& p) 
    {
        os<<"("<<p.x<<","<<p.y<<","<<p.z<<","<<p.w<<")";
        return os;
    }


    //Dot product
    double dot(const Point& p) const 
    {
        return x*p.x+y*p.y+z*p.z;
    }

    //Cross product
    Point cross(const Point& p) const 
    {
        return Point(
            y*p.z-z*p.y,
            z*p.x-x*p.z,
            x*p.y-y*p.x,
            0.0
        );
    }
};


class Matrix 
{
public:
    double matrix[4][4];

    //Default constructor initialized to identity matrix
    Matrix() 
    {
        for(int i=0;i<4;i++) 
        {
            for(int j=0;j<4;j++) 
            {
                matrix[i][j]=(i==j)?1.0:0.0;
            }
        }
    }

    //Matrix multiplication
    Matrix operator*(const Matrix& mat) const 
    {
        Matrix result;
        for(int i=0;i<4;i++) 
        {
            for(int j=0;j<4;j++) 
            {
                result.matrix[i][j]=0;
                for (int k=0;k<4;k++) 
                {
                    result.matrix[i][j]+=matrix[i][k]*mat.matrix[k][j];
                }
            }
        }
        return result;
    }

    //Matrix transformation to a point
    Point transform(const Point& p) const 
    {
        Point result;
        result.x=matrix[0][0]*p.x+matrix[0][1]*p.y+matrix[0][2]*p.z+matrix[0][3]*p.w;
        result.y=matrix[1][0]*p.x+matrix[1][1]*p.y+matrix[1][2]*p.z+matrix[1][3]*p.w;
        result.z=matrix[2][0]*p.x+matrix[2][1]*p.y+matrix[2][2]*p.z+matrix[2][3]*p.w;
        result.w=matrix[3][0]*p.x+matrix[3][1]*p.y+matrix[3][2]*p.z+matrix[3][3]*p.w;
        return result;
    }
};


//Translation Matrix
Matrix createTranslationMatrix(double tx,double ty,double tz) 
{
    Matrix matrix;
    matrix.matrix[0][3]=tx;
    matrix.matrix[1][3]=ty;
    matrix.matrix[2][3]=tz;
    return matrix;
}

//Matrix Scaling
Matrix createScalingMatrix(double sx,double sy,double sz) 
{
    Matrix matrix;
    matrix.matrix[0][0]=sx;
    matrix.matrix[1][1]=sy;
    matrix.matrix[2][2]=sz;
    return matrix;
}


//Rodrigues rotation formula generation
Point rodriguesRotate(const Point& v,const Point& axis,double angle)
{

    const Point& a=axis; 
    
    double angleRad=angle*PI/180.0;

    
    //1st term=cos(theta)*v
    double cosTheta=cos(angleRad);
    Point term1(v.x*cosTheta,v.y*cosTheta,v.z*cosTheta);

    //2nd term=(1-cos(theta))*(a·v)*a 
    double dotProduct=a.dot(v);
    double oneMinusCosTheta=1.0-cosTheta;
    Point term2(a.x*dotProduct*oneMinusCosTheta,a.y*dotProduct*oneMinusCosTheta,a.z*dotProduct*oneMinusCosTheta);

    //3rd term=sin(theta)*(axv)
    Point crossProduct=a.cross(v);
    double sinTheta=sin(angleRad);
    Point term3(crossProduct.x*sinTheta,crossProduct.y*sinTheta,crossProduct.z*sinTheta);

    //R
    return Point(
        term1.x+term2.x+term3.x,
        term1.y+term2.y+term3.y,
        term1.z+term2.z+term3.z
    );
}

//Rotation Matrix
Matrix createRotationMatrix(double angle,double ax,double ay,double az) 
{
    
    //i,j,k unit vectors creation
    Point i(1,0,0);
    Point j(0,1,0);
    Point k(0,0,1);
    
    Point axis(ax,ay,az);
    //Normalize axis
    axis.normalize();
    
    Point c1=rodriguesRotate(i,axis,angle);
    Point c2=rodriguesRotate(j,axis,angle);
    Point c3=rodriguesRotate(k,axis,angle);

    Matrix matrix;
    
    matrix.matrix[0][0]=c1.x; 
    matrix.matrix[0][1]=c2.x; 
    matrix.matrix[0][2]=c3.x;

    matrix.matrix[1][0]=c1.y; 
    matrix.matrix[1][1]=c2.y; 
    matrix.matrix[1][2]=c3.y;

    matrix.matrix[2][0]=c1.z;
    matrix.matrix[2][1]=c2.z;
    matrix.matrix[2][2]=c3.z;

    //Already 4*4 identity matrix,so no need to initialize the last row and column
    
    return matrix;
}


//Stage-1
void stageOne() 
{
    ifstream inputFile("scene.txt");
    ofstream outputFile("stage1.txt");
    
    if(!inputFile.is_open()) 
    {
        cout<<"Can't Open 'scene.txt!"<<endl;
        return;
    }
    
    if(!outputFile.is_open()) 
    {
        cout<<"Can't Open 'stage1.txt!"<<endl;
        inputFile.close();
        return;
    }
    
    //Set precision for output
    outputFile<<fixed<<setprecision(7);
    
    //Skip the first 4 lines
    string line;
    for(int i=0;i<4;i++) 
    {
        getline(inputFile,line);
    }
    
    //Initialize the transformation matrix stack with identity matrix
    //Handling opengl's glPushMatrix & glPopMatrix
    stack<Matrix>matrixStack;
    matrixStack.push(Matrix());
    
    while(getline(inputFile,line)) 
    {
        if(line=="triangle") 
        {
            //Read points of the triangle
            Point p[3];
            for(int i=0;i<3;i++) 
            {
                double x,y,z;
                inputFile>>x>>y>>z;
                p[i]=Point(x,y,z);
            }
            
            //Transform each point
            for(int i=0;i<3;i++) 
            {
                Point transformedPoint=matrixStack.top().transform(p[i]);
                transformedPoint.homogenize();//Make w=1
                outputFile<<transformedPoint.x<<" "<<transformedPoint.y<<" "<<transformedPoint.z<<endl;
            }
            outputFile << endl;
        } 

        else if(line=="translate") 
        {
            //Read parameters
            double tx,ty,tz;
            inputFile>>tx>>ty>>tz;
            
            Matrix translationMatrix=createTranslationMatrix(tx,ty,tz);
            
            //Update the top of the stack
            matrixStack.top()=matrixStack.top()*translationMatrix;
            //cout<<matrixStack.size()<<endl;
        } 

        else if(line=="scale") 
        {
            //Read parameters
            double sx,sy,sz;
            inputFile>>sx>>sy>>sz;
            
            Matrix scalingMatrix=createScalingMatrix(sx,sy,sz);
            
            //Update the top of the stack
            matrixStack.top()=matrixStack.top()*scalingMatrix;
        } 
        
        else if(line=="rotate") 
        {
            //Read parameters
            double angle,ax,ay,az;
            inputFile>>angle>>ax>>ay>>az;
            
            Matrix rotationMatrix=createRotationMatrix(angle,ax,ay,az);
            
            //Update the top of the stack
            matrixStack.top()=matrixStack.top()*rotationMatrix;
        } 
        
        else if(line=="push") 
        {
            //Duplicate the matrix of top of the stack and push it
            matrixStack.push(matrixStack.top());
            //cout<<matrixStack.size()<<endl;
        }
        
        else if(line=="pop") 
        {
            //Remove the top matrix
            if(matrixStack.size()>1) 
            {
                matrixStack.pop();
            } 
            else 
            {
                cout<<"There is only one matrix in the stack(the initial one)"<<endl;
            }

            //cout<<matrixStack.size()<<endl;
        } 
        
        else if(line=="end")
        {
            //cout<<matrixStack.size()<<endl;
            break;
        }
    }
    
    inputFile.close();
    outputFile.close();
    
    cout<<"Stage 1 done.Output written to stage1.txt"<<endl;
}

//Stage-2
void stageTwo() 
{    
    ifstream inputFile("stage1.txt");
    ofstream outputFile("stage2.txt");
    
    if(!inputFile.is_open()) 
    {
        cout<<"Can't Open 'stage1.txt!"<<endl;
        return;
    }
    
    if(!outputFile.is_open()) 
    {
        cout<<"Can't Open 'stage2.txt!"<<endl;
        inputFile.close();
        return;
    }
    
    //Set precision for output
    outputFile<<fixed<<setprecision(7);
    
    //Read camera parameters from scene.txt
    ifstream sceneFile("scene.txt");
    if(!sceneFile.is_open()) 
    {
        cout<<"Can't Open 'scene.txt!"<<endl;
        inputFile.close();
        outputFile.close();
        return;
    }
    
    //Read eye,look,up
    Point eye,look,up;
    sceneFile>>eye.x>>eye.y>>eye.z;
    sceneFile>>look.x>>look.y>>look.z;
    sceneFile>>up.x>>up.y>>up.z;
    
    sceneFile.close();
    
    //View Transformation
    //Calculate l,r,u
    Point l=look-eye;
    l.normalize();
    
    Point r=l.cross(up);
    r.normalize();
    
    Point u=r.cross(l);
    //u is already normalized since r and l are normalized
    
    //Create translation matrix T to move eye to origin
    Matrix T;
    T.matrix[0][3]=-eye.x;
    T.matrix[1][3]=-eye.y;
    T.matrix[2][3]=-eye.z;
    
    //Create rotation matrix R
    Matrix R;
    R.matrix[0][0]=r.x;
    R.matrix[0][1]=r.y; 
    R.matrix[0][2]=r.z;

    R.matrix[1][0]=u.x; 
    R.matrix[1][1]=u.y; 
    R.matrix[1][2]=u.z;

    R.matrix[2][0]=-l.x; 
    R.matrix[2][1]=-l.y; 
    R.matrix[2][2]=-l.z;
    
    //Calculate view transformation matrix V
    Matrix V=R*T;
    
    //Read each point from stage1.txt and transformed by V
    string line;
    while(getline(inputFile,line))
    {
        if(line.empty()) 
        {
            outputFile<<endl;//Empty lines between triangles
            continue;
        }
        
        
        istringstream iss(line);
        double x,y,z;
        iss>>x>>y>>z;
        Point p(x,y,z);
        
        //Apply view transformation
        Point transformedPoint=V.transform(p);
        transformedPoint.homogenize();
        
        //Write in the file
        outputFile<<transformedPoint.x<<" "<<transformedPoint.y<<" "<<transformedPoint.z <<endl;
    }
    
    inputFile.close();
    outputFile.close();
    
    cout<<"Stage 2 done.Output written to stage2.txt"<<endl;
}

//Stage-3
void stageThree() 
{
    ifstream inputFile("stage2.txt");
    ofstream outputFile("stage3.txt");
    
    if(!inputFile.is_open()) 
    {
        cout<<"Can't Open 'stage2.txt!"<<endl;
        return;
    }
    
    if(!outputFile.is_open()) 
    {
        cout<<"Can't Open 'stage3.txt!"<<endl;
        inputFile.close();
        return;
    }
    
    //Set precision for output
    outputFile<<fixed<<setprecision(7);
    
    //Read perspective parameters from scene.txt
    ifstream sceneFile("scene.txt");
    if(!sceneFile.is_open()) 
    {
        cout<<"Can't Open 'scene.txt!"<<endl;
        inputFile.close();
        outputFile.close();
        return;
    }
    
    //Skip first 3 lines
    string line;
    for(int i=0;i<3;i++) {
        getline(sceneFile, line);
    }
    
    //Read fovY,aspectRatio,near,far
    double fovY,aspectRatio,near,far;
    sceneFile>>fovY>>aspectRatio>>near>>far;
    
    sceneFile.close();
    
    //Projection Matrix calculation

    //Degree to radian
    double fovYRad=fovY*PI/180.0;
    double fovXRad=fovYRad*aspectRatio;
    
    //Calculate r and t
    double t=near*tan(fovYRad/2.0);
    double r=near*tan(fovXRad/2.0);
    
    //Create Projection Matrix P
    Matrix P;

    P.matrix[0][0]=near/r;
    P.matrix[1][1]=near/t;

    P.matrix[2][2]=-(far+near)/(far-near);
    P.matrix[2][3]=-(2.0*far*near)/(far-near);

    P.matrix[3][2]=-1.0;
    P.matrix[3][3]=0.0;
    
    //Read each point from stage2.txt and transformed by P
    while(getline(inputFile,line)) 
    {
        if(line.empty())
        {
            outputFile<<endl;//Empty lines between triangles
            continue;
        }
        
        istringstream iss(line);
        double x,y,z;
        iss>>x>>y>>z;
        Point p(x,y,z);
        
        //Apply projection transformation
        Point transformedPoint=P.transform(p);
        transformedPoint.homogenize();
        
        //Write in the file
        outputFile<<transformedPoint.x<<" "<<transformedPoint.y<<" "<<transformedPoint.z <<endl;
    }
    
    inputFile.close();
    outputFile.close();
    
    cout<<"Stage 3 done.Output written to stage3.txt"<<endl;
}

//Stage-4
struct Color 
{
    unsigned char r,g,b;
    Color() 
    {
        this->r=0;
        this->g=0;
        this->b=0;
    }

    Color(unsigned char r,unsigned char g,unsigned char b) 
    {
        this->r=r;
        this->g=g;
        this->b=b;
    }

    //Generate random color
    static Color random() 
    {
        return Color
        (
            rand()%256,
            rand()%256, 
            rand()%256 
        );
    }
};


struct Triangle
{
    Point vertices[3];
    Color color;
    
    Triangle(){}
    
    Triangle(const Point& v1,const Point& v2,const Point& v3)
    {
        vertices[0]=v1;
        vertices[1]=v2;
        vertices[2]=v3;
        color=Color::random();//random color assigned
    }
};


vector<Triangle>readTrianglesFromStage3()
{
    vector<Triangle>triangles;
    ifstream inputFile("stage3.txt");
    
    if(!inputFile.is_open())
    {
        cout<<"Can't open 'stage3.txt'"<<endl;
        return triangles;
    }
    
    string line;
    Point vertices[3];
    int vertexIndex=0;
    
    while(getline(inputFile,line))
    {
        if(line.empty())
        {
            //Ends of Triangle block,resets index 0
            if(vertexIndex==3)
            {
                triangles.push_back(Triangle(vertices[0],vertices[1],vertices[2]));
            }
            vertexIndex=0;
            continue;
        }
        
        istringstream iss(line);
        double x,y,z;
        iss>>x>>y>>z;
        
        if(vertexIndex<3)
        {
            vertices[vertexIndex++]=Point(x,y,z);
        }
    }
    
    //Handle the last triangle if file doesn't end with an empty line
    if(vertexIndex==3)
    {
        triangles.push_back(Triangle(vertices[0], vertices[1], vertices[2]));
    }
    
    inputFile.close();
    return triangles;
}


//Config file reading
void readConfig(int& screenWidth,int& screenHeight,double& leftX,double& rightX,
                double& bottomY,double& topY,double& frontZ,double& rearZ)
{
    ifstream configFile("config.txt");
    
    if(!configFile.is_open())
    {
        cout<<"Can't open 'config.txt'!"<<endl;
        return;
    }
    
    //Read Window Size
    configFile>>screenWidth>>screenHeight;
    
    //Read left and right limit of X
    configFile>>leftX;
    rightX=-leftX; 
    
    //Read bottom and top limit of Y
    configFile>>bottomY;
    topY=-bottomY;  
    
    //Read front and rear limit of Z
    configFile>>frontZ>>rearZ;
    configFile.close();
}


double interpolateX(double y,const Point& p1,const Point& p2)
{
    if(fabs(p1.y-p2.y)<EPSILON) return p1.x;//If value is too small,it may cause divide by 0 problme,so return only p1.x
    return p1.x+(y-p1.y)*(p2.x-p1.x)/(p2.y-p1.y);
}



double interpolateZ(double x,double y,const Triangle& triangle)
{
    const Point& A=triangle.vertices[0];
    const Point& B=triangle.vertices[1]; 
    const Point& C=triangle.vertices[2];
    
    //Calculate plane using cross product
    Point AB=B-A;
    Point AC=C-A;
    Point normal=AB.cross(AC);
    
    // Plane equation coefficients
    double a=normal.x;
    double b=normal.y;
    double c=normal.z;
    double d=-(a*A.x+b*A.y+c*A.z);
    
    //Calculate z
    if(fabs(c)>EPSILON)
    {
        return -(a*x+b*y+d)/c;
    }
    
    return (A.z+B.z+C.z)/3.0;
}

//Check the area of triangle.If it is too small,then skip
bool isValidTriangle(const Triangle& triangle)
{
    const Point& A=triangle.vertices[0];
    const Point& B=triangle.vertices[1];
    const Point& C=triangle.vertices[2];
    
    Point AB=B-A;
    Point AC=C-A;
    Point cross=AB.cross(AC);
    
    return cross.dot(cross)>EPSILON*EPSILON;
}

bool compareByYDescending(const Point& a,const Point& b)
{
    return a.y>b.y;
}

void stageFour()
{
    //Read parameters
    int screenWidth,screenHeight;
    double leftX,rightX,bottomY,topY,frontZ,rearZ;
    readConfig(screenWidth,screenHeight,leftX,rightX,bottomY,topY,frontZ,rearZ);

    vector<Triangle>triangles=readTrianglesFromStage3();

    //Initialize buffers
    vector<vector<double>>zBuffer(screenHeight,vector<double>(screenWidth,rearZ));
    vector<vector<Color>>frameBuffer(screenHeight,vector<Color>(screenWidth,Color(0,0,0)));//storing the color of each pixel

    //Pixel dimensions
    double dx=(rightX-leftX)/screenWidth;
    double dy=(topY-bottomY)/screenHeight;
    double topLeftY=topY-dy/2.0;
    double topLeftX=leftX+dx/2.0;

    //Apply procedure
    for(const Triangle& triangle:triangles)
    {
        if(!isValidTriangle(triangle))
            continue;


        
        const Point& A=triangle.vertices[0];
        const Point& B=triangle.vertices[1];
        const Point& C=triangle.vertices[2];

        Point vertices[3]={A,B,C};
        //Sort vertices by Y co-ordinate
        std::sort(vertices,vertices+3,compareByYDescending);

        Point& topVertex=vertices[0];
        Point& middleVertex=vertices[1];
        Point& bottomVertex=vertices[2];

        //Calculate scanline range
        //Clipping is done here
        int topScanline=max(0,(int)((topY-topVertex.y)/dy));
        int bottomScanline=min(screenHeight-1,(int)((topY-bottomVertex.y)/dy));


        for (int row=topScanline;row<=bottomScanline;row++)
        {
            double scanlineY=topLeftY-row*dy;
            vector<double>intersections;

            //Check top-middle
            if(topVertex.y!=middleVertex.y && scanlineY>=min(topVertex.y,middleVertex.y)&&scanlineY<=max(topVertex.y,middleVertex.y))
            {
                intersections.push_back(interpolateX(scanlineY,topVertex,middleVertex));
            }

            //Check middle-bottom
            if(middleVertex.y!=bottomVertex.y && scanlineY>=min(middleVertex.y,bottomVertex.y)&&scanlineY<=max(middleVertex.y,bottomVertex.y))
            {
                intersections.push_back(interpolateX(scanlineY,middleVertex,bottomVertex));
            }

            //Check top-bottom
            if (topVertex.y!=bottomVertex.y && scanlineY>=min(topVertex.y,bottomVertex.y)&&scanlineY<=max(topVertex.y,bottomVertex.y))
            {
                intersections.push_back(interpolateX(scanlineY,topVertex,bottomVertex));
            }

            //cout<<intersections.size()<<endl;
            if(intersections.size()>=2)
            {
                sort(intersections.begin(),intersections.end());
                double leftXPoint=intersections.front();
                double rightXPoint=intersections.back();

                int leftCol=max(0,(int)((leftXPoint-leftX)/dx));
                int rightCol=min(screenWidth-1,(int)((rightXPoint-leftX)/dx));

                //Color pixels
                for(int col=leftCol;col<=rightCol;col++)
                {
                    double x=topLeftX+col*dx;
                    double z=interpolateZ(x,scanlineY,triangle);
                    
                    //If the current value is closer than previous,update it
                    if(z>=frontZ&&z<zBuffer[row][col]-EPSILON)
                    {
                        zBuffer[row][col]=z;
                        frameBuffer[row][col]=triangle.color;
                    }
                }
            }
        }
    }

    //Write into z-buffer.txt
    ofstream zBufferFile("z-buffer.txt");
    if(zBufferFile.is_open())
    {
        zBufferFile<<fixed<<setprecision(6);
        for(int row=0;row<screenHeight;row++)
        {
            for(int col=0;col<screenWidth;col++)
            {
                if(zBuffer[row][col]<rearZ)
                {
                    zBufferFile<<zBuffer[row][col]<<"\t";
                }
            }
            zBufferFile<<endl;
        }
        zBufferFile.close();
    }

    //Save image
    bitmap_image image(screenWidth,screenHeight);
    for(int row=0;row<screenHeight;row++)
    {
        for(int col=0;col<screenWidth;col++)
        {
            const Color& color=frameBuffer[row][col];
            image.set_pixel(col,row,color.r,color.g,color.b);
        }
    }
    image.save_image("out.bmp");
    cout<<"Stage 4 done.Output written to z-buffer.txt"<<endl;
}


int main()
{
    srand(time(0));

    stageOne();
    stageTwo();
    stageThree();
    stageFour();
    cout<<"All stages done"<<endl;
    return 0;
}