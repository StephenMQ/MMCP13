#pragma once

#include <iostream>
#include <cmath>
#include"templateFunctions.h"
#include "namesp.h"


using namespace std;
using namespace consts;

class vec3
{
public:
    double x;
    double y;
    double z;


    vec3();
    vec3(double v_x, double v_y, double v_z);
    ~vec3();
    
   
  // void set(double v_x, double v_y, double v_z);
   
    vec3 operator=(vec3 v);
    bool operator==(const vec3& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

    bool almostEqual(const vec3& other) const {
        return fabs(x - other.x) < eps && fabs(y - other.y) < eps && fabs(z - other.z) < eps;
	}

    vec3 operator+(vec3 v);
    vec3 operator+(const vec3& v) const { return vec3(x + v.x, y + v.y, z + v.z); }
    vec3 operator-(const vec3& v) const;
   
    vec3 operator*(const double n) const;//time operation from right
    vec3 operator-() const;

    friend vec3 operator*(double n, vec3 v);//time operation from left; Not member function
 
    vec3 operator/(double n);

    double operator*(vec3 v);//dot product       
   double operator*(const vec3& v) const {return double(x * v.x + y * v.y + z * v.z);
   }

    vec3 crossPdt(const vec3& v2) const;//cross product
   

    double getLength() const;
    
    vec3 normalization();
    vec3 normalization() const 
    {
        return vec3(x / getLength(), y / getLength(), z / getLength());
    }
    double getDist(const vec3& v) const;
    
   
    void show() const;
};
double getAngle(const vec3& v0,const vec3 &v1);
double TriplePdt(const vec3& v1, const vec3& v2, const vec3& v3);
struct triVec3 //(x,y,z)-->(v,u,e) 	
{
    vec3 e_x, e_y, e_z;
};


vec3 rotateVec3_ZYZ(const vec3& vertex, double gamma, double beta, double alpha);
vec3 rotateVec3_ZYZ02(const vec3& vertex, double gamma, double beta, double alpha);
vec3 rotateVec3_YZ(const vec3& v, double alpha, double beta);
//class point3:public vec3{};
//class direction:public vec3 {};

