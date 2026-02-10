#pragma once
#include"templateFunctions.h"
#include"vec3.h"

class line
{
public:
	vec3 point_p;
	vec3 direction_d;
	//double x, y, z;
 
    line();

    
	line(vec3 point_P0, vec3 direction_D0);

    
    line(const line& other) : point_p(other.point_p), direction_d(other.direction_d) {}

  
    line& operator=(const line& other) {
        if (this != &other) { 
            this->point_p = other.point_p;
            this->direction_d = other.direction_d;
        }
        return *this;
    }
	line update(const vec3& point_p1, const vec3& direction_d1);

};
double DistPointLine(const vec3& p, const line& l);///this function returns a double value - distance between a point and a line

/// 3D Plane Class, Plane3D(const vec3 & point, vec3 & normal)
class Plane3D {
private:
    vec3 point;   // A point on the plane
    vec3 normal;   // Unit normal vector of the plane

public:
    // Constructor 1: From a point and normal vector
    Plane3D(const vec3 & point, const vec3 & normal)
        : point(point), normal(normal.normalization()) {
    }

    // Constructor 2: From three points (computes normal)
    Plane3D(const vec3& p1, const vec3& p2, const vec3& p3) {
        vec3 v1 = p2 - p1;
        vec3 v2 = p3 - p1;
        vec3 n = v1.crossPdt(v2).normalization(); // Compute and normalization normal vector
        point = p1;
        normal = n;
    }

    // Calculate distance from point to plane
    double distanceToPoint(const vec3& p) const {
        vec3 diff = p - point;
        return fabs(diff*normal); // Absolute value gives distance
    }

    // Check if point lies on the plane (with tolerance)
    bool containsPoint(const vec3& p, double epsilon = 1e-6) const {
        return distanceToPoint(p) < epsilon;
    }

    // Get the plane's normal vector
    vec3 getNormal() const { return normal; }

    // Get a reference point on the plane
    vec3 getPoint() const { return point; }
};

vec3 IntersecPointPlaneLine(const Plane3D & P, const line& l);

class sphere //sphere class in 3D space
{
public:
	vec3 center;
	double radius;
	sphere();
	sphere(vec3 ct, double rd);
};
/// <summary>
/// class ellipsoid,  ellipsoid(double a = 1.0, double b = 1.0, double c = 1.0) 
/// </summary>
class ellipsoid
{
private:
    double a;  // semi-major axis (x-axis radius)
    double b;  // semi-middle axis (y-axis radius)
    double c;  // semi-minor axis (z-axis radius)

public:
    // Constructor
    ellipsoid(double a = 1.0, double b = 1.0, double c = 1.0) {
        if (a <= 0 || b <= 0 || c <= 0) {
            throw std::invalid_argument("All axes must be positive");
        }
        this->a = a;
        this->b = b;
        this->c = c;
    }

    // Getter methods
    double getA() const { return a; }
    double getB() const { return b; }
    double getC() const { return c; }

    // Setter methods with validation
    void setA(double newA) {
        if (newA <= 0) throw std::invalid_argument("Axis must be positive");
        a = newA;
    }

    void setB(double newB) {
        if (newB <= 0) throw std::invalid_argument("Axis must be positive");
        b = newB;
    }

    void setC(double newC) {
        if (newC <= 0) throw std::invalid_argument("Axis must be positive");
        c = newC;
    }

};
/* enum class relativePositionOfLineAndSphere
{
	notIntersect, //d>r
    Intersect
};
 */
bool  CheckLineIntersectSphere(const sphere& sph, const line& l);
bool CheckPointOnSphere(const sphere& sph, const vec3& p);
vec3 findIntersecPointOfLineAndSphere(const sphere& sph, const line& l);

bool CheckPointOnEllipsoid(const ellipsoid& ellips, const vec3& p);
bool  CheckLineIntersectEllipsoid(const ellipsoid& ellips, const line& l);
vec3 findIntersecPointOfLineAndEllipsoid(const ellipsoid& ellips, const line& l);
vec3 normalEllipsoid(const ellipsoid& ellips, const vec3& p);
vec3 rotateVec3(const vec3& vertex, double gamma, double beta, double alpha);
vec3 rotateVec3_02(const vec3& vertex, double gamma, double beta, double alpha);

