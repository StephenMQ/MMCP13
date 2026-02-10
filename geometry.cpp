#include "geometry.h"
#include "namesp.h"
#include "functionList.h"
using namespace consts;

line::line(void){}

line::line(vec3 point_P0, vec3 direction_D0)
{
    point_p = point_P0;
    if (fabs(direction_D0.getLength() - 1.0) > 1000.0*eps) {
        cout << "fabs(direction_D0.getLength() - 1.0)= "<< fabs(direction_D0.getLength() - 1.0) <<endl;
        throw std::runtime_error("direction_D0.getLength() should be normalizationd to 1");    
    }
    else{ direction_D0 = direction_D0.normalization(); }
    direction_d = direction_D0;
}
line line::update(const vec3& point_p1, const vec3& direction_d1)
{      
    return line(point_p1, direction_d1);
}
double DistPointLine(const vec3& p, const line& l)///this function returns a double value - distance between a point and a line
{
    vec3 tmp_v;
    tmp_v = p - l.point_p;    
    double angle = getAngle(tmp_v, l.direction_d);
    double dist = 0.0;
    if (fabs(angle) < eps)  throw std::runtime_error("Error, the point is on the line!");
    else
    {
        dist = tmp_v.getLength() * sin(angle);
    }
    if(dist<0.0) throw std::runtime_error("Error, check function**DistPointLine()**!");
    return dist;
}

//==========================================================================


vec3 IntersecPointPlaneLine(const Plane3D& P, const line& l)
{    
    vec3 IntersecPoint_p;
    if (fabs(P.getNormal() * l.direction_d) < 1e-6)
        throw std::runtime_error("The line is parallel to the plane or the line is in the plane");
    else {
        double t = (P.getPoint() * P.getNormal() - l.point_p * P.getNormal()) / (l.direction_d * P.getNormal());
        IntersecPoint_p = l.point_p + l.direction_d * t;
    }
    return IntersecPoint_p;
}

sphere::sphere(void)
{
    vec3 v0(0.0, 0.0, 0.0);
    center = v0;
    radius = 0.0;
}
sphere::sphere(vec3 ct, double rd)
{
    center = ct;
    radius = rd;
}

bool CheckLineIntersectSphere(const sphere& sph, const line& l)
{
    vec3 v_tmp = sph.center - l.point_p;    
    double distanceOfCenterAndLine = DistPointLine(sph.center, l);
    //cout << "distanceOfCenterAndLine=" << distanceOfCenterAndLine<< "    sph.radius="<< sph.radius<<endl;
    if (distanceOfCenterAndLine > sph.radius) return false;//the line does not intersect the sphere at all
    else return true; //intersect(1 or 2 intersection points)
}
bool CheckPointOnSphere(const sphere& sph, const vec3& p)
{
    double dist = sph.center.getDist(p);
   // cout << "fabs(dist - sph.radius)=" << fabs(dist - sph.radius) << endl;
    //cout << "-----------------------------" << endl;
    if (fabs(dist - sph.radius) > 1000.0*eps) return false;//the point is not on the sphere
    else return true; //the point is on the sphere
}

/// <summary>
/// this function use vectors to calaulate the geometric points and relations.
/// </summary>
/// <param name="sph">the sphere</param>
/// <param name="l">the line</param>
/// <returns>intersection point in format vec3</returns>

vec3 findIntersecPointOfLineAndSphere(const sphere& sph, const line& l) 
{  
    vec3 IntersecPoint_p;

    if (!CheckLineIntersectSphere(sph, l)) throw std::runtime_error("Error, the line does not intersect the sphere at all!");
    vec3 tmp_v = sph.center - l.point_p;
    double dist_t;
    double hypotenuseTimesCosAngle = sph.center.getDist(l.point_p) * cos(getAngle(tmp_v, l.direction_d)); //hypotenuse- longest side of a right triangle 
  //  cout << "ray's initial point is:"; l.point_p.show();
    if (CheckPointOnSphere(sph, l.point_p))   //ray's initial point is on sphere!
    {        
        dist_t = 2.0 * hypotenuseTimesCosAngle;
      // cout << "ray's initial point is on sphere!" << endl;
    }
    else
    {
       // cout << "ray's initial point is NOT on sphere!" << endl;
        dist_t = hypotenuseTimesCosAngle - sqrt(SQR(sph.radius) - SQR(DistPointLine(sph.center, l)));

    }
    IntersecPoint_p = l.point_p + dist_t * l.direction_d;
    //recheck
    if(fabs(IntersecPoint_p.getLength() - sph.radius)>1000.0*eps)
    {
        cout << "fabs(IntersecPoint_p.getLength() - sph.radius)"<< fabs(IntersecPoint_p.getLength() - sph.radius) <<endl;
        throw std::runtime_error("Error, the IntersecPoint_p is not properly calculated!");
    }
        
    return IntersecPoint_p;
}

bool CheckPointOnEllipsoid(const ellipsoid& ellips, const vec3& p) {
    // Get ellipsoid parameters
    double a = ellips.getA();
    double b = ellips.getB();
    double c = ellips.getC();

    // Calculate left side of ellipsoid equation
    double left_side = (p.x * p.x) / (a * a) +
        (p.y * p.y) / (b * b) +
        (p.z * p.z) / (c * c);

    // Due to floating-point precision, we need to check if the value is
    // approximately equal to 1.0 within some small epsilon
    const double epsilon = 1e-6; // Adjust this value based on your precision needs

    return fabs(left_side - 1.0) < epsilon;
}




bool CheckLineIntersectEllipsoid(const ellipsoid& ellips, const line& l) {
    // Get ellipsoid parameters
    double a = ellips.getA();
    double b = ellips.getB();
    double c = ellips.getC();

    // Get line parameters
    vec3 p = l.point_p;
    vec3 d = l.direction_d;

    // normalization the ellipsoid to a unit sphere by scaling coordinates
    // This transforms the problem into checking intersection of a transformed line with a unit sphere

    // Scale factors
    double a2 = a * a;
    double b2 = b * b;
    double c2 = c * c;

    // Compute coefficients of the quadratic equation
    // Equation: A*t? + B*t + C = 0
    double A = (d.x * d.x) / a2 + (d.y * d.y) / b2 + (d.z * d.z) / c2;
    double B = 2 * ((p.x * d.x) / a2 + (p.y * d.y) / b2 + (p.z * d.z) / c2);
    double C = (p.x * p.x) / a2 + (p.y * p.y) / b2 + (p.z * p.z) / c2 - 1;

    // Compute discriminant
    double discriminant = B * B - 4 * A * C;

    // If discriminant is negative or zero, no real roots -> no intersection
    if (discriminant <= 1e-6) {
        return false;
    }

    // If discriminant is zero or positive, there's at least one intersection point
    return true;
}


vec3 findIntersecPointOfLineAndEllipsoid(const ellipsoid& ellips, const line& l)
{
    // Get ellipsoid parameters
    double a = ellips.getA();
    double b = ellips.getB();
    double c = ellips.getC();

    // Get line parameters
    vec3 p = l.point_p;
    vec3 d = l.direction_d;

    // normalization the ellipsoid to a unit sphere by scaling coordinates
    // This transforms the problem into checking intersection of a transformed line with a unit sphere

    // Scale factors
    double a2 = a * a;
    double b2 = b * b;
    double c2 = c * c;

    // Compute coefficients of the quadratic equation
    // Equation: A*t? + B*t + C = 0
    double A = (d.x * d.x) / a2 + (d.y * d.y) / b2 + (d.z * d.z) / c2;
    double B = 2 * ((p.x * d.x) / a2 + (p.y * d.y) / b2 + (p.z * d.z) / c2);
    double C = (p.x * p.x) / a2 + (p.y * p.y) / b2 + (p.z * p.z) / c2 - 1;

    // Compute discriminant
    double discriminant = B * B - 4 * A * C;
    vec3 intersection(0.0,0.0,0.0);
    // If discriminant is negative or zero, no real roots -> no intersection
    if (discriminant <= eps) {
        cout << "discriminant=" << discriminant << endl;
        throw std::runtime_error("error, please recheck!");
    }
    else {       
        double sqrt_discriminant = sqrt(discriminant);
        double t1 = (-B + sqrt_discriminant) / (2 * A);
        double t2 = (-B - sqrt_discriminant) / (2 * A);
        NaNCheck(t1);
        NaNCheck(t2);
        double t;

        double t_larger, t_smaller;
        if (fabs(t1) > fabs(t2)) { t_larger = t1; t_smaller = t2; }
        else {
            t_larger = t2; t_smaller = t1;
        }      

        if (CheckPointOnEllipsoid(ellips, l.point_p))   //ray's initial point is on ellipsoid!
        {
            t = t_larger;
        }
        else {
            t = t_smaller;
        }
        vec3 intsP(
            p.x + t * d.x,
            p.y + t * d.y,
            p.z + t * d.z);

        intersection = intsP;

    }

    // If discriminant is zero or positive, there's at least one intersection point
    return intersection;
}
/// <summary>
/// note: the given normal points outward from the ellipsoid's surface.
/// </summary>
/// <param name="ellips"></param>
/// <param name="p"></param>
/// <returns></returns>
vec3 normalEllipsoid(const ellipsoid& ellips, const vec3& p)
{
    // Get ellipsoid parameters
    double a = ellips.getA();
    double b = ellips.getB();
    double c = ellips.getC();

    // Compute the normal vector using the gradient
    vec3 normal(
        2.0 * p.x / (a * a),
        2.0 * p.y / (b * b),
        2.0 * p.z / (c * c)
    );

    // normalization the vector to get a unit normal (optional but often useful)
    normal = normal.normalization();
    return normal;

}



vec3 rotateVec3(const vec3& vertex, double gamma, double beta, double alpha)
{
    // Rotate vertex using Euler angles
    double cosGamma = cos(gamma);
    double sinGamma = sin(gamma);
    double cosBeta = cos(beta);
    double sinBeta = sin(beta);
    double cosAlpha = cos(alpha);
    double sinAlpha = sin(alpha);
    double x_new = vertex.x * (cosGamma * cosAlpha - sinGamma * cosBeta * sinAlpha) -
        vertex.y * (sinGamma * cosAlpha + cosGamma * cosBeta * sinAlpha) +
        vertex.z * (sinBeta * sinAlpha);
    double y_new = vertex.x * (cosGamma * sinAlpha + sinGamma * cosBeta * cosAlpha) +
        vertex.y * (-sinGamma * sinAlpha + cosGamma * cosBeta * cosAlpha) -
        vertex.z * (sinBeta * cosAlpha);
    double z_new = vertex.x * (sinGamma * sinBeta) +
        vertex.y * (cosGamma * sinBeta) +
        vertex.z * (cosBeta);
    return vec3(x_new, y_new, z_new);
}

vec3 rotateVec3_02(const vec3& vertex, double gamma, double beta, double alpha)
{
    // Rotate vertex using Euler angles Z(gamma) Y(beta) Z(alpha)!!! Note:middle axis is Y, not X£¡£¡£¡
    matrix Z_alpha = Z_sita(alpha);
    matrix Y_beta = Y_sita(beta);
    matrix Z_gamma = Z_sita(gamma);


    vec3 vertex_new = Z_gamma * Y_beta * Z_alpha * vertex;

    return vertex_new;
}
