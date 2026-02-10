#include "vec3.h"
#include "namesp.h"
using namespace consts;
vec3::vec3(void)
{
    x = 0.0;
    y = 0.0;
    z = 0.0;
}
vec3::vec3(double v_x, double v_y, double v_z)
{
    x = v_x;
    y = v_y;
    z = v_z;
}
vec3::~vec3(){}

vec3 vec3::operator=(vec3 v)
{
    return vec3(x = v.x, y = v.y, z = v.z);
}
vec3 vec3::operator+(vec3 v)
{
    return vec3(x + v.x, y+ v.y,  z + v.z);
}
vec3 vec3::operator-(const vec3& v) const
{
    return vec3(x - v.x, y - v.y, z - v.z);
}

vec3 operator*(double n, vec3 v) //friend non member function
{
    //return vec3(v.x * n, v.y * n, v.z * n);
    return v * n;
}

vec3 vec3::operator*(const double n) const
{
    return vec3(x * n, y * n, z * n);
}
vec3 vec3::operator-() const //-v 
{
    return vec3(-x, -y, -z);
}
double vec3::operator*(vec3 v)
{
    return double(x * v.x+ y * v.y+ z * v.z);
}


vec3 vec3::operator/(double n)
{
    return vec3(x / n, y / n, z / n);
}

vec3 vec3::crossPdt(const vec3& v2) const
{    
    return vec3((y * v2.z) - (z * v2.y), (z * v2.x) - (x * v2.z), (x * v2.y) - (y * v2.x));
}

double vec3::getLength() const
{
    return  sqrt(x * x + y * y + z * z);
}
vec3 vec3::normalization()
{
    double length = getLength();
    x = x / length;
    y = y / length;
    z = z / length;
    return vec3(x, y, z);
}
void vec3::show() const
{
    cout << "x:" << x << "  y:" << y << "  z:" << z << endl;
}

double vec3::getDist(const vec3& v) const
{
    double tmp = SQR(x - v.x) + SQR(y - v.y) + SQR(z - v.z);
    return sqrt(tmp);
}

double getAngle(const vec3& v0, const vec3& v1) {
    double angle=0.0;
    double tmp;
    tmp = (v0 * v1) / (v0.getLength() * v1.getLength());
    angle = acos(tmp);
    if (isnan(angle)) {
        cout << "cos(sita)=" << tmp<< endl;
        throw std::runtime_error("Error! fabs(cos(sita)) should not greater then 1!");
    }
    else   angle = acos(tmp);
    return angle;
}

double TriplePdt(const vec3& v1, const vec3& v2, const vec3& v3)
{
    return  v1.crossPdt(v2) * v3;
}

//rotate about fixed axes.Rotation Order: first Z_alpha, then Y_beta, finally Z_gamma
//this matrix has been checked by myself, it is correct.
// this matrix has been tested and verified to be correct for ZYZ Euler angles.
//in our code, we use the ZYZ Euler angles to rotate the vector.
vec3 rotateVec3_ZYZ(const vec3& v, double alpha, double beta, double gamma)
{
    double ca = cos(alpha), sa = sin(alpha);
    double cb = cos(beta), sb = sin(beta);
    double cg = cos(gamma), sg = sin(gamma);

    double m00 = cg * cb * ca - sg * sa;
    double m01 = -cg * cb * sa - sg * ca;
    double m02 = cg * sb;

    double m10 = sg * cb * ca + cg * sa;
    double m11 = -sg * cb * sa + cg * ca;
    double m12 = sg * sb;

    double m20 = -sb * ca;
    double m21 = sb * sa;
    double m22 = cb;

    return vec3(
        v.x * m00 + v.y * m01 + v.z * m02,
        v.x * m10 + v.y * m11 + v.z * m12,
        v.x * m20 + v.y * m21 + v.z * m22
    );
}
vec3 rotateVec3_ZYZ02(const vec3& v, double a, double b, double g)//Macke's version
{
    double c1 = cos(a), s1 = sin(a);
    double c2 = cos(b), s2 = sin(b);
    double c3 = cos(g), s3 = sin(g);

    double m00 = -c2 * s1 * s3 + c1 * c3;
    double m01 = -c2 * s1 * c3 - c1 * s3;
    double m02 = s2 * s1;

    double m10 = c2 * c1 * s3 + s1 * c3;
    double m11 = c2 * c1 * c3 - s1 * s3;
    double m12 = -s2 * c1;

    double m20 = s2 * s3;
    double m21 = s2 * c3;
    double m22 = c2;

    return vec3(
        v.x * m00 + v.y * m01 + v.z * m02,
        v.x * m10 + v.y * m11 + v.z * m12,
        v.x * m20 + v.y * m21 + v.z * m22
    );
}

vec3 rotateVec3_YZ(const vec3& v, double alpha, double beta)
{
    // alpha: rotation around the Z axis
    // beta:  rotation around the Y axis

    double ca = cos(alpha), sa = sin(alpha);
    double cb = cos(beta), sb = sin(beta);

    // First rotate around the Z axis (alpha), then around the Y axis (beta)
    // Both are with respect to fixed (global) axes
    // The total rotation matrix is: Ry(beta) * Rz(alpha)

    double m00 = cb * ca;
    double m01 = -cb * sa;
    double m02 = sb;

    double m10 = sa;
    double m11 = ca;
    double m12 = 0.0;

    double m20 = -sb * ca;
    double m21 = sb * sa;
    double m22 = cb;

    return vec3(
        v.x * m00 + v.y * m01 + v.z * m02,
        v.x * m10 + v.y * m11 + v.z * m12,
        v.x * m20 + v.y * m21 + v.z * m22
    );
}