/*
In this header file, all function declarations are listed except for class member functions.
The function definitions can be found in the following different source files:
otherFunctions.cpp
rayTracing.cpp

*/
#pragma once

#include"matrix.h"
#include <cmath>
#include"vec3.h"
#include"geometry.h"
#include"prism.h"
#include"complex.h"
#include"templateFunctions.h"
#include"namesp.h"
#include <vector>

double rn(void);
double rn1(void);
double RNG1(void);
double RNG2(void);

vec3 isotropicUnitDirVec();
vec3 isotropicUnitDirVec(vec3& v_hat, vec3& u_hat);
vec3 isotropicUnitDirVec02(vec3& v_hat, vec3& u_hat);//(x,y,z)-->(v,u,e) 
vec3 isotropicUnitDirVec03(int positionNumber, vec3& v_hat, vec3& u_hat);
vec3 initialPhotonPosition(const vec3& v_hat, const vec3& u_hat, const vec3& initialPhotonDirection, double R);
vec3 initialPhotonPosition_02(double gamma, double beta, double r_base, double R_max);
vec3 initialPhotonPosition_03(double beta, double r_base, double R_max);

vec3 isotropicUnitDirVec04(int positionNumber, double& sita, double& phi, vec3& v_hat, vec3& u_hat);
vec3 isotropicUnitDirVec05(vec3& v_hat, vec3& u_hat);
vec3 initialPhotonPosition01(double R, const double& sita, const double& phi, vec3& isoUniVec);



void showProgressBar(int current, int total, int barWidth, int step);
double maxOf3Num(const double a, const double b, const double c);
vec3 reflect(const char* caller, const vec3& normal, const vec3& incident);
bool TIR_check(const double m, const vec3& normal, const vec3& incident);
vec3 refract(const char* caller, const vec3& normal, const vec3& incident, const double n21);
matrix GammaMtx(const vec3& u0, const vec3& v0, const vec3& u1, const vec3& v1);
bool CheckGammaMtx(const matrix& G); //this rotation matrix should be orthogonal and det=1
matrix Z_sita(double sita);
matrix Y_sita(double sita);
matrix X_sita(double sita);
matrix R_abg(double gamma, double beta, double alpha);
matrix R_abg_T(double s3, double s2, double s1);  //s3-gamma; s2-beta; s1-alpha


matrixC reflectionMatrix(const complex m, const double sita_i, const double sita_t);//m-refractive index; sita_i-incident angle; sita_t-refraction angle
matrixC totalReflectionMatrix(complex m, double sita_i);
matrixC totalReflectionMatrix02(complex m, double sita_i);
matrixC refractionMatrix01(complex m, double sita_i, double sita_t);//m-refractive index; sita_i-incident angle; sita_t-refraction angle
matrixC refractionMatrix02(const complex m, const double sita_i, const double sita_t);

void normalUpdateForSphere(vec3& normal, line& ray, vec3& intersP, const vec3& e_r, const sphere& sphericalParticle);
void normalUpdateForEllipsoid(vec3& normal, line& ray, vec3& intersP, const vec3& e_r, const ellipsoid& ellips);
void normalUpdateForPrism(vec3& normal, line& ray, IntersectionInfo& hit, const vec3& e_r, const NPrism& prism);

void ReflThenRefrGoIn(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp, matrixC& Jones_refract_go_in, matrixC& Jones_ExtenalReflection_go_out, const double m);
void ReflInsideThenRefrGoOut(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp_1, matrixC& Jones_i, matrixC& Jones_refract_go_out, const double& m);

void ReflThenRefrGoIn02(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp, matrixC& Jones_refract_go_in, matrixC& Jones_ExtenalReflection_go_out, const double m);
void ReflInsideThenRefrGoOut02(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp_1, matrixC& Jones_i, matrixC& Jones_refract_go_out, const double& m);

void ReflThenRefrGoIn03(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp, matrixC& Jones_refract_go_in, matrixC& Jones_ExtenalReflection_go_out, const complex m);
void ReflInsideThenRefrGoOut03(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp_1, matrixC& Jones_i, matrixC& Jones_refract_go_out, const complex& m);


vec3 Event_Reflection_external(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_1, const vec3& normal, vec3& vp_1, matrixC& Jones_R, const double m); //Gamma_i included
vec3 Event_Refraction_go_in(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_1, const vec3& normal, vec3& vp_1, matrixC& Jones_T, const double& m);
vec3 Event_Reflection_inside(const vec3& ep_1, const vec3& normal, vec3& vp_1, matrixC& Jones_R, const double& m, bool& totalRefl);
vec3 Event_Refraction_go_out(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_i, const vec3& normal, vec3& vp_1, const matrixC& Jones_R, matrixC& Jones_T, const double m);
matrix Mueller(const matrixC& in);
#define vec3UnityCheck(v) vec3UnityCheckImpl(v, #v)
void vec3UnityCheckImpl(vec3& v, const char* varName);
#define NaNCheck(v) NaNCheckImpl(v, #v)
void NaNCheckImpl(double& num, const char* varName);

#define matrixJonesNonZeroCheck(M) matrixJonesNonZeroCheckImpl(M, #M)
bool matrixJonesNonZeroCheckImpl(matrixC& M, const char* varName);

//void vec3UnityCheck(vec3& v);
void LawConsEnergyVerify(const char* caller, double sita_i, double sita_t);
//matrixC realToComplexMatrix(matrix& M_real);
void bin_and_collect(vector<matrix>& MlrMtxVector, const vec3& e_out, const vec3& PhotonDirection, const matrix& Mueller_Mtx, const int binNumber);
void bin_and_collect02(vector<matrix>& MlrMtxVector, const vec3& e_out, const vec3& PhotonDirection, const matrix& Mueller_Mtx, const int binNumber);

// Clamp function for double values
inline double clamp(double value, double minVal, double maxVal)
{
    return std::max(minVal, std::min(value, maxVal));
}






