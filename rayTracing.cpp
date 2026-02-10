#include"functionList.h"


using namespace consts;
vec3 reflect(const char* caller, const vec3& normal, const vec3& incident)
{
	const double cosI = -normal * incident;
	if (cosI < 0.0)
	{
		normal.show(); cout << "|normal|=" << normal.getLength() << endl;
		incident.show(); cout << "|incident|=" << incident.getLength() << endl;
		cout << "CosSita_i=" << cosI << endl;
		cout << "Called by: " << caller << endl;
		throw std::runtime_error("Error found in function \"reflect\"! CosSita_i must be positive.");

	}
	return incident + 2.0 * cosI * normal;
}
bool TIR_check(const double m, const vec3& normal, const vec3& incident)
{
	const double sita_i = acos(-normal * incident);
	if (m * sin(sita_i) > 1.0) return true;//TIR case
	else return false;
}
vec3 refract(const char* caller, const vec3& normal, const vec3& incident, const double m)
{
	const double cosI = -normal * incident;
	
	if (cosI < 0.0)
	{
		cout << "normal: "; normal.show();
		cout << "|normal|=" << normal.getLength() << endl;
		cout << "incident: "; incident.show();
		cout << "|incident|=" << incident.getLength() << endl;
		
		cout << "CosSita_i=" << cosI << endl;
		
		cout << "Called by: " << caller << endl;
		throw std::runtime_error("Error in function \"refract()\"! CosSita_i must be positive.");
	}
	const double sinT2 = SQR(1.0/ m) * (1.0 - cosI * cosI);
	if (sinT2 > 1.0) throw std::runtime_error("Error in function \"refract()\": case:Total Internal Reflection."); //TIR
	const double cosT = sqrt(1.0 - sinT2); 
	return (1.0/ m) * incident + ((1.0 / m) * cosI - cosT)* normal;
}

matrixC reflectionMatrix(const complex m, const double sita_i, const double sita_t)//m-refractive index; sita_i-incident angle; sita_t-refraction angle
{
	matrixC R(2, 2);
	complex R_1u, R_1v;
	R_1u = (m * cos(sita_i) - cos(sita_t)) / (m * cos(sita_i) + cos(sita_t));
	R_1v = (cos(sita_i) - m * cos(sita_t)) / (cos(sita_i) + m * cos(sita_t));
	R[0][0] = R_1u;	R[0][1] = 0.0;
	R[1][0] = 0.0;	R[1][1] = R_1v;

	return R;
}
matrixC totalReflectionMatrix(complex m, double sita_i)//m-refractive index; sita_i-incident angle; sita_t-refraction angle
{
	matrixC R(2, 2);
	complex R_u, R_v,CosSita_t(0.0,0.0);
	double tmp = 1.0 - SQR((1.0/m.re) * sin(sita_i));
	if(tmp>0.0) throw std::runtime_error("Error in function \"totalReflectionMatrix()\": double check, \"tmp<0.0\" should be met.");
	CosSita_t.im = sqrt(-tmp);
	
	R_u = (m * cos(sita_i) - CosSita_t) / (m * cos(sita_i) + CosSita_t);
	R_v = (cos(sita_i) - m * CosSita_t) / (cos(sita_i) + m * CosSita_t);
	R[0][0] = R_u;	R[0][1] = zeroComplex();
	R[1][0] = zeroComplex();	R[1][1] = R_v;

	return R;
}

matrixC totalReflectionMatrix02(complex m, double theta_i)
{
	constexpr double eps = 1e-12;

	matrixC R(2, 2);

	double ci = cos(theta_i);
	ci = clamp(ci, -1.0, 1.0);

	double si2 = max(0.0, 1.0 - ci * ci);

	// 用复数形式算 sinθt
	complex sin_t = (sqrt(si2)) / m;

	complex cos_t2 = 1.0 - sin_t * sin_t;

	// 强制选取虚部正方向 (TIR)
	complex cos_t = complex(0.0, sqrt(max(0.0, -cos_t2.re)));

	complex R_u = (m * ci - cos_t) / (m * ci + cos_t);
	complex R_v = (ci - m * cos_t) / (ci + m * cos_t);

	R[0][0] = R_u;
	R[0][1] = zeroComplex();
	R[1][0] = zeroComplex();
	R[1][1] = R_v;

	return R;
}



//The following code is based on eqs(2.27) in chaper 2, Yang2006. Notation!: may be wrong!!!
matrixC refractionMatrix01(complex m, double sita_i, double sita_t)//m-refractive index; sita_i-incident angle; sita_t-refraction angle
{
	matrixC T(2, 2);
	complex R_1u, R_1v, T_1u, T_1v;
	R_1u = (m * cos(sita_i) - cos(sita_t)) / (m * cos(sita_i) + cos(sita_t));
	R_1v = (cos(sita_i) - m * cos(sita_t)) / (cos(sita_i) + m * cos(sita_t));
	
	T_1u = sqrt(1.0 - R_1u * R_1u);
	T_1v = sqrt(1.0 - R_1v * R_1v);
	if(isnan(T_1u.re)|| isnan(T_1v.re))
	{
		cout << "testing..."<< endl;
		cout << "R_1u=" << R_1u.re << endl;
		//cout << "R_1u*R_1u=" << R_1u * R_1u << endl;
		//throw std::runtime_error("Error in function \"refractionMatrix()\": check input data in function sqrt().");
	}	

	T[0][0] = T_1u;	T[0][1] = 0.0;
	T[1][0] = 0.0;	T[1][1] = T_1v;

	return T;
}


matrixC refractionMatrix02(const complex m, const double sita_i, const double sita_t)//m-refractive index; sita_i-incident angle; sita_t-refraction angle
{
	matrixC T(2, 2);
	complex  T_1u, T_1v;
	

	T_1u = (2.0 * cos(sita_i)) / (m * cos(sita_i) + cos(sita_t));
	T_1v = (2.0 * cos(sita_i)) / (cos(sita_i) + m*cos(sita_t));
	
	T[0][0] = T_1u;	T[0][1] = 0.0;
	T[1][0] = 0.0;	T[1][1] = T_1v;

	return T;
}
matrix GammaMtx(const vec3& u0, const vec3& v0, const vec3& u1, const vec3& v1)//following notation in [Yang Ping,2006]
{
	matrix Gamma(2, 2);
	Gamma[0][0] = u1 * u0;	Gamma[0][1] = u1 * v0;
	Gamma[1][0] = v1 * u0;	Gamma[1][1] = v1 * v0;

	return Gamma;
}
bool CheckGammaMtx(const matrix& G) //this rotation matrix should be orthogonal and det=1
{
	//cout << "testing..." << endl;
	double g00 = G[0][0], g01 = G[0][1];
	double g10 = G[1][0], g11 = G[1][1];

	// ---- 1. 行正交 ----
	double r0_norm = g00 * g00 + g01 * g01;
	double r1_norm = g10 * g10 + g11 * g11;
	double dot01 = g00 * g10 + g01 * g11;

	// ---- 2. determinant ----
	double det = g00 * g11 - g01 * g10;

	bool ok =
		std::abs(r0_norm - 1.0) < eps &&
		std::abs(r1_norm - 1.0) < eps &&
		std::abs(dot01) < eps &&
		std::abs(det - 1.0) < eps;

	if (!ok)
	{
		std::cout << "\n[Gamma matrix INVALID]\n";
		std::cout << g00 << "  " << g01 << "\n";
		std::cout << g10 << "  " << g11 << "\n";

		std::cout << "row0 norm = " << r0_norm << "\n";
		std::cout << "row1 norm = " << r1_norm << "\n";
		std::cout << "dot       = " << dot01 << "\n";
		std::cout << "det       = " << det << "\n\n";
	}

	return ok;
}

//#define vec3UnityCheck(v) vec3UnityCheckImpl(v, #v)
void vec3UnityCheckImpl(vec3& v, const char* varName) {

	if ((v.getLength() - 1.0) > 100*eps) {
		cout << scientific << "Error in vector                                                              : " << varName << ", |" << varName << "| = " << v.getLength() << endl;
		
		throw std::runtime_error("Error in function \"vec3UnityCheck()\": vector length should be 1.");
	}
}
//#define NaNCheck(v) NaNCheckImpl(v, #v)
void NaNCheckImpl(double& num, const char* varName) {
	if (isnan(num)) {
		cout << "Error in number: " << varName << ", |" << varName << "| = " << num << endl;
		throw std::runtime_error("Error found. Ckeck NaN");
	}
}
void LawConsEnergyVerify(const char* caller, double sita_i, double sita_t)
{
	if ((sita_i < 0.0) || (sita_i * 180.0 / pi > 90.0) || (sita_t < 0.0)||(sita_t * 180.0 / pi > 90.0))
	{   
		cout << "sita_i= " << sita_i * 180.0 / pi << "(in degree), sita_t=" << sita_t * 180.0 / pi << "(in degree)" << endl;
		cout << "Called by: " << caller << endl;
		throw std::runtime_error("Error found. sita_i and sita_t both should be positive and less than 90 degree (0<=sita<=90)! Please Check!");

	}
	double Rp, Rs, Tp, Ts ;//p-parallel, s- perpendicular 
	Rp = SQR(tan(sita_i - sita_t))/ SQR(tan(sita_i + sita_t));
	Rs = SQR(sin(sita_i - sita_t)) / SQR(sin(sita_i + sita_t));
	Tp = (sin(2.0 * sita_i) * sin(2.0 * sita_t)) / SQR(sin(sita_i + sita_t) * cos(sita_i - sita_t));
	Ts = (sin(2.0 * sita_i) * sin(2.0 * sita_t)) / SQR(sin(sita_i + sita_t));
	//cout << "testing..." << endl;
	if (abs(Rp + Tp - 1.0) > 1.0*eps|| abs(Rs + Ts - 1.0) > 1.0*eps) {
		
		cout << "Rp= " << Rp << ", Tp=" << Tp << endl;
		cout << "Rs= " << Rs << ", Ts=" << Ts << endl;
		cout << "abs(Rp + Tp - 1.0)= " << abs(Rp + Tp - 1.0) << ", abs(Rs + Ts - 1.0)=" << abs(Rs + Ts - 1.0) << endl;
		throw std::runtime_error("Error found in function LawConsEnergyVerify. Please Check!");
	}



}

//=================================================================================
/// <summary>
/// this function simulate the reflection event and update the photo state via &param, which includes (position, direction, Jones matrix)
/// </summary>
/// <param name="e_i">photon incident direction-unit vector </param>
/// <param name="e_r">reflected direction-unit vector </param>
/// <param name="e_t">refracted or transmitted direction-unit vector </param> 
/// <param name="n">crystal surface normal-unit vector</param>
/// <param name="Jones_M">complex matrix 2*2</param>
/// <param name="m">refractive index-complex number</param>
//(x,y,z)-->(v,u,e) 	
vec3 Event_Reflection_external(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_1, const vec3& normal, vec3& vp_1, matrixC& Jones_R, const double m) //Gamma_i included
{
	vec3 v1,u1;
	vec3 vr, ur,e_r;
	vec3 vs, us, e_s;
	double CosSita_i = -normal * e_1;
	
	double sita_i = acos(CosSita_i);
	if (isnan(sita_i)) throw std::runtime_error("Error #1! The domain of the \"acos()\" function is [-1,1].");	
	if (abs(sin(sita_i)) < eps) { v1=v_0; }//throw std::runtime_error("Error #1.2! Division by zero!");
	else v1 = e_1.crossPdt(normal) / sin(sita_i);

	vp_1 = v1;//tracing vector vp_1
	u1 = e_1.crossPdt(v1);
	matrix Gamma1 = GammaMtx(u_0, v_0, u1, v1);

	double tmp3 = sin(sita_i) / m;
	if (abs(tmp3) > 1.0 + eps) throw std::runtime_error("Error #1.1! The domain of the \"acos()\" function is [-1,1].");
	double sita_t = asin(tmp3);
	//e_r = e_1 + 2.0 * cos(sita_i) * normal;	
	e_r = reflect(__FUNCTION__, normal, e_0);
	vr = v1; 
	ur = e_r.crossPdt(vr);
	
	matrixC Refl1 = reflectionMatrix(m, sita_i, sita_t);

	double tmp1 = e_1 * e_r;
	if (abs(tmp1) > 1.0 + eps) throw std::runtime_error("Error #2! The domain of the \"acos()\" function is [-1,1].");
	double sita_s = acos(tmp1);

	double tmp2 = sin(sita_s);	
	if (abs(tmp2) < eps) vs = v_0;
	else vs = e_1.crossPdt(e_r) / tmp2;
	us = e_r.crossPdt(vs);
	matrix Gamma_s = GammaMtx(ur, vr, us, vs);

	vec3 u_tmp = e_0.crossPdt(vs);
	matrix Gamma_1i = GammaMtx(u_tmp, vs, u_0, v_0);
	Jones_R = Gamma_s * Refl1 * Gamma1 * Gamma_1i;
	//Jones_R = Gamma_s * Refl1 * Gamma1;//this is maybe wrong, just for test
	vec3UnityCheck(v1);
	vec3UnityCheck(u1);
	LawConsEnergyVerify(__FUNCTION__, sita_i, sita_t);
	return e_r;//reflected unit vector
}


vec3 Event_Refraction_go_in(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_1, const vec3& normal, vec3& vp_1, matrixC& Jones_T, const double& m)
{
	//m = 1.0 / m; this information has been considered in formulas;

	vec3 v1, u1;
	vec3 e_t;// vt, ut,
	double CosAngle = -normal * e_1;
	double sita_i = acos(CosAngle);
	if (isnan(sita_i)) throw std::runtime_error("Error #3! The domain of the \"acos()\" function is [-1,1].");
	

	if (abs(sin(sita_i)) < eps) { v1 = v_0; }//throw std::runtime_error("Error #1.2! Division by zero!");
	else v1 = e_1.crossPdt(normal) / sin(sita_i);
	vp_1 = v1;//tracing vector vp_1
	
	u1 = e_1.crossPdt(v1);

	double sita_t = asin(sin(sita_i) / m);
	if(isnan(sita_t)) throw std::runtime_error("Error #3.1! The domain of the \"asin()\" function is [-1,1].");
	//e_t = e_1 / m + (cos(sita_i) / m - cos(sita_t)) * normal;
	e_t= refract(__FUNCTION__, normal,e_0, m);
	//vt = v1;
	//ut = e_t.crossPdt(vt);
	matrix Gamma1 = GammaMtx(u_0, v_0, u1, v1);
	matrixC Trmt1 = refractionMatrix02(m, sita_i, sita_t);

	Jones_T = Trmt1 * Gamma1;
	vec3UnityCheck(v1);
	vec3UnityCheck(u1);
	LawConsEnergyVerify(__FUNCTION__,sita_i, sita_t);
	return e_t;//note! carefully check! here we use e_i to return the e_r, so before using this function, the e_i MUST be copied to another vec3;
}




vec3 Event_Reflection_inside(const vec3& ep_1, const vec3& normal,vec3& vp_1, matrixC& Jones_R, const double& m, bool& totalRefl)
{
	vec3 vp_1_i = vp_1;
	vec3 up_1 = ep_1.crossPdt(vp_1);	
	vec3 v_p, u_p, e_p; e_p = ep_1;
	vec3 vr, ur, e_r;
	vec3 minusNormal = -normal;//note: all other formulas assumed that all normal vector direction is facing the incoming rays!!!
	
	double tmp = -minusNormal * ep_1;
	double sita_i = acos(tmp);
	NaNCheck(sita_i);
		
	if (abs(sin(sita_i)) < eps)
	{
		v_p = vp_1;
		
	}
		
	else { v_p = ep_1.crossPdt(minusNormal) / sin(sita_i);
	//cout << "|v_p|=" << v_p.getLength()<<"  "; v_p.show();
	}
	vp_1 = v_p;
	u_p = e_p.crossPdt(v_p);

	//e_r = ep_1 + 2.0 * cos(sita_i) * minusNormal;
	e_r = reflect(__FUNCTION__, -normal, ep_1);
	vr = v_p;
	ur = e_r.crossPdt(vr);

	matrix Gamma = GammaMtx(up_1, vp_1_i, u_p, v_p);

	double tmp1 = sin(sita_i) * m;
	matrixC Refl_p(2, 2);

	if (fabs(tmp1) > 1.0) //
	{ totalRefl = true;  	
	Refl_p = totalReflectionMatrix(1.0 / m, sita_i);
	//cout << "TIR case!" << endl;
	//Refl_p.Identity(); simply wrong!!!
	}
	else{ 
		double sita_t = asin(tmp1); 
		NaNCheck(sita_t);
	    Refl_p = reflectionMatrix(1.0/m, sita_i, sita_t);
	}

	Jones_R = Refl_p * Gamma * Jones_R;
	vec3UnityCheck(v_p);
	vec3UnityCheck(u_p);
	return e_r;
}



vec3 Event_Refraction_go_out(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_i, const vec3& normal, vec3& vp_1, const matrixC& Jones_R, matrixC& Jones_T, const double m)
{
	vec3 vi,ui,e_t,vt,ut;
	vec3 ep_1 = e_i;
	vec3 up_1= ep_1.crossPdt(vp_1);
	vec3 vp_1_i = vp_1;

	vec3 minusNormal = -normal;//note: all other formulas assumed that all normal vector direction is facing the incoming rays!!!	
	double sita_i = acos(-minusNormal * e_i);
	NaNCheck(sita_i);
	if (abs(sin(sita_i)) < eps)  vi = vp_1;
	//{vec3 v0 = e_i.crossPdt(n_i) / sin(sita_i); v0.show(); cout << "sin(sita_i)=" << sin(sita_i) << ";    |v0|=" << v0.getLength() << endl; throw std::runtime_error("Error #5.1! Division by zero!");}
	else { vi = e_i.crossPdt(minusNormal) / sin(sita_i); }
	vp_1 = vi;//tracing vector vp_1
	ui = e_i.crossPdt(vi);

	double sita_t = asin(sin(sita_i) * m);
	if (isnan(sita_t)) {
		cout << "sita_i=" <<sita_i<<"; in degree:"<< sita_i * 180.0 / pi << endl;
		cout << "sin(sita_i) * m=" << sin(sita_i) * m << endl;
		throw std::runtime_error("Error #3.1! The domain of the \"asin()\" function is [-1,1]. case:Total Internal Reflection.");
	}
	
	//e_t = e_i * m + (cos(sita_i) * m - cos(sita_t)) * minusNormal;
	e_t = refract(__FUNCTION__, -normal, e_i, 1.0/m);
	vt = vi;
	ut = e_t.crossPdt(vt);
	matrix Gamma = GammaMtx(up_1, vp_1_i, ui, vi);
	matrixC Trmt = refractionMatrix02(1.0/m, sita_i, sita_t);



	double tmp1 = e_0 * e_t;
	double sita_s = acos(tmp1);
	if (isnan(sita_s)) throw std::runtime_error("Error #6! The domain of the \"acos()\" function is [-1,1].");
	

	double tmp2 = sin(sita_s);
	vec3 vs;
	if (abs(tmp2) < eps) vs = v_0;
	else vs = e_0.crossPdt(e_t) / tmp2;
	vec3 us = e_t.crossPdt(vs);
	matrix Gamma_s = GammaMtx(ut, vt, us, vs);

	vec3 u_tmp = e_0.crossPdt(vs);
	matrix Gamma_i = GammaMtx(u_tmp, vs, u_0, v_0);
	
	Jones_T = Gamma_s * Trmt * Gamma * Jones_R* Gamma_i;
	//Jones_T = Gamma_s * Trmt * Gamma * Jones_R;//this is maybe wrong, just for test
	if (isnan(norm(Jones_T)))
	{
		cout << "-nan(ind), NaN,Not a Number,indeterminate detected!" << endl;
		Gamma_s.show();
		Trmt.show();
		Gamma.show();
		Gamma_i.show();

	}

	vec3UnityCheck(vi);
	vec3UnityCheck(ui);
	vec3UnityCheck(vs);
	vec3UnityCheck(us);
	LawConsEnergyVerify(__FUNCTION__, sita_i, sita_t);

	return e_t;//note! carefully check! here we use e_i to return the e_r, so before using this function, the e_i MUST be copied to another vec3;

}
/* matrixC realToComplexMatrix(matrix& M_real)//convert a real value 2*2 matrix to complex value by adding i*0.0;
{
	matrixC M_complex(2,2);
	M_complex[0][0].re = M_real[0][0], M_complex[0][0].im = 0.0;
	M_complex[0][1].re = M_real[0][1]; M_complex[0][1].im = 0.0;
	M_complex[1][0].re = M_real[1][0]; M_complex[1][0].im = 0.0;
	M_complex[1][1].re = M_real[1][1]; M_complex[1][1].im = 0.0;
	return M_complex;
	
}
*/


//==============================================================================
// deriving 4*4 double Mueller matrix from a given 2*2 complex Jones matrix
// 
matrix Mueller(const matrixC& in)
{
	matrix M(4, 4);
	const double a11 = norm(in[0][0]), a12 = norm(in[0][1]),
		a21 = norm(in[1][0]), a22 = norm(in[1][1]);
	//----------------------------------------------------------------------------
	// first & second lines
	double A1 = a11 + a21, A2 = a12 + a22;
	M[0][0] = (A1 + A2) / 2.0; 	M[0][1] = (A1 - A2) / 2.0;
	//----------------------------------------------------------------------------
	A1 = a11 - a21; 			A2 = a12 - a22;
	M[1][0] = (A1 + A2) / 2.0; 	M[1][1] = (A1 - A2) / 2.0;
	complex C1 = in[0][0] * conj(in[0][1]), C2 = in[1][1] * conj(in[1][0]);
	M[0][2] = -real(C1) - real(C2); 			M[0][3] = imag(C2) - imag(C1);
	M[1][2] = real(C2) - real(C1); 			M[1][3] = -imag(C1) - imag(C2);
	//----------------------------------------------------------------------------
	// forth & third lines
	C1 = in[0][0] * conj(in[1][0]); 			C2 = in[1][1] * conj(in[0][1]);
	M[2][0] = -real(C1) - real(C2); 			M[2][1] = real(C2) - real(C1);
	M[3][0] = imag(C1) - imag(C2); 			M[3][1] = imag(C2) + imag(C1);
	C1 = in[0][0] * conj(in[1][1]); 			C2 = in[0][1] * conj(in[1][0]);
	M[2][2] = real(C1) + real(C2); 			M[2][3] = imag(C1) - imag(C2);
	M[3][2] = -imag(C1) - imag(C2); 			M[3][3] = real(C1) - real(C2);
	return M;
}

//================================================================binning and photon collection
void bin_and_collect(vector<matrix>& MlrMtxVector, const vec3& e_out, const vec3& PhotonDirection, const matrix& Mueller_Mtx, const int binNumber)
{
	double CosAng = e_out * PhotonDirection;

	//*
	if ((fabs(CosAng) - 1.0) < eps)
	{
		if (CosAng > 1.0) CosAng = 1.0;
		if (CosAng < -1.0) CosAng = -1.0;
	}
	//*/
	
	double scattAngle = (180.0 / pi) * acos(CosAng);
	
	if (isnan(scattAngle)&& isinf(scattAngle))
	{
		cout << "Error in function bin_and_collect, CosAng=" << CosAng<< endl;
		throw std::runtime_error("Error #7! The domain of the \"acos()\" function is [-1,1].");
	}
		
	//if (scattAngle < 80.0) cout << "scattAngle=" << scattAngle << endl;//testing code 
	for (int bn = 0; bn < binNumber; bn++)
	{
		if (scattAngle <= (bn + 1) * 180.0 / binNumber)
		{
			MlrMtxVector[bn] = MlrMtxVector[bn] + Mueller_Mtx;
			break;//very important!!!
		}
		if (scattAngle > 180.0) cout << "Energy loss! Error in function bin_and_collect, scattAngle=" << scattAngle << endl;
	}
}

void bin_and_collect02(vector<matrix>& MlrMtxVector, const vec3& e_out, const vec3& PhotonDirection, const matrix& Mueller_Mtx, const int binNumber)
{
	double CosAng = e_out * PhotonDirection;

	//*
	if ((fabs(CosAng) - 1.0) < eps)
	{
		if (CosAng > 1.0) CosAng = 1.0;
		if (CosAng < -1.0) CosAng = -1.0;
	}
	//*/

	double scattAngle = (180.0 / pi) * acos(CosAng);

	if (isnan(scattAngle) && isinf(scattAngle))
	{
		cout << "Error in function bin_and_collect, CosAng=" << CosAng << endl;
		throw std::runtime_error("Error #7! The domain of the \"acos()\" function is [-1,1].");
	}

	//if (scattAngle < 80.0) cout << "scattAngle=" << scattAngle << endl;//testing code 
	for (int bn = 0; bn < binNumber+1; bn++)
	{
		if (scattAngle <= (bn + 0.5) * 180.0 / binNumber)
		{
			MlrMtxVector[bn] = MlrMtxVector[bn] + Mueller_Mtx;
			break;//very important!!!
		}
	}
}