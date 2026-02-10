#include"functionList.h"
int nonTotRefCounter, totReflCounter; //for main.cpp
using namespace consts;

vec3 isotropicUnitDirVec()
{
	vec3 isoUniVec;
	double sita = acos(2.0 * rn() - 1.0);//sita---[0,pi]
	double phi = 2.0 * pi * rn();//phi---[0,2pi]
	isoUniVec.x = sin(sita) * cos(phi);
	isoUniVec.y = sin(sita) * sin(phi);
	isoUniVec.z = cos(sita);
	return isoUniVec;
}


vec3 isotropicUnitDirVec(vec3& v_hat, vec3& u_hat)//(x,y,z)-->(v,u,e) 
{
	vec3 isoUniVec;
	double sita = acos(2.0 * rn() - 1.0);//sita---[0,pi]
	double phi = 2.0 * pi * rn();//phi---[0,2pi]
	isoUniVec.x = sin(sita) * cos(phi);
	isoUniVec.y = sin(sita) * sin(phi);
	isoUniVec.z = cos(sita);

	vec3 k(1.0, 0.0, 0.0);
	if (fabs(fabs(k * isoUniVec) - 1.0) < 1.0e-5)
	{
		vec3 k0(0.0, 1.0, 0.0);
		k = k0;
	}
	u_hat = k.crossPdt(isoUniVec).normalization();
	v_hat = -isoUniVec.crossPdt(u_hat);

	return isoUniVec;
}

vec3 isotropicUnitDirVec05(vec3& v_hat, vec3& u_hat)//(x,y,z)-->(v,u,e) 
{
	vec3 isoUniVec;
	double sita = acos(2.0 * rn() - 1.0);//sita---[0,pi]
	NaNCheck(sita);

	double phi = 2.0 * pi * rn();//phi---[0,2pi]
	isoUniVec.x = sin(sita) * cos(phi);
	isoUniVec.y = sin(sita) * sin(phi);
	isoUniVec.z = cos(sita);

	vec3 k = isotropicUnitDirVec();
	while (fabs(fabs(k * isoUniVec) - 1.0) < eps1)
	{
		k = isotropicUnitDirVec();
	}
	u_hat = (k.crossPdt(isoUniVec)).normalization();
	v_hat = -isoUniVec.crossPdt(u_hat);

	vec3UnityCheck(u_hat);
	vec3UnityCheck(v_hat);
	return isoUniVec;
}

vec3 isotropicUnitDirVec01(vec3& v_hat, vec3& u_hat)//(x,y,z)-->(v,u,e) 
{
	vec3 isoUniVec;
	double sita = acos(2.0 * rn() - 1.0);//sita---[0,pi]
	double phi = 2.0 * pi * rn();//phi---[0,2pi]
	isoUniVec.x = sin(sita) * cos(phi);
	isoUniVec.y = sin(sita) * sin(phi);
	isoUniVec.z = cos(sita);

	vec3 k(1.0, 0.0, 0.0);

	v_hat = Z_sita(phi) * Y_sita(sita) *k;
	u_hat = isoUniVec.crossPdt(v_hat);
	
	return isoUniVec;
}
//Mikhailov, Gennadii Alekseevich, and Anton Vatslavovich Voitisek. Numerical Statistical Modeling: Monte Carlo Methods. 2006.
//You can find the relevant information on pages 91¨C92.
vec3 isotropicUnitDirVec02(vec3& v_hat, vec3& u_hat)//(x,y,z)-->(v,u,e) 
{
	vec3 isoUniVec;
	//double alpha1 = rn();
	//double alpha2 = rn();
	double alpha1 = RNG1();
	double alpha2 = RNG1();
	double tmp;
	isoUniVec.x = 1.0 - 2.0 * alpha1; tmp = sqrt(1.0 - SQR(isoUniVec.x));
	isoUniVec.y = tmp * cos(2.0 * pi * alpha2);
	isoUniVec.z = tmp * sin(2.0 * pi * alpha2);

	vec3 k(1.0, 0.0, 0.0);
	if (fabs(fabs(k * isoUniVec) - 1.0) < 1.0e-5)
	{
		vec3 k0(0.0, 1.0, 0.0);
		k = k0;
	}
	u_hat = k.crossPdt(isoUniVec).normalization();
	v_hat = -isoUniVec.crossPdt(u_hat);

	return isoUniVec;
}
//non random, absolute uniformly distributed 
vec3 isotropicUnitDirVec03(int positionNumber, vec3& v_hat, vec3& u_hat)//(x,y,z)-->(v,u,e) 
{
	int SqrtPositionNumber = static_cast<int>(sqrt(positionNumber));
	vec3 isoUniVec;
	double sita, phi;
	for (static int counter_i = 0; counter_i < SqrtPositionNumber; counter_i++)
	{
		for (static int counter_j = 0; counter_j < SqrtPositionNumber; counter_j++)
		{
		    sita = acos(1.0 - 2.0 * (counter_i + 0.5) / SqrtPositionNumber);//sita---[0,pi]
		    phi = 2.0 * pi * ((counter_j + 0.5) / SqrtPositionNumber);//phi---[0,2pi]
		}
	}	

	isoUniVec.x = sin(sita) * cos(phi);
	isoUniVec.y = sin(sita) * sin(phi);
	isoUniVec.z = cos(sita);

	vec3 k(1.0, 0.0, 0.0);
	v_hat = Z_sita(phi) * Y_sita(sita) * k;
	u_hat = isoUniVec.crossPdt(v_hat);

	/*
	vec3 k= isotropicUnitDirVec();
	while (fabs(fabs(k * isoUniVec) - 1.0) < 1.0e-5)
	{
		k = isotropicUnitDirVec();
	}
	u_hat = k.crossPdt(isoUniVec).normalization();
	v_hat = -isoUniVec.crossPdt(u_hat);
	//*/

	return isoUniVec;
}


vec3 isotropicUnitDirVec04(int positionNumber,double &sita, double& phi, vec3& v_hat, vec3& u_hat)//(x,y,z)-->(v,u,e) 
{	
	int SqrtPositionNumber = static_cast<int>(sqrt(positionNumber));
	vec3 isoUniVec;
	for (static int counter_i = 0; counter_i < SqrtPositionNumber; counter_i++)
	{
		for (static int counter_j = 0; counter_j < SqrtPositionNumber; counter_j++)
		{
			sita = acos(1.0 - 2.0 * (counter_i + 0.5) / SqrtPositionNumber);//sita---[0,pi]
			phi = 2.0 * pi * ((counter_j + 0.5) / SqrtPositionNumber);//phi---[0,2pi]
		}
	}

	isoUniVec.x = sin(sita) * cos(phi);
	isoUniVec.y = sin(sita) * sin(phi);
	isoUniVec.z = cos(sita);

	vec3 k(1.0, 0.0, 0.0);
	v_hat = Z_sita(phi) * Y_sita(sita) * k;
	u_hat = isoUniVec.crossPdt(v_hat);
	
	return isoUniVec;
}

vec3 initialPhotonPosition01(double R, const double& sita, const double& phi,  vec3& isoUniVec)//(x,y,z)-->(v,u,e) 
{
	
	vec3 p(-R + 2.0 * R * rn(), -R + 2.0 * R * rn(), -2.0 * R);
	while (p.getLength() > R)
	{
		p.x = -R + 2.0 * R * rn();
		p.y = -R + 2.0 * R * rn();
		
	}
	vec3 initialPhotonPosition = Z_sita(phi) *Y_sita(sita) * p;

	Plane3D Pl(-2.0 * R * isoUniVec, isoUniVec);
	if (Pl.containsPoint(initialPhotonPosition) == false)	
	{		
		cout << "distanceToPoint(p)=" << Pl.distanceToPoint(initialPhotonPosition) << endl;
		throw std::runtime_error("error in function initialPhotonPosition01!");
	}
	return initialPhotonPosition;
}




vec3 initialPhotonPosition(const vec3& v_hat, const vec3& u_hat, const vec3& initialPhotonDirection, double R)//(x,y,z)-->(v,u,e) 
{

	double s = -R + 2.0 * R * rn();
	double t = -R + 2.0 * R * rn();
	//double s = -R + 2.0 * R * RNG1();
	//double t = -R + 2.0 * R * RNG1();
	///in plane equation Ax+By+Cz+D=0, in our case, we set the center of ellipsoid is the origin point (0,0,0),  
	///and we set the distance between origin and initial plane, from which we sample initial position of photon, to be 2*R, 
	/// where R=ConHull_max_radius = maxOf3Num(semi_axis_a, semi_axis_b, semi_axis_c);
	vec3 P_0 = -2.0 * R * initialPhotonDirection;
	vec3 iPP = P_0 + s * u_hat + t * v_hat;
	//cout << fabs((iPP - P_0) * initialPhotonDirection) << endl;

	while ((iPP - P_0).getLength() > R)
		//Only points inside the circle will be selected;
		//points outside the circle will be rejected and require resampling
	{
		s = -R + 2.0 * R * rn();
		t = -R + 2.0 * R * rn();
		iPP = P_0 + s * u_hat + t * v_hat;
	}
	return iPP;

}

//this fuction is ONLY used to sample initial position of photon in the case of hexagonal crystal!!!
//Here, the initial position sampling of photons can be adjusted based on the orientation of ice crystal particles, 
// thereby narrowing the sampling range and improving the program's execution speed.
vec3 initialPhotonPosition_02(double gamma, double beta, double r_base, double R_max)
{
	double h = R_max * sin(beta);
	double r = r_base;
	double s = -r + 2.0 * r * rn();
	double t = -h + 2.0 * h * rn();
	vec3 e_axis(cos(gamma), sin(gamma), 0.0);//main axis of the hexagonal crystal
	vec3 e_perpendicular(-sin(gamma), cos(gamma), 0.0);//perpendicular to the main axis of the hexagonal crystal


	
	vec3 iPP = (s+t) * e_axis + s * e_perpendicular;
	
	return iPP;

}

vec3 initialPhotonPosition_03(double beta, double r_base, double R_max)
{
	double h = R_max * sin(beta);
	double r = r_base;
	double s = -r + 2.0 * r * rn();
	double t = -h + 2.0 * h * rn();
	vec3 e_axis(1.0, 0.0, 0.0);//main axis of the hexagonal crystal
	vec3 e_perpendicular(0.0, 1.0, 0.0);//perpendicular to the main axis of the hexagonal crystal



	vec3 iPP = (s + t) * e_axis + s * e_perpendicular;

	return iPP;

}




bool matrixJonesNonZeroCheckImpl(matrixC& M, const char* varName)
{
	matrixC Jones0(2, 2); Jones0.Fill(zeroComplex());
	if (M == Jones0) {
		cout << "Error in matrix: " << varName << ", |" << varName << "| = ";
		M.show();				
		return true;
	}
	return false;
}


void normalUpdateForSphere(vec3& normal, line& ray, vec3& intersP, const vec3& e_r, const sphere& sphericalParticle)
{
	//cout << "ray:" << endl;
	//ray.point_p.show();
	//ray.direction_d.show();

	ray=ray.update(intersP, e_r);

	//cout << "ray after update:" << endl;
	//ray.point_p.show();
	//ray.direction_d.show();

	intersP = findIntersecPointOfLineAndSphere(sphericalParticle, ray);	
	normal = intersP / sphericalParticle.radius;//This is correct just for spherical model!	
}
void normalUpdateForEllipsoid(vec3& normal, line& ray, vec3& intersP, const vec3& e_r, const ellipsoid& ellipsoidParticle)
{
	//cout << "ray:" << endl;
	//ray.point_p.show();
	//ray.direction_d.show();

	ray = ray.update(intersP, e_r);

	//cout << "ray after update:" << endl;
	//ray.point_p.show();
	//ray.direction_d.show();

	intersP = findIntersecPointOfLineAndEllipsoid(ellipsoidParticle, ray);
	normal = normalEllipsoid(ellipsoidParticle, intersP);
}


void normalUpdateForPrism(vec3& normal, line& ray, IntersectionInfo & hit, const vec3& e_r, const NPrism& prism)
{	
	
	ray = ray.update(hit.point, e_r);
	//auto nextIntersections = prism.findNextLineIntersections(ray,hit.face);
	auto nextIntersections = prism.findNextLineIntersections2004(ray, hit.face);
	if (nextIntersections.size() < 1)
	{
		cout << "error found in function normalUpdateForPrism" << endl;
	}
	hit = nextIntersections[0];	
	normal = hit.face->getNormal();	
}



void ReflThenRefrGoIn(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp, matrixC& Jones_refract_go_in, matrixC& Jones_ExtenalReflection_go_out, const double m)
{
	vec3 v1, u1;
	vec3 vr, ur;
	vec3 vs, us, e_s;
	double CosSita_i = -normal * e_incident;

	double sita_i = acos(CosSita_i);
	if (isnan(sita_i)) throw std::runtime_error("Error #1! The domain of the \"acos()\" function is [-1,1].");
	if (abs(sin(sita_i)) < eps) { v1 = v_0; }//throw std::runtime_error("Error #1.2! Division by zero!");
	else v1 = e_incident.crossPdt(normal) / sin(sita_i);


	vp = v1;//tracing vector vp_1
	u1 = e_incident.crossPdt(v1).normalization();
	matrix Gamma1 = GammaMtx(u_0, v_0, u1, v1);

	double tmp3 = sin(sita_i) / m;
	if (abs(tmp3) > 1.0 + eps) throw std::runtime_error("Error #1.1! The domain of the \"acos()\" function is [-1,1].");
	double sita_t = asin(tmp3);
	//e_r = e_incident + 2.0 * cos(sita_i) * normal;
	e_r = reflect(__FUNCTION__, normal, e_incident);
	//e_t = e_1 / m + (cos(sita_i) / m - cos(sita_t)) * normal;
	e_t = refract(__FUNCTION__, normal, e_incident, m);
	vr = v1; ur = e_r.crossPdt(vr);
	matrixC Refl1 = reflectionMatrix(m, sita_i, sita_t);

	double tmp1 = e_incident * e_r;
	if (abs(tmp1) > 1.0 + eps) throw std::runtime_error("Error #2! The domain of the \"acos()\" function is [-1,1].");
	double sita_s = acos(tmp1);

	double tmp2 = sin(sita_s);
	if (abs(tmp2) < eps1) vs = v_0;
	else vs = e_incident.crossPdt(e_r) / tmp2;
	us = e_r.crossPdt(vs).normalization();
	matrix Gamma_s = GammaMtx(ur, vr, us, vs);

	vec3 u_tmp = e_0.crossPdt(vs).normalization();
	matrix Gamma_1i = GammaMtx(u_tmp, vs, u_0, v_0);
	Jones_ExtenalReflection_go_out = Gamma_s * Refl1 * Gamma1 * Gamma_1i;

	matrixC Trmt1 = refractionMatrix01(m, sita_i, sita_t);
	Jones_refract_go_in = Trmt1 * Gamma1;

	/*
	vec3UnityCheck(v1);
	vec3UnityCheck(u1);
	vec3UnityCheck(us);
	vec3UnityCheck(vs);
	vec3UnityCheck(u_tmp);
	vec3UnityCheck(vs);
	*/
	matrixJonesNonZeroCheck(Jones_ExtenalReflection_go_out);
	matrixJonesNonZeroCheck(Jones_refract_go_in);
	LawConsEnergyVerify(__FUNCTION__, sita_i, sita_t);
}

void ReflInsideThenRefrGoOut(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp_1, matrixC& Jones_i, matrixC& Jones_refract_go_out, const double& m)
{
	
	vec3 vp_1_i = vp_1;
	vec3 up_1 = e_incident.crossPdt(vp_1);
	vec3 v_p, u_p, e_p; e_p = e_incident;
	vec3 vr, ur;
	vec3 minusNormal = -normal;//note: all other formulas assumed that all normal vector direction is facing the incoming rays!!!

	double tmp = -minusNormal * e_incident;
	double sita_i = acos(tmp);
	NaNCheck(sita_i);

	if (fabs(sin(sita_i)) < eps1)
	{
		v_p = vp_1;

	}

	else {
		v_p = e_incident.crossPdt(minusNormal) / sin(sita_i);
		//cout << "|v_p|=" << v_p.getLength()<<"  "; v_p.show();
	}
	vp_1 = v_p;
	u_p = e_p.crossPdt(v_p).normalization();
	//vec3UnityCheck(v_p);
	//vec3UnityCheck(u_p);
	//e_r = ep_1 + 2.0 * cos(sita_i) * minusNormal;

	//cout << "e_incident: "; e_incident.show();
	e_r = reflect(__FUNCTION__, -normal, e_incident);
	//cout << "e_incident after reflect(): "; e_incident.show();


	vr = v_p;
	ur = e_r.crossPdt(vr);
	
	matrix Gamma = GammaMtx(up_1, vp_1_i, u_p, v_p);

	double tmp1 = sin(sita_i) * m;
	matrixC Refl_p(2, 2);

	if (fabs(tmp1) > 1.0) //TIR case! refracted ray should not be expected!
	{
		//totalRefl = true;
		totReflCounter++;//added 20260203
		Refl_p = totalReflectionMatrix(1.0 / m, sita_i);
		Jones_i = Refl_p * Gamma * Jones_i;
		matrixC zeroMatrix(2, 2); zeroMatrix.Fill(zeroComplex());
		Jones_refract_go_out = zeroMatrix;
		//cout << "TIR case!" << endl;
		//Refl_p.Identity(); simply wrong!!!
	}
	else {
		nonTotRefCounter++;//added 20260203
		double sita_t = asin(tmp1);
		NaNCheck(sita_t);
		Refl_p = reflectionMatrix(1.0 / m, sita_i, sita_t);
		matrixC Trmt = refractionMatrix01(1.0 / m, sita_i, sita_t);
		LawConsEnergyVerify(__FUNCTION__, sita_i, sita_t);
		
		//-------------------------------------------------------------reflection ^
		vec3 vt, ut;
		e_t = refract(__FUNCTION__, -normal, e_incident, 1.0 / m);
		vt = v_p;
		ut = e_t.crossPdt(vt);
		double tmp = e_0 * e_t;

		if ((fabs(tmp) - 1.0) < eps)
		{
			if (tmp > 1.0) tmp = 1.0;
			if (tmp < -1.0) tmp = -1.0;
		}

		double sita_s = acos(tmp);
		if (isnan(sita_s))
		{
			cout << "e_0 * e_t=" << e_0 * e_t<<endl;
			throw std::runtime_error("Error #6.1! The domain of the \"acos()\" function is [-1,1].");
		}
		double tmp2 = sin(sita_s);
		vec3 vs;
		if (fabs(tmp2) < eps1) vs = v_0;
		else vs = (e_0.crossPdt(e_t) / tmp2).normalization();

		vec3 us = e_t.crossPdt(vs).normalization();
		matrix Gamma_s = GammaMtx(ut, vt, us, vs);
		vec3 u_tmp = e_0.crossPdt(vs).normalization();
		//vec3UnityCheck(u_tmp);
		//vec3UnityCheck(vs);
		matrix Gamma_i = GammaMtx(u_tmp, vs, u_0, v_0);


		Jones_refract_go_out = Gamma_s * Trmt * Gamma * Jones_i * Gamma_i;
		if (matrixJonesNonZeroCheck(Jones_refract_go_out)) 
		{
			cout << "tmp2=" << tmp2 << endl;
			cout << "v_0:  ";  v_0.show();
			cout << "e_0:  ";  e_0.show();
			cout << "e_t:  ";  e_t.show();
			cout << "e_0.crossPdt(e_t):  ";  (e_0.crossPdt(e_t)).show();
			cout << "vs:  ";  vs.show();
			Gamma_s.show();
			Trmt.show();
			Jones_i.show();
			Gamma.show();
			Gamma_i.show();
			throw std::runtime_error("Error in function \"matrixJonesNonZeroCheck()\": this Jones matrix should not be zero-matrix.");
		};
		
		Jones_i = Refl_p * Gamma * Jones_i;//Caution: this line must be put after "Jones_T = Gamma_s * Trmt * Gamma * Jones_i * Gamma_i;"!
		//Jones_T = Gamma_s * Trmt * Gamma * Jones_R;//this is maybe wrong, just for test
		if (isnan(norm(Jones_refract_go_out)))
		{
			cout << "-nan(ind), NaN,Not a Number,indeterminate detected!" << endl;
			Gamma_s.show();
			Trmt.show();
			Gamma.show();
			Gamma_i.show();

		}

		//vec3UnityCheck(vs);
		//vec3UnityCheck(us);
	}
	matrixJonesNonZeroCheck(Jones_i);
	
}

void ReflThenRefrGoIn02(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp, matrixC& Jones_refract_go_in, matrixC& Jones_ExtenalReflection_go_out, const double m)
{
	vec3 v1, u1;
	vec3 vr, ur;
	vec3 vs, us, e_s;
	double CosSita_i = -normal * e_incident;

	double sita_i = acos(CosSita_i);
	if (isnan(sita_i)) throw std::runtime_error("Error #1! The domain of the \"acos()\" function is [-1,1].");
	if (abs(sin(sita_i)) < eps) { v1 = v_0; }//throw std::runtime_error("Error #1.2! Division by zero!");
	else v1 = e_incident.crossPdt(normal).normalization();


	vp = v1;//tracing vector vp_1
	u1 = e_incident.crossPdt(v1).normalization();
	matrix Gamma1 = GammaMtx(u_0, v_0, u1, v1);

	double tmp3 = sin(sita_i) / m;
	if (abs(tmp3) > 1.0 + eps) throw std::runtime_error("Error #1.1! The domain of the \"acos()\" function is [-1,1].");
	double sita_t = asin(tmp3);
	//e_r = e_incident + 2.0 * cos(sita_i) * normal;
	e_r = reflect(__FUNCTION__, normal, e_incident);
	//e_t = e_1 / m + (cos(sita_i) / m - cos(sita_t)) * normal;
	e_t = refract(__FUNCTION__, normal, e_incident, m);
	vr = v1; ur = e_r.crossPdt(vr).normalization();
	matrixC Refl1 = reflectionMatrix(m, sita_i, sita_t);

	double tmp1 = e_incident * e_r;
	if (abs(tmp1) > 1.0 + eps) throw std::runtime_error("Error #2! The domain of the \"acos()\" function is [-1,1].");
	double sita_s = acos(tmp1);

	double tmp2 = sin(sita_s);
	if (abs(tmp2) < eps) vs = v_0;
	else vs = e_incident.crossPdt(e_r).normalization();
	us = e_r.crossPdt(vs).normalization();
	matrix Gamma_s = GammaMtx(ur, vr, us, vs);

	vec3 u_tmp = e_0.crossPdt(vs).normalization();
	matrix Gamma_1i = GammaMtx(u_tmp, vs, u_0, v_0);
	Jones_ExtenalReflection_go_out = Gamma_s * Refl1 * Gamma1 * Gamma_1i;

	matrixC Trmt1 = refractionMatrix02(m, sita_i, sita_t);
	Jones_refract_go_in = Trmt1 * Gamma1;

	//vec3UnityCheck(v1);
	//vec3UnityCheck(u1);
	//matrixJonesNonZeroCheck(Jones_ExtenalReflection_go_out);
	//matrixJonesNonZeroCheck(Jones_refract_go_in);
	LawConsEnergyVerify(__FUNCTION__, sita_i, sita_t);
}

void ReflInsideThenRefrGoOut02(const vec3& v_0, const vec3& u_0, const vec3& e_0, const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal, vec3& vp_1, matrixC& Jones_i, matrixC& Jones_refract_go_out, const double& m)
{

	vec3 vp_1_i = vp_1;
	vec3 up_1 = e_incident.crossPdt(vp_1).normalization();
	vec3 v_p, u_p, e_p; e_p = e_incident;
	vec3 vr, ur;
	vec3 minusNormal = -normal;//note: all other formulas assumed that all normal vector direction is facing the incoming rays!!!

	double tmp = -minusNormal * e_incident;
	double sita_i = acos(tmp);
	NaNCheck(sita_i);

	if (fabs(sin(sita_i)) < eps)
	{
		v_p = vp_1;

	}

	else {
		v_p = e_incident.crossPdt(minusNormal).normalization();
		//cout << "|v_p|=" << v_p.getLength()<<"  "; v_p.show();
	}
	vp_1 = v_p;
	u_p = e_p.crossPdt(v_p).normalization();
	//vec3UnityCheck(v_p);
	//vec3UnityCheck(u_p);
	//e_r = ep_1 + 2.0 * cos(sita_i) * minusNormal;

	//cout << "e_incident: "; e_incident.show();
	e_r = reflect(__FUNCTION__, -normal, e_incident);
	//cout << "e_incident after reflect(): "; e_incident.show();


	vr = v_p;
	ur = e_r.crossPdt(vr).normalization();

	matrix Gamma = GammaMtx(up_1, vp_1_i, u_p, v_p);

	double tmp1 = sin(sita_i) * m;
	matrixC Refl_p(2, 2);

	if (fabs(tmp1) > 1.0) //TIR case! refracted ray should not be expected!
	{
		//totalRefl = true;
		totReflCounter++;//added 20260203
		Refl_p = totalReflectionMatrix(1.0 / m, sita_i);
		Jones_i = Refl_p * Gamma * Jones_i;
		matrixC zeroMatrix(2, 2); zeroMatrix.Fill(zeroComplex());
		Jones_refract_go_out = zeroMatrix;
		//cout << "TIR case!" << endl;
		//Refl_p.Identity(); simply wrong!!!
	}
	else {
		nonTotRefCounter++;//added 20260203
		double sita_t = asin(tmp1);
		NaNCheck(sita_t);
		Refl_p = reflectionMatrix(1.0 / m, sita_i, sita_t);
		matrixC Trmt = refractionMatrix02(1.0 / m, sita_i, sita_t);
		LawConsEnergyVerify(__FUNCTION__, sita_i, sita_t);

		//-------------------------------------------------------------reflection ^
		vec3 vt, ut;
		e_t = refract(__FUNCTION__, -normal, e_incident, 1.0 / m);
		vt = v_p;
		ut = e_t.crossPdt(vt).normalization();
		double tmp = e_0 * e_t;

		if ((fabs(tmp) - 1.0) < eps)
		{
			if (tmp > 1.0) tmp = 1.0;
			if (tmp < -1.0) tmp = -1.0;
		}

		double sita_s = acos(tmp);
		if (isnan(sita_s))
		{
			cout << "e_0 * e_t=" << e_0 * e_t << endl;
			throw std::runtime_error("Error #6.1! The domain of the \"acos()\" function is [-1,1].");
		}
		double tmp2 = sin(sita_s);
		vec3 vs;
		if (fabs(tmp2) < eps) vs = v_0;
		else vs = e_0.crossPdt(e_t).normalization();

		vec3 us = e_t.crossPdt(vs).normalization();
		matrix Gamma_s = GammaMtx(ut, vt, us, vs);
		vec3 u_tmp = e_0.crossPdt(vs).normalization();
		matrix Gamma_i = GammaMtx(u_tmp, vs, u_0, v_0);
		Jones_refract_go_out = Gamma_s * Trmt * Gamma * Jones_i * Gamma_i;
		if (matrixJonesNonZeroCheck(Jones_refract_go_out))
		{
			cout << "tmp2=" << tmp2 << endl;
			cout << "v_0:  ";  v_0.show();
			cout << "e_0:  ";  e_0.show();
			cout << "e_t:  ";  e_t.show();
			cout << "e_0.crossPdt(e_t):  ";  (e_0.crossPdt(e_t)).show();
			cout << "vs:  ";  vs.show();
			Gamma_s.show();
			Trmt.show();
			Jones_i.show();
			Gamma.show();
			Gamma_i.show();
			throw std::runtime_error("Error in function \"matrixJonesNonZeroCheck()\": this Jones matrix should not be zero-matrix.");
		};

		Jones_i = Refl_p * Gamma * Jones_i;//Caution: this line must be put after "Jones_T = Gamma_s * Trmt * Gamma * Jones_i * Gamma_i;"!
		//Jones_T = Gamma_s * Trmt * Gamma * Jones_R;//this is maybe wrong, just for test
		if (isnan(norm(Jones_refract_go_out)))
		{
			cout << "-nan(ind), NaN,Not a Number,indeterminate detected!" << endl;
			Gamma_s.show();
			Trmt.show();
			Gamma.show();
			Gamma_i.show();

		}

		//vec3UnityCheck(vs);
		//vec3UnityCheck(us);
		LawConsEnergyVerify(__FUNCTION__, sita_i, sita_t);
	}
	//matrixJonesNonZeroCheck(Jones_i);

}

void ReflThenRefrGoIn03(const vec3& v_0, const vec3& u_0, const vec3& e_0,
	const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal,
	vec3& vp, matrixC& Jones_refract_go_in,
	matrixC& Jones_ExtenalReflection_go_out, const complex m)
{
	vec3 v1, u1;
	vec3 vr, ur;
	vec3 vs, us;

	// ======= using only cos¦È =======
	double cosI = -normal * e_incident;
	cosI = clamp(cosI, -1.0, 1.0);

	double sin2I = std::max(0.0, 1.0 - cosI * cosI);

	if (sin2I < eps)
		v1 = v_0;
	else
		v1 = e_incident.crossPdt(normal).normalization();

	vp = v1;
	u1 = e_incident.crossPdt(v1).normalization();

	matrix Gamma1 = GammaMtx(u_0, v_0, u1, v1);

	
	double sin2T = sin2I / (m.re * m.re);
	double cosT2 = sqrt(std::max(0.0, 1.0 - sin2T));

	
	e_r = reflect(__FUNCTION__, normal, e_incident);
	e_t = refract(__FUNCTION__, normal, e_incident, m.re);

	vr = v1;
	ur = e_r.crossPdt(vr).normalization();

	// ======= Fresnel£¨using only cos£©=======
	//matrixC Refl1 = reflectionMatrix(m, acos(cosI), acos(cosT)); //£¡£¡£¡wrong
	matrixC Refl1 = reflectionMatrix(m, acos(cosI), acos(cosT2));

	// ======= £¨no acos/sin£©=======
	vec3 cross_ir = e_incident.crossPdt(e_r);
	if (cross_ir.getLength() < eps)
		vs = v_0;
	else
		vs = cross_ir.normalization();

	us = e_r.crossPdt(vs).normalization();

	matrix Gamma_s = GammaMtx(ur, vr, us, vs);

	vec3 u_tmp = e_0.crossPdt(vs).normalization();
	matrix Gamma_1i = GammaMtx(u_tmp, vs, u_0, v_0);

	Jones_ExtenalReflection_go_out =
		Gamma_s * Refl1 * Gamma1 * Gamma_1i;

	matrixC Trmt1 = refractionMatrix02(m, acos(cosI), acos(cosT2));
	//matrixC Trmt1 = refractionMatrix01(m, acos(cosI), acos(cosT)); //see 2006 Yang eq.(2.27) may not work
	Jones_refract_go_in = Trmt1 * Gamma1;
	LawConsEnergyVerify(__FUNCTION__, sqrt(sin2I), asin(sqrt(sin2I)/m.re));
	CheckGammaMtx(Gamma1);
	CheckGammaMtx(Gamma_s);
	CheckGammaMtx(Gamma_1i);
}


void ReflInsideThenRefrGoOut03(const vec3& v_0, const vec3& u_0, const vec3& e_0,
	const vec3& e_incident, vec3& e_r, vec3& e_t, const vec3& normal,
	vec3& vp_1, matrixC& Jones_i,
	matrixC& Jones_refract_go_out, const complex& m)
{
	vec3 vp_1_i = vp_1;
	vec3 up_1 = e_incident.crossPdt(vp_1).normalization();

	vec3 v_p, u_p;
	vec3 vr, ur;

	vec3 minusNormal = -normal;

	// ======= using only cos¦È =======
	double cosI = -minusNormal * e_incident;
	cosI = clamp(cosI, -1.0, 1.0);

	double sin2I = std::max(0.0, 1.0 - cosI * cosI);

	if (sin2I < eps)
		v_p = vp_1;
	else
		v_p = e_incident.crossPdt(minusNormal).normalization();

	vp_1 = v_p;
	u_p = e_incident.crossPdt(v_p).normalization();

	e_r = reflect(__FUNCTION__, -normal, e_incident);

	vr = v_p;
	ur = e_r.crossPdt(vr).normalization();

	matrix Gamma = GammaMtx(up_1, vp_1_i, u_p, v_p);
	CheckGammaMtx(Gamma);
	matrixC Refl_p(2, 2);

	// ======= TIR =======
	if (m.re * m.re * sin2I > 1.0-eps)
	{
		totReflCounter++;

		Refl_p = totalReflectionMatrix02(1.0 / m, acos(cosI));
		Jones_i = Refl_p * Gamma * Jones_i;

		matrixC zeroMatrix(2, 2);
		zeroMatrix.Fill(zeroComplex());
		Jones_refract_go_out = zeroMatrix;
		
	}
	else
	{
		nonTotRefCounter++;

		double sin2T = m.re * m.re * sin2I;
		double cosT2 = sqrt(max(0.0, 1.0 - sin2T));

		Refl_p = reflectionMatrix(1.0 / m, acos(cosI), acos(cosT2));
		matrixC Trmt = refractionMatrix02(1.0 / m, acos(cosI), acos(cosT2));

		vec3 vt, ut;

		e_t = refract(__FUNCTION__, -normal, e_incident, 1.0 / m.re);

		vt = v_p;
		ut = e_t.crossPdt(vt).normalization();

		vec3 cross_et = e_0.crossPdt(e_t);

		vec3 vs;
		if (cross_et.getLength() < eps)
			vs = v_0;
		else
			vs = cross_et.normalization();

		vec3 us = e_t.crossPdt(vs).normalization();

		matrix Gamma_s = GammaMtx(ut, vt, us, vs);
		vec3 u_tmp = e_0.crossPdt(vs).normalization();
		matrix Gamma_i = GammaMtx(u_tmp, vs, u_0, v_0);

		Jones_refract_go_out =
			Gamma_s * Trmt * Gamma * Jones_i * Gamma_i;

		Jones_i = Refl_p * Gamma * Jones_i;
		LawConsEnergyVerify(__FUNCTION__, sqrt(sin2I), asin(sqrt(sin2I) * m.re));
		CheckGammaMtx(Gamma_s);
		CheckGammaMtx(Gamma_i);
	}
}
