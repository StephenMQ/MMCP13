/*reference:
1.Beam-splitting algorithm by Konoshonkin, Russia
2.Yang, Ping, and Kuo-Nan Liou. "Light scattering and absorption by nonspherical ice crystals." 
  Light Scattering Reviews: Single and Multiple Light Scattering. Berlin, Heidelberg: Springer Berlin Heidelberg, 2006. 31-71.
3.Zhang Z, Yang P, Kattawar GW, Tsay SC, Baum BA, Hu Y, Heymsfield AJ, Reichardt J. Geometrical-optics solution to light scattering by droxtal ice crystals. 
  Appl Opt. 2004 Apr 20;43(12):2490-9. doi: 10.1364/ao.43.002490. PMID: 15119619.
  
 

*/
#include <iostream>
#include <string>
#include <fstream>
#include <cctype>
#include <math.h>
#include <time.h>
#include <random> //for random_device
#include <functional>// usage of std::bind requires #include <functional>
#include "vec3.h"
#include "namesp.h"
#include "convexhull.h"
#include "functionList.h"

#include "CHVizer08.h" //In this version, users can zoom in/out the visualization by scrolling the mouse wheel.
#include "convex_polyhedron.h"

#include "ConvexPolyhedronVisualizer03.h" //In this version, Added multi-threading functionality to 
    //enable simultaneous execution of convex hull visualization and scattering matrix computation.                 

#include <stdexcept>  //for Exception handling: runtime_error
#include"matrix.h"
#include"geometry.h"
#include"prism.h"
#include"ColorConsole.h"
#include"TimeEstimator.h"
#include"templateFunctions.h"
#include <iomanip>
#include<cmath>
#include<vector>


string Title = "*************************************************************************\
             \n Light Scattering by Non-Spherical Particles + Polarization.                 \
             \n Computation of Mueller matrix for convex polyhedron. \
	         \n (Convex Hull Construction + Merge If Triangular Faces On Same Plane.)\
             \n Geometric optics + Monte Carlo, Diffraction is not considered.       \
		     \n 2025.06 \
		     \n Mu Quan,  Email: mu.quan@foxmail.com \
		     \n*************************************************************************";


using namespace std;
using namespace consts;//this namespace is defined by myself



//=====================================================================
int
positionNumber,
hits_counter,//added 20260202
batchSize,
batchNumber,
totalNumPhoton,
binNumber,
reflInsideControl,
totReflControl

;
extern int nonTotRefCounter,totReflCounter;
int innerExceptionCounter = 0;

extern double semi_axis_a, semi_axis_b, semi_axis_c; //semi-axis 
extern double ConHull_max_radius;

double
scattAngle=0.0, //scattering angle in degree (0-180)
refractiveIndex_r, refractiveIndex_i,
gmNumber , //number of gamma angles
btNumber , //number of beta angles
alNumber  //number of alpha angles
;
complex refractiveIndex(0.0, 0.0); //complex refractive index of the particle


vec3 ct(0.0, 0.0, 0.0);

matrix Mueller_Mtx(4,4);

void printData(ostream& out0, const complex& refractiveIndex, int positionNumber, int batchSize, int batchNumber, int binNumber, int reflInsideControl, int totReflControl) {
	out0 << "The data read from the file is as follows:" << endl;
	out0 << "refractiveIndex: " << refractiveIndex.re << "+ i "<<refractiveIndex.im <<endl;
	out0 << "positionNumber: " << positionNumber << endl;
	out0 << "batchSize: " << batchSize << endl;
	out0 << "batchNumber: " << batchNumber << endl;
	out0 << "binNumber: " << binNumber << endl;
	out0 << "reflInsideControl: " << reflInsideControl << endl;
	out0 << "totReflControl: " << totReflControl << endl;
	out0 << "========================================================================== " << endl;
}

int main(int argc, char** argv)
{

	//----------------------convex hull construction and visualization begin---------------------
	vertex_index_vector.clear();
	convex_hull_vertices.clear();
	faces_normal.clear();
	final_faces.clear();

	//-------------------------------------case 1: generate points
	//generatePointsOnSphere(p);
	//generatePointsOnSphere02(p);
	//generateRandomPoints(p);
	//generatePointsOnEllipsoid(p);
	//generatePointsOnEllipsoidWithOnePole(p);
	//generatePointsOnEllipsoidWithoutPoles(p);
	generatePointsOnEllipticCylinder(p);

	//--------------------------------------case 2: read points from file
	//readIn(p);  outputReadInDataForCheck(p);

	//----------------------------------begin convex hull construction
	hull.make();
	//output information about the convex hull
	ofstream outfile02;
	outfile02.flags(ios::scientific);
	outfile02.open("Convex Hull.info");
	if (!outfile02.is_open()) {
		cerr << "Error opening output file!" << endl;
		return 1;
	}
	//Raw String Literal R"(...)"
	outfile02 << R"(=== Convex Hull Construction Result ===
NIP: Number of Initial points;
NFF: Number of Final Faces;
NFV: Number of Final Vertices;)" << endl;
	outfile02 << "NIP=" << hull.n << "   NFF=" << hull.final_faces_number << "   NFV=" << hull.final_vertices_number << endl;
	outfile02 << "semi_axis_a=" << semi_axis_a << "   semi_axis_b=" << semi_axis_b << "   semi_axis_c=" << semi_axis_c << endl;


	outfile02.close();


	//Convex_polyhedron ConvexPldron01;
	//ConvexHullVisualizer::run(argc, argv,convex_hull_vertices);
	//upper two lines are needed for before merging triangular faces on the same plane.
	
	//-------------------begin convex polyhedron construction and visualization
	Convex_polyhedron ConvexPldron01;
	ConvexPldron01 = MergeIfTriangularFacesOnSamePlane(convex_hull_vertices);
	//output information about the convex polyhedron
	ofstream outfile05;
	outfile05.flags(ios::scientific);
	outfile05.open("Convex Polyhedron.info");
	if (!outfile05.is_open()) {
		cerr << "Error opening output file!" << endl;
		return 1;
	}
	
	outfile05 << "NF(Number of Faces)=" << ConvexPldron01.CP_faces.size()
		<< ",   NV(Number of Vertices)=" << ConvexPldron01.CP_vertices.size() << endl;
	outfile05.close();
	
	
	ConvexPolyhedronVisualizer::run(argc, argv, ConvexPldron01.CP_faces);


	//----------------------convex hull construction and visualization end---------------------




	try{
		//vec3 test(0.0, 0.0, 1.0); vec3 PhotonDirectionTest(0.0, 0.0, 1.0); cout << "dot product=" << test * PhotonDirectionTest << endl;
		//cout << "acos=" << acos(test * PhotonDirectionTest) << endl;
	clock_t start = clock();
	//random_device rd;
	//mt19937 mt(rd());
	//uniform_real_distribution<double> rn(0.0, 1.0);
	cout << Title << endl;

	//=================================================================read parameters from file 'params.in'
	cout << endl;
	cout << "The particle size data can be modified in file <convexhull.cpp>, while other data is in the file <params.in>. " << endl;
	cout << "Reading parameters from file 'params.in'......" << endl;
	ifstream in;
	in.open("params.in");	
	in >> refractiveIndex.re;
	in >> refractiveIndex.im;
	in >> positionNumber;
	in >> batchSize;
	in >> batchNumber;
	in >> binNumber;
	in >> reflInsideControl;
	in >> totReflControl;
	in.close();
	
	//positionNumber = gmNumber* btNumber * alNumber; //total number of positions, this is used to control the number of Euler angles for the orientation of ice crystal particles.

	//positionNumber = (gmNumber + 1) * (btNumber + 1) * (alNumber + 1); //for start from 0
	
	totalNumPhoton = batchSize * batchNumber* positionNumber;
	TimeEstimator estimator(positionNumber);

	//For simplicity, we will continue to use the name "prism," though in this context it actually refers to a convex polyhedron.
	NPrism prism(ConvexPldron01.CP_faces);
	//prism.printPrismFacesInfo();//for test

	
	
	ofstream out; out<< scientific << setprecision(6);
	out.open("Info.out");
	out<< Title << endl;

	ofstream out02; out02 << scientific << setprecision(6);
	out02.open("Info02.out");
	out02 << Title << endl;

	ofstream out03; out03 << scientific << setprecision(6);
	out03.open("Info03.out");
	out03 << Title << endl;

	ofstream logFile; logFile.open("error_log.txt"); 
	logFile << Title << endl;
	/*
	//random number test
	ofstream outest; outest << scientific << setprecision(6);
	outest.open("Test.out");
	for (int i = 0; i < 100; i++)
	{
		outest << isotropicUnitDirVec().x << "   " << isotropicUnitDirVec().y << "   " << isotropicUnitDirVec().z  << endl;
	}
	
	//*/	

	printData(cout, refractiveIndex,  positionNumber, batchSize, batchNumber, binNumber, reflInsideControl, totReflControl);//output data to console
	printData(out, refractiveIndex, positionNumber, batchSize, batchNumber, binNumber, reflInsideControl, totReflControl);//output data to file
	printData(out02, refractiveIndex, positionNumber, batchSize, batchNumber, binNumber, reflInsideControl, totReflControl);//output data to file
	
	//matrixC polarizationState(2, 1);
	vector<matrix> MlrMtxVector(binNumber+1);//this matrix vector is needed to storage the Final Scattering Matrix.
	for (size_t i = 0; i < MlrMtxVector.size(); ++i) 
	{
		    MlrMtxVector[i] = Mueller_Mtx;//use matrix Mul to set the columns(m) and rows(n);
			MlrMtxVector[i].Fill(0.0); //initialize all the Mueller Matrix to Zero matrix;
	}
	
	//-------------------------------------------------------------ray tracing begin
	//here + for loop to control ray direction, if need;
	vec3 v_00, u_00;
	//vec3 PhotonDirection_00 = isotropicUnitDirVec05(v_00, u_00); //editted 20260202
	
	hits_counter = 0; //added 20260202
	for (int posiNum = 0; posiNum < positionNumber; posiNum++)
	//for (int gm = 1; gm <= gmNumber; gm++)
		//for (int gm = 0; gm <= gmNumber; gm++)
	{
			double gamma = 2.0 * pi * rn();
		//double gamma = 2.0 * pi * gm / gmNumber;
		//double gamma =  pi/6.0 * gm / gmNumber;
		//double gamma = pi * gm / gmNumber;
		//double gamma = 0.0;
		
		//for (int bt = 1; bt <= btNumber; bt++)
			//for (int bt = 0; bt <= btNumber; bt++)
	//	{
				double beta = acos(2.0 * rn() - 1.0);
			//double beta = acos(2.0 * bt / btNumber - 1.0); //beta is the angle between the incident light direction and the z-axis, beta in [0, pi]
			//double beta = acos(bt / btNumber);
			//for (int al = 1; al <= alNumber; al++)
				//for (int al = 0; al <= alNumber; al++)
			//{
					double alpha = 2.0 * pi * rn();
				//double alpha = 2.0 * pi * al / alNumber; //alpha is the angle between the x-axis and the projection of the incident light direction on the xy-plane, alpha in [0, 2pi]
				//double alpha = pi * al / alNumber;
				//double alpha = pi/3.0 * al / alNumber;

				static int pN = 0;
				pN++;

				//for (int pN = 1; pN <= positionNumber; pN++)
				{
					///rotate ray, not the crystal; the direction vector should be a Isotropic Random Vector!!!
					vec3 v_00(1.0, 0.0, 0.0), u_00(0.0, 1.0, 0.0);//(x,y,z)-->(v,u,e) 	
					vec3 PhotonDirection_00(0.0, 0.0, 1.0); //for fixed incident light direction, the orientation of ice crystal particles is controlled using Euler angles

					//vec3 PhotonDirection = isotropicUnitDirVec(v_0,u_0);
					//vec3 PhotonDirection = isotropicUnitDirVec02(v_0, u_0);
					//vec3 PhotonDirection = isotropicUnitDirVec05(v_0, u_0);
					//vec3 PhotonDirection = isotropicUnitDirVec03(positionNumber, v_0, u_0);
					//double sita, phi;
					//vec3 PhotonDirection = isotropicUnitDirVec04(positionNumber, sita, phi, v_0, u_0);

					//double gamma = 2.0 * pi * rn();
					//double beta = acos(2.0 * rn() - 1.0); //beta is the angle between the incident light direction and the z-axis, beta in [0, pi]
					//double alpha = 2.0 * pi * rn(); //alpha is the angle between the x-axis and the projection of the incident light direction on the xy-plane, alpha in [0, 2pi]

					prism.rotateCrystal(gamma, beta, alpha);//rotate the crystal using Euler angles, the incident light direction is fixed.

					//prism.rotateCrystal_BA(beta, alpha);

					vec3 e0 = PhotonDirection_00;
					vec3 PhotonDirection = PhotonDirection_00;
					vec3 v0 = v_00;
					vec3 u0 = u_00;
					for (int bN = 1; bN <= batchNumber; bN++)
					{
						for (int bs = 1; bs <= batchSize; bs++)
						{
							vec3 v_0 = v0;
							vec3 u_0 = u0;
							vec3 e_0 = e0;

							vec3 PhotonPosition = initialPhotonPosition(v_0, u_0, PhotonDirection, ConHull_max_radius);
							//vec3 PhotonPosition = initialPhotonPosition_03(beta, semi_axis_a, ConHull_max_radius);//note:Only works for hegagonal prism, not for other convex polyhedrons.
							//vec3 PhotonPosition = initialPhotonPosition_02(gamma, beta, semi_axis_a, ConHull_max_radius);//note:Only works for hegagonal prism, not for other convex polyhedrons.
							//vec3 PhotonPosition = initialPhotonPosition01(ConHull_max_radius, sita, phi, PhotonDirection);
							//polarizationState[0][0].re = 1.0; polarizationState[0][0].im = 0.0;
							//polarizationState[1][0].re = 0.0; polarizationState[1][0].im = 0.0;				

							line ray(PhotonPosition, PhotonDirection);
							//int counter = 0;
							//while (prism.doesLineIntersect(ray) == false)// re-sample photon
							//while (prism.findLineIntersections2004(ray).empty())
							if (!prism.findLineIntersections2004(ray).empty()) //2026.02.02 editted
							{
								hits_counter++;
								nonTotRefCounter = totReflCounter = 0; //added 20260202
								//vec3 PhotonPosition = initialPhotonPosition(v_0, u_0, PhotonDirection, ConHull_max_radius);//2026.02.02 editted
								//vec3 PhotonPosition = initialPhotonPosition01(ConHull_max_radius, sita, phi, PhotonDirection);
								//vec3 PhotonPosition = initialPhotonPosition_02(gamma, beta, semi_axis_a, ConHull_max_radius);
								//vec3 PhotonPosition = initialPhotonPosition_03(beta, semi_axis_a, ConHull_max_radius);

								//ray.point_p = PhotonPosition;

								/* //2026.02.02 editted
								counter++;
								if (counter >= 1000)
								{
									throw std::runtime_error("error, please check the \" while loop \"!");
									break;
								}
								*/


								//}

								//if (prism.doesLineIntersect(ray)) //if intersection happens, do the ray tracing process, else do nothing and re-sample photon
								try {
									matrixC Jones0(2, 2); Jones0.Fill(zeroComplex());
									matrixC Jones_in = Jones0;
									matrixC Jones_out = Jones0;

									//auto firstIntersections = prism.findFirstLineIntersections(ray);
									auto firstIntersections = prism.findLineIntersections2004(ray);

									if (firstIntersections.size() < 1)
									{
										cout << "error found in function firstIntersections" << endl;

									}
									IntersectionInfo Hit = firstIntersections[firstIntersections.size() - 1];//take max of all possible distances
									//vec3 firstintersP = firstHit.point;				
									vec3 n = Hit.face->getNormal();

									vec3 e_r, e_t, vp;

									ReflThenRefrGoIn(v_0, u_0, e_0, PhotonDirection, e_r, e_t, n, vp, Jones_in, Jones_out, refractiveIndex.re);

									Mueller_Mtx = Mueller(Jones_out);
									bin_and_collect02(MlrMtxVector, e_r, PhotonDirection, Mueller_Mtx, binNumber);

									normalUpdateForPrism(n, ray, Hit, e_t, prism);


									ReflInsideThenRefrGoOut(v_0, u_0, e_0, e_t, e_r, e_t, n, vp, Jones_in, Jones_out, refractiveIndex.re);
									Mueller_Mtx = Mueller(Jones_out);
									bin_and_collect02(MlrMtxVector, e_t, PhotonDirection, Mueller_Mtx, binNumber);

									//	normalUpdateForEllipsoid(n, ray, intersP, e_r, ellipsoidParticle);
									normalUpdateForPrism(n, ray, Hit, e_r, prism);
									//for (int reflectionInsideCount = 0; reflectionInsideCount < reflInsideControl; reflectionInsideCount++) //consider reflection times
									while (nonTotRefCounter < reflInsideControl && totReflCounter < totReflControl) //added 20260202
									{
										vec3 e_incident = e_r;

										ReflInsideThenRefrGoOut(v_0, u_0, e_0, e_incident, e_r, e_t, n, vp, Jones_in, Jones_out, refractiveIndex.re);

										Mueller_Mtx = Mueller(Jones_out);
										bin_and_collect02(MlrMtxVector, e_t, PhotonDirection, Mueller_Mtx, binNumber);

										//normalUpdateForEllipsoid(n, ray, intersP, e_r, ellipsoidParticle);
										normalUpdateForPrism(n, ray, Hit, e_r, prism);
									}

									//TotalPhotonNumberCounter++;

								}


								catch (const std::exception& e)
								{
									innerExceptionCounter++;
									logFile << "Inner exception caught: " << e.what() << endl;
								}
							}

						}
						
						cout << "\rSimulating the " << GREEN << setw(3) << bN << RESET << " -th batch (out of " << GREEN << batchNumber << RESET << ") for "
							<< GREEN << setw(5) << pN << RESET << " - th position (out of " << GREEN << positionNumber << RESET
							<< ") ......" << endl;


						cout << "Total simulation progress : " << GREEN << 100.0 * ((1.0 * (pN - 1.0) / positionNumber) + (1.0 * bN / batchNumber / positionNumber))
							<< RESET << "%  ......   " << endl;

						//if (pN % 100 == 0)
						estimator.update(pN);
						int seconds = estimator.get_remaining_time() / 1000.0;
						cout << "Remaining:  " << floor(seconds / 3600) << "  hours  " << floor((seconds % 3600) / 60) << "  minutes  ......   " << endl;
						cout << "\033[3F"; // Move the cursor up 3 lines
					}
					//showProgressBar(pN, positionNumber, 50, positionNumber / 100);



				}
			
			//}
		//}
	}

	
	
	//=================================================================normalization and output data
	
	/*
	double sum_Mlr_M11 = 0.0;
	for (int bn = 0; bn < binNumber; bn++)
	{
		sum_Mlr_M11 = sum_Mlr_M11 + MlrMtxVector[bn][0][0];// =dE
	}
	for (int bn = 0; bn < binNumber; bn++)
	{
		
		double SolAng = 2.0 * pi * (cos(bn * pi / binNumber) - cos((bn + 1) * pi / binNumber));//solid angle
		

		//MlrMtxVector[bn] = MlrMtxVector[bn] / totalNumPhoton / SolAng;
		// MlrMtxVector[bn] = MlrMtxVector[bn]/totalNumPhoton/SolAng; //
		//MlrMtxVector[bn] = MlrMtxVector[bn] * (4 * pi) / hits_counter / SolAng; //added 20260202
		MlrMtxVector[bn] = MlrMtxVector[bn] * (4 * pi) / sum_Mlr_M11 / SolAng; //added 20260202
		//MlrMtxVector[bn] = MlrMtxVector[bn] / sum_Mlr_M11 / SolAng*2; //added 20260205
	}
	*/

	//*
	//second normalization form from 20260205 to campare with Macke
	double sum_Mlr_M11 = 0.0;
	for (int bn = 0; bn < binNumber + 1; bn++)
	{
		sum_Mlr_M11 = sum_Mlr_M11 + MlrMtxVector[bn][0][0];// =dE
	}

	double help1 = pi / (binNumber);//??
	//double help2= pi / (binNumber)/2.0;
	for (int bn = 0; bn < binNumber+1; bn++)
	{
		
		double SolAng = 2.0 * pi * 2 * sin(bn * help1) * sin(help1 / 2.0); //solid angle interval
			//double SolAng = 2.0 * pi * (cos(bn * pi / binNumber) - cos((bn + 1) * pi / binNumber));//solid angle
		if (bn == 0)  SolAng = 2.0 * pi * (cos(0.0) - cos(0.5 * pi / binNumber));
		if (bn == binNumber) SolAng = 2.0 * pi * (cos(179.5 * pi / binNumber) - cos(180.0 * pi / binNumber));

		//MlrMtxVector[bn] = MlrMtxVector[bn] / totalNumPhoton / SolAng;
		// MlrMtxVector[bn] = MlrMtxVector[bn]/totalNumPhoton/SolAng; //
		//MlrMtxVector[bn] = MlrMtxVector[bn] * (4 * pi) / hits_counter / SolAng; //added 20260202
		MlrMtxVector[bn] = MlrMtxVector[bn] * (4 * pi) / sum_Mlr_M11 / SolAng; //added 20260205 from Macke
	}


	//*/
	

	

	clock_t end01 = clock();
	double run_time01 = static_cast<double>(end01 - start) / CLOCKS_PER_SEC;
	out << "Hits counter: " << hits_counter << "; Total incoming NumPhoton:"<< totalNumPhoton<<" (about "<< hits_counter * 100.0 / totalNumPhoton<<"%)." << endl; //added 20260202
	out << "sum_Mlr_M11: " << sum_Mlr_M11 << "; Missing:" << hits_counter - sum_Mlr_M11 << " (about " << (hits_counter - sum_Mlr_M11) * 100.0 / hits_counter << "%)." << endl;
	out << "Total run time: " << run_time01 << " seconds."<< endl;
	out << "========================================================================== " << endl;
	out << "Scattering Angle      F11              F12              F21              F22              F33              F34              F43              F44" << endl;//6+14 spaces
	for (int bn = 0; bn < binNumber; bn++)
	{
		out << showpos << bn * 180.0 / binNumber			
			<< "    " << MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][0][1]
			<< "    " << MlrMtxVector[bn][1][0]
			<< "    " << MlrMtxVector[bn][1][1]
			<< "    " << MlrMtxVector[bn][2][2]
			<< "    " << MlrMtxVector[bn][2][3]
			<< "    " << MlrMtxVector[bn][3][2]
			<< "    " << MlrMtxVector[bn][3][3]
			<< endl;

	}


	out02 << "Hits counter: " << hits_counter << "; Total incoming NumPhoton:" << totalNumPhoton << " (about " << hits_counter * 100.0 / totalNumPhoton  << "%)." << endl; //added 20260202
	out02 << "sum_Mlr_M11: " << sum_Mlr_M11 << "; Missing:" << hits_counter - sum_Mlr_M11 << " (about " << (hits_counter - sum_Mlr_M11) * 100.0 / hits_counter << "%)." << endl;
	out02 << "Total run time: " << run_time01 << " seconds." << endl;
	out02 << "========================================================================== " << endl;
	out02 << "Scattering Angle      F11           -F12/F11         -F21/F11          F22/F11          F33/F11          F34/F11          F43/F11          F44/F11" << endl;//6+12+10 spaces
	
	/*
	for (int bn = 0; bn < binNumber; bn++)
	{
		out02 << showpos << bn * 180.0 / binNumber
			<< "    " << MlrMtxVector[bn][0][0]
			<< "    " << -MlrMtxVector[bn][0][1] / MlrMtxVector[bn][0][0]
			<< "    " << -MlrMtxVector[bn][1][0] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][1][1] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][2][2] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][2][3] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][3][2] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][3][3] / MlrMtxVector[bn][0][0]
			<< endl;

	}
	*/

	//*
		for (int bn = 0; bn < binNumber+1; bn++)
	{
		if (bn == 0)

			out02 << showpos << (0.25) * 180.0 / binNumber
			<< "    " << MlrMtxVector[bn][0][0]
			<< "    " << -MlrMtxVector[bn][0][1]/ MlrMtxVector[bn][0][0]
			<< "    " << -MlrMtxVector[bn][1][0]/ MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][1][1]/ MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][2][2]/ MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][2][3]/ MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][3][2]/ MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][3][3]/ MlrMtxVector[bn][0][0]
			<< endl;
		else
		{
			if (bn < binNumber)
				out02 << showpos << (bn) * 180.0 / binNumber
				<< "    " << MlrMtxVector[bn][0][0]
				<< "    " << -MlrMtxVector[bn][0][1] / MlrMtxVector[bn][0][0]
				<< "    " << -MlrMtxVector[bn][1][0] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][1][1] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][2][2] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][2][3] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][3][2] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][3][3] / MlrMtxVector[bn][0][0]
				<< endl;
			else
				out02 << showpos << (180.0 - 0.25) * 180.0 / binNumber
				<< "    " << MlrMtxVector[bn][0][0]
				<< "    " << -MlrMtxVector[bn][0][1] / MlrMtxVector[bn][0][0]
				<< "    " << -MlrMtxVector[bn][1][0] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][1][1] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][2][2] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][2][3] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][3][2] / MlrMtxVector[bn][0][0]
				<< "    " << MlrMtxVector[bn][3][3] / MlrMtxVector[bn][0][0]
				<< endl;

		}
			

	}
	//*/

	out03 << "Hits counter: " << hits_counter << "; Total incoming NumPhoton:" << totalNumPhoton << " (about " << hits_counter * 100.0 / totalNumPhoton  << "%)." << endl; //added 20260202
	out03 << "sum_Mlr_M11: " << sum_Mlr_M11 << "; Missing:" << hits_counter - sum_Mlr_M11 << " (about " << (hits_counter - sum_Mlr_M11) * 100.0 / hits_counter << "%)." << endl;
	out03 << "Total run time: " << run_time01 << " seconds." << endl;
	out03 << "========================================================================== " << endl;
	out03 << "Scattering Angle      F11           -F12/F11         -F13/F11         -F14/F11         -F21/F11          F22/F11          F23/F11          F24/F11          F31/F11          F32/F11          F33/F11          F34/F11          F41/F11          F42/F11          F43/F11          F44/F11" << endl;//6+12+10 spaces
	for (int bn = 0; bn < binNumber; bn++)
	{
		out03 << showpos << bn * 180.0 / binNumber
			//-----------------------------------------------------------row 1
			<< "    " <<  MlrMtxVector[bn][0][0]
			<< "    " << -MlrMtxVector[bn][0][1] / MlrMtxVector[bn][0][0]
			<< "    " << -MlrMtxVector[bn][0][2] / MlrMtxVector[bn][0][0]
			<< "    " << -MlrMtxVector[bn][0][3] / MlrMtxVector[bn][0][0]

			//-----------------------------------------------------------row 2
			<< "    " << -MlrMtxVector[bn][1][0] / MlrMtxVector[bn][0][0]
			<< "    " <<  MlrMtxVector[bn][1][1] / MlrMtxVector[bn][0][0]
			<< "    " <<  MlrMtxVector[bn][1][2] / MlrMtxVector[bn][0][0]
			<< "    " <<  MlrMtxVector[bn][1][3] / MlrMtxVector[bn][0][0]

			//-----------------------------------------------------------row 3
			<< "    " << MlrMtxVector[bn][2][0] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][2][1] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][2][2] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][2][3] / MlrMtxVector[bn][0][0]

			//-----------------------------------------------------------row 4
			<< "    " << MlrMtxVector[bn][3][0] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][3][1] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][3][2] / MlrMtxVector[bn][0][0]
			<< "    " << MlrMtxVector[bn][3][3] / MlrMtxVector[bn][0][0]
			<< endl;

	}

	
	out.close();
	out02.close();
	out03.close();
	logFile.close();

	clock_t end = clock();
	double run_time = static_cast<double>(end - start) / CLOCKS_PER_SEC;

	cout << endl << endl << endl; // Move the cursor down 3 lines

	cout << "InnerExceptionCounter=" << innerExceptionCounter <<". See the log file (error_log.txt) for more details."<< endl;
	cout << GREEN<< "Done! Total run time: " << floor(run_time / 3600.0) 
		<< " hours " << floor((run_time - floor(run_time / 3600.0)*3600.0)/60.0) << " minutes." <<RESET<< endl;
	} catch (const std::exception& e) {
		std::cerr << endl << endl << "Exception caught: " << e.what() << std::endl;
	}

	cin.get();
	return 0;
}