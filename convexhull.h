#pragma once

#include <iostream>
#include <fstream>
#include <iomanip>  // for formatted output
#include <vector>
#include "vec3.h" 
#include"namesp.h"
#include "functionList.h"
#include <unordered_set>
#include <string>  
#include <sstream>  // must include, otherwise std::istringstream cannot be used

using namespace consts;
using namespace std;

const int M0 = 20000; // maximum number of points
void readInPointCoords(vec3 p[]);
void generateRandomPoints(vec3 p[]);
void generatePointsOnSphere(vec3 p[]);
void generatePointsOnSphere02(vec3 p[]);
void generatePointsOnEllipsoid(vec3 p[]);
void generatePointsOnEllipsoidWithOnePole(vec3 p[]);
void generatePointsOnEllipsoidWithoutPoles(vec3 p[]);
void generatePointsOnEllipticCylinder(vec3 p[]);
void readIn(vec3 p[]);

void outputReadInDataForCheck(vec3 p[]);
bool isPointInConvexHull(vec3& p);

struct face
{
	int a, b, c;  // indices of the three vertices of one convex hull face
	bool is_final_face;     // whether this face belongs to the final convex hull
};

extern vec3 p[M0];// all initial input points
extern face facesNumber[M0 * 8];// convex hull triangles
extern vector<int> vertex_index_vector; // indices of convex hull triangle vertices, store the three vertices of each face
extern vector<vec3> convex_hull_vertices; // convex hull vertex coordinates
extern vector<vec3> faces_normal; // normals of convex hull faces
extern vector<face> final_faces; // faces after hull construction
extern vec3 geometric_center;

struct Convex_hull//3D convex hull
{
	int n;// number of initial points	
	int conhull_triangle_index;// convex hull triangle index      
	int face_containing_edge[M0][M0];// which face contains edge i->j
	int final_faces_number; // number of final convex hull faces
	int final_vertices_number; // number of final convex hull vertices


	// tetrahedron volume * 6, i.e., scalar triple product
	double volume(vec3 face_a, vec3 face_b, vec3 face_c, vec3 point_p)
	{
		return (face_b - face_a).crossPdt(face_c - face_a) * (point_p - face_a);
	}

	double pointToFace(vec3 point_p, face f)     // positive: point and face have same orientation
	{
		vec3 m = p[f.b] - p[f.a];
		vec3 n = p[f.c] - p[f.a];
		vec3 t = point_p - p[f.a];
		return TriplePdt(m, n, t);
	}

	void deal(int q, int a, int b)
	{
		int edge_shared_face_index = face_containing_edge[a][b];      // the other face sharing edge ab with current face
		face add;
		if (facesNumber[edge_shared_face_index].is_final_face)         // if this face is still on the hull surface
		{
			if (pointToFace(p[q], facesNumber[edge_shared_face_index]) > eps1)
				deepdeal(q, edge_shared_face_index);// if point q sees this face, continue deep search on its three edges


			else            // if q sees current face but not this shared face, create a new triangle (p,a,b)
			{
				add.a = b;
				add.b = a;
				add.c = q;
				add.is_final_face = 1;
				face_containing_edge[b][a] = face_containing_edge[a][q] = face_containing_edge[q][b] = conhull_triangle_index;
				facesNumber[conhull_triangle_index++] = add;
			}
		}
	}

	// each triangle records edges counterclockwise, reversing the edge finds the adjacent face
	void deepdeal(int q, int face_current)// maintain hull, update if point q is outside
	{
		facesNumber[face_current].is_final_face = 0;// remove current face since it becomes internal

		// reverse edges to find adjacent faces
		deal(q, facesNumber[face_current].b, facesNumber[face_current].a);
		deal(q, facesNumber[face_current].c, facesNumber[face_current].b);
		deal(q, facesNumber[face_current].a, facesNumber[face_current].c);
	}



	void make()// construct 3D convex hull
	{
		/*Part1*/
		conhull_triangle_index = 0;
		if (n < 4)
		{
			cout << "Error! The number of points is less than 4, cannot construct a convex hull." << endl;
			return;// void function, exit early without return value
		}

		int indicator = 1;              // flag, bad case if indicator = 1;

		for (int i = 1; i < n; i++)   // ensure first two points are not identical
		{
			static int concurrentCounter = 1;
			if ((p[0] - p[i]).getLength() > eps1)
			{
				swap(p[1], p[i]);
				indicator = 0;
				break;
			}
			else { concurrentCounter++; }

			if (concurrentCounter == n) // all points are identical
			{
				cout << "Error! All points are the same, cannot construct a convex hull." << endl;
				return;
			}
		}
		if (indicator) return;

		indicator = 1;

		for (int i = 2; i < n; i++)   // ensure first three points are not collinear
		{
			static int collinearCounter = 2;
			if (fabs((p[1] - p[0]) * (p[i] - p[0]) - 1.0) > eps1)
			{
				swap(p[2], p[i]);
				indicator = 0;
				break;
			}
			else { collinearCounter++; }

			if (collinearCounter == n) // all points are collinear
			{
				cout << "Error! All points are collinear, cannot construct a convex hull." << endl;
				return;
			}
		}
		if (indicator) return;

		indicator = 1;

		for (int i = 3; i < n; i++)// ensure first four points are not coplanar
		{
			static int coplanarCounter = 3;
			if (fabs(volume(p[0], p[1], p[2], p[i])) > eps1)
			{
				swap(p[3], p[i]);
				indicator = 0;
				break;
			}
			else { coplanarCounter++; }

			if (coplanarCounter == n) // all points are coplanar
			{
				cout << "Error! All points are coplanar, cannot construct a convex hull." << endl;
				return;
			}
		}
		if (indicator) return;

		/*Part2*/
		face add;
		for (int i = 0; i < 4; i++)// construct initial tetrahedron
		{
			add.a = (i + 1) % 4;
			add.b = (i + 2) % 4;
			add.c = (i + 3) % 4;
			add.is_final_face = 1;

			if (pointToFace(p[i], add) > eps1)   swap(add.c, add.b);    // ensure outward normal

			face_containing_edge[add.a][add.b]
				= face_containing_edge[add.b][add.c]
				= face_containing_edge[add.c][add.a]
				= conhull_triangle_index;

			facesNumber[conhull_triangle_index++] = add;
		}

		/*Part3*/
		for (int i = 4; i < n; i++)// incremental update of hull
		{
			for (int j = 0; j < conhull_triangle_index; j++)
			{
				if (facesNumber[j].is_final_face && pointToFace(p[i], facesNumber[j]) > eps1)
				{
					deepdeal(i, j);
					break;
				}
			}
		}

		/*Part4*/
		int tmp = conhull_triangle_index; // remove non-final faces
		conhull_triangle_index = 0;

		for (int i = 0; i < tmp; i++)
			if (facesNumber[i].is_final_face)
			{
				facesNumber[conhull_triangle_index++] = facesNumber[i];

				vertex_index_vector.push_back(facesNumber[i].a);
				vertex_index_vector.push_back(facesNumber[i].b);
				vertex_index_vector.push_back(facesNumber[i].c);
				final_faces.push_back(facesNumber[i]);
			}

		final_faces_number = final_faces.size();

		ofstream outfile;
		outfile << std::scientific << std::setprecision(6) << std::showpos;
		outfile.open("convex_hull_vertices.out");

		outfile << vertex_index_vector.size() << endl;

		for (int i = 0; i < vertex_index_vector.size(); i++)
		{
			outfile << p[vertex_index_vector[i]].x << "   "
				<< p[vertex_index_vector[i]].y << "   "
				<< p[vertex_index_vector[i]].z << endl;

			convex_hull_vertices.push_back(p[vertex_index_vector[i]]);
		}
		outfile.close();

		cout << "<Convex hull vertices saved to file \"convex_hull_vertices.out\"." << endl;

		struct Vec3Hash {
			size_t operator()(const vec3& v) const {
				return std::hash<float>()(v.x) ^ std::hash<float>()(v.y) ^ std::hash<float>()(v.z);
			}
		};

		struct Vec3Equal {
			bool operator()(const vec3& a, const vec3& b) const {
				return a.x == b.x && a.y == b.y && a.z == b.z;
			}
		};

		std::vector<vec3> tmp_ch_vertices = convex_hull_vertices;

		std::unordered_set<vec3, Vec3Hash, Vec3Equal> unique_vertices(
			tmp_ch_vertices.begin(), tmp_ch_vertices.end()
		);

		tmp_ch_vertices.assign(unique_vertices.begin(), unique_vertices.end());

		final_vertices_number = tmp_ch_vertices.size(); // final number of vertices

		for (const auto& vertex : tmp_ch_vertices)
		{
			geometric_center = geometric_center + vertex;
		}
		geometric_center = geometric_center / tmp_ch_vertices.size(); // geometric center

		if (isPointInConvexHull(geometric_center))
		{
			cout << "<Geometric center is inside the convex hull." << endl;
		}
		else
		{
			cout << "Error! Geometric center is NOT inside the convex hull!" << endl;
		}

		for (int i = 0; i < vertex_index_vector.size() / 3; i++)
		{
			vec3 point_a, point_b, point_c, ab, ac;
			point_a = convex_hull_vertices[3 * i];
			point_b = convex_hull_vertices[3 * i + 1];
			point_c = convex_hull_vertices[3 * i + 2];

			ab = point_b - point_a;
			ac = point_c - point_a;

			faces_normal.push_back((ab.crossPdt(ac)).normalization()); // compute face normal

			if (faces_normal[i] * (geometric_center - point_a) > 0.0) // check if outward
			{
				cout << "Error! Face normal should be outward!" << endl;
			}
		}
	}
};

extern Convex_hull hull;
