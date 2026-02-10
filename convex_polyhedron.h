#pragma once

#include "prism.h" //for Face class 
#include "convexhull.h"

//warning: Pay attention to distinguishing between "class Face" and "struct face"!!!


class Convex_polyhedron
{
public:
    vector<Face> CP_faces;  // CP stands for Convex_Polyhedron
	vector<vec3> CP_vertices; // Vertices of the convex polyhedron
  //  size_t CP_faces_number;                 // Number of faces

public:
    // Default constructor
    Convex_polyhedron() {
        // Initialize empty vectors for faces and vertices
        CP_faces = {};
        CP_vertices = {};
    }

    // Constructor - Create prism from constructed convex hull faces
    Convex_polyhedron(const std::vector<Face>& convexHullFaces) : CP_faces(convexHullFaces) {
        // Validate prism geometry
        if (CP_faces.size() < 4) {
            throw std::invalid_argument("Error warning: A convex polyhedron must have at least 4 faces.");
        }

    }
    // Constructor - Create prism from constructed convex hull faces
    Convex_polyhedron(const std::vector<vec3>& convexHullVertices) : CP_vertices(convexHullVertices) {
        // Validate prism geometry
        if (CP_vertices.size() < 4) {
            throw std::invalid_argument("Error warning: A convex polyhedron must have at least 4 vertices.");
        }

    }

};

vector<vec3> removeDuplicateVertices(const std::vector<vec3>& verticesIn);
Convex_polyhedron MergeIfTriangularFacesOnSamePlane(vector<vec3>& verticesFromConHull);
void ensureCCW(Face& face);