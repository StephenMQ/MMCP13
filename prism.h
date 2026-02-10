#pragma once
#include <vector>
#include <algorithm> 
#include <cmath>

#include"vec3.h"
#include "namesp.h"
#include "geometry.h"

//#include <stdexcept>
using namespace std;
using namespace consts;
const double kEpsilon = 1e-6; ///< Tolerance for floating-point comparisons

class Face {
private:
    vector<vec3> vertices;  // Vertices of the face
    vec3 normal;                 // Normal vector of the face
    bool normalCalculated;       // Flag indicating if normal has been calculated

public:
    // Constructor - Initialize face with vertices
    Face(const vector<vec3>& verts) : vertices(verts), normalCalculated(false) {
        if (verts.size() < 3) {
            throw invalid_argument("A face must have at least 3 vertices");
        }
    }

    // Get number of vertices
    size_t getVertexCount() const {
        return vertices.size();
    }

    vector<vec3>&  setVertices(const vector<vec3>& verts) {
        vertices = verts;
        normalCalculated = false; // Reset normal calculation flag
        return vertices;
	}

    // Get list of vertices
    const std::vector<vec3>& getVertices() const {
        return vertices;
    }

	//non const version
    vector<vec3>& getVertices() { return vertices; }

    // Get normal vector (calculates if not already done)
    vec3 getNormal() {
        if (!normalCalculated) {
            calculateNormal();
            normalCalculated = true;
        }
        return normal;
    }

    //const case
    vec3 getNormal() const {       
        return normal;
    }

    // Set normal vector (manual specification)
    void setNormal(const vec3& n) {
        normal = n.normalization();
        normalCalculated = true;
    }

	bool getNormalCalculated() const {
		return normalCalculated;
	}
	bool setNormalCalculated(bool value) {
		normalCalculated = value;
		return normalCalculated;
	}

public:
    // Calculate face normal using cross product of edges
    void calculateNormal() {
        // Use first three vertices to compute normal
        vec3 v1 = vertices[1] - vertices[0];
        vec3 v2 = vertices[2] - vertices[0];
        normal = v1.crossPdt(v2).normalization();
    }
};

/**
 * @struct IntersectionInfo
 * @brief Stores information about line-prism intersection
 */
struct IntersectionInfo 
{
    vec3 point;             ///< Intersection point coordinates
    double distance;        ///< Distance from line origin to intersection
    const Face* face = nullptr; ///< Pointer to intersected face (null if edge/vertex)
    bool isValid = false;   ///< Flag indicating valid intersection

    /// Comparison operator for sorting
    bool operator<(const IntersectionInfo& other) const {
        return distance < other.distance;
    }
};

/// <summary>
/// The center of the starting prism must coincide with the coordinate origin and be placed vertically on the xy-plane.
/// </summary>
class NPrism {
private:
    vector<Face> faces;  // All faces of the prism
    vector<vec3> vertices;  // All vertices of the prism
    size_t n;                 // Number of sides of the prism

public:
    // Constructor01 - Create regular n-sided prism
    NPrism(size_t numberSides, double diameter, double height) : n(numberSides) {
        if (numberSides < 3) {
            throw std::invalid_argument("A prism must have at least 3 sides");
        }

        // Create vertices for bottom and top faces
       vector<vec3> bottomVertices;
       vector<vec3> topVertices;

        // Calculate vertices in circular pattern
        for (size_t i = 0; i < numberSides; ++i) {
            double angle = 2 * pi * i / numberSides;
            double radius = diameter / 2.0;
            double x = radius * cos(angle);
            double y = radius * sin(angle);

            bottomVertices.emplace_back(x, y, -height / 2.0);
            topVertices.emplace_back(x, y, height/2.0);
        }

        // Create bottom and top faces
        Face bottomFace(bottomVertices);       
        bottomFace.setNormal(vec3(0.0, 0.0, -1.0));  // Bottom face normal points downward

        Face topFace(topVertices);       
        topFace.setNormal(vec3(0.0, 0.0, 1.0));     // Top face normal points upward

        faces.push_back(bottomFace);
        faces.push_back(topFace);

        // Create rectangular side faces
        for (size_t i = 0; i < numberSides; ++i) 
        {
            size_t next = (i + 1) % numberSides;//careful!!!
            std::vector<vec3> sideVertices = 
            {
                bottomVertices[i],
                bottomVertices[next],
                topVertices[next],
                topVertices[i]
            };
            Face sideFace(sideVertices);          

            //Calculate normals in circular pattern
            double angle = 2 * pi * i / numberSides + 2 * pi / numberSides / 2;
            double radius = diameter / 2.0;
            double x = cos(angle);
            double y = sin(angle);
            sideFace.setNormal(vec3(x, y, 0));//sides face normal points outward

            faces.push_back(sideFace);

        }

    }

    // Constructor02 - Create prism from custom face collection
    NPrism(const std::vector<Face>& customFaces) : faces(customFaces) {
        // Validate prism geometry
        if (faces.size() < 4) {  // Minimum 4 faces
            throw std::invalid_argument("A convex polyhedron must have at least 4 faces.");
        }

    }

    // Constructor03 - Create prism from vertices
    NPrism(const std::vector<vec3>& vertices_Vector) : vertices(vertices_Vector) {
        // Validate prism geometry
        if (vertices.size() < 4) {  // Minimum 4 vertices
            throw std::invalid_argument("A convex polyhedron must have at least 4 vertices.");
        }

    }

    void rotateCrystal(double gamma, double beta, double alpha)
    {
        for (auto& face : faces)
        {
			for (auto& vertex : face.getVertices())
			{
				// Rotate each vertex using Euler angles
				//vertex = rotateVec3(vertex, gamma, beta, alpha);
               // vertex = rotateVec3_02(vertex, gamma, beta, alpha);
                vertex = rotateVec3_ZYZ02(vertex, gamma, beta, alpha);
			}
            face.calculateNormal();
            face.setNormalCalculated(true);
        }

    }

    void rotateCrystal_BA(double beta, double alpha)
    {
        for (auto& face : faces)
        {
            for (auto& vertex : face.getVertices())
            {
                // Rotate each vertex using Euler angles
                //vertex = rotateVec3(vertex, gamma, beta, alpha);
               // vertex = rotateVec3_02(vertex, gamma, beta, alpha);
               // vertex = rotateVec3_ZYZ(vertex, gamma, beta, alpha);
                vertex = rotateVec3_YZ(vertex, beta, alpha);
            }
            face.calculateNormal();
            face.setNormalCalculated(true);
        }

    }

    // Get number of sides
    size_t getSideCount() const {
        return n;
    }

    // Get total number of faces
    size_t getFaceCount() const {
        return faces.size();
    }

    // Get all faces
    const std::vector<Face>& getFaces() const {
        return faces;
    }

    // Get specific face by index
    const Face& getFace(size_t index) const {
        if (index >= faces.size()) {
            throw std::out_of_range("Face index out of range");
        }
        return faces[index];
    }

    /**
 * @brief Pre-filters potentially intersecting faces using normal dot product
 * @param lineDir normalizationd direction vector of the line
 * @return Vector of face pointers that may intersect with the line
 *
 * @details This implements fast preliminary culling by checking the angle between
 * the line direction and each face normal. Only faces where the line enters from
 * the front (dot product < -ε) are considered for actual intersection tests.
 */
    std::vector<const Face*> prefilterFaces(const line& ray) const {
        std::vector<const Face*> candidateFaces;

        for (const auto& face : faces) 
        {
            vec3 faceNormal = face.getNormal();
            double dot = ray.direction_d * faceNormal;  // Dot product = |D||N|cosθ

            // Mathematical insight:
            // - dot ≈ 0: Line parallel to face (θ ≈ 90°)
            // - dot > 0: Line faces back side (θ < 90°)
            // - dot < 0: Line faces front side (θ > 90°)
            if (dot < 0.0) {
                candidateFaces.push_back(&face);
            }
        }

        return candidateFaces;
    }

    /**
     * @brief Finds all intersections between line and prism with optimization
     * @param linePoint Origin point of the line
     * @param lineDirection Direction vector of the line (will be normalizationd)
     * @return Sorted vector of intersection info (by distance)
     *
     * @note This optimized version:
     * 1. normalizations direction vector
     * 2. Pre-filters faces using normal culling
     * 3. Only performs precise intersection tests on candidate faces
     * 4. Returns sorted intersections
     */
   vector<IntersectionInfo> findLineIntersections( const line& ray) const
    {
       vector<IntersectionInfo> intersections;       
        // Optimization Stage 1: Face Normal Culling
        auto candidateFaces = prefilterFaces(ray);
        // Optimization Stage 2: Precise Intersection Tests
        for (const auto& face : candidateFaces) {
            checkFaceIntersection(*face, ray, intersections);
        }

        // Sort results by distance from line origin
       sort(intersections.begin(), intersections.end());
        return intersections;
    }

   void findDistance2004(const Face& face, const line& ray, vector<IntersectionInfo>& results) const
   {
       Plane3D facePlane(face.getVertices()[0], face.getNormal());
       vec3 intersectionPoint;

       if (intersectLinePlane(ray, facePlane, intersectionPoint))
       {
           
               IntersectionInfo info;
               info.point = intersectionPoint;
               double tmp = (intersectionPoint - ray.point_p) * ray.direction_d;
               if (tmp > 0.0)//note: important!!!
                   info.distance = (intersectionPoint - ray.point_p).getLength();
               else info.distance = -((intersectionPoint - ray.point_p).getLength());
               info.face = &face;
               info.isValid = true;
               results.push_back(info);
           
       }

   }

   /// <summary>
   /// ***2004 - for method used from paper Zhang2004.
   /// </summary>
   /// <param name="ray"></param>
   /// <returns></returns>
  
   vector<IntersectionInfo> findLineIntersections2004(const line& ray) const
   {
       vector<IntersectionInfo> possibleIntersections;
       vector<IntersectionInfo> impossibleIntersections;
       vector<IntersectionInfo> tmp;

       vector<const Face*> possibleFaces;
       vector<const Face*> impossibleFaces;

       for (const auto& face : faces)
       {
           vec3 faceNormal = face.getNormal();
           double dot = ray.direction_d * faceNormal; // Dot product = |D||N|cosθ
           if (dot < 0.0) 
           {
               possibleFaces.push_back(&face);
           }
           else impossibleFaces.push_back(&face);
       }

       for (const auto& face : possibleFaces) {
           findDistance2004(*face, ray, possibleIntersections);
       }
       for (const auto& face : impossibleFaces) {
           findDistance2004(*face, ray, impossibleIntersections);
       }

       // Sort results by distance from line origin
       sort(possibleIntersections.begin(), possibleIntersections.end());
       sort(impossibleIntersections.begin(), impossibleIntersections.end());
       if (possibleIntersections.size() < 1 || impossibleIntersections.size() < 1)
       {
           cout << "error found in function findLineIntersections2004" << endl;
           
       }
       if (possibleIntersections[possibleIntersections.size() - 1].distance < impossibleIntersections[0].distance) //max(d_i)<min(D_i) intesects!
           return possibleIntersections;
       else
           return tmp;
   }

   //following Zhang2004 method note: must careful with the start point face!!! need to be added!!!
   vector<IntersectionInfo> findNextLineIntersections2004(const line& ray, const Face* faceToSkip = nullptr) const
   {
       vector<IntersectionInfo> intersections;
      // vector<const Face*> possibleFaces;
      
       for (const auto& face : faces)
       {
           vec3 faceNormal = face.getNormal();
           double dot = ray.direction_d * faceNormal; // Dot product = |D||N|cosθ

           if (faceToSkip != &face && dot > 0.0) findNextFaceIntersection(face, ray, intersections);
       }        
      
       sort(intersections.begin(), intersections.end());

       return intersections;
   }



    /**
    * Checks if a line intersects with the prism by testing against each face
    * @param linePoint A point on the line
    * @param lineDirection Direction vector of the line
    * @return True if the line intersects the prism, false otherwise
    */
    bool doesLineIntersect(const line& ray) const {
        // 1. Check intersection with bottom and top faces
        const Face& bottomFace = faces[0];
        const Face& topFace = faces[1];

        Plane3D bottomPlane(bottomFace.getVertices()[0], bottomFace.getNormal());
        Plane3D topPlane(topFace.getVertices()[0], topFace.getNormal());

        // Store potential intersection points
        std::vector<vec3> intersectionPoints;

        vec3 bottomIntersection;
        if (intersectLinePlane(ray, bottomPlane, bottomIntersection)) {
            if (isPointInFace(bottomFace, bottomIntersection)) {
                return true; // Intersects bottom face
            }
        }

        vec3 topIntersection;
        if (intersectLinePlane(ray, topPlane, topIntersection)) {
            if (isPointInFace(topFace, topIntersection)) {
                return true; // Intersects top face
            }
        }

        // 2. Check intersection with side faces
       
        vec3 sideIntersection;
        for (size_t i = 2; i < faces.size(); ++i) 
        {
            const Face& sideFace = faces[i];
            Plane3D sidePlane(sideFace.getVertices()[0], sideFace.getNormal());
            const auto& verts = sideFace.getVertices();           
            if (intersectLinePlane(ray, sidePlane, sideIntersection))
            {
                    if (isPointInFace(sideFace, sideIntersection))
                    {
                        return true; // Intersects side face
                    }
            }
           
        }

        return false; // No intersections found
    }

   vector<IntersectionInfo> findFirstLineIntersections(const line& ray) const
    {
       vector<IntersectionInfo> intersections;
       //0. Check initial point of line "on or not on" prism faces
       
        // 1. Check intersections with bottom and top faces
         checkFaceIntersection(faces[0], ray, intersections);
         checkFaceIntersection(faces[1], ray, intersections);

        // 2. Check intersections with side faces
        for (size_t i = 2; i < faces.size(); ++i) 
        {
            checkFaceIntersection(faces[i], ray, intersections);
        }

        // 3. Sort intersections by distance
        std::sort(intersections.begin(), intersections.end());

        if (intersections.size() > 2)
        {
            cout << "intersections.size()="<< intersections.size() <<endl;
            for (size_t i = 0; i < intersections.size(); ++i)
            {
               cout << "intersections[i].distance" << intersections[i].distance << endl;
               intersections[i].point.show();
            }
            cout << "Error in function .findLineIntersections, intersections.size should be 1 or 2." << endl;
        }
        
        return intersections;
    }

   vector<IntersectionInfo> findNextLineIntersections(const line& ray, const Face* faceToSkip = nullptr) const
   {
       vector<IntersectionInfo> intersections;
       //0. Check initial point of line "on or not on" prism faces

        // 1. Check intersections with bottom and top faces
       if (faceToSkip != &faces[0])  findNextFaceIntersection(faces[0], ray, intersections);
       if (faceToSkip != &faces[1])  findNextFaceIntersection(faces[1], ray, intersections);

       // 2. Check intersections with side faces
       for (size_t i = 2; i < faces.size(); ++i)
       {
           if (faceToSkip != &faces[i]) findNextFaceIntersection(faces[i], ray, intersections);
       }

       // 3. Sort intersections by distance
       std::sort(intersections.begin(), intersections.end());

       return intersections;
   }


    /**
    * Checks intersection between line and plane
    * @param linePoint Point on line
    * @param lineDir Direction of line
    * @param plane The plane to test against
    * @param intersection Output intersection point if exists
    * @return True if intersection exists
    */
    bool intersectLinePlane(const line& ray, const Plane3D& plane, vec3& intersection) const 
    {
        double denom = plane.getNormal() * ray.direction_d;

        // Line is parallel to plane
        if (fabs(denom) < 1e-6) {
            return false;
        }

        // Calculate intersection parameter
        double t = (plane.getNormal() * (plane.getPoint() - ray.point_p)) / denom;
        intersection = ray.point_p + ray.direction_d * t;
        return true;
    }

    /**
     * Checks if point is inside a face polygon
     * @param face The face to test against
     * @param point The point to check
     * @return True if point is inside face
     */
    bool isPointInFace(const Face& face, const vec3& point) const 
    {
        const auto& verts = face.getVertices();
        size_t n = verts.size();
        double totalAngle = 0.0;

        // Sum angles between point and each edge
        for (size_t i = 0; i < n; ++i) {
            vec3 v1 = verts[i] - point;
            vec3 v2 = verts[(i + 1) % n] - point;
            if(v1.getLength() <1e-6|| v2.getLength()<1e-6)  return false;//special case 1: The point coincides with the vertex.
            double tmp = (v1 * v2) / (v1.getLength() * v2.getLength());
            if ((fabs(tmp) - 1.0) < eps)
            {
                if (tmp > 1.0) tmp = 1.0;
                if (tmp < -1.0) tmp = -1.0;
            }
            double angle = acos(tmp);
            //if (isnan(angle) || isinf(angle)) return false;//special case 1: The point coincides with the vertex.
            if (fabs(angle - pi) < 1e-3) return false;//special case 2: The point coincides with the edge.
            totalAngle += angle;
        }
       
        // Point is inside if total angle is ~2π
        return fabs(totalAngle - 2 * pi) < 1e-3;
    }

    

    /**
   * @brief Checks intersection with a single face and stores results
   * @param face Face to test against
   * @param origin Line origin point
   * @param dir Line direction vector
   * @param[out] results Collection to store found intersections
   */
    void checkFaceIntersection(const Face& face,  const line& ray, vector<IntersectionInfo>& results) const
    {
        Plane3D facePlane(face.getVertices()[0], face.getNormal());
        vec3 intersectionPoint;

        if (intersectLinePlane(ray, facePlane, intersectionPoint)) 
        {
            if (isPointInFace(face, intersectionPoint)) 
            {
                IntersectionInfo info;
                info.point = intersectionPoint;
                double tmp = (intersectionPoint - ray.point_p) * ray.direction_d;
                if (tmp>0.0) 
                info.distance = (intersectionPoint - ray.point_p).getLength();
                else info.distance = -(intersectionPoint - ray.point_p).getLength();
                info.face = &face;
                info.isValid = true;

                if(info.distance>0.0)   results.push_back(info);
            }
        }

    }

    void findNextFaceIntersection(const Face& face, const line& ray, vector<IntersectionInfo>& results) const
    {
        Plane3D facePlane(face.getVertices()[0], face.getNormal());
        vec3 intersectionPoint;
        if (intersectLinePlane(ray, facePlane, intersectionPoint))
        {            
                IntersectionInfo info;
                info.point = intersectionPoint;
                double tmp = (intersectionPoint - ray.point_p) * ray.direction_d;
                info.distance = (intersectionPoint - ray.point_p).getLength();

                if (tmp > 0.0&& info.distance>0.0)//note: very important!!! 
                {                   
                    info.face = &face;
                    info.isValid = true;
                    results.push_back(info);
                }             

            
        }

    }

    /**
 * @brief Prints detailed information about all faces of the prism
 *
 * Output includes:
 * - Total number of faces
 * - For each face:
 *   - Number of vertices
 *   - Coordinates of each vertex
 *   - Face normal vector
 */
    void printPrismFacesInfo() const {
        // Print total face count
        std::cout << "Prism has " << faces.size() << " faces:\n";

        // Print information for each face
        for (size_t i = 0; i < faces.size(); ++i) {
            const Face& face = faces[i];
            const auto& vertices = face.getVertices();

            // Face header
            std::cout << "\nFace " << i + 1 << ":\n";
            std::cout << "  Vertex count: " << vertices.size() << "\n";

            // Print all vertex coordinates
            std::cout << "  Vertex coordinates:\n";
            for (size_t j = 0; j < vertices.size(); ++j) {
                std::cout << "    Vertex " << j + 1 << ": ("
                    << vertices[j].x << ", "
                    << vertices[j].y << ", "
                    << vertices[j].z << ")\n";
            }

            // Print face normal
            vec3 normal = face.getNormal();
            std::cout << "  Face normal: ("
                << normal.x << ", "
                << normal.y << ", "
                << normal.z << ")\n";
        }
    }

};
