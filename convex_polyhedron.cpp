#include "convex_polyhedron.h"




vector<vec3> removeDuplicateVertices(const std::vector<vec3>& verticesIn) {
    std::vector<vec3> result;

    for (const auto& vertex : verticesIn) {
        // Check whether the current vertex already exists in the result vector
        auto it = std::find_if(result.begin(), result.end(),
            [&vertex](const vec3& existingVertex) {
                return vertex == existingVertex;
            });

        // If no duplicate is found, add it to the result
        if (it == result.end()) {
            result.push_back(vertex);
        }
    }

    return result;
}

//ensure that the face vertices are in Counter-ClockWise (CCW) order

void ensureCCW(Face& face) {
    auto& verts = face.getVertices();
    if (verts.size() < 3) return;

    // 1. Compute the face center (centroid)
    vec3 center(0.0f, 0.0f, 0.0f);
    for (const auto& v : verts) {
        center = center + v;
    }
    center = center / static_cast<float>(verts.size());

    // 2. Get the face normal and normalize it
    vec3 normal = face.getNormal();
    normal.normalization();

    // 3. Choose a reference direction (from center to the first vertex)
    vec3 refDir = verts[0] - center;
    refDir.normalization();

    // 4. Compute the polar angle of each vertex relative to refDir (on the plane perpendicular to the normal)
    std::vector<std::pair<float, size_t>> angleIndices;
    for (size_t i = 0; i < verts.size(); ++i) {
        vec3 dir = verts[i] - center;
        dir.normalization();

        // Compute the polar angle: use atan2 to calculate the angle on the plane perpendicular to the normal
        // Use cross product and dot product to determine the angle
        float dot = refDir * dir;
        vec3 cross = refDir.crossPdt(dir);
        float angle = atan2(cross * normal, dot);
        angleIndices.emplace_back(angle, i);
    }

    // 5. Sort by polar angle (from small to large, i.e., counterclockwise)
    std::sort(angleIndices.begin(), angleIndices.end(),
        [](const auto& a, const auto& b) {
            return a.first < b.first;
        });

    // 6. Rearrange vertices
    std::vector<vec3> sortedVerts;
    for (const auto& pair : angleIndices) {
        sortedVerts.push_back(verts[pair.second]);
    }

    // 7. Update the face vertex order
    face.setVertices(sortedVerts);

    // Optional: verify whether the normal points outward (may require additional computation)
}

Convex_polyhedron MergeIfTriangularFacesOnSamePlane(vector<vec3>& verticesFromConHull)
{
    Convex_polyhedron ConvexPldron(verticesFromConHull);
    for (size_t i = 0; i < ConvexPldron.CP_vertices.size() / 3; ++i) {
        vector<vec3> tmp_vertices;

        tmp_vertices.push_back(ConvexPldron.CP_vertices[i * 3]);
        tmp_vertices.push_back(ConvexPldron.CP_vertices[i * 3 + 1]);
        tmp_vertices.push_back(ConvexPldron.CP_vertices[i * 3 + 2]);
        ConvexPldron.CP_faces.push_back(Face(tmp_vertices));
        ConvexPldron.CP_faces[i].calculateNormal(); // Compute face normal for each face

    }

    for (size_t i = 0; i < ConvexPldron.CP_faces.size(); ) {
        bool erased = false;
        for (size_t j = i + 1; j < ConvexPldron.CP_faces.size(); ) {
            if (ConvexPldron.CP_faces[i].getNormal().almostEqual(ConvexPldron.CP_faces[j].getNormal())) {
                ConvexPldron.CP_faces[i].getVertices().insert(
                    ConvexPldron.CP_faces[i].getVertices().end(),
                    ConvexPldron.CP_faces[j].getVertices().begin(),
                    ConvexPldron.CP_faces[j].getVertices().end()
                );
                ConvexPldron.CP_faces.erase(ConvexPldron.CP_faces.begin() + j);
                erased = true;
            }
            else {
                ++j;
            }
        }


        if (!erased) {
            ++i;
        }
    }
    // Remove duplicate vertices from each face
    ConvexPldron.CP_vertices.clear();
    for (auto& face : ConvexPldron.CP_faces) {
        face.getVertices() = removeDuplicateVertices(face.getVertices());
        ensureCCW(face); // Ensure the vertices are in CCW order

        ConvexPldron.CP_vertices.insert(ConvexPldron.CP_vertices.end(), face.getVertices().begin(), face.getVertices().end());

    }

    //recheck if the face normals are outward facing
    for (auto& face : ConvexPldron.CP_faces) {
        vec3 faceNormal = face.getNormal();
        if (faceNormal * (geometric_center - face.getVertices()[0]) > 0.0) // recheck whether the normal points outward
        {
            cout << "Error warning! Convex Polyhedron's face normal should be outward!" << endl;
        }
    }


    return ConvexPldron;
}
