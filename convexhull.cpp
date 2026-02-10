#include "convexhull.h"


vec3 p[M0];// All initially given points
face facesNumber[M0 * 8];// Convex hull triangles


vector<int> vertex_index_vector;
vector<vec3> convex_hull_vertices; // Convex hull vertex coordinates
vector<vec3> faces_normal; // Normals of convex hull faces
vector<face> final_faces; // Faces of the completed convex hull
Convex_hull hull;
vec3 geometric_center;


double semi_axis_a, semi_axis_b, semi_axis_c;
//note: It's radius, NOT diameter!
double ConHull_max_radius; //for main.cpp, the maximum radius of the convex hull. for sample ray's initial points.



void readIn(vec3 point[]) {
    const string& filename = "convex_hull_vertices.in";
    cout << "<Reading parameters from \" " << filename << " \" " << endl;
    ifstream ReadIn;
    ReadIn.open(filename);

    if (!ReadIn.is_open()) {
        cerr << "Error: Could not open file " << filename << ", Please check input data file name." << endl;
        return;
    }
    ReadIn >> hull.n;
    for (int i = 0; i < hull.n; i++)
    {
        ReadIn >> point[i].x >> point[i].y >> point[i].z;
    }
    ReadIn.close();
}


void outputReadInDataForCheck(vec3 point[]) {
    const string& filename = "convex_hull_vertices_ReadInCheck.out";
    cout << "<Printing the number of points loaded and their coordinates to file \" " << filename << " \" " << endl;
    ofstream outfile03;
    outfile03 << std::scientific << std::setprecision(6) << std::showpos;
    outfile03.open(filename);

    if (!outfile03.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        return;
    }

    outfile03 << hull.n << endl;
    for (int i = 0; i < hull.n; i++) {
        outfile03 << point[i].x << "   "
            << point[i].y << "   "
            << point[i].z << endl;
    }
    outfile03.close();
    // outfile03.flush();
}





bool isPointInConvexHull(vec3& p)
{

    int i;
    for (i = 0; i < final_faces.size(); i++)
    {
        if (hull.pointToFace(p, final_faces[i]) > eps1) {

            cout << "<The point is outside the convex hull." << endl;
            return false;
        }
        else {
            cout << "<The point is inside the convex hull." << endl;
            return true;
        }
    }


}


void readInPointCoords(vec3 p[]) {
    cout << "<Reading Point number and coordinates from file\"pointCoords.in\" " << endl;
    ifstream ReadIn;
    ReadIn.open("pointCoords.in");

    if (!ReadIn) {
        cerr << "Error: Could not open \"pointCoords.in\"" << endl;
        return;
    }

    ReadIn >> hull.n;
    for (int i = 0; i < hull.n; i++) {
        ReadIn >> p[i].x >> p[i].y >> p[i].z;
    }

    ReadIn.close();
}

void generateRandomPoints(vec3 p[])
{
    //hull.n = ceil(4.0 + rn() * 190.0);//  input points number setup
    hull.n = 25;

    semi_axis_a = 1.0;
    semi_axis_b = 1.0;
    semi_axis_c = 1.0;

    for (int i = 0; i < hull.n; i++)
    {
        p[i].x = -semi_axis_a + 2.0 * semi_axis_a * rn();
        p[i].y = -semi_axis_b + 2.0 * semi_axis_b * rn();
        p[i].z = -semi_axis_c + 2.0 * semi_axis_c * rn();

    }
}

void generatePointsOnSphere(vec3 p[])
{
    double sphere_radius = 1.0; // Sphere radius
    //hull.n = ceil(4.0 + rn() * 190.0);//  input points number setup
    hull.n = 100;
    double sita, phi;
    for (int i = 0; i < hull.n; i++)
    {
        //sita = pi * rn();
        sita = acos(1.0 - 2.0 * rn()); //sita ~ [0, pi]
        phi = 2.0 * pi * rn();

        p[i].x = sphere_radius * sin(sita) * cos(phi);
        p[i].y = sphere_radius * sin(sita) * sin(phi);
        p[i].z = sphere_radius * cos(sita);

    }
}

void generatePointsOnSphere02(vec3 p[])
{
    double sphere_radius = 1.0; // Sphere radius
    int sita_N = 10, phi_N = 10, hull_n_counter;

    double sita, phi;
    hull.n = sita_N * phi_N + 2;
    //cout << "hull.n = " << hull.n << endl;
    //two poles
    p[0].x = 0.0;
    p[0].y = 0.0;
    p[0].z = sphere_radius;

    p[1].x = 0.0;
    p[1].y = 0.0;
    p[1].z = -sphere_radius;
    hull_n_counter = 2; //two poles already assigned

    for (int sita_i = 0; sita_i < sita_N; sita_i++)
    {
        for (int phi_j = 0; phi_j < phi_N; phi_j++)
        {
            //sita = pi * (sita_i + 0.5) / sita_N; //sita ~ [0, pi]  
            sita = acos(1.0 - 2.0 * (sita_i + 0.5) / sita_N); // Better ¦È distribution
            phi = 2.0 * pi * phi_j / phi_N; //phi ~ [0, 2pi]

            p[hull_n_counter].x = sphere_radius * sin(sita) * cos(phi);
            p[hull_n_counter].y = sphere_radius * sin(sita) * sin(phi);
            p[hull_n_counter].z = sphere_radius * cos(sita);
            hull_n_counter++;
        }
    }

}

//a, b, c are the semi-axis lengths of the ellipsoid
void generatePointsOnEllipsoid(vec3 p[])
{
    semi_axis_a = 0.2, semi_axis_b = 0.5, semi_axis_c = 1.0; // Semi-axis lengths of the ellipsoid
    int sita_N = 20, phi_N = 10, hull_n_counter;
    double sita, phi;

    // Total number of points = two poles + sita_N * phi_N
    hull.n = sita_N * phi_N + 2;


    // Two poles (top and bottom)
    p[0].x = 0.0;
    p[0].y = 0.0;
    p[0].z = semi_axis_c;  // Top pole, z-axis = c

    p[1].x = 0.0;
    p[1].y = 0.0;
    p[1].z = -semi_axis_c; // Bottom pole, z-axis = -c
    hull_n_counter = 2; // Two poles already assigned

    for (int sita_i = 0; sita_i < sita_N; sita_i++)
    {
        for (int phi_j = 0; phi_j < phi_N; phi_j++)
        {
            // Improved ¦È distribution to avoid dense points near poles
            sita = acos(1.0 - 2.0 * (sita_i + 0.5) / sita_N); // ¦È ¡Ê [0, ¦Ð]
            phi = 2.0 * pi * phi_j / phi_N;                      // ¦Õ ¡Ê [0, 2¦Ð]

            // Ellipsoid parametric equation: x = a*sin¦È*cos¦Õ, y = b*sin¦È*sin¦Õ, z = c*cos¦È
            p[hull_n_counter].x = semi_axis_a * sin(sita) * cos(phi);
            p[hull_n_counter].y = semi_axis_b * sin(sita) * sin(phi);
            p[hull_n_counter].z = semi_axis_c * cos(sita);

            hull_n_counter++;
        }
    }
}

//a, b, c are the semi-axis lengths of the ellipsoid
void generatePointsOnEllipsoidWithOnePole(vec3 p[])
{
    semi_axis_a = 0.5, semi_axis_b = 0.5, semi_axis_c = 1.0; // Semi-axis lengths of the ellipsoid
    int sita_N = 2, phi_N = 6, hull_n_counter;
    double sita, phi;


    hull.n = sita_N * phi_N + 1;



    p[0].x = 0.0;
    p[0].y = 0.0;
    p[0].z = semi_axis_c;  // Top pole, z-axis = c


    hull_n_counter = 1; // One pole already assigned

    for (int sita_i = 0; sita_i < sita_N; sita_i++)
    {
        for (int phi_j = 0; phi_j < phi_N; phi_j++)
        {
            // Improved ¦È distribution to avoid dense points near poles
            sita = acos(1.0 - 2.0 * (sita_i + 0.5) / sita_N); // ¦È ¡Ê [0, ¦Ð]
            phi = 2.0 * pi * phi_j / phi_N;                      // ¦Õ ¡Ê [0, 2¦Ð]

            // Ellipsoid parametric equation: x = a*sin¦È*cos¦Õ, y = b*sin¦È*sin¦Õ, z = c*cos¦È
            p[hull_n_counter].x = semi_axis_a * sin(sita) * cos(phi);
            p[hull_n_counter].y = semi_axis_b * sin(sita) * sin(phi);
            p[hull_n_counter].z = semi_axis_c * cos(sita);

            hull_n_counter++;
        }
    }
}

//a, b, c are the semi-axis lengths of the ellipsoid
void generatePointsOnEllipsoidWithoutPoles(vec3 p[])
{
    semi_axis_a = 1, semi_axis_b = 3, semi_axis_c = 5; // Semi-axis lengths of the ellipsoid
    int sita_N = 20, phi_N = 10, hull_n_counter = 0;
    double sita, phi;


    hull.n = sita_N * phi_N;


    for (int sita_i = 0; sita_i < sita_N; sita_i++)
    {
        for (int phi_j = 0; phi_j < phi_N; phi_j++)
        {
            // Improved ¦È distribution to avoid dense points near poles
           // sita = acos(1.0 - 2.0 * (sita_i + 0.5) / sita_N); // ¦È ¡Ê [0, ¦Ð]
            sita = pi * (sita_i + 0.5) / sita_N; //sita ~ [0, pi]  
            phi = 2.0 * pi * phi_j / phi_N;                      // ¦Õ ¡Ê [0, 2¦Ð]

            // Ellipsoid parametric equation: x = a*sin¦È*cos¦Õ, y = b*sin¦È*sin¦Õ, z = c*cos¦È
            p[hull_n_counter].x = semi_axis_a * sin(sita) * cos(phi);
            p[hull_n_counter].y = semi_axis_b * sin(sita) * sin(phi);
            p[hull_n_counter].z = semi_axis_c * cos(sita);

            hull_n_counter++;
        }
    }
}

//a, b are the semi-axis lengths of the base ellipse, h is the semi-height of the cylinder 
void generatePointsOnEllipticCylinder(vec3 p[])
{
    semi_axis_a = 2.0, semi_axis_b = 2.0, semi_axis_c = 5.0; // Semi-axis lengths of the elliptical cylinder, h is the semi-height
    ConHull_max_radius = sqrt(semi_axis_a * semi_axis_a + semi_axis_b * semi_axis_b + semi_axis_c * semi_axis_c);  // Maximum radius of the elliptical cylinder

    int height_N = 0, phi_N = 6, hull_n_counter = 0;
    double height, phi;


    hull.n = (height_N + 2) * phi_N;


    for (int height_i = 0; height_i < height_N + 2; height_i++)
    {
        for (int phi_j = 0; phi_j < phi_N; phi_j++)
        {
            height = -semi_axis_c + 2.0 * semi_axis_c * height_i / (height_N + 1.0); // height ~ [-h, h]
            phi = 2.0 * pi * phi_j / phi_N;            // ¦Õ ¡Ê [0, 2¦Ð]

            p[hull_n_counter].x = semi_axis_a * cos(phi);
            p[hull_n_counter].y = semi_axis_b * sin(phi);
            p[hull_n_counter].z = height;

            hull_n_counter++;
        }
    }
}
