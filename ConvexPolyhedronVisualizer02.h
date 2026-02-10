/*
*This version is modified to support arbitrary convex polyhedron faces, not limited to triangles.
*20250613
*/
#pragma once
#include <gl/glut.h>
#include <vector>
#include <algorithm>
#include "vec3.h"
#include "namesp.h"
#include "prism.h"

class ConvexPolyhedronVisualizer {
private:
    float angleX, angleY;
    float scale;
    bool mouseLeftDown, mouseRightDown;
    float mouseX, mouseY;

    std::vector<Face> faces; // Support for arbitrary polygonal faces

    bool showCube;
    float transparency;
    float colorR, colorG, colorB;
    bool showColorControls;
    bool showEdges;

    int colorMapWidth, colorMapHeight;
    int colorMapStartX, colorMapStartY;
    bool pickingColor;
    int pickerX, pickerY;

public:
    ConvexPolyhedronVisualizer()
        : angleX(0), angleY(0), scale(1.0f),
        mouseLeftDown(false), mouseRightDown(false),
        showCube(false), transparency(1.0f),
        colorR(0.85f), colorG(0.92f), colorB(1.0f),
        showColorControls(true), showEdges(true),
        colorMapWidth(180), colorMapHeight(100),
        colorMapStartX(30), colorMapStartY(40),
        pickingColor(false), pickerX(-1), pickerY(-1)
    {
    }

    void initialize(const std::vector<Face>& convex_poly_faces) {
        faces = convex_poly_faces;
    }

    void pickColorFromMap(int x, int y) {
        int localX = std::min(std::max(0, x - colorMapStartX), colorMapWidth - 1);
        int localY = std::min(std::max(0, y - colorMapStartY), colorMapHeight - 1);

        pickerX = localX;
        pickerY = localY;

        float hue = localX / float(colorMapWidth);
        float sat = 1.0f - (localY / float(colorMapHeight));
        float val = 1.0f;

        float r, g, b;
        int i = int(hue * 6);
        float f = hue * 6 - i;
        float p = val * (1 - sat);
        float q = val * (1 - f * sat);
        float t = val * (1 - (1 - f) * sat);
        switch (i % 6) {
        case 0: r = val, g = t, b = p; break;
        case 1: r = q, g = val, b = p; break;
        case 2: r = p, g = val, b = t; break;
        case 3: r = p, g = q, b = val; break;
        case 4: r = t, g = p, b = val; break;
        case 5: r = val, g = p, b = q; break;
        }
        colorR = r;
        colorG = g;
        colorB = b;
    }

    void display() {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glLoadIdentity();

        gluLookAt(0, 0, 5, 0, 0, 0, 0, 1, 0);

        glRotatef(angleX, 1, 0, 0);
        glRotatef(angleY, 0, 1, 0);
        glScalef(scale, scale, scale);

        if (showCube) drawCoordinateCube();

        GLfloat ambient[] = { 0.8f, 0.8f, 0.8f, 1.0f };
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);

        GLfloat diffuse[] = { 0.0f, 0.0f, 0.0f, 1.0f };
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
        GLfloat specular[] = { 0.0f, 0.0f, 0.0f, 1.0f };
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

        // --------- Draw Faces (no blend) ---------
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
        glEnable(GL_LIGHTING);

        for (const auto& face : faces) {
            const auto& verts = face.getVertices();
            if (verts.size() < 3) continue;

            // Compute normal for this face
            vec3 v1 = verts[1] - verts[0];
            vec3 v2 = verts[2] - verts[0];
            vec3 normal = v1.crossPdt(v2).normalization();
            glNormal3f(normal.x, normal.y, normal.z);
            glColor4f(colorR, colorG, colorB, 1.0f);

            glBegin(GL_POLYGON);
            for (const auto& v : verts) {
                glVertex3f(v.x, v.y, v.z);
            }
            glEnd();
        }

        // --------- Draw Edges ---------
        if (showEdges) {
            glEnable(GL_BLEND);
            glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

            // Pass 1: draw all edges ignoring depth, dashed, alpha*0.5
            glDisable(GL_DEPTH_TEST);
            glDisable(GL_LIGHTING);
            glEnable(GL_LINE_STIPPLE);
            glLineStipple(1, 0xFFFF); // dashed
            glColor4f(0.0f, 0.0f, 0.0f, transparency * 0.2f);
            glLineWidth(0.3f);

            glBegin(GL_LINES);
            for (const auto& face : faces) {
                const auto& verts = face.getVertices();
                for (size_t i = 0; i < verts.size(); ++i) {
                    const vec3& v1 = verts[i];
                    const vec3& v2 = verts[(i + 1) % verts.size()];
                    glVertex3f(v1.x, v1.y, v1.z);
                    glVertex3f(v2.x, v2.y, v2.z);
                }
            }
            glEnd();
            glDisable(GL_LINE_STIPPLE);

            // Pass 2: draw edges with depth, solid, alpha
            glEnable(GL_DEPTH_TEST);
            glLineWidth(2.0f);
            glColor4f(0.0f, 0.0f, 0.0f, transparency);

            glBegin(GL_LINES);
            for (const auto& face : faces) {
                const auto& verts = face.getVertices();
                for (size_t i = 0; i < verts.size(); ++i) {
                    const vec3& v1 = verts[i];
                    const vec3& v2 = verts[(i + 1) % verts.size()];
                    glVertex3f(v1.x, v1.y, v1.z);
                    glVertex3f(v2.x, v2.y, v2.z);
                }
            }
            glEnd();
            glEnable(GL_LIGHTING);
            glDepthMask(GL_TRUE);
            glDisable(GL_BLEND);
        }

        if (showColorControls) drawColorControls();

        glutSwapBuffers();
    }

    void drawColorControls() {
        glMatrixMode(GL_PROJECTION);
        glPushMatrix();
        glLoadIdentity();
        gluOrtho2D(0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), 0);

        glMatrixMode(GL_MODELVIEW);
        glPushMatrix();
        glLoadIdentity();

        glDisable(GL_LIGHTING);
        glDisable(GL_DEPTH_TEST);

        glColor4f(0.2f, 0.2f, 0.2f, 0.7f);
        glBegin(GL_QUADS);
        glVertex2f(colorMapStartX - 10, colorMapStartY - 10);
        glVertex2f(colorMapStartX + colorMapWidth + 10, colorMapStartY - 10);
        glVertex2f(colorMapStartX + colorMapWidth + 10, colorMapStartY + colorMapHeight + 10);
        glVertex2f(colorMapStartX - 10, colorMapStartY + colorMapHeight + 10);
        glEnd();

        glPointSize(1.0f);
        glBegin(GL_POINTS);
        for (int i = 0; i < colorMapWidth; ++i) {
            for (int j = 0; j < colorMapHeight; ++j) {
                float hue = i / float(colorMapWidth);
                float sat = 1.0f - (j / float(colorMapHeight));
                float val = 1.0f;
                float r, g, b;
                int hi = int(hue * 6);
                float f = hue * 6 - hi;
                float p = val * (1 - sat);
                float q = val * (1 - f * sat);
                float t = val * (1 - (1 - f) * sat);
                switch (hi % 6) {
                case 0: r = val, g = t, b = p; break;
                case 1: r = q, g = val, b = p; break;
                case 2: r = p, g = val, b = t; break;
                case 3: r = p, g = q, b = val; break;
                case 4: r = t, g = p, b = val; break;
                case 5: r = val, g = p, b = q; break;
                }
                glColor3f(r, g, b);
                glVertex2f(colorMapStartX + i, colorMapStartY + j);
            }
        }
        glEnd();

        if (pickerX >= 0 && pickerY >= 0) {
            float hue = pickerX / float(colorMapWidth);
            float sat = 1.0f - (pickerY / float(colorMapHeight));
            float val = 1.0f;
            float r, g, b;
            int hi = int(hue * 6);
            float f = hue * 6 - hi;
            float p = val * (1 - sat);
            float q = val * (1 - f * sat);
            float t = val * (1 - (1 - f) * sat);
            switch (hi % 6) {
            case 0: r = val, g = t, b = p; break;
            case 1: r = q, g = val, b = p; break;
            case 2: r = p, g = val, b = t; break;
            case 3: r = p, g = q, b = val; break;
            case 4: r = t, g = p, b = val; break;
            case 5: r = val, g = p, b = q; break;
            }
            glColor3f(r, g, b);
            glLineWidth(2.0f);
            int px = colorMapStartX + pickerX;
            int py = colorMapStartY + pickerY;
            glBegin(GL_LINE_LOOP);
            glVertex2f(px - 6, py - 6);
            glVertex2f(px + 6, py - 6);
            glVertex2f(px + 6, py + 6);
            glVertex2f(px - 6, py + 6);
            glEnd();
        }

        int transBarY = colorMapStartY + colorMapHeight + 25;
        int transBarH = 15, transBarW = colorMapWidth;
        glColor3f(0.7f, 0.7f, 0.7f);
        glBegin(GL_QUADS);
        glVertex2f(colorMapStartX, transBarY);
        glVertex2f(colorMapStartX + transBarW, transBarY);
        glVertex2f(colorMapStartX + transBarW, transBarY + transBarH);
        glVertex2f(colorMapStartX, transBarY + transBarH);
        glEnd();

        float knobX = colorMapStartX + transparency * transBarW;
        glColor3f(1.0f, 1.0f, 1.0f);
        glBegin(GL_QUADS);
        glVertex2f(knobX - 5, transBarY - 3);
        glVertex2f(knobX + 5, transBarY - 3);
        glVertex2f(knobX + 5, transBarY + transBarH + 3);
        glVertex2f(knobX - 5, transBarY + transBarH + 3);
        glEnd();

        glEnable(GL_DEPTH_TEST);
        glEnable(GL_LIGHTING);

        glMatrixMode(GL_PROJECTION);
        glPopMatrix();
        glMatrixMode(GL_MODELVIEW);
        glPopMatrix();
    }

    void drawCoordinateCube() {
        glDisable(GL_LIGHTING);
        glLineWidth(1.0f);

        glBegin(GL_LINES);
        glColor3f(1.0f, 0.0f, 0.0f);
        glVertex3f(-1.0f, 0.0f, 0.0f);
        glVertex3f(1.0f, 0.0f, 0.0f);

        glColor3f(0.0f, 1.0f, 0.0f);
        glVertex3f(0.0f, -1.0f, 0.0f);
        glVertex3f(0.0f, 1.0f, 0.0f);

        glColor3f(0.0f, 0.0f, 1.0f);
        glVertex3f(0.0f, 0.0f, -1.0f);
        glVertex3f(0.0f, 0.0f, 1.0f);
        glEnd();

        glColor3f(0.5f, 0.5f, 0.5f);
        glutWireCube(2.0f);

        glEnable(GL_LIGHTING);
    }

    void reshape(int w, int h) {
        glViewport(0, 0, w, h);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0, (double)w / (double)h, 0.1, 100.0);
        glMatrixMode(GL_MODELVIEW);
    }

    void mouse(int button, int state, int x, int y) {
        // 鼠标滚轮处理（GLUT环境下）
        if (state == GLUT_DOWN && (button == 3 || button == 4)) {
            if (button == 3)
                scale *= 1.05f;
            else
                scale /= 1.05f;
            if (scale < 0.01f) scale = 0.01f;
            if (scale > 100.0f) scale = 100.0f;
            glutPostRedisplay();
            return;
        }

        if (showColorControls && state == GLUT_DOWN) {
            if (x >= colorMapStartX && x <= colorMapStartX + colorMapWidth &&
                y >= colorMapStartY && y <= colorMapStartY + colorMapHeight) {
                pickColorFromMap(x, y);
                pickingColor = true;
                glutPostRedisplay();
                return;
            }
            int transBarY = colorMapStartY + colorMapHeight + 25;
            int transBarH = 15;
            if (x >= colorMapStartX && x <= colorMapStartX + colorMapWidth &&
                y >= transBarY && y <= transBarY + transBarH) {
                transparency = (x - colorMapStartX) / float(colorMapWidth);
                transparency = std::max(0.05f, std::min(1.0f, transparency));
                glutPostRedisplay();
                return;
            }
        }

        mouseX = x;
        mouseY = y;

        if (button == GLUT_LEFT_BUTTON) {
            mouseLeftDown = (state == GLUT_DOWN);
            if (state == GLUT_UP) pickingColor = false;
        }
        else if (button == GLUT_RIGHT_BUTTON) {
            mouseRightDown = (state == GLUT_DOWN);
        }
    }

    void motion(int x, int y) {
        if (showColorControls && pickingColor) {
            if (x >= colorMapStartX && x <= colorMapStartX + colorMapWidth &&
                y >= colorMapStartY && y <= colorMapStartY + colorMapHeight) {
                pickColorFromMap(x, y);
                glutPostRedisplay();
                return;
            }
        }

        if (mouseLeftDown && !pickingColor) {
            angleY += (x - mouseX);
            angleX += (y - mouseY);
            mouseX = x;
            mouseY = y;
        }
        if (mouseRightDown) {
            scale += (y - mouseY) * 0.01f;
            mouseY = y;
        }

        glutPostRedisplay();
    }

    static void run(int argc, char** argv,
        const std::vector<Face>& convex_poly_faces) {
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
        glutInitWindowSize(800, 600);
        glutCreateWindow("3D Convex Polyhedron Visualizer");

        static ConvexPolyhedronVisualizer visualizer;
        visualizer.initialize(convex_poly_faces);

        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

        GLfloat light_pos[] = { 0.0, 0.0, 1.0, 0.0 };
        glLightfv(GL_LIGHT0, GL_POSITION, light_pos);

        GLfloat ambient[] = { 0.8f, 0.8f, 0.8f, 1.0f };
        glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);

        GLfloat diffuse[] = { 0.0f, 0.0f, 0.0f, 1.0f };
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
        GLfloat specular[] = { 0.0f, 0.0f, 0.0f, 1.0f };
        glLightfv(GL_LIGHT0, GL_SPECULAR, specular);

        glutDisplayFunc([]() { visualizer.display(); });
        glutReshapeFunc([](int w, int h) { visualizer.reshape(w, h); });
        glutMouseFunc([](int b, int s, int x, int y) { visualizer.mouse(b, s, x, y); });
        glutMotionFunc([](int x, int y) { visualizer.motion(x, y); });

        glutMainLoop();
    }
};