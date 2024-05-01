#include <GL/glut.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cstdio>
#include <cmath>

int gNumVertices = 0;  
int gNumTriangles = 0; 
int*    gIndexBuffer = NULL; 
const double M_PI = 3.14159265358979323846;
float* vertices;

float ka[3] = { 0.0, 1.0, 0.0 }; 
float kd[3] = { 0.0, 0.5, 0.0 };
float ks[3] = { 0.5, 0.5, 0.5 }; 
float p = 32.0; 
float* colors;
float l=-0.1;
float r = 0.1;
float b = -0.1;
float t = 0.1;
float n = -0.1;
float f = -1000.0;

int nx = 512;
int ny = 512;
float Ia = 0.2;  
float lightPosition[3] = { -4.0, 4.0, -3.0 }; 
float lightColor[3] = { 1.0, 1.0, 1.0 };

float clearcolor[4] = { 0.0, 0.0,0.0,1.0 };



void create_scene()
{
    int width = nx;
    int height = ny;

    float theta, phi;
    int t;
    
    gNumVertices = (height - 2) * width + 2;

    vertices = new float[3 * gNumVertices];
    gNumTriangles = (height - 2) * (width - 1) * 2;
    colors = new float[3 * gNumVertices];

    gIndexBuffer = new int[3 * gNumTriangles];

    t = 0;
    for (int j = 1; j < height - 1; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            theta = (float)j / (height - 1) * M_PI;
            phi = (float)i / (width - 1) * M_PI * 2;

            float   x = sinf(theta) * cosf(phi);
            float   y = cosf(theta);
            float   z = -sinf(theta) * sinf(phi);

            vertices[3 * t] = x;
            vertices[3 * t + 1] = y;
            vertices[3 * t + 2] = z;
            t++;
        }
    }

    vertices[3 * t] = 0;
    vertices[3 * t + 1] = 1;
    vertices[3 * t + 2] = 0;

    t++;

    vertices[3 * t] = 0;
    vertices[3 * t + 1] = -1;
    vertices[3 * t + 2] = 0;

    t++;
    
    t = 0;
    for (int j = 0; j < height - 3; ++j)
    {
        for (int i = 0; i < width - 1; ++i)
        {
            gIndexBuffer[t++] = j * width + i;
            gIndexBuffer[t++] = (j + 1) * width + (i + 1);
            gIndexBuffer[t++] = j * width + (i + 1);
            gIndexBuffer[t++] = j * width + i;
            gIndexBuffer[t++] = (j + 1) * width + i;
            gIndexBuffer[t++] = (j + 1) * width + (i + 1);

        }
    }
    for (int i = 0; i < width - 1; ++i)
    {
        gIndexBuffer[t++] = (height - 2) * width;
        gIndexBuffer[t++] = i;
        gIndexBuffer[t++] = i + 1;
        gIndexBuffer[t++] = (height - 2) * width + 1;
        gIndexBuffer[t++] = (height - 3) * width + (i + 1);
        gIndexBuffer[t++] = (height - 3) * width + i;
    }

}
void normalize(float* v)
{
    float length = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] /= length;
    v[1] /= length;
    v[2] /= length;
}
void calculateNormal(int triangleIndex, float* normal)
{
    int i1 = gIndexBuffer[3 * triangleIndex];
    int i2 = gIndexBuffer[3 * triangleIndex + 1];
    int i3 = gIndexBuffer[3 * triangleIndex + 2];

    float v1[3], v2[3];
    v1[0] = vertices[3 * i2] - vertices[3 * i1];
    v1[1] = vertices[3 * i2 + 1] - vertices[3 * i1 + 1];
    v1[2] = vertices[3 * i2 + 2] - vertices[3 * i1 + 2];

    v2[0] = vertices[3 * i3] - vertices[3 * i1];
    v2[1] = vertices[3 * i3 + 1] - vertices[3 * i1 + 1];
    v2[2] = vertices[3 * i3 + 2] - vertices[3 * i1 + 2];

    normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
    normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
    normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

    float length = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);
    normal[0] /= length;
    normal[1] /= length;
    normal[2] /= length;
}
void calculateVertexNormals(float* normals)
{
    // Initialize all normals to zero
    for (int i = 0; i < gNumVertices; ++i)
    {
        normals[3 * i] = 0.0;
        normals[3 * i + 1] = 0.0;
        normals[3 * i + 2] = 0.0;
    }

    // Accumulate face normals to each vertex
    for (int i = 0; i < gNumTriangles; ++i)
    {
        float faceNormal[3];
        calculateNormal(i, faceNormal);

        for (int j = 0; j < 3; ++j)
        {
            int index = gIndexBuffer[3 * i + j];
            normals[3 * index] += faceNormal[0];
            normals[3 * index + 1] += faceNormal[1];
            normals[3 * index + 2] += faceNormal[2];
        }
    }

    // Normalize vertex normals
    for (int i = 0; i < gNumVertices; ++i)
    {
        float* normal = &normals[3 * i];
        normalize(normal);
    }
}
void calculateTriangleNormals(int triangleIndex, float* n1, float* n2, float* n3)
{
    int i1 = gIndexBuffer[3 * triangleIndex];
    int i2 = gIndexBuffer[3 * triangleIndex + 1];
    int i3 = gIndexBuffer[3 * triangleIndex + 2];

    float v1[3], v2[3];
    v1[0] = vertices[3 * i2] - vertices[3 * i1];
    v1[1] = vertices[3 * i2 + 1] - vertices[3 * i1 + 1];
    v1[2] = vertices[3 * i2 + 2] - vertices[3 * i1 + 2];

    v2[0] = vertices[3 * i3] - vertices[3 * i1];
    v2[1] = vertices[3 * i3 + 1] - vertices[3 * i1 + 1];
    v2[2] = vertices[3 * i3 + 2] - vertices[3 * i1 + 2];

    n1[0] = n2[0] = n3[0] = v1[1] * v2[2] - v1[2] * v2[1];
    n1[1] = n2[1] = n3[1] = v1[2] * v2[0] - v1[0] * v2[2];
    n1[2] = n2[2] = n3[2] = v1[0] * v2[1] - v1[1] * v2[0];

    normalize(n1);
    normalize(n2);
    normalize(n3);
}
void shadeVertex(float* vertex, float* normal, float* color)
{
    float L[3], V[3], R[3];

    // Light vector
    L[0] = lightPosition[0] - vertex[0];
    L[1] = lightPosition[1] - vertex[1];
    L[2] = lightPosition[2] - vertex[2];
    normalize(L);

    // View vector
    V[0] = -vertex[0];
    V[1] = -vertex[1];
    V[2] = -vertex[2];
    normalize(V);

    // Ambient component
    color[0] = Ia * ka[0];
    color[1] = Ia * ka[1];
    color[2] = Ia * ka[2];

    // Diffuse component
    float NdotL = normal[0] * L[0] + normal[1] * L[1] + normal[2] * L[2];
    if (NdotL > 0)
    {
        color[0] += lightColor[0] * kd[0] * NdotL;
        color[1] += lightColor[1] * kd[1] * NdotL;
        color[2] += lightColor[2] * kd[2] * NdotL;
    }

    // Specular component
    R[0] = 2 * NdotL * normal[0] - L[0];
    R[1] = 2 * NdotL * normal[1] - L[1];
    R[2] = 2 * NdotL * normal[2] - L[2];
    normalize(R);

    float RdotV = R[0] * V[0] + R[1] * V[1] + R[2] * V[2];
    if (RdotV > 0)
    {
        float spec = pow(RdotV, p);
        color[0] += lightColor[0] * ks[0] * spec;
        color[1] += lightColor[1] * ks[1] * spec;
        color[2] += lightColor[2] * ks[2] * spec;
    }
    int index = (vertex - vertices) / 3;
    if (index >= 0 && index < gNumVertices)
    {
        colors[3 * index] = color[0];
        colors[3 * index + 1] = color[1];
        colors[3 * index + 2] = color[2];
    }
}
void flatShading()
{
    //glEnable(GL_DEPTH_TEST);
    //glClear(GL_DEPTH_BUFFER_BIT);

    glShadeModel(GL_FLAT);
    for (int i = 0; i < gNumTriangles; ++i)
    {
        float centroid[3] = { 0.0, 0.0, 0.0 };
        for (int j = 0; j < 3; ++j)
        {
            int index = gIndexBuffer[3 * i + j];
            centroid[0] += vertices[3 * index];
            centroid[1] += vertices[3 * index + 1];
            centroid[2] += vertices[3 * index + 2];
        }
        centroid[0] /= 3.0;
        centroid[1] /= 3.0;
        centroid[2] /= 3.0;

        float normal[3];
        calculateNormal(i, normal);

        float color[3];
        shadeVertex(centroid, normal, color);

        for (int j = 0; j < 3; ++j)
        {
            int index = gIndexBuffer[3 * i + j];
            colors[3 * index] = color[0];
            colors[3 * index + 1] = color[1];
            colors[3 * index + 2] = color[2];
        }
    }
}

void gouraudShading()
{
    float* normals = new float[3 * gNumVertices];
    calculateVertexNormals(normals);

    for (int i = 0; i < gNumVertices; ++i)
    {
        float* vertex = &vertices[3 * i];
        float* normal = &normals[3 * i];

        float color[3];
        shadeVertex(vertex, normal, color);

        colors[3 * i] = color[0];
        colors[3 * i + 1] = color[1];
        colors[3 * i + 2] = color[2];
    }

    delete[] normals;
}


void phongShading()
{
    float* normals = new float[3 * gNumVertices];
    calculateVertexNormals(normals);

    for (int i = 0; i < gNumTriangles; ++i)
    {
        float n1[3], n2[3], n3[3];
        calculateTriangleNormals(i, n1, n2, n3);

        for (int j = 0; j < 3; ++j)
        {
            int index = gIndexBuffer[3 * i + j];
            float* vertex = &vertices[3 * index];
            float* normal;

            if (j == 0)
                normal = n1;
            else if (j == 1)
                normal = n2;
            else
                normal = n3;

            float color[3];
            shadeVertex(vertex, normal, color);

            colors[3 * index] = color[0];
            colors[3 * index + 1] = color[1];
            colors[3 * index + 2] = color[2];
        }
    }

    delete[] normals;
}


void gammaCorrection(float gamma)
{
    float gammaInv = 1.0 / gamma;

    // Apply gamma correction to each color component
    for (int i = 0; i < gNumVertices; ++i)
    {
        float* color = &colors[3 * i];  // colors 배열 사용

        color[0] = pow(color[0], gammaInv);
        color[1] = pow(color[1], gammaInv);
        color[2] = pow(color[2], gammaInv);
    }
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glFrustum(l, r, b, t, -n, -f);
    gluLookAt(0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0, 1.0, 0.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    glTranslatef(0.0f, 0.0f, -7.0f);
    glScalef(2.0f, 2.0f, 2.0f);

    glEnableClientState(GL_VERTEX_ARRAY);
    glVertexPointer(3, GL_FLOAT, 0, vertices);

    glEnableClientState(GL_COLOR_ARRAY);
    glColorPointer(3, GL_FLOAT, 0, colors);

    
    glEnable(GL_DEPTH_TEST);

    glDrawElements(GL_TRIANGLES, gNumTriangles * 3, GL_UNSIGNED_INT, gIndexBuffer);
    glDisableClientState(GL_VERTEX_ARRAY);
    flatShading();
    //gouraudShading();
    //phongShading();

    gammaCorrection(2.2);
    glutSwapBuffers();
}
int main(int argc, char** argv)
{
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(512, 512);
    glutCreateWindow("Sphere");

    create_scene();

    glutDisplayFunc(display);

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);

    glutMainLoop();
    delete[] gIndexBuffer;

    delete[] colors;

    return 0;
}