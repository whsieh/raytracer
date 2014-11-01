
#ifndef __PARSING_H__
#define __PARSING_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <math.h>
#include <Eigen/Dense>

#include "util.h"

using Eigen::Affine;
using Eigen::Vector3f;
using Eigen::Quaternionf;
using Eigen::Transform;
using Eigen::Translation;
using Eigen::AngleAxisf;
using namespace std;

bool parseObjFile(string objFileName, vector<Object*>& objects, const Transform<float, 3, Affine>& currentTransform, const Material& currentMaterial, bool flipNormals = false)
{
    ifstream objFile;
    objFile.open(objFileName + string(".obj"), ios::binary | ios::in);
    if (!objFile.is_open()) {
        cout << "Error: could not open \"" << objFileName << ".obj\". Make sure the file exists." << endl;
        return false;
    }
    vector<Vector3f> vertices;
    vector<Vector3f> normals;
    int lineNo = 0;
    string currentLine;
    while (!objFile.eof()) {
        getline(objFile, currentLine);
        lineNo++;
        trim(currentLine);
        if (!currentLine.size())
            continue;

        if (currentLine[0] == '#')
            continue;

        vector<string> tokens = split(currentLine);
        string identifier = tokens[0];
        tokens.erase(tokens.begin());
        if (identifier == "v") {
            float x, y, z;
            if (tokens.size() < 3) {
                printf("Warning: not enough parameters for vertex on line %d of %s.obj.\n", lineNo, objFileName.c_str());
                continue;
            }
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z)) {
                printf("Warning: could not parse vertex on line %d of %s.obj.\n", lineNo, objFileName.c_str());
                continue;
            }
            vertices.push_back(Vector3f(x, y, z));

        } else if (identifier == "vn") {
            float x, y, z;
            if (tokens.size() < 3) {
                printf("Warning: not enough parameters for vertex normal on line %d of %s.obj.\n", lineNo, objFileName.c_str());
                continue;
            }
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z)) {
                printf("Warning: could not parse vertex normal on line %d of %s.obj.\n", lineNo, objFileName.c_str());
                continue;
            }
            normals.push_back(Vector3f(x, y, z));

        } else if (identifier == "f") {
            int v1, v2, v3, vn1 = 0, vn2 = 0, vn3 = 0;
            if (tokens.size() < 3) {
                printf("Warning: not enough parameters for face on line %d of %s.obj.\n", lineNo, objFileName.c_str());
                continue;
            }
            if (!parseIntPair(tokens[0], v1, vn1) || !parseIntPair(tokens[1], v2, vn2) || !parseIntPair(tokens[2], v3, vn3)) {
                printf("Warning: could not parse face on line %d of %s.obj.\n", lineNo, objFileName.c_str());
                continue;
            }
            v1--;
            v2--;
            v3--;
            vn1 = vn1 ? (vn1 - 1) : v1;
            vn2 = vn2 ? (vn2 - 1) : v2;
            vn3 = vn3 ? (vn3 - 1) : v3;
            if (v1 >= vertices.size() || v1 < 0) {
                printf("Warning: invalid vertex %d on line %d of %s.obj.\n", v1, lineNo, objFileName.c_str());
                continue;
            }
            if (v2 >= vertices.size() || v2 < 0) {
                printf("Warning: invalid vertex %d on line %d of %s.obj.\n", v2, lineNo, objFileName.c_str());
                continue;
            }
            if (v3 >= vertices.size() || v3 < 0) {
                printf("Warning: invalid vertex %d on line %d of %s.obj.\n", v3, lineNo, objFileName.c_str());
                continue;
            }
            if (normals.size() && (vn1 >= normals.size() || vn1 < 0)) {
                printf("Warning: invalid normal %d on line %d of %s.obj.\n", vn1, lineNo, objFileName.c_str());
                continue;
            }
            if (normals.size() && (vn2 >= normals.size() || vn2 < 0)) {
                printf("Warning: invalid normal %d on line %d of %s.obj.\n", vn2, lineNo, objFileName.c_str());
                continue;
            }
            if (normals.size() && (vn3 >= normals.size() || vn3 < 0)) {
                printf("Warning: invalid normal %d on line %d of %s.obj.\n", vn3, lineNo, objFileName.c_str());
                continue;
            }
            Triangle* triangle;
            if (normals.size() == 0)
                triangle = flipNormals ? new Triangle(vertices[v3], vertices[v2], vertices[v1]) : new Triangle(vertices[v1], vertices[v2], vertices[v3]);
            else
                triangle = flipNormals ? new Triangle(vertices[v3], vertices[v2], vertices[v1], normals[vn3], normals[vn2], normals[vn1])
                    : new Triangle(vertices[v1], vertices[v2], vertices[v3], normals[vn1], normals[vn2], normals[vn3]);

            triangle->setTransform(currentTransform);
            triangle->setMaterial(currentMaterial);
            objects.push_back(triangle);
        }
    }
    return true;
}

bool parseSceneFile(int argc, char* argv[], Camera& camera, vector<Light*>& lights, vector<Object*>& objects, SkyBox& environment, int& width, int& height, int& antiAliasing, int& reflectionDepth)
{
    string inputFileName = sceneFileNameFromArgs(argc, argv);
    if (!inputFileName.size())
        return false;

    ifstream sceneFile;
    sceneFile.open(inputFileName, ios::binary | ios::in);
    if (!sceneFile.is_open()) {
        cout << "Error: could not open \"" << inputFileName << "\". Make sure the file exists." << endl;
        return false;
    }
    string currentLine;
    unsigned int lineNo = 0;
    Transform<float, 3, Affine> currentTransform = Transform<float, 3, Affine>();
    resetTransformToIdentity(currentTransform);
    bool isTransformSet = false;

    Material currentMaterial;
    while (!sceneFile.eof()) {
        lineNo++;
        getline(sceneFile, currentLine);
        trim(currentLine);
        if (!currentLine.size())
            continue;

        if (currentLine[0] == '#')
            continue;

        vector<string> tokens = split(currentLine);
        string identifier = tokens[0];
        tokens.erase(tokens.begin());
        if (identifier == "dim") { // Image dimensions.
            int w, h;
            if (!parseInt(tokens[0], w) || !parseInt(tokens[1], h)) {
                printf("Warning: failed to parse image dimensions on line %d.\n", lineNo);
                continue;
            }
            width = w;
            height = h;

        } else if (identifier == "aa") { // Anti-aliasing.
            int a;
            if (!parseInt(tokens[0], a)) {
                printf("Warning: failed to parse anti-aliasing sample count on line %d.\n", lineNo);
                continue;
            }
            antiAliasing = a;

        } else if (identifier == "rd") { // Maximum reflection depth.
            int d;
            if (!parseInt(tokens[0], d)) {
                printf("Warning: failed to parse max reflection depth on line %d.\n", lineNo);
                continue;
            }
            reflectionDepth = d;

        } else if (identifier == "cam") {
            if (tokens.size() != 15) {
                printf("Warning: camera on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float x, y, z;
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z)) {
                printf("Warning: could not parse eye position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3f eye(x, y, z);
            if (!parseFloat(tokens[3], x) || !parseFloat(tokens[4], y) || !parseFloat(tokens[5], z)) {
                printf("Warning: could not parse LL position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3f ll(x, y, z);
            if (!parseFloat(tokens[6], x) || !parseFloat(tokens[7], y) || !parseFloat(tokens[8], z)) {
                printf("Warning: could not parse LR position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3f lr(x, y, z);
            if (!parseFloat(tokens[9], x) || !parseFloat(tokens[10], y) || !parseFloat(tokens[11], z)) {
                printf("Warning: could not parse UL position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3f ul(x, y, z);
            if (!parseFloat(tokens[12], x) || !parseFloat(tokens[13], y) || !parseFloat(tokens[14], z)) {
                printf("Warning: could not parse UR position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3f ur(x, y, z);
            camera = Camera(eye, ll, lr, ul, ur);

        } else if (identifier == "tri") {
            if (tokens.size() < 9) {
                printf("Warning: triangle on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float x, y, z;
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z)) {
                printf("Warning: could not parse vertex A for triangle on line %d.\n", lineNo);
                continue;
            }
            Vector3f vertexA(x, y, z);
            if (!parseFloat(tokens[3], x) || !parseFloat(tokens[4], y) || !parseFloat(tokens[5], z)) {
                printf("Warning: could not parse vertex B for triangle on line %d.\n", lineNo);
                continue;
            }
            Vector3f vertexB(x, y, z);
            if (!parseFloat(tokens[6], x) || !parseFloat(tokens[7], y) || !parseFloat(tokens[8], z)) {
                printf("Warning: could not parse vertex C for triangle on line %d.\n", lineNo);
                continue;
            }
            Vector3f vertexC(x, y, z);
            Triangle* triangle = new Triangle(vertexA, vertexB, vertexC);
            triangle->setTransform(currentTransform);
            triangle->setMaterial(currentMaterial);
            objects.push_back(triangle);

        } else if (identifier == "sph") {
            if (tokens.size() < 4) {
                printf("Warning: sphere on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float x, y, z, r;
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z) || !parseFloat(tokens[3], r)) {
                printf("Warning: could not parse sphere on line %d.\n", lineNo);
                continue;
            }
            Sphere* sphere = new Sphere(Vector3f(x, y, z), r);
            sphere->setTransform(currentTransform);
            sphere->setMaterial(currentMaterial);
            objects.push_back(sphere);

        } else if (identifier == "ltp") {
            if (tokens.size() < 6) {
                printf("Warning: point light on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float x, y, z, r, g, b;
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z) ||
                !parseFloat(tokens[3], r) || !parseFloat(tokens[4], g) || !parseFloat(tokens[5], b)) {
                printf("Warning: could not parse point light on line %d.\n", lineNo);
                continue;
            }
            int falloff = 0;
            if (tokens.size() >= 7) {
                if (!parseInt(tokens[6], falloff)) {
                    printf("Warning: could not parse falloff for point light on line %d.\n", lineNo);
                    continue;
                }
            }
            lights.push_back(new PointLight(Vector3f(x, y, z), Color(r, g, b), falloff));

        } else if (identifier == "ltd") {
            if (tokens.size() < 6) {
                printf("Warning: directional light on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float x, y, z, r, g, b;
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z) ||
                !parseFloat(tokens[3], r) || !parseFloat(tokens[4], g) || !parseFloat(tokens[5], b)) {
                printf("Warning: could not parse directional light on line %d.\n", lineNo);
                continue;
            }
            lights.push_back(new DirectionalLight(Vector3f(x, y, z), Color(r, g, b)));

        } else if (identifier == "lta") {
            if (tokens.size() < 3) {
                printf("Warning: ambient light on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float r, g, b;
            if (!parseFloat(tokens[0], r) || !parseFloat(tokens[1], g) || !parseFloat(tokens[2], b)) {
                printf("Warning: could not parse ambient light on line %d.\n", lineNo);
                continue;
            }
            lights.push_back(new AmbientLight(Color(r, g, b)));

        } else if (identifier == "mat") {
            if (tokens.size() < 13) {
                printf("Warning: material on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float r, g, b;
            if (!parseFloat(tokens[0], r) || !parseFloat(tokens[1], g) || !parseFloat(tokens[2], b)) {
                printf("Warning: could not parse ambient color for material on line %d.\n", lineNo);
                continue;
            }
            Color ambient(r, g, b);

            if (!parseFloat(tokens[3], r) || !parseFloat(tokens[4], g) || !parseFloat(tokens[5], b)) {
                printf("Warning: could not parse diffuse color for material on line %d.\n", lineNo);
                continue;
            }
            Color diffuse(r, g, b);

            if (!parseFloat(tokens[6], r) || !parseFloat(tokens[7], g) || !parseFloat(tokens[8], b)) {
                printf("Warning: could not parse specular color for material on line %d.\n", lineNo);
                continue;
            }
            Color specular(r, g, b);

            float specularCoefficient;
            if (!parseFloat(tokens[9], specularCoefficient)) {
                printf("Warning: could not parse specular coefficient for material on line %d.\n", lineNo);
                continue;
            }

            if (!parseFloat(tokens[10], r) || !parseFloat(tokens[11], g) || !parseFloat(tokens[12], b)) {
                printf("Warning: could not parse reflective color for material on line %d.\n", lineNo);
                continue;
            }
            Color reflective(r, g, b);
            currentMaterial = Material(ambient, diffuse, specular, reflective, specularCoefficient);

        } else if (identifier == "xft") {
            if (tokens.size() < 3) {
                printf("Warning: translation transformation on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float x, y, z;
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z)) {
                printf("Warning: could not parse translation transformation on line %d.\n", lineNo);
                continue;
            }
            if (isTransformSet)
                currentTransform = currentTransform * Translation<float, 3>(Vector3f(x, y, z));
            else {
                currentTransform = Translation<float, 3>(Vector3f(x, y, z));
                isTransformSet = true;
            }

        } else if (identifier == "xfs") {
            if (tokens.size() < 3) {
                printf("Warning: scale transformation on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float x, y, z;
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z)) {
                printf("Warning: could not parse scale transformation on line %d.\n", lineNo);
                continue;
            }
            if (isTransformSet)
                currentTransform = currentTransform * Scaling(Vector3f(x, y, z));
            else {
                currentTransform = Scaling(Vector3f(x, y, z));
                isTransformSet = true;
            }

        } else if (identifier == "xfr") {
            if (tokens.size() < 3) {
                printf("Warning: rotation transformation on line %d is missing parameters.\n", lineNo);
                continue;
            }
            float x, y, z;
            if (!parseFloat(tokens[0], x) || !parseFloat(tokens[1], y) || !parseFloat(tokens[2], z)) {
                printf("Warning: could not parse rotation transformation on line %d.\n", lineNo);
                continue;
            }
            if (isTransformSet)
                currentTransform = currentTransform * rotationTransformFromAxisAngle(x, y, z);
            else {
                currentTransform = rotationTransformFromAxisAngle(x, y, z);
                isTransformSet = true;
            }

        } else if (identifier == "xfz") {
            resetTransformToIdentity(currentTransform);
            isTransformSet = false;

        } else if (identifier == "obj") {
            if (tokens.size() < 1) {
                printf("Warning: obj specified without file name on line %d.\n", lineNo);
                continue;
            }
            if (!parseObjFile(tokens[0], objects, currentTransform, currentMaterial, tokens.size() >= 2 && tokens[2] == "-f"))
                printf("Warning: error occurred while parsing \"%s.obj\" on line %d.\n", tokens[0].c_str() , lineNo);

        } else if (identifier == "bg") {
            if (tokens.size() < 3) {
                printf("Warning: background on line %d is missing parameters.\n", lineNo);
                continue;
            }
            if (tokens.size() == 3) {
                float red, green, blue;
                if (!parseFloat(tokens[0], red) || !parseFloat(tokens[1], green) || !parseFloat(tokens[2], blue)) {
                    printf("Warning: failed to parse background color on line %d.\n", lineNo);
                    continue;
                }
                environment = SkyBox(Color(red, green, blue));

            } else {
                float minX, maxX, minY, maxY, minZ, maxZ;
                if (!parseFloat(tokens[1], minX) || !parseFloat(tokens[2], maxX) || !parseFloat(tokens[3], minY)
                    || !parseFloat(tokens[4], maxY) || !parseFloat(tokens[5], minZ) || !parseFloat(tokens[6], maxZ)) {
                    printf("Warning: failed to parse min/max for axes on line %d.\n", lineNo);
                    continue;
                }
                int bgWidth, bgHeight;
                if (!parseInt(tokens[7], bgWidth) || !parseInt(tokens[8], bgHeight)) {
                    printf("Warning: failed to parse image dimensions on line %d.\n", lineNo);
                    continue;
                }
                environment = SkyBox(tokens[0], minX, maxX, minY, maxY, minZ, maxZ, bgWidth, bgHeight);
            }
        }
    }
    sceneFile.close();
    return true;
}

#endif /** __PARSING_H__ */
