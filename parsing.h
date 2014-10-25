
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

bool parseObjFile(string objFileName, vector<Object*>& objects, const Transform<float, 3, Affine>& currentTransform, const Material& currentMaterial)
{
    ifstream objFile;
    objFile.open(objFileName + string(".obj"), ios::binary | ios::in);
    if (!objFile.is_open()) {
        cout << "Error: could not open \"" << objFileName << ".obj\". Make sure the file exists." << endl;
        return false;
    }
    vector<Vector3f> vertices;
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

        } else if (identifier == "f") {
            int v1, v2, v3;
            if (tokens.size() < 3) {
                printf("Warning: not enough parameters for face on line %d of %s.obj.\n", lineNo, objFileName.c_str());
                continue;
            }
            if (!parseInt(tokens[0], v1) || !parseInt(tokens[1], v2) || !parseInt(tokens[2], v3)) {
                printf("Warning: could not parse face on line %d of %s.obj.\n", lineNo, objFileName.c_str());
                continue;
            }
            v1--;
            v2--;
            v3--;
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
            Triangle* triangle = new Triangle(vertices[v1], vertices[v2], vertices[v3]);
            triangle->transform = currentTransform;
            triangle->material = currentMaterial;
            objects.push_back(triangle);
        }
    }
    return true;
}

bool parseSceneFile(int argc, char* argv[], Camera& camera, vector<Light*>& lights, vector<Object*>& objects)
{
    string inputFileName;
    for (int i = 0; i < argc; i++)
        if (string(argv[i]).find(".scene") != -1)
            inputFileName = argv[i];

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

        if (currentLine[0] == '/' && currentLine[1] == '/')
            continue;

        vector<string> tokens = split(currentLine);
        string identifier = tokens[0];
        tokens.erase(tokens.begin());

        if (identifier == "cam") {
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
            triangle->transform = currentTransform;
            triangle->material = currentMaterial;
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
            sphere->transform = currentTransform;
            sphere->material = currentMaterial;
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
            lights.push_back(new PointLight(Vector3f(x, y, z), Color(r, g, b)));

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

        } else if (identifier == "xfz") {
            resetTransformToIdentity(currentTransform);
            isTransformSet = false;

        } else if (identifier == "obj") {
            if (tokens.size() < 1) {
                printf("Warning: obj specified without file name on line %d.\n", lineNo);
                continue;
            }
            if (!parseObjFile(tokens[0], objects, currentTransform, currentMaterial))
                printf("Warning: error occurred while parsing \"%s.obj\" on line %d.\n", tokens[0].c_str() , lineNo);
        }
    }

    sceneFile.close();
    return true;
}

#endif /** __PARSING_H__ */
