
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

using Eigen::Vector3d;
using Eigen::Vector4d;
using namespace std;

bool parseSceneFile(int argc, char* argv[], Camera& camera, vector<Light*>& lights, vector<Object*>& objects)
{
    string inputFileName;
    for (int i = 0; i < argc; i++)
        if (string(argv[i]).find(".scene") != -1)
            inputFileName = argv[i];

    if (!inputFileName.size())
        return false;

    std::ifstream sceneFile;
    sceneFile.open(inputFileName, ios::binary | ios::in);
    if (!sceneFile.is_open()) {
        std::cout << "Error: could not open \"" << inputFileName << "\". Make sure the file exists." << std::endl;
        return false;
    }
    string currentLine;
    unsigned int lineNo = 0;
    Transformation currentTransformation;
    Material currentMaterial;
    while (!sceneFile.eof()) {
        lineNo++;
        getline(sceneFile, currentLine);
        trim(currentLine);
        if (!currentLine.size())
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
            Vector3d eye(x, y, z);
            if (!parseFloat(tokens[3], x) || !parseFloat(tokens[4], y) || !parseFloat(tokens[5], z)) {
                printf("Warning: could not parse LL position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3d ll(x, y, z);
            if (!parseFloat(tokens[6], x) || !parseFloat(tokens[7], y) || !parseFloat(tokens[8], z)) {
                printf("Warning: could not parse LR position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3d lr(x, y, z);
            if (!parseFloat(tokens[9], x) || !parseFloat(tokens[10], y) || !parseFloat(tokens[11], z)) {
                printf("Warning: could not parse UL position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3d ul(x, y, z);
            if (!parseFloat(tokens[12], x) || !parseFloat(tokens[13], y) || !parseFloat(tokens[14], z)) {
                printf("Warning: could not parse UR position for camera on line %d.\n", lineNo);
                continue;
            }
            Vector3d ur(x, y, z);
            camera = Camera(eye, ll, lr, ul, ur);

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
            Sphere* sphere = new Sphere(Vector3d(x, y, z), r);
            sphere->transformation = currentTransformation;
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
            lights.push_back(new PointLight(Vector3d(x, y, z), Color(r, g, b)));

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
            lights.push_back(new DirectionalLight(Vector3d(x, y, z), Color(r, g, b)));

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
            Translation translation(x, y, z);
            // FIXME: Apply transforms in reverse order!
            currentTransformation.push(translation);

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
            Scale scale(x, y, z);
            // FIXME: Apply transforms in reverse order!
            currentTransformation.push(scale);

        } else if (identifier == "xfz") {
            currentTransformation.reset();
        }
    }

    sceneFile.close();
    return true;
}

#endif /** __PARSING_H__ */
