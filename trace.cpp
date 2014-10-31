
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <math.h>
#include <Eigen/Dense>

#include "util.h"
#include "parsing.h"

using Eigen::Affine;
using Eigen::Vector3f;
using Eigen::Quaternionf;
using Eigen::Transform;
using Eigen::Translation;
using Eigen::AngleAxisf;
using namespace std;

static int IMAGE_WIDTH = 400;
static int IMAGE_HEIGHT = 400;
static int MAX_REFLECTION_DEPTH = 3;
static int ANTI_ALIASING_SAMPLES = 0;

// Assumes that the shadow ray to the given light is not occluded.
Color computeDiffuseComponent(const Color& diffuseColor, Vector3f intersection, const Ray& normal, Light* light)
{
    float dotTerm = 0;
    if (light->type() == DIRECTIONAL) {
        DirectionalLight* directionalLight = (DirectionalLight*) light;
        dotTerm = normal.direction.normalized().dot(-directionalLight->direction().normalized());

    } else if (light->type() == POINT) {
        PointLight* pointLight = (PointLight*) light;
        dotTerm = normal.direction.normalized().dot((pointLight->position() - intersection).normalized());
    }
    // If dotTerm < 0, color clamping will naturally take care of it.
    return (diffuseColor * light->color()) * dotTerm;
}

Color computeSpecularComponent(const Color& specularColor, float specularCoefficient, Vector3f intersection, const Ray& reflection, Light* light)
{
    float dotTerm = 0;
    if (light->type() == DIRECTIONAL) {
        DirectionalLight* directionalLight = (DirectionalLight*) light;
        dotTerm = reflection.direction.normalized().dot(-directionalLight->direction().normalized());

    } else if (light->type() == POINT) {
        PointLight* pointLight = (PointLight*) light;
        dotTerm = reflection.direction.normalized().dot((pointLight->position() - intersection).normalized());
    }
    if (dotTerm < 0) // If this isn't here, taking dotTerm to the specularCoefficient might result in positive number.
        return Color(0, 0, 0);
    // If dotTerm < 0, color clamping will naturally take care of it.
    return (specularColor * light->color()) * pow(dotTerm, specularCoefficient);
}

void setImagePixelFromColor(vector<unsigned char>& image, unsigned int i, unsigned int j, unsigned int width, unsigned int height, Color color)
{
    unsigned int startingIndex = 4 * ((width * j) + i);
    image[startingIndex] = round(color.red * 255);
    image[startingIndex + 1] = round(color.green * 255);
    image[startingIndex + 2] = round(color.blue * 255);
    image[startingIndex + 3] = 255;
}

Color colorFromRay(const Ray& ray, const vector<Light*>& lights, const vector<Object*>& objects, const SkyBox& environment, const AABBNode* aabb, int depth = 0)
{
    float closestTIntersection = FLT_MAX;
    Vector3f intersection, currentIntersection;
    Ray bounce, normal, currentBounce, currentNormal;
    Object* intersectionObject;
    vector<Object*> potentialIntersectingObjects;
    aabb->collectObjectsForRayIntersection(ray, potentialIntersectingObjects);
    // cout << "Narrowed search to " << potentialIntersectingObjects.size() << " objects." << endl;
    for (Object* object : potentialIntersectingObjects) {
        if (object->intersects(ray, currentIntersection, currentNormal, currentBounce)) {
            float currentTIntersection = ray.timeAtPosition(currentIntersection);
            if (currentTIntersection < closestTIntersection) {
                closestTIntersection = currentTIntersection;
                intersection = currentIntersection;
                bounce = currentBounce;
                normal = currentNormal;
                intersectionObject = object;
            }
        }
    }
    // Return early if ray did not intersect with any objects.
    // Instead of defaulting to black, compute the background color using the SkyBox.
    if (closestTIntersection == FLT_MAX)
        return environment.colorFromRay(ray);

    Color pixelColor;
    Ray shadowRay;
    for (Light* light : lights) {
        pixelColor = pixelColor + (intersectionObject->material.ambient * light->color());
        if (light->type() == AMBIENT)
            continue;

        float maxTime = FLT_MAX;
        if (light->type() == DIRECTIONAL) {
            DirectionalLight* directionalLight = (DirectionalLight*) light;
            shadowRay.start = intersection;
            shadowRay.direction = -directionalLight->direction();

        } else if (light->type() == POINT) {
            PointLight* pointLight = (PointLight*) light;
            shadowRay.start = intersection;
            shadowRay.direction = pointLight->position() - intersection;
            maxTime = 1;
        }

        bool isOccluded = false;
        vector<Object*> potentialShadowRayIntersectingObjects;
        aabb->collectObjectsForRayIntersection(shadowRay, potentialShadowRayIntersectingObjects);
        for (Object* object : potentialShadowRayIntersectingObjects) {
            if (object->intersects(shadowRay, maxTime)) {
                isOccluded = true;
                break;
            }
        }
        if (isOccluded)
            continue;

        pixelColor = pixelColor + computeSpecularComponent(intersectionObject->material.specular, intersectionObject->material.specularCoefficient, intersection, bounce, light)
            + computeDiffuseComponent(intersectionObject->material.diffuse, intersection, normal, light);
    }
    Color reflectiveColor = intersectionObject->material.reflective;
    if (depth < MAX_REFLECTION_DEPTH && (reflectiveColor.red > 0 || reflectiveColor.green > 0 || reflectiveColor.blue > 0))
        pixelColor = pixelColor + (reflectiveColor * colorFromRay(bounce, lights, objects, environment, aabb, depth + 1));

    return pixelColor;
}

int main(int argc, char* argv[])
{
    long startTime = time(NULL);
    Camera camera;
    vector<Light*> lights;
    vector<Object*> objects;
    SkyBox environment;
    parseSceneFile(argc, argv, camera, lights, objects, environment, IMAGE_WIDTH, IMAGE_HEIGHT, ANTI_ALIASING_SAMPLES, MAX_REFLECTION_DEPTH);
    AABBNode* aabb = makeAABBNode(objects);
    cout << "Preparing to render scene. (" << objects.size() << " objects)" << endl;
    unsigned width = IMAGE_WIDTH, height = IMAGE_HEIGHT;
    std::vector<unsigned char> image;
    image.resize(4 * width * height);
    double antiAliasingWidthMargin = 0.4 / width;
    double antiAliasingHeightMargin = 0.4 / height;
    double u, v;
    for (int j = 0; j < IMAGE_HEIGHT; j++) {
        v = (j + 0.5) / (double) height;
        for (int i = 0; i < IMAGE_WIDTH; i++) {
            u = (i + 0.5) / (double) width;

            Vector3f cameraPlanePosition = camera.viewPlanePositionFrom2D(u, v);
            Ray viewRay = Ray(camera.eye, cameraPlanePosition - camera.eye);
            Color pixelColor = colorFromRay(viewRay, lights, objects, environment, aabb);
            float totalRed = pixelColor.red;
            float totalGreen = pixelColor.green;
            float totalBlue = pixelColor.blue;
            for (int k = 0; k < ANTI_ALIASING_SAMPLES; k++) {
                cameraPlanePosition = camera.viewPlanePositionFrom2D(u + randomDouble() * antiAliasingWidthMargin, v + randomDouble() * antiAliasingHeightMargin);
                viewRay = Ray(camera.eye, cameraPlanePosition - camera.eye);
                pixelColor = colorFromRay(viewRay, lights, objects, environment, aabb);
                totalRed += pixelColor.red;
                totalGreen += pixelColor.green;
                totalBlue += pixelColor.blue;
            }
            Color finalColor(totalRed / (1 + ANTI_ALIASING_SAMPLES), totalGreen / (1 + ANTI_ALIASING_SAMPLES), totalBlue / (1 + ANTI_ALIASING_SAMPLES));
            setImagePixelFromColor(image, i, j, width, height, finalColor);
        }
    }

    cout << "Rendering finished. Time: " << setprecision(3) << (time(NULL) - startTime) << " seconds." << endl;

    string label = sceneFileNameFromArgs(argc, argv, false);
    if (lodepng::encode(string(label + ".png").c_str(), image, width, height))
        cout << "Error: Failed to save image." << endl;

    // Tear down and deallocate data structures.
    delete aabb;
    for (Light* light : lights)
        free(light);

    for (Object* object : objects)
        free(object);
}
