
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <math.h>
#include <Eigen/Dense>

#include "util.h"
#include "parsing.h"

using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix4d;
using namespace std;

int main(int argc, char* argv[])
{
    Camera camera;
    vector<Light*> lights;
    vector<Object*> objects;
    parseSceneFile(argc, argv, camera, lights, objects);

    cout << objects[0]->toString() << endl;

    // Vector3d cameraPoint(50, 50, 0);
    // Ray viewRay(camera.eye, cameraPoint - camera.eye);
    // Ray bounceRay(Vector3d(0, 0, 0), Vector3d(0, 0, 0));

    // if (objects[0]->intersects(viewRay, bounceRay))
    //     cout << "INTERSECTION." << endl;
    // else
    //     cout << "Did not intersect." << endl;





    // for (int i = 0; i < 400; i++) {
    //     for (int j = 0; j < 400; j++) {
    //         // FIXME: Hard-coded value for testing. Compute using LL, LR, UL, UR later.
    //         Vector3d cameraPlanePosition(i + 0.5, j + 0.5, 0);
    //         Ray r(camera.eye, cameraPlanePosition - camera.eye);
    //     }
    // }

    // Tear down and deallocate data structures.
    for (Light* light : lights)
        free(light);

    for (Object* object : objects)
        free(object);
}
