
#ifndef __UTIL_H__
#define __UTIL_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
#include <math.h>
#include <Eigen/Dense>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include <cfloat>
#include <iomanip>

#include "lib/lodepng.cpp"

using Eigen::Affine;
using Eigen::Vector3f;
using Eigen::Matrix3f;
using Eigen::Quaternionf;
using Eigen::Transform;
using Eigen::Translation;
using Eigen::AngleAxisf;
using namespace std;

#define PI 3.14159265
#define EPSILON 0.01
#define MAX_AABB_DEPTH 20
#define AABB_LEAF_SIZE 10

static int GLOBAL_OBJECT_COUNT = 0;

enum LightSourceType { DIRECTIONAL, POINT, AMBIENT, UNSPECIFIED };
enum TransformationType { ROTATION_ONLY, TRANSLATION_ONLY, SCALE_ONLY, COMBINATION };

double randomDouble()
{
    return ((double)rand()) / RAND_MAX;
}

unsigned int randomInt(unsigned int lower, unsigned int upper)
{
    if (upper > lower)
        return (rand() % (upper - lower)) + lower;
    return -1;
}

void resetTransformToIdentity(Transform<float, 3, Affine>& transform)
{
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            transform.matrix()(i, j) = i == j ? 1 : 0;
        }
    }
}

string sceneFileNameFromArgs(int argc, char* argv[], bool includeFormat = true)
{
    string inputFileName;
    bool isFinalInput = false;
    for (int i = 0; i < argc; i++) {
        isFinalInput = string(argv[i]).find("input-") != -1;
        if (string(argv[i]).find(".scene") != -1 || isFinalInput) {
            inputFileName = argv[i];
            break;
        }
    }
    if (!includeFormat && !isFinalInput)
        inputFileName = inputFileName.substr(0, inputFileName.size() - 6);

    return inputFileName;
}

void loadImageToVector(vector<unsigned char>& image, const char* fileName, unsigned width, unsigned height)
{
    if (lodepng::decode(image, width, height, fileName)) {
        cout << "Could not load image " << string(fileName) << endl;
        exit(0);
    }
}

float square(float x) { return x * x; }

template<typename T>
void clamp(T& value, T lowerBound, T upperBound)
{
    if (value < lowerBound)
        value = lowerBound;

    else if (value > upperBound)
        value = upperBound;
}

bool parseFloat(string str, float& value)
{
    try {
        stringstream ss(str);
        ss >> value;
        if (!ss.eof())
            return false;
        return true;
    } catch (int e) {
        return false;
    }
}

bool parseInt(string str, int& value)
{
    try {
        stringstream ss(str);
        ss >> value;
        if (!ss.eof())
            return false;
        return true;
    } catch (int e) {
        return false;
    }
}

Transform<float, 3, Affine> rotationTransformFromAxisAngle(float x, float y, float z)
{
    float rotationAmount = PI * sqrt(square(x) + square(y) + square(z)) / 180.0;
    return Eigen::Transform<float, 3, Affine>(AngleAxisf(rotationAmount, Vector3f(x, y, z).normalized()).toRotationMatrix());
}

/** Boilerplate C++ string manipulation routines from the following sources:
 *      http://stackoverflow.com/questions/216823/whats-the-best-way-to-trim-stdstring
 *      http://stackoverflow.com/questions/236129/how-to-split-a-string-in-c
 */
static inline string &ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
    return s;
}

static inline string &rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
    return s;
}

static inline string &trim(string &s) {
    return ltrim(rtrim(s));
}

vector<string> split(const string& str, const string& delimiter = " ") {
    vector <string> tokens;
    string::size_type lastPos = 0;
    string::size_type pos = str.find(delimiter, lastPos);
    while (string::npos != pos) {
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        lastPos = pos + delimiter.size();
        pos = str.find(delimiter, lastPos);
    }
    tokens.push_back(str.substr(lastPos, str.size() - lastPos));
    return tokens;
}
/** End code that I didn't write. */

bool parseIntPair(string str, int& first, int& second, const string& delimiter = "//")
{
    vector<string> list = split(str, delimiter);
    if (list.size() == 0)
        return false;
    else if (list.size() == 1)
        return parseInt(list[0], first);
    return parseInt(list[0], first) && parseInt(list[1], second);
}

string vector3fAsString(Vector3f vector)
{
    ostringstream out;
    out << "[" << vector(0) << ", " << vector(1) << ", " << vector(2) << "]";
    return out.str();
}

/** The Color class is taken from my implementation of AS1. */
class Color
{
public:
    Color()
        : red(0)
        , green(0)
        , blue(0)
    {
    }

    Color(float r, float g, float b)
        : red(r)
        , green(g)
        , blue(b)
    {
        clamp<float>(red, 0, 1);
        clamp<float>(green, 0, 1);
        clamp<float>(blue, 0, 1);
    }

    string toString() const
    {
        return "<" + to_string(red) + ", " + to_string(green) + ", " + to_string(blue) + ">";
    }

    float red;
    float green;
    float blue;
};

inline Color operator*(const Color& c1, const Color& c2)
{
    return Color(c1.red * c2.red, c1.green * c2.green, c1.blue * c2.blue);
}

inline Color operator*(const Color& c, float scale)
{
    return Color(c.red * scale, c.green * scale, c.blue * scale);
}

inline Color operator*(float scale, const Color& c)
{
    return c * scale;
}

inline Color operator+(const Color& c1, const Color& c2)
{
    return Color(c1.red + c2.red, c1.green + c2.green, c1.blue + c2.blue);
}

Color colorFromImageVector(const vector<unsigned char>& image,  int width, int height, float u, float v)
{
    clamp<float>(u, 0, 1);
    clamp<float>(v, 0, 1);
    int x = min<int>(width - 1, u * width),
        y = min<int>(height - 1, v * height);
    int startIndex = 4 * (width * y + x);
    return Color(image[startIndex] / 255.0,
        image[startIndex + 1] / 255.0,
        image[startIndex + 2] / 255.0);
}

class Ray {
public:
    Ray(Vector3f _start, Vector3f _direction)
        : start(_start(0), _start(1), _start(2))
        , direction(_direction(0), _direction(1), _direction(2))
    {
    }

    Ray()
        : start(0, 0, 0)
        , direction(0, 0, 0)
    {
    }

    // Assumes that position lies somewhere on this ray. Averages to counteract floating point error.
    // Used to determine the "earliest" point of intersection with some set of objects.
    float timeAtPosition(Vector3f position) const
    {
        float t = 0;
        if (direction(0))
            t += (position(0) - start(0)) / direction(0);

        if (direction(1))
            t += (position(1) - start(1)) / direction(1);

        if (direction(2))
            t += (position(2) - start(2)) / direction(2);

        return t / 3;
    }

    Vector3f positionAtTime(float t) const
    {
        return start + (t * direction);
    }

    string toString() const
    {
        ostringstream out;
        out << "Ray(start=[" << setprecision(3) << start(0) << ", "
            << setprecision(3) << start(1) << ", "
            << setprecision(3) << start(2) << "], direction=["
            << setprecision(3) << direction(0) << ", "
            << setprecision(3) << direction(1) << ", "
            << setprecision(3) << direction(2) << "])";
        return out.str();
    }

    Vector3f start;
    Vector3f direction;
};

class AABB {
public:
    AABB(float _xmin, float _xmax, float _ymin, float _ymax, float _zmin, float _zmax)
        : xmin(_xmin)
        , xmax(_xmax)
        , ymin(_ymin)
        , ymax(_ymax)
        , zmin(_zmin)
        , zmax(_zmax)
    {
    }

    AABB()
        : xmin(-FLT_MAX)
        , xmax(FLT_MAX)
        , ymin(-FLT_MAX)
        , ymax(FLT_MAX)
        , zmin(-FLT_MAX)
        , zmax(FLT_MAX)
    {
    }

    string toString() const
    {
        ostringstream out;
        out << "AABB(x=[" << setprecision(3) << xmin << ", " << setprecision(3) << xmax << "], "
            << "y=[" << setprecision(3) << ymin << ", " << setprecision(3) << ymax << "], "
            << "z=[" << setprecision(3) << zmin << ", " << setprecision(3) << zmax << "])";
        return out.str();
    }

    Vector3f midpoint() const
    {
        return Vector3f((xmin + xmax) / 2, (ymin + ymax) / 2, (zmin + zmax) / 2);
    }

    bool intersects(const Ray& ray) const
    {
        Vector3f intersection;
        float tmin = (xmin - ray.start(0)) / ray.direction(0);
        float tmax = (xmax - ray.start(0)) / ray.direction(0);
        if (!isnan(tmin) && tmin > 0) {
            intersection = ray.positionAtTime(tmin);
            if (ymin <= intersection(1) && intersection(1) <= ymax && zmin <= intersection(2) && intersection(2) <= zmax)
                return true;
        }
        if (!isnan(tmax) && tmax > 0) {
            intersection = ray.positionAtTime(tmax);
            if (ymin <= intersection(1) && intersection(1) <= ymax && zmin <= intersection(2) && intersection(2) <= zmax)
                return true;
        }
        tmin = (ymin - ray.start(1)) / ray.direction(1);
        tmax = (ymax - ray.start(1)) / ray.direction(1);
        if (!isnan(tmin) && tmin > 0) {
            intersection = ray.positionAtTime(tmin);
            if (xmin <= intersection(0) && intersection(0) <= xmax && zmin <= intersection(2) && intersection(2) <= zmax)
                return true;
        }
        if (!isnan(tmax) && tmax > 0) {
            intersection = ray.positionAtTime(tmax);
            if (xmin <= intersection(0) && intersection(0) <= xmax && zmin <= intersection(2) && intersection(2) <= zmax)
                return true;
        }
        tmin = (zmin - ray.start(2)) / ray.direction(2);
        tmax = (zmax - ray.start(2)) / ray.direction(2);
        if (!isnan(tmin) && tmin > 0) {
            intersection = ray.positionAtTime(tmin);
            if (xmin <= intersection(0) && intersection(0) <= xmax && ymin <= intersection(1) && intersection(1) <= ymax)
                return true;
        }
        if (!isnan(tmax) && tmax > 0) {
            intersection = ray.positionAtTime(tmax);
            if (xmin <= intersection(0) && intersection(0) <= xmax && ymin <= intersection(1) && intersection(1) <= ymax)
                return true;
        }
        return false;
    }

    float xmin;
    float xmax;
    float ymin;
    float ymax;
    float zmin;
    float zmax;
};

class Camera {
public:
    Camera()
        : eye(0, 0, -1)
        , ll(-1, -1, 0)
        , lr(1, -1, 0)
        , ul(-1, 1, 0)
        , ur(1, 1, 0)
    {
    }

    Camera(Vector3f _eye, Vector3f _ll, Vector3f _lr, Vector3f _ul, Vector3f _ur)
        : eye(_eye)
        , ll(_ll)
        , lr(_lr)
        , ul(_ul)
        , ur(_ur)
    {
    }

    Vector3f viewPlanePositionFrom2D(float u, float v)
    {
        return u * (v * lr + (1 - v) * ur) + (1 - u) * (v * ll + (1 - v) * ul);
    }

    string toString() const
    {
        ostringstream out;
        out << "Camera(eye=[" << setprecision(3) << eye(0) << ", "
            << setprecision(3) << eye(1) << ", "
            << setprecision(3) << eye(2) << "], ll=["
            << setprecision(3) << ll(0) << ", "
            << setprecision(3) << ll(1) << ", "
            << setprecision(3) << ll(2) << "], lr=["
            << setprecision(3) << lr(0) << ", "
            << setprecision(3) << lr(1) << ", "
            << setprecision(3) << lr(2) << "], ul=["
            << setprecision(3) << ul(0) << ", "
            << setprecision(3) << ul(1) << ", "
            << setprecision(3) << ul(2) << "], ur=["
            << setprecision(3) << ur(0) << ", "
            << setprecision(3) << ur(1) << ", "
            << setprecision(3) << ur(2) << "])";
        return out.str();
    }

    Vector3f eye;
    Vector3f ll;
    Vector3f lr;
    Vector3f ul;
    Vector3f ur;
};

class Light
{
public:
    Light(Vector3f vector, Color color)
        : m_vector(vector)
        , m_color(color)
    {
    }

    virtual const LightSourceType type() const
    {
        return m_type;
    }

    const Color color() const
    {
        return m_color;
    }

protected:
    Vector3f m_vector;
    LightSourceType m_type;

private:
    Color m_color;
};

class DirectionalLight : public Light
{
public:
    DirectionalLight(Vector3f vector, Color color)
        : Light(vector, color)
    {
        m_type = DIRECTIONAL;
        m_vector = m_vector / sqrt(m_vector.dot(m_vector));
    }

    const Vector3f direction() const
    {
        return m_vector;
    }
};

class PointLight : public Light
{
public:
    PointLight(Vector3f vector, Color color)
        : Light(vector, color)
    {
        m_type = POINT;
    }

    const Vector3f position() const
    {
        return m_vector;
    }
};

class AmbientLight : public Light
{
public:
    AmbientLight(Color color)
        : Light(Vector3f(0, 0, 0), color)
    {
        m_type = AMBIENT;
    }
};

class Material {
public:
    Material()
        : ambient(Color(0, 0, 0))
        , diffuse(Color(0, 0, 0))
        , specular(Color(0, 0, 0))
        , reflective(Color(0, 0, 0))
        , specularCoefficient(1)
    {
    }

    Material(Color _ambient, Color _diffuse, Color _specular, Color _reflective, int _specularCoefficient)
        : ambient(_ambient)
        , diffuse(_diffuse)
        , specular(_specular)
        , reflective(_reflective)
        , specularCoefficient(_specularCoefficient)
    {
    }

    Color ambient;
    Color diffuse;
    Color specular;
    Color reflective;
    int specularCoefficient;
};

class Object {
public:
    Object(int _id)
        : id(_id)
    {
        resetTransformToIdentity(transform);
    }

    Ray rayInObjectSpace(const Ray& viewRay) const
    {
        Transform<float, 3, Affine> inverseTransform = transform.inverse();
        Vector3f transformedStart = inverseTransform * viewRay.start;
        Vector3f transformedOffsetPosition = inverseTransform * (viewRay.start + viewRay.direction);
        return Ray(transformedStart, transformedOffsetPosition - transformedStart);
    }

    virtual string toString() const = 0;
    virtual bool intersects(const Ray&, float) const = 0;
    virtual bool intersects(const Ray&, Vector3f&, Ray&, Ray&) const = 0;
    virtual void computeWorldAABB() = 0;

    void setTransform(const Transform<float, 3, Affine>& _transform)
    {
        transform = _transform;
        computeWorldAABB();
        midpoint = worldAABB.midpoint();
    }

    void setMaterial(const Material& _material)
    {
        material = _material;
    }

    AABB worldAABB;
    Vector3f midpoint;
    Transform<float, 3, Affine> transform;
    Material material;
    int id;
};

class Sphere : public Object {
public:
    Sphere(Vector3f _center, float _radius)
        : Object(GLOBAL_OBJECT_COUNT++)
        , center(_center)
        , radius(_radius)
    {
    }

    void computeWorldAABB()
    {
        Vector3f corners[8];
        // Compute the corners of the bounding box in object space and transform them into world space.
        corners[0] = transform * (center + Vector3f(-radius, -radius, -radius));
        corners[1] = transform * (center + Vector3f(-radius, -radius, radius));
        corners[2] = transform * (center + Vector3f(-radius, radius, -radius));
        corners[3] = transform * (center + Vector3f(-radius, radius, radius));
        corners[4] = transform * (center + Vector3f(radius, -radius, -radius));
        corners[5] = transform * (center + Vector3f(radius, -radius, radius));
        corners[6] = transform * (center + Vector3f(radius, radius, -radius));
        corners[7] = transform * (center + Vector3f(radius, radius, radius));
        // Compute a new bounding box around the transformed corners.
        float xmin = FLT_MAX, xmax = -FLT_MAX, ymin = FLT_MAX, ymax = -FLT_MAX, zmin = FLT_MAX, zmax = -FLT_MAX;
        for (int i = 0; i < 8; i++) {
            float x = corners[i](0),
                y = corners[i](1),
                z = corners[i](2);
            if (x < xmin)
                xmin = x;
            if (x > xmax)
                xmax = x;
            if (y < ymin)
                ymin = y;
            if (y > ymax)
                ymax = y;
            if (z < zmin)
                zmin = z;
            if (z > zmax)
                zmax = z;
        }
        if (xmin == FLT_MAX)
            xmin = xmax - 0.001;
        if (xmax == -FLT_MAX)
            xmax = xmin + 0.001;
        if (ymin == FLT_MAX)
            ymin = ymax - 0.001;
        if (ymax == -FLT_MAX)
            ymax = ymin + 0.001;
        if (zmin == FLT_MAX)
            zmin = zmax - 0.001;
        if (zmax == -FLT_MAX)
            zmax = zmin + 0.001;
        worldAABB = AABB(xmin, xmax, ymin, ymax, zmin, zmax);
    }

    bool intersects(const Ray& originalViewRay, float maxTime = FLT_MAX) const
    {
        Ray viewRayInObjectSpace = rayInObjectSpace(originalViewRay);
        float A = square(viewRayInObjectSpace.direction(0)) + square(viewRayInObjectSpace.direction(1)) + square(viewRayInObjectSpace.direction(2));
        float B = 2 * ((viewRayInObjectSpace.direction(0) * (viewRayInObjectSpace.start(0) - center(0))) +
                    (viewRayInObjectSpace.direction(1) * (viewRayInObjectSpace.start(1) - center(1))) +
                    (viewRayInObjectSpace.direction(2) * (viewRayInObjectSpace.start(2) - center(2))));
        float C = square(viewRayInObjectSpace.start(0) - center(0)) +
            square(viewRayInObjectSpace.start(1) - center(1)) +
            square(viewRayInObjectSpace.start(2) - center(2)) - square(radius);
        float discriminant = square(B) - (4 * A * C);
        float t1 = -1, t2 = -1;
        if (discriminant >= 0) {
            float rootTerm = sqrt(square(B) - (4 * A * C));
            t1 = (-B - rootTerm) / (2 * A);
            t2 = (-B + rootTerm) / (2 * A);
        }
        if (discriminant < 0 || (t1 <= EPSILON && t2 <= EPSILON))
            return false;

        return (t1 > EPSILON ? t1 : t2) < maxTime;
    }

    // Arguments intersection, normal and bounce are set in terms of world coordinates.
    bool intersects(const Ray& originalViewRay, Vector3f& intersection, Ray& normal, Ray& bounce) const
    {
        Ray viewRayInObjectSpace = rayInObjectSpace(originalViewRay);
        float A = square(viewRayInObjectSpace.direction(0)) + square(viewRayInObjectSpace.direction(1)) + square(viewRayInObjectSpace.direction(2));
        float B = 2 * ((viewRayInObjectSpace.direction(0) * (viewRayInObjectSpace.start(0) - center(0))) +
                    (viewRayInObjectSpace.direction(1) * (viewRayInObjectSpace.start(1) - center(1))) +
                    (viewRayInObjectSpace.direction(2) * (viewRayInObjectSpace.start(2) - center(2))));
        float C = square(viewRayInObjectSpace.start(0) - center(0)) +
            square(viewRayInObjectSpace.start(1) - center(1)) +
            square(viewRayInObjectSpace.start(2) - center(2)) - square(radius);
        float discriminant = square(B) - (4 * A * C);
        float t1 = -1, t2 = -1;
        if (discriminant >= 0) {
            float rootTerm = sqrt(square(B) - (4 * A * C));
            t1 = (-B - rootTerm) / (2 * A);
            t2 = (-B + rootTerm) / (2 * A);
        }
        if (discriminant < 0 || (t1 < EPSILON && t2 < EPSILON))
            return false; // The object does not intersect.

        float tIntersection = t1 > EPSILON ? t1 : t2;
        Vector3f intersectionPoint = viewRayInObjectSpace.positionAtTime(tIntersection);
        Vector3f incidenceNormal = (intersectionPoint - center).normalized();
        Vector3f incomingIncidenceVector = viewRayInObjectSpace.direction.normalized();
        Vector3f outgoingReflectionVector = (-2 * incomingIncidenceVector.dot(incidenceNormal)) * incidenceNormal + incomingIncidenceVector;

        intersection = transform * intersectionPoint;
        Vector3f bounceOffset = transform * (intersectionPoint + outgoingReflectionVector);
        Vector3f normalOffset = transform * (intersectionPoint + incidenceNormal);
        bounce.start = intersection;
        bounce.direction = bounceOffset - intersection;
        normal.start = intersection;
        normal.direction = normalOffset - intersection;
        return true;
    }

    string toString() const
    {
        ostringstream out;
        out << "Sphere(center=[" << setprecision(3) << center(0) << ", " << setprecision(3) << center(1) << ", "
            << setprecision(3) << center(2) << "], radius=" << setprecision(3) << radius << ")";
        return out.str();
    }

    Vector3f center;
    float radius;
};

class Triangle : public Object {
public:
    Triangle(Vector3f _a, Vector3f _b, Vector3f _c, Vector3f _normalA, Vector3f _normalB, Vector3f _normalC)
        : Object(GLOBAL_OBJECT_COUNT++)
        , a(_a)
        , b(_b)
        , c(_c)
        , normalA(_normalA)
        , normalB(_normalB)
        , normalC(_normalC)
    {
        cachedIntersectionMatrix(0, 0) = b(0) - a(0);
        cachedIntersectionMatrix(1, 0) = b(1) - a(1);
        cachedIntersectionMatrix(2, 0) = b(2) - a(2);
        cachedIntersectionMatrix(0, 1) = c(0) - a(0);
        cachedIntersectionMatrix(1, 1) = c(1) - a(1);
        cachedIntersectionMatrix(2, 1) = c(2) - a(2);
    }

    Triangle(Vector3f _a, Vector3f _b, Vector3f _c)
        : Triangle(_a, _b, _c, (_b - _a).cross(_c - _a).normalized(), (_b - _a).cross(_c - _a).normalized(), (_b - _a).cross(_c - _a).normalized())
    {
    }

    void computeWorldAABB()
    {
        // Compute the three corners of the triangle in world coordinates.
        Vector3f corners[3];
        corners[0] = transform * a;
        corners[1] = transform * b;
        corners[2] = transform * c;
        // Compute a new bounding box around the transformed corners.
        float xmin = FLT_MAX, xmax = -FLT_MAX, ymin = FLT_MAX, ymax = -FLT_MAX, zmin = FLT_MAX, zmax = -FLT_MAX;
        for (int i = 0; i < 3; i++) {
            float x = corners[i](0),
                y = corners[i](1),
                z = corners[i](2);
            if (x < xmin)
                xmin = x;
            if (x > xmax)
                xmax = x;
            if (y < ymin)
                ymin = y;
            if (y > ymax)
                ymax = y;
            if (z < zmin)
                zmin = z;
            if (z > zmax)
                zmax = z;
        }
        if (xmin == FLT_MAX)
            xmin = xmax - 0.025;
        if (xmax == -FLT_MAX)
            xmax = xmin + 0.025;
        if (ymin == FLT_MAX)
            ymin = ymax - 0.025;
        if (ymax == -FLT_MAX)
            ymax = ymin + 0.025;
        if (zmin == FLT_MAX)
            zmin = zmax - 0.025;
        if (zmax == -FLT_MAX)
            zmax = zmin + 0.025;
        worldAABB = AABB(xmin, xmax, ymin, ymax, zmin, zmax);
    }

    bool intersects(const Ray& originalViewRay, float maxTime = FLT_MAX) const
    {
        Ray viewRayInObjectSpace = rayInObjectSpace(originalViewRay);
        Matrix3f intersectionMatrix = cachedIntersectionMatrix;
        intersectionMatrix(0, 2) = -viewRayInObjectSpace.direction(0);
        intersectionMatrix(1, 2) = -viewRayInObjectSpace.direction(1);
        intersectionMatrix(2, 2) = -viewRayInObjectSpace.direction(2);
        Vector3f intersectionVector(viewRayInObjectSpace.start(0) - a(0), viewRayInObjectSpace.start(1) - a(1), viewRayInObjectSpace.start(2) - a(2));
        Vector3f solution = intersectionMatrix.inverse() * intersectionVector;
        float beta = solution(0), gamma = solution(1), t = solution(2);
        return t > EPSILON && t < maxTime && beta > 0 && gamma > 0 && (beta + gamma) < 1;
    }

    // Returns the normal at a given intersection point given by beta and gamma.
    Vector3f interpolatedNormal(float beta, float gamma) const
    {
        return (normalA + (beta * (normalB - normalA)) + (gamma * (normalC - normalA))).normalized();
    }

    // Arguments intersection, normal and bounce are set in terms of world coordinates.
    bool intersects(const Ray& originalViewRay, Vector3f& intersection, Ray& normal, Ray& bounce) const
    {
        Ray viewRayInObjectSpace = rayInObjectSpace(originalViewRay);
        Matrix3f intersectionMatrix = cachedIntersectionMatrix;
        intersectionMatrix(0, 2) = -viewRayInObjectSpace.direction(0);
        intersectionMatrix(1, 2) = -viewRayInObjectSpace.direction(1);
        intersectionMatrix(2, 2) = -viewRayInObjectSpace.direction(2);
        Vector3f intersectionVector(viewRayInObjectSpace.start(0) - a(0), viewRayInObjectSpace.start(1) - a(1), viewRayInObjectSpace.start(2) - a(2));
        Vector3f solution = intersectionMatrix.inverse() * intersectionVector;
        float beta = solution(0), gamma = solution(1), t = solution(2);
        if (beta < 0 || gamma < 0 || (beta + gamma) > 1 || t <= EPSILON)
            return false;

        Vector3f intersectionInObjectSpace = viewRayInObjectSpace.positionAtTime(t);
        Vector3f incomingIncidenceVector = viewRayInObjectSpace.direction.normalized();
        Vector3f intersectionNormal = interpolatedNormal(beta, gamma);
        Vector3f outgoingReflectionVector = (-2 * incomingIncidenceVector.dot(intersectionNormal)) * intersectionNormal + incomingIncidenceVector;
        intersection = transform * intersectionInObjectSpace;
        Vector3f bounceOffset = transform * (intersectionInObjectSpace + outgoingReflectionVector);
        Vector3f normalOffset = transform * (intersectionInObjectSpace + intersectionNormal);
        bounce.start = intersection;
        bounce.direction = bounceOffset - intersection;
        normal.start = intersection;
        normal.direction = normalOffset - intersection;
        return true;
    }

    string toString() const
    {
        ostringstream out;
        out << "Triangle(a=[" << setprecision(3) << a(0) << ", " << setprecision(3) << a(1) << ", " << setprecision(3) << a(2) << "], "
            << "b=[" << setprecision(3) << b(0) << ", " << setprecision(3) << b(1) << ", " << setprecision(3) << b(2) << "], "
            << "c=[" << setprecision(3) << c(0) << ", " << setprecision(3) << c(1) << ", " << setprecision(3) << c(2) << "])";
        return out.str();
    }

    Vector3f a;
    Vector3f b;
    Vector3f c;
    Vector3f normalA;
    Vector3f normalB;
    Vector3f normalC;
    Matrix3f cachedIntersectionMatrix;
};


AABB computeAABBFromObjects(const vector<Object*>& objects)
{
    float xmin = FLT_MAX, xmax = -FLT_MAX, ymin = FLT_MAX, ymax = -FLT_MAX, zmin = FLT_MAX, zmax = -FLT_MAX;
    for (Object* object : objects) {
        if (object->worldAABB.xmin < xmin)
            xmin = object->worldAABB.xmin;
        if (object->worldAABB.xmax > xmax)
            xmax = object->worldAABB.xmax;
        if (object->worldAABB.ymin < ymin)
            ymin = object->worldAABB.ymin;
        if (object->worldAABB.ymax > ymax)
            ymax = object->worldAABB.ymax;
        if (object->worldAABB.zmin < zmin)
            zmin = object->worldAABB.zmin;
        if (object->worldAABB.zmax > zmax)
            zmax = object->worldAABB.zmax;
    }
    return AABB(xmin, xmax, ymin, ymax, zmin, zmax);
}

class AABBNode {
public:
    AABBNode(const vector<Object*> _objects, int _depth)
        : objects(_objects)
        , depth(_depth)
        , pos(NULL)
        , neg(NULL)
    {
        box = computeAABBFromObjects(objects);
    }

    AABBNode(AABB _box, int _depth, AABBNode* _pos = NULL, AABBNode* _neg = NULL)
        : box(_box)
        , depth(_depth)
        , pos(_pos)
        , neg(_neg)
    {
    }

    ~AABBNode()
    {
        if (pos) delete pos;
        if (neg) delete neg;
    }

    bool isLeaf() const
    {
        return !pos && !neg;
    }

    void collectObjectsForRayIntersection(const Ray& ray, vector<Object*>& result, int depth = 0) const
    {
        if (!box.intersects(ray))
            return;

        // From here on out, assume that the ray intersects my AABB.
        if (isLeaf()) {
            for (Object* object : objects)
                result.push_back(object);
        } else {
            pos->collectObjectsForRayIntersection(ray, result, depth + 1);
            neg->collectObjectsForRayIntersection(ray, result, depth + 1);
        }
    }

    void nodeCount(int& numTreeNodes, int& numLeafNodes) const
    {
        if (isLeaf())
            numLeafNodes++;
        else {
            numTreeNodes++;
            pos->nodeCount(numTreeNodes, numLeafNodes);
            neg->nodeCount(numTreeNodes, numLeafNodes);
        }
    }

    string toString() const
    {
        ostringstream out;
        if (isLeaf()) {
            out << "AABBNode(objects=[ ";
            for (Object* object : objects)
                out << object->id << ", ";
            out << "])";
            return out.str();
        }
        out << "AABBNode(box=" << box.toString() << ",\n";
        for (int _ = 0; _ < depth; _++)
            out << "    ";
        out << "    pos=" << pos->toString() << ",\n";
        for (int _ = 0; _ < depth; _++)
            out << "    ";
        out << "    neg=" << neg->toString() << "\n";
        for (int _ = 0; _ < depth; _++)
            out << "    ";
        out << ")";
        return out.str();
    }

    vector<Object*> objects;
    AABB box;
    int depth;
    AABBNode* pos;
    AABBNode* neg;
};

void splitObjectsByMidpointAlongAxis(const vector<Object*>& objects, const Vector3f& midpoint, int axis, vector<Object*>& pos, vector<Object*>& neg)
{
    bool discriminator = objects.size() % 2 == 0;
    for (Object* object : objects) {
        if (object->midpoint(axis) < midpoint(axis))
            neg.push_back(object);
        else if (object->midpoint(axis) > midpoint(axis))
            pos.push_back(object);
        else {
            if (discriminator)
                pos.push_back(object);
            else
                neg.push_back(object);
            discriminator = !discriminator;
        }
    }
}

AABBNode* makeAABBNode(const vector<Object*>& objects, int depth = 0)
{
    if (depth > MAX_AABB_DEPTH || objects.size() <= AABB_LEAF_SIZE)
        return new AABBNode(objects, depth);

    AABB rootBoundingBox = computeAABBFromObjects(objects);
    Vector3f midpoint = rootBoundingBox.midpoint();
    vector<Object*> pos, neg;
    splitObjectsByMidpointAlongAxis(objects, midpoint, depth % 3, pos, neg);
    return new AABBNode(rootBoundingBox, depth, makeAABBNode(pos, depth + 1), makeAABBNode(neg, depth + 1));
}

class SkyBox {
public:
    SkyBox()
        : backgroundColor(Color(0, 0, 0))
        , isSolidColorSkyBox(true)
    {
    }

    SkyBox(Color _backgroundColor)
        : backgroundColor(_backgroundColor)
        , isSolidColorSkyBox(true)
    {
    }

    SkyBox(string fileNamePrefix, float _xmin, float _xmax, float _ymin, float _ymax, float _zmin, float _zmax, int _width, int _height)
        : xmin(_xmin)
        , xmax(_xmax)
        , ymin(_ymin)
        , ymax(_ymax)
        , zmin(_zmin)
        , zmax(_zmax)
        , sizeX(_xmax - _xmin)
        , sizeY(_ymax - _ymin)
        , sizeZ(_zmax - _zmin)
        , width(_width)
        , height(_height)
        , backgroundColor(Color(0, 0, 0))
        , isSolidColorSkyBox(false)
    {
        loadImageToVector(faces[0], (fileNamePrefix + "1.png").c_str(), width, height);
        loadImageToVector(faces[1], (fileNamePrefix + "2.png").c_str(), width, height);
        loadImageToVector(faces[2], (fileNamePrefix + "3.png").c_str(), width, height);
        loadImageToVector(faces[3], (fileNamePrefix + "4.png").c_str(), width, height);
        loadImageToVector(faces[4], (fileNamePrefix + "5.png").c_str(), width, height);
        loadImageToVector(faces[5], (fileNamePrefix + "6.png").c_str(), width, height);
    }

    float faceRayIntersectionTime(int face, const Ray& ray) const
    {
        switch (face) {
            case 1:
            return (zmin - ray.start(2)) / ray.direction(2);
            case 2:
            return (xmax- ray.start(0)) / ray.direction(0);
            case 3:
            return (zmax - ray.start(2)) / ray.direction(2);
            case 4:
            return (xmin - ray.start(0)) / ray.direction(0);
            case 5:
            return (ymax - ray.start(1)) / ray.direction(1);
            case 6:
            return (ymin - ray.start(1)) / ray.direction(1);
        }
        return -1;
    }

    Color colorForCoordinateOnFace(int face, float u, float v) const
    {
        return colorFromImageVector(faces[face - 1], width, height, u, v);
    }

    Color colorFromRay(const Ray& ray) const
    {
        if (isSolidColorSkyBox)
            return backgroundColor;

        float intersectionTime = FLT_MAX;
        int intersectionFace;
        for (int face = 1; face <= 6; face++) {
            float tFace = faceRayIntersectionTime(face, ray);
            if (tFace > 0 && tFace < intersectionTime) {
                intersectionTime = tFace;
                intersectionFace = face;
            }
        }
        if (intersectionTime == FLT_MAX)
            return Color(1, 1, 1); // Should never happen.

        Vector3f intersection = ray.positionAtTime(intersectionTime);
        // Compute u and v, which are parameters on the interval [0, 1].
        float u, v;
        float intersectionX = intersection(0),
            intersectionY = intersection(1),
            intersectionZ = intersection(2);
        switch (intersectionFace) {
            case 1:
            u = (intersectionX - xmin) / sizeX;
            v = (ymax - intersectionY) / sizeY;
            break;

            case 2:
            u = (intersectionZ - zmin) / sizeZ;
            v = (ymax - intersectionY) / sizeY;
            break;

            case 3:
            u = (xmax - intersectionX) / sizeX;
            v = (ymax - intersectionY) / sizeY;
            break;

            case 4:
            u = (zmax - intersectionZ) / sizeZ;
            v = (ymax - intersectionY) / sizeY;
            break;

            case 5:
            u = (intersectionX - xmin) / sizeX;
            v = (zmax - intersectionZ) / sizeZ;
            break;

            case 6:
            u = (intersectionX - xmin) / sizeX;
            v = (intersectionZ - zmin) / sizeZ;
            break;
        }
        clamp<float>(u, 0, 1);
        clamp<float>(v, 0, 1);
        return colorForCoordinateOnFace(intersectionFace, u, v);
    }

    float xmin;
    float xmax;
    float ymin;
    float ymax;
    float zmin;
    float zmax;
    float sizeX;
    float sizeY;
    float sizeZ;
    int width;
    int height;
    vector<unsigned char> faces[6];

    Color backgroundColor;
    bool isSolidColorSkyBox;
};

#endif /* __UTIL_H__ */
