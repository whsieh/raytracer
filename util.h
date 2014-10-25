
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

#include "lodepng.cpp"

using Eigen::Affine;
using Eigen::Vector3f;
using Eigen::Matrix3f;
using Eigen::Quaternionf;
using Eigen::Transform;
using Eigen::Translation;
using Eigen::AngleAxisf;
using namespace std;

// THIS IS HOW YOU DO THE ROTATON THING
// Transform<float, 3, Affine> testTransform = Eigen::Transform<float, 3, Affine>(AngleAxisf(PI / 2, Vector3f(1, 0, 0)).toRotationMatrix());

#define PI 3.14159265
#define EPSILON 0.01

static int GLOBAL_OBJECT_COUNT = 0;

enum LightSourceType { DIRECTIONAL, POINT, AMBIENT, UNSPECIFIED };

enum TransformationType { ROTATION_ONLY, TRANSLATION_ONLY, SCALE_ONLY, COMBINATION };

unsigned int randint(unsigned int lower, unsigned int upper)
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

class Ray {
public:
    Ray(Vector3f _start, Vector3f _direction, float _tMin = 0, float _tMax = FLT_MAX)
        : start(_start(0), _start(1), _start(2))
        , direction(_direction(0), _direction(1), _direction(2))
        , tMax(_tMax)
        , tMin(_tMin)
    {
    }

    Ray()
        : start(0, 0, 0)
        , direction(0, 0, 0)
        , tMax(FLT_MAX)
        , tMin(0)
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

    Vector3f positionAtTime(float t)
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
    float tMax;
    float tMin;
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
    }

    Ray rayInObjectSpace(const Ray& viewRay) const
    {
        Transform<float, 3, Affine> inverseTransform = transform.inverse();
        Vector3f transformedStart = inverseTransform * viewRay.start;
        Vector3f transformedOffsetPosition = inverseTransform * (viewRay.start + viewRay.direction);
        return Ray(transformedStart, transformedOffsetPosition - transformedStart);
    }

    virtual string toString() const = 0;
    virtual bool intersects(const Ray&) const = 0;
    virtual bool intersects(const Ray&, Vector3f&, Ray&, Ray&) const = 0;

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

    bool intersects(const Ray& originalViewRay) const
    {
        Ray viewRayInObjectSpace = rayInObjectSpace(originalViewRay);
        // cout << "The view ray in object coordinates is " << viewRayInObjectSpace.toString() << endl;
        float A = square(viewRayInObjectSpace.direction(0)) + square(viewRayInObjectSpace.direction(1)) + square(viewRayInObjectSpace.direction(2));
        float B = 2 * ((viewRayInObjectSpace.direction(0) * (viewRayInObjectSpace.start(0) - center(0))) +
                    (viewRayInObjectSpace.direction(1) * (viewRayInObjectSpace.start(1) - center(1))) +
                    (viewRayInObjectSpace.direction(2) * (viewRayInObjectSpace.start(2) - center(2))));
        float C = square(viewRayInObjectSpace.start(0) - center(0)) +
            square(viewRayInObjectSpace.start(1) - center(1)) +
            square(viewRayInObjectSpace.start(2) - center(2)) - square(radius);
        // cout << "A: " << A << " B: " << B << "C: " << C << endl;
        float discriminant = square(B) - (4 * A * C);
        float t1 = -1, t2 = -1;
        if (discriminant >= 0) {
            float rootTerm = sqrt(square(B) - (4 * A * C));
            t1 = (-B - rootTerm) / (2 * A);
            t2 = (-B + rootTerm) / (2 * A);
        }
        return discriminant >= 0 && (t1 > EPSILON || t2 > EPSILON);
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
        out << "Circle(center=[" << setprecision(3) << center(0) << ", " << setprecision(3) << center(1) << ", "
            << setprecision(3) << center(2) << "], radius=" << setprecision(3) << radius << ")";
        return out.str();
    }

    Vector3f center;
    float radius;
};

class Triangle : public Object {
public:
    Triangle(Vector3f _a, Vector3f _b, Vector3f _c)
        : Object(GLOBAL_OBJECT_COUNT++)
        , a(_a)
        , b(_b)
        , c(_c)
    {
        triangleNormal = (b - a).cross(c - a).normalized();
        cachedIntersectionMatrix(0, 0) = b(0) - a(0);
        cachedIntersectionMatrix(1, 0) = b(1) - a(1);
        cachedIntersectionMatrix(2, 0) = b(2) - a(2);
        cachedIntersectionMatrix(0, 1) = c(0) - a(0);
        cachedIntersectionMatrix(1, 1) = c(1) - a(1);
        cachedIntersectionMatrix(2, 1) = c(2) - a(2);
    }

    bool intersects(const Ray& originalViewRay) const
    {
        Ray viewRayInObjectSpace = rayInObjectSpace(originalViewRay);
        Matrix3f intersectionMatrix = cachedIntersectionMatrix;
        intersectionMatrix(0, 2) = -viewRayInObjectSpace.direction(0);
        intersectionMatrix(1, 2) = -viewRayInObjectSpace.direction(1);
        intersectionMatrix(2, 2) = -viewRayInObjectSpace.direction(2);
        Vector3f intersectionVector(viewRayInObjectSpace.start(0) - a(0), viewRayInObjectSpace.start(1) - a(1), viewRayInObjectSpace.start(2) - a(2));
        Vector3f solution = intersectionMatrix.inverse() * intersectionVector;
        return solution(2) > EPSILON && solution(0) > 0 && solution(1) > 0 && (solution(0) + solution(1)) < 1;
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
        Vector3f outgoingReflectionVector = (-2 * incomingIncidenceVector.dot(triangleNormal)) * triangleNormal + incomingIncidenceVector;
        intersection = transform * intersectionInObjectSpace;
        Vector3f bounceOffset = transform * (intersectionInObjectSpace + outgoingReflectionVector);
        Vector3f normalOffset = transform * (intersectionInObjectSpace + triangleNormal);
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
    Vector3f triangleNormal;
    Matrix3f cachedIntersectionMatrix;
};

#endif /* __UTIL_H__ */
