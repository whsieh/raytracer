
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

#include "lodepng.cpp"

using Eigen::Vector3d;
using namespace std;

#define PI 3.14159265

enum LightSourceType { DIRECTIONAL, POINT, AMBIENT, UNSPECIFIED };

enum TransformationType { ROTATION_ONLY, TRANSLATION_ONLY, SCALE_ONLY, COMBINATION };

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
        clamp<float>(r, 0, 1);
        clamp<float>(g, 0, 1);
        clamp<float>(b, 0, 1);
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
    Ray(Vector3d _start, Vector3d _direction, float _tMin = 0, float _tMax = FLT_MAX)
        : start(_start(0), _start(1), _start(2), 1)
        , direction(_direction(0), _direction(1), _direction(2), 0)
        , tMax(_tMax)
        , tMin(_tMin)
    {
    }

    string toString() const
    {
        return "Ray(start=<" + to_string(start(0)) + ","
            + to_string(start(1)) + ","
            + to_string(start(2)) + ">, direction=<"
            + to_string(direction(0)) + ","
            + to_string(direction(1)) + ","
            + to_string(direction(2)) + ">)";
    }

    Vector4d start;
    Vector4d direction;
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

    Camera(Vector3d _eye, Vector3d _ll, Vector3d _lr, Vector3d _ul, Vector3d _ur)
        : eye(_eye)
        , ll(_ll)
        , lr(_lr)
        , ul(_ul)
        , ur(_ur)
    {
    }

    Vector3d eye;
    Vector3d ll;
    Vector3d lr;
    Vector3d ul;
    Vector3d ur;
};

class Light
{
public:
    Light(Vector3d vector, Color color)
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
    Vector3d m_vector;
    LightSourceType m_type;

private:
    Color m_color;
};

class DirectionalLight : public Light
{
public:
    DirectionalLight(Vector3d vector, Color color)
        : Light(vector, color)
    {
        m_type = DIRECTIONAL;
        m_vector = m_vector / sqrt(m_vector.dot(m_vector));
    }

    const Vector3d direction() const
    {
        return m_vector;
    }
};

class PointLight : public Light
{
public:
    PointLight(Vector3d vector, Color color)
        : Light(vector, color)
    {
        m_type = POINT;
    }

    const Vector3d position() const
    {
        return m_vector;
    }
};

class AmbientLight : public Light
{
public:
    AmbientLight(Color color)
        : Light(Vector3d(0, 0, 0), color)
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

class Rotation;
class Scale;
class Translation;

class Transformation {
public:
    Transformation()
    {
        reset();
        needsUpdate = false;
    }

    virtual Matrix4d matrix() const
    {
        return tMatrix;
    }

    void push(Transformation& transformation)
    {
        tMatrix = transformation.matrix() * tMatrix;
        needsUpdate = true;
    }

    void reset()
    {
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                tMatrix(i, j) = i == j;
    }

    Vector4d apply(Vector4d v)
    {
        return matrix() * v;
    }

    Vector4d apply(Vector3d v, bool indicatesPosition = true)
    {
        return matrix() * (indicatesPosition ? Vector4d(v(0), v(1), v(2), 1) : Vector4d(v(0), v(1), v(2), 0));
    }

    Matrix4d tMatrix;
    bool needsUpdate;
};

class Rotation : public Transformation {
public:
    Rotation(float degX, float degY, float degZ)
        : axis(Vector3d(degX, degY, degZ))
    {
        needsUpdate = true;
        magnitude = sqrt(square(degX) + square(degY) + square(degZ));
        axis = axis / magnitude;
        magnitude = PI / 180;
        for (int i = 0; i < 4; i++)
            for (int j = 0; j < 4; j++)
                tMatrix(i, j) = i == j;
    }

    void push(Transformation& transformation) { }

    void reset() {}

    Vector3d axis;
    float magnitude; // In radians.
};

class Scale : public Transformation {
public:
    Scale(float _x, float _y, float _z)
        : x(_x)
        , y(_y)
        , z(_z)
    {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i != j)
                    tMatrix(i, j) = 0;
                else if (i == 0)
                    tMatrix(i, j) = x;
                else if (i == 1)
                    tMatrix(i, j) = y;
                else if (i == 2)
                    tMatrix(i, j) = z;
                else
                    tMatrix(i, j) = 1;
            }
        }
    }

    void push(Transformation& transformation) { }

    void reset() {}

    float x;
    float y;
    float z;
};

class Translation : public Transformation {
public:
    Translation(float _x, float _y, float _z)
        : x(_x)
        , y(_y)
        , z(_z)
    {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (i == j)
                    tMatrix(i, j) = 1;
                else if (i == 0 && j == 3)
                    tMatrix(i, j) = x;
                else if (i == 1 && j == 3)
                    tMatrix(i, j) = y;
                else if (i == 2 && j == 3)
                    tMatrix(i, j) = z;
                else
                    tMatrix(i, j) = 0;
            }
        }
    }

    void push(Transformation& transformation) { }

    void reset() {}

    float x;
    float y;
    float z;
};

class Object {
public:
    Object()
    {
    }

    virtual Ray transform(const Ray& ray) const = 0;
    virtual string toString() const = 0;
    virtual bool intersects(const Ray&, Ray&) const = 0;

    Transformation transformation;
    Material material;
};

class Sphere : public Object {
public:
    Sphere(Vector3d _center, float _radius)
        : center(_center)
        , radius(_radius)
    {
    }

    Ray transform(const Ray& ray) const
    {
        return Ray(transformation.matrix().transpose() * ray.start, transformation.matrix().transpose() * ray.direction);
        // cout << "ray.start: " << endl << ray.start << endl;
        // cout << "ray.direction: " << endl << ray.direction << endl;
        // cout << "M.T: " << endl << transformation.matrix().transpose() << endl;
        // return Ray(Vector3d(0, 0, 0), Vector3d(0, 0, 0));
    }

    bool intersects(const Ray& ray, Ray& bounce) const
    {
        Ray viewRay = transform(ray);
        cout << "The ray is: " << ray.toString() << endl;
        float a = viewRay.direction(0),
            b = viewRay.direction(1),
            c = viewRay.direction(2),
            d = viewRay.start(0) - center(0),
            e = viewRay.start(1) - center(1),
            f = viewRay.start(2) - center(2);
        float A = square(a) + square(b) + square(c);
        float B = 2 * (a * d + b * e + c * f);
        float C = square(d) + square(e) + square(f) - square(radius);
        cout << "The discriminant is: " << square(B) - (4 * A * C) << endl;
        return false;
    }

    string toString() const
    {
        string s = "Circle(";
        s += "center=[" + to_string(center(0)) + ", "
            + to_string(center(1)) + ", "
            + to_string(center(2)) + "], radius=" + to_string(radius);
        s += ")";
        return s;
    }

    Vector3d center;
    float radius;
};

class Triangle : public Object {
public:
    Triangle(Vector3d _a, Vector3d _b, Vector3d _c)
        : a(_a)
        , b(_b)
        , c(_c)
    {
    }

    Vector3d a;
    Vector3d b;
    Vector3d c;
};

#endif /* __UTIL_H__ */
