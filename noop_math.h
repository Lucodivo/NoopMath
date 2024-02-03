#pragma once

#include <math.h>
#include <stdint.h>

// Note: If degrees/radians is not specified for angle, angle is assumed to be in radians

#define Min(x, y) (x < y ? x : y)
#define Max(x, y) (x > y ? x : y)

#define Pi32 3.14159265359f
#define PiOverTwo32 1.57079632679f
#define SqrtTwoOverTwo32 0.70710678119f
#define Tau32 6.28318530717958647692f
#define RadiansPerDegree (Pi32 / 180.0f)

#define COMPARISON_EPSILON 0.001f

namespace noop {
  struct vec2 {
    float values[2];
    inline float operator[](uint32_t i) const { return values[i]; }
    inline float& operator[](uint32_t i) { return values[i]; }
  };

  struct vec3 {
    union {
      float values[3];
      vec2 xy;
    };
    inline float operator[](uint32_t i) const { return values[i]; }
    inline float& operator[](uint32_t i) { return values[i]; }
  };

  struct vec4 {
    union {
      float values[4];
      vec3 xyz;
    };
    inline float operator[](uint32_t i) const { return values[i]; }
    inline float& operator[](uint32_t i) { return values[i]; }
  };

// NOTE: column-major
  struct mat2 {
    float values[4];
    inline float operator[](uint32_t i) const { return values[i]; }
    inline float& operator[](uint32_t i) { return values[i]; }
  };

// NOTE: column-major
  struct mat3 {
    float values[9];
    inline float operator[](uint32_t i) const { return values[i]; }
    inline float& operator[](uint32_t i) { return values[i]; }
    inline mat2 toMat2(){ return mat2{
      values[0], values[1],
      values[3], values[4]
    };}
  };

// NOTE: column-major
  struct mat4 {
    float values[16];
    inline float operator[](uint32_t i) const { return values[i]; }
    inline float& operator[](uint32_t i) { return values[i]; }
    inline mat3 toMat3(){ return mat3{
          values[0], values[1], values[2],
          values[4], values[5], values[6],
          values[8], values[9], values[10]
      };
    }
    inline static mat4 fromMat3(const mat3& M) {
      return mat4{
        M[0], M[1], M[2], 0.0f,
        M[3], M[4], M[5], 0.0f,
        M[6], M[7], M[8], 0.0f,
        0.0f, 0.0f, 0.0f, 1.0f,
      };
    }
  };

  struct quaternion {
    float r;
    union {
      struct {
        float i; float j; float k;
      };
      vec3 ijk;
    };
  };

  struct complex {
    float r;
    float i;
  };

  struct BoundingRect {
    vec2 min;
    vec2 diagonal;
  };

  struct BoundingBox {
    vec3 min;
    vec3 diagonal;
  };

// floating point
  bool epsilonComparison(float a, float b, float epsilon = COMPARISON_EPSILON);

  bool epsilonComparison(double a, double b, double epsilon = COMPARISON_EPSILON);

  float step(float edge, float x);

  float clamp(float minVal, float maxVal, float x);

  float smoothStep(float edge1, float edge2, float x);

  float lerp(float a, float b, float t);

  float sign(float x);

  float radians(float degrees);

// real-time rendering 4.7.2
// ex: screenWidth = 20.0f, screenDist = 30.0f will provide the horizontal field of view
// for a person sitting 30 units away from a 20 unit screen, assuming the screen is
// perpendicular to the line of sight.
// NOTE: Any units work as long as they are the same. Works for vertical and horizontal.
  float fieldOfView(float screenWidth, float screenDist);

// == vec2 ==
  bool operator==(const vec2 &v1, const vec2 &v2);

  float dot(vec2 xy1, vec2 xy2);

  float magnitudeSquared(vec2 xy);

  float magnitude(vec2 xy);

  vec2 normalize(const vec2 &xy);

  vec2 normalize(float x, float y);

// v2 exists in the same half circle that centers around v1. Acute (<90º) angle between the two vectors.
  bool similarDirection(vec2 v1, vec2 v2);

  vec2 operator-(const vec2 &xy);

  vec2 operator-(const vec2 &xy1, const vec2 &xy2);

  vec2 operator+(const vec2 &xy1, const vec2 &xy2);

  void operator-=(vec2 &xy1, const vec2 &xy2);

  void operator+=(vec2 &xy1, const vec2 &xy2);

  vec2 operator*(float s, vec2 xy);

  vec2 operator*(vec2 xy, float s);

  void operator*=(vec2 &xy, float s);

  vec2 operator/(const vec2 &xy, const float s);

  vec2 operator/(const vec2 &xy1, const vec2 &xy2);

  vec2 lerp(const vec2 &a, const vec2 &b, float t);

  bool lineSegmentsIntersection(vec2 A1, vec2 A2, vec2 B1, vec2 B2, vec2* intersection = nullptr);

// == vec3 ==
  vec3 Vec3(vec2 xy, float z);

  vec3 Vec3(float value);

  bool operator==(const vec3 &v1, const vec3 &v2);

  vec3 operator-(const vec3 &xyz);

  vec3 operator-(const vec3 &xyz1, const vec3 &xyz2);

  vec3 operator+(const vec3 &xyz1, const vec3 &xyz2);

  void operator-=(vec3 &xyz1, const vec3 &xyz2);

  void operator+=(vec3 &xyz1, const vec3 &xyz2);

  vec3 operator*(float s, vec3 xyz);

  vec3 operator*(vec3 xyz, float s);

  void operator*=(vec3 &xyz, float s);

  vec3 operator/(const vec3 &xyz, const float s);

  vec3 operator/(const vec3 &xyz1, const vec3 &xyz2);

  float dot(vec3 xyz1, vec3 xyz2);

  vec3 hadamard(const vec3 &xyz1, const vec3 &xyz2);

  vec3 cross(const vec3 &xyz1, const vec3 &xyz2);

  float magnitudeSquared(vec3 xyz);

  float magnitude(vec3 xyz);

  vec3 normalize(const vec3 &xyz);

  vec3 normalize(float x, float y, float z);

  bool degenerate(const vec3 &v);

  vec3 projection(const vec3 &v1, const vec3 &ontoV2 /* assumed to be normalized */);

  vec3 perpendicularTo(const vec3 &v1, const vec3 &ontoV2 /* assumed to be normalized */);

// v2 exists in the same hemisphere that centers around v1. Acute (<90º) angle between the two vectors.
  bool similarDirection(const vec3 &v1, const vec3 &v2);

  vec3 lerp(const vec3 &a, const vec3 &b, float t);

// vec4
  bool operator==(const vec4 &v1, const vec4 &v2);

  vec4 Vec4(vec2 xy, float z,  float w);

  vec4 Vec4(vec3 xyz, float w);

  vec4 Vec4(float x, vec3 yzw);

  vec4 Vec4(vec2 xy, vec2 zw);

  float dot(vec4 xyzw1, vec4 xyzw2);

  float magnitudeSquared(vec4 xyzw);

  float magnitude(vec4 xyzw);

  vec4 normalize(const vec4 &xyzw);

  vec4 normalize(float x, float y, float z, float w);

  vec4 operator-(const vec4 &xyzw);

  vec4 operator-(const vec4 &xyzw1, const vec4 &xyzw2);

  vec4 operator+(const vec4 &xyzw1, const vec4 &xyzw2);

  void operator-=(vec4 &xyzw1, const vec4 &xyzw2);

  void operator+=(vec4 &xyzw1, const vec4 &xyzw2);

  vec4 operator*(float s, vec4 xyzw);

  vec4 operator*(vec4 xyzw, float s);

  void operator*=(vec4 &xyzw, float s);

  vec4 operator/(const vec4 &xyzw, const float s);

  vec4 operator/(const vec4 &xyzw1, const vec4 &xyzw2);

  vec4 lerp(const vec4 &a, const vec4 &b, float t);

  vec4 min(const vec4 &xyzw1, const vec4 &xyzw2);

  vec4 max(const vec4 &xyzw1, const vec4 &xyzw2);

// Complex
// This angle represents a counter-clockwise rotation
  complex Complex(float angle);

  bool operator==(const complex &c1, const complex &c2);

  float magnitudeSquared(complex c);

  float magnitude(complex c);

  vec2 operator*(const complex &ri, vec2 xy /* treated as complex number*/);

  void operator*=(vec2 &xy /* treated as complex number*/, const complex &ri);

// Quaternions
  bool operator==(const quaternion &q1, const quaternion &q2);

  quaternion Quaternion(vec3 v, float angle);

  quaternion identity_quaternion();

  float magnitudeSquared(quaternion rijk);

  float magnitude(quaternion q);

// NOTE: Conjugate is equal to the inverse when the quaternion is unit length
  quaternion conjugate(quaternion q);

  float dot(const quaternion &q1, const quaternion &q2);

  quaternion normalize(const quaternion &q);

  vec3 rotate(const vec3 &v, const quaternion &q);

  quaternion operator*(const quaternion &q1, const quaternion &q2);

  quaternion operator+(const quaternion &q1, const quaternion &q2);

  quaternion operator-(const quaternion &q1, const quaternion &q2);

  quaternion operator*(const quaternion &q1, float s);

  quaternion operator*(float s, const quaternion &q1);

  quaternion operator/(const quaternion &q1, float s);

  quaternion lerp(quaternion a, quaternion b, float t);

  quaternion slerp(quaternion a, quaternion b, float t);

  quaternion orient(const vec3 &startOrientation, const vec3 &endOrientation);

// mat2
  bool operator==(const mat2 &A, const mat2 &B);

  mat2 identity_mat2();

  mat2 scale_mat2(float scale);

  mat2 scale_mat2(vec3 scale);

  mat2 transpose(const mat2 &A);

  mat2 rotate(float radians);

  float determinant(mat2 m);

// mat3
  bool operator==(const mat3 &A, const mat3 &B);

  mat3 identity_mat3();

  mat3 scale_mat3(float scale);

  mat3 scale_mat3(vec3 scale);

  mat3 transpose(const mat3 &A);

  mat3 rotate_mat3(float radians, vec3 v);

  mat3 rotate_mat3(quaternion q);

  // vector treated as column vector
  vec3 operator*(const mat3 &M, const vec3 &v);

  // vector treated as row vector
  vec3 operator*(const vec3 &v, const mat3 &M);

  mat3 operator*(const mat3 &A, const mat3 &B);

// mat4
  mat4 Mat4(mat3 M);

  bool operator==(const mat4 &A, const mat4 &B);

  mat4 identity_mat4();

  mat4 scale_mat4(float scale);

  mat4 scale_mat4(vec3 scale);

  mat4 translate_mat4(vec3 translation);

  mat4 transpose(const mat4 &A);

  mat4 rotate_xyPlane_mat4(float radians);

  mat4 rotate_mat4(float radians, vec3 v);

  mat4 scaleTrans_mat4(const float scale, const vec3 &translation);

  mat4 scaleTrans_mat4(const vec3 &scale, const vec3 &translation);

  mat4 scaleRotTrans_mat4(const vec3 &scale, const vec3 &rotAxis, const float angle,
                          const vec3 &translation);

  mat4 scaleRotTrans_mat4(const vec3 &scale, const quaternion &q, const vec3 &translation);

  mat4 rotate_mat4(quaternion q);

  vec4 operator*(const mat4 &M, const vec4 &v);

  mat4 operator*(const mat4 &A, const mat4 &B);

// real-time rendering 4.7.1
// This projection is for a canonical view volume goes from <-1,1>
  mat4 orthographic(float l, float r, float b, float t, float n, float f);

// real-time rendering 4.7.2
  mat4 perspective(float l, float r, float b, float t, float n, float f);

  mat4 perspectiveInverse(float l, float r, float b, float t, float n, float f);

// real-time rendering 4.7.2
// aspect ratio is equivalent to width / height
  mat4 perspective(float fovVert, float aspect, float n, float f);
  mat4 perspective_fovHorz(float fovHorz, float aspect, float n, float f);

  mat4 perspectiveInverse(float fovVert, float aspect, float n, float f);

/*
 * // NOTE: Oblique View Frustum Depth Projection and Clipping by Eric Lengyel (Terathon Software)
 * Arguments:
 * mat4 perspectiveMat: Perspective matrix that is being used for the rest of the scene
 * vec3 planeNormal_viewSpace: Plane normal in view space. MUST be normalized. This is a normal that points INTO the frustum, NOT one that is generally pointing towards the camera.
 * vec3 planePos_viewSpace: Any position on the plane in view space.
 */
  mat4 obliquePerspective(const mat4 &perspectiveMat, vec3 planeNormal_viewSpace,
                          vec3 planePos_viewSpace, float farPlane);
  mat4 obliquePerspective_fovHorz(float fovHorz, float aspect, float nearPlane, float farPlane, vec3 planeNormal_viewSpace, vec3 planePos_viewSpace);

/*
 * // NOTE: Oblique View Frustum Depth Projection and Clipping by Eric Lengyel (Terathon Software)
 * Arguments:
 * vec3 planeNormal_viewSpace: Plane normal in view space. MUST be normalized. This is a normal that points INTO the frustum, NOT one that is generally pointing towards the camera.
 * vec3 planePos_viewSpace: Any position on the plane in view space.
 */
  mat4 obliquePerspective(float fovVert, float aspect, float nearPlane, float farPlane,
                          vec3 planeNormal_viewSpace, vec3 planePos_viewSpace);

  void adjustAspectPerspProj(mat4 *projectionMatrix, float fovVert, float aspect);

  void adjustNearFarPerspProj(mat4 *projectionMatrix, float n, float f);

// etc
  bool insideRect(BoundingRect boundingRect, vec2 position);

  bool insideBox(BoundingBox boundingBox, vec3 position);

  bool overlap(BoundingRect bbA, BoundingRect bbB);

  bool overlap(BoundingBox bbA, BoundingBox bbB);
}

#undef COMPARISON_EPSILON