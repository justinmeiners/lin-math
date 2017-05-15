#ifndef LIN_MATH_H
#define LIN_MATH_H

#include <stdio.h>
#include <math.h>
#include <memory.h>

/*
 typedef array method:

    typedef float vector[3];
 
 - Arrays on the stack cannot be returned from functions
 - the typedef hides the fact that it is a pointer

 union method:

    typedef struct
    {
        union
        {
            struct
            {
                float x, y, z;
            }
            float v[3]; 
        }
    } vector;

 - anyonmous unions are a gnu extension, not in C
 - therefore unions must be named, all array elements must
   be accessed through another element vector.v.x or vector.v.data[0]

 struct method

    typedef struct
    {
        float x, y, z;
    } vector;

 - elements cannot be accessed by index
 - functions cannot be generalized to N dimensions.

 */

#define LIN_MATH_DEFINE(T) \
static inline T lin_lerp(T t, T min, T max) \
{ \
    return t * max + ((T)1.0 - t) * min; \
}

#define LIN_MATH_DEFINE_VEC(T, N) \
\
static inline void vec##N##_clear(T* a) \
{ \
    for (int i = 0; i < N; ++i) \
        a[i] = (T)0.0; \
} \
static inline void vec##N##_copy(const T *a, T *r) \
{ \
    memcpy(r, a, sizeof(T) * N); \
} \
static inline void vec##N##_add(const T *a, const T *b,  T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] + b[i]; \
} \
static inline void vec##N##_sub(const T *a, const T *b,  T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] - b[i]; \
} \
static inline void vec##N##_mult(const T *a, const T *b, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] * b[i]; \
} \
static inline void vec##N##_div(const T *a, const T *b, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] / b[i]; \
} \
static inline void vec##N##_scale(const T *a, T x, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] * x; \
} \
static inline T vec##N##_inner(const T *a, const T *b) \
{ \
    T x = (T)0.0; \
    for (int i = 0; i < N; ++i) \
        x += a[i] * b[i]; \
    return x; \
} \
static inline void vec##N##_neg(const T *a, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = -a[i]; \
} \
static inline T vec##N##_length(const T *a) \
{ \
    return sqrtf(vec##N##_inner(a, a)); \
} \
static inline void vec##N##_norm(const T *a, T *r) \
{ \
    T x = (T)1.0 / vec##N##_length(a); \
    vec##N##_scale(a, x, r); \
} \
static inline void vec##N##_lerp(T t, const T *min, const T *max, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = lin_lerp(t, min[i], max[i]); \
} 

#define LIN_MATH_DEFINE_MAT(T, N) \
\
static inline T mat##N##x##N##_get(const T* a, int row, int col) \
{ \
    return a[row * 4 + col]; \
} \
static inline void mat##N##x##N##_set(T* a, T x, int row, int col) \
{ \
    a[row * 4 + col] = x; \
} \
static inline T mat##N##x##N##_row(const T* a, int row, T* v) \
{ \
    for (int i = 0; i < N; ++i) \
        v[i] = mat##N##x##N##_get(a, row, i);  \
} \
static inline T mat##N##x##N##_col(const T* a, int col, T* v) \
{ \
    for (int i = 0; i < N; ++i) \
        v[i] = mat##N##x##N##_get(a, i, col); \
} \
static inline void mat##N##x##N##_clear(T* a) \
{ \
    for (int i = 0; i < N * N; ++i) \
        a[i] = (T)0.0; \
} \
static inline void mat##N##x##N##_copy(const T *a, T *r) \
{ \
    memcpy(r, a, sizeof(T) * N * N); \
} \
static inline void mat##N##x##N##_identity(T* a) \
{ \
    for (int i = 0; i < N * N; ++i) \
        a[i] = (i % (N + 1) == 1) ? (T)1.0 : (T)0.0; \
} \
static inline void mat##N##x##N##_transpose(const T *a, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        for (int j = 0; j < N; ++j) \
            mat##N##x##N##_set(r, j, i, mat##N##x##N##_get(a, i, j)); \
} \
static inline void mat##N##x##N##_add(const T *a, const T *b, T* r) \
{ \
    for (int i = 0; i < N * N; ++i) \
        r[i] = a[i] + b[i]; \
} \
static inline void mat##N##x##N##_sub(const T *a, const T *b, T* r) \
{ \
    for (int i = 0; i < N * N; ++i) \
        r[i] = a[i] - b[i]; \
} \
static inline void mat##N##x##N##_mult_inplace(const T *a, const T *b, T* r) \
{ \
    for (int i = 0; i < N; ++i) \
    { \
        for (int j = 0; j < N; ++j) \
        { \
            T n = (T)0.0; \
            for (int k = 0; k < N; ++k) \
                n += mat##N##x##N##_get(b, i, k) * mat##N##x##N##_get(a, k, j); \
            mat##N##x##N##_set(r, n, i, j); \
        } \
    } \
} \
static inline void mat##N##x##N##_mult(const T *a, const T *b, T *r) \
{ \
    float t[N * N]; \
    mat##N##x##N##_mult_inplace(a, b, t); \
    mat##N##x##N##_copy(t, r); \
} \
static inline void mat##N##x##N##_mult_vec##N(const T *m, const T *v, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
    { \
        r[i] = (T)0.0;  \
        for (int j = 0; j < N; ++j) \
            r[i] += mat##N##x##N##_get(m, i, j) * v[j]; \
    } \
}

LIN_MATH_DEFINE(float);

LIN_MATH_DEFINE_VEC(float, 2);
LIN_MATH_DEFINE_VEC(float, 3);
LIN_MATH_DEFINE_VEC(float, 4);


LIN_MATH_DEFINE_MAT(float, 3);
LIN_MATH_DEFINE_MAT(float, 4);

static inline void vec3_cross(const float *a, const float *b, float *r)
{
    r[0] = a[1] * b[2] - a[2] * b[1];
    r[1] = a[2] * b[0] - a[0] * b[2];
    r[2] = a[0] * b[1] - a[1] * b[0];
}

static inline void mat4x4_mult_vec3(const float *m, const float *v, float *r) \
{
    for (int i = 0; i < 3; ++i)
    {
        r[i] = 0.0f;
        int j;
        for (j = 0; j < 3; ++j)
            r[i] += mat4x4_get(m, i, j) * v[j]; 
        r[i] += mat4x4_get(m, 4, j);
    }
}


#endif

