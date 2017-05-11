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

#define LIN_MATH(T) \
static inline T lin_lerp(T t, T min, T max) \
{ \
    return t * max + ((T)1.0 - t) * min; \
}

#define LIN_MATH_VEC(N, T) \
\
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
static inline T vec##N##_length(const T *a) \
{ \
    return sqrtf(vec##N##_inner(a, a)); \
} \
static inline void vec##N##_norm(const T *a,  T*r) \
{ \
    T x = (T)1.0 / vec##N##_length(a); \
    vec##N##_scale(a, x, r); \
} \
static inline void vec##N##_lerp(T t, const T *min, const T *max, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = lin_lerp(t, min[i], max[i]); \
} 

LIN_MATH(float);
LIN_MATH_VEC(2, float);
LIN_MATH_VEC(3, float);
LIN_MATH_VEC(4, float);


#endif

