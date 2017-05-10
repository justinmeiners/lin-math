#ifndef LIN_MATH_H
#define LIN_MATH_H

#include <stdio.h>
#include <math.h>
#include <memory.h>


#define LIN_MATH(T) \
static inline T float_lerp(T t, T min, T max) \
{ \
    return t * max + ((T)1.0 - t) * min; \
}


#define LIN_MATH_VEC(N, T) \
\
static inline void float##N##_add(const T *a, const T *b,  T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] + b[i]; \
} \
\
static inline void float##N##_sub(const T *a, const T *b,  T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] - b[i]; \
} \
static inline void float##N##_mult(const T *a, const T *b, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] * b[i]; \
} \
static inline void float##N##_div(const T *a, const T *b, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] / b[i]; \
} \
static inline void float##N##_scale(const T *a, float x, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = a[i] * x; \
} \
static inline T float##N##_inner(const T *a, const T *b) \
{ \
    T x = (T)0.0; \
    for (int i = 0; i < N; ++i) \
        x += a[i] * b[i]; \
    return x; \
} \
static inline T float##N##_length(const T *a) \
{ \
    return sqrtf(float##N##_inner(a, a)); \
} \
static inline void float##N##_norm(const T *a, T* r) \
{ \
    T x = (T)1.0 / float##N##_length(a); \
    float##N##_scale(a, x, r); \
} \
static inline void float##N##_lerp(T t, const T *min, const T *max, T *r) \
{ \
    for (int i = 0; i < N; ++i) \
        r[i] = float_lerp(t, min[i], max[i]); \
} 

LIN_MATH(float);
LIN_MATH_VEC(3, float);


#endif

