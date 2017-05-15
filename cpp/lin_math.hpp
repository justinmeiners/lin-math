#include <cmath>

namespace lin_math
{
    template <typename T>
    inline T lerp(T t, T min, T max)
    {
        return t * max + ((T)1.0 - t) * min;
    }

    template <typename T=float, int N=3>
    struct vec
    {
        typedef vec<T, N> vec_type;
        T v[N];

        vec() {}

        vec(T x, T y)
        {
            v[0] = x;
            v[0] = y;
        }

        vec(T x, T y, T z) 
        {
            v[0] = x;
            v[1] = y;
            v[2] = z;
        }

        vec(T x, T y, T z, T w)
        {
            v[0] = x;
            v[1] = y;
            v[2] = z;
            v[3] = w;
        }

        vec(const vec_type& a)
        {
            memcpy(v, a.v, sizeof(T) * N);
        }

        vec_type& operator=(const vec_type& a)
        {
            memcpy(v, a.v, sizeof(T) * N);
            return *this;
        }

        T& operator[](int i) 
        {
            return v[i];
        }
        const T& operator[](int i) const
        {
            return v[i];
        }
        vec_type& operator+=(const vec_type& rhs)
        {
            for (int i = 0; i < N; ++i)
                v[i] += rhs[i];
            return *this;
        }
        friend vec_type& operator+(vec_type lhs, const vec_type& rhs)
        {
            lhs += rhs;
            return lhs;
        }
        vec_type& operator-=(const vec_type& rhs)
        {
            for (int i = 0; i < N; ++i)
                v[i] -= rhs[i];
            return *this;
        }
        friend vec_type& operator-(vec_type lhs, const vec_type& rhs)
        {
            lhs -= rhs;
            return lhs;
        }
        vec_type& operator*=(const vec_type& rhs)
        {
            for (int i = 0; i < N; ++i)
                v[i] *= rhs[i];
            return *this;
        }
        friend vec_type& operator*(vec_type lhs, const vec_type& rhs)
        {
            lhs *= rhs;
            return lhs;
        }
        vec_type& operator/=(const vec_type& rhs)
        {
            for (int i = 0; i < N; ++i)
                v[i] /= rhs[i];
            return *this;
        }
        friend vec_type& operator/(vec_type lhs, const vec_type& rhs)
        {
            lhs /= rhs;
            return lhs;
        }
        vec_type& operator+=(const T& rhs)
        {
            for (int i = 0; i < N; ++i)
                v[i] += rhs;
            return *this;
        }
        friend vec_type& operator+(vec_type lhs, const T& rhs)
        {
            lhs += rhs;
            return lhs;
        }
        vec_type& operator*=(const T& rhs)
        {
            for (int i = 0; i < N; ++i)
                v[i] *= rhs;
            return *this;
        }
        friend vec_type& operator*(vec_type lhs, const T& rhs)
        {
            lhs *= rhs;
            return lhs;
        }
 
        T inner(const vec_type& x) const
        {
            T p = 0.0;
            for (int i = 0; i < N; ++i)
                p += v[i] * x[i];
            return p;
        }
        T length() const
        {
            return std::sqrt(inner(*this)); 
        }
        vec_type norm() const
        {
            return *this / length();
        }
        void normalize()
        {
            *this /= length();
        }   
        friend vec_type lerp(T t, const vec_type& min, const vec_type& max)
        {
            vec_type r;
            for (int i = 0; i < N; ++i)
                r.v[i] = lerp(min[i], max[i]);
            return r;
        }
    };

    template <typename T=float>
    using vec2 = vec<T, 2>;

    template <typename T=float>
    using vec3 = vec<T, 3>;

    template <typename T=float>
    using vec4 = vec<T, 4>;
}



