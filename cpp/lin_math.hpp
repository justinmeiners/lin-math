#include <cmath>
#include <cstring>
#include <cassert>
#include <ostream>

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

        // constructors
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

        // conversion
        operator T*() 
        {
            return v;
        }

        operator const T*() const 
        {
            return v;
        }

        T& operator[](int i) 
        {
            assert(i >= 0 && i < N);
            return v[i];
        }

        const T& operator[](int i) const
        {
            assert(i >= 0 && i < N);
            return v[i];
        }

        // vector operations
        vec_type& operator=(const vec_type& a)
        {
            memcpy(v, a.v, sizeof(T) * N);
            return *this;
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

        // scalar operations
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

        vec_type& operator/=(const T& rhs)
        {
            for (int i = 0; i < N; ++i)
                v[i] /= rhs;
            return *this;
        }

        friend vec_type& operator/(vec_type lhs, const T& rhs)
        {
            lhs /= rhs;
            return lhs;
        }

        void clear()
        {
            for (int i = 0; i < N; ++i)
                v[i] = T();
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

        vec_type normalized() const
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

        vec<T, 3> cross(const vec<T, 3>& rhs)
        {
            vec<T, 3> r;
            r[0] = v[1] * rhs[2] - v[2] * rhs[1];
            r[1] = v[2] * rhs[0] - v[0] * rhs[2];
            r[2] = v[0] * rhs[1] - v[1] * rhs[0];
            return r;
        }   
 
        static vec_type zero()
        {
            vec_type a;
            for (int i = 0; i < N; ++i)
                a[i] = T();
            return a;
        }

        static vec_type one()
        {
            vec_type a;
            for (int i = 0; i < N; ++i)
                a[i] = T(1);
            return a;
        }

        friend std::ostream& operator<<(std::ostream& lhs, const vec_type& rhs)
        {
            lhs << "[" << rhs[0];
            for (int i = 1; i < N; ++i)
                lhs << ", " << rhs[i];
            lhs << "]";
            return lhs;
        }
    };

    template <typename T, int N>
    struct mat
    {
        typedef mat<T, N> mat_type;
        typedef vec<T, N> vec_type;

        T m[N * N];

        mat() {}

        mat(const mat_type& rhs)
        {
            memcpy(m, rhs.m, sizeof(T) * N * N);
        }

        mat_type& operator=(const mat& rhs)
        {
            memcpy(m, rhs.m, sizeof(T) * N * N);
            return *this;
        }

        // acccessors
        T& operator[](int i)
        {
            assert(i >= 0 && i < N * N);
            return m[i];
        }

        const T& operator[](int i) const
        {
            assert(i >= 0 && i < N * N);
            return m[i];
        }

        T& at(int row, int column)
        {
            assert(row >= 0 && row < N);
            assert(column >= 0 && column < N);
            return m[row * N + column];
        }

        const T& at(int row, int column) const
        {
            assert(row >= 0 && row < N);
            assert(column >= 0 && column < N);
            return m[row * N + column];
        }

        // conversion
        operator T*() 
        {
            return m;
        }

        operator const T*() const
        {
            return m;
        }

        // vector operators
        friend vec_type operator*(const mat_type& lhs, const vec_type& v)
        {
            vec_type r;
            for (int i = 0; i < N; ++i)
            {
                r[i] = T();
                for (int j = 0; j < N; ++j)
                    r[i] += lhs.at(i, j) * v[j];
            }
            return r;
        }

        // matrix operators
        mat_type& operator+=(const mat_type& rhs)
        {
            for (int i = 0; i < N; ++i)
                m[i] += rhs[i];
            return *this;
        }

        friend mat_type operator+(mat_type lhs, const mat_type& rhs)
        {
            lhs += rhs;
            return lhs;
        }

        mat_type& operator-=(const mat_type& rhs)
        {
            for (int i = 0; i < N; ++i)
                m[i] -= rhs[i];
            return *this;
        }

        friend mat_type operator-(mat_type lhs, const mat_type& rhs)
        {
            lhs -= rhs;
            return lhs;
        }

        vec_type row(int row) const
        {
            assert(row >= 0 && row < N);
            vec_type r;
            for (int i = 0; i < N; ++i)
                r[i] = at(row, i);
            return r;
        }

        vec_type column(int column) const
        {
            assert(column >= 0 && column < N);
            vec_type r;
            for (int i = 0; i < N; ++i)
                r[i] = at(i, column);
            return r;
        }

        void clear()
        {
            for (int i = 0; i < N * N; ++i)
                m[i] = T();
        }

        mat_type transposed() const
        {
            mat_type x;
            for (int r = 0; r < N; ++r)
                for (int c = 0; c < N; ++c)
                    x.at(r, c) = at(c, r);
            return x;
        }

        friend std::ostream& operator<<(std::ostream& lhs, const mat_type& rhs)
        {
            for (int c = 0; c < N; ++c)
            {
                lhs << rhs.at(0, c);
                for (int r = 1; r < N; ++r)
                    lhs << ", " << rhs.at(r, c);
                if (c < N -1)
                    lhs << std::endl;
            }

            return lhs;
        }
 
        static mat_type zero()
        {
            mat_type m;
            for (int i = 0; i < N * N; ++i)
                m[i] = T();
            return m;
        }

        static mat_type identity()
        {
            mat_type m; 
            for (int i = 0; i < N * N; ++i) 
                m[i] = (i % (N + 1) == 0) ? T(1) : T();
            return m;
        }

        static mat_type translation(T x, T y, T z)
        {
            mat_type m = identity();
            m.at(3, 0) = x;
            m.at(3, 1) = y; 
            m.at(3, 2) = z;
            return m;
        }
    };

    template<typename T>
    struct quat : public vec<T, 4>
    {
        typedef quat<T> quat_type;

        static quat_type axis_angle(T angle, T x, T y, T z)
        {
            T halfAngle = angle / T(2);
	        T sinA = std::sin(halfAngle);

            quat_type q;
            q[3] = cosf(halfAngle);
            q[0] = x * sinA;
            q[1] = y * sinA;
            q[2] = z * sinA;
	        return q;
        }
    };

    template <typename T=float>
    using vec2 = vec<T, 2>;

    template <typename T=float>
    using vec3 = vec<T, 3>;

    template <typename T=float>
    using vec4 = vec<T, 4>;

    template <typename T=float>
    using mat3 = mat<T, 3>;

    template <typename T=float>
    using mat4 = mat<T, 4>;
}



