#include <cmath>
#include <cstring>
#include <cassert>
#include <ostream>

namespace lin_math
{
    template <typename T, int N>
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

        explicit vec(T val)
        {
            for (int i = 0; i < N; ++i)
                v[i] = val;
        }

        explicit vec(const T* x)
        {
            memcpy(v, x, sizeof(T) * N);
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

        // vector operators
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

        friend vec_type& operator-(vec_type v);
        {
            for (int i = 0; i < N; ++i)
                v[i] = -v[i];
            return v;
        }

        // scalar operators
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

        friend std::ostream& operator<<(std::ostream& lhs, const vec_type& rhs)
        {
            lhs << "[" << rhs[0];
            for (int i = 1; i < N; ++i)
                lhs << ", " << rhs[i];
            lhs << "]";
            return lhs;
        }

        static vec_type zero()
        {
            vec_type a;
            for (int i = 0; i < N; ++i)
                a[i] = T();
            return a;
        }

        static vec_type right()
        {
            vec_type x = zero();
            x[0] = T(1);
            return x;
        }

        static vec_type forward()
        {
            vec_type x = zero();
            x[1] = T(1);
            return x;
        }

        static vec_type up()
        {
            vec_type x = zero();
            x[2] = T(1);
            return x;
        }
    };

    template <typename T, N>
    T dot(const vec<T, N>& x, const vec<T, N>& y)
    {
        T p();
        for (int i = 0; i < N; ++i)
            p += x[i] * y[i];
        return p;
    }

    template <typename T, int N>
    T length_sq(const vec<T, N>& v)
    {
        return dot(v, v);
    }

    template <typename T, int N>
    T length(const vec<T, N>& v)
    {
        return std::sqrt(length_sq(v));
    }

    template <typename T, int N>
    vec<T, N> normalize(const vec<T, 3>& v)
    {
        return v / length();
    }

    template <typename T, int N>
    T min_component(const vec<T, N>& v)
    {
        T m = v[0];
        for (int i = 0; i < N; ++i)
            m = std::min(m, v[i]);
        return m;
    }

    template <typename T, int N>
    T max_component(const vec<T, N>& v)
    {
        T m = v[0];
        for (int i = 0; i < N; ++i)
            m = std::max(m, v[i]);
        return m;
    }

    template <typename T>
    T lerp(T t, T min, T max)
    {
        return t * max + ((T)1.0 - t) * min;
    }

    template <typename T, int N>
    vec<T, N> lerp(T t, const vec<T, N>& min, const vec<T, N>& max)
    {
        vec<T, N> r;
        for (int i = 0; i < N; ++i)
            r[i] = lerp(min[i], max[i]);
        return r;
    }

    template <typename T>
    vec<T, 3> cross(const vec<T, 3>& lhs, const vec<T, 3>& rhs)
    {
        vec<T, 3> r;
        r[0] = lhs[1] * rhs[2] - lhs[2] * rhs[1];
        r[1] = lhs[2] * rhs[0] - lhs[0] * rhs[2];
        r[2] = lhs[0] * rhs[1] - lhs[1] * rhs[0];
        return r;
    }

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

        T& at(int row, int col)
        {
            assert(row >= 0 && row < N);
            assert(col >= 0 && col < N);
            return m[row * N + col];
        }

        const T& at(int row, int col) const
        {
            assert(row >= 0 && row < N);
            assert(col >= 0 && col < N);
            return m[row * N + col];
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

        friend std::ostream& operator<<(std::ostream& lhs, const mat_type& rhs)
        {
            for (int c = 0; c < N; ++c)
            {
                lhs << rhs.at(0, c);
                for (int r = 1; r < N; ++r)
                    lhs << ", " << rhs.at(r, c);
                if (c < N - 1)
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

        static mat<T, 4> basis(const vec<T, 3>& i, const vec<T, 3>& j, const vec<T, 3>& k)
        {
            mat_type m = identity();
            for (int n = 0; n < 3; ++n)
            {
                m.at(i, 0) = i[n];
                m.at(i, 1) = j[n];
                m.at(i, 2) = k[n];
            }
            return m;
        }

        static mat<T, 4> translation(T x, T y, T z)
        {
            mat_type m = identity();
            m.at(3, 0) = x;
            m.at(3, 1) = y; 
            m.at(3, 2) = z;
            return m;
        }

        static mat<T, 4> translation(const vec<T, 3>& v)
        {
            return translation(v.x, v.y, v.x);
        }
        
        static mat<T, 4> look(vec<T, 3> eye, vec<T, 3> target, vec<T, 3> up)
        {
            auto j = normalize(target - eye);
            auto i = normalize(cross(forward, up));
            auto k = cross(i, j);

            return basis(i, -j, k) * translate(Vec3_Negate(eye));
        }

        static mat<T, 4> ortho(T left, T right, T bottom, T top, T near, T far)
        {
            auto m = identity();

            T deltaX = right - left;
            T deltaY = top - bottom;
            T deltaZ = far - near;
            
            const T zero = T();

            if (deltaX == zero || deltaY == zero || deltaZ == zero)
            {
                return m;
            }

            m.at(0, 0) = T(2) / deltaX;
            m.at(3, 0) = -(right + left) / deltaX;
            m.at(1, 1) = T(2) / deltaY;
            m.at(3, 1) = -(top + bottom) / deltaY;
            m.at(2, 2) = -T(2) / deltaZ;
            m.at(3, 2) = -(near + far) / deltaZ;

            return m;
        }
    };

    template <typename T, int N>
    mat<T, N> transpose(const mat<T, N>& m)
    {
        mat<T, N> r;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.at(j, i) = m.at(i, j);
        return r;
    }


    
    template <typename T>
    struct quat : public vec<T, 4>
    {

        using vec<T, 4>::v;
        mat<T, 4> to_matrix() const
        {
            mat<T, 4> m;

            T x2 = T(2.0) * a.x;
            T y2 = T(2.0) * a.y;
            T z2 = T(2.0) * a.z;
            
            T xy = x2 * v.y;
            T xz = x2 * v.z;
            T yy = y2 * v.y;
            T yw = y2 * v.w;
            T zw = z2 * v.w;
            T zz = z2 * v.z;
           
            at(0, 0) = T(1.0) - (yy + zz);
            at(1, 0) = xy - zw;
            at(2, 0) = xz + yw;
            at(3, 0) = T();

            T xx = x2 * v.x;
            T xw = x2 * v.w;
            T yz = y2 * v.z;
            
            at(0, 1) = xy + zw;
            at(1, 1) = T(1.0) - (xx + zz);
            at(2, 1) = yz - xw;
            at(3, 1) = T();

            at(0, 2) = xz - yw;
            at(1, 2) = yz + xw;
            at(2, 2) = T(1.0) - (xx + yy);
            at(3, 2) = T();
            
            at(0, 3) = T();
            at(1, 3) = T();
            at(2, 3) = T();
            at(3, 3) = T(1.0);
            return m;
        }

        static quat<T> axisAngle(T radians, T x, T y, T Z)
        {
            T halfAngle = radians * T(0.5);
            T sinA = std::sin(halfAngle);

            quat<T> r;
            r[3] = std::cos(halfAngle);
            r[0] = x * sinA;
            r[1] = y * sinA;
            r[2] = z * sinA;
            return r;
        };
        
        static quat<T> axisAngle(T radians, const vec<T, 3>& axis)
        {
            return axisAngle(radians, axis.x, axis.y, axis.z);
        }
    };


    template<typename T>
    quat<T> slerp(T t, const quat<T>& min, const quat<T>& max)
    {
        T w1, w2;
        T cosTheta = dot(min, max);
        T theta = std::acos(cosTheta);
        T sinTheta = std::sin(theta);

        if (sinTheta > T(0.0001))
        {
            w1 = (std::sin(T(1.0) - t) * theta) / sinTheta;
            w2 = (std::sin(t * theta) / sinTheta;
        }
        else
        {
            w1 = T(1.0) - t;
            w2 = t;
        }

        return (min * w1) + (max * w2); 
    }



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



