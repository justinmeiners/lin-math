#include "lin_math.hpp"
#include <iostream>
#include <vector>

using namespace lin_math;
typedef vec3<float> vec3f;

float rand_float()
{
    return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

int main(int argc, const char* argv[])
{
    srand(time(NULL));

    std::cout << "inserting data" << std::endl;

    int sample_count = 10000;

    auto samples = std::vector<vec3f>(sample_count);
    
    for (int i = 0; i < sample_count; ++i)
    {
        samples[i] = vec3f(rand_float(), rand_float(), rand_float());
    }

    std::cout << "adding" << std::endl;

    for (int i = 0; i < sample_count; ++i)
    {
        for (int j = i; j < sample_count; ++j)
        {
            if (rand() % 2 == 1)
            {
                samples[i] += samples[j];
            }
            else
            {
                samples[i] -= samples[j];
            }
        }
    }

    std::cout << "multiplying" << std::endl;

    for (int i = 1; i < sample_count; ++i)
    {
        for (int j = i; j < sample_count; ++j)
        {
            samples[i] *= samples[j];
            samples[i] *= rand_float();
        } 
    }

    std::cout << "length" << std::endl;

    float sum = 0.0f;
    for (int i = 0; i < sample_count; ++i)
    {
        sum += samples[i].length();
    }

    std::cout << "done: " << sum << std::endl; 


    mat3<float> m = mat3<float>::identity();
    std::cout << m;

    quat<float> q = quat<float>::axis_angle(45.0f, 1.0f, 0.0f, 0.0f);
    std::cout << q;
}
