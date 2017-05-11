#include "lin_math.h"

int main(int argc, const char* argv[])
{

    float mvp[16];
    mat4x4_identity(mvp);

    float model[16];
    mat4x4_identity(model);

    mat4x4_mult(mvp, model, mvp);

    float v[] = {0.0f, 0.0f, 1.0f};
    float b[3];
    mat4x4_mult_vec3(mvp, v, b);
}
