#version 410 core

layout(points) in;
layout(points, max_vertices = 512) out;

// buffer texture for the positions
uniform samplerBuffer u_tbo_tex;

// stage of this calculation
uniform int stage;

in float work_group;

vec2 muliply_complex(vec2d a, vec2d b)
{
    return vec2(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
}

void main()
{
    //go over the groups
    for (int group = 0; group < 512 / pow(2, stage); group++)
    {
        for (int butterfly = 0; butterfly < pow(2, stage - 1); butterfly++)
        {
            // get the two indices
            int index1 = group * pow(2, stage) + butterfly;
            int index2 = group * pow(2, stage) + pow(2, stage - 1) + butterfly;

            // get the two positions
            vec4 pos1 = texelFetch(u_tbo_tex, index1);
            vec4 pos2 = texelFetch(u_tbo_tex, index2);

            int index = group * pow(2, stage) + butterfly;
            float w = index * 2 * 3.14159265359 / pow(2, stage);
            vec2 twiddle = vec2(cos(w), sin(w));

            vec2 p = pos1 + muliply_complex(twiddle, pos2);
            gl_Position = vec4(p.x, p.y, 0.0, 1.0);
            EmitVertex();

            p = pos1 - muliply_complex(twiddle, pos2);
            gl_Position = vec4(p.x, p.y, 0.0, 1.0);
            EmitVertex();
        }
    }
}
