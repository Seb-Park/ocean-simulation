#version 410 core

layout(location = 0) in float n;

out float work_group;

void main()
{
    // pass into geometry shader
    work_group = n;
}