#version 410 core

// follow steps below to update the number of vertices

layout(binding = 0) uniform FreqBufferObject {
    float freq[65536];
} freqBuffer;

layout(vertices = 256) out;

void main()
{

}