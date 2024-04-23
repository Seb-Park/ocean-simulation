#version 330 core

layout (points) in;
layout (triangle_strip, max_vertices = 36) out;
uniform mat4 view, projection;
uniform vec3 worldSpace_camPos;
uniform float u_time;
uniform sampler2D wind_texture;
uniform sampler2D meadow_texture;
uniform vec3 playerPos;

uniform vec4 plane;



// texture atlasing
//uniform float numRows;
//uniform vec2 atlas_offset;

in VS_OUT {
   float visibility;
} gs_in[];

out GS_OUT {
    vec2 texCoord;
    float visibility;
} gs_out;

// grass shader referenced from: https://vulpinii.github.io/tutorials/grass-modelisation/en/
/////////
float PI = 3.141592653589793;
mat4 windModel = mat4(1.f);
mat4 trampleModel = mat4(1.f);
int numRows = 1;

vec2 atlas_offset = vec2(0.f);

vec2 getOffset(int index, int numRows){
    int column = int(mod(index, numRows));
    float xoffset = float(column)/float(numRows);

    int row = index/numRows;
    float yOffset = float(row)/float(numRows);
    return vec2(xoffset, yOffset);
}

mat4 rotX(float a){
    mat4 rx = mat4(1.f);
    rx[1] = vec4(0.f, cos(a), -sin(a), 0.f);
    rx[2] = vec4(0.f, sin(a), cos(a), 0.f);

    return rx;
}

mat4 rotY(float a){
    mat4 ry = mat4(1.f);
    ry[0] = vec4(cos(a), 0.f, sin(a), 0.f);
    ry[2] = vec4(-sin(a), 0.f, cos(a), 0.f);

    return ry;
}

mat4 rotZ(float a){
    mat4 rz = mat4(1.f);
    rz[0] = vec4(cos(a), sin(a), 0.f, 0.f);
    rz[1] = vec4(-sin(a), cos(a), 0.f, 0.f);

    return rz;
}

float randomize(vec2 st){
    return fract(sin(dot(st.xy,vec2(12.9898,78.233)))*43758.5453123);
}


void makeQuad(vec4 grass_pos, mat4 crossmodel){
    vec4 vertexPos[4];
    vertexPos[0] = vec4(-.9f, 0.f, 0.f, 0.f); // bottom left
    vertexPos[1] = vec4(0.9f, 0.f, 0.f, 0.f); // bottom right
    vertexPos[2] = vec4(-0.9f, 1.f, 0.f, 0.f); // upper left
    vertexPos[3] = vec4(0.9f, 1.f, 0.f, 0.f); // upper right

    vec2 texCoords[4];
    texCoords[0] = vec2(0.f, 0.f);
    texCoords[1] = vec2(1.f, 0.f);
    texCoords[2] = vec2(0.f, 1.f);
    texCoords[3] = vec2(1.f, 1.f);

    mat4 randomY = rotY(randomize(grass_pos.zx)*PI);

    // will apply wind only to top two corners of quad
    mat4 defaultWind = mat4(1.f);
    mat4 defaultTrample = mat4(1.f);

    gl_ClipDistance[0]=dot(vec4(grass_pos), plane);


    vec2 clampedOffset = atlas_offset;
    for (int i=0; i<4; i++){
        if (i>=2) {
            defaultWind = windModel;
            defaultTrample = trampleModel;

        }

        gl_Position = projection*view*(grass_pos + (defaultTrample*defaultWind*randomY*crossmodel*vertexPos[i]));
        gs_out.texCoord = (texCoords[i]/numRows) + clampedOffset;
        EmitVertex();
    }

    EndPrimitive();

}

struct TrampleInfo {
    vec3 trampleOffset;
    float windMultiplier;

};


// trample referenced from NedMakesGames: https://www.youtube.com/watch?v=AmO7k-Lr0XM
TrampleInfo calculateTrample(vec3 entityWorldPos, float maxDistance, float falloff, float pushAwayStrength, float pushDownStrength){
    vec3 offset = vec3(0.f);
    float windMultiplier = 1.f;
    vec3 distanceVec = gl_in[0].gl_Position.xyz-playerPos;
    float distance = length(distanceVec);

    // convert to trample strength
    float strength = 1 - pow(clamp(distance / maxDistance, 0.f, 1.f), falloff);

    // apply pushAway offest in xz plane
    vec3 xzDistance = vec3(distanceVec.x, 0.f, distanceVec.z);
    vec3 pushAwayOffset = normalize(xzDistance) * pushAwayStrength * strength;

    // pushDown offset always points downwards
    vec3 squishOffset = vec3(0.f, -1.f, 0.f) * pushDownStrength * strength;

    offset += pushAwayOffset + squishOffset;

    // supress wind when this grass is being trampled

    windMultiplier = min(windMultiplier, 1-strength);

    TrampleInfo info;
    info.trampleOffset = offset;
    info.windMultiplier = windMultiplier;

    return info;
}

void makeGrass(int numQuads){
    mat4 model0, model45, modelm45;
    model0 = mat4(1.f);
    model45 = rotY(radians(45));
    modelm45 = rotY(-radians(45));


    // wind calculation
    vec2 windDir = vec2(1.f);
    float windStrength = .05f;

    // uv coordinates of wind
    vec2 uv = gl_in[0].gl_Position.xz/10.f + windDir * u_time;
    uv.x = mod(uv.x, 1.0);
    uv.y = mod(uv.y, 1.0);
    vec4 windSpot = texture(wind_texture, uv);

    // get index in meadow texture atlas depending on where uv is

    vec2 meadow_uv = gl_in[0].gl_Position.xz * .03f;
    meadow_uv.x = mod(meadow_uv.x, 1.0);
    meadow_uv.y = mod(meadow_uv.y, 1.0);
    vec4 meadowSpot = texture(meadow_texture, meadow_uv);
    float meadowSpotAccumulate = meadowSpot.r + meadowSpot.g + meadowSpot.b;

    // convert meadowSpot to range 0-maxindex

    numRows = 4;
    int index = int(((meadowSpotAccumulate) / (3.f) ) * (16.f));
    if (index < 0) index = 0;
    if (index > 15) index = 15;
    atlas_offset = getOffset(index, numRows);


    // matrix that rotates top of grass away from player
    TrampleInfo info = calculateTrample(playerPos, 2.f, 5.f, .52f, .51f);
    trampleModel = (rotX(info.trampleOffset.x*PI*1.f) * rotZ(info.trampleOffset.z*PI*1.f));

    // matrix that tilts quad in x and z dir, accordiing to wind dir and force
    windModel = (rotX(windSpot.x*PI*.75f*info.windMultiplier - PI*.25f) * rotZ(windSpot.z*PI*.75f*info.windMultiplier - PI*.25f));


    // draw number of quads based on level of detail
    switch(numQuads){
        case 1: {
            makeQuad(gl_in[0].gl_Position, model0);
            break;
        }
        case 2: {
            makeQuad(gl_in[0].gl_Position, model45);
            makeQuad(gl_in[0].gl_Position, modelm45);
            break;
        }
        default: {
            makeQuad(gl_in[0].gl_Position, model0);
            makeQuad(gl_in[0].gl_Position, model45);
            makeQuad(gl_in[0].gl_Position, modelm45);
            break;
        }
    }

}


const float LOD1 = 5.f;
const float LOD2 = 10.f;
const float LOD3 = 30.f;


void main()
{
    gs_out.visibility = gs_in[0].visibility;

    vec3 dist_to_camera = gl_in[0].gl_Position.xyz - worldSpace_camPos;
    float dist_length = length(dist_to_camera);

    // add transition for smooth levels
    float t = 6.f;
    if (dist_length > LOD2) t*=1.5f;
    dist_length += (randomize(gl_in[0].gl_Position.xz)*t - t/2.f);

    // change depending on distance
    int detailLevel = 3;
    if (dist_length > LOD1) detailLevel = 2;
    if (dist_length > LOD2) detailLevel = 1;
    if (dist_length > LOD3) detailLevel = 0;

    // make grass with transition levels
    if (detailLevel != 1
      || (detailLevel == 1 && (int(gl_in[0].gl_Position.x * 10) % 1) == 0
      || (int(gl_in[0].gl_Position.z * 10) % 1) == 0)
      || (detailLevel == 2 && (int(gl_in[0].gl_Position.x * 5) % 1) == 0
      || (int(gl_in[0].gl_Position.z * 5) % 1) == 0)
    ) {
        makeGrass(detailLevel);
    }
}


