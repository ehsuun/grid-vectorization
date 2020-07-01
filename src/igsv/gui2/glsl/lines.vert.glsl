#version 330

uniform mat4 model, view, proj;

layout(location=0) in vec2 position;
layout(location=1) in vec3 color;

out vec3 vcolor;
out vec4 position_world;

void main() {
    position_world = model * vec4(position.x, position.y, 0.0, 1.0);
    vcolor = color;
}
