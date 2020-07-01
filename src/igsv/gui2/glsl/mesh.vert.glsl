#version 330
uniform mat4 model, view, proj;
layout(location=0) in vec2 position;
layout(location=1) in vec3 color;
layout(location=2) in vec2 uvcoords;
out vec3 vcolor;
out vec2 uv;
void main() {
  gl_Position = proj * view * model * vec4(position.x, position.y, 0.0, 1.0);
  vcolor = color;
  uv = uvcoords;
}
