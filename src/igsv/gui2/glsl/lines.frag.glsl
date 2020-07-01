#version 330
in vec3 gcolor;
out vec4 fcolor;
void main() {
  fcolor = vec4(gcolor, 1.0);
}
