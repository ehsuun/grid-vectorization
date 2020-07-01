#version 330
in vec3 vcolor;
out vec4 fcolor;
void main() {
  //  fcolor = vec4(1.0, 1.0, 0.0, 1.0);
  fcolor = vec4(vcolor, 1.0);
}
