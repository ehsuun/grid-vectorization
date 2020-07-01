#version 330

in vec3 vcolor;
in vec2 uv;

out vec4 fcolor;

const float alpha = 0.5;

uniform sampler2D imageTexture;
uniform sampler2D maskTexture;
uniform sampler2D gridTexture;

uniform bool enableImageTexture;
uniform bool enableMaskTexture;
uniform bool enableGridTexture;

void main() {

  if (enableGridTexture){
    vec4 grid_color = texture(gridTexture, uv);
    if (grid_color.r < 1 || grid_color.g < 1 || grid_color.b < 1)
      fcolor = vec4(grid_color.r, grid_color.g, grid_color.b, 0.9);
    else
      fcolor = vec4(0.9, 0.9, 0.9, 0.5);
    return;
  }

  vec4 image_color = texture(imageTexture, uv);
  vec4 mask_color  = texture(maskTexture, uv);

  if (enableImageTexture && enableMaskTexture) {
    if (mask_color.a > 0)
      fcolor = (1.0 - alpha) * mask_color * image_color + alpha * mask_color; // simple linear blending
    else
      fcolor = vec4(vcolor,1.0) * image_color;

  } else if (enableImageTexture) {
    fcolor = image_color;

  } else if (enableMaskTexture) {
    fcolor = mask_color;

  } else {
    fcolor = vec4(1.0, 1.0, 1.0, 1.0);
  }
}
