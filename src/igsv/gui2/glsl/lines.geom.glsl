#version 330

layout (lines) in;
layout (triangle_strip,max_vertices=4) out;

uniform mat4 model, view, proj;

uniform float lineThickness;

in vec3 vcolor[];
in vec4 position_world[];

out vec3 gcolor;

void main() {
    vec2 p0 = position_world[0].xy;
    vec2 p1 = position_world[1].xy;
    vec2 tangent = normalize(p1 - p0);
    vec2 normal = vec2(-tangent.y, tangent.x);

    vec2 A = p0 - 0.5 * lineThickness * normal;
    vec2 B = p0 + 0.5 * lineThickness * normal;
    vec2 C = p1 + 0.5 * lineThickness * normal;
    vec2 D = p1 - 0.5 * lineThickness * normal;

    gcolor = vcolor[0];

    gl_Position = proj * view * vec4(A.x, A.y, 0.0, 1.0);
    EmitVertex();
    gl_Position = proj * view * vec4(B.x, B.y, 0.0, 1.0);
    EmitVertex();
    gl_Position = proj * view * vec4(D.x, D.y, 0.0, 1.0);
    EmitVertex();
    gl_Position = proj * view * vec4(C.x, C.y, 0.0, 1.0);
    EmitVertex();
    EndPrimitive();
}
