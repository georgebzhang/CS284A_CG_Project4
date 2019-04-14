#version 330

uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

uniform vec4 u_color;

uniform sampler2D u_texture_2;
uniform vec2 u_texture_2_size;

uniform float u_normal_scaling;
uniform float u_height_scaling;

in vec4 v_position;
in vec4 v_normal;
in vec4 v_tangent;
in vec2 v_uv;

out vec4 out_color;

float h(vec2 uv) {
  // You may want to use this helper function...
  return texture(u_texture_2, uv).x;
}

void main() {
  // YOUR CODE HERE
  
  // (Placeholder code. You will want to replace it.)
  vec3 b = cross(v_normal.xyz, v_tangent.xyz);
  mat3 TBN = mat3(v_tangent.xyz, b, v_normal.xyz);
  float kh = 1.0;
  float kn = u_normal_scaling;
  float dU = kh*kn*(h(v_uv + vec2(1/u_texture_2_size.x, 0)) - h(v_uv));
  float dV = kh*kn*(h(v_uv + vec2(0, 1/u_texture_2_size.y)) - h(v_uv));
  vec3 ls_normal = vec3(-dU, -dV, 1.0);
  vec3 dms_normal = TBN * ls_normal;

  float ka = 0.15;
  float kd = 0.5;
  float ks = 1.0;
  vec4 Ia = vec4(1, 0, 1, 0);
  vec4 intensity_falloff = vec4(u_light_intensity, 1.0)/pow(length(vec4(u_light_pos, 1.0) - v_position), 2);
  vec4 l = vec4(u_light_pos, 1.0) - v_position;
  vec4 v = vec4(u_cam_pos, 1.0) - v_position;
  vec4 h = (v + l)/length(v + l);
  int p = 5;

  vec4 ambient = ka*Ia;
  vec4 diffuse = kd * intensity_falloff * max(0, dot(vec4(dms_normal, 1.0), normalize(l)));
  vec4 bp = ks * intensity_falloff * pow(max(0, dot(vec4(dms_normal, 1.0), normalize(h))), p);
  out_color = ambient + diffuse + bp;
  out_color.a = 1;
}

