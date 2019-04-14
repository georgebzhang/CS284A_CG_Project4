#version 330

uniform vec4 u_color;
uniform vec3 u_cam_pos;
uniform vec3 u_light_pos;
uniform vec3 u_light_intensity;

in vec4 v_position;
in vec4 v_normal;
in vec2 v_uv;

out vec4 out_color;

void main() {
  // YOUR CODE HERE
  
  // (Placeholder code. You will want to replace it.)
  float ka = 0.15;
  float kd = 0.5;
  float ks = 1.0;
  vec4 Ia = vec4(0, 1, 1, 0);
  vec4 intensity_falloff = vec4(u_light_intensity, 1.0)/pow(length(vec4(u_light_pos, 1.0) - v_position), 2);
  vec4 l = vec4(u_light_pos, 1.0) - v_position;
  vec4 v = vec4(u_cam_pos, 1.0) - v_position;
  vec4 h = (v + l)/length(v + l);
  int p = 5;

  vec4 ambient = ka*Ia;
  vec4 diffuse = kd * intensity_falloff * max(0, dot(v_normal, normalize(l)));
  vec4 bp = ks * intensity_falloff * pow(max(0, dot(v_normal, normalize(h))), p);
  out_color = ambient + diffuse + bp;
  out_color.a = 1;
}

