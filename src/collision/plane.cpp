#include "iostream"
#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../clothSimulator.h"
#include "../leak_fix.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide(PointMass &pm) {
	// TODO (Part 3): Handle collisions with planes.
	Vector3D o = pm.last_position;
	Vector3D d = -normal;
	// pm.last_position may be on the side of the plane with normal or it may be on the opposite side, thus we should test two rays with opposite directions (d, -d)
	double t1 = dot((point - o), normal) / dot(d, normal);
	double t2 = dot((point - o), normal) / dot(-d, normal);
	double t;
	// if either t1 or t2 is a valid time, store that in t, else return
	if (t1 >= 0 && t1 < INF_D) t = t1;
	else if (t2 >= 0 && t2 < INF_D) t = t2;
	else return;

	Vector3D tan_pt = o + t * d; // tangent point
	Vector3D correctv = tan_pt - pm.last_position; // correction vector
	if (dot(correctv, normal) < 0) correctv += (normal*SURFACE_OFFSET); // if pm.last_position is on same side of plane as normal is pointing
	else correctv -= (normal*SURFACE_OFFSET);
	// sign is defined as positive if the point is on the side of the plane with normal
	int sign_position = dot(tan_pt - pm.position, normal) < 0 ? 1 : -1;
	int sign_last_position = dot(tan_pt - pm.last_position, normal) < 0 ? 1 : -1;
	// if signs are unequal, then position and last_position are on opposite sides of the plane, thus we collide
	if (sign_position != sign_last_position) pm.position = pm.last_position + (1.0 - friction)*correctv;
}

void Plane::render(GLShader &shader) {
  nanogui::Color color(0.7f, 0.7f, 0.7f, 1.0f);

  Vector3f sPoint(point.x, point.y, point.z);
  Vector3f sNormal(normal.x, normal.y, normal.z);
  Vector3f sParallel(normal.y - normal.z, normal.z - normal.x,
                     normal.x - normal.y);
  sParallel.normalize();
  Vector3f sCross = sNormal.cross(sParallel);

  MatrixXf positions(3, 4);
  MatrixXf normals(3, 4);

  positions.col(0) << sPoint + 2 * (sCross + sParallel);
  positions.col(1) << sPoint + 2 * (sCross - sParallel);
  positions.col(2) << sPoint + 2 * (-sCross + sParallel);
  positions.col(3) << sPoint + 2 * (-sCross - sParallel);

  normals.col(0) << sNormal;
  normals.col(1) << sNormal;
  normals.col(2) << sNormal;
  normals.col(3) << sNormal;

  if (shader.uniform("u_color", false) != -1) {
    shader.setUniform("u_color", color);
  }
  shader.uploadAttrib("in_position", positions);
  if (shader.attrib("in_normal", false) != -1) {
    shader.uploadAttrib("in_normal", normals);
  }

  shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);
#ifdef LEAK_PATCH_ON
  shader.freeAttrib("in_position");
  if (shader.attrib("in_normal", false) != -1) {
    shader.freeAttrib("in_normal");
  }
#endif
}
