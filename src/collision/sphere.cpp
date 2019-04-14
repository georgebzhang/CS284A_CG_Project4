#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
	// TODO (Part 3): Handle collisions with spheres.
	Vector3D o = pm.position;
	Vector3D sphere_to_pm = o - origin;
	if (sphere_to_pm.norm() < radius) { // if PointMass' position is within sphere
		Vector3D tan_pt = origin + radius * sphere_to_pm.unit(); // tangent point
		Vector3D correctv = tan_pt - pm.last_position; // correction vector
		pm.position = pm.last_position + (1.0 - friction)*correctv;
	}
	//Vector3D d = Vector3D(0, -1, 0);
	//double a = dot(d, d);
	//double b = 2.0 * dot(sphere_to_pm, d);
	//double c = dot(sphere_to_pm, sphere_to_pm) - radius2;
	//double t1, t2;
	//double discr = b * b - 4.0 * a * c;
	//if (discr < 0) return;
	//else {
	//	t1 = (-b + sqrt(discr)) / (2 * a);
	//	t2 = (-b - sqrt(discr)) / (2 * a);
	//}

	//Vector3D tan_pt = o + t2 * d; // tangent point
	//Vector3D correctv = tan_pt - o; // correction vector
	//pm.position = o + (1.0 - friction)*correctv;

	//std::cout << "t1: " << t1 << ", t2: " << t2 << std::endl;
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  m_sphere_mesh.draw_sphere(shader, origin, radius * 0.92);
}
