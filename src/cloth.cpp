#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Cloth::buildGrid() {
	// TODO (Part 1): Build a grid of masses and springs.
	double d_width = width / (num_width_points - 1);
	double d_height = height / (num_height_points - 1);

	//std::cout << "width: " << width << ", height: " << height << std::endl;
	//std::cout << "num_width_points: " << num_width_points << ", num_height_points: " << num_height_points << std::endl;
	//std::cout << "d_width: " << d_width << ", d_height: " << d_height << std::endl;

	//orientation = HORIZONTAL;
	//point_masses = *(new vector<PointMass>());
	//springs = *(new vector<Spring>());

	if (orientation == HORIZONTAL) {
		std::cout << "HORIZONTAL" << std::endl;
		double z = 0;
		for (int j = 0; j < num_height_points; ++j) {
			double x = 0;
			for (int i = 0; i < num_width_points; ++i) {
				//PointMass pm = PointMass(Vector3D(x, 1, z), false);
				//point_masses.push_back(pm); // pinned false for all point masses
				point_masses.emplace_back(Vector3D(x, 1, z), false); // pinned false for all point masses
				//std::cout << x << ", " << z << std::endl;
				x += d_width;
			}
			z += d_height;
		}

		//std::cout << point_masses.size() << std::endl;
	}
	else if (orientation == VERTICAL) {
		std::cout << "VERTICAL" << std::endl;
		double fmin = -1.0 / 1000;
		double fmax = 1.0 / 1000;
		double y = 0;
		for (int j = 0; j < num_height_points; ++j) {
			double x = 0;
			for (int i = 0; i < num_width_points; ++i) {
				double f = (double)rand() / RAND_MAX;
				double z = fmin + f * (fmax - fmin);
				//point_masses.emplace_back(PointMass(Vector3D(x, y, z), false)); // assume pinned false for all point masses, will update pinned true for pinned point masses later
				point_masses.emplace_back(Vector3D(x, y, z), false);
				x += d_width;
			}
			y += d_height;
		}

	}
	else throw "Orientation not valid!";

	std::cout << "point_masses.size(): " << point_masses.size() << std::endl;
	// update pinned true for pinned point masses
	for (vector<int> coords : pinned) {
		//std::cout << "coords[0]: " << coords[0] << ", coords[1]: " << coords[1] << std::endl;
		point_masses[coords[1] * num_width_points + coords[0]].pinned = true;
	}

	// springs
	// structural springs
	int num_structural_springs = 0;
	for (int j = 1; j < num_height_points; ++j) { // leftmost column
		PointMass* curr = &point_masses[j * num_width_points];
		PointMass* top = &point_masses[(j - 1) * num_width_points];
		//springs.emplace_back(*(new Spring(&top, &curr, STRUCTURAL)));
		//std::cout << top.position << std::endl;
		springs.emplace_back(top, curr, STRUCTURAL);
		//springs.emplace_back(new PointMass(top), new PointMass(curr), STRUCTURAL);
		num_structural_springs += 1;
	}
	for (int i = 1; i < num_width_points; ++i) { // topmost row
		PointMass* curr = &point_masses[i];
		PointMass* left = &point_masses[i - 1];
		//springs.emplace_back(*(new Spring(&left, &curr, STRUCTURAL)));
		springs.emplace_back(left, curr, STRUCTURAL);
		//springs.emplace_back(new PointMass(left), new PointMass(curr), STRUCTURAL);
		num_structural_springs += 1;
	}
	for (int j = 1; j < num_height_points; ++j) { // rest
		for (int i = 1; i < num_width_points; ++i) {
			PointMass* curr = &point_masses[j * num_width_points + i];
			PointMass* left = &point_masses[j * num_width_points + i - 1];
			PointMass* top = &point_masses[(j - 1) * num_width_points + i];
			//springs.emplace_back(*(new Spring(&left, &curr, STRUCTURAL)));
			//springs.emplace_back(*(new Spring(&top, &curr, STRUCTURAL)));
			springs.emplace_back(left, curr, STRUCTURAL);
			springs.emplace_back(top, curr, STRUCTURAL);
			//springs.emplace_back(new PointMass(left), new PointMass(curr), STRUCTURAL);
			//springs.emplace_back(new PointMass(top), new PointMass(curr), STRUCTURAL);
			num_structural_springs += 2;
		}
	}
	// shearing springs
	int num_shear_springs = 0;
	for (int j = 1; j < num_height_points; ++j) { // leftmost column
		PointMass* curr = &point_masses[j * num_width_points];
		PointMass* topright = &point_masses[(j - 1) * num_width_points + 1];
		//springs.emplace_back(*(new Spring(&topright, &curr, SHEARING)));
		springs.emplace_back(topright, curr, SHEARING);
		//springs.emplace_back(new PointMass(topright), new PointMass(curr), SHEARING);
		num_shear_springs += 1;
	}
	for (int j = 1; j < num_height_points; ++j) { // rightmost column
		int i = num_width_points - 1;
		PointMass* curr = &point_masses[j * num_width_points + i];
		PointMass* topleft = &point_masses[(j - 1) * num_width_points + i - 1];
		//springs.emplace_back(*(new Spring(&topleft, &curr, SHEARING)));
		springs.emplace_back(topleft, curr, SHEARING);
		//springs.emplace_back(new PointMass(topleft), new PointMass(curr), SHEARING);
		num_shear_springs += 1;
	}
	for (int j = 1; j < num_height_points; ++j) { // rest
		for (int i = 1; i < num_width_points - 1; ++i) {
			PointMass* curr = &point_masses[j * num_width_points + i];
			PointMass* topleft = &point_masses[(j - 1) * num_width_points + i - 1];
			PointMass* topright = &point_masses[(j - 1) * num_width_points + i + 1];
			//springs.emplace_back(*(new Spring(&topleft, &curr, SHEARING)));
			//springs.emplace_back(*(new Spring(&topright, &curr, SHEARING)));
			springs.emplace_back(topleft, curr, SHEARING);
			springs.emplace_back(topright, curr, SHEARING);
			//springs.emplace_back(new PointMass(topleft), new PointMass(curr), SHEARING);
			//springs.emplace_back(new PointMass(topright), new PointMass(curr), SHEARING);
			num_shear_springs += 2;
		}
	}
	// bending springs
	int num_bending_springs = 0;
	for (int j = 2; j < num_height_points; ++j) { // bottom left
		for (int i = 0; i < 2; ++i) {
			PointMass* curr = &point_masses[j * num_width_points + i];
			PointMass* top = &point_masses[(j - 2) * num_width_points + i];
			//springs.emplace_back(*(new Spring(&top, &curr, BENDING)));
			springs.emplace_back(top, curr, BENDING);
			//springs.emplace_back(new PointMass(top), new PointMass(curr), BENDING);
			num_bending_springs += 1;
		}
	}
	for (int j = 0; j < 2; ++j) { // top right
		for (int i = 2; i < num_width_points; ++i) {
			PointMass* curr = &point_masses[j * num_width_points + i];
			PointMass* left = &point_masses[j * num_width_points + i - 2];
			//springs.emplace_back(*(new Spring(&left, &curr, BENDING)));
			springs.emplace_back(left, curr, BENDING);
			//springs.emplace_back(new PointMass(left), new PointMass(curr), BENDING);
			num_bending_springs += 1;
		}
	}
	for (int j = 2; j < num_height_points; ++j) { // rest
		for (int i = 2; i < num_width_points; ++i) {
			PointMass* curr = &point_masses[j * num_width_points + i];
			PointMass* left = &point_masses[j * num_width_points + i - 2];
			PointMass* top = &point_masses[(j - 2) * num_width_points + i];
			//springs.emplace_back(*(new Spring(&left, &curr, BENDING)));
			//springs.emplace_back(*(new Spring(&top, &curr, BENDING)));
			springs.emplace_back(left, curr, BENDING);
			springs.emplace_back(top, curr, BENDING);
			//springs.emplace_back(new PointMass(left), new PointMass(curr), BENDING);
			//springs.emplace_back(new PointMass(top), new PointMass(curr), BENDING);
			num_bending_springs += 2;
		}
	}
	//std::cout << "num_structural_springs: " << num_structural_springs << std::endl;
	//std::cout << "num_shear_springs: " << num_shear_springs << std::endl;
	//std::cout << "num_bending_springs: " << num_bending_springs << std::endl;
	//std::cout << "in cloth.cpp: " << springs.size() << std::endl;
	//std::cout << "in cloth.cpp: " << point_masses.size() << std::endl;
	//std::cout << "in cloth.cpp: " << &(springs[0].pm_a) << std::endl;
	//std::cout << "in cloth.cpp: " << &(springs[0].pm_a->position) << std::endl;
	//std::cout << "in cloth.cpp: " << springs[0].pm_a->position << std::endl;
	//std::cout << "in cloth.cpp: " << point_masses[0].position << std::endl;
}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
	double mass = width * height * cp->density / num_width_points / num_height_points;
	double delta_t = 1.0f / frames_per_sec / simulation_steps;
	//std::cout << "simulate" << std::endl;

	// TODO (Part 2): Compute total force acting on each point mass.
	//std::cout << "external_accelerations size: " << external_accelerations.size() << std::endl;
	//std::cout << "external_accelerations: " << external_accelerations[0] << std::endl;
	//std::cout << "point_masses size: " << point_masses.size() << std::endl;
	//for (PointMass &pm : point_masses) {
	//	pm.forces = Vector3D();
	//}
	for (PointMass &pm : point_masses) {
		pm.forces = mass * external_accelerations[0];
	}
	//for (Vector3D v : external_accelerations) {
	//	std::cout << "external_accelerations: " << v << std::endl;
	//}
	//std::cout << "test" << std::endl;
	for (Spring &s : springs) {
		if (cp->enable_bending_constraints || cp->enable_shearing_constraints || cp->enable_structural_constraints) {
			double Fs = cp->ks * ((s.pm_a->position - s.pm_b->position).norm() - s.rest_length);
			if (s.spring_type == BENDING) Fs *= 0.2;
			Vector3D Fsv = Fs * ((s.pm_b->position - s.pm_a->position).unit());
			s.pm_a->forces += Fsv;
			s.pm_b->forces -= Fsv;
		}
	}

	// TODO (Part 2): Use Verlet integration to compute new point mass positions
	//std::cout << point_masses[1].pinned << std::endl;
	//std::cout << point_masses[1].position << std::endl;
	for (PointMass &pm : point_masses) {
		if (!pm.pinned) {
			Vector3D xt = pm.position;
			//std::cout << xt << std::endl;
			Vector3D xtmdt = pm.last_position;
			Vector3D at = pm.forces / mass;
			pm.position = xt + (1.0 - cp->damping/100.0)*(xt - xtmdt) + at * pow(delta_t, 2);
			pm.last_position = xt;
		}
	}

	// TODO (Part 4): Handle self-collisions.
	build_spatial_map();
	for (PointMass &pm : point_masses) {
		self_collide(pm, simulation_steps);
	}

	// TODO (Part 3): Handle collisions with other primitives.
	//int i = 0;
	for (PointMass &pm : point_masses) {
		//std::cout << i << "th pm" << std::endl;
		for (CollisionObject* co : *collision_objects) {
			co->collide(pm);
		}
		//++i;
	}

	// TODO (Part 2): Constrain the changes to be such that the spring does not change
	// in length more than 10% per timestep [Provot 1995].
	for (Spring s : springs) {
		//if (cp->enable_bending_constraints || cp->enable_shearing_constraints || cp->enable_structural_constraints) {
		if (s.pm_a->pinned && s.pm_b->pinned) continue; // both PointMasses pinned
		double diff = (s.pm_b->position - s.pm_a->position).norm() - 1.1*s.rest_length;
		if (diff <= 0) continue; // Spring's length less than 110% of its rest length
		Vector3D diffv = diff * (s.pm_b->position - s.pm_a->position).unit();
		if (s.pm_a->pinned) {
			s.pm_b->position -= diffv;
		}
		else if (s.pm_b->pinned) {
			s.pm_a->position += diffv;
		}
		else {
			diffv /= 2.0;
			s.pm_b->position -= diffv;
			s.pm_a->position += diffv;
		}
		//}
	}
}

void Cloth::build_spatial_map() {
	for (const auto &entry : map) {
		delete(entry.second);
	}
	map.clear();

	// TODO (Part 4): Build a spatial map out of all of the point masses.
	for (PointMass &pm : point_masses) {
		float hash = hash_position(pm.position);
		if (map[hash] == nullptr) {
			map[hash] = new vector<PointMass*>;
		}
		map[hash]->push_back(&pm);
	}

	//int nbuckets = map.bucket_count();

	//std::cout << "nbuckets: " << nbuckets << std::endl;

	//for (int i = 0; i < nbuckets; ++i) {
	//	std::cout << "bucket #" << i << " has size: " << map.bucket_size(i) << std::endl;
	//}

	//std::cout << "begin" << std::endl;
	//for (auto itr = map.begin(); itr != map.end(); itr++) {
	//	std::cout << itr->first << "  " << itr->second->size() << std::endl;
	//}
	//std::cout << "end" << std::endl;
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
	// TODO (Part 4): Handle self-collision for a given point mass.
	//std::cout << "thickness: " << thickness << std::endl;
	Vector3D correctv = Vector3D();
	float hash = hash_position(pm.position);
	int num = 0;
	for (PointMass* pm_cand : *map[hash]) { // for each candidate PointMass
		if (&pm == pm_cand) continue; // do not collide a PointMass with itself
		double diff = 2.0*thickness - ((pm.position - pm_cand->position).norm());
		if (diff < 0) continue;
		correctv += diff * ((pm.position - pm_cand->position).unit());
		++num;
	}
	//std::cout << "num: " << num << std::endl;
	if (num == 0) return;
	//std::cout << "correctv: " << correctv << std::endl;
	correctv /= num;
	//std::cout << "simulation_steps: " << simulation_steps << std::endl;
	correctv /= simulation_steps;
	pm.position += correctv;
}

float Cloth::hash_position(Vector3D pos) {
	// TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
	//std::cout << "pos.x: " << pos.x << ", pos.y: " << pos.y << ", pos.z: " << pos.z << std::endl;
	double w = 3.0 * width / (num_width_points - 1);
	double h = 3.0 * height / (num_height_points - 1);
	double t = max(w, h);
	//std::cout << "w: " << w << ", h: " << h << ", t: " << t << std::endl;
	//std::cout << fmod(pos.x, w) << std::endl;
	//float trunc_x = (pos.x - fmod(pos.x, w))/w;
	//float trunc_y = (pos.y - fmod(pos.y, h))/h;
	//float trunc_z = (pos.z - fmod(pos.z, t))/t;
	float trunc_x = pos.x / w;
	float trunc_y = pos.y / h;
	float trunc_z = pos.z / t;
	//std::cout << "trunc_x: " << trunc_x << ", trunc_y: " << trunc_y << ", trunc_z: " << trunc_z << std::endl;
	int trunc[3] = { trunc_x, trunc_y, trunc_z };
	float hash = 7.0;
	for (int i = 0; i < 3; ++i) {
		hash = hash * 31.0 + trunc[i];
	}
	return hash; 
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;
  std::cout << "Cloth::buildClothMesh" << std::endl;
  //std::cin.get();

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  // num_XXX_points - 1 because of image below
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
	  //std::cout << pm->position << std::endl;
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
