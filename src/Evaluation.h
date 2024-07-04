#pragma once
#include <algorithm>
#include <cinolib/meshes/meshes.h>
#include <cinolib/ARAP.h>
#include <cinolib/gl/glcanvas.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <omp.h>
#include <cinolib/gl/surface_mesh_controls.h>
namespace Eva {
	cinolib::vec3d rotate_p(cinolib::vec3d p, double angle);
	cinolib::vec2d map_coordinate(cinolib::vec2d p, cinolib::vec2d p_min, cinolib::vec2d p_max);
	double compute_area(const cinolib::vec2d v1, const cinolib::vec2d v2, const cinolib::vec2d v3);
	bool is_point_in_triangle1(const cinolib::vec2d v0, const cinolib::vec2d v1, const cinolib::vec2d v2, const cinolib::vec2d v3);
	double evaluate_similarity(const cinolib::DrawableTrimesh<> m0, const cinolib::DrawableTrimesh<> m, uint evaluate_angle_num = 20);
}

