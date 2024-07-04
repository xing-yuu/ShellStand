#pragma once
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/drawable_octree.h>
#include <cinolib/drawable_segment_soup.h>
#include <cmath>
#include <vector>

extern cinolib::vec3d target_center_g;

namespace unv {
	enum data_type {
		FACE_DATA, POINT_DATA
	};
	cinolib::vec3d get_projection(const cinolib::vec3d p);
	std::vector <std::vector<double>>* normalize(std::vector <std::vector<double>>* f,double min_f=0, double max_f = 1);
	std::vector<double>* normalize(std::vector<double>* f, double min_f = 0, double max_f = 1);
	std::vector<double>* face2point(const std::vector<double> facedata, const cinolib::DrawableTrimesh<> m);
	std::vector<double>* point2face(const std::vector<double> pointdata, const cinolib::DrawableTrimesh<> m) ;
}


