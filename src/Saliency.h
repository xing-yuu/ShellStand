#pragma once


#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/drawable_octree.h>
#include <cinolib/drawable_segment_soup.h>

#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/drawable_octree.h>
#include <cinolib/drawable_segment_soup.h>

#include <cmath>



#include <vector>
#include <queue>
#include <universal.h>
#include <maths_funcs.h>
namespace Env {

	std::vector<double>* calculate_mean_curvature(cinolib::DrawableTrimesh<> m);
	std::vector<double>* get_saliency_score(cinolib::DrawableTrimesh<> m, unv::data_type type = unv::POINT_DATA);
}