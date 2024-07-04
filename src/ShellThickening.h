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
#include <ThickeningRange.h>
#include <SelfIntersect.h>
#include <Saliency.h>
#include <Visibility.h>
#include <Deformation.h>
#include <MOSEKOptimize.h>
#include <universal.h>
#include <ModelOperation.h>
#include <universal.h>
#include <io.h>
#include <direct.h>
#include <string>

extern cinolib::vec3d target_center_g;
namespace Env {
	double get_thickening_score(std::vector<double> RHO, std::vector<double>delta_m);
	std::pair<std::vector<double>, double> shell_THICKENING(const cinolib::DrawableTrimesh<> m, std::pair<std::vector<double>, std::vector<double>> RHO_RANGE);
	std::pair<std::vector<double>, double> shell_THICKENING(const cinolib::DrawableTrimesh<> m);
	void output_volume_model(const std::string output, const cinolib::DrawableTrimesh<> m, const std::vector<double> offset);

	void Mosek_test();
}

