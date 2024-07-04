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
#include <universal.h>

namespace Env {

    std::vector<cinolib::vec3d>* generate_viewpoints(uint cnt = 200);
    std::vector<cinolib::vec3d>* scale_viewpoints(std::vector<cinolib::vec3d>* viewpoints, cinolib::vec3d center, double scale);
    cinolib::vec3d get_bottom_center(cinolib::DrawableTrimesh<> m);
    double get_radius(cinolib::DrawableTrimesh<> m);
    std::vector <std::vector<double>>* get_visibility_score(cinolib::DrawableTrimesh<> m, unv::data_type type= unv::POINT_DATA, int viewpoints_number = 200, bool debug_mode = true);
}