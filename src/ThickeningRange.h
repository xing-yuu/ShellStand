#pragma once

#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/drawable_octree.h>
#include <cinolib/drawable_segment_soup.h>

#include <cinolib/gl/surface_mesh_controls.h>

#include <cinolib/gl/glcanvas.h>
#include <cinolib/meshes/meshes.h>
#include <cinolib/drawable_octree.h>

#include <cmath>

#include <vector>
#include <queue>
#include <universal.h>
#include <stdio.h>
#include <SelfIntersect.h>

namespace Env {
    void output_range(const cinolib::DrawableTrimesh<> m, std::vector <std::vector<double>> thickness, std::string debug_path, std::string model_name);
    void normal_smooth(const cinolib::DrawableTrimesh<> m, std::vector<cinolib::vec3d>& normals, uint it = 20);
    std::vector<double> thickening(const cinolib::DrawableTrimesh<> m, const std::vector<cinolib::vec3d> normals, int type = 1);
    std::vector <std::vector<double>>* get_thickening_range(const cinolib::DrawableTrimesh<> m);
}