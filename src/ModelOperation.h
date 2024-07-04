#pragma once

#include <cinolib/meshes/meshes.h>
#include <cinolib/ARAP.h>
#include <cinolib/gl/glcanvas.h>
#include <time.h>
#include <limits.h>
#include <omp.h>
#include <yaml-cpp/yaml.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <filesystem>
#include <algorithm>
#include <Evaluation.h>
#include <SelfIntersect.h>
#include <Deformation.h>
using namespace cinolib;
#include <BayesianOptimization.h>
#include <iostream>
#include <iostream>
#include <string>
#include <filesystem>
#include <direct.h>
#include <ShellThickening.h>
#include <mosek.h>




#include <cmath>
#include <vector>
#include <universal.h>
namespace Ope {
	cinolib::DrawableTrimesh<> scale(const cinolib::DrawableTrimesh<> m, double scale_size = 100);
	cinolib::DrawableTrimesh<> add_base(const cinolib::DrawableTrimesh<> m, cinolib::vec3d center, double radius = 3);
}