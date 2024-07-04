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
#include <stdio.h>
#include "mosek.h"
/* This function prints log output from MOSEK to the terminal. */
static void MSKAPI printstr(void       *handle,
                            const char str[])
{
  //printf("%s", str);
} /* printstr */


namespace msk {


    std::vector < cinolib::vec3d>* get_poly_center(cinolib::DrawableTrimesh<> m);

    std::vector<double>* MOSEK_deltam(std::vector<double> m0, std::vector<double> rho, cinolib::vec3d target_center, std::vector<double> upper_bound, double lower_bound, cinolib::DrawableTrimesh<> m);

    std::vector < std::vector < uint>>* poly_link_poly(const cinolib::DrawableTrimesh<> m);

    void Mosek_test();
}