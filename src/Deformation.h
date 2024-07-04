#pragma once

#include <cinolib/meshes/meshes.h>
#include <cinolib/ARAP.h>
#include <cinolib/gl/glcanvas.h>
#include <cinolib/gl/surface_mesh_controls.h>
#include <time.h>
#include <limits.h>
#include <omp.h>
#include <algorithm>

#include <universal.h>
#include <Evaluation.h>
#include <SelfIntersect.h>
extern std::vector<std::pair<cinolib::DrawableTrimesh<>, double>> meshlist;
namespace Deform {



    struct handle_select {
        handle_select(uint ui, double vi) {
            handle_id = ui;
            handle_score = vi;
        }
        uint handle_id;
        double handle_score;
    };


    bool greater_cmp(handle_select u, handle_select v);
    //::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    //cinolib::vec3d get_projection(const cinolib::vec3d p);
    void get_deform_database(const cinolib::DrawableTrimesh<> m);
    std::vector<uint>* get_static_handles(const cinolib::DrawableTrimesh<> m);
    std::vector<uint>* get_random_handles(uint handle_range, uint handle_number, std::vector<uint>* avoid_list);
    std::vector<handle_select>* score_handle(std::vector<cinolib::vec3d> allhandles, std::vector<cinolib::vec3d> normals, cinolib::vec3d c_current, cinolib::vec3d c_target);
    std::vector<uint>* get_random_handles(const cinolib::DrawableTrimesh<> m, std::vector<uint>* avoid_list, cinolib::vec3d c_current, cinolib::vec3d target_center, int top_num = INT_MAX, uint handle_num = 1);
    cinolib::vec3d get_gravity_center(const cinolib::DrawableTrimesh<> m);
    cinolib::vec3d get_target_center(const cinolib::DrawableTrimesh<> m, std::vector<uint>* handles);
    std::vector<cinolib::vec3d>* get_positive_normals(const cinolib::DrawableTrimesh<> m, cinolib::vec3d gravity_center, cinolib::vec3d target_center, bool IF_TO_CENTER = true);
    std::pair<cinolib::DrawableTrimesh<>, double> deform_ARAP_DATABASE(const cinolib::DrawableTrimesh<> m, double alpha);
    std::vector<cinolib::vec3d>* get_shift_pos(const std::vector<cinolib::vec3d>* ori_vertexs, const std::vector<cinolib::vec3d>* normals, const double distance);
    std::pair<cinolib::DrawableTrimesh<>, double> deform_ARAP(const cinolib::DrawableTrimesh<> m, double alpha);
    bool cmp(std::pair<cinolib::DrawableTrimesh<>, double> u, std::pair<cinolib::DrawableTrimesh<>, double> v);
}