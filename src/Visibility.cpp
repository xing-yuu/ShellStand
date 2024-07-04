#include "Visibility.h"




std::vector<cinolib::vec3d>* Env::generate_viewpoints(uint cnt){
    std::vector<cinolib::vec3d>* viewpoint= new std::vector<cinolib::vec3d>;
    cnt *= 2;
    for (uint i=0;i<cnt;i++){
        double phi=acos(-1.0 + (2.0 * i - 1.0)/cnt);
        double theta = sqrt(cnt * M_PI ) * phi;
        double x = cos(theta) * sin(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(phi);
        if (z >= 0)
            viewpoint->push_back(cinolib::vec3d(x, y, z));
    }
    return viewpoint;
}
std::vector<cinolib::vec3d>* Env::scale_viewpoints(std::vector<cinolib::vec3d>* viewpoints, cinolib::vec3d center,double scale) {
    for (uint i = 0; i < viewpoints->size(); i++) {
        viewpoints->at(i) = viewpoints->at(i) * scale * 1.2;
        viewpoints->at(i) = viewpoints->at(i) + center;
    }
    return viewpoints;
}
cinolib::vec3d Env::get_bottom_center(cinolib::DrawableTrimesh<> m) {
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    double max_x = INT_MIN;
    double min_x = INT_MAX;
    double max_y = INT_MIN;
    double min_y = INT_MAX;
    double min_z = INT_MAX;
    for (uint pid = 0; pid < vertexs.size(); pid++) {
        max_x = std::max(max_x, vertexs[pid][0]);
        min_x = std::min(min_x, vertexs[pid][0]);
        max_y = std::max(max_y, vertexs[pid][1]);
        min_y = std::min(min_y, vertexs[pid][1]);
        min_z = std::min(min_z, vertexs[pid][2]);
    }
    return cinolib::vec3d((max_x + min_x) / 2, (max_y + min_y) / 2, min_z);
}
double Env::get_radius(cinolib::DrawableTrimesh<> m) {
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    double max_x = INT_MIN;
    double min_x = INT_MAX;
    double max_y = INT_MIN;
    double min_y = INT_MAX;
    double max_z = INT_MIN;
    double min_z = INT_MAX;
    for (uint pid = 0; pid < vertexs.size(); pid++) {
        max_x = std::max(max_x, vertexs[pid][0]);
        min_x = std::min(min_x, vertexs[pid][0]);
        max_y = std::max(max_y, vertexs[pid][1]);
        min_y = std::min(min_y, vertexs[pid][1]);
        max_z = std::max(max_z, vertexs[pid][2]);
        min_z = std::min(min_z, vertexs[pid][2]);
    }
    return (cinolib::vec3d((max_x + min_x) / 2, (max_y + min_y) / 2, min_z) - cinolib::vec3d(max_x, max_y, min_z)).norm();
}
std::vector <std::vector<double>>* Env::get_visibility_score(cinolib::DrawableTrimesh<> m, unv::data_type type,int viewpoints_number,bool debug_mode) {
    std::vector<std::vector<uint>> faces = m.vector_polys();
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    std::vector<cinolib::vec3d> normals = m.vector_vert_normals();
    
    std::vector <std::vector<double>>* res_point_ab = new std::vector <std::vector<double>>;
    std::vector<double> di(vertexs.size());
    res_point_ab->push_back(di);
    res_point_ab->push_back(di);
    
    

    uint max_depth = 7;
    uint items_per_leaf = 50;

    //view points
    std::vector<cinolib::vec3d>* viewpoints = generate_viewpoints(viewpoints_number);
    cinolib::vec3d bottem_center = Env::get_bottom_center(m);
    double radius = Env::get_radius(m);
    viewpoints = Env::scale_viewpoints(viewpoints, bottem_center, radius);


    cinolib::DrawableOctree octree(max_depth, items_per_leaf);
    octree.debug_mode(debug_mode); // dump times and statistics
    octree.build_from_mesh_polys(m);
    uint all_block_num = viewpoints->size();
#pragma omp parallel
    {
#pragma omp for
        for (int vid = 0; vid < viewpoints->size(); vid++) {
            printf("\r[%d%%]>", vid * 100 / (all_block_num - 1));
            for (int j = 1; j <= vid * 40 / all_block_num; j++)
                std::cout << "▉";
            for (uint pid = 0; pid < vertexs.size(); pid++) {
                cinolib::vec3d  dir = vertexs[pid] - viewpoints->at(vid);
                double t;
                uint   fid;
                if (octree.intersects_ray(viewpoints->at(vid), dir, t, fid)) {
                    if (faces[fid][0] == pid || faces[fid][1] == pid || faces[fid][2] == pid) {
                        //pid可以被看见，命中，计算法向
                        if (normals[pid].dot(dir) > 0) {
                            res_point_ab->at(0)[pid] = res_point_ab->at(0)[pid] + 1;
                        }
                        else {
                            res_point_ab->at(1)[pid] = res_point_ab->at(1)[pid] + 1;
                        }
                    }
                }
            }
        }
    }
    printf("\n");
    if (type == unv::POINT_DATA) {
        return unv::normalize(res_point_ab);
    }
    else {
        std::vector <std::vector<double>>* res_face_ab = new std::vector <std::vector<double>>;
        res_face_ab->push_back(*unv::point2face(res_point_ab->at(0), m));
        res_face_ab->push_back(*unv::point2face(res_point_ab->at(1), m));
        delete res_point_ab;
        return unv::normalize(res_face_ab);
    }
}