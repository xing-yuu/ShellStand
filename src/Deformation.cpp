#include "Deformation.h"

//using namespace Deform;
std::vector<std::pair<cinolib::DrawableTrimesh<>, double>> meshlist;
bool Deform::greater_cmp(handle_select u, handle_select v) {
    return u.handle_score > v.handle_score;
}
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

std::vector<uint>* Deform::get_static_handles(const cinolib::DrawableTrimesh<> m) {

    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    double min_z = vertexs[0][2];
    std::vector<uint>* min_index = new std::vector<uint>;
    for (uint vid = 0; vid < m.num_verts(); ++vid)
    {
        if (min_z > vertexs[vid][2])
            min_z = vertexs[vid][2];
    }
    for (uint vid = 0; vid < m.num_verts(); ++vid)
    {
        if (min_z == vertexs[vid][2]) {
            min_index->push_back(vid);
        }

    }
    //uint vid = m.pick_vert(vi);
    return min_index;
}
std::vector<uint>* Deform::get_random_handles(uint handle_range, uint handle_number, std::vector<uint>* avoid_list) {
    std::vector<uint>* random_handles = new std::vector<uint>;
    std::set<uint> in;
    for (int i = 0; i < avoid_list->size(); i++) {
        in.insert(avoid_list->at(i));
    }
    for (int i = 0; i < handle_number; i++) {
        uint rand_handle = rand() % handle_range;
        while (in.find(rand_handle) != in.end()) rand_handle = rand() % handle_range;
        random_handles->push_back(rand_handle);
    }
    return random_handles;
}
std::vector<Deform::handle_select>* Deform::score_handle(std::vector<cinolib::vec3d> allhandles, std::vector<cinolib::vec3d> normals, cinolib::vec3d c_current, cinolib::vec3d c_target) {
    std::vector<Deform::handle_select>* hle = new std::vector<Deform::handle_select>;//result
    std::vector<double> s_dis;
    std::vector<double> s_ang;
    double dis_max = INT_MIN;
    double dis_min = INT_MAX;
    double ang_max = INT_MIN;
    double ang_min = INT_MAX;

    for (uint i = 0; i < allhandles.size(); i++) {
        cinolib::vec3d std(unv::get_projection((c_target - c_current) / (c_current - c_target).norm()));
        s_dis.push_back(fabs(unv::get_projection(allhandles[i] - c_current).dot(std)));
        s_ang.push_back(fabs(normals[i].dot(std)));
        dis_max = std::max(dis_max, s_dis[i]);
        dis_min = std::min(dis_min, s_dis[i]);
        ang_max = std::max(ang_max, s_ang[i]);
        ang_min = std::max(ang_min, s_ang[i]);
    }
    for (uint i = 0; i < allhandles.size(); i++) {
        s_dis[i] = (s_dis[i] - dis_min) / (dis_max - dis_min);
        s_ang[i] = (s_ang[i] - ang_min) / (ang_max - ang_min);
        s_ang[i] = 0;
        Deform::handle_select vi(i, (s_dis[i] + s_ang[i]) / 2);
        hle->push_back(vi);
    }
    std::sort(hle->begin(), hle->end(), greater_cmp);
    return hle;
}
std::vector<uint>* Deform::get_random_handles(const cinolib::DrawableTrimesh<> m, std::vector<uint>* avoid_list, cinolib::vec3d c_current, cinolib::vec3d target_center, int top_num , uint handle_num ) {
    std::vector<cinolib::vec3d> normals = m.vector_vert_normals();
    std::vector<uint>* random_handles = new std::vector<uint>;
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    std::vector<Deform::handle_select>* hle = score_handle(vertexs, normals, c_current, target_center);
    uint handle_range = vertexs.size();
    uint max_handle = 0;
    std::set<uint> in;
    if (INT_MAX == top_num) {
        top_num = handle_range * 0.2;
    }
    double max_length = -1;
    for (int i = 0; i < avoid_list->size(); i++) {
        in.insert(avoid_list->at(i));
    }
    for (uint i = 0; i < handle_num; i++) {
        uint rand_handle = rand() % top_num;
        while (in.find(hle->at(rand_handle).handle_id) != in.end()) rand_handle = rand() % top_num;
        random_handles->push_back(hle->at(rand_handle).handle_id);
        in.insert(hle->at(rand_handle).handle_id);
    }
    delete hle;
    return random_handles;
}

cinolib::vec3d Deform::get_gravity_center(const cinolib::DrawableTrimesh<> m) {
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    std::vector<std::vector<uint>> faces = m.vector_polys();
    double mesh_area = 0;
    cinolib::vec3d center = cinolib::vec3d(0, 0, 0);
    for (int i = 0; i < m.num_polys(); i++) {
        cinolib::vec3d p_i = vertexs[faces[i][0]];
        cinolib::vec3d p_j = vertexs[faces[i][1]];
        cinolib::vec3d p_k = vertexs[faces[i][2]];
        cinolib::vec3d c = (p_i + p_j + p_k) / 3;
        cinolib::vec3d u = p_j - p_i;
        cinolib::vec3d v = p_k - p_i;
        cinolib::vec3d area_temp = u.cross(v);
        double area = area_temp.norm() * 0.5;
        mesh_area += area;
        center += c * area;
    }
    return center / mesh_area;
}
cinolib::vec3d Deform::get_target_center(const cinolib::DrawableTrimesh<> m, std::vector<uint>* handles) {
    cinolib::vec3d center = cinolib::vec3d(0, 0, 0);
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    for (int i = 0; i < handles->size(); i++) {
        center += vertexs[handles->at(i)];
    }
    return center /= handles->size();
} 
std::vector<cinolib::vec3d>* Deform::get_positive_normals(const cinolib::DrawableTrimesh<> m, cinolib::vec3d gravity_center, cinolib::vec3d target_center, bool IF_TO_CENTER) {
    std::vector<cinolib::vec3d> normals = m.vector_vert_normals();

    //to center
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    if (IF_TO_CENTER)
        for (int i = 0; i < m.num_verts(); ++i) {
            cinolib::vec3d vi = target_center - vertexs.at(i);
            vi = vi / vi.norm();
            normals[i] = vi;
        }
    std::vector<cinolib::vec3d>* positive_normals = new std::vector<cinolib::vec3d>;
    cinolib::vec3d n1 = unv::get_projection(target_center - gravity_center);
    for (int i = 0; i < m.num_verts(); ++i) {
        if (n1.dot(normals[i]) < 0) positive_normals->push_back(normals[i] * (-1.0));
        else positive_normals->push_back(unv::get_projection(normals[i]));
    }
    return positive_normals;
}

std::vector<cinolib::vec3d>* Deform::get_shift_pos(const std::vector<cinolib::vec3d>* ori_vertexs, const std::vector<cinolib::vec3d>* normals, const double distance) {
    std::vector<cinolib::vec3d>* shift_vertexs = new std::vector<cinolib::vec3d>;
    for (int i = 0; i < ori_vertexs->size(); i++) {
        shift_vertexs->push_back(ori_vertexs->at(i) + normals->at(i) * distance);
    }
    return shift_vertexs;
}
std::pair<cinolib::DrawableTrimesh<>, double> Deform::deform_ARAP_DATABASE(const cinolib::DrawableTrimesh<> m0,double alpha) {
    alpha = std::max(alpha, 0.05);
    int min_Id = -1;
    double min_score = MAXINT;
    for (uint i = 0; i < meshlist.size(); i++) {
        if (meshlist[i].second<alpha)
        {
            double per_score = Eva::evaluate_similarity(m0, meshlist[i].first);
            if (per_score < min_score) {
                min_score = per_score;
                min_Id = i;
            }
        }
        if (meshlist[i].second < alpha * 0.9 && min_Id != -1) {
            break;
        }
    }
    return std::make_pair(meshlist[min_Id].first, min_score);
}

std::pair<cinolib::DrawableTrimesh<>, double> Deform::deform_ARAP(const cinolib::DrawableTrimesh<> m, double alpha) {
    alpha = std::max(alpha, 0.05);
    std::vector<cinolib::vec3d>vertexs = m.vector_verts();
    double xMin = INT_MAX, yMin = INT_MAX, zMin = INT_MAX;
    double xMax = INT_MIN, yMax = INT_MIN, zMax = INT_MIN;
    for (uint i = 0; i < vertexs.size(); i++) {
        xMin = std::min(xMin, vertexs[i][0]);
        xMax = std::max(xMax, vertexs[i][0]);
        yMin = std::min(yMin, vertexs[i][1]);
        yMax = std::max(yMax, vertexs[i][1]);
        zMin = std::min(zMin, vertexs[i][2]);
        zMax = std::max(zMax, vertexs[i][2]);
    }
    float diagonalLength = sqrt((xMax - xMin) * (xMax - xMin) + (yMax - yMin) * (yMax - yMin) + (zMax - zMin) * (zMax - zMin));
    std::vector<uint>* static_handles = get_static_handles(m);
    int n_init = 20;
    cinolib::vec3d ori_gravity_center = get_gravity_center(m);
    cinolib::vec3d ori_target_center = get_target_center(m, static_handles);
    double per_step = (ori_target_center - ori_gravity_center).norm() * 0.4;
    double tar_dis = unv::get_projection((ori_target_center - ori_gravity_center)).norm() * alpha;
    std::vector<cinolib::DrawableTrimesh<>> res_model_list(n_init);
    std::vector < double>  score(n_init + 2);
    for (int i = 0; i < n_init; i++) {
        score[i] = INT_MAX;
    }
    clock_t start, finish;
    start = clock();
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < n_init; i++) {
            double per_step = unv::get_projection((ori_target_center - ori_gravity_center)).norm() * 0.01;
            cinolib::DrawableTrimesh<> model(m);
            cinolib::DrawableTrimesh<> last_model(m);
            cinolib::vec3d gravity_center = ori_gravity_center;
            int handle_num = 3;
            uint deform_num = 0;
            clock_t start_per, finish_per;
            start_per = clock();
            int tip = 1;
            double last_dis = unv::get_projection(ori_target_center - gravity_center).norm();
            bool IF_VALID = true;
            int MAX_DEFORM = 50;

            while (unv::get_projection(ori_target_center - gravity_center).norm() > tar_dis && deform_num < MAX_DEFORM) {

                per_step = diagonalLength * 0.05;
                cinolib::ARAP_data data;
                data.n_iters = 5;
               // data.n_iters = 20;
                std::vector<cinolib::vec3d> all_vertexs = model.vector_verts();
                //std::vector<uint>* random_handles = get_random_handles(model.num_verts(), handle_num, static_handles);
                std::vector<uint>* random_handles = get_random_handles(model, static_handles, gravity_center, ori_target_center);
                std::vector<cinolib::vec3d>* normals = get_positive_normals(model, gravity_center, ori_target_center, false);
                for (int i = 0; i < static_handles->size(); i++) {
                    data.bcs[static_handles->at(i)] = all_vertexs[static_handles->at(i)];
                }
                //per_step = get_projection(ori_target_center - gravity_center).norm() * 0.2;
                for (int i = 0; i < random_handles->size(); i++) {
                    cinolib::vec3d normal_plane = unv::get_projection(normals->at(random_handles->at(i)));
                    data.bcs[random_handles->at(i)] = all_vertexs[random_handles->at(i)] + normal_plane * per_step;
                }
                clock_t s_per, f_per;
                s_per = clock();

                ARAP(model, data);
                gravity_center = get_gravity_center(model);
                deform_num++;
                if (Eva::SelfIntersection(model.vector_verts(), model.vector_polys()) || unv::get_projection(ori_target_center - gravity_center).norm() > last_dis) {
                    model = last_model;
                    for (int ti = 1; ti <= 5; ti++) {
//                        std::cout << "EEEOR\t" << ti << std::endl;
                        per_step /= 2;
                        for (int i = 0; i < random_handles->size(); i++) {
                            cinolib::vec3d normal_plane = unv::get_projection(normals->at(random_handles->at(i)));
                            data.bcs[random_handles->at(i)] = all_vertexs[random_handles->at(i)] + normal_plane * per_step;
                        }
                        ARAP(model, data);
                        cinolib::vec3d temp_gravity_center= get_gravity_center(model);
                        if (!(Eva::SelfIntersection(model.vector_verts(), model.vector_polys()) || unv::get_projection(ori_target_center - temp_gravity_center).norm() > last_dis)) {
                            break;
                        }
                        else {
                            model = last_model;
                        }
                    }
                }
                gravity_center = get_gravity_center(model);
                last_dis = unv::get_projection(ori_target_center - gravity_center).norm();
                last_model = model;
                f_per = clock();
                //std::cout << "Per deformation:" << double(f_per - s_per) / CLOCKS_PER_SEC <<"target dis:\t"<< tar_dis<<"\t now dis:\t"<<(ori_target_center - gravity_center).norm() << std::endl;

                //model.updateGL();
                //gui.pop_all_markers();
                //for (uint vid : *random_handles) gui.push_marker(model.vert(vid), "", cinolib::Color::BLUE(), 10);
                //for (uint vid : *static_handles) gui.push_marker(model.vert(vid), "", cinolib::Color::GREEN(), 25);



                std::string s = std::to_string(tip) + ".obj";
//                model.save(s.data());
                tip++;
                delete random_handles;
                delete normals;

               // std::cout << "Now dis:\t" << unv::get_projection(ori_target_center - gravity_center).norm() << "\t\tTarget dis:\t" << tar_dis << std::endl;
            }
            double per_score = MAXINT;
            if (deform_num < MAX_DEFORM) {
                per_score = Eva::evaluate_similarity(m, model);

                finish_per = clock();
                //std::cout << i << "\tFinished\t" << deform_num << "\t the time cost is:" << double(finish_per - start_per) / CLOCKS_PER_SEC << std::endl;
                std::string s = std::to_string(i) + ".obj";
                model.save(s.data());
            }
#pragma omp critical
            {
                score[i] = per_score;
                res_model_list[i] = model;
            }
            printf("%d\n", i);
        }
    }
    double minele = MAXINT;
    int min_num_index = 0;
    for (int i = 0; i < n_init; i++) {
        if (minele> score[i])
        {
            minele = score[i];
            min_num_index = i;
        }
        cout << score[i] << endl;
    }
    cout << "MIN SCORE INDEX:\t" << min_num_index << "\tMIN SCORE:\t" << score[min_num_index] << std::endl;
    finish = clock();
    res_model_list[min_num_index].save("Deform_result.obj");
    std::cout << "the time cost is:" << double(finish - start) / CLOCKS_PER_SEC << std::endl;
    std::pair<cinolib::DrawableTrimesh<>, double> ans = std::make_pair(res_model_list[min_num_index], score[min_num_index]);
    return ans;
}


void Deform::get_deform_database(const cinolib::DrawableTrimesh<> m) {
    double alpha = 0.05;
    std::cout << "========================" << std::endl;
    std::cout << "Deform Database " << std::endl;
    std::cout << "========================" << std::endl;
    std::vector<cinolib::vec3d>vertexs = m.vector_verts();
    double xMin = INT_MAX, yMin = INT_MAX, zMin = INT_MAX;
    double xMax = INT_MIN, yMax = INT_MIN, zMax = INT_MIN;
    for (uint i = 0; i < vertexs.size(); i++) {
        xMin = std::min(xMin, vertexs[i][0]);
        xMax = std::max(xMax, vertexs[i][0]);
        yMin = std::min(yMin, vertexs[i][1]);
        yMax = std::max(yMax, vertexs[i][1]);
        zMin = std::min(zMin, vertexs[i][2]);
        zMax = std::max(zMax, vertexs[i][2]);
    }
    float diagonalLength = sqrt((xMax - xMin) * (xMax - xMin) + (yMax - yMin) * (yMax - yMin) + (zMax - zMin) * (zMax - zMin));
    std::vector<uint>* static_handles = get_static_handles(m);
    int n_init = 20;
    cinolib::vec3d ori_gravity_center = get_gravity_center(m);
    cinolib::vec3d ori_target_center = get_target_center(m, static_handles);
    target_center_g = ori_target_center;
    double per_step = (ori_target_center - ori_gravity_center).norm() * 0.4;
    double tar_dis = unv::get_projection((ori_target_center - ori_gravity_center)).norm() * alpha;
    std::vector<cinolib::DrawableTrimesh<>> res_model_list(n_init);
    std::vector < double>  score(n_init + 2);
    std::vector < std::vector<std::pair<cinolib::DrawableTrimesh<>, double>>> all_mesh(n_init);
    for (int i = 0; i < n_init; i++) {
        score[i] = INT_MAX;
    }
    clock_t start, finish;
    start = clock();
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < n_init; i++) {
            printf("\r 1/2 [%.2f%%]\t>", float((i + 1) * 100.0 / n_init));
            for (int j = 1; j <= i * 40 / n_init; j++)
                std::cout << "¨¨";
            double per_step = unv::get_projection((ori_target_center - ori_gravity_center)).norm() * 0.005;
            cinolib::DrawableTrimesh<> model(m);
            cinolib::DrawableTrimesh<> last_model(m);
            cinolib::vec3d gravity_center = ori_gravity_center;
            int handle_num = 3;
            uint deform_num = 0;
            clock_t start_per, finish_per;
            start_per = clock();
            int tip = 1;
            /*         cinolib::GLcanvas gui;
                       gui.push(&model);
                       gui.push(new cinolib::SurfaceMeshControls<cinolib::DrawableTrimesh<>>(&model, &gui));
           */
            double last_dis = unv::get_projection(ori_target_center - gravity_center).norm();
            bool IF_VALID = true;
            int MAX_DEFORM = 50;

            std::vector<std::pair<cinolib::DrawableTrimesh<>, double>> temp_meshlist;
            while (unv::get_projection(ori_target_center - gravity_center).norm() > tar_dis && deform_num < MAX_DEFORM) {
                //                
                //                cout << "DEFORM\t" << deform_num << "\tNOW DIS:\t" << unv::get_projection(ori_target_center - gravity_center).norm() << endl;

                                //std::cout << "deform_num\t" << deform_num << std::endl;
                                //per_step = std::max(unv::get_projection((ori_target_center - gravity_center)).norm() * 0.1, diagonalLength * 0.01);
                per_step = diagonalLength * 0.05;
                cinolib::ARAP_data data;
                //data.n_iters = 5;
                data.n_iters = 50;
                std::vector<cinolib::vec3d> all_vertexs = model.vector_verts();
                //std::vector<uint>* random_handles = get_random_handles(model.num_verts(), handle_num, static_handles);
                std::vector<uint>* random_handles = get_random_handles(model, static_handles, gravity_center, ori_target_center);
                std::vector<cinolib::vec3d>* normals = get_positive_normals(model, gravity_center, ori_target_center, false);
                //                cout << "NOW STEP\t" << per_step << "\tNOW HANDLE:\t" << random_handles->at(0) << endl;
                for (int i = 0; i < static_handles->size(); i++) {
                    data.bcs[static_handles->at(i)] = all_vertexs[static_handles->at(i)];
                }
                //per_step = get_projection(ori_target_center - gravity_center).norm() * 0.2;
                for (int i = 0; i < random_handles->size(); i++) {
                    cinolib::vec3d normal_plane = unv::get_projection(normals->at(random_handles->at(i)));
                    data.bcs[random_handles->at(i)] = all_vertexs[random_handles->at(i)] + normal_plane * per_step;
                }
                clock_t s_per, f_per;
                s_per = clock();

                ARAP(model, data);
                gravity_center = get_gravity_center(model);
                deform_num++;
                if (Eva::SelfIntersection(model.vector_verts(), model.vector_polys()) || unv::get_projection(ori_target_center - gravity_center).norm() > last_dis) {
                    model = last_model;
                    for (int ti = 1; ti <= 5; ti++) {
                        //                        std::cout << "EEEOR\t" << ti << std::endl;
                        per_step /= 2;
                        for (int i = 0; i < random_handles->size(); i++) {
                            cinolib::vec3d normal_plane = unv::get_projection(normals->at(random_handles->at(i)));
                            data.bcs[random_handles->at(i)] = all_vertexs[random_handles->at(i)] + normal_plane * per_step;
                        }
                        ARAP(model, data);
                        cinolib::vec3d temp_gravity_center = get_gravity_center(model);
                        if (!(Eva::SelfIntersection(model.vector_verts(), model.vector_polys()) || unv::get_projection(ori_target_center - temp_gravity_center).norm() > last_dis)) {
                            break;
                        }
                        else {
                            model = last_model;
                        }
                    }
                }
                gravity_center = get_gravity_center(model);
                last_dis = unv::get_projection(ori_target_center - gravity_center).norm();
                double now_alpha = last_dis / unv::get_projection((ori_target_center - ori_gravity_center)).norm();
                temp_meshlist.push_back(std::make_pair(model, now_alpha));
                last_model = model;
                f_per = clock();
                //std::cout << "Per deformation:" << double(f_per - s_per) / CLOCKS_PER_SEC <<"target dis:\t"<< tar_dis<<"\t now dis:\t"<<(ori_target_center - gravity_center).norm() << std::endl;


                delete random_handles;
                delete normals;

                // std::cout << "Now dis:\t" << unv::get_projection(ori_target_center - gravity_center).norm() << "\t\tTarget dis:\t" << tar_dis << std::endl;
            }
            all_mesh.at(i) = temp_meshlist;
        }
    }
    finish = clock();
    for (uint i = 0; i < all_mesh.size(); i++) {
        for (uint j = 0; j < all_mesh.at(i).size(); j++) {
            meshlist.push_back(all_mesh.at(i).at(j));
        }
    }
    std::cout << "Deformation Time cost: " << (double)(finish - start) / CLOCKS_PER_SEC << std::endl;
    sort(meshlist.begin(), meshlist.end(), cmp);
    std::cout << "DataSet Size: " << meshlist.size() << std::endl;
    for (uint i = 0; i < meshlist.size(); i++) {
        cinolib::DrawableTrimesh<> model = meshlist[i].first;
        double alpha = meshlist[i].second;
        std::string str = std::to_string(alpha)+"_deform.obj";
        model.save(str.data());
    }
}
bool Deform::cmp(std::pair<cinolib::DrawableTrimesh<>, double> u, std::pair<cinolib::DrawableTrimesh<>, double> v) {
    return u.second > v.second;
}