#include "ThickeningRange.h"


void Env::output_range(const cinolib::DrawableTrimesh<> m, std::vector <std::vector<double>> thickness,std::string debug_path, std::string model_name) {
    std::vector<std::vector<uint>> faces = m.vector_polys();
    std::vector<cinolib::vec3d> vertexs_a = m.vector_verts();
    std::vector<cinolib::vec3d> vertexs_b = m.vector_verts();
    std::vector<cinolib::vec3d> normals = m.vector_vert_normals();
    for (uint vid = 0; vid < vertexs_a.size(); vid++) {
        vertexs_a[vid] += thickness[0][vid] * normals[vid];
        vertexs_b[vid] -= thickness[1][vid] * normals[vid];
    }
    cinolib::DrawableTrimesh<> m_a(vertexs_a, faces);
    cinolib::DrawableTrimesh<> m_b(vertexs_b, faces);
    m_a += m_b;
    std::string save_path = debug_path + "\\" + model_name + "_boundary_test.obj";
    m_a.save(save_path.data());
}

void Env::normal_smooth(const cinolib::DrawableTrimesh<> m, std::vector<cinolib::vec3d>& normals, uint it) {
    uint numvertexs = m.num_verts();
    for (uint k = 0; k < it; k++) {
        std::vector<cinolib::vec3d> temp_normals;
        for (uint vid = 0; vid < numvertexs; vid++) {
            std::vector<uint> link = m.vert_verts_link(vid);
            cinolib::vec3d avg_normal = normals[vid];
            for (uint lid = 0; lid < link.size(); lid++) {
                avg_normal += normals[link[lid]];
            }
            avg_normal /= (link.size() + 1);
            temp_normals.push_back(avg_normal);
        }
        normals = temp_normals;
    }
}
std::vector<double> Env::thickening(const cinolib::DrawableTrimesh<> m, const std::vector<cinolib::vec3d> normals, int type) {
    std::vector<std::vector<uint>> faces = m.vector_polys();
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    std::vector<double> a_side;
    std::vector<double> last_side;
    std::vector<double> basic_step_size;
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

    double step_size = sqrt((xMax - xMin) * (xMax - xMin) + (yMax - yMin) * (yMax - yMin) + (zMax - zMin) * (zMax - zMin)) * 0.02;
    double max_size = sqrt((xMax - xMin) * (xMax - xMin) + (yMax - yMin) * (yMax - yMin) + (zMax - zMin) * (zMax - zMin)) * 0.4;
    for (uint i = 0; i < vertexs.size(); i++) {
        basic_step_size.push_back(step_size);
        a_side.push_back(0);
        last_side.push_back(0);
    }
    uint idnum = 10;


    //printf("GL0BAL %d\n", intersection_facets_list_ori->size());

    for (uint numstp = 0; numstp < idnum; numstp++) {
        printf("\r 1/2 [%.2f%%]\t>", float(numstp * 100.0 / (idnum - 1)));
        for (int j = 1; j <= numstp * 40 / idnum; j++)
            std::cout << "¨";
        std::vector<cinolib::vec3d> test_vertexs;
        for (uint i = 0; i < vertexs.size(); i++) {
            basic_step_size[i] = step_size;
            test_vertexs.push_back(vertexs[i] + double(type) * normals[i] * (a_side[i] + basic_step_size[i]));
        }
        std::set<uint>* intersection_facets_list = Eva::IntersectionList(test_vertexs, faces);//facets
        while (!intersection_facets_list->empty()) {
            //printf("Size of intersection_facets_list: %d\n", intersection_facets_list->size());
            std::set<uint> intersection_vertexs_list;
            for (std::set<uint>::iterator i = intersection_facets_list->begin(); i != intersection_facets_list->end(); i++) {
                std::vector<uint> fi = faces[*i];
                intersection_vertexs_list.insert(fi[0]);
                intersection_vertexs_list.insert(fi[1]);
                intersection_vertexs_list.insert(fi[2]);
            }
            for (std::set<uint>::iterator i = intersection_vertexs_list.begin(); i != intersection_vertexs_list.end(); i++) {
                basic_step_size[*i] /= 2;
            }
            for (uint i = 0; i < vertexs.size(); i++) {
                test_vertexs[i]=(vertexs[i] + double(type) * normals[i] * (a_side[i] + basic_step_size[i]));
            }
            delete intersection_facets_list;
            intersection_facets_list = Eva::IntersectionList(test_vertexs, faces);
        }
        
        for (uint i = 0; i < vertexs.size(); i++) {
            a_side[i] = a_side[i] + basic_step_size[i];
        }
    }
    printf("\n");
    return a_side;
}

std::vector <std::vector<double>>* Env::get_thickening_range(const cinolib::DrawableTrimesh<> m) {

    
    std::vector<cinolib::vec3d> normals = m.vector_vert_normals();
    std::vector <std::vector<double>>* two_side_range = new std::vector <std::vector<double>>;
    
    Env::normal_smooth(m, normals);
    two_side_range->push_back(Env::thickening(m, normals, 1));
    two_side_range->push_back(Env::thickening(m, normals, -1));
    return two_side_range;
}