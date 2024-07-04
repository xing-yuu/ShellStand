#include "ModelOperation.h"
cinolib::DrawableTrimesh<> Ope::scale(const cinolib::DrawableTrimesh<> m, double scale_size ) {
    std::vector<std::vector<uint>> faces = m.vector_polys();
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();
    std::vector<cinolib::vec3d> vertexs_new;
    double min_x = vertexs[0][0], max_x = vertexs[0][0];
    for (uint vid = 0; vid < m.num_verts(); vid++) {
        max_x = std::max(max_x, vertexs[vid][0]);
        min_x = std::min(min_x, vertexs[vid][0]);
    }
    double factor = scale_size / (max_x - min_x);
    for (uint vid = 0; vid < m.num_verts(); vid++) {
        vertexs_new.push_back(vertexs[vid] * factor);
    }
    return cinolib::DrawableTrimesh<>(vertexs_new, faces);
}

cinolib::DrawableTrimesh<> Ope::add_base(const cinolib::DrawableTrimesh<> m, cinolib::vec3d center, double radius ) {
    std::vector<cinolib::vec3d> vertexs = m.vector_verts();

    std::vector< std::vector<uint>> facet_a, facet_b;
    std::vector<cinolib::vec3d> ring_a, ring_b;
    int sampling = 20;
    double height = 2.5;
    
    ring_a.push_back(cinolib::vec3d(center[0], center[1], center[2] - height / 2));
    ring_b.push_back(cinolib::vec3d(center[0], center[1], center[2] + height / 2));
    
    for (int i = 0; i < sampling; i++) {
        cinolib::vec3d vi(cos(double(i * 1.0 / sampling * 2 * M_PI)), sin(double(i * 1.0 / sampling * 2 * M_PI)), 0);
        vi = vi * radius;
        vi += center;
        vi[2] -= height / 2;
        ring_a.push_back(vi);
        vi[2] += height;
        ring_b.push_back(vi);
        std::vector<uint> fi;
        
        fi.push_back(0);
        fi.push_back(i+1);
        if (i == sampling - 1) 
            fi.push_back(1);
        
        else 
            fi.push_back(i + 2);
        facet_a.push_back(fi);
        
    }
    for (int i = 0; i < sampling; i++) {
        std::vector<uint> fi = facet_a[i];
        for (int j = 0; j < fi.size(); j++) {
            fi[j] += sampling + 1;
        }
        facet_b.push_back(fi);
    }
    for (int i = 1; i <= sampling; i++) {
        uint v0 = i, v1, v2, v3;
        if (i == sampling) {
            v1 = 1;
        }
        else {
            v1 = i + 1;
        }
        v2 = v0 + sampling + 1;
        v3 = v1 + sampling + 1;
        std::vector<uint> fa = { v0,v1,v2 };
        std::vector<uint> fb = { v2,v1,v3 };
        facet_b.push_back(fa);
        facet_b.push_back(fb);
    }
    facet_a.insert(facet_a.end(), facet_b.begin(), facet_b.end());
    ring_a.insert(ring_a.end(), ring_b.begin(), ring_b.end());
    cinolib::DrawableTrimesh<> ans(ring_a, facet_a);
    ans += m;
    return ans;
}