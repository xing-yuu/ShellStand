#include "SelfIntersect.h"
bool Eva::SelfIntersection(const std::vector<cinolib::vec3d> vertexs, const std::vector<std::vector<uint>> faces) {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Surface_mesh<K::Point_3>             Mesh;
    typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
    typedef Mesh::Vertex_index vertex_descriptor;
    namespace PMP = CGAL::Polygon_mesh_processing;
    std::vector<vertex_descriptor> face_cgal_list;
    Mesh mesh;
    for (uint i = 0; i < vertexs.size(); i++) {
        vertex_descriptor u = mesh.add_vertex(K::Point_3(vertexs[i][0], vertexs[i][1], vertexs[i][2]));
        face_cgal_list.push_back(u);
    }
    for (uint i = 0; i < faces.size(); i++) {
        mesh.add_face(face_cgal_list[faces[i][0]], face_cgal_list[faces[i][1]], face_cgal_list[faces[i][2]]);
    }
    bool intersecting = PMP::does_self_intersect(mesh,
        PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));
    return intersecting;
}

std::set<uint>* Eva::IntersectionList(const std::vector<cinolib::vec3d> vertexs, const std::vector<std::vector<uint>> faces) {
    typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
    typedef CGAL::Surface_mesh<K::Point_3>             Mesh;
    typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
    typedef Mesh::Vertex_index vertex_descriptor;
    namespace PMP = CGAL::Polygon_mesh_processing;
    std::vector<vertex_descriptor> face_cgal_list;
    std::set<uint>* intersection_list = new std::set<uint>;
    Mesh mesh;
    double max_a = 1, min_a = 100000;
    for (uint i = 0; i < vertexs.size(); i++) {
        vertex_descriptor u = mesh.add_vertex(K::Point_3(vertexs[i][0], vertexs[i][1], vertexs[i][2]));
        max_a = std::max(max_a, std::max(std::max(vertexs[i][0], vertexs[i][1]), vertexs[i][2]));
        min_a = std::min(min_a, std::min(std::min(vertexs[i][0], vertexs[i][1]), vertexs[i][2]));
        face_cgal_list.push_back(u);
    }
    //cout << "MAX " << max_a << "MIN " << min_a << endl;
    for (uint i = 0; i < faces.size(); i++) {
        mesh.add_face(face_cgal_list[faces[i][0]], face_cgal_list[faces[i][1]], face_cgal_list[faces[i][2]]);
    }
    
 /*   bool intersecting = PMP::does_self_intersect(mesh,
        PMP::parameters::vertex_point_map(get(CGAL::vertex_point, mesh)));
    if (!intersecting) return intersection_list;
*/
    std::vector<std::pair<face_descriptor, face_descriptor> > intersected_tris;
    //std::cout << "TESTBv" << std::endl;
    PMP::self_intersections(mesh, back_inserter(intersected_tris));
    //std::cout << "TESTB" << std::endl;
    //cout << intersected_tris.size() << "对三角形相交" << endl;
    std::vector<std::pair<face_descriptor, face_descriptor> >::iterator iter;
    for (iter = intersected_tris.begin(); iter != intersected_tris.end(); iter++) {
        intersection_list->insert(iter->first);
        intersection_list->insert(iter->second);
    }
    //std::cout << "TESTC" << std::endl;
    return intersection_list;
}