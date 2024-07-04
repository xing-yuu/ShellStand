#include "ShellThickening.h"


double Env::get_thickening_score(std::vector<double> RHO, std::vector<double>delta_m) {
	double res = 0;
	for (uint i = 0; i < delta_m.size(); i++) {
		res += (RHO[i] * delta_m[i]);
	}
	return res / delta_m.size();
}

std::pair<std::vector<double>, double> Env::shell_THICKENING(const cinolib::DrawableTrimesh<> m, std::pair<std::vector<double>, std::vector<double>> RHO_RANGE) {
	uint numfaces = m.num_polys();
	std::vector<double> RHO = RHO_RANGE.first;
	std::vector<double> range = RHO_RANGE.second;


	//TRAGET CENTER
	std::vector<uint>* base = Deform::get_static_handles(m);
	cinolib::vec3d target_center = Deform::get_target_center(m, base);


	//Min thickening
	std::vector<double> m0;
	for (uint i = 0; i < numfaces; i++) {
		m0.push_back(1.2);
	}
	std::vector<double>* thickening = msk::MOSEK_deltam(m0, RHO, target_center, range, 1.2, m);

	double score = get_thickening_score(RHO, *thickening);
	//facet to vertex

	//	delete saliency_score;
	//	delete visibility_score;
		//delete thickening_range;
		//delete face_range_a;
		//delete face_range_b;

	delete base;
	std::vector<double>* point_thickening = unv::face2point(*thickening, m);
	std::pair<std::vector<double>, double> ans(*point_thickening, score);
	delete point_thickening;
	delete thickening;
	return ans;

}

std::pair<std::vector<double>, double> Env::shell_THICKENING(const cinolib::DrawableTrimesh<> m) {
	uint numfaces = m.num_polys();
	std::vector<double>* RHO = new std::vector<double>;

	std::cout << "++++++++++++++++++++" << std::endl;
	std::cout << "Start calculate saliency." << std::endl;
	clock_t start = clock();
	//RHO
	std::vector<double>* saliency_score = Env::get_saliency_score(m, unv::FACE_DATA);
	clock_t end = clock();
	std::cout << "End saliency." << std::endl;
	std::cout << "Time cost: " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;
	std::cout << "++++++++++++++++++++" << std::endl << std::endl;



	std::cout << "++++++++++++++++++++" << std::endl;
	std::cout << "Start calculate visibility." << std::endl;
	clock_t start_2 = clock();
	std::vector <std::vector<double>>* visibility_score = Env::get_visibility_score(m, unv::FACE_DATA, 200, false);
	clock_t end_2 = clock();
	std::cout << "End visibility." << std::endl;
	std::cout << "Time cost: " << (double)(end_2 - start_2) / CLOCKS_PER_SEC << std::endl;
	std::cout << "++++++++++++++++++++" << std::endl << std::endl;

	for (uint i = 0; i < saliency_score->size(); i++) {
		double vi = (saliency_score->at(i) + std::max(visibility_score->at(0)[i], visibility_score->at(1)[i])) / 2;
		RHO->push_back(vi);
	}

	//RANGE
	std::cout << "++++++++++++++++++++" << std::endl;
	std::cout << "Start calculate upper bound." << std::endl;
	clock_t start_3 = clock();
	std::vector <std::vector<double>>* thickening_range = Env::get_thickening_range(m);
	std::vector<double>* face_range_a = unv::point2face(thickening_range->at(0), m);
	std::vector<double>* face_range_b = unv::point2face(thickening_range->at(1), m);
	std::vector<double>range;
	for (uint i = 0; i < face_range_b->size(); i++) {
		range.push_back(std::max(1.2, std::max(face_range_a->at(i), face_range_b->at(i))));
	}
	clock_t end_3 = clock();
	std::cout << "End upper bound.." << std::endl;
	std::cout << "Time cost: " << (double)(end_3 - start_3) / CLOCKS_PER_SEC << std::endl;
	std::cout << "++++++++++++++++++++" << std::endl << std::endl;



	//TRAGET CENTER
	std::vector<uint>* base = Deform::get_static_handles(m);
	cinolib::vec3d target_center = Deform::get_target_center(m, base);


	//Min thickening
	std::vector<double> m0;
	for (uint i = 0; i < numfaces; i++) {
		m0.push_back(1.2);
	}

	std::vector<double>* thickening = msk::MOSEK_deltam(m0, *RHO, target_center, range, 1.2, m);

	double score = get_thickening_score(*RHO, *thickening);
	double maxn = -1, minn = 1000000;
	for (uint pid = 0; pid < thickening->size(); pid++) {
		maxn = std::max(maxn, thickening->at(pid));
		minn = std::min(minn, thickening->at(pid));
		thickening->at(pid) += m0.at(pid);
	}
	std::cout << "MAX RANGE: " << maxn << std::endl;
	std::cout << "MIN RANGE: " << minn << std::endl;
	//	delete saliency_score;
	//	delete visibility_score;
		//delete thickening_range;
		//delete face_range_a;
		//delete face_range_b;
	delete RHO;
	delete base;
	
	std::vector<double> *point_thickening = unv::face2point(*thickening,m);
	std::pair<std::vector<double>, double> ans(*point_thickening, score);
	delete point_thickening;
	delete thickening;
	return ans;

}


void Env::output_volume_model(const std::string output, const cinolib::DrawableTrimesh<> m, const std::vector<double> offset) {
	std::vector<std::vector<uint>> faces_a = m.vector_polys();
	std::vector<std::vector<uint>> faces_b = m.vector_polys();
	std::vector<cinolib::vec3d> vertexs_a = m.vector_verts();
	std::vector<cinolib::vec3d> vertexs_b = m.vector_verts();
	std::vector<cinolib::ipair> boundary = m.get_boundary_edges();
	std::vector<cinolib::vec3d> normals = m.vector_vert_normals();
	uint num_vertex = m.num_verts();
	uint num_face = m.num_polys();
	for (uint vid = 0; vid < num_vertex; vid++) {
		vertexs_a[vid] = vertexs_a[vid] - offset[vid] * 0.5 * normals[vid];
		vertexs_b[vid] = vertexs_b[vid] + offset[vid] * 0.5 * normals[vid];
	}
	for (uint pid = 0; pid < num_face; pid++) {
		for (uint vid = 0; vid < faces_b.at(pid).size(); vid++) {
			faces_b[pid][vid] += num_vertex;
		}
	}

	for (uint vid = 0; vid < boundary.size(); vid++) {
		cinolib::ipair edge_i = boundary[vid];
		uint vid0 = edge_i.first;
		uint vid1 = edge_i.second;
		uint vid2 = vid0 + num_vertex;
		uint vid3 = vid1 + num_vertex;
		std::vector<uint> face0 = { vid0 ,vid1,vid2 };
		std::vector<uint> face1 = { vid2 ,vid1,vid3 };
		faces_b.push_back(face0);
		faces_b.push_back(face1);
	}
	faces_a.insert(faces_a.end(), faces_b.begin(), faces_b.end());
	vertexs_a.insert(vertexs_a.end(), vertexs_b.begin(), vertexs_b.end());
	cinolib::DrawableTrimesh<> ans(vertexs_a, faces_a);


	ans = Ope::add_base(ans, target_center_g);

	ans.save(output.data());
}

