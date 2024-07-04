#include<universal.h>

cinolib::vec3d target_center_g;

	cinolib::vec3d unv::get_projection(const cinolib::vec3d p) {
		return cinolib::vec3d(p[0], p[1], 0);
	}
	std::vector <std::vector<double>>* unv::normalize(std::vector <std::vector<double>>* f,double min_f, double max_f) {
		double max_t = f->at(0)[0];
		double min_t = f->at(0)[0];
		for (uint rid = 0; rid < f->size(); rid++) {
			for (uint cid = 0; cid < f->at(rid).size(); cid++) {
				max_t = std::max(max_t, f->at(rid)[cid]);
				min_t = std::min(min_t, f->at(rid)[cid]);
			}
		}
		for (uint rid = 0; rid < f->size(); rid++) {
			for (uint cid = 0; cid < f->at(rid).size(); cid++) {
				f->at(rid)[cid] = (f->at(rid)[cid] - min_t) / (max_t - min_t) * (max_f - min_f) + min_f;
			}
		}
		return f;
	}
	std::vector<double>* unv::normalize(std::vector<double>* f, double min_f, double max_f  ){
		double max_t = f->at(0);
		double min_t = f->at(0);
		for (uint rid = 0; rid < f->size(); rid++) {
			max_t = std::max(max_t, f->at(rid));
			min_t = std::min(min_t, f->at(rid));
		}
		for (uint rid = 0; rid < f->size(); rid++) {
			f->at(rid) = (f->at(rid) - min_t) / (max_t - min_t) * (max_f - min_f) + min_f;
			
		}
		return f;
	}
	std::vector<double>* unv::face2point(const std::vector<double> facedata, const cinolib::DrawableTrimesh<> m) {
		std::vector<std::vector<uint>> faces = m.vector_polys();
		std::vector<cinolib::vec3d> vertexs = m.vector_verts();
		std::vector<double> arealist(vertexs.size());
		std::vector<double>* pointdata = new std::vector<double>;
		for (uint i = 0; i < vertexs.size(); i++) {
			std::vector<uint> link = m.vert_ordered_polys_star(i);
			double area = 0;
			double score = 0;
			for (uint j = 0; j < link.size(); j++) {
				double temp_area = m.poly_area(link[j]);
				score += (temp_area * facedata.at(link[j]));
				area += temp_area;
			}
			score /= area;
			pointdata->push_back(score);
		}
		return pointdata;
	}
	std::vector<double>* unv::point2face(const std::vector<double> pointdata, const cinolib::DrawableTrimesh<> m) {
		std::vector<std::vector<uint>> faces = m.vector_polys();
		std::vector<double>* facedata = new std::vector<double>;
		for (uint i = 0; i < faces.size(); i++) {
			facedata->push_back((pointdata.at(faces[i][0]) + pointdata.at(faces[i][1]) + pointdata.at(faces[i][2])) / 3);
		}
		return facedata;
	}
