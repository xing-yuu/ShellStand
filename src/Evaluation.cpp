#include <Evaluation.h>




cinolib::vec3d Eva::rotate_p(cinolib::vec3d p, double angle) {
	cinolib::vec3d line0(cos(angle), sin(angle), 0);
	cinolib::vec3d line1(-sin(angle), cos(angle), 0);
	cinolib::vec3d line2(0, 0, 1);
	return cinolib::vec3d(line0.dot(p), line1.dot(p), line2.dot(p));
}
cinolib::vec2d Eva::map_coordinate(cinolib::vec2d p, cinolib::vec2d p_min, cinolib::vec2d p_max) {
	p -= p_min;
	cinolib::vec2d p0 = p_max - p_min;
	return cinolib::vec2d(p[0] / p0[0], p[1] / p0[1]);
}
double Eva::compute_area(const cinolib::vec2d v1, const cinolib::vec2d v2, const cinolib::vec2d v3) {
	return 0.5 * fabs(v1[0] * (v2[1] - v3[1]) + v2[0] * (v3[1] - v1[1]) + v3[0] * (v1[1] - v2[1]));
}
bool Eva::is_point_in_triangle1(const cinolib::vec2d v0, const cinolib::vec2d v1, const cinolib::vec2d v2, const cinolib::vec2d v3) {
	double area_ABC = compute_area(v1, v2, v3);
	double area_PAB = compute_area(v0, v2, v3);
	double area_PAC = compute_area(v0, v1, v3);
	double area_PBC = compute_area(v0, v1, v2);

	if (fabs(area_PAB + area_PBC + area_PAC - area_ABC) < 0.000001)
		return true;
	else return false;
}
double Eva::evaluate_similarity(const cinolib::DrawableTrimesh<> m0, const cinolib::DrawableTrimesh<> m, uint evaluate_angle_num) {
	const int resolution = 256;
	double sumn = 0;
	std::vector<cinolib::vec3d> vertexs0 = m0.vector_verts();
	std::vector<cinolib::vec3d> vertexs1 = m.vector_verts();
	std::vector<std::vector<uint>> faces = m0.vector_polys();

	cinolib::vec3d rotation_axis(0, 0, 0), max_p(vertexs0[0]), min_p(vertexs0[0]);
	for (int i = 0; i < vertexs0.size(); i++) {
		max_p[0] = std::max(vertexs0[i][0], max_p[0]);
		max_p[1] = std::max(vertexs0[i][1], max_p[1]);
		max_p[2] = std::min(vertexs0[i][2], max_p[2]);
		min_p[0] = std::min(vertexs0[i][0], min_p[0]);
		min_p[1] = std::min(vertexs0[i][1], min_p[1]);
		min_p[2] = std::min(vertexs0[i][2], min_p[2]);
	}

	for (int i = 0; i < vertexs0.size(); i++) {
		vertexs0[i] -= rotation_axis;
		vertexs1[i] -= rotation_axis;
	}
	rotation_axis = (max_p + min_p) / 2;
	for (uint i = 0; i < evaluate_angle_num; i++)
	{
		bool in0[resolution][resolution];
		bool in1[resolution][resolution];
		memset(in0, 0, sizeof(in0));
		memset(in1, 0, sizeof(in1));
		std::vector<cinolib::vec3d> v0(vertexs0);
		std::vector<cinolib::vec3d> v1(vertexs1);
		cinolib::vec2d max_2d(INT_MIN, INT_MIN), min_2d(INT_MAX, INT_MAX);
		double rotate_angle = double(i) / double(evaluate_angle_num) * 2.0 * M_PI;
		for (uint j = 0; j < vertexs0.size(); j++) {
			vertexs0[j] = rotate_p(vertexs0[j], rotate_angle);
			vertexs1[j] = rotate_p(vertexs1[j], rotate_angle);
			max_2d[0] = std::max(max_2d[0], vertexs0[j][0]);
			max_2d[0] = std::max(max_2d[0], vertexs1[j][0]);
			max_2d[1] = std::max(max_2d[1], vertexs0[j][2]);
			max_2d[1] = std::max(max_2d[1], vertexs1[j][2]);
			min_2d[0] = std::min(min_2d[0], vertexs0[j][0]);
			min_2d[0] = std::min(min_2d[0], vertexs1[j][0]);
			min_2d[1] = std::min(min_2d[1], vertexs0[j][2]);
			min_2d[1] = std::min(min_2d[1], vertexs1[j][2]);
		}

		for (uint j = 0; j < faces.size(); j++) {
			cinolib::vec2d v1 = map_coordinate(cinolib::vec2d(vertexs0[faces[j][0]][0], vertexs0[faces[j][0]][2]), min_2d, max_2d) * (resolution - 1);
			cinolib::vec2d v2 = map_coordinate(cinolib::vec2d(vertexs0[faces[j][1]][0], vertexs0[faces[j][1]][2]), min_2d, max_2d) * (resolution - 1);
			cinolib::vec2d v3 = map_coordinate(cinolib::vec2d(vertexs0[faces[j][2]][0], vertexs0[faces[j][2]][2]), min_2d, max_2d) * (resolution - 1);
			int minx = std::max(0, int(std::min(std::min(v1[0], v2[0]), v3[0]) * 1.0));
			int miny = std::max(0, int(std::min(std::min(v1[1], v2[1]), v3[1]) * 1.0));
			int maxx = std::min((resolution - 1), int(std::max(std::max(v1[0], v2[0]), v3[0]) * 1.0));
			int maxy = std::min((resolution - 1), int(std::max(std::max(v1[1], v2[1]), v3[1]) * 1.0));
			for (int k = minx; k <= maxx; k++) {
				for (int l = miny; l <= maxy; l++)
					if (is_point_in_triangle1(cinolib::vec2d(k, l), v1, v2, v3))
						in0[k][l] = true;
			}
		}
		for (uint j = 0; j < faces.size(); j++) {
			cinolib::vec2d v1 = map_coordinate(cinolib::vec2d(vertexs1[faces[j][0]][0], vertexs1[faces[j][0]][2]), min_2d, max_2d) * (resolution - 1);
			cinolib::vec2d v2 = map_coordinate(cinolib::vec2d(vertexs1[faces[j][1]][0], vertexs1[faces[j][1]][2]), min_2d, max_2d) * (resolution - 1);
			cinolib::vec2d v3 = map_coordinate(cinolib::vec2d(vertexs1[faces[j][2]][0], vertexs1[faces[j][2]][2]), min_2d, max_2d) * (resolution - 1);
			int minx = std::max(0, int(std::min(std::min(v1[0], v2[0]), v3[0]) * 1.0));
			int miny = std::max(0, int(std::min(std::min(v1[1], v2[1]), v3[1]) * 1.0));
			int maxx = std::min((resolution - 1), int(std::max(std::max(v1[0], v2[0]), v3[0]) * 1.0));
			int maxy = std::min((resolution - 1), int(std::max(std::max(v1[1], v2[1]), v3[1]) * 1.0));
			for (int k = minx; k <= maxx; k++) {
				for (int l = miny; l <= maxy; l++)
					if (is_point_in_triangle1(cinolib::vec2d(k, l), v1, v2, v3))
						in1[k][l] = true;
			}
		}
		for (uint j = 0; j < resolution; j++)
			for (uint k = 0; k < resolution; k++)
				if (in0[j][k] ^ in1[j][k])
					sumn += 1;

	}
	sumn = sumn / (evaluate_angle_num * resolution * resolution);
	return sumn;
}