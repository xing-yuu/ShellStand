/*
 * author: Yupan Liu
 * date: Dec 27, 2015
 * brief: vector and matrix library
 * comment: modified by "Anton's OpenGL 4 Tutorial"
 * */
#include "maths_funcs.h"
//using namespace MFUN;

/*--------------------------------CONSTRUCTORS--------------------------------*/
MFUN::vec2::vec2 () {}

MFUN::vec2::vec2 (float x, float y) {
	v[0] = x;
	v[1] = y;
}

MFUN::vec3::vec3 () {}

MFUN::vec3::vec3 (float x, float y, float z) {
	v[0] = x;
	v[1] = y;
	v[2] = z;
}

MFUN::vec3::vec3 (const MFUN::vec2& vv, float z) {
	v[0] = vv.v[0];
	v[1] = vv.v[1];
	v[2] = z;
}

MFUN::vec3::vec3 (const MFUN::vec4& vv) {
	v[0] = vv.v[0];
	v[1] = vv.v[1];
	v[2] = vv.v[2];
}

MFUN::vec4::vec4 () {}

MFUN::vec4::vec4 (float x, float y, float z, float w) {
	v[0] = x;
	v[1] = y;
	v[2] = z;
	v[3] = w;
}

MFUN::vec4::vec4 (const MFUN::vec2& vv, float z, float w) {
	v[0] = vv.v[0];
	v[1] = vv.v[1];
	v[2] = z;
	v[3] = w;
}

MFUN::vec4::vec4 (const MFUN::vec3& vv, float w) {
	v[0] = vv.v[0];
	v[1] = vv.v[1];
	v[2] = vv.v[2];
	v[3] = w;
}

MFUN::mat3::mat3 () {}

/* note: entered in COLUMNS */
MFUN::mat3::mat3 (float a, float b, float c,
						float d, float e, float f,
						float g, float h, float i) {
	m[0] = a;
	m[1] = b;
	m[2] = c;
	m[3] = d;
	m[4] = e;
	m[5] = f;
	m[6] = g;
	m[7] = h;
	m[8] = i;
}

MFUN::mat4::mat4 () {}

/* note: entered in COLUMNS */
MFUN::mat4::mat4 (float a, float b, float c, float d,
						float e, float f, float g, float h,
						float i, float j, float k, float l,
						float mm, float n, float o, float p) {
	m[0] = a;
	m[1] = b;
	m[2] = c;
	m[3] = d;
	m[4] = e;
	m[5] = f;
	m[6] = g;
	m[7] = h;
	m[8] = i;
	m[9] = j;
	m[10] = k;
	m[11] = l;
	m[12] = mm;
	m[13] = n;
	m[14] = o;
	m[15] = p;
}

/*-----------------------------PRINT FUNCTIONS--------------------------------*/
void MFUN::print (const MFUN::vec2& v) {
	printf ("[%.2f, %.2f]\n", v.v[0], v.v[1]);
}

void MFUN::print (const MFUN::vec3& v) {
	printf ("[%.2f, %.2f, %.2f]\n", v.v[0], v.v[1], v.v[2]);
}

void MFUN::print (const MFUN::vec4& v) {
	printf ("[%.2f, %.2f, %.2f, %.2f]\n", v.v[0], v.v[1], v.v[2], v.v[3]);
}

void MFUN::print (const MFUN::mat3& m) {
	printf("\n");
	printf ("[%.2f][%.2f][%.2f]\n", m.m[0], m.m[3], m.m[6]);
	printf ("[%.2f][%.2f][%.2f]\n", m.m[1], m.m[4], m.m[7]);
	printf ("[%.2f][%.2f][%.2f]\n", m.m[2], m.m[5], m.m[8]);
}

void MFUN::print (const MFUN::mat4& m) {
	printf("\n");
	printf ("[%.2f][%.2f][%.2f][%.2f]\n", m.m[0], m.m[4], m.m[8], m.m[12]);
	printf ("[%.2f][%.2f][%.2f][%.2f]\n", m.m[1], m.m[5], m.m[9], m.m[13]);
	printf ("[%.2f][%.2f][%.2f][%.2f]\n", m.m[2], m.m[6], m.m[10], m.m[14]);
	printf ("[%.2f][%.2f][%.2f][%.2f]\n", m.m[3], m.m[7], m.m[11], m.m[15]);
}

/*------------------------------VECTOR FUNCTIONS------------------------------*/
float MFUN::length (const MFUN::vec3& v) {
	return sqrt (v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2]);
}

// squared length
float MFUN::length2 (const MFUN::vec3& v) {
	return v.v[0] * v.v[0] + v.v[1] * v.v[1] + v.v[2] * v.v[2];
}

// note: proper spelling (hehe)
MFUN::vec3 MFUN::normalise (const MFUN::vec3& v) {
	MFUN::vec3 vb;
	float l = MFUN::length (v);
	if (0.0f == l) {
		return MFUN::vec3 (0.0f, 0.0f, 0.0f);
	}
	vb.v[0] = v.v[0] / l;
	vb.v[1] = v.v[1] / l;
	vb.v[2] = v.v[2] / l;
	return vb;
}

MFUN::vec3 MFUN::vec3::operator+ (const MFUN::vec3& rhs) {
	MFUN::vec3 vc;
	vc.v[0] = v[0] + rhs.v[0];
	vc.v[1] = v[1] + rhs.v[1];
	vc.v[2] = v[2] + rhs.v[2];
	return vc;
}

MFUN::vec3& MFUN::vec3::operator+= (const MFUN::vec3& rhs) {
	v[0] += rhs.v[0];
	v[1] += rhs.v[1];
	v[2] += rhs.v[2];
	return *this; // return self
}

MFUN::vec3 MFUN::vec3::operator- (const MFUN::vec3& rhs) {
	MFUN::vec3 vc;
	vc.v[0] = v[0] - rhs.v[0];
	vc.v[1] = v[1] - rhs.v[1];
	vc.v[2] = v[2] - rhs.v[2];
	return vc;
}

MFUN::vec3& MFUN::vec3::operator-= (const MFUN::vec3& rhs) {
	v[0] -= rhs.v[0];
	v[1] -= rhs.v[1];
	v[2] -= rhs.v[2];
	return *this;
}

MFUN::vec3 MFUN::vec3::operator+ (float rhs) {
	MFUN::vec3 vc;
	vc.v[0] = v[0] + rhs;
	vc.v[1] = v[1] + rhs;
	vc.v[2] = v[2] + rhs;
	return vc;
}

MFUN::vec3 MFUN::vec3::operator- (float rhs) {
	MFUN::vec3 vc;
	vc.v[0] = v[0] - rhs;
	vc.v[1] = v[1] - rhs;
	vc.v[2] = v[2] - rhs;
	return vc;
}

MFUN::vec3 MFUN::vec3::operator* (float rhs) {
	MFUN::vec3 vc;
	vc.v[0] = v[0] * rhs;
	vc.v[1] = v[1] * rhs;
	vc.v[2] = v[2] * rhs;
	return vc;
}

MFUN::vec3 MFUN::vec3::operator/ (float rhs) {
	MFUN::vec3 vc;
	vc.v[0] = v[0] / rhs;
	vc.v[1] = v[1] / rhs;
	vc.v[2] = v[2] / rhs;
	return vc;
}

MFUN::vec3& MFUN::vec3::operator*= (float rhs) {
	v[0] = v[0] * rhs;
	v[1] = v[1] * rhs;
	v[2] = v[2] * rhs;
	return *this;
}

MFUN::vec3& MFUN::vec3::operator= (const MFUN::vec3& rhs) {
	v[0] = rhs.v[0];
	v[1] = rhs.v[1];
	v[2] = rhs.v[2];
	return *this;
}

float MFUN::dot (const MFUN::vec3& a, const MFUN::vec3& b) {
	return a.v[0] * b.v[0] + a.v[1] * b.v[1] + a.v[2] * b.v[2];
}

MFUN::vec3 MFUN::cross (const MFUN::vec3& a, const MFUN::vec3& b) {
	float x = a.v[1] * b.v[2] - a.v[2] * b.v[1];
	float y = a.v[2] * b.v[0] - a.v[0] * b.v[2];
	float z = a.v[0] * b.v[1] - a.v[1] * b.v[0];
	return MFUN::vec3 (x, y, z);
}

//new
MFUN::mat3 MFUN::wedge (const MFUN::vec3& a, const MFUN::vec3& b) {
	MFUN::mat3 r = MFUN::zero_mat3();
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++)
			r.m[3*i+j] = a.v[i] * b.v[j];
	return r;
}

float MFUN::get_squared_dist (MFUN::vec3 from, MFUN::vec3 to) {
	float x = (to.v[0] - from.v[0]) * (to.v[0] - from.v[0]);
	float y = (to.v[1] - from.v[1]) * (to.v[1] - from.v[1]);
	float z = (to.v[2] - from.v[2]) * (to.v[2] - from.v[2]);
	return x + y + z;
}

/* converts an un-normalised direction into a heading in degrees
NB i suspect that the z is backwards here but i've used in in
several places like this. d'oh! */
float MFUN::direction_to_heading (MFUN::vec3 d) {
	return atan2 (-d.v[0], -d.v[2]) * ONE_RAD_IN_DEG;
}

MFUN::vec3 MFUN::heading_to_direction (float degrees) {
	float rad = degrees * ONE_DEG_IN_RAD;
	return MFUN::vec3 (-sinf (rad), 0.0f, -cosf (rad));
}

/*-----------------------------MATRIX FUNCTIONS-------------------------------*/
MFUN::mat3 MFUN::zero_mat3 () {
	return MFUN::mat3 (
		0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f
	);
}

MFUN::mat3 MFUN::identity_mat3 () {
	return MFUN::mat3 (
		1.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 1.0f
	);
}

MFUN::mat4 MFUN::zero_mat4 () {
	return MFUN::mat4 (
		0.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 0.0f
	);
}

MFUN::mat4 MFUN::identity_mat4 () {
	return MFUN::mat4 (
		1.0f, 0.0f, 0.0f, 0.0f,
		0.0f, 1.0f, 0.0f, 0.0f,
		0.0f, 0.0f, 1.0f, 0.0f,
		0.0f, 0.0f, 0.0f, 1.0f
	);
}

/* mat4 array layout
 0  4  8 12
 1  5  9 13
 2  6 10 14
 3  7 11 15
*/

//new
MFUN::vec3 MFUN::mat3::operator* (const MFUN::vec3& rhs) {
	// 0x + 3y + 6z
	float x = m[0] * rhs.v[0] +
		m[3] * rhs.v[1] +
		m[6] * rhs.v[2];
	// 1x + 4y + 7z
	float y = m[1] * rhs.v[0] +
		m[4] * rhs.v[1] +
		m[7] * rhs.v[2];
	// 2x + 5y + 8z
	float z = m[2] * rhs.v[0] +
		m[5] * rhs.v[1] +
		m[8] * rhs.v[2];
	return MFUN::vec3 (x, y, z);
}

MFUN::vec4 MFUN::mat4::operator* (const MFUN::vec4& rhs) {
	// 0x + 4y + 8z + 12w
	float x =
		m[0] * rhs.v[0] +
		m[4] * rhs.v[1] +
		m[8] * rhs.v[2] +
		m[12] * rhs.v[3];
	// 1x + 5y + 9z + 13w
	float y = m[1] * rhs.v[0] +
		m[5] * rhs.v[1] +
		m[9] * rhs.v[2] +
		m[13] * rhs.v[3];
	// 2x + 6y + 10z + 14w
	float z = m[2] * rhs.v[0] +
		m[6] * rhs.v[1] +
		m[10] * rhs.v[2] +
		m[14] * rhs.v[3];
	// 3x + 7y + 11z + 15w
	float w = m[3] * rhs.v[0] +
		m[7] * rhs.v[1] +
		m[11] * rhs.v[2] +
		m[15] * rhs.v[3];
	return MFUN::vec4 (x, y, z, w);
}

//new
MFUN::mat3 MFUN::mat3::operator* (const MFUN::mat3& rhs) {
	MFUN::mat3 r = MFUN::zero_mat3 ();
	int r_index = 0;
	for (int col = 0; col < 3; col++) {
		for (int row = 0; row < 3; row++) {
			float sum = 0.0f;
			for (int i = 0; i < 3; i++) {
				sum += rhs.m[i + col * 3] * m[row + i * 3];
			}
			r.m[r_index] = sum;
			r_index++;
		}
	}
	return r;
}

MFUN::mat4 MFUN::mat4::operator* (const MFUN::mat4& rhs) {
	MFUN::mat4 r = MFUN::zero_mat4 ();
	int r_index = 0;
	for (int col = 0; col < 4; col++) {
		for (int row = 0; row < 4; row++) {
			float sum = 0.0f;
			for (int i = 0; i < 4; i++) {
				sum += rhs.m[i + col * 4] * m[row + i * 4];
			}
			r.m[r_index] = sum;
			r_index++;
		}
	}
	return r;
}

//new
MFUN::mat3 MFUN::mat3::operator* (const float& rhs) {
	MFUN::mat3 r = MFUN::zero_mat3();
	for(int i = 0; i < 9; i++)
		r.m[i] = rhs * m[i];
	return r;
}

//new
MFUN::mat3 MFUN::mat3::operator+ (const MFUN::mat3& rhs) {
	MFUN::mat3 r = MFUN::zero_mat3 ();
	for(int i = 0; i < 9; i++)
		r.m[i] = m[i] + rhs.m[i];
	return r;
}

//new
MFUN::mat3 MFUN::mat3::operator- (const MFUN::mat3& rhs) {
	MFUN::mat3 r = MFUN::zero_mat3 ();
	for(int i = 0; i < 9; i++)
		r.m[i] = m[i] - rhs.m[i];
	return r;
}

//new
MFUN::mat3& MFUN::mat3::operator= (const MFUN::mat3& rhs) {
	for (int i = 0; i < 9; i++) {
		m[i] = rhs.m[i];
	}
	return *this;
}

MFUN::mat4& MFUN::mat4::operator= (const MFUN::mat4& rhs) {
	for (int i = 0; i < 16; i++) {
		m[i] = rhs.m[i];
	}
	return *this;
}

// returns a scalar value with the determinant for a 4x4 matrix
// see http://www.euclideanspace.com/maths/algebra/matrix/functions/determinant/fourD/index.htm
float MFUN::determinant (const MFUN::mat4& mm) {
	return
		mm.m[12] * mm.m[9] * mm.m[6] * mm.m[3] -
		mm.m[8] * mm.m[13] * mm.m[6] * mm.m[3] -
		mm.m[12] * mm.m[5] * mm.m[10] * mm.m[3] +
		mm.m[4] * mm.m[13] * mm.m[10] * mm.m[3] +
		mm.m[8] * mm.m[5] * mm.m[14] * mm.m[3] -
		mm.m[4] * mm.m[9] * mm.m[14] * mm.m[3] -
		mm.m[12] * mm.m[9] * mm.m[2] * mm.m[7] +
		mm.m[8] * mm.m[13] * mm.m[2] * mm.m[7] +
		mm.m[12] * mm.m[1] * mm.m[10] * mm.m[7] -
		mm.m[0] * mm.m[13] * mm.m[10] * mm.m[7] -
		mm.m[8] * mm.m[1] * mm.m[14] * mm.m[7] +
		mm.m[0] * mm.m[9] * mm.m[14] * mm.m[7] +
		mm.m[12] * mm.m[5] * mm.m[2] * mm.m[11] -
		mm.m[4] * mm.m[13] * mm.m[2] * mm.m[11] -
		mm.m[12] * mm.m[1] * mm.m[6] * mm.m[11] +
		mm.m[0] * mm.m[13] * mm.m[6] * mm.m[11] +
		mm.m[4] * mm.m[1] * mm.m[14] * mm.m[11] -
		mm.m[0] * mm.m[5] * mm.m[14] * mm.m[11] -
		mm.m[8] * mm.m[5] * mm.m[2] * mm.m[15] +
		mm.m[4] * mm.m[9] * mm.m[2] * mm.m[15] +
		mm.m[8] * mm.m[1] * mm.m[6] * mm.m[15] -
		mm.m[0] * mm.m[9] * mm.m[6] * mm.m[15] -
		mm.m[4] * mm.m[1] * mm.m[10] * mm.m[15] +
		mm.m[0] * mm.m[5] * mm.m[10] * mm.m[15];
}

/* returns a 16-element array that is the inverse of a 16-element array (4x4
matrix). see http://www.euclideanspace.com/maths/algebra/matrix/functions/inverse/fourD/index.htm */
MFUN::mat4 inverse (const MFUN::mat4& mm) {
	float det = MFUN::determinant (mm);
	/* there is no inverse if determinant is zero (not likely unless scale is
	broken) */
	if (0.0f == det) {
		fprintf (stderr, "WARNING. matrix has no determinant. can not invert\n");
		return mm;
	}
	float inv_det = 1.0f / det;
	
	return MFUN::mat4 (
		inv_det * (
			mm.m[9] * mm.m[14] * mm.m[7] - mm.m[13] * mm.m[10] * mm.m[7] +
			mm.m[13] * mm.m[6] * mm.m[11] - mm.m[5] * mm.m[14] * mm.m[11] -
			mm.m[9] * mm.m[6] * mm.m[15] + mm.m[5] * mm.m[10] * mm.m[15]
		),
		inv_det * (
			mm.m[13] * mm.m[10] * mm.m[3] - mm.m[9] * mm.m[14] * mm.m[3] -
			mm.m[13] * mm.m[2] * mm.m[11] + mm.m[1] * mm.m[14] * mm.m[11] +
			mm.m[9] * mm.m[2] * mm.m[15] - mm.m[1] * mm.m[10] * mm.m[15]
		),
		inv_det * (
			mm.m[5] * mm.m[14] * mm.m[3] - mm.m[13] * mm.m[6] * mm.m[3] +
			mm.m[13] * mm.m[2] * mm.m[7] - mm.m[1] * mm.m[14] * mm.m[7] -
			mm.m[5] * mm.m[2] * mm.m[15] + mm.m[1] * mm.m[6] * mm.m[15]
		),
		inv_det * (
			mm.m[9] * mm.m[6] * mm.m[3] - mm.m[5] * mm.m[10] * mm.m[3] -
			mm.m[9] * mm.m[2] * mm.m[7] + mm.m[1] * mm.m[10] * mm.m[7] +
			mm.m[5] * mm.m[2] * mm.m[11] - mm.m[1] * mm.m[6] * mm.m[11]
		),
		inv_det * (
			mm.m[12] * mm.m[10] * mm.m[7] - mm.m[8] * mm.m[14] * mm.m[7] -
			mm.m[12] * mm.m[6] * mm.m[11] + mm.m[4] * mm.m[14] * mm.m[11] +
			mm.m[8] * mm.m[6] * mm.m[15] - mm.m[4] * mm.m[10] * mm.m[15]
		),
		inv_det * (
			mm.m[8] * mm.m[14] * mm.m[3] - mm.m[12] * mm.m[10] * mm.m[3] +
			mm.m[12] * mm.m[2] * mm.m[11] - mm.m[0] * mm.m[14] * mm.m[11] -
			mm.m[8] * mm.m[2] * mm.m[15] + mm.m[0] * mm.m[10] * mm.m[15]
		),
		inv_det * (
			mm.m[12] * mm.m[6] * mm.m[3] - mm.m[4] * mm.m[14] * mm.m[3] -
			mm.m[12] * mm.m[2] * mm.m[7] + mm.m[0] * mm.m[14] * mm.m[7] +
			mm.m[4] * mm.m[2] * mm.m[15] - mm.m[0] * mm.m[6] * mm.m[15]
		),
		inv_det * (
			mm.m[4] * mm.m[10] * mm.m[3] - mm.m[8] * mm.m[6] * mm.m[3] +
			mm.m[8] * mm.m[2] * mm.m[7] - mm.m[0] * mm.m[10] * mm.m[7] -
			mm.m[4] * mm.m[2] * mm.m[11] + mm.m[0] * mm.m[6] * mm.m[11]
		),
		inv_det * (
			mm.m[8] * mm.m[13] * mm.m[7] - mm.m[12] * mm.m[9] * mm.m[7] +
			mm.m[12] * mm.m[5] * mm.m[11] - mm.m[4] * mm.m[13] * mm.m[11] -
			mm.m[8] * mm.m[5] * mm.m[15] + mm.m[4] * mm.m[9] * mm.m[15]
		),
		inv_det * (
			mm.m[12] * mm.m[9] * mm.m[3] - mm.m[8] * mm.m[13] * mm.m[3] -
			mm.m[12] * mm.m[1] * mm.m[11] + mm.m[0] * mm.m[13] * mm.m[11] +
			mm.m[8] * mm.m[1] * mm.m[15] - mm.m[0] * mm.m[9] * mm.m[15]
		),
		inv_det * (
			mm.m[4] * mm.m[13] * mm.m[3] - mm.m[12] * mm.m[5] * mm.m[3] +
			mm.m[12] * mm.m[1] * mm.m[7] - mm.m[0] * mm.m[13] * mm.m[7] -
			mm.m[4] * mm.m[1] * mm.m[15] + mm.m[0] * mm.m[5] * mm.m[15]
		),
		inv_det * (
			mm.m[8] * mm.m[5] * mm.m[3] - mm.m[4] * mm.m[9] * mm.m[3] -
			mm.m[8] * mm.m[1] * mm.m[7] + mm.m[0] * mm.m[9] * mm.m[7] +
			mm.m[4] * mm.m[1] * mm.m[11] - mm.m[0] * mm.m[5] * mm.m[11]
		),
		inv_det * (
			mm.m[12] * mm.m[9] * mm.m[6] - mm.m[8] * mm.m[13] * mm.m[6] -
			mm.m[12] * mm.m[5] * mm.m[10] + mm.m[4] * mm.m[13] * mm.m[10] +
			mm.m[8] * mm.m[5] * mm.m[14] - mm.m[4] * mm.m[9] * mm.m[14]
		),
		inv_det * (
			mm.m[8] * mm.m[13] * mm.m[2] - mm.m[12] * mm.m[9] * mm.m[2] +
			mm.m[12] * mm.m[1] * mm.m[10] - mm.m[0] * mm.m[13] * mm.m[10] -
			mm.m[8] * mm.m[1] * mm.m[14] + mm.m[0] * mm.m[9] * mm.m[14]
		),
		inv_det * (
			mm.m[12] * mm.m[5] * mm.m[2] - mm.m[4] * mm.m[13] * mm.m[2] -
			mm.m[12] * mm.m[1] * mm.m[6] + mm.m[0] * mm.m[13] * mm.m[6] +
			mm.m[4] * mm.m[1] * mm.m[14] - mm.m[0] * mm.m[5] * mm.m[14]
		),
		inv_det * (
			mm.m[4] * mm.m[9] * mm.m[2] - mm.m[8] * mm.m[5] * mm.m[2] +
			mm.m[8] * mm.m[1] * mm.m[6] - mm.m[0] * mm.m[9] * mm.m[6] -
			mm.m[4] * mm.m[1] * mm.m[10] + mm.m[0] * mm.m[5] * mm.m[10]
		)
	);
}

//new
MFUN::mat3 MFUN::transpose (const MFUN::mat3& mm) {
	return MFUN::mat3 (
		mm.m[0], mm.m[3], mm.m[6],
		mm.m[1], mm.m[4], mm.m[7],
		mm.m[2], mm.m[5], mm.m[8]
	);
}

// returns a 16-element array flipped on the main diagonal
MFUN::mat4 MFUN::transpose (const MFUN::mat4& mm) {
	return MFUN::mat4 (
		mm.m[0], mm.m[4], mm.m[8], mm.m[12],
		mm.m[1], mm.m[5], mm.m[9], mm.m[13],
		mm.m[2], mm.m[6], mm.m[10], mm.m[14],
		mm.m[3], mm.m[7], mm.m[11], mm.m[15]
	);
}

/*--------------------------AFFINE MATRIX FUNCTIONS---------------------------*/
// translate a 4d matrix with xyz array
MFUN::mat4 MFUN::translate (const MFUN::mat4& m, const MFUN::vec3& v) {
	MFUN::mat4 m_t = MFUN::identity_mat4 ();
	m_t.m[12] = v.v[0];
	m_t.m[13] = v.v[1];
	m_t.m[14] = v.v[2];
	return m_t * m;
}

// rotate around x axis by an angle in degrees
MFUN::mat4 MFUN::rotate_x_deg (const MFUN::mat4& m, float deg) {
	// convert to radians
	float rad = deg * ONE_DEG_IN_RAD;
	MFUN::mat4 m_r = MFUN::identity_mat4 ();
	m_r.m[5] = cos (rad);
	m_r.m[9] = -sin (rad);
	m_r.m[6] = sin (rad);
	m_r.m[10] = cos (rad);
	return m_r * m;
}

// rotate around y axis by an angle in degrees
MFUN::mat4 MFUN::rotate_y_deg (const MFUN::mat4& m, float deg) {
	// convert to radians
	float rad = deg * ONE_DEG_IN_RAD;
	MFUN::mat4 m_r = MFUN::identity_mat4 ();
	m_r.m[0] = cos (rad);
	m_r.m[8] = sin (rad);
	m_r.m[2] = -sin (rad);
	m_r.m[10] = cos (rad);
	return m_r * m;
}

// rotate around z axis by an angle in degrees
MFUN::mat4 MFUN::rotate_z_deg (const MFUN::mat4& m, float deg) {
	// convert to radians
	float rad = deg * ONE_DEG_IN_RAD;
	MFUN::mat4 m_r = MFUN::identity_mat4 ();
	m_r.m[0] = cos (rad);
	m_r.m[4] = -sin (rad);
	m_r.m[1] = sin (rad);
	m_r.m[5] = cos (rad);
	return m_r * m;
}

// scale a matrix by [x, y, z]
MFUN::mat4 MFUN::scale (const MFUN::mat4& m, const MFUN::vec3& v) {
	MFUN::mat4 a = MFUN::identity_mat4 ();
	a.m[0] = v.v[0];
	a.m[5] = v.v[1];
	a.m[10] = v.v[2];
	return a * m;
}

/*-----------------------VIRTUAL CAMERA MATRIX FUNCTIONS----------------------*/
// returns a view matrix using the opengl lookAt style. COLUMN ORDER.
MFUN::mat4 MFUN::look_at (const MFUN::vec3& cam_pos, MFUN::vec3 targ_pos, const MFUN::vec3& up) {
	// inverse translation
	MFUN::mat4 p = MFUN::identity_mat4 ();
	p = MFUN::translate (p, MFUN::vec3 (-cam_pos.v[0], -cam_pos.v[1], -cam_pos.v[2]));
	// distance vector
	MFUN::vec3 d = targ_pos - cam_pos;
	// forward vector
	MFUN::vec3 f = MFUN::normalise (d);
	// right vector
	MFUN::vec3 r = MFUN::normalise (MFUN::cross (f, up));
	// real up vector
	MFUN::vec3 u = MFUN::normalise (MFUN::cross (r, f));
	MFUN::mat4 ori = MFUN::identity_mat4 ();
	ori.m[0] = r.v[0];
	ori.m[4] = r.v[1];
	ori.m[8] = r.v[2];
	ori.m[1] = u.v[0];
	ori.m[5] = u.v[1];
	ori.m[9] = u.v[2];
	ori.m[2] = -f.v[0];
	ori.m[6] = -f.v[1];
	ori.m[10] = -f.v[2];
	
	return ori * p;//p * ori;
}

// returns a perspective function mimicking the opengl projection style.
MFUN::mat4 MFUN::perspective (float fovy, float aspect, float near, float far) {
	float fov_rad = fovy * ONE_DEG_IN_RAD;
	float range = tan (fov_rad / 2.0f) * near;
	float sx = (2.0f * near) / (range * aspect + range * aspect);
	float sy = near / range;
	float sz = -(far + near) / (far - near);
	float pz = -(2.0f * far * near) / (far - near);
	MFUN::mat4 m = MFUN::zero_mat4 (); // make sure bottom-right corner is zero
	m.m[0] = sx;
	m.m[5] = sy;
	m.m[10] = sz;
	m.m[14] = pz;
	m.m[11] = -1.0f;
	return m;
}

/*----------------------------HAMILTON IN DA HOUSE!---------------------------*/
MFUN::versor::versor () { }

MFUN::versor MFUN::versor::operator/ (float rhs) {
	MFUN::versor result;
	result.q[0] = q[0] / rhs;
	result.q[1] = q[1] / rhs;
	result.q[2] = q[2] / rhs;
	result.q[3] = q[3] / rhs;
	return result;
}

MFUN::versor MFUN::versor::operator* (float rhs) {
	MFUN::versor result;
	result.q[0] = q[0] * rhs;
	result.q[1] = q[1] * rhs;
	result.q[2] = q[2] * rhs;
	result.q[3] = q[3] * rhs;
	return result;
}

void MFUN::print (const MFUN::versor& q) {
	printf ("[%.2f ,%.2f, %.2f, %.2f]\n", q.q[0], q.q[1], q.q[2], q.q[3]);
}

MFUN::versor MFUN::versor::operator* (const MFUN::versor& rhs) {
	MFUN::versor result;
	result.q[0] = rhs.q[0] * q[0] - rhs.q[1] * q[1] -
		rhs.q[2] * q[2] - rhs.q[3] * q[3];
	result.q[1] = rhs.q[0] * q[1] + rhs.q[1] * q[0] -
		rhs.q[2] * q[3] + rhs.q[3] * q[2];
	result.q[2] = rhs.q[0] * q[2] + rhs.q[1] * q[3] +
		rhs.q[2] * q[0] - rhs.q[3] * q[1];
	result.q[3] = rhs.q[0] * q[3] - rhs.q[1] * q[2] +
		rhs.q[2] * q[1] + rhs.q[3] * q[0];
	// re-normalise in case of mangling
	return MFUN::normalise (result);
}

MFUN::versor MFUN::versor::operator+ (const MFUN::versor& rhs) {
	MFUN::versor result;
	result.q[0] = rhs.q[0] + q[0];
	result.q[1] = rhs.q[1] + q[1];
	result.q[2] = rhs.q[2] + q[2];
	result.q[3] = rhs.q[3] + q[3];
	// re-normalise in case of mangling
	return MFUN::normalise (result);
}

MFUN::versor MFUN::quat_from_axis_rad (float radians, float x, float y, float z) {
	MFUN::versor result;
	result.q[0] = cos (radians / 2.0);
	result.q[1] = sin (radians / 2.0) * x;
	result.q[2] = sin (radians / 2.0) * y;
	result.q[3] = sin (radians / 2.0) * z;
	return result;
}

MFUN::versor MFUN::quat_from_axis_deg (float degrees, float x, float y, float z) {
	return MFUN::quat_from_axis_rad (ONE_DEG_IN_RAD * degrees, x, y, z);
}

MFUN::mat4 MFUN::quat_to_mat4 (const MFUN::versor& q) {
	float w = q.q[0];
	float x = q.q[1];
	float y = q.q[2];
	float z = q.q[3];
	return MFUN::mat4 (
		1.0f - 2.0f * y * y - 2.0f * z * z,
		2.0f * x * y + 2.0f * w * z,
		2.0f * x * z - 2.0f * w * y,
		0.0f,
		2.0f * x * y - 2.0f * w * z,
		1.0f - 2.0f * x * x - 2.0f * z * z,
		2.0f * y * z + 2.0f * w * x,
		0.0f,
		2.0f * x * z + 2.0f * w * y,
		2.0f * y * z - 2.0f * w * x,
		1.0f - 2.0f * x * x - 2.0f * y * y,
		0.0f,
		0.0f,
		0.0f,
		0.0f,
		1.0f
	);
}

MFUN::versor MFUN::normalise (MFUN::versor& q) {
	// norm(q) = q / magnitude (q)
	// magnitude (q) = sqrt (w*w + x*x...)
	// only compute sqrt if interior sum != 1.0
	float sum =
		q.q[0] * q.q[0] + q.q[1] * q.q[1] +
		q.q[2] * q.q[2] + q.q[3] * q.q[3];
	// NB: floats have min 6 digits of precision
	const float thresh = 0.0001f;
	if (fabs (1.0f - sum) < thresh) {
		return q;
	}
	float mag = sqrt (sum);
	return q / mag;
}

float MFUN::dot (const MFUN::versor& q, const MFUN::versor& r) {
	return q.q[0] * r.q[0] + q.q[1] * r.q[1] + q.q[2] * r.q[2] + q.q[3] * r.q[3];
}

MFUN::versor slerp (MFUN::versor& q, MFUN::versor& r, float t) {
	// angle between q0-q1
	float cos_half_theta = MFUN::dot (q, r);
	// as found here http://stackoverflow.com/questions/2886606/flipping-issue-when-interpolating-rotations-using-quaternions
	// if dot product is negative then one quaternion should be negated, to make
	// it take the short way around, rather than the long way
	// yeah! and furthermore Susan, I had to recalculate the d.p. after this
	if (cos_half_theta < 0.0f) {
		for (int i = 0; i < 4; i++) {
			q.q[i] *= -1.0f;
		}
		cos_half_theta = MFUN::dot (q, r);
	}
	// if qa=qb or qa=-qb then theta = 0 and we can return qa
	if (fabs (cos_half_theta) >= 1.0f) {
		return q;
	}
	// Calculate temporary values
	float sin_half_theta = sqrt (1.0f - cos_half_theta * cos_half_theta);
	// if theta = 180 degrees then result is not fully defined
	// we could rotate around any axis normal to qa or qb
	MFUN::versor result;
	if (fabs (sin_half_theta) < 0.001f) {
		for (int i = 0; i < 4; i++) {
			result.q[i] = (1.0f - t) * q.q[i] + t * r.q[i];
		}
		return result;
	}
	float half_theta = acos (cos_half_theta);
	float a = sin ((1.0f - t) * half_theta) / sin_half_theta;
	float b = sin (t * half_theta) / sin_half_theta;
	for (int i = 0; i < 4; i++) {
		result.q[i] = q.q[i] * a + r.q[i] * b;
	}
	return result;
}
