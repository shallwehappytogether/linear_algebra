
#include <linear_algebra/matrix.h>
#include <linear_algebra/vector.h>
#include <linear_algebra/graphic/transform.h>
#include <linear_algebra/trigonometry.h>
#include <iostream>
#include <glm\glm.hpp>
#include <glm\gtx\transform.hpp>
#include <glm\gtx\euler_angles.hpp>
#include "print.h"
using namespace lin;

template <typename LINTy, typename GLMTy>
void compare(const char *content, const LINTy &linval, const GLMTy &glmval)
{
	std::cout << "===== " << content << " =====\n"
		"lin\n" << linval << "\n" <<
		"glm\n" << glmval << "\n\n";
}

int main()
{
	// vector construction
	vector_4d v0; // zero vector
	vector_4d v1{ 1, 2 }; // with components (1.0, 2.0, 0.0, 0.0)
	std::cout << v1.y << " " << v1[1] << "\n";
	vector_4d v2{ 1, 2, 3, 4 }; // with components (1.0, 2.0, 3.0, 4.0)
	vector_4d v3(1, 2); // same as v1
	if (v1 == v3)
		std::cout << "v1 == v3\n";
	if (+v1 != -v3)
		std::cout << "+v1 != -v3\n";
	std::cout << sizeof(v0) << "\n";
	auto negv2 = -v2;

	vector_4f v00 = v0;

	auto dv = v1 - v0;

	// vector normalize
	vector_4d v4;
	auto v5 = v4.get_normalized(); // v4 itself is not affected
	v4.normalize(); // normalize v4 itself, so that v4 == v5 then

	// vector dot product
	auto dotv = std::inner_product(v1, v2);

	// vector cross product, only defined for 3d vectors
	vector_3d v6 = cross_product(v1.reduce(), v2.reduce());

	// make axis basis
	auto xaxis = vector_3d::basis<x_coord>();
	auto yaxis = vector_3d::basis<y_coord>();
	auto zaxis = cross_product(xaxis, yaxis);

	// matrix construction
	matrix_4d m0; // zero matrix

	// get an identity matrix
	matrix_4d m1 = identity_matrix<double, 4>();

	// matrix transpose
	matrix_4d m2;
	matrix_4d m3 = m2.get_transpose(); // m2 itself is not affected
	m2.transpose(); // transpose m2 itself

	// matrix inverse
	matrix_4d m4 = identity_matrix<double, 4>();
	matrix_4d m5 = m4.inverse(); // m4 itself is not affected
	m4.invert(); // inverse m4 itself

	// no inverse or inverse_to_assign or transpose_to_assign is
	// defined for non-sqare matrix
	basic_matrix<double, 3, 4> nonsqrm;
	// the following code will static_assert with false
	//nonsqrm.inverse();
	//nonsqrm.inverse_to_assign();
	//nonsqrm.transpose_to_assign();
	// but tranpose is allowed, it return a new matrix
	auto nonsqrmT = nonsqrm.transpose();

	// affine transform, treat vector as column vector
	matrix_4d translMatrix, rotMatrix, scaleMatrix;
	vector_3d localPosition;
	vector_3d worldPosition = (translMatrix * rotMatrix * scaleMatrix * localPosition.homogeneous()).reduce();

	// affine transform, treat vector as row vector
	matrix_4d translMatrixT, rotMatrixT, scaleMatrixT;
	vector_3d localPositionT;
	vector_3d worldPositionT = (localPositionT.homogeneous() * scaleMatrixT * rotMatrixT * translMatrixT).reduce();

	using GLM_VEC_SYS = lin::c_v;

	lin::transform<double, GLM_VEC_SYS> parentTrans;
	parentTrans.
		then<lin::rotating>(radius(45.0), vector_3d(1.0, 1.0, 1.0)).
		then<lin::scalling>(1.0, 1.0, 1.0).
		then<lin::translating>(1.0, 2.0, 3.0);
	lin::transform<double, GLM_VEC_SYS> myTrans;
	myTrans.
		then<lin::rotating>(radius(45.0), vector_3d(1.0, 1.0, 1.0)).
		then<lin::scalling>(1.0, 1.0, 1.0).
		then<lin::translating>(1.0, 2.0, 3.0).
		then(parentTrans);
	vector_3d vertex;
	auto vertexTransformed = myTrans.apply(vertex.homogeneous()).reduce();

	glm::mat4 parentTransGLM;
	parentTransGLM = glm::translate(parentTransGLM, glm::vec3(1.0, 2.0, 3.0));
	parentTransGLM = glm::scale(parentTransGLM, glm::vec3(1.0, 1.0, 1.0));
	parentTransGLM = glm::rotate(parentTransGLM, glm::radians(45.0f), glm::normalize(glm::vec3(1.0, 1.0, 1.0)));

	auto randomAxisLin = vector_3d(1., 2., 3.).get_normalized();
	auto randomAngleLin = lin::radius(76.);
	auto randomQuatLin = lin::rotation_quaternion(randomAngleLin, randomAxisLin);
	auto randomAxis2Lin = vector_3d(45., -84., 23.).get_normalized();
	auto randomAngle2Lin = lin::radius(-26.);
	auto randomQuat2Lin = lin::rotation_quaternion(randomAngle2Lin, randomAxis2Lin);
	auto randomScaleLin = lin::vector_3d(2.0, 2.0, 2.0);
	auto randomMoveLin = lin::vector_3d(4.0, 5.0, 6.0);

	auto randomAxisGlm = glm::normalize(glm::vec3(1.f, 2.f, 3.f));
	auto randomAngleGlm = glm::radians(76.f);
	auto randomQuatGlm = glm::angleAxis(randomAngleGlm, randomAxisGlm);
	auto randomAxis2Glm = glm::normalize(glm::vec3(45.f, -84.f, 23.f));
	auto randomAngle2Glm = glm::radians(-26.f);
	auto randomQuat2Glm = glm::angleAxis(randomAngle2Glm, randomAxis2Glm);
	auto randomScaleGlm = glm::vec3(2.0f, 2.0f, 2.0f);
	auto randomMoveGlm = glm::vec3(4.0f, 5.0f, 6.0f);

	compare("Translation matrix",
		lin::translating<GLM_VEC_SYS>::get(randomMoveLin),
		glm::translate(glm::mat4(), randomMoveGlm));

	compare("Scalling matrix",
		lin::scalling<GLM_VEC_SYS>::get(randomScaleLin),
		glm::scale(glm::mat4(), randomScaleGlm));

	compare("Rotation matrix(Axis-Angle to matrix)",
		lin::rotating<GLM_VEC_SYS>::get(randomAngleLin, randomAxisLin),
		glm::rotate(glm::mat4(), randomAngleGlm, randomAxisGlm));

	compare("TRS matrix multiply",
		parentTrans.matrix(),
		parentTransGLM);

	compare("matrix inverse",
		parentTrans.matrix().inverse(),
		glm::inverse(parentTransGLM));

	compare("Apply quaternion to vector",
		randomQuatLin * vector_3d(10., 9., 8.),
		randomQuatGlm * glm::vec3(10.f, 9.f, 8.f));

	compare("Axis-Angle to Quaternion",
		randomQuatLin,
		randomQuatGlm);

	compare("Quaternion to matrix",
		lin::rotating<GLM_VEC_SYS>::get(randomQuatLin),
		glm::mat4_cast(randomQuatGlm));

	compare("Lookat matrix",
		lin::lookat_rhs<GLM_VEC_SYS>::get(vector_3d{ 400, 300, 200 }, vector_3d{ 1, 2, 3 }, vector_3d{ 0, 1, 0 }),
		glm::lookAtRH(glm::vec3{ 400, 300, 200 }, glm::vec3{ 400, 300, 200 } +glm::vec3{ 1, 2, 3 }, { 0, 1, 0 }));

	compare("Eular angle to quaternion",
		lin::to_quaternion(lin::eular_angle_xyz(lin::radius(45.0), lin::radius(45.0), lin::radius(45.0))),
		glm::toQuat(glm::eulerAngleXYZ(glm::radians(45.0), glm::radians(45.0), glm::radians(45.0))));

	compare("Quaternion multiply", 
		randomQuatLin * randomQuat2Lin,
		randomQuatGlm * randomQuat2Glm);

	compare("Quaternion multiply(reversed)",
		randomQuat2Lin * randomQuatLin,
		randomQuat2Glm * randomQuatGlm);
}