
#include <linear_algebra/matrix.h>
#include <linear_algebra/vector.h>
#include <linear_algebra/graphic/transform.h>
#include <iostream>
using namespace lin;

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

	// vector normalize
	vector_4d v4;
	auto v5 = v4.normalized(); // v4 itself is not affected
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
	square_matrix_4d m0; // zero matrix

	// get an identity matrix
	square_matrix_4d m1 = make_identity_matrix<double, 4>();

	// matrix transpose
	square_matrix_4d m2;
	square_matrix_4d m3 = m2.get_transpose(); // m2 itself is not affected
	m2.transpose(); // transpose m2 itself

	// matrix inverse
	square_matrix_4d m4 = make_identity_matrix<double, 4>();
	square_matrix_4d m5 = m4.inverse(); // m4 itself is not affected
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
	square_matrix_4d translMatrix, rotMatrix, scaleMatrix;
	vector_3d localPosition;
	vector_3d worldPosition = (translMatrix * rotMatrix * scaleMatrix * localPosition.homogeneous()).reduce();

	// affine transform, treat vector as row vector
	square_matrix_4d translMatrixT, rotMatrixT, scaleMatrixT;
	vector_3d localPositionT;
	vector_3d worldPositionT = (localPositionT.homogeneous() * scaleMatrixT * rotMatrixT * translMatrixT).reduce();

	auto m = lin::translating<>::get(vector_3d(1.0, 2.0, 3.0));

	lin::transform<> parentTrans;
	parentTrans.
		then<lin::rotating>(radius(45.0), vector_3d(1.0, 1.0, 1.0)).
		then<lin::scalling>(1.0, 1.0, 1.0).
		then<lin::translating>(1.0, 2.0, 3.0);
	lin::transform<> myTrans;
	myTrans.
		then<lin::rotating>(radius(45.0), vector_3d(1.0, 1.0, 1.0)).
		then<lin::scalling>(1.0, 1.0, 1.0).
		then<lin::translating>(1.0, 2.0, 3.0).
		then(parentTrans);
	vector_3d vertex;
	auto vertexTransformed = myTrans.apply(vertex.homogeneous()).reduce();


	std::cout << "Bye!\n";
}