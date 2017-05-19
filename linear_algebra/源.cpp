
#include "linear_algebra.h"
#include <iostream>
using namespace linear_algebra;

int main()
{
	basic_vector_4d<double> v1, v2;
	basic_matrix<double, 4, 4> mat;

	auto dotv = std::inner_product(v1, v2);
	std::cout << dotv << "\n";

	v1 *= mat;
	mat *= v1;

	basic_vector_4d<double> v3{1, 2};
	basic_vector_4d<double> v4{ 1, 2, 3, 4 };

	basic_vector_4d<double> v3_1(1, 2);

	std::cout << "pause\n";
}