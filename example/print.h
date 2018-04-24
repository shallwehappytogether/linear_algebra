#pragma once

#include <iostream>
#include <glm\glm.hpp>
#include <linear_algebra\matrix.h>

std::ostream& operator<<(std::ostream &os, const lin::square_matrix_4d &mat)
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
			os << mat.element_at(i, j) << " ";
		if (i != 3)
			os << "\n";
	}
	return os;
}

std::ostream& operator<<(std::ostream &os, const glm::mat4 &mat)
{
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
			os << mat[j][i] << " ";
		if (i != 3)
			os << "\n";
	}
	return os;
}