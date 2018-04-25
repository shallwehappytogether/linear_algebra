#pragma once

#include <iostream>
#include <glm\glm.hpp>
#include <linear_algebra\matrix.h>
#include <glm\gtc\type_ptr.hpp>
#include <glm\gtx\quaternion.hpp>
#include <linear_algebra\quaternion.h>

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
			os << glm::value_ptr(mat)[j * 4 + i] << " ";
		if (i != 3)
			os << "\n";
	}
	return os;
}

std::ostream& operator<<(std::ostream &os, const lin::quaternion &q)
{
	for (int j = 0; j < 4; ++j)
		os << q[j] << " ";
	return os;
}

std::ostream& operator<<(std::ostream &os, const glm::quat &q)
{
	for (int j = 0; j < 4; ++j)
		os << q[j] << " ";
	return os;
}

std::ostream& operator<<(std::ostream &os, const lin::vector_3d &q)
{
	for (int j = 0; j < 3; ++j)
		os << q[j] << " ";
	return os;
}

std::ostream& operator<<(std::ostream &os, const glm::vec3 &q)
{
	for (int j = 0; j < 3; ++j)
		os << q[j] << " ";
	return os;
}