#pragma once

#include <linear_algebra/matrix.h>
#include <linear_algebra/vector.h>

namespace lin
{
	namespace column_vector
	{
		template <typename Ty>
		void translate(
			basic_square_matrix_4d<Ty> &dest_,
			const basic_vector_3d<Ty> &pos_)
		{
			dest_.element_at(3, 0) = pos_.x;
			dest_.element_at(3, 1) = pos_.y;
			dest_.element_at(3, 2) = pos_.z;
		}

		template <typename Ty>
		basic_square_matrix_4d<Ty> translation(
			const basic_vector_3d<Ty> &pos_)
		{
			auto result = make_identity_matrix<Ty, 4>();
			translate(result, pos_);
			return result;
		}
	}

	namespace row_vector
	{
		template <typename Ty>
		void translate(
			basic_square_matrix_4d<Ty> &dest_,
			const basic_vector_3d<Ty> &pos_)
		{
			dest_.element_at(0, 3) = pos_.x;
			dest_.element_at(1, 3) = pos_.y;
			dest_.element_at(2, 3) = pos_.z;
		}

		template <typename Ty>
		basic_square_matrix_4d<Ty> translation(
			const basic_vector_3d<Ty> &pos_)
		{
			auto result = make_identity_matrix<Ty, 4>();
			translate(result, pos_);
			return result;
		}
	}
}