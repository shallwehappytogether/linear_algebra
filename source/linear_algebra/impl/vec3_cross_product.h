#pragma once

#include <type_traits>

namespace lin
{
	namespace impl
	{
		template <typename ResultTy, typename LeftTy, typename RightTy>
		void vec3_cross_product(
			ResultTy *result,
			const LeftTy *left,
			const RightTy *right)
		{
			result[0] = left[1] * right[2] - left[2] * right[1];
			result[1] = left[2] * right[0] - left[0] * right[2];
			result[2] = left[0] * right[1] - left[1] * right[0];
		}
	}
}