#pragma once

#include <cstddef>

namespace lin
{
	inline namespace convention
	{
		/* CONSTEXPR VARIANBLES x_coord y_coord z_coord w_coord
		Convention signs for x, y, z, w denotes of vector components.
		*/
		constexpr std::size_t x_coord = 0, y_coord = 1, z_coord = 2, w_coord = 3;
	}
}