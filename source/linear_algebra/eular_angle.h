#pragma once

#include <linear_algebra/convention.h>
#include <linear_algebra/impl/_number_array.h>

namespace lin
{
	template <typename Ty, std::size_t ...Basis>
	class basic_eular_angle
		:public impl::_number_array<Ty, sizeof...(Basis)>
	{
		using MyBase = impl::_number_array<Ty, sizeof...(Basis)>;
	public:
		using MyBase::_number_array;
	};

	template <typename Ty>
	using basic_eular_angle_xyz = basic_eular_angle<Ty, x_coord, y_coord, z_coord>;

	using eular_angle_xyz = basic_eular_angle_xyz<double>;

	using eular_angle_xyz_f = basic_eular_angle_xyz<float>;
}