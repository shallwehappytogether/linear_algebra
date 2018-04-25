#pragma once

#include <linear_algebra/matrix.h>
#include <linear_algebra/vector.h>
#include <linear_algebra/quaternion.h>

namespace lin
{
	struct r_v {};

	struct c_v {};

	using default_vector_system = c_v;

	namespace impl
	{
		template <typename VectorSystem, typename Ty>
		static inline void set_c_v(
			basic_square_matrix_4d<Ty> &m,
			typename basic_square_matrix_4d<Ty>::dimension_type r,
			typename basic_square_matrix_4d<Ty>::dimension_type c,
			const Ty &val_)
		{
			if constexpr (std::is_same_v<VectorSystem, c_v>)
				m.element_at(r, c) = val_;
			else
				m.element_at(c, r) = val_;
		}
	}

	template <typename VectorSystem = default_vector_system>
	struct translating
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const basic_vector_3d<Ty> &pos_)
		{
			return get(pos_.x, pos_.y, pos_.z);
		}

		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &x_,
			const Ty &y_,
			const Ty &z_)
		{
			auto i = identity_matrix<Ty, 4>();
			fill(i, x_, y_, z_);
			return i;
		}

		template <typename Ty>
		static void fill(
			basic_square_matrix_4d<Ty> &dest_,
			const Ty &x_,
			const Ty &y_,
			const Ty &z_)
		{
			impl::set_c_v<VectorSystem>(dest_, 0, 3, x_);
			impl::set_c_v<VectorSystem>(dest_, 1, 3, y_);
			impl::set_c_v<VectorSystem>(dest_, 2, 3, z_);
		}
	};

	template <typename = default_vector_system>
	struct scalling
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const basic_vector_3d<Ty> &scale_)
		{
			return get(scale_.x, scale_.y, scale_.z);
		}

		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &x_,
			const Ty &y_,
			const Ty &z_)
		{
			auto i = identity_matrix<Ty, 4>();
			i.element_at(0, 0) = x_;
			i.element_at(1, 1) = y_;
			i.element_at(2, 2) = z_;
			return i;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct rotating_x
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &rangle_)
		{
			auto i = identity_matrix<Ty, 4>();
			Ty cosa = std::cos(rangle_), sina = std::sin(rangle_);
			impl::set_c_v<VectorSystem>(i, 1, 1, cosa);
			impl::set_c_v<VectorSystem>(i, 2, 2, cosa);
			impl::set_c_v<VectorSystem>(i, 1, 2, -sina);
			impl::set_c_v<VectorSystem>(i, 2, 1, sina);
			return i;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct rotating_y
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &rangle_)
		{
			auto i = identity_matrix<Ty, 4>();
			Ty cosa = std::cos(rangle_), sina = std::sin(rangle_);
			impl::set_c_v<VectorSystem>(i, 0, 0, cosa);
			impl::set_c_v<VectorSystem>(i, 2, 2, cosa);
			impl::set_c_v<VectorSystem>(i, 0, 2, sina);
			impl::set_c_v<VectorSystem>(i, 2, 0, -sina);
			return i;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct rotating_z
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &rangle_)
		{
			auto i = identity_matrix<Ty, 4>();
			Ty cosa = std::cos(rangle_), sina = std::sin(rangle_);
			impl::set_c_v<VectorSystem>(i, 0, 0, cosa);
			impl::set_c_v<VectorSystem>(i, 1, 1, cosa);
			impl::set_c_v<VectorSystem>(i, 0, 1, -sina);
			impl::set_c_v<VectorSystem>(i, 1, 0, sina);
			return i;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct rotating
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &rangle_, const basic_vector_3d<Ty> &axis_)
		{
			auto u = axis_.get_normalized();
			auto i = identity_matrix<Ty, 4>();
			Ty cosa = std::cos(rangle_), sina = std::sin(rangle_);
			Ty xcosa = 1 - cosa;
			impl::set_c_v<VectorSystem>(i, 0, 0, cosa + u.x * u.x * xcosa);
			impl::set_c_v<VectorSystem>(i, 0, 1, u.x * u.y * xcosa - u.z * sina);
			impl::set_c_v<VectorSystem>(i, 0, 2, u.x * u.z * xcosa + u.y * sina);
			impl::set_c_v<VectorSystem>(i, 1, 0, u.y * u.x * xcosa + u.z * sina);
			impl::set_c_v<VectorSystem>(i, 1, 1, cosa + u.y * u.y * xcosa);
			impl::set_c_v<VectorSystem>(i, 1, 2, u.y * u.z * xcosa - u.x * sina);
			impl::set_c_v<VectorSystem>(i, 2, 0, u.z * u.x * xcosa - u.y * sina);
			impl::set_c_v<VectorSystem>(i, 2, 1, u.z * u.y * xcosa + u.x * sina);
			impl::set_c_v<VectorSystem>(i, 2, 2, cosa + u.z * u.z * xcosa);
			return i;
		}

		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const basic_quaternion<Ty> &quat_)
		{
			auto q = quat_.get_normalized();
			auto qx = q.x;
			auto qy = q.y;
			auto qz = q.z;
			auto qw = q.w;
			auto i = identity_matrix<Ty, 4>();
			impl::set_c_v<VectorSystem>(i, 0, 0, 1 - 2 * (qy * qy + qz * qz));
			impl::set_c_v<VectorSystem>(i, 0, 1, 2 * (qx * qy - qz * qw));
			impl::set_c_v<VectorSystem>(i, 0, 2, 2 * (qx * qz + qy * qw));
			impl::set_c_v<VectorSystem>(i, 1, 0, 2 * (qx * qy + qz * qw));
			impl::set_c_v<VectorSystem>(i, 1, 1, 1 - 2 * (qx * qx + qz * qz));
			impl::set_c_v<VectorSystem>(i, 1, 2, 2 * (qy * qz - qx * qw));
			impl::set_c_v<VectorSystem>(i, 2, 0, 2 * (qx * qz - qy * qw));
			impl::set_c_v<VectorSystem>(i, 2, 1, 2 * (qy * qz + qx * qw));
			impl::set_c_v<VectorSystem>(i, 2, 2, 1 - 2 * (qx * qx + qy * qy));
			return i;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct lookat_rhs
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const basic_vector_3d<Ty> &pos_,
			const basic_vector_3d<Ty> &expectedz_,
			const basic_vector_3d<Ty> &expectedy_)
		{
			basic_square_matrix_4d<Ty> result;
			auto orthoZ = expectedz_.get_normalized();
			auto orthoX = cross_product(expectedz_, expectedy_).get_normalized();
			auto orthoY = cross_product(orthoX, orthoZ);
			impl::set_c_v<VectorSystem>(result, 0, 0, orthoX.x);
			impl::set_c_v<VectorSystem>(result, 0, 1, orthoX.y);
			impl::set_c_v<VectorSystem>(result, 0, 2, orthoX.z);
			impl::set_c_v<VectorSystem>(result, 0, 3, -std::inner_product(pos_, orthoX));
			impl::set_c_v<VectorSystem>(result, 1, 0, orthoY.x);
			impl::set_c_v<VectorSystem>(result, 1, 1, orthoY.y);
			impl::set_c_v<VectorSystem>(result, 1, 2, orthoY.z);
			impl::set_c_v<VectorSystem>(result, 1, 3, -std::inner_product(pos_, orthoY));
			// 因为是右手坐标系，z轴朝向应和视线方向相反。
			impl::set_c_v<VectorSystem>(result, 2, 0, -orthoZ.x);
			impl::set_c_v<VectorSystem>(result, 2, 1, -orthoZ.y);
			impl::set_c_v<VectorSystem>(result, 2, 2, -orthoZ.z);
			impl::set_c_v<VectorSystem>(result, 2, 3, std::inner_product(pos_, orthoZ));
			impl::set_c_v<VectorSystem>(result, 3, 0, static_cast<Ty>(0));
			impl::set_c_v<VectorSystem>(result, 3, 1, static_cast<Ty>(0));
			impl::set_c_v<VectorSystem>(result, 3, 2, static_cast<Ty>(0));
			impl::set_c_v<VectorSystem>(result, 3, 3, static_cast<Ty>(1));
			return result;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct lookat_lhs
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const basic_vector_3d<Ty> &pos_,
			const basic_vector_3d<Ty> &expectedz_,
			const basic_vector_3d<Ty> &expectedy_)
		{
			basic_square_matrix_4d<Ty> result;
			auto orthoZ = expectedz_.get_normalized();
			auto orthoX = cross_product(expectedy_, expectedz_).get_normalized();
			auto orthoY = cross_product(orthoZ, orthoX);
			impl::set_c_v<VectorSystem>(result, 0, 0, orthoX.x);
			impl::set_c_v<VectorSystem>(result, 0, 1, orthoX.y);
			impl::set_c_v<VectorSystem>(result, 0, 2, orthoX.z);
			impl::set_c_v<VectorSystem>(result, 0, 3, -std::inner_product(pos_, orthoX));
			impl::set_c_v<VectorSystem>(result, 1, 0, orthoY.x);
			impl::set_c_v<VectorSystem>(result, 1, 1, orthoY.y);
			impl::set_c_v<VectorSystem>(result, 1, 2, orthoY.z);
			impl::set_c_v<VectorSystem>(result, 1, 3, -std::inner_product(pos_, orthoY));
			impl::set_c_v<VectorSystem>(result, 2, 0, orthoZ.x);
			impl::set_c_v<VectorSystem>(result, 2, 1, orthoZ.y);
			impl::set_c_v<VectorSystem>(result, 2, 2, orthoZ.z);
			impl::set_c_v<VectorSystem>(result, 2, 3, -std::inner_product(pos_, orthoZ));
			impl::set_c_v<VectorSystem>(result, 3, 0, static_cast<Ty>(0));
			impl::set_c_v<VectorSystem>(result, 3, 1, static_cast<Ty>(0));
			impl::set_c_v<VectorSystem>(result, 3, 2, static_cast<Ty>(0));
			impl::set_c_v<VectorSystem>(result, 3, 3, static_cast<Ty>(1));
			return result;
		}
	};

	template <typename VectorSystem = default_vector_system, typename Ty = double>
	class transform
	{
	public:
		transform()
			:_mat(identity_matrix<Ty, 4>())
		{

		}

		void clear()
		{
			*this = transform();
		}

		transform& then(const transform &other_)
		{
			_multi(other_._mat);
			return *this;
		}

		template <template <typename> class Transform, typename ...Args>
		transform& then(Args&& ...args)
		{
			auto m = Transform<VectorSystem>::get(std::forward<Args>(args)...);
			_multi(m);
			return *this;
		}

		basic_vector_4d<Ty> apply(const basic_vector_4d<Ty> &v_) const
		{
			if constexpr (std::is_same_v<VectorSystem, c_v>)
				return _mat * v_;
			else
				return v_ * _mat;
		}

		const basic_square_matrix_4d<Ty>& matrix() const
		{
			return _mat;
		}
	private:
		basic_square_matrix_4d<Ty> _mat;

		void _multi(const basic_square_matrix_4d<Ty> &m)
		{
			if constexpr (std::is_same_v<VectorSystem, c_v>)
				_mat = m * _mat;
			else
				_mat = _mat * m;
		}
	};
}