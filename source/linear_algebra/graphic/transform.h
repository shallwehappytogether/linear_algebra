#pragma once

#include <linear_algebra/matrix.h>
#include <linear_algebra/vector.h>
#include <linear_algebra/quaternion.h>
#include <linear_algebra/eular_angle.h>
#include <linear_algebra/impl/utility.h>

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
			typename basic_square_matrix_4d<Ty>::row_dimension_type r,
			typename basic_square_matrix_4d<Ty>::column_dimension_type c,
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
			impl::set_c_v<VectorSystem>(result, 2, 3, -std::inner_product(pos_, -orthoZ));
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

	template <typename VectorSystem = default_vector_system>
	struct perspective_rhs
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &near_,
			const Ty &far_,
			const Ty &fovy_,
			const Ty &aspectratio_)
		{
			auto tanHalfFOVY = std::tan(fovy_ / static_cast<Ty>(2));
			basic_square_matrix_4d<Ty> result;
			impl::set_c_v<VectorSystem>(result, 0, 0, static_cast<Ty>(1) / (aspectratio_ * tanHalfFOVY));
			impl::set_c_v<VectorSystem>(result, 1, 1, static_cast<Ty>(1) / tanHalfFOVY);
			impl::set_c_v<VectorSystem>(result, 2, 3, (static_cast<Ty>(2) * near_ * far_) / (near_ - far_));

			impl::set_c_v<VectorSystem>(result, 2, 2, -(far_ + near_) / (far_ - near_));
			impl::set_c_v<VectorSystem>(result, 3, 2, -static_cast<Ty>(1));
			return result;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct perspective_lhs
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &near_,
			const Ty &far_,
			const Ty &fovy_,
			const Ty &aspectratio_)
		{
			auto tanHalfFOVY = std::tan(fovy_ / static_cast<Ty>(2));
			basic_square_matrix_4d<Ty> result;
			impl::set_c_v<VectorSystem>(result, 0, 0, static_cast<Ty>(1) / (aspectratio_ * tanHalfFOVY));
			impl::set_c_v<VectorSystem>(result, 1, 1, static_cast<Ty>(1) / tanHalfFOVY);
			impl::set_c_v<VectorSystem>(result, 2, 3, (static_cast<Ty>(2) * near_ * far_) / (near_ - far_));

			impl::set_c_v<VectorSystem>(result, 2, 2, (far_ + near_) / (far_ - near_));
			impl::set_c_v<VectorSystem>(result, 3, 2, static_cast<Ty>(1));
			return result;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct ortho_projection_rhs
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &left_,
			const Ty &right_,
			const Ty &top_,
			const Ty &bottom_,
			const Ty &near_,
			const Ty &far_)
		{
			basic_square_matrix_4d<Ty> result;
			impl::set_c_v<VectorSystem>(result, 0, 0, static_cast<Ty>(2) / (right_ - left_));
			impl::set_c_v<VectorSystem>(result, 0, 3, -(right_ + left_) / (right_ - left_));
			impl::set_c_v<VectorSystem>(result, 1, 1, static_cast<Ty>(2) / (top_ - bottom_));
			impl::set_c_v<VectorSystem>(result, 1, 3, -(top_ + bottom_) / (top_ - bottom_));
			impl::set_c_v<VectorSystem>(result, 2, 3, -(far_ + near_) / (far_ - near_));
			impl::set_c_v<VectorSystem>(result, 3, 3, static_cast<Ty>(1));

			impl::set_c_v<VectorSystem>(result, 2, 2, -static_cast<Ty>(2) / (far_ - near_));
			return result;
		}
	};

	template <typename VectorSystem = default_vector_system>
	struct ortho_projection_lhs
	{
		template <typename Ty>
		static basic_square_matrix_4d<Ty> get(
			const Ty &left_,
			const Ty &right_,
			const Ty &top_,
			const Ty &bottom_,
			const Ty &near_,
			const Ty &far_)
		{
			basic_square_matrix_4d<Ty> result;
			impl::set_c_v<VectorSystem>(result, 0, 0, static_cast<Ty>(2) / (right_ - left_));
			impl::set_c_v<VectorSystem>(result, 0, 3, -(right_ + left_) / (right_ - left_));
			impl::set_c_v<VectorSystem>(result, 1, 1, static_cast<Ty>(2) / (top_ - bottom_));
			impl::set_c_v<VectorSystem>(result, 1, 3, -(top_ + bottom_) / (top_ - bottom_));
			impl::set_c_v<VectorSystem>(result, 2, 3, -(far_ + near_) / (far_ - near_));
			impl::set_c_v<VectorSystem>(result, 3, 3, static_cast<Ty>(1));

			impl::set_c_v<VectorSystem>(result, 2, 2, static_cast<Ty>(2) / (far_ - near_));
			return result;
		}
	};

	template <typename Ty = double, typename VectorSystem = default_vector_system>
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

		template <typename Ty>
		transform& then(const basic_square_matrix_4d<Ty> &mat_)
		{
			_multi(mat_);
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

	namespace impl
	{
		template <typename Ty, std::size_t ...EularAngleBasis>
		struct to_quaternion_helper
		{
		private:
			using Quat = basic_quaternion<Ty>;

			using EularAngle = basic_eular_angle<Ty, EularAngleBasis...>;

			constexpr static std::size_t EularAngleBasisSize = sizeof...(EularAngleBasis);

			template <std::size_t N>
			struct Progress
			{
				static Quat get(const Quat &lastquat, const EularAngle &eularangle_)
				{
					return Progress<N + 1>::get(
						lastquat * get_cur_quat(eularangle_), eularangle_);
				}
			private:
				static Quat get_cur_quat(const EularAngle &eularangle_)
				{
					constexpr std::size_t coord =
						get_nth_param_value<N, std::size_t, EularAngleBasis...>::value;
					return rotation_quaternion(
						eularangle_[N], basic_vector_3d<Ty>::basis<coord>());
				}
			};

			template <>
			struct Progress<EularAngleBasisSize>
			{
				static Quat get(const Quat &lastquat, const EularAngle &eularangle_)
				{
					return lastquat;
				}
			};
		public:
			static Quat get(const EularAngle &eularangle_)
			{
				return Progress<0>::get(rotation_quaternion<Ty>(), eularangle_);
			}
		};
	}

	template <typename Ty, std::size_t ...EularAngleBasis>
	basic_quaternion<Ty> to_quaternion(
		const basic_eular_angle<Ty, EularAngleBasis...> &eularangle_)
	{
		return impl::to_quaternion_helper<Ty, EularAngleBasis...>::get(eularangle_);
	}
}