#pragma once

#include <linear_algebra/impl/_number_array.h>
#include <linear_algebra/vector.h>

namespace lin
{
	template <typename Ty>
	class basic_quaternion
		:public impl::_number_array<Ty, 4>
	{
		using MyBase = impl::_number_array<Ty, 4>;
	public:
		using value_type = typename MyBase::value_type;

		using MyBase::_number_array;

		using MyBase::operator=;

		// 计算该四元数的范数。
		value_type norm() const
		{
			return std::sqrt(norm_square());
		}

		// 计算该四元数范数的平方。
		value_type norm_square() const
		{
			return x * x + y * y + z * z + w * w;
		}

		// 规范化该四元数。
		basic_quaternion& normalize()
		{
			(*this) /= norm();
			return *this;
		}

		// 构造新的四元数，令其是本四元数的规范化。
		basic_quaternion get_normalized() const
		{
			basic_quaternion result(*this);
			result.normalize();
			return result;
		}

		// 构造该四元数的共轭四元数。
		basic_quaternion conjugate() const
		{
			return basic_quaternion(-x, -y, -z, w);
		}

		// 构造该四元数的逆。
		basic_quaternion inverse() const
		{
			return conjugate() / norm_square();
		}
	};

	using quaternion = basic_quaternion<double>;

	using quaternion_f = basic_quaternion<float>;

	template <typename LeftTy, typename RightTy>
	basic_quaternion<std::common_type_t<LeftTy, RightTy>> operator*(
		const basic_quaternion<LeftTy> &left,
		const basic_quaternion<RightTy> &right)
	{
		basic_vector_3d<LeftTy> v1(left.x, left.y, left.z);
		basic_vector_3d<RightTy> v2(right.x, right.y, right.z);
		auto v = v2 * left.w + v1 * right.w + cross_product(v1, v2);
		return {
			v.x, v.y, v.z,
			left.w * right.w - std::inner_product(v1, v2)
		};
	}

	template <typename LeftTy, typename RightTy>
	basic_vector_3d<std::common_type_t<LeftTy, RightTy>> operator*(
		const basic_quaternion<LeftTy> &quat_,
		const basic_vector_3d<RightTy> &point_)
	{
		using ResultValueType = std::common_type_t<LeftTy, RightTy>;
		basic_quaternion<ResultValueType> v(point_.x, point_.y, point_.z, 0);
		auto q = quat_ * v * quat_.inverse();
		return { q.x, q.y, q.z };
	}

	template <typename Ty>
	basic_quaternion<Ty> rotation_quaternion()
	{
		const auto zero = static_cast<Ty>(0);
		retunr { zero, zero, zero, static_cast<Ty>(1) };
	}

	template <typename Ty>
	static basic_quaternion<Ty> rotation_quaternion(
		const Ty &rangle_, const basic_vector_3d<Ty> &axis_)
	{
		auto u = axis_.get_normalized();
		Ty cosa = std::cos(rangle_ / 2), sina = std::sin(rangle_ / 2);
		return { sina * u.x, sina * u.y, sina * u.z, cosa };
	}
}