#pragma once

#include <linear_algebra/impl/_number_array.h>
#include <cassert>
#include <type_traits>

namespace lin
{
	/* TEMPLATE CLASS basic_vector
	The basic_vector represents a mathematical [Dimension] dimension vector with elements of type [Ty].
	It underlying use a plain array with [Dimension] elements of type [Ty].
	*/
	template <typename Ty, std::size_t Dimension>
	class basic_vector
		:public impl::_number_array<Ty, Dimension>
	{
		using MyBase = impl::_number_array<Ty, Dimension>;
	public:
		using MyBase::_number_array;

		// Large enough to hold dimension of this vector.
		using dimension_type = std::size_t;

		using value_type = typename MyBase::value_type;

		// Equal to [Dimension].
		constexpr static dimension_type dimension()
		{
			return Dimension;
		}

		// Return the square of length of this vector.
		value_type length_square() const
		{
			value_type result = 0;
			for (auto &comp : *this)
				result += comp * comp;
			return result;
		}

		// Return the length of this vector.
		value_type length() const
		{
			return std::sqrt(length_square());
		}

		// Let this vector be normalized vector of original one and return this vector.
		basic_vector& normalize()
		{
			return static_cast<basic_vector&>(*this /= length());
		}

		// Return a vector which is normalized vector of this vector.
		basic_vector get_normalized() const
		{
			return basic_vector(*this).normalize();
		}

		// Return a [Dimension - 1] dimensions vector
		// as if the last component is removed from this vector.
		basic_vector<Ty, Dimension - 1> reduce() const
		{
			static_assert(Dimension > 1,
				"Only vector with 2 or higher dimensions may perform the reduce function.");

			basic_vector<Ty, Dimension - 1> result;
			for (decltype(result.dimension()) i = 0; i < result.dimension(); ++i)
				result[i] = (*this)[i];
			return result;
		}

		// Return a [Dimension + 1] dimensions vector
		// as if an component of value [lastcomp] is append to this vector.
		basic_vector<Ty, Dimension + 1> homogeneous(Ty lastcomp = 1) const
		{
			static_assert(Dimension < std::numeric_limits<std::size_t>::max(),
				"Dimension overflow.");

			basic_vector<Ty, Dimension + 1> result;
			for (decltype(result.dimension()) i = 0; i < result.dimension() - 1; ++i)
				result[i] = (*this)[i];
			result[result.dimension() - 1] = lastcomp;
			return result;
		}

		basic_vector operator+() const
		{
			return MyBase::operator+();
		}

		basic_vector operator-() const
		{
			return MyBase::operator-();
		}

		template <typename RightTy>
		basic_vector& operator+=(const basic_vector<RightTy, Dimension> &right)
		{
			return static_cast<basic_vector&>(MyBase::operator+=(right));
		}

		template <typename RightTy>
		basic_vector& operator-=(const basic_vector<RightTy, Dimension> &right)
		{
			return static_cast<basic_vector&>(MyBase::operator-=(right));
		}

		template <typename RightTy>
		basic_vector<std::common_type_t<Ty, RightTy>, Dimension>
			operator+(const basic_vector<RightTy, Dimension> &right) const
		{
			return MyBase::operator+(right);
		}

		template <typename RightTy>
		basic_vector<std::common_type_t<Ty, RightTy>, Dimension>
			operator-(const basic_vector<RightTy, Dimension> &right) const
		{
			return MyBase::operator-(right);
		}

		template <typename RightTy>
		basic_vector& operator*=(const basic_vector<RightTy, Dimension> &right)
		{
			return static_cast<basic_vector&>(MyBase::operator*=(right));
		}

		template <typename RightTy>
		basic_vector& operator/=(const basic_vector<RightTy, Dimension> &right)
		{
			return static_cast<basic_vector&>(MyBase::operator/=(right));
		}

		template <typename RightTy>
		basic_vector<std::common_type_t<Ty, RightTy>, Dimension>
			operator*(const basic_vector<RightTy, Dimension> &right) const
		{
			return MyBase::operator*(right);
		}

		template <typename RightTy>
		basic_vector<std::common_type_t<Ty, RightTy>, Dimension>
			operator/(const basic_vector<RightTy, Dimension> &right) const
		{
			return MyBase::operator/(right);
		}

		basic_vector& operator*=(const value_type &right)
		{
			return static_cast<basic_vector&>(MyBase::operator*=(right));
		}

		basic_vector operator*(const value_type &right) const
		{
			return MyBase::operator*(right);
		}

		basic_vector& operator/=(const value_type &right)
		{
			return static_cast<basic_vector&>(MyBase::operator/=(right));
		}

		basic_vector operator/(const value_type &right) const
		{
			return MyBase::operator/(right);
		}

		// Return a vector whose [Comp]-th component has value 1 and other components are all 0.
		template <dimension_type Comp>
		static basic_vector basis()
		{
			basic_vector result{};
			result[Comp] = 1;
			return result;
		}
	};

	/* TEMPLATE FUNCTION cross_product
	Return a vector which is cross product(also called vector product) of
	3-dimension vectors [left] and [right].
	*/
	template <typename LeftTy, typename RightTy>
	basic_vector<std::common_type_t<LeftTy, RightTy>, 3> cross_product(
		const basic_vector<LeftTy, 3> &left,
		const basic_vector<RightTy, 3> &right)
	{
		return
		{
			left[1] * right[2] - left[2] * right[1],
			left[2] * right[0] - left[0] * right[2],
			left[0] * right[1] - left[1] * right[0],
		};
	}

	// d means dimension not double....

	template <typename Ty>
	using basic_vector_2d = basic_vector<Ty, 2>;

	template <typename Ty>
	using basic_vector_3d = basic_vector<Ty, 3>;

	template <typename Ty>
	using basic_vector_4d = basic_vector<Ty, 4>;

	using vector_2d = basic_vector_2d<double>;

	using vector_3d = basic_vector_3d<double>;

	using vector_4d = basic_vector_4d<double>;

	using vector_2f = basic_vector_2d<float>;

	using vector_3f = basic_vector_3d<float>;

	using vector_4f = basic_vector_4d<float>;

	inline namespace convention
	{
		/* CONSTEXPR VARIANBLES x_coord y_coord z_coord w_coord
		Convention signs for x, y, z, w denotes of vector components.
		*/
		constexpr std::size_t x_coord = 0, y_coord = 1, z_coord = 2, w_coord = 3;
	}
}

namespace std
{
	/* TEMPLATE FUNCTION inner_product
	Return the inner product(also called dot product) of vectors [left] and [right].
	*/
	template <typename LeftTy, typename RightTy, std::size_t Dimension>
	std::common_type_t<LeftTy, RightTy> inner_product(
		const lin::basic_vector<LeftTy, Dimension> &left,
		const lin::basic_vector<RightTy, Dimension> &right)
	{
		std::common_type_t<LeftTy, RightTy> result = 0;
		for (std::size_t i = 0; i != Dimension; ++i)
			result += (left[i] * right[i]);
		return result;
	}
}

namespace lin
{
	/* TEMPLATE FUNCTION angle_between
	Return the angle, in radius, of vectors [left] and [right].
	*/
	template <typename LeftTy, typename RightTy, std::size_t Dimension>
	std::common_type_t<LeftTy, RightTy> angle_between(
		const lin::basic_vector<LeftTy, Dimension> &left,
		const lin::basic_vector<RightTy, Dimension> &right)
	{
		using ResultTy = std::common_type_t<LeftTy, RightTy>;
		auto lenprod = left.length() * right.length();
		assert(lenprod != static_cast<ResultTy>(0) &&
			"Neither the first or the second operand of angle_between function can be zero vector.");
		return std::acos(std::inner_product(left, right) / lenprod);
	}
}