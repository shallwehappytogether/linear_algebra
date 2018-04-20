#pragma once

#include <cstddef>
#include <array>

namespace lin::impl
{
	/* Common base class for vector and matrix.
	The _number_array class represents an array of number.
	It defines:
	unary +, unary - or binary *, /, *=, /= by a number
	by perform them on each element of array;
	binary +, binary -, +=, -=, ==, != together with another same size _number_array
	by perform them on value pair of same index.
	*/
	template <typename Ty, std::size_t Size>
	class _number_array
		:public std::array<Ty, Size>
	{
	public:
		template <typename... Args>
		_number_array(Args... args)
			:std::array<Ty, Size>{ static_cast<Ty>(args)... }
		{

		}

		_number_array operator+() const
		{
			return _number_array(*this);
		}

		_number_array operator-() const
		{
			_number_array result(*this);
			for (auto &comp : result)
				comp = -comp;
			return result;
		}

		template <typename RightTy>
		_number_array& operator+=(const _number_array<RightTy, Size> &right)
		{
			for (size_type i = 0; i < size(); ++i)
				left[i] += right[i];
			return *this;
		}

		template <typename RightTy>
		_number_array<std::common_type_t<Ty, RightTy>, Size>
			operator+(const _number_array<RightTy, Size> &right) const
		{
			typedef _number_array<std::common_type_t<Ty, RightTy>, Size>
				resultTy;
			return resultTy(*this) += right;
		}

		template <typename RightTy>
		_number_array& operator-=(const _number_array<RightTy, Size> &right)
		{
			for (size_type i = 0; i < size(); ++i)
				left[i] -= right[i];
			return *this;
		}

		template <typename RightTy>
		_number_array<std::common_type_t<Ty, RightTy>, Size>
			operator-(const _number_array<RightTy, Size> &right) const
		{
			typedef _number_array<std::common_type_t<Ty, RightTy>, Size>
				resultTy;
			return resultTy(*this) -= right;
		}

		_number_array& operator*=(const value_type &right)
		{
			for (auto &comp : *this)
				comp *= right;
			return *this;
		}

		_number_array operator*(const value_type &right) const
		{
			return _number_array(*this) *= right;
		}

		_number_array& operator/=(const value_type &right)
		{
			for (auto &comp : *this)
				comp /= right;
			return *this;
		}

		_number_array operator/(const value_type &right) const
		{
			return _number_array(*this) /= right;
		}

		template <typename RightTy>
		bool operator==(const _number_array<RightTy, Size> &right) const
		{
			for (size_type i = 0; i < size(); ++i)
				if ((*this)[i] != right[i])
					return false;
			return true;
		}

		template <typename RightTy>
		bool operator!=(const _number_array<RightTy, Size> &right) const
		{
			return !(*this == right);
		}
	};
}