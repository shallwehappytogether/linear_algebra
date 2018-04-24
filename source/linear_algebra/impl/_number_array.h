#pragma once

#include <cstddef>
#include <array>

namespace lin::impl
{
	template <typename Ty, std::size_t Size>
	class _number_array_storage
		:public std::array<Ty, Size>
	{
	public:
		using std::array<Ty, Size>::array;
	};

	template <typename Ty>
	class _number_array_storage<Ty, 1>
	{
	public:
		_number_array_storage(const Ty &v = 0)
			:x(v)
		{
		};

		Ty* data()
		{
			return &x;
		}

		const Ty* data() const
		{
			return &x;
		}

		Ty x;
	};

	template <typename Ty>
	class _number_array_storage<Ty, 2>
	{
	public:
		_number_array_storage(const Ty &x_ = 0, const Ty &y_ = 0)
			:x(x_), y(y_)
		{
		};

		Ty* data()
		{
			return _array;
		}

		const Ty* data() const
		{
			return _array;
		}

		union
		{
			struct { Ty x, y; };
			Ty _array[2];
		};
	};

	template <typename Ty>
	class _number_array_storage<Ty, 3>
	{
	public:
		_number_array_storage(const Ty &x_ = 0, const Ty &y_ = 0, const Ty &z_ = 0)
			:x(x_), y(y_), z(z_)
		{
		};

		Ty* data()
		{
			return _array;
		}

		const Ty* data() const
		{
			return _array;
		}

		union
		{
			struct { Ty x, y, z; };
			Ty _array[3];
		};
	};

	template <typename Ty>
	class _number_array_storage<Ty, 4>
	{
	public:
		_number_array_storage(const Ty &x_ = 0, const Ty &y_ = 0, const Ty &z_ = 0, const Ty &w_ = 0)
			:x(x_), y(y_), z(z_), w(w_)
		{
		};

		Ty* data()
		{
			return _array;
		}

		const Ty* data() const
		{
			return _array;
		}

		union
		{
			struct { Ty x, y, z, w; };
			Ty _array[4];
		};
	};

	template <typename Ty, std::size_t Size>
	class _number_array_base
		:public _number_array_storage<Ty, Size>
	{
		using MyBase = _number_array_storage<Ty, Size>;
	public:
		typedef Ty value_type;

		typedef std::size_t size_type;

		using _number_array_storage<Ty, Size>::_number_array_storage;

		Ty & operator[](std::size_t i)
		{
			return MyBase::data()[i];
		}

		const Ty& operator[](std::size_t i) const
		{
			return MyBase::data()[i];
		}

		Ty* begin()
		{
			return MyBase::data();
		}

		const Ty* begin() const
		{
			return MyBase::data();
		}

		Ty* end()
		{
			return MyBase::data() + Size;
		}

		const Ty* end() const
		{
			return MyBase::data() + Size;
		}

		std::size_t size() const
		{
			return Size;
		}
	};

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
		:public _number_array_base<Ty, Size>
	{
		using MyBase = _number_array_base<Ty, Size>;
	public:
		using value_type = typename MyBase::value_type;

		using size_type = typename MyBase::size_type;

		template <typename... Args>
		_number_array(Args... args)
			:MyBase{ static_cast<Ty>(args)... }
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
			for (size_type i = 0; i < MyBase::size(); ++i)
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