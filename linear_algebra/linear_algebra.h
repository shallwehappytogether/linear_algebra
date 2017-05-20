#pragma once

#include <cstddef>
#include <array>
#include <cmath>
#include <numeric>
#include <cassert>
#include <memory>

namespace linear_algebra
{
	namespace impl
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
				_number_array(*this) += right;
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
				_number_array(*this) -= right;
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

		constexpr bool _column_first_storage = false;
	}

	template <typename Ty, std::size_t Dimension>
	class basic_vector
		:public impl::_number_array<Ty, Dimension>
	{
	public:
		using impl::_number_array<Ty, Dimension>::_number_array;

		typedef std::size_t dimension_type;

		constexpr dimension_type dimension() const
		{
			return Dimension;
		}

		value_type length_square() const
		{
			value_type result = 0;
			for (auto &comp : *this)
				result += comp * comp;
			return result;
		}

		value_type length() const
		{
			return std::sqrt(length_square());
		}

		basic_vector& normalize_to_assign()
		{
			return static_cast<basic_vector&>(*this /= length());
		}

		basic_vector normalize() const
		{
			return basic_vector(*this).normalize_to_assign();
		}

		basic_vector<Ty, Dimension - 1> reduce() const
		{
			static_assert(Dimension > 1,
				"Only vector with 2 or higher dimensions may perform the reduce function.");

			basic_vector<Ty, Dimension - 1> result;
			for (decltype(result.dimension()) i = 0; i < result.dimension(); ++i)
				result[i] = (*this)[i];
			return result;
		}

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
	};

	template <typename Ty, std::size_t Dimension>
	Ty angle_between(
		const linear_algebra::basic_vector<Ty, Dimension> &left,
		const linear_algebra::basic_vector<Ty, Dimension> &right)
	{
		auto lenprod = left.length() * right.length();
		assert(lenprod != static_cast<Ty>(0) &&
			"Neither the first or the second operand of angle_between function can be zero vector.");
		return std::acos(std::inner_product(left, right) / lenprod);
	}

	template <typename Ty>
	basic_vector<Ty, 3> cross_product(
		const basic_vector<Ty, 3> &left,
		const basic_vector<Ty, 3> &right)
	{
		return
		{
			left[1] * right[2] - left[2] * right[1],
			left[2] * right[0] - left[0] * right[2],
			left[0] * right[1] - left[1] * right[0],
		};
	}

	template <typename Ty>
	using basic_vector_2d = basic_vector<Ty, 2>;

	template <typename Ty>
	using basic_vector_3d = basic_vector<Ty, 3>;

	template <typename Ty>
	using basic_vector_4d = basic_vector<Ty, 4>;

	using vector_2d = basic_vector_2d<double>;

	using vector_3d = basic_vector_3d<double>;

	using vector_4d = basic_vector_4d<double>;

	template <typename Ty, std::size_t RowDimension, std::size_t ColumnDimension>
	class basic_matrix
		:public impl::_number_array<Ty, RowDimension * ColumnDimension>
	{
		typedef impl::_number_array<Ty, RowDimension * ColumnDimension> MyBase;
	protected:
		using MyBase::operator[];

		using MyBase::at;

		using MyBase::begin;

		using MyBase::end;
	public:
		using MyBase::_number_array;

		typedef std::size_t row_dimension_type;

		typedef std::size_t column_dimension_type;

		constexpr row_dimension_type row_dimension() const
		{
			return RowDimension;
		}

		constexpr column_dimension_type column_dimension() const
		{
			return ColumnDimension;
		}

		reference element_at(row_dimension_type rowdimension, column_dimension_type columndimension)
		{
			return MyBase::operator[](_elemIndexAt(rowdimension, columndimension));
		}

		const_reference element_at(row_dimension_type rowdimension, column_dimension_type columndimension) const
		{
			return MyBase::operator[](_elemIndexAt(rowdimension, columndimension));
		}

		reference at(row_dimension_type rowdimension, column_dimension_type columndimension)
		{
			return MyBase::at(_elemIndexAt(rowdimension, columndimension));
		}

		const_reference at(row_dimension_type rowdimension, column_dimension_type columndimension) const
		{
			return MyBase::at(_elemIndexAt(rowdimension, columndimension));
		}

		basic_matrix<Ty, ColumnDimension, RowDimension> transpose() const
		{
			basic_matrix<Ty, ColumnDimension, RowDimension> result;
			for (row_dimension_type r = 0; r < row_dimension(); ++r)
				for (column_dimension_type c = 0; c < column_dimension(); ++c)
					result.element_at(c, r) = this->element_at(r, c);
			return result;
		}

		template <typename RightTy, std::size_t RightColumnDimension>
		basic_matrix<std::common_type_t<Ty, RightTy>, ColumnDimension, RightColumnDimension>
			operator*(const basic_matrix<RightTy, RowDimension, RightColumnDimension> &right) const
		{
			typedef basic_matrix<std::common_type_t<Ty, RightTy>, ColumnDimension, RightColumnDimension>
				resultTy;
			resultTy result;
			for (decltype(result.row_dimension()) r = 0; r < result.row_dimension(); ++r)
				for (decltype(result.column_dimension()) c = 0; c < result.column_dimension(); ++c)
				{
					typename resultTy::value_type v = 0;
					for (column_dimension_type i = 0; i < column_dimension(); ++i)
						v += this->element_at(r, i) * right.element_at(i, c);
					result.element_at(r, c) = v;
				}
			return result;
		}
	private:
		template <typename RowFirstStorage = void>
		size_type _elemIndexAt(row_dimension_type rowdimension, column_dimension_type columndimension) const
		{
			return rowdimension * ColumnDimension + columndimension;
		}

		template <>
		size_type _elemIndexAt<std::enable_if_t<impl::_column_first_storage>>(
			row_dimension_type rowdimension, column_dimension_type columndimension) const
		{
			return columndimension * RowDimension + rowdimension;
		}
	};

	template <typename Ty, std::size_t Dimension>
	class basic_square_matrix
		:public basic_matrix<Ty, Dimension, Dimension>
	{
		typedef basic_matrix<Ty, Dimension, Dimension> MyBase;
	public:
		using MyBase::basic_matrix;

		typedef row_dimension_type dimension_type;

		constexpr static dimension_type dimension()
		{
			return Dimension;
		}

		basic_square_matrix& transpose_to_assign()
		{
			return *this = transpose();
		}

		basic_square_matrix& adjoint_to_assign()
		{
			return *this = adjoint();
		}

		template <typename Enable = void>
		basic_square_matrix adjoint() const
		{
			basic_square_matrix result;
			for (dimension_type r = 0; r < dimension(); ++r)
				for (dimension_type c = 0; c < dimension(); ++c)
				{
					auto detFac = algebraic_cofactor(r, c).determinant();
					result.element_at(c, r) = detFac;
				}
			return result;
		}

		basic_square_matrix& inverse_to_assign()
		{
			return *this = inverse();
		}

		basic_square_matrix inverse() const
		{
			auto det = determinant();
			if (det == 0)
				return *this;
			auto result = adjoint();
			result /= det;
			return result;
		}

		basic_square_matrix<Ty, Dimension - 1>
			algebraic_cofactor(dimension_type r, dimension_type c) const
		{
			static_assert(Dimension > 0,
				"Only square matrix with non-zero dimension(s) may get the algebraic_cofactor.");
			basic_square_matrix<Ty, Dimension - 1> result;
			_saveCofactor(data(), dimension(), r, c, result.data());
			return result;
		}

		template <typename HigherDimension = void>
		value_type determinant() const
		{
			return _det_helper<Dimension>::get(*this);
		}

	private:
		template <typename RowFirstStorage = void>
		static void _saveCofactor(
			const Ty *from,
			size_type fromdim,
			size_type r,
			size_type c,
			Ty *result)
		{
			size_type iOutput = 0;
			for (size_type i = 0; i < r; ++i)
			{
				for (size_type j = 0; j < c; ++j)
					result[iOutput++] = from[i * fromdim + j];
				if (c != std::numeric_limits<size_type>::max() - 1)
					for (size_type j = c + 1; j < fromdim; ++j)
						result[iOutput++] = from[i * fromdim + j];
			}
			if (r != std::numeric_limits<size_type>::max())
				for (size_type i = r + 1; i < fromdim; ++i)
				{
					for (size_type j = 0; j < c; ++j)
						result[iOutput++] = from[i * fromdim + j];
					if (c != std::numeric_limits<size_type>::max() - 1)
						for (size_type j = c + 1; j < fromdim; ++j)
							result[iOutput++] = from[i * fromdim + j];
				}
		}

		template <>
		static void _saveCofactor<std::enable_if_t<impl::_column_first_storage>>(
			const Ty *from,
			size_type fromdim,
			size_type r,
			size_type c,
			Ty *result)
		{
			size_type iOutput = 0;
			for (size_type i = 0; i < c; ++i)
			{
				for (size_type j = 0; j < r; ++j)
					result[iOutput++] = from[j * fromdim + i];
				if (r != std::numeric_limits<size_type>::max() - 1)
					for (size_type j = r + 1; j < fromdim; ++j)
						result[iOutput++] = from[j * fromdim + i];
			}
			if (c != std::numeric_limits<size_type>::max())
				for (size_type i = c + 1; i < fromdim; ++i)
				{
					for (size_type j = 0; j < r; ++j)
						result[iOutput++] = from[j * fromdim + i];
					if (r != std::numeric_limits<size_type>::max() - 1)
						for (size_type j = r + 1; j < fromdim; ++j)
							result[iOutput++] = from[j * fromdim + i];
				}
		}

		template <typename RowFirstStorage = void>
		static inline Ty _getElment(const Ty *arr, dimension_type dim, dimension_type r, dimension_type c)
		{
			return arr[r * dim + c];
		}

		template <>
		static inline Ty _getElment<std::enable_if_t<impl::_column_first_storage>>(
			const Ty *arr, dimension_type dim, dimension_type r, dimension_type c)
		{
			return arr[c * dim + r];
		}

		static inline value_type _5dOrHigerDet(
			Ty *arr, dimension_type dim)
		{
			assert(dim >= 4);

			if (dim == 4)
				return _4dDet(arr);

			auto elemAt = [=](auto r, auto c) { return _getElment(arr, dim, r, c); };
			Ty result = 0;
			for (dimension_type i = 0; i < dimension(); ++i)
			{
				auto tmp = std::make_unique<Ty[]>((dim - 1) * (dim - 1));
				_saveCofactor(arr, dim, 0, i, tmp.data());
				auto v = _5dOrHigerDet(tmp.data());
				result += (i % 2 == 0 ? v : -v);
			}

			return result;
		}

		static inline value_type _2dDet(
			const Ty *arr)
		{
			auto elemAt = [=](auto r, auto c) { return _getElment(arr, 2, r, c); };

			return elemAt(0, 0) * elemAt(1, 1) -
				elemAt(0, 1) * elemAt(1, 0);
		}

		static inline value_type _3dDetHelper(
			value_type a, value_type b, value_type c,
			value_type d, value_type e, value_type f,
			value_type g, value_type h, value_type i)
		{
			return a * e * i + b * f * g + c * d * h -
				c * e * g - b * d * i - a * f * h;
		}

		static inline value_type _3dDet(
			const Ty *arr)
		{
			auto elemAt = [=](auto r, auto c) { return _getElment(arr, 3, r, c); };

			auto a = elemAt(0, 0);
			auto b = elemAt(0, 1);
			auto c = elemAt(0, 2);
			auto d = elemAt(1, 0);
			auto e = elemAt(1, 1);
			auto f = elemAt(1, 2);
			auto g = elemAt(2, 0);
			auto h = elemAt(2, 1);
			auto i = elemAt(2, 2);

			return _3dDetHelper(a, b, c, d, e, f, g, h, i);
		}

		static inline value_type _4dDet(
			const Ty *arr)
		{
			auto elemAt = [=] (auto r, auto c) { return _getElment(arr, 4, r, c); };

			auto a = elemAt(0, 0);
			auto b = elemAt(0, 1);
			auto c = elemAt(0, 2);
			auto d = elemAt(0, 3);
			auto e = elemAt(1, 0);
			auto f = elemAt(1, 1);
			auto g = elemAt(1, 2);
			auto h = elemAt(1, 3);
			auto i = elemAt(2, 0);
			auto j = elemAt(2, 1);
			auto k = elemAt(2, 2);
			auto l = elemAt(2, 3);
			auto m = elemAt(3, 0);
			auto n = elemAt(3, 1);
			auto o = elemAt(3, 2);
			auto p = elemAt(3, 3);

			auto aterm = _3dDetHelper(f, g, h, j, k, l, n, o, p);
			auto bterm = _3dDetHelper(e, g, h, i, k, l, m, o, p);
			auto cterm = _3dDetHelper(e, f, h, i, j, l, m, n, p);
			auto dterm = _3dDetHelper(e, f, g, i, j, k, m, n, o);

			return a * aterm - b * bterm + c * cterm - d * dterm;
		}

		template <std::size_t Dimen>
		struct _det_helper
		{
			static value_type get(
				const basic_square_matrix &m)
			{
				return _5dOrHigerDet(m.data(), m.dimension());
			}
		};

		template <>
		struct _det_helper<1>
		{
			static value_type get(
				const basic_square_matrix &m)
			{
				return m.data()[0];
			}
		};

		template <>
		struct _det_helper<2>
		{
			static value_type get(
				const basic_square_matrix &m)
			{
				return _2dDet(m.data());
			}
		};

		template <>
		struct _det_helper<3>
		{
			static value_type get(
				const basic_square_matrix &m)
			{
				return _3dDet(m.data());
			}
		};

		template <>
		struct _det_helper<4>
		{
			static value_type get(
				const basic_square_matrix &m)
			{
				return _4dDet(m.data());
			}
		};
	};

	template <typename Ty>
	using basic_square_matrix_2d = basic_square_matrix<Ty, 2>;

	template <typename Ty>
	using basic_square_matrix_3d = basic_square_matrix<Ty, 3>;

	template <typename Ty>
	using basic_square_matrix_4d = basic_square_matrix<Ty, 4>;

	template <std::size_t Dimension>
	using square_matrix = basic_square_matrix<double, Dimension>;

	using square_matrix_2d = basic_square_matrix_2d<double>;

	using square_matrix_3d = basic_square_matrix_3d<double>;

	using square_matrix_4d = basic_square_matrix_4d<double>;

	template <typename Ty, std::size_t Dimension>
	basic_square_matrix<Ty, Dimension>
		make_identity_matrix()
	{
		basic_square_matrix<Ty, Dimension> result;
		for (decltype(result.row_dimension()) i = 0; i < result.row_dimension(); ++i)
			result.element_at(i, i) = static_cast<Ty>(1);
		return result;
	}

	template <typename VectorValueTy, typename MatrixTy,
		std::size_t VectorDimension, std::size_t MatrixColumnDimension>
	basic_vector<VectorValueTy, VectorDimension>&
		operator*=(
			const basic_matrix<MatrixTy, VectorDimension, MatrixColumnDimension> &mat,
			basic_vector<VectorValueTy, VectorDimension> &vec)
	{
		for (decltype(mat.row_dimension()) r = 0; r < VectorDimension; ++r)
		{
			std::common_type_t<VectorValueTy, MatrixTy> v = 0;
			for (decltype(mat.column_dimension()) c = 0; c < MatrixColumnDimension; ++c)
				v += (vec[c] * mat.element_at(r, c));
			vec[r] = v;
		}
		return vec;
	}

	template <typename VectorValueTy, typename MatrixTy,
		std::size_t VectorDimension, std::size_t MatrixColumnDimension>
	basic_vector<std::common_type_t<VectorValueTy, MatrixTy>, VectorDimension>
		operator*(
			const basic_matrix<MatrixTy, VectorDimension, MatrixColumnDimension> &mat,
			const basic_vector<VectorValueTy, VectorDimension> &vec)
	{
		typedef
			basic_vector<std::common_type_t<VectorValueTy, MatrixTy>, VectorDimension>
			resultType;
		return mat *= resultType(vec);
	}

	template <typename VectorValueTy, typename MatrixTy,
		std::size_t VectorDimension, std::size_t MatrixRowDimension>
	basic_vector<VectorValueTy, VectorDimension>
		operator*=(
			basic_vector<VectorValueTy, VectorDimension> &vec,
			const basic_matrix<MatrixTy, MatrixRowDimension, VectorDimension> &mat)
	{
		for (decltype(mat.column_dimension()) c = 0; c < VectorDimension; ++c)
		{
			std::common_type_t<VectorValueTy, MatrixTy> v = 0;
			for (decltype(mat.row_dimension()) r = 0; r < MatrixRowDimension; ++r)
				v += (vec[r] * mat.element_at(r, c));
			vec[c] = v;
		}
		return vec;
	}

	template <typename VectorValueTy, typename MatrixTy,
		std::size_t VectorDimension, std::size_t MatrixRowDimension>
	basic_vector<std::common_type_t<VectorValueTy, MatrixTy>, VectorDimension>
		operator*(
			const basic_vector<VectorValueTy, VectorDimension> &vec,
			const basic_matrix<MatrixTy, MatrixRowDimension, VectorDimension> &mat)
	{
		typedef
			basic_vector<std::common_type_t<VectorValueTy, MatrixTy>, VectorDimension>
			resultType;
		return resultType(vec) *= mat;
	}
}

namespace std
{
	template <typename Ty, std::size_t Dimension>
	Ty inner_product(
		const linear_algebra::basic_vector<Ty, Dimension> &left,
		const linear_algebra::basic_vector<Ty, Dimension> &right)
	{
		Ty result = 0;
		for (std::size_t i = 0; i != Dimension; ++i)
			result += (left[i] * right[i]);
		return result;
	}
}