#pragma once

#include <cstddef>
#include <array>
#include <cmath>
#include <numeric>
#include <cassert>
#include <memory>

namespace lin
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

		constexpr bool _column_first_storage = false;


	}

	/* TEMPLATE CLASS basic_vector
		The basic_vector represents a mathematical [Dimension] dimension vector with elements of type [Ty].
		It underlying use a plain array with [Dimension] elements of type [Ty].
	*/
	template <typename Ty, std::size_t Dimension>
	class basic_vector
		:public impl::_number_array<Ty, Dimension>
	{
	public:
		using impl::_number_array<Ty, Dimension>::_number_array;

		// Large enough to hold dimension of this vector.
		typedef std::size_t dimension_type;

		// Equal to [Dimension].
		constexpr dimension_type dimension() const
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
		basic_vector& normalize_to_assign()
		{
			return static_cast<basic_vector&>(*this /= length());
		}

		// Return a vector which is normalized vector of this vector.
		basic_vector normalize() const
		{
			return basic_vector(*this).normalize_to_assign();
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

		class
		{
		public:
			Ty& operator= (const Ty &v) { return (*this)[0] = v; }
			operator Ty() const { return (*this)[0]; }
		} x;

		class
		{
		public:
			Ty& operator= (const Ty &v) { return (*this)[1] = v; }
			operator Ty() const { return (*this)[1]; }
		} y;

		class
		{
		public:
			Ty& operator= (const Ty &v) { return (*this)[2] = v; }
			operator Ty() const { return (*this)[2]; }
		} z;

		class
		{
		public:
			Ty& operator= (const Ty &v) { return (*this)[3] = v; }
			operator Ty() const { return (*this)[3]; }
		} w;

		// Return a vector whose [Comp]-th component has value 1 and other components are all 0.
		template <dimension_type Comp>
		static basic_vector basis()
		{
			basic_vector result{};
			result[Comp] = 1;
			return result;
		}
	};

	/* TEMPLATE FUNCTION angle_between
		Return the angle, in radius, of vectors [left] and [right].
	*/
	template <typename LeftTy, typename RightTy, std::size_t Dimension>
	std::common_type_t<LeftTy, RightTy> angle_between(
		const lin::basic_vector<LeftTy, Dimension> &left,
		const lin::basic_vector<RightTy, Dimension> &right)
	{
		auto lenprod = left.length() * right.length();
		assert(lenprod != static_cast<Ty>(0) &&
			"Neither the first or the second operand of angle_between function can be zero vector.");
		return std::acos(std::inner_product(left, right) / lenprod);
	}

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

	template <typename Ty>
	using basic_vector_2d = basic_vector<Ty, 2>;

	template <typename Ty>
	using basic_vector_3d = basic_vector<Ty, 3>;

	template <typename Ty>
	using basic_vector_4d = basic_vector<Ty, 4>;

	using vector_2d = basic_vector_2d<double>;

	using vector_3d = basic_vector_3d<double>;

	using vector_4d = basic_vector_4d<double>;

	inline namespace convention
	{
		/* CONSTEXPR VARIANBLES x_coord y_coord z_coord w_coord
			Convention signs for x, y, z, w denotes of vector components.
		*/
		constexpr std::size_t x_coord = 0, y_coord = 1, z_coord = 2, w_coord = 3;
	}

	/* TEMPLATE CLASS basic_matrix
		The basic_matrix represents a mathematical matrix with [RowDimension] row and [ColumnDimension] column elements of type [Ty].
		It underlying use a plain array with [RowDimension] * [ColumnDimension] elements of type [Ty].

		Matrix's elements are stored sequenced, no padding,
		either row by row or column by column, which is implement-defined
		and user shall not concern.

		An element at row 'r' and column 'c' can also be located
		using a pair of numeric value ('f', 's').
		The 'f'('s') is defined to be 'r'('c'), if
		elements in matrix are stored row by row, or be 'c'('r') otherwise.
		The 'f' and 's' is called the first and secondnary seq-index of this element.

		An Element with seq-index 'f' and 's' is the (f * secondnar_dimension() + s)-th element of underlying array.
	*/
	template <typename Ty, std::size_t RowDimension, std::size_t ColumnDimension>
	class basic_matrix
		:public impl::_number_array<Ty, RowDimension * ColumnDimension>
	{
		typedef impl::_number_array<Ty, RowDimension * ColumnDimension> MyBase;
		constexpr static bool _rowFirstStorage = true;
	public:
		using MyBase::_number_array;

		// Large enough to hold row dimension of this matrix.
		typedef std::size_t row_dimension_type;

		// Large enough to hold column dimension of this matrix.
		typedef std::size_t column_dimension_type;

		// Equal to [RowDimension].
		constexpr static row_dimension_type row_dimension()
		{
			return RowDimension;
		}

		// Equal to [ColumnDimension].
		constexpr static column_dimension_type column_dimension()
		{
			return ColumnDimension;
		}

		// Large enough to hold any first seq-index,
		// or the one beyond the past of the largest first seq-index of this matrix.
		typedef std::conditional_t<_rowFirstStorage, row_dimension_type, column_dimension_type>
			first_dimension_type;

		// Large enough to hold any secondary seq-index,
		// or the one beyond the past of the largest secondary seq-index of this matrix.
		typedef std::conditional_t<_rowFirstStorage, column_dimension_type, row_dimension_type>
			secondary_dimension_type;

		// The largest first seq-index.
		constexpr static row_dimension_type first_dimension()
		{
			return _rowFirstStorage ?
				row_dimension() : column_dimension();
		}

		// The largest secondary seq-index.
		constexpr static column_dimension_type secondary_dimension()
		{
			return _rowFirstStorage ?
				column_dimension() : row_dimension();
		}

		// Return the seq-index of element at row [r] and column [c].
		static std::pair<first_dimension_type, secondary_dimension_type>
			select_priority(row_dimension_type r, column_dimension_type c)
		{
			return _rowFirstStorage ? { r, c } : { c, r };
		}

		// Return the reference to element at row [r] and column [c].
		constexpr reference element_at(row_dimension_type r, column_dimension_type c)
		{
			return MyBase::operator[](_elemIndexAt(r, c));
		}

		// Return the const reference to element at row [r] and column [c].
		constexpr const_reference element_at(row_dimension_type r, column_dimension_type c) const
		{
			return MyBase::operator[](_elemIndexAt(r, c));
		}

		// Return the reference to element of seq-index [firstdim] and [secdim].
		constexpr reference element_at_seq(first_dimension_type f, secondary_dimension_type s)
		{
			return (*this)[f * secondary_dimension() + s];
		}

		// Return the const reference to element of seq-index [firstdim] and [secdim].
		constexpr const_reference element_at_seq(first_dimension_type f, secondary_dimension_type s) const
		{
			return (*this)[f * secondary_dimension() + s];
		}

		// Same as [element_at], but with boundary check.
		constexpr reference at(row_dimension_type r, column_dimension_type c)
		{
			return MyBase::at(_elemIndexAt(r, c));
		}

		// Same as [element_at], but with boundary check.
		constexpr const_reference at(row_dimension_type r, column_dimension_type c) const
		{
			return MyBase::at(_elemIndexAt(r, c));
		}

		// Return a matrix which is transpose matrix of this matrix.
		basic_matrix<Ty, ColumnDimension, RowDimension> transpose() const
		{
			basic_matrix<Ty, ColumnDimension, RowDimension> result;
			for (row_dimension_type r = 0; r < row_dimension(); ++r)
				for (column_dimension_type c = 0; c < column_dimension(); ++c)
					result.element_at(c, r) = this->element_at(r, c);
			return result;
		}

		// Return a matrix which is the result of multiply this matrix by matrix [right];
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

	/* TEMPLATE CLASS basic_square_matrix
		The basic_square_matrix represents a mathematical square matrix with elements of type [Ty].
		It is essential a derived class of basic_matrix
		with [RowDimension] and [ColumnDimension] both equal to [Dimension].
	*/
	template <typename Ty, std::size_t Dimension>
	class basic_square_matrix
		:public basic_matrix<Ty, Dimension, Dimension>
	{
		typedef basic_matrix<Ty, Dimension, Dimension> MyBase;
	public:
		using MyBase::basic_matrix;

		// Large enough to hold row dimension or column dimension of this matrix.
		typedef row_dimension_type dimension_type;

		// Equal to [Dimension].
		constexpr static dimension_type dimension()
		{
			return Dimension;
		}

		// Let this matrix be transpose matrix of the orignal one and return this matrix.
		basic_square_matrix& transpose_to_assign()
		{
			return *this = transpose();
		}

		// Let this matrix be adjoint matrix of the orignal one and return this matrix.
		basic_square_matrix& adjoint_to_assign()
		{
			return *this = adjoint();
		}

		// Return a matrix which is adjoint matrix of this matrix.
		basic_square_matrix adjoint() const
		{
			basic_square_matrix result;
			for (dimension_type r = 0; r < dimension(); ++r)
				for (dimension_type c = 0; c < dimension(); ++c)
				{
					auto detFac = algebraic_cofactor_seq(r, c).determinant();
					result.element_at_seq(c, r) = detFac;
				}
			return result;
		}

		// Let this matrix be inverse matrix of the orignal one and return this matrix.
		basic_square_matrix& inverse_to_assign()
		{
			return *this = inverse();
		}

		// Return a matrix which is inverse matrix of this matrix.
		basic_square_matrix inverse() const
		{
			auto det = determinant();
			if (det == 0)
				return *this;
			auto result = adjoint();
			result /= det;
			return result;
		}

		// Return a matrix which is algebraic cofactor at row [r] and column [c] of the matrix.
		basic_square_matrix<Ty, Dimension - 1>
			algebraic_cofactor(dimension_type r, dimension_type c) const
		{
			static_assert(Dimension > 0,
				"Only square matrix with non-zero dimension(s) may get the algebraic_cofactor.");

			typedef basic_square_matrix<Ty, Dimension - 1>
				reusultTy;
			auto seqdim = reusultTy::select_priority(r, c);
			return algebraic_cofactor_seq(seqdim.first, seqdim.second);
		}

		// Return a matrix which is algebraic cofactor at seq-index [f] and [s] of the matrix.
		basic_square_matrix<Ty, Dimension - 1>
			algebraic_cofactor_seq(dimension_type f, dimension_type s) const
		{
			static_assert(Dimension > 0,
				"Only square matrix with non-zero dimension(s) may get the algebraic_cofactor.");

			basic_square_matrix<Ty, Dimension - 1> result;
			_saveCofactorSeq(data(), dimension(), f, s, result.data());
			return result;
		}

		// Return the determinant of the matrix.
		value_type determinant() const
		{
			return _det_helper<Dimension>::get(*this);
		}
	private:
		static void _saveCofactorSeq(
			const Ty *from,
			size_type fromdim,
			size_type first,
			size_type second,
			Ty *result)
		{
			auto elemAt = [=](auto r, auto c) { return from[r * fromdim + c]; };

			size_type iOutput = 0;
			for (size_type i = 0; i < first; ++i)
			{
				for (size_type j = 0; j < second; ++j)
					result[iOutput++] = elemAt(i, j);
				if (second != std::numeric_limits<size_type>::max() - 1)
					for (size_type j = second + 1; j < fromdim; ++j)
						result[iOutput++] = elemAt(i, j);
			}
			if (first != std::numeric_limits<size_type>::max())
				for (size_type i = first + 1; i < fromdim; ++i)
				{
					for (size_type j = 0; j < second; ++j)
						result[iOutput++] = elemAt(i, j);
					if (second != std::numeric_limits<size_type>::max() - 1)
						for (size_type j = second + 1; j < fromdim; ++j)
							result[iOutput++] = elemAt(i, j);
				}
		}

		static inline value_type _5dOrHigerDet(
			Ty *arr, dimension_type dim)
		{
			assert(dim >= 4);

			if (dim == 4)
				return _4dDet(arr);

			auto elemAt = [=](auto r, auto c) { return arr[r * dim + c]; };
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
			auto elemAt = [=](auto r, auto c) { return arr[r * 2 + c]; };

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
			auto elemAt = [=](auto r, auto c) { return arr[r * 3 + c]; };

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
			auto elemAt = [=](auto r, auto c) { return arr[r * 4 + c]; };

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

	/* TEMPLATE FUNCTION make_identity_matrix
		Return an identity [Dimension]x[Dimension] matrix with elements of [Ty].
	*/
	template <typename Ty, std::size_t Dimension>
	basic_square_matrix<Ty, Dimension>
		make_identity_matrix()
	{
		basic_square_matrix<Ty, Dimension> result;
		for (decltype(result.dimension()) i = 0; i < result.dimension(); ++i)
			result.element_at(i, i) = 1;
		return result;
	}

	/* TEMPLATE FUNCTION matrix *= vector
		Treate [vec] as a row vector.
		Let [vec] be the result row vector of multiply [mat] by original [vec].
	*/
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

	/* TEMPLATE FUNCTION matrix * vector
		Treate [vec] as a row vector.
		Return the result row vector of multiply [mat] by [vec].
	*/
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

	/* TEMPLATE FUNCTION vector *= matrix
		Treate [vec] as a column vector.
		Let [vec] be the result column vector of multiply original [vec] by [mat].
	*/
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

	/* TEMPLATE FUNCTION vector * matrix
		Treate [vec] as a column vector.
		Return the result column vector of multiply [vec] by [mat].
	*/
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