#pragma once

#include "impl/_number_array.h"
#include "vector.h"
#include <type_traits>
#include <cassert>

namespace lin
{
	/* TEMPLATE CLASS basic_matrix
	The basic_matrix represents a mathematical matrix with [RowDimension] row and [ColumnDimension] column elements of type [Ty].
	It underlying use a plain array with [RowDimension] * [ColumnDimension] elements of type [Ty].

	Matrix's elements are stored sequenced, no padding, column by column.
	*/
	template <typename Ty, std::size_t RowDimension, std::size_t ColumnDimension>
	class basic_matrix
		:public impl::_number_array<Ty, RowDimension * ColumnDimension>
	{
		typedef impl::_number_array<Ty, RowDimension * ColumnDimension> MyBase;
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

		// Return the reference to element at row [r] and column [c].
		constexpr reference element_at(row_dimension_type r, column_dimension_type c)
		{
			return MyBase::operator[](_elemIndex(r, c));
		}

		// Return the const reference to element at row [r] and column [c].
		constexpr const_reference element_at(row_dimension_type r, column_dimension_type c) const
		{
			return MyBase::operator[](_elemIndex(r, c));
		}

		// Same as [element_at], but with boundary check.
		constexpr reference at(row_dimension_type r, column_dimension_type c)
		{
			return MyBase::at(_elemIndex(r, c));
		}

		// Same as [element_at], but with boundary check.
		constexpr const_reference at(row_dimension_type r, column_dimension_type c) const
		{
			return MyBase::at(_elemIndex(r, c));
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
		size_type _elemIndex(row_dimension_type rowdimension, column_dimension_type columndimension) const
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

		// 将本矩阵转置并返回。
		basic_square_matrix& transpose()
		{
			return *this = get_transpose();
		}

		// 返回本矩阵的转置矩阵。
		basic_square_matrix get_transpose()
		{
			basic_square_matrix result;
			for (dimension_type r = 0; r < dimension(); ++r)
				for (dimension_type c = 0; c < dimension(); ++c)
				{
					result.element_at(c, r) = element_at(r, c);
				}
			return result;
		}

		// 返回本矩阵的伴随矩阵。
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

		// 将本矩阵求逆并返回。
		basic_square_matrix& invert()
		{
			return *this = inverse();
		}

		// 返回本矩阵的逆矩阵。
		basic_square_matrix inverse() const
		{
			auto det = determinant();
			if (det == 0)
				return *this;
			auto result = adjoint();
			result /= det;
			return result;
		}

		// 返回一个新的矩阵，这个新矩阵对应的行列式是本矩阵对应的行列式的第[r]行第[c]列元素的代数余子式。
		basic_square_matrix<Ty, Dimension - 1>
			algebraic_cofactor(dimension_type r, dimension_type c) const
		{
			static_assert(Dimension > 0,
				"Only square matrix with non-zero dimension(s) may get the algebraic_cofactor.");

			basic_square_matrix<Ty, Dimension - 1> result;
			_cofactor(data(), dimension(), r, c, result.data());
			return result;
		}

		// 计算本矩阵对应的行列式的值。
		value_type determinant() const
		{
			if constexpr (Dimension == 1)
				return data()[0];
			else if constexpr (Dimension == 2)
				return _det2(data());
			else if constexpr (Dimension == 3)
				return _det3(data());
			else if constexpr (Dimension == 4)
				return _det4(data());
			else
				return _det(data(), Dimension);
		}
	private:
		
		// 存储[fromdim]阶行列式[from]的元素[dstr, dstc]的代数余子式到[result]中。
		static void _cofactor(
			const Ty *from,
			size_type fromdim,
			size_type dstr,
			size_type dstc,
			Ty *result)
		{
			assert(dstr < fromdim && dstc < fromdim);

			size_type r = 0, c = 0, iOutput = 0;

			value_type fac = (dstr + dstr) % 2 == 0 ? 1 : -1;

			auto writeElem = [&]() { result[iOutput++] = from[c * fromdim + r] * fac; };

			auto processColumn = [&]()
			{
				for (r = 0; r < dstr; ++r)
					writeElem();
				for (r = dstr + 1; r < fromdim; ++r)
					writeElem();
			};

			for (c = 0; c < dstc; ++c)
				processColumn();
			for (c = dstc + 1; c < fromdim; ++c)
				processColumn();
		}

		// 计算[dim]阶行列式[arr]的值。
		// 算法：行列式等于它的任一行（列）的各元素与其对应的代数余子式乘积之和。
		static inline value_type _det(
			Ty *arr, dimension_type dim)
		{
			assert(dim >= 4);

			if (dim == 4)
				return _det4(arr);

			auto elemAt = [=](auto r, auto c) { return arr[c * dim + r]; };

			Ty result = 0;
			for (dimension_type i = 0; i < dimension(); ++i)
			{ // process first row
				auto tmp = std::make_unique<Ty[]>((dim - 1) * (dim - 1));
				_cofactor(arr, dim, 0, i, tmp.data());
				auto v = _det(tmp.data());
				result += (i % 2 == 0 ? v : -v);
			}
			return result;
		}

		// 计算2阶行列式[arr]的值。
		static inline value_type _det2(
			const Ty *arr)
		{
			auto elemAt = [=](auto r, auto c) { return arr[c * 2 + r]; };
			return elemAt(0, 0) * elemAt(1, 1) - elemAt(0, 1) * elemAt(1, 0);
		}

		static inline value_type _det3Helper(
			value_type a, value_type b, value_type c,
			value_type d, value_type e, value_type f,
			value_type g, value_type h, value_type i)
		{
			return a * e * i + b * f * g + c * d * h -
				c * e * g - b * d * i - a * f * h;
		}

		// 计算3阶行列式[arr]的值。
		static inline value_type _det3(
			const Ty *arr)
		{
			auto elemAt = [=](auto r, auto c) { return arr[c * 3 + r]; };

			auto a = elemAt(0, 0);
			auto b = elemAt(0, 1);
			auto c = elemAt(0, 2);
			auto d = elemAt(1, 0);
			auto e = elemAt(1, 1);
			auto f = elemAt(1, 2);
			auto g = elemAt(2, 0);
			auto h = elemAt(2, 1);
			auto i = elemAt(2, 2);

			return _det3Helper(a, b, c, d, e, f, g, h, i);
		}

		// 计算4阶行列式[arr]的值。
		static inline value_type _det4(
			const Ty *arr)
		{
			auto elemAt = [=](auto r, auto c) { return arr[c * 4 + r]; };

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

			auto aterm = _det3Helper(f, g, h, j, k, l, n, o, p);
			auto bterm = _det3Helper(e, g, h, i, k, l, m, o, p);
			auto cterm = _det3Helper(e, f, h, i, j, l, m, n, p);
			auto dterm = _det3Helper(e, f, g, i, j, k, m, n, o);

			return a * aterm - b * bterm + c * cterm - d * dterm;
		}
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