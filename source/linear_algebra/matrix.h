#pragma once

#include "impl/_number_array.h"
#include "vector.h"
#include <type_traits>
#include <cassert>
#include <memory>

namespace lin
{
	/* TEMPLATE CLASS basic_matrix_base
	The basic_matrix_base represents a mathematical matrix with [RowDimension] row and [ColumnDimension] column elements of type [Ty].
	It underlying use a plain array with [RowDimension] * [ColumnDimension] elements of type [Ty].

	Matrix's elements are stored sequenced, no padding, column by column.
	*/
	template <typename Ty, std::size_t RowDimension, std::size_t ColumnDimension>
	class basic_matrix
		:public impl::_number_array<Ty, RowDimension * ColumnDimension>
	{
		typedef impl::_number_array<Ty, RowDimension * ColumnDimension> MyBase;
	protected:
		using size_type = MyBase::size_type;
	public:
		using MyBase::_number_array;

		using MyBase::operator=;

		// Large enough to hold row dimension of this matrix.
		using row_dimension_type = std::size_t;

		// Large enough to hold column dimension of this matrix.
		using column_dimension_type = std::size_t;

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

		constexpr static bool is_squared()
		{
			return RowDimension == ColumnDimension;
		}

		// Return the reference to element at row [r] and column [c].
		constexpr Ty& element_at(row_dimension_type r, column_dimension_type c)
		{
			return (*this)[_elemIndex(r, c)];
		}

		// Return the const reference to element at row [r] and column [c].
		constexpr const Ty& element_at(row_dimension_type r, column_dimension_type c) const
		{
			return (*this)[_elemIndex(r, c)];
		}

		// Same as [element_at], but with boundary check.
		constexpr Ty& at(row_dimension_type r, column_dimension_type c)
		{
			return MyBase::at(_elemIndex(r, c));
		}

		// Same as [element_at], but with boundary check.
		constexpr const Ty& at(row_dimension_type r, column_dimension_type c) const
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

		// 返回本矩阵的转置矩阵。
		basic_matrix<Ty, ColumnDimension, RowDimension> get_transpose()
		{
			basic_matrix<Ty, ColumnDimension, RowDimension> result;
			for (row_dimension_type r = 0; r < row_dimension(); ++r)
				for (column_dimension_type c = 0; c < column_dimension(); ++c)
					result.element_at(c, r) = element_at(r, c);
			return result;
		}

		// 将本矩阵转置并返回。
		basic_matrix& transpose()
		{
			static_assert(is_squared(), "Worked only for squared matrix.");

			return *this = get_transpose();
		}

		// 返回本矩阵的伴随矩阵。
		basic_matrix adjoint() const
		{
			static_assert(is_squared(), "Worked only for squared matrix.");

			basic_matrix result;
			for (row_dimension_type r = 0; r < row_dimension(); ++r)
				for (row_dimension_type c = 0; c < column_dimension(); ++c)
				{
					auto detFac = algebraic_cofactor(r, c).determinant();
					result.element_at(c, r) = detFac;
				}
			return result;
		}

		// 将本矩阵求逆并返回。
		basic_matrix& invert()
		{
			static_assert(is_squared(), "Worked only for squared matrix.");

			return *this = inverse();
		}

		// 返回本矩阵的逆矩阵。
		basic_matrix inverse() const
		{
			static_assert(is_squared(), "Worked only for squared matrix.");

			auto det = determinant();
			if (det == 0)
				return *this;
			auto result = adjoint();
			result /= det;
			return result;
		}

		// 返回一个新的矩阵，这个新矩阵对应的行列式是本矩阵对应的行列式的第[r]行第[c]列元素的代数余子式。
		basic_matrix<Ty, RowDimension - 1, RowDimension - 1>
			algebraic_cofactor(row_dimension_type r, row_dimension_type c) const
		{
			static_assert(is_squared(), "Worked only for squared matrix.");

			static_assert(RowDimension > 0,
				"Only square matrix with non-zero dimension(s) may get the algebraic_cofactor.");

			basic_matrix<Ty, RowDimension - 1, RowDimension - 1> result;
			_cofactor(MyBase::data(), row_dimension(), r, c, result.data());
			return result;
		}

		// 计算本矩阵对应的行列式的值。
		value_type determinant() const
		{
			static_assert(is_squared(), "Worked only for squared matrix.");

			if constexpr (RowDimension == 1)
				return data()[0];
			else if constexpr (RowDimension == 2)
				return _det2(MyBase::data());
			else if constexpr (RowDimension == 3)
				return _det3(MyBase::data());
			else if constexpr (RowDimension == 4)
				return _det4(MyBase::data());
			else
				return _det(MyBase::data(), RowDimension);
		}
	private:
		inline MyBase::size_type _elemIndex(row_dimension_type rowdimension, column_dimension_type columndimension) const
		{
			return columndimension * RowDimension + rowdimension;
		}

		// 存储[fromdim]阶行列式[from]的元素[dstr, dstc]的代数余子式到[result]中。
		static void _cofactor(
			const Ty *from,
			MyBase::size_type fromdim,
			MyBase::size_type dstr,
			MyBase::size_type dstc,
			Ty *result)
		{
			assert(dstr < fromdim && dstc < fromdim);

			typename MyBase::size_type r = 0, c = 0, iOutput = 0;

			value_type fac = static_cast<value_type>((dstr + dstc) % 2 == 0 ? 1 : -1);

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
			Ty *arr, row_dimension_type dim)
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

	template <typename Ty, std::size_t Dimension>
	using basic_square_matrix = basic_matrix<Ty, Dimension, Dimension>;

	template <typename Ty>
	using basic_square_matrix_2d = basic_square_matrix<Ty, 2>;

	template <typename Ty>
	using basic_square_matrix_3d = basic_square_matrix<Ty, 3>;

	template <typename Ty>
	using basic_square_matrix_4d = basic_square_matrix<Ty, 4>;

	template <std::size_t Dimension>
	using square_matrix = basic_square_matrix<double, Dimension>;

	using matrix_2d = basic_square_matrix_2d<double>;

	using matrix_3d = basic_square_matrix_3d<double>;

	using matrix_4d = basic_square_matrix_4d<double>;

	template <std::size_t Dimension>
	using square_matrix_f = basic_square_matrix<float, Dimension>;

	using matrix_2f = basic_square_matrix_2d<float>;

	using matrix_3f = basic_square_matrix_3d<float>;

	using matrix_4f = basic_square_matrix_4d<float>;

	/* TEMPLATE FUNCTION identity_matrix
	Return an identity [Dimension]x[Dimension] matrix with elements of [Ty].
	*/
	template <typename Ty, std::size_t Dimension>
	basic_square_matrix<Ty, Dimension>
		identity_matrix()
	{
		basic_square_matrix<Ty, Dimension> result;
		for (decltype(result.row_dimension()) i = 0; i < result.row_dimension(); ++i)
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
		resultType result(vec);
		return mat *= result;
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
			ResultType;
		ResultType result(vec);
		return result *= mat;
	}
}