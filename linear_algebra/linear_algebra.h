#pragma once

#include <cstddef>
#include <array>
#include <cmath>
#include <numeric>
#include <cassert>

namespace linear_algebra
{
	template <typename Ty, std::size_t Dimension>
	class basic_vector
		:public std::array<Ty, Dimension>
	{
	public:
		typedef std::size_t dimension_type;

		template <typename... Args>
		basic_vector(Args... args)
			:std::array<Ty, Dimension>{ static_cast<Ty>(args)... }
		{

		}

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

		basic_vector operator+() const
		{
			return basic_vector(*this);
		}

		basic_vector operator-() const
		{
			basic_vector result(*this);
			for (auto &comp : result)
				comp = -comp;
			return result;
		}

		basic_vector& normalize_to_assign()
		{
			return *this /= length();
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

		template <typename RightTy>
		basic_vector& operator+=(const basic_vector<RightTy, Dimension> &right)
		{
			for (dimension_type i = 0; i < dimension(); ++i)
				left[i] += right[i];
			return *this;
		}

		template <typename RightTy>
		basic_vector<std::common_type_t<Ty, RightTy>, Dimension>
			operator+(const basic_vector<RightTy, Dimension> &right) const
		{
			basic_vector(*this) += right;
		}

		template <typename RightTy>
		basic_vector& operator-=(const basic_vector<RightTy, Dimension> &right)
		{
			for (dimension_type i = 0; i < dimension(); ++i)
				left[i] -= right[i];
			return *this;
		}

		template <typename RightTy>
		basic_vector<std::common_type_t<Ty, RightTy>, Dimension>
			operator-(const basic_vector<RightTy, Dimension> &right) const
		{
			basic_vector(*this) -= right;
		}

		basic_vector& operator*=(const value_type &right)
		{
			for (auto &comp : *this)
				comp *= right;
			return *this;
		}

		basic_vector operator*(const value_type &right) const
		{
			return basic_vector(*this) *= right;
		}

		basic_vector& operator/=(const value_type &right)
		{
			for (auto &comp : *this)
				comp /= right;
			return *this;
		}

		basic_vector operator/(const value_type &right) const
		{
			return basic_vector(*this) /= right;
		}

		template <typename RightTy>
		bool operator==(const basic_vector<RightTy, Dimension> &right) const
		{
			for (dimension_type i = 0; i < dimension(); ++i)
				if ((*this)[i] != right[i])
					return false;
			return true;
		}

		template <typename RightTy>
		bool operator!=(const basic_vector<RightTy, Dimension> &right) const
		{
			return !(*this == right);
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
		:public std::array<Ty, RowDimension * ColumnDimension>
	{
		typedef std::array<Ty, RowDimension * ColumnDimension> MyBase;
	private:
		using MyBase::operator[];

		using MyBase::at;

		using MyBase::begin;

		using MyBase::end;

		constexpr static bool _isSquare = RowDimension == ColumnDimension;
	public:
		basic_matrix()
			:MyBase{}
		{

		}

		typedef std::size_t row_dimension_type;

		typedef std::size_t column_dimension_type;

		typedef std::pair<row_dimension_type, column_dimension_type> location_type;

		constexpr row_dimension_type row_dimension() const
		{
			return RowDimension;
		}

		constexpr column_dimension_type column_dimension() const
		{
			return ColumnDimension;
		}

		reference operator[](const location_type &idimensions)
		{
			return MyBase::operator[](_elemIndexAt(idimensions.first, idimensions.second));
		}

		const_reference operator[](const location_type &idimensions) const
		{
			return MyBase::operator[](_elemIndexAt(idimensions.first, idimensions.second));
		}

		reference at(row_dimension_type rowdimension, column_dimension_type columndimension)
		{
			return MyBase::at(_elemIndexAt(rowdimension, columndimension));
		}

		const_reference at(row_dimension_type rowdimension, column_dimension_type columndimension) const
		{
			return MyBase::at(_elemIndexAt(rowdimension, columndimension));
		}

		template <class Enabled = void>
		basic_matrix& transpose_to_assign()
		{
			static_assert(false,
				"Only square matrix may perform the transpose_to_assign function.");
		}

		template <>
		basic_matrix& transpose_to_assign<std::enable_if_t<_isSquare>>()
		{
			_fillTransposeTo(*this);
			return *this;
		}

		basic_matrix<Ty, ColumnDimension, RowDimension> transpose() const
		{
			basic_matrix<Ty, ColumnDimension, RowDimension> result;
			_fillTransposeTo(result);
			return result;
		}

		template <class Enabled = void>
		basic_matrix& inverse_to_assign()
		{
			static_assert(false,
				"Only square matrix may perform the inverse or inverse_to_assign function.");
		}

		template <>
		basic_matrix& inverse_to_assign<std::enable_if_t<_isSquare>>()
		{
			return *this;
		}

		template <class Enabled = void>
		basic_matrix inverse() const
		{
			return basic_matrix(*this).inverse_to_assign();
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
						v += (*this)[{ r, i }] * right[{ i, c }];
					result[{ r, c }] = v;
				}
			return result;
		}
	private:
		size_type _elemIndexAt(row_dimension_type rowdimension, column_dimension_type columndimension) const
		{
			return rowdimension * ColumnDimension + columndimension;
		}

		void _fillTransposeTo(basic_matrix<Ty, ColumnDimension, RowDimension> &dest) const
		{
			for (row_dimension_type r = 0; r < row_dimension(); ++r)
				for (column_dimension_type c = 0; c < column_dimension(); ++c)
					dest[{ r, c }] = (*this)[{ c, r }];
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

	using square_matrix_2d = basic_square_matrix_2d<double>;

	using square_matrix_3d = basic_square_matrix_3d<double>;

	using square_matrix_4d = basic_square_matrix_4d<double>;

	template <typename Ty, std::size_t Dimension>
	basic_square_matrix<Ty, Dimension>
		make_identity_matrix()
	{
		basic_square_matrix<Ty, Dimension> result;
		for (decltype(result.row_dimension()) i = 0; i < result.row_dimension(); ++i)
			result[{ i, i }] = static_cast<Ty>(1);
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
				v += (vec[c] * mat[{ r, c }]);
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
				v += (vec[r] * mat[{ r, c }]);
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