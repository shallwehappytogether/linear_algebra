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

		template <typename RightTy>
		basic_vector& operator+=(const basic_vector<RightTy, Dimension> &right)
		{
			for (dimension_type i = 0; i != Dimension; ++i)
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
			for (dimension_type i = 0; i != Dimension; ++i)
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

		basic_vector& normalize()
		{
			return *this /= length();
		}

		basic_vector get_normalized() const
		{
			return basic_vector(*this).normalize();
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

	template <typename Ty, std::size_t RowDimension, std::size_t ColumnDimension>
	class basic_matrix
		:public std::array<Ty, RowDimension * ColumnDimension>
	{
		typedef std::array<Ty, RowDimension * ColumnDimension> MyBase;
	private:
		using MyBase::operator[];

		using MyBase::at;
	public:
		typedef std::size_t dimension_type;

		constexpr dimension_type row_dimension() const
		{
			return RowDimension;
		}

		constexpr dimension_type column_dimension() const
		{
			return ColumnDimension;
		}

		reference operator[](const std::pair<dimension_type, dimension_type> &idimensions)
		{
			return MyBase::operator[](_elemIndexAt(idimensions.first, idimensions.second));
		}

		const_reference operator[](const std::pair<dimension_type, dimension_type> &idimensions) const
		{
			return MyBase::operator[](_elemIndexAt(idimensions.first, idimensions.second));
		}

		reference at(dimension_type rowdimension, dimension_type columndimension)
		{
			return MyBase::at(_elemIndexAt(rowdimension, columndimension));
		}

		const_reference at(dimension_type rowdimension, dimension_type columndimension) const
		{
			return MyBase::at(_elemIndexAt(rowdimension, columndimension));
		}
	private:
		size_type _elemIndexAt(dimension_type rowdimension, dimension_type columndimension) const
		{
			return rowdimension * ColumnDimension + columndimension;
		}
	};

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