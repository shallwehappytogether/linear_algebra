#pragma once

namespace lin
{
	template <typename Ty>
	constexpr Ty pi = static_cast<Ty>(3.14159265358979323846264338327950288l);

	template <typename Ty>
	Ty radius(const Ty &deg_)
	{
		return static_cast<Ty>(deg_ / 180.0 * pi<Ty>);
	}

	template <typename Ty>
	Ty degrees(const Ty &rad_)
	{
		return static_cast<Ty>(rad_ / pi<Ty> * 180.0);
	}
}