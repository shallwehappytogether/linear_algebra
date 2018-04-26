#pragma once

namespace lin
{
	namespace impl
	{
		template <std::size_t N, typename Ty, Ty ...Values>
		struct get_nth_param_value;

		template <typename Ty, Ty HeadValue, Ty ...RemainValues>
		struct get_nth_param_value<0, Ty, HeadValue, RemainValues...>
		{
			constexpr static Ty value = HeadValue;
		};

		template <std::size_t N, typename Ty, Ty HeadValue, Ty ...RemainValues>
		struct get_nth_param_value<N, Ty, HeadValue, RemainValues...>
		{
			constexpr static Ty value = get_nth_param_value<N - 1, Ty, RemainValues...>::value;
		};
	}
}