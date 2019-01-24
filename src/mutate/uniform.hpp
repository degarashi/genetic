#pragma once
#include <random>

namespace gene::mutate {
	// 突然変異操作： 一様突然変異
	template <class T>
	class Uniform {
		private:
			using value_t = T;
			value_t		_min,
						_max;
		public:
			Uniform(const value_t& min, const value_t& max):
				_min(min),
				_max(max)
			{
				assert(min <= max);
			}
			template <class RAND, class Gene>
			void operator()(RAND& rd, Gene& g) const {
				const auto randI = [&rd](auto... arg){
					using UD = std::uniform_int_distribution<size_t>;
					return UD(arg...)(rd);
				};
				const auto len = g.length();
				assert(len > 1);
				const size_t idx = randI(0, len-1);
				std::uniform_real_distribution<value_t> dist(_min, _max);
				g[idx] = dist(rd);
			}
	};
}
