#pragma once
#include <random>

namespace gene::mutate {
	// 突然変異操作： 要素を一箇所スワップ
	class Swap {
		public:
			template <class RAND, class Gene>
			void operator()(RAND& rd, Gene& g) const {
				const auto randI = [&rd](auto... arg){
					using UD = std::uniform_int_distribution<size_t>;
					return UD(arg...)(rd);
				};
				const auto len = g.length();
				assert(len > 1);
				const size_t i0 = randI(0, len-2),
							i1 = randI(i0+1, len-1);
				std::swap(g[i0], g[i1]);
			}
	};
}
