#pragma once
#include <cstdint>

namespace gene::test {
	class Fit_Ascend {
		public:
			template <class Gene>
			double operator()(const Gene& g) const {
				const auto len = g.length();
				uint64_t count = 1;
				for(std::size_t i=1 ; i<len ; i++) {
					// 前の数値より低ければその差分だけ引く
					if(g[i-1] > g[i])
						count += (g[i-1] - g[i]) * 4;
					// 大き過ぎも駄目
					else if(g[i] != g[i-1]+1)
						count += g[i] - g[i-1];
				}
				return 1.0 / count;
			}
	};
}
