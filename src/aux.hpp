#pragma once
#include <vector>
#include <random>
#include <cassert>

namespace gene {
	namespace detail {
		template <class T>
		using Dist = std::uniform_int_distribution<T>;
	}
	//! 重複のないランダムなN個のインデックス(0からmax-1まで)を生成
	template <class T, class RAND>
	auto PickRandomIndices(RAND& rd, const T n, const T max) {
		assert(max >= n);

		std::vector<T> ret(n);
		if(max > 0) {
			T r = max - n,
			  l = 0;
			auto* dst = ret.data();
			for(T i=0 ; i<n ; i++) {
				const auto idx = detail::Dist<T>{l, r}(rd);
				++r;
				*dst++ = idx;
				l = idx+1;
				assert(l<=r && r<=max);
			}
			assert(dst == ret.data()+n);
		}
		return ret;
	}
}
