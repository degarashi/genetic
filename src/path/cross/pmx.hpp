#pragma once
#include <cassert>
#include <cstdint>
#include <vector>
#include <random>

namespace gene::path::cross {
	// 部分写像交叉 (PMX) 二点交叉
	class PartiallyMapped {
		private:
			template <class Gene>
			Gene _make(const size_t i0, const size_t i1,
						const Gene& inner, const Gene& outer) const
			{
				assert(i0 < i1);
				const auto len = outer.length();
				std::vector<bool> used(len);

				Gene ret(len);
				// i0とi1に挟まれた所はinnerから、それ以外はouterから取る
				for(size_t i=i0 ; i<i1 ; i++) {
					ret[i] = inner[i];
					// 使用済みの数字を記録
					used[ret[i]] = true;
				}

				// 空の遺伝子座リスト
				std::vector<uint8_t> empty;
				const auto putUniqueNum = [&](const auto i){
					const auto c = outer[i];
					if(!used[c]) {
						used[c] = true;
						ret[i] = c;
					} else
						empty.emplace_back(i);
				};
				// 規定の物と衝突を起こさない分についてはそのまま代入
				for(size_t i=0 ; i<i0 ; i++)
					putUniqueNum(i);
				for(size_t i=i1 ; i<len ; i++)
					putUniqueNum(i);

				// usedリストをどこまで使ったか
				size_t cur = 0;
				const auto getNextNum = [&used, &cur](){
					for(;;) {
						assert(cur < used.size());
						if(!used[cur]) {
							used[cur] = true;
							return cur++;
						}
						++cur;
					}
				};
				// 残った数字を順番に入れていく
				const auto empLen = empty.size();
				for(size_t i=0 ; i<empLen ; i++)
					ret[empty[i]] = getNextNum();
				return ret;
			}
		public:
			template <class RAND, class Gene>
			std::pair<Gene, Gene> operator()(RAND& rd, const Gene& src0, const Gene& src1) const {
				// 同じ長さの遺伝子(2以上)が対象
				const auto len = src0.length();
				assert(len == src1.length());
				using UD = std::uniform_int_distribution<size_t>;
				// 交叉する地点を2点、決める
				const size_t i0 = UD(0, len-1)(rd),
							   i1 = UD(i0+1, len)(rd);
				assert(i0 != i1);

				return {
					_make(i0, i1, src0, src1),
					_make(i0, i1, src1, src0)
				};
			}
	};
}
