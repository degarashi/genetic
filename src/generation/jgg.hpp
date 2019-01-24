#pragma once
#include "../aux.hpp"
#include <cstddef>

namespace gene {
	// 世代交代モデル：Just Generation Gap
	class JustGenerationGap {
		private:
			size_t		_nParent,
						_nChild;
		public:
			// 集団からnParent分の個体を抽出して、nChildの個体を作り、nParent数だけ戻す
			JustGenerationGap(const size_t nParent, const size_t nChild):
				_nParent(nParent),
				_nChild(nChild)
			{
				assert(nParent >= 2);
				assert(nChild >= nParent);
			}
			template <class RAND, class Pool, class Cross, class Mutate>
			void operator()(RAND& rd, Pool& pool, const Cross& cross, Mutate& mutate) const {
				assert(pool.nGene() > 1);
				// 交叉に使う親個体をランダムに抽出
				const auto parent = pool.extractRandom(rd, _nParent);
				using Gene = typename Pool::Gene_t;
				std::vector<Gene> child;

				// 交叉で選んだ親個体で子個体を生成
				const auto nc = _nChild;
				std::vector<const Gene*> ptr;
				while(child.size() < nc) {
					const auto nParent = cross.prepare();
					const auto idx = PickRandomIndices<size_t>(
										rd,
										nParent,
										parent.length
									);
					ptr.resize(nParent);
					for(size_t i=0 ; i<nParent ; i++)
						ptr[i] = &(parent.data[idx[i]].gene);
					const auto c = cross.crossover(rd, ptr.data());
					for(auto&& ch : c)
						child.emplace_back(std::move(ch));
				}
				// 子個体数の上限補正
				child.resize(std::min(child.size(), nc));

				// 交叉に使った親個体は取り除く
				pool.popBack(parent.length);
				// 突然変異
				for(auto& c : child)
					mutate(rd, c);
				// 集団の中へ戻す
				pool.putPrime(child.begin(), child.end(), _nParent);
			}
	};
}
