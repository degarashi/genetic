#pragma once
#include <cstddef>
#include <cassert>
#include <random>

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
				const auto randI = [&rd](auto... arg){
					using UD = std::uniform_int_distribution<size_t>;
					return UD(arg...)(rd);
				};
				// 交叉に使う親個体をランダムに抽出
				const auto parent = pool.extractRandom(rd, _nParent);
				using Gene = typename Pool::Gene_t;
				std::vector<Gene> child;

				// 交叉で選んだ親個体で子個体を生成
				const auto nc = _nChild;
				while(child.size() < nc) {
					const size_t i0 = randI(0, parent.length-2),
								i1 = randI(i0+1, parent.length-1);
					auto [c0, c1] = cross(rd, parent.data[i0].gene, parent.data[i1].gene);
					child.emplace_back(std::move(c0));
					child.emplace_back(std::move(c1));
				}
				if(child.size() > nc)
					child.pop_back();

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
