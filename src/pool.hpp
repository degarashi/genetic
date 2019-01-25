#pragma once
#include "ptrlen.hpp"
#include <vector>
#include <iterator>
#include <random>
#include <algorithm>
#include <cassert>

namespace gene {
	// 遺伝子プール
	template <class Gene, class Fit>
	class Pool {
		public:
			using Score = decltype(std::declval<Fit>()(std::declval<Gene>()));
			using Gene_t = Gene;
			struct Ent {
				Gene	gene;
				Score	score;

				Ent(Gene&& g, const Fit& fit):
					gene(std::move(g)),
					score(fit(gene))
				{}
				bool operator < (const Ent& ent) const noexcept {
					return score < ent.score;
				}
			};
		private:
			using EntV = std::vector<Ent>;

			Fit		_fit;
			EntV	_gene;

		public:
			template <
				class RAND,
				class FitA,
				class... Args
			>
			Pool(RAND& rd, FitA&& fit, const size_t population, const Args&... args):
				_fit(std::forward<FitA>(fit))
			{
				for(size_t i=0 ; i<population ; i++) {
					_gene.emplace_back(Gene::MakeRandom(rd, args...), _fit);
				}
			}

			const Ent& getBest() {
				assert(!_gene.empty());
				Score best = std::numeric_limits<Score>::lowest();
				const Ent* e = &_gene.front();
				for(auto& ent : _gene) {
					if(best < ent.score) {
						best = ent.score;
						e = &ent;
					}
				}
				return *e;
			}
			// 集団からN個体をランダムに抽出(末尾に寄せてそのポインタを返す)
			template <class RAND>
			PtrLen<Ent> extractRandom(RAND& rd, const size_t n) {
				const auto ng = nGene();
				if(ng<=1 || n>=ng)
					return {_gene.data(), ng};

				size_t cR = ng;
				for(size_t i=0 ; i<n ; i++)
					extractSingle(rd, --cR);
				return {_gene.data() + cR, ng-cR};
			}
			// Indexよりちいさい場所の個体をひとつ選んでIndexの個体と取り替え
			template <class RAND>
			void extractSingle(RAND& rd, const size_t idx) {
				assert(idx < nGene());
				if(idx == 0)
					return;
				const auto target = std::uniform_int_distribution<size_t>(0, idx)(rd);
				if(target == idx)
					return;
				std::swap(_gene[target], _gene[idx]);
			}
			// 末尾からN個体を削除
			size_t popBack(size_t n) {
				n = std::min(nGene(), n);
				_gene.erase(_gene.begin()+_gene.size()-n, _gene.end());
				return n;
			}
			// 集団の中へ個体を戻す(move)
			template <class Itr>
			void put(Itr itr, const Itr itrE) {
				while(itr != itrE) {
					_gene.emplace_back(std::move(*itr++), _fit);
				}
			}
			// 優秀な遺伝子をN個、集団へ加える(move)
			template <class Itr>
			void putPrime(Itr itr, const Itr itrE, const size_t nPut) {
				assert(std::ptrdiff_t(nPut) <= (itrE-itr));
				EntV tmp;
				while(itr != itrE) {
					tmp.emplace_back(std::move(*itr), _fit);
					++itr;
				}
				std::sort(tmp.begin(), tmp.end());
				std::move(tmp.end() - nPut, tmp.end(), std::back_inserter(_gene));
			}
			// ソート済みの個体群を取得
			PtrLen<Ent> getSorted() {
				std::sort(_gene.begin(), _gene.end());
				return {_gene.data(), nGene()};
			}
			size_t nGene() const noexcept {
				return _gene.size();
			}
	};
}
