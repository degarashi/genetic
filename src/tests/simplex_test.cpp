#include "lubee/src/tests/test.hpp"
#include "../gene.hpp"
#include "../real/cross/simplex.hpp"
#include "../mutate/uniform.hpp"
#include "../environment.hpp"
#include "../generation/jgg.hpp"
#include "../bernoulli.hpp"

namespace gene::test {
	using Crossover = lubee::test::Random;
	TEST_F(Crossover, Simplex) {
		auto& rd = this->mt();
		auto& mt = rd.refMt();
		const size_t geneLen = rd.getUniform<size_t>({2, 32});
		const size_t nGene = geneLen+1;
		using vec_t = Vec<double>;
		using Gene = VariableGene<double, vec_t>;

		std::vector<Gene> gene(nGene);
		std::generate(gene.begin(), gene.end(), [&](){
			return Gene::MakeRandom(mt, geneLen, -1e2, 1e2);
		});

		{
			constexpr const double Eps = 1.0;
			real::cross::Simplex spx(geneLen, Eps);
			ASSERT_GE(gene.size(), spx.prepare());
			std::vector<const Gene*> parent(nGene);
			for(size_t i=0 ; i<nGene ; i++)
				parent[i] = &gene[i];
			const auto child = spx.crossover(mt, parent.data());
			ASSERT_EQ(parent.size(), child.size());

			// 重心
			vec_t g(geneLen, 0);
			for(auto* p : parent) {
				g += *p;
			}
			g /= parent.size();

			// Epsilonが1以下の場合、子は親同士の最大頂点距離以下
			// (面倒なので雑な判定)
			for(size_t i=0 ; i<nGene ; i++) {
				double md = 0;
				const auto& from = gene[i];
				for(size_t j=0 ; j<nGene ; j++) {
					md = std::max(md, from.distance(gene[j]));
				}
				constexpr const double Threshold = 1e-9;
				for(size_t j=0 ; j<nGene ; j++) {
					ASSERT_LE(from.distance(child[j]), md + Threshold);
				}
			}
		}
	}
	namespace {
		// 1, 2, 4, 8, 16, ... となるような数列
		class Fit_Factorial {
			private:
				mutable std::vector<double>		_cache;
				double _getNum(const size_t n) const {
					auto cur = _cache.size();
					if(n >= cur) {
						_cache.resize(n+1);
						while(cur < n+1) {
							_cache[cur] = _cache[cur-1] * 2;
							++cur;
						}
					}
					return _cache[n];
				}
			public:
				Fit_Factorial():
					_cache{1}
				{}
				template <class Gene>
				double operator()(const Gene& g) const {
					const auto len = g.length();
					double score = 0;
					for(size_t i=0 ; i<len ; i++) {
						const auto s = std::abs(g[i] - _getNum(i));
						score -= s;
					}
					return score;
				}
				static double Ideal() noexcept {
					return 0;
				}
		};
	}
	using Env = lubee::test::Random;
	TEST_F(Env, Simplex) {
		auto& mt = this->mt().refMt();
		using Gene = VariableGene<double, Vec<double>>;
		using Simplex = real::cross::Simplex;
		using Mutate = Bernoulli<mutate::Uniform<double>>;
		using Env_t = Environment<std::mt19937, Gene, Fit_Factorial, Simplex, Mutate, JustGenerationGap>;
		constexpr const size_t GeneLen = 4,
						NParent = GeneLen*2,
						NChild = NParent*16,
						Population = NChild*2;
		constexpr double MutateP = 0.01;
		Fit_Factorial fit;
		Env_t env(
			mt,
			fit,
			Simplex(GeneLen),
			Mutate(MutateP, mutate::Uniform<double>(0, 1e2)),
			JustGenerationGap(NParent, NChild),
			Population,
			GeneLen,
			0, 1e2
		);
		constexpr const size_t MaxIteration = 0x10000;
		constexpr const double Threshold = 1e-1;
		for(size_t i=0 ; i<MaxIteration ; i++) {
			env.advance();
			const auto score = fit(env.getBest().gene);
			// 理想スコアに達したら終了
			if(std::abs(score - fit.Ideal()) < Threshold) {
				return;
			}
		}
		// MaxIteration回数処理しても理想スコアに達しない場合はアルゴリズムに問題がある可能性が非常に高い
		FAIL();
	}
}
