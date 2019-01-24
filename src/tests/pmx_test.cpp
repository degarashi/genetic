#include "lubee/src/tests/test.hpp"
#include "fit_ascend.hpp"
#include "../path/gene.hpp"
#include "../environment.hpp"
#include "../path/cross/pmx.hpp"
#include "../generation/jgg.hpp"
#include "../bernoulli.hpp"
#include "../mutate/swap.hpp"

namespace gene::test {
	using VariableGene = lubee::test::Random;
	TEST_F(VariableGene, Random) {
		auto& mt = this->mt().refMt();
		using Gene = path::VariableGene<int>;
		const auto randI = [&mt=this->mt()](auto... arg){
			return mt.getUniform<size_t>({arg...});
		};
		const auto gene = Gene::MakeRandom(mt, randI(1, 64));
		ASSERT_TRUE(gene.checkValidness());
	}
	using Crossover = lubee::test::Random;
	TEST_F(Crossover, PMX) {
		auto& mt = this->mt().refMt();
		using Gene = path::VariableGene<int>;
		path::cross::PartiallyMapped pmx;

		const auto randI = [&mt=this->mt()](auto... arg){
			return mt.getUniform<size_t>({arg...});
		};
		const size_t geneLen = randI(1, 64);
		const auto g0 = Gene::MakeRandom(mt, geneLen),
					g1 = Gene::MakeRandom(mt, geneLen);
		ASSERT_EQ(2, pmx.prepare());
		const Gene* parent[2]{&g0, &g1};
		const auto ng = pmx.crossover(mt, parent);
		ASSERT_EQ(2, ng.size());

		// 交叉させた遺伝子の番号が重複してないか
		ASSERT_TRUE(ng[0].checkValidness());
		ASSERT_TRUE(ng[1].checkValidness());
	}
	TEST_F(Crossover, Equal) {
		auto& mt = this->mt().refMt();
		using Gene = path::VariableGene<int>;
		path::cross::PartiallyMapped pmx;

		const auto g0 = Gene::MakeRandom(mt, 32);
		ASSERT_EQ(2, pmx.prepare());
		const Gene* parent[2]{&g0, &g0};
		const auto ng = pmx.crossover(mt, parent);
		ASSERT_EQ(2, ng.size());

		// 自身と交叉したら同一の遺伝子が生成される
		ASSERT_EQ(g0, ng[0]);
		ASSERT_EQ(g0, ng[1]);
	}
	using Env = lubee::test::Random;
	TEST_F(Env, Regular) {
		auto& mt = this->mt().refMt();
		using Gene = path::VariableGene<int>;
		using PMX = path::cross::PartiallyMapped;
		using Mutate = Bernoulli<mutate::Swap>;
		using Env_t = Environment<std::mt19937, Gene, Fit_Ascend, PMX, Mutate, JustGenerationGap>;
		constexpr size_t GeneLen = 8,
						Population = 128,
						NParent = 32,
						NChild = 64;
		constexpr double MutateP = 0.01;
		Fit_Ascend fit;
		Env_t env(
				mt,
				fit,
				PMX(),
				Mutate(MutateP, mutate::Swap()),
				JustGenerationGap(NParent, NChild),
				Population,
				GeneLen
		);
		// 理想スコア
		double ideal;
		{
			Gene idealG(GeneLen);
			std::iota(idealG.begin(), idealG.end(), 0);
			ideal = fit(idealG);
		}
		constexpr const size_t MaxIteration = 0x10000;
		for(size_t i=0 ; i<MaxIteration ; i++) {
			env.advance();
			// 理想スコアに達したら終了
			if(fit(env.getBest().gene) == ideal) {
				return;
			}
		}
		// MaxIteration回数処理しても理想スコアに達しない場合はアルゴリズムに問題がある可能性が非常に高い
		FAIL();
	}
}
