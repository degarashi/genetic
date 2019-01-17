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
		using Gene = order::path::VariableGene<int>;
		const auto randI = [&mt=this->mt()](auto... arg){
			return mt.getUniform<size_t>({arg...});
		};
		const auto gene = Gene::MakeRandom(mt, randI(1, 64));
		ASSERT_TRUE(gene.checkValidness());
	}
	using Crossover = lubee::test::Random;
	TEST_F(Crossover, Regular) {
		auto& mt = this->mt().refMt();
		using Gene = order::path::VariableGene<int>;
		order::cross::PartiallyMapped pmx;

		const auto randI = [&mt=this->mt()](auto... arg){
			return mt.getUniform<size_t>({arg...});
		};
		const size_t geneLen = randI(1, 64);
		const auto g0 = Gene::MakeRandom(mt, geneLen),
					g1 = Gene::MakeRandom(mt, geneLen);
		const auto [ng0, ng1] = pmx(mt, g0, g1);

		// 交叉させた遺伝子の番号が重複してないか
		ASSERT_TRUE(ng0.checkValidness());
		ASSERT_TRUE(ng1.checkValidness());
	}
	TEST_F(Crossover, Equal) {
		auto& mt = this->mt().refMt();
		using Gene = order::path::VariableGene<int>;
		order::cross::PartiallyMapped pmx;

		const auto g0 = Gene::MakeRandom(mt, 32);
		const auto [ng0, ng1] = pmx(mt, g0, g0);

		// 自身と交叉したら同一の遺伝子が生成される
		ASSERT_EQ(g0, ng0);
		ASSERT_EQ(g0, ng1);
	}
	using Env = lubee::test::Random;
	TEST_F(Env, Regular) {
		auto& mt = this->mt().refMt();
		using Gene = order::path::VariableGene<int>;
		using PMX = order::cross::PartiallyMapped;
		using Mutate = Bernoulli<mutate::Swap>;
		using Env_t = Environment<std::mt19937, Gene, Fit_Ascend, PMX, Mutate, JustGenerationGap>;
		constexpr size_t GeneLen = 8,
						Population = 128,
						NParent = 32,
						NChild = 64;
		constexpr double MutateP = 0.01;
		Env_t env(
				mt,
				Fit_Ascend(),
				PMX(),
				Mutate(MutateP, mutate::Swap()),
				JustGenerationGap(NParent, NChild),
				Population,
				GeneLen
		);
		for(size_t i=0 ; i<16; i++) {
			env.advance();
		}
	}
}
