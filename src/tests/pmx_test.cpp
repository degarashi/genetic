#include "lubee/src/tests/test.hpp"
#include "../gene_order.hpp"
#include "../environment.hpp"
#include "../pmx.hpp"
#include "../jgg.hpp"
#include "../bernoulli.hpp"
#include "../mutate/swap.hpp"

class Fit_Ascend {
	public:
		template <class Gene>
		double operator()(const Gene& g) const {
			const auto len = g.length();
			uint64_t count = 1;
			for(size_t i=1 ; i<len ; i++) {
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
	using Pool = lubee::test::Random;
	TEST_F(Pool, Regular) {
		auto& mt = this->mt().refMt();
		using Gene = order::path::VariableGene<int>;
		using Pool_t = ::gene::Pool<Gene, Fit_Ascend>;

		const auto randI = [&mt=this->mt()](auto... arg){
			return mt.getUniform<size_t>({arg...});
		};
		size_t population = randI(0, 32);
		const size_t geneLength = randI(1, 32);
		Pool_t pool(mt, Fit_Ascend(), population, geneLength);
		ASSERT_EQ(population, pool.nGene());

		size_t nAct = randI(1, 32);
		while(nAct != 0) {
			--nAct;
			enum Action {
				Act_Insert,
				Act_Pop,
				Act_Extract,
				Act_Sort,
				N_Act
			};
			switch(randI(0, N_Act)) {
				case Act_Extract:
				{
					// 範囲外の値を入れてもOKなので、そのテスト
					const size_t n = randI(0, population * 2);
					pool.extractRandom(mt, n);
					// 末尾に寄せるだけなので数は変わらない
					ASSERT_EQ(population, pool.nGene());
				}
				break;
				case Act_Pop:
				{
					// 範囲外の値を入れてもOKなので、そのテスト
					size_t n = randI(0, population * 2);
					n = pool.popBack(n);
					population -= n;
					// ここで初めて数が変わる
					ASSERT_EQ(population, pool.nGene());
				}
				break;
				case Act_Insert:
				{
					// ランダムに個体を作って挿入
					const size_t n = randI(0, 32);
					std::vector<Gene> gene;
					for(size_t i=0 ; i<n ;  i++)
						gene.emplace_back(Gene::MakeRandom(mt, geneLength));
					pool.put(gene.begin(), gene.end());
					population += n;
					ASSERT_EQ(population, pool.nGene());
				}
				break;
				case Act_Sort:
				{
					// ソート機能のチェック
					auto cur = std::numeric_limits<Pool_t::Score>::lowest();
					const auto pl = pool.getSorted();
					for(size_t i=0 ; i<pl.length ; i++) {
						const auto score = pl.data[i].score;
						ASSERT_LE(cur, score);
						cur = score;
					}
				}
				break;
			}
		}
	}
	using Env = lubee::test::Random;
	TEST_F(Env, Regular) {
		auto& mt = this->mt().refMt();
		using Gene = order::path::VariableGene<int>;
		using PMX = order::cross::PartiallyMapped;
		using Mutate = order::Bernoulli<mutate::Swap>;
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
