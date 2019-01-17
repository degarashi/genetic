#include "lubee/src/tests/test.hpp"
#include "fit_ascend.hpp"
#include "../path/gene.hpp"
#include "../pool.hpp"

namespace gene::test {
	using Pool = lubee::test::Random;
	TEST_F(Pool, Regular) {
		auto& mt = this->mt().refMt();
		using Gene = path::VariableGene<int>;
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
}
