#include "lubee/src/tests/test.hpp"
#include "fit_ascend.hpp"
#include "../path/gene.hpp"
#include "../pool.hpp"
#include  "lubee/src/compare.hpp"

namespace gene::test {
	using Pool = lubee::test::Random;
	namespace {
		template <class T>
		struct Fit_Zero {
			using range_t = lubee::Range<T>;
			range_t		_range;
			Fit_Zero(const range_t r):
				_range(r)
			{}

			template <class Gene>
			double operator()(const Gene& g) const {
				// 評価の時点で値が適正範囲に収まっている事を確認
				for(auto&& gv : g) {
					if(!_range.hit(gv))
						throw "Assertion failed";
				}
				return 0;
			}
		};
		template <class T>
		struct ClipVal {
			using range_t = lubee::Range<T>;
			range_t		_range;
			ClipVal(const range_t r):
				_range(r)
			{}
			template <class Gene>
			void operator()(Gene& g) const noexcept {
				for(auto&& gv : g)
					gv = std::min(_range.to, std::max(gv, _range.from));
			}
		};
	}
	TEST_F(Pool, Clip) {
		auto& mt = this->mt().refMt();
		using value_t = int;
		using Gene = VariableGene<value_t>;
		using L = std::numeric_limits<value_t>;
		using Clip_t = ClipVal<value_t>;
		using Pool_t = ::gene::Pool<Gene, Fit_Zero<value_t>, Clip_t>;
		const auto randI = [&mt=this->mt()](auto... arg){
			return mt.getUniform<size_t>({arg...});
		};
		const size_t population = randI(1, 32);
		const size_t geneLength = randI(1, 32);
		const lubee::Range<value_t> range(L::lowest()/2, L::max()/2);
		Pool_t pool(mt, Fit_Zero(range), Clip_t(range), population, geneLength, L::lowest(), L::max());

		const auto chk = [&pool, range, geneLength](){
			const auto gene = pool.getSorted();
			for(size_t i=0 ; i<gene.length ; i++) {
				for(size_t j=0 ; j<geneLength ; j++)
					ASSERT_TRUE(lubee::IsInRange(gene.data[i].gene[j], range.from, range.to));
			}
		};
		// Init直後の状態をチェック
		chk();

		const auto makePut = [&randI, &mt, geneLength](){
			const size_t np = randI(1, 32);
			std::vector<Gene> toPut(np);
			for(size_t i=0 ; i<np ; i++) {
				toPut[i] = Gene::MakeRandom(mt, geneLength, L::lowest(), L::max());
			}
			return toPut;
		};
		// put後のチェック
		{
			auto p = makePut();
			pool.put(p.begin(), p.end());
			chk();
		}
		// // putPrime後のチェック
		{
			auto p = makePut();
			pool.putPrime(p.begin(), p.end(), randI(1, p.size()));
			chk();
		}
	}
	TEST_F(Pool, Regular) {
		auto& mt = this->mt().refMt();
		using Gene = path::VariableGene<int>;
		using Pool_t = ::gene::Pool<Gene, Fit_Ascend>;

		const auto randI = [&mt=this->mt()](auto... arg){
			return mt.getUniform<size_t>({arg...});
		};
		size_t population = randI(0, 32);
		const size_t geneLength = randI(1, 32);
		Pool_t pool(mt, Fit_Ascend(), NoClip{}, population, geneLength);
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
