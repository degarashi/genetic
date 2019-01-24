#include "lubee/src/tests/test.hpp"
#include "../gene.hpp"
#include "../real/cross/simplex.hpp"

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
}
