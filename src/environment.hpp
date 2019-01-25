#pragma once
#include "gene.hpp"
#include "pool.hpp"

namespace gene {
	// GAフレームワーク
	template <
		class RAND,
		class Gene,
		class Pool,
		class Cross,
		class Mutate,
		class Generation
	>
	class Environment {
		private:
			RAND&		_rand;
			Pool		_pool;
			Cross		_cross;
			Mutate		_mutate;
			Generation	_generation;

		public:
			template <
				class PoolA,
				class CrossA,
				class MutateA,
				class GenerationA
			>
			Environment(
				RAND& rand,
				PoolA&& pool,
				CrossA&& cross,
				MutateA&& mutate,
				GenerationA&& generation
			):
				_rand(rand),
				_pool(std::forward<PoolA>(pool)),
				_cross(std::forward<CrossA>(cross)),
				_mutate(std::forward<MutateA>(mutate)),
				_generation(std::forward<GenerationA>(generation))
			{}

			// 一世代進める
			void advance() {
				const auto nG = _pool.nGene();
				_generation(_rand, _pool, _cross, _mutate);
				assert(nG == _pool.nGene());
			}
			// 一番成績の良かった遺伝子を返す
			decltype(auto) getBest() {
				return _pool.getBest();
			}
	};
}
