#pragma once
#include "gene.hpp"
#include "pool.hpp"

namespace gene {
	// GAフレームワーク
	template <
		class RAND,
		class Gene,
		class Fit,
		class Cross,
		class Mutate,
		class Generation
	>
	class Environment {
		private:
			using Pool_t = Pool<Gene,Fit>;
			using Score = typename Pool_t::Score;
			using GeneEnt_t = typename Pool_t::Ent;

			RAND&		_rand;
			Pool_t		_pool;
			Cross		_cross;
			Mutate		_mutate;
			Generation	_generation;

		public:
			template <class... Ts>
			Environment(
				RAND& rand,
				const Fit& fit,
				const Cross& cross,
				const Mutate& mutate,
				const Generation& generation,
				const size_t population,
				const Ts&... ts
			):
				_rand(rand),
				_pool(rand, fit, population, ts...),
				_cross(cross),
				_mutate(mutate),
				_generation(generation)
			{}

			// 一世代進める
			void advance() {
				_generation(_rand, _pool, _cross, _mutate);
			}
			// 一番成績の良かった遺伝子を返す
			const GeneEnt_t& getBest() {
				return _pool.getBest();
			}
	};
}
