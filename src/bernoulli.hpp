#pragma once
#include <random>

namespace gene::order {
	// ベルヌーイ施行 (一定確率でMutateを実行)
	template <class Mutate>
	class Bernoulli {
		private:
			std::bernoulli_distribution	_dist;
			Mutate						_mutate;

		public:
			Bernoulli(const double probability, const Mutate& mutate):
				_dist(probability),
				_mutate(mutate)
			{}
			template <class RAND, class Gene>
			void operator()(RAND& rd, Gene& g) {
				if(_dist(rd))
					_mutate(rd, g);
			}
	};
}
