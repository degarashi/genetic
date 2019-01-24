#pragma once
#include "../../vector.hpp"
#include <cstddef>
#include <utility>
#include <cassert>
#include <random>

namespace gene::real::cross {
	class Simplex {
		private:
			using value_t = double;
			using vec_t = Vec<value_t>;
			static_assert(std::is_floating_point_v<value_t>);

			size_t		_dim;
			value_t		_eps;
		public:
			Simplex(const size_t dim, const value_t eps):
				_dim(dim),
				_eps(eps)
			{
				assert(dim > 1);
			}
			Simplex(const size_t dim):
				Simplex(dim, std::sqrt(dim+2))
			{}
			size_t prepare() const noexcept {
				return _dim + 1;
			}
			template <class RAND, class Gene>
			std::vector<Gene> crossover(RAND& rd, const std::vector<const Gene*>& src) const {
				const auto size = src.size();
				assert(size == _dim+1);
				assert(src.at(0)->length() == _dim);

				vec_t g(_dim, 0);
				for(auto&& s : src) {
					g += *s;
				}
				g /= size;

				std::vector<vec_t>		s(_dim+1);
				std::vector<value_t>	r(_dim+1);
				std::uniform_real_distribution<value_t> dist(0, 1);
				for(size_t i=0 ; i<size ; i++) {
					s[i] = g + (*src[i] - g) * _eps;
					r[i] = std::pow(dist(rd), 1/value_t(i+1));
				}

				std::vector<Gene> c(size);
				c[0] = vec_t(_dim, 0);
				for(size_t i=1 ; i<size ; i++)
					c[i] = (s[i-1] - s[i] + c[i-1]) * r[i-1];

				std::vector<Gene> ret(size);
				for(size_t i=0 ; i<size ; i++)
					ret[i] = s[i] + c[i];
				return ret;
			}
	};
}
