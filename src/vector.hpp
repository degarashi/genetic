#pragma once
#include <vector>
#include <cmath>

namespace gene {
	template <class V>
	struct Scalar {
		V	value;

		Scalar() = default;
		Scalar(const V v):
			value(v) {}
		decltype(auto) operator [](const size_t) noexcept {
			return value;
		}
		decltype(auto) operator [](const size_t) const noexcept {
			return value;
		}
	};
	template <class T, class Container=std::vector<T>>
	class Vec {
		private:
			using value_t = T;
			using container_t = Container;

			container_t		_m;

			template <class Op, class V>
			Vec _op(Op&& op, const V& v) const noexcept {
				const auto s = size();
				Vec ret(s);
				for(size_t i=0 ; i<s ; i++)
					ret[i] = op(_m[i], v[i]);
				return ret;
			}
			template <class Op, class V>
			Vec& _op2(Op&& op, const V& v) noexcept {
				const auto s = size();
				for(size_t i=0 ; i<s ; i++)
					_m[i] = op(_m[i], v[i]);
				return *this;
			}

		public:
			auto begin() const noexcept { return _m.begin(); }
			auto begin() noexcept { return _m.begin(); }
			auto cbegin() const noexcept { return _m.cbegin(); }
			auto end() const noexcept { return _m.end(); }
			auto end() noexcept { return _m.end(); }
			auto cend() const noexcept { return _m.cend(); }

			size_t size() const noexcept {
				return _m.size();
			}
			Vec() = default;
			explicit Vec(const size_t dim):
				_m(dim)
			{}
			Vec(const size_t dim, const value_t v):
				Vec(dim)
			{
				std::fill(_m.begin(), _m.end(), v);
			}
			explicit Vec(const container_t& data):
				_m(data)
			{}
			void resize(const size_t n) {
				_m.resize(n);
			}
			Vec operator * (const value_t r) const noexcept {
				return _op(std::multiplies{}, Scalar{r});
			}
			Vec operator / (const value_t r) const noexcept {
				return _op(std::divides{}, Scalar{r});
			}
			Vec operator + (const Vec& v) const noexcept {
				return _op(std::plus{}, v);
			}
			Vec operator - (const Vec& v) const noexcept {
				return _op(std::minus{}, v);
			}
			Vec& operator += (const Vec& v) noexcept {
				return _op2(std::plus{}, v);
			}
			Vec& operator -= (const Vec& v) noexcept {
				return _op2(std::minus{}, v);
			}
			Vec& operator *= (const value_t r) noexcept {
				return _op2(std::multiplies{}, Scalar{r});
			}
			Vec& operator /= (const value_t r) noexcept {
				return _op2(std::divides{}, Scalar{r});
			}
			decltype(auto) operator [](const std::size_t i) noexcept { return _m[i]; }
			decltype(auto) operator [](const std::size_t i) const noexcept { return _m[i]; }

			value_t dot(const Vec& v) const noexcept {
				value_t sum = 0;
				const auto s = size();
				for(size_t i=0 ; i<s ; i++)
					sum += _m[i] * v[i];
				return sum;
			}
			value_t len_sq() const noexcept {
				return dot(*this);
			}
			value_t length() const noexcept {
				return std::sqrt(len_sq());
			}
			value_t distance(const Vec& v) const noexcept {
				return (v - *this).length();
			}
			bool operator == (const Vec& v) const noexcept {
				const auto s = size();
				for(size_t i=0 ; i<s ; i++)
					if(_m[i] != v[i])
						return false;
				return true;
			}
			bool operator != (const Vec& v) const noexcept {
				return !this->operator == (v);
			}
			void normalize() noexcept {
				const value_t r = 1 / length();
				*this *= r;
			}
	};
}
