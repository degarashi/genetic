#pragma once
#include "explode.hpp"
#include <vector>
#include <iostream>
#include <random>

namespace gene {
	// 長さ可変の数値配列
	template <class T, class Container=std::vector<T>>
	struct VariableGene :
		Container
	{
		using base_t = Container;
		using value_t = T;

		VariableGene() = default;
		VariableGene(const size_t len):
			base_t(len)
		{}
		size_t length() const noexcept {
			return base_t::size();
		}
		VariableGene& operator = (const base_t& v) {
			data() = v;
			return *this;
		}
		bool operator == (const VariableGene& g) const noexcept {
			return data() == g.data();
		}
		base_t& data() noexcept {
			return static_cast<base_t&>(*this);
		}
		const base_t& data() const noexcept {
			return static_cast<const base_t&>(*this);
		}
		template <class RAND>
		static VariableGene MakeRandom(RAND& rd, const size_t len, const value_t min, const value_t max) {
			VariableGene ret(len);
			const auto gen = [&rd, &ret, len](auto&& dist){
				for(size_t i=0 ; i<len ; i++)
					ret[i] = dist(rd);
			};
			if constexpr (std::is_floating_point_v<value_t>) {
				gen(std::uniform_real_distribution{min, max});
			} else {
				gen(std::uniform_int_distribution{min, max});
			}
			return ret;
		}
	};
	template <class T, class C>
	std::ostream& operator << (std::ostream& os, const VariableGene<T,C>& g) {
		os << "[";
		Explode(os, g.begin(), g.end(), ", ");
		return os << "]";
	}
	// 固定長の数値配列
	// 長さ可変のビット配列
	// 固定長のビット配列
}
