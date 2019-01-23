#pragma once
#include "explode.hpp"
#include <vector>
#include <iostream>

namespace gene {
	// 長さ可変の数値配列
	template <class T>
	struct VariableGene {
		using value_t = T;
		using Ar = std::vector<value_t>;
		Ar	array;

		VariableGene(const size_t len):
			array(len)
		{}
		size_t length() const noexcept {
			return array.size();
		}
		value_t& operator [](const size_t n) noexcept {
			return array[n];
		}
		const value_t& operator [](const size_t n) const noexcept {
			return array[n];
		}
		bool operator == (const VariableGene& g) const noexcept {
			return this->array == g.array;
		}
	};
	template <class T>
	std::ostream& operator << (std::ostream& os, const VariableGene<T>& g) {
		os << "[";
		Explode(os, g.array.begin(), g.array.end(), ", ");
		return os << "]";
	}
	// 固定長の数値配列
	// 長さ可変のビット配列
	// 固定長のビット配列
}
