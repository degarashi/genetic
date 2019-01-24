#pragma once
#include "../gene.hpp"
#include "../ptrlen.hpp"

namespace gene {
	namespace path {
		// パス表現
		template <class T>
		struct VariableGene : ::gene::VariableGene<T> {
			using base_t = ::gene::VariableGene<T>;
			using base_t::base_t;
			using value_t = typename base_t::value_t;

			PtrLen<const value_t> getPath() const {
				return {this->array.data(), this->length()};
			}
			// 同じ番号が重複していないか、順列になっているかを確認
			bool checkValidness() const noexcept {
				const auto len = this->length();
				typename base_t::base_t tmp = *this;
				std::sort(tmp.begin(), tmp.end());
				for(size_t i=0 ; i<len ; i++) {
					if(tmp[i] != value_t(i))
						return false;
				}
				return true;
			}
			// 主に初期値生成に使われる、ランダム遺伝子生成
			template <class RAND>
			static VariableGene MakeRandom(RAND& rd, const size_t len) {
				VariableGene ret(len);
				std::iota(ret.begin(), ret.end(), 0);
				std::shuffle(ret.begin(), ret.end(), rd);
				return ret;
			}
		};
	}
}
