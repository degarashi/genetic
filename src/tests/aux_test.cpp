#include "lubee/src/tests/test.hpp"
#include "../aux.hpp"

namespace gene {
	using AuxFuncs = lubee::test::Random;
	TEST_F(AuxFuncs, PickRandomIndices) {
		auto& mt = this->mt();
		using typ = size_t;
		const auto len = mt.getUniform<typ>({0, 128}),
					pick = mt.getUniform<typ>({0, len});
		auto res = PickRandomIndices<typ>(mt.refMt(), pick, len);
		// pickで指定した長さの配列
		ASSERT_EQ(pick, res.size());
		// インデックスの値は範囲内
		for(auto&& idx : res) {
			ASSERT_GE(idx, 0);
			ASSERT_LE(idx, len-1);
		}
		// インデックスの値に重複がない
		const auto itr = std::unique(res.begin(), res.end());
		ASSERT_EQ(res.end(), itr);
	}
}
