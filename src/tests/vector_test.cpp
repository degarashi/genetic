#include "lubee/src/tests/test.hpp"
#include "../vector.hpp"

namespace gene {
	template <class RDF>
	auto MakeRandomArray(RDF&& rdf, const size_t n) {
		std::vector<double>	ar(n);
		std::generate(ar.begin(), ar.end(), [&](){
			return rdf();
		});
		return ar;
	}
	template <class T0, class T1>
	bool CompareArray(const T0& t0, const T1& t1) {
		const auto size = t0.size();
		if(size != t1.size())
			return false;
		for(size_t i=0 ; i<size ; i++) {
			if(t0[i] != t1[i])
				return false;
		}
		return true;
	}
	template <class Ar0, class Ar1, class Op>
	auto ProcArray(const Ar0& ar0, const Ar1& ar1, Op&& op) {
		const auto s = ar0.size();
		std::vector<double> ret(s);
		for(size_t i=0 ; i<s ; i++)
			ret[i] = op(ar0[i], ar1[i]);
		return ret;
	}

	struct Vector :
		lubee::test::Random
	{
		using vec_t = Vec<double>;
	};
	TEST_F(Vector, Constructor) {
		using vec_t = Vector::vec_t;
		auto& mt = this->mt();
		const auto n = mt.getUniform<size_t>({1, 64});

		{
			// 1要素での初期化
			const auto num = mt.getUniform<double>();
			const vec_t v(n, num);
			for(size_t i=0 ; i<n ; i++)
				ASSERT_EQ(v[i], num);
		}
		{
			// 配列を指定しての初期化
			const auto ar = MakeRandomArray(mt.getUniformF<double>(), n);
			const vec_t v(ar);
			ASSERT_TRUE(CompareArray(ar, v));
		}
	}
	TEST_F(Vector, Dot) {
		using vec_t = Vector::vec_t;
		auto& mt = this->mt();
		const auto ang = mt.getUniform<double>({-10, 10});
		ASSERT_NEAR(
			vec_t(std::vector<double>{std::sin(ang), std::cos(ang)}).length(),
			1.0,
			1e-3
		);
	}
	TEST_F(Vector, Normalize) {
		using vec_t = Vector::vec_t;
		auto& mt = this->mt();
		const auto n = mt.getUniform<size_t>({1, 64});
		auto ar = MakeRandomArray(mt.getUniformF<double>({-1e2, 1e2}), n);
		if(std::abs(ar[0]) < 1e-2)
			ar[0] = 1;
		vec_t v(ar);
		v.normalize();
		ASSERT_NEAR(v.length(), 1.0, 1e-3);
	}
	TEST_F(Vector, Compare) {
		using vec_t = Vector::vec_t;
		auto& mt = this->mt();
		const auto n = mt.getUniform<size_t>({1, 64});
		const auto ar = MakeRandomArray(mt.getUniformF<double>({-1e2, 1e2}), n);
		vec_t v0(ar),
			  v1(ar);
		ASSERT_EQ(v0, v0);
		ASSERT_EQ(v0, v1);

		const auto idx = mt.getUniform<size_t>({0, n-1});
		v1[idx] += mt.getUniform<double>({1e-2, 1e2});
		ASSERT_NE(v0, v1);
		for(size_t i=0 ; i<n ; i++)
			ASSERT_EQ(i!=idx, v0[i]==v1[i]);
	}
	TEST_F(Vector, BinaryOperator) {
		using vec_t = Vector::vec_t;
		auto& mt = this->mt();
		const auto n = mt.getUniform<size_t>({1, 64});
		const auto ar0 = MakeRandomArray(mt.getUniformF<double>({-1e2, 1e2}), n),
					ar1 = MakeRandomArray(mt.getUniformF<double>({-1e2, 1e2}), n);
		const vec_t v0(ar0),
					v1(ar1);

		#define OpTest(op, obj) \
			{ \
				const auto ar2 = ProcArray(ar0, ar1, obj{}); \
				ASSERT_TRUE(CompareArray(v0 op v1, ar2)); \
				auto v2 = v0; \
				v2 op##= v1; \
				ASSERT_TRUE(CompareArray(v2, ar2)); \
			}
		OpTest(+, std::plus)
		OpTest(-, std::minus)

		const auto s = mt.getUniform<double>({-1e2, 1e2});
		#define OpTest2(op, obj) \
		{ \
			const auto ar2 = ProcArray(ar0, Scalar{s}, obj{}); \
			ASSERT_TRUE(CompareArray(v0 op s, ar2)); \
			auto v2 = v0; \
			v2 op##= s; \
			ASSERT_TRUE(CompareArray(v2, ar2)); \
		}
		OpTest2(*, std::multiplies)
		OpTest2(/, std::divides)
	}
}
