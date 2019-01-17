#pragma once

namespace gene::mutate {
	// 突然変異操作を何もしない時に使用するダミークラス
	class Nothing {
		public:
			template <class RAND, class Gene>
			void operator()(RAND&, Gene&) const noexcept {}
	};
}
