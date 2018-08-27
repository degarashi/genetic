#pragma once
#include <cstddef>

namespace gene {
	template <class T>
	struct PtrLen {
		T*		data;
		size_t	length;
	};
}
