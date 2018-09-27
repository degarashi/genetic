#pragma once

namespace gene {
	template <class OS, class Itr, class Sep>
	OS& Explode(OS& os, Itr itr, const Itr itrE, const Sep& sep) {
		bool first = true;
		while(itr != itrE) {
			if(first)
				first = false;
			else
				os << sep;
			os << *itr;
			++itr;
		}
		return os;
	}
}
