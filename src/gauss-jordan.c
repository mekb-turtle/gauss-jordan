#include "gauss-jordan.h"

int gauss_jordan(get_func get, set_func set, m_func sub, m_func mul, m_func div, int (*nonzero)(void *a), int m, int n, void *data) {
	if (n < m || n < 1 || m < 1) return -1;
	for (int i = 0; i < m; ++i) { // loop pivot rows
		// find pivot (make sure it's non-zero)
		void *pivot = get(i, i, data);
		if (nonzero(pivot) == 0)
			for (int k = i + 1; k < m; ++k)
				if (nonzero(get(k, i, data)) != 0) {
					// swap rows
					for (int j = 0; j < n; ++j) {
						void *tmp = get(i, j, data);
						set(i, j, get(k, j, data), data);
						set(k, j, tmp, data);
					}
					pivot = get(i, i, data);
					break;
				}

		if (nonzero(pivot) == 0) continue; // skip zero pivot (singular matrix)

		// normalise
		for (int j = 0; j < n; ++j)
			set(i, j, div(get(i, j, data), pivot), data);

		// eliminate
		for (int k = 0; k < m; ++k) {
			if (k == i) continue; // skip pivot row
			void *factor = get(k, i, data);
			for (int j = 0; j < n; ++j) {
				void *value = get(k, j, data);
				void *product = mul(factor, get(i, j, data));
				value = sub(value, product);
				set(k, j, value, data);
			}
		}
	}
	return 0;
}
