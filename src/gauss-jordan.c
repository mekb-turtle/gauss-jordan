#include "gauss-jordan.h"

void gauss_jordan(double (*get)(int x, int y, void *data), void (*set)(int x, int y, double value, void *data), int m, int n, void *data) {
	if (n < m || n < 1 || m < 1) return;
	for (int i = 0; i < m; ++i) { // loop pivot rows
		// find pivot (make sure it's non-zero)
		double pivot = get(i, i, data);
		if (pivot == 0)
			for (int k = i + 1; k < m; ++k)
				if (get(k, i, data) != 0) {
					// swap rows
					for (int j = 0; j < n; ++j) {
						double tmp = get(i, j, data);
						set(i, j, get(k, j, data), data);
						set(k, j, tmp, data);
					}
					pivot = get(i, i, data);
					break;
				}

		if (pivot == 0) continue; // skip zero pivot (singular matrix)

		// normalise
		for (int j = 0; j < n; ++j)
			set(i, j, get(i, j, data) / pivot, data);

		// eliminate
		for (int k = 0; k < m; ++k) {
			if (k == i) continue; // skip pivot row
			double factor = get(k, i, data);
			for (int j = 0; j < n; ++j) {
				double v = get(k, j, data);
				v -= factor * get(i, j, data);
				set(k, j, v, data);
			}
		}
	}
}
