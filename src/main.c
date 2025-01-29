#include "gauss-jordan.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdint.h>
#include <stdbool.h>
#include <signal.h>
#include <limits.h>
#include <math.h>
#define eprintf(...) fprintf(stderr, __VA_ARGS__)

struct matrix {
	double **matrix;
	int m, n;
};

double get(int x, int y, void *data) {
	struct matrix *matrix = (struct matrix *) data;
	if (x < 0 || y < 0 || x >= matrix->m || y >= matrix->n) raise(SIGABRT);
	return matrix->matrix[x][y];
}
void set(int x, int y, double value, void *data) {
	struct matrix *matrix = (struct matrix *) data;
	if (x < 0 || y < 0 || x >= matrix->m || y >= matrix->n) raise(SIGABRT);
	matrix->matrix[x][y] = value;
}

char *letters[] = {"x", "y", "z", "w", "v", "u", "t", "s", "r"};

static int frac_places(double x) {
	if (!isfinite(x)) return 0;
	x = fabs(x);
	x = fmod(x, 1); // get fractional part

	int places;
	for (places = 0; x > 1e-12 && places < 10; ++places) {
		x *= 10.0;      // shift decimal point
		x = fmod(x, 1); // get fractional part
	}
	return places;
}

#define FRAC(x) (frac_places(x)), (x)

static void print_d(double d, bool first, bool one) {
	bool minus = d < 0;
	if (minus) d = -d;

	if (!first) printf(" ");

	// print sign
	if (minus)
		printf("-");
	else if (!first)
		printf("+");

	// print space
	if (!first) printf(" ");

	if (isnan(d)) {
		printf("nan");
		return;
	}
	if (isinf(d)) {
		printf("inf");
		return;
	}

	d = fabs(d);
	if (!one && d == 1) return;
	printf("%.*f", FRAC(d));
}

void print_matrix(struct matrix mat) {
	for (int i = 0; i < mat.m; ++i) {
		for (int j = 0; j < mat.n; ++j) {
			double v = get(i, j, &mat);
			if (j == mat.m) printf(" | ");
			print_d(v, j == 0 || j == mat.m, true);
		}
		printf("\n");
	}
	printf("\n");
}

void print_system(struct matrix mat) {
	for (int i = 0; i < mat.m; ++i) {
		bool first = true;
		for (int j = 0, l = 0; j < mat.n; ++j) {
			double v = get(i, j, &mat);
			if (j == mat.m) {
				first = true;
				printf(" = ");
			}
			if (v == 0) goto inc_l;

			// print coefficient
			print_d(v, first, false);
			first = false;

			if (j != mat.m) {
				printf("%s", letters[l]);
			inc_l:
				++l;
			}
		}
		printf("\n");
	}
	printf("\n");
}

#define ROW(...) ((double[]) {__VA_ARGS__})

int main() {
	struct matrix mat = {
	        (double *[]) {
	                      ROW(2, -3, 4, 6),
	                      ROW(3, 4, -5, 7),
	                      ROW(4, -5, 6, 8)},
	        .m = 3, .n = 4
    };

	if (mat.n < mat.m) return 1;

	print_system(mat);
	print_matrix(mat);

	gauss_jordan(get, set, mat.m, mat.n, &mat);

	print_matrix(mat);
	print_system(mat);

	bool identity = true;

	for (int i = 0; i < mat.m; ++i)
		for (int j = 0; j < mat.m; ++j)
			// check if the matrix is the identity matrix
			if (get(i, j, &mat) != (i == j ? 1.0 : 0.0))
				identity = false;

	printf("%solution found\n", identity ? "S" : "No s");

	return 0;
}
