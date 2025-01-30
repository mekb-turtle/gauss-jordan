#include "gauss-jordan.h"
#include "rational.h"

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <stdint.h>
#include <stdbool.h>
#include <signal.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#define eprintf(...) fprintf(stderr, __VA_ARGS__)

struct matrix {
	rational **matrix;
	int m, n;
};

static struct numlist {
	rational value;
	struct numlist *next;
} *tmplist = NULL, *tmplist_tail = NULL;

static void *get(int x, int y, void *data) {
	struct matrix *matrix = (struct matrix *) data;
	if (x < 0 || y < 0 || x >= matrix->m || y >= matrix->n) raise(SIGABRT);

	// create copy of value to avoid old values in pointers
	struct numlist *num = malloc(sizeof(struct numlist));
	if (!num) raise(SIGABRT);
	num->value = matrix->matrix[x][y];
	if (tmplist_tail) {
		tmplist_tail->next = num;
		tmplist_tail = num;
	} else {
		tmplist = tmplist_tail = num;
	}
	return &num->value;
}

static void free_tmplist() {
	for (; tmplist;) {
		void *t = tmplist->next;
		free(tmplist);
		tmplist = t;
	}
}

static void set(int x, int y, void *value_, void *data) {
	rational value = *(rational *) value_;
	struct matrix *matrix = (struct matrix *) data;
	if (x < 0 || y < 0 || x >= matrix->m || y >= matrix->n) raise(SIGABRT);
	matrix->matrix[x][y] = value;
}

// void* math functions
static void *m_sub(void *a_, void *b_) {
	rational a = *(rational *) a_;
	rational b = *(rational *) b_;
	static rational v;
	v = r_subtract(a, b);
	return &v;
}
static void *m_mul(void *a_, void *b_) {
	rational a = *(rational *) a_;
	rational b = *(rational *) b_;
	static rational v;
	v = r_multiply(a, b);
	return &v;
}
static void *m_div(void *a_, void *b_) {
	rational a = *(rational *) a_;
	rational b = *(rational *) b_;
	static rational v;
	v = r_divide(a, b);
	return &v;
}
static int m_nonzero(void *a_) {
	// returns 0 iff the value is 0
	rational a = *(rational *) a_;
	return a.num;
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

struct display_settings {
	bool frac;
};

static void print_d(rational d, struct display_settings settings, bool first, bool one, bool space, int pad, FILE *out, int *len) {
	d = r_simplify(d);

	bool is_negative = (d.num < 0) != (d.den < 0);

	if (d.num < 0) d.num = -d.num;
	if (d.den < 0) d.den = -d.den;

	char fmt[64] = {0};
	int r;

#define PRINT(...) r = snprintf(fmt, sizeof(fmt), __VA_ARGS__)
	if (d.den == 1) {
		if (!one) {
			fmt[0] = '\0';
			r = 0;
		} else {
			// print as integer
			PRINT("%d", d.num);
		}
	} else if (d.den == 0)
		PRINT("inf");
	else if (d.num == 0)
		PRINT("0");
	else if (settings.frac) {
		// print as fraction
		PRINT("%d/%d", d.num, d.den);
	} else {
		// print with fractional places
		double f = (double) d.num / d.den;
		PRINT("%.*f", FRAC(f));
	}
#undef PRINT

	if (r < 0 || (size_t) r >= sizeof(fmt)) {
		if (out) fprintf(out, "?");
		if (len) *len = 1;
		return;
	}

	// calculate padding
	if (space && !first) ++r;

	if (len) *len = r;
	if (!out) return;

	if (r > pad) pad = 0;
	else
		pad -= r;

	// print padding
	for (int i = 0; i < pad; ++i)
		fprintf(out, " ");

	// print sign
	if (is_negative)
		fprintf(out, "-");
	else if (!first)
		fprintf(out, one ? " " : "+");

	if (space && !first)
		fprintf(out, " ");

	fprintf(out, "%s", fmt);
}

void print_matrix(struct matrix mat, struct display_settings settings) {
	int lens[mat.n];
	for (int j = 0; j < mat.n; ++j) {
		lens[j] = 0;
		for (int i = 0; i < mat.m; ++i) {
			// calculate padding required for each column
			int cur = 0;
			rational *v = get(i, j, &mat);
			print_d(*v, settings, false, true, false, lens[j], NULL, &cur);
			if (cur > lens[j]) lens[j] = cur;
		}
	}
	for (int i = 0; i < mat.m; ++i) {
		bool top = i == 0, bottom = i == mat.m - 1;
#define BORDER(top_, middle_, bottom_) (top ? top_ : bottom ? bottom_ \
	                                                        : middle_)
		printf("%s", BORDER("┌", "│", "└"));
		for (int j = 0; j < mat.n; ++j) {
			rational *v = get(i, j, &mat);
			if (j == mat.m) printf(" %s", BORDER("╷", "│", "╵"));
			if (j != 0) printf(" ");
			print_d(*v, settings, false, true, false, lens[j], stdout, NULL);
		}
		printf(" %s", BORDER("┐", "│", "┘"));
		printf("\n");
	}
	printf("\n");
}

void print_system(struct matrix mat, struct display_settings settings) {
	for (int i = 0; i < mat.m; ++i) {
		bool first = true;
		bool term = false;
		for (int j = 0, l = 0; j < mat.n; ++j) {
			rational *v = get(i, j, &mat);

			if (j == mat.m) {
				first = true;
				if (!term) {
					printf("0 "); // if no terms
				}
				printf("= ");
			}
			if (v->num == 0) goto inc_l;

			// print coefficient
			print_d(*v, settings, first, j == mat.m, true, 0, stdout, NULL);
			first = false;

			if (j != mat.m) {
				if ((size_t) l >= sizeof(letters) / sizeof(letters[0])) {
					// prevent out of bounds
					printf("?");
					continue;
				}
				printf("%s ", letters[l]);
				if (j < mat.m) term = true;
			inc_l:
				++l;
			} else
				printf(" ");
		}
		printf("\n");
	}
	printf("\n");
}

#define ROW(...) ((rational[]) {__VA_ARGS__})
#define F(x, y) ((rational) {(x), (y)})

int main(int argc, char *argv[]) {
	int opt;
	bool invalid = false;
	struct display_settings settings = {.frac = true};

	// argument handling
	while ((opt = getopt_long(argc, argv, ":hVf", (struct option[]) {
	                                                       {"help",        no_argument, 0, 'h'},
	                                                       {"version",     no_argument, 0, 'V'},
	                                                       {"no-fraction", no_argument, 0, 'f'},
	                                                       {0,             0,           0, 0  }
    },
	                          NULL)) != -1) {
		switch (opt) {
			case 'h':
				printf("Usage: %s [OPTION]... <main size> <augment cols> <main matrix>... <augment matrix>...\n", PROJECT_NAME);
				printf("Solves a system of linear equations using Gauss-Jordan elimination\n");
				printf("The main matrix must be m x m, where m = main size\n");
				printf("The augment matrix must be m x n, where n = augment cols\n");
				printf("-h --help: Shows help text\n");
				printf("-V --version: Shows the version\n");
				printf("-f --no-fraction: Print decimal numbers instead of fractions\n");
				return 0;
			case 'V':
				printf("%s %s\n", PROJECT_NAME, PROJECT_VERSION);
#ifdef PROJECT_URL
				printf("See more at %s\n", PROJECT_URL);
#endif
				return 0;
			default:
				if (invalid) break;
				switch (opt) {
					case 'f':
						settings.frac = false;
						break;
					default:
						invalid = true;
						break;
				}
		}
	}

	if (invalid) {
		eprintf("Invalid syntax, use --help for help\n");
		return 1;
	}

	struct matrix mat = {
	        (rational *[]) {
	                        ROW(F(2, 1), F(-3, 1), F(4, 1), F(6, 1), F(2, 1)),
	                        ROW(F(3, 1), F(4, 1), F(-5, 1), F(7, 1), F(1, 1)),
	                        ROW(F(4, 1), F(-5, 1), F(6, 1), F(8, 1), F(-9, 1))},
	        .m = 3, .n = 4
    };

	if (mat.n < mat.m) return 1;

	print_system(mat, settings);
	print_matrix(mat, settings);

	int r = gauss_jordan(get, set, m_sub, m_mul, m_div, m_nonzero, mat.m, mat.n, &mat);
	if (r < 0) {
		eprintf("Invalid matrix\n");
		free_tmplist();
		return 1;
	}

	print_matrix(mat, settings);
	print_system(mat, settings);

	bool identity = true;

	for (int i = 0; i < mat.m; ++i)
		for (int j = 0; j < mat.m; ++j) {
			rational *g = get(i, j, &mat);
			// check if the matrix is the identity matrix
			if (g->num != (i == j ? g->den : 0.0))
				identity = false;
		}

	printf("%solution found\n", identity ? "S" : "No s");

	free_tmplist();

	return 0;
}
