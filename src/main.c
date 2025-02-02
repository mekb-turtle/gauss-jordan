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
	rational *matrix;
	int m, n;
};

static struct numlist {
	rational value;
	struct numlist *next;
} *tmplist = NULL, *tmplist_tail = NULL;

rational *matrix_get(struct matrix matrix, int row, int col) {
	if (col < 0 || row < 0 || col >= matrix.n || row >= matrix.m) return NULL;
	return &matrix.matrix[row * matrix.n + col];
}

static void *get(int row, int col, void *data) {
	struct matrix *matrix = (struct matrix *) data;
	rational *v = matrix_get(*matrix, row, col);
	if (!v) raise(SIGABRT);

	// create copy of value to avoid old values in pointers
	struct numlist *num = malloc(sizeof(struct numlist));
	if (!num) raise(SIGABRT);
	num->value = *v;
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

static void set(int row, int col, void *value_, void *data) {
	rational value = *(rational *) value_;
	struct matrix *matrix = (struct matrix *) data;

	rational *v = matrix_get(*matrix, row, col);
	if (!v) raise(SIGABRT);

	*v = value;
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

struct display_settings {
	bool frac, latex;
};

struct digit_add_data {
	char *fmt;
	size_t len;
	int *fmt_len;
	bool first, reoccuring;
	struct display_settings settings;
};

static bool write_digit(int digit, bool is_reoccuring, void *data) {
	int res;
	size_t max_len_tmp;
	struct digit_add_data *d = (struct digit_add_data *) data;

#define PRINT(...) (max_len_tmp = d->len - *d->fmt_len, res = snprintf(d->fmt + *d->fmt_len, max_len_tmp, __VA_ARGS__), *d->fmt_len += res)
#define RES() \
	if (!(res >= 0 && (size_t) res < max_len_tmp)) return false

	if (d->first) {
		// write decimal point
		PRINT(".");
		RES();
		d->first = false;
	}
	if (!d->reoccuring && is_reoccuring) {
		if (d->settings.latex) PRINT("\\overline{"); // start overline
		else
			PRINT("("); // opening bracket
		RES();
		d->reoccuring = true;
	}
	PRINT("%d", digit);
	RES();
	return true;
#undef PRINT
#undef RES
}

static void print_d(rational d, struct display_settings settings, bool first, bool one, bool space, int pad, FILE *out, int *len) {
	if (settings.latex) {
		pad = 0;
		space = false;
	}

	d = r_simplify(d);
	bool is_neg = r_is_negative(d);
	d = r_abs(d);

	char fmt[64] = {0};
	int fmt_len = 0, res = 0;
	size_t max_len_tmp = 0;

#define PRINT(...) (max_len_tmp = sizeof(fmt) - fmt_len, res = snprintf(fmt + fmt_len, max_len_tmp, __VA_ARGS__), fmt_len += res)
#define RES() \
	if (!(res >= 0 && (size_t) res < max_len_tmp)) goto invalid
	if (d.den == 1) {
		if (d.num == 1 && !one) {
			fmt[0] = '\0';
			res = fmt_len = 0;
			max_len_tmp = 1;
		} else {
			// print as integer
			PRINT("%d", d.num);
		}
	} else if (d.den == 0)
		PRINT(settings.latex ? "\\infty" : "inf");
	else if (d.num == 0)
		PRINT("0");
	else if (settings.frac)
		// print as fraction
		PRINT(settings.latex ? "\\tfrac{%d}{%d}" : "%d/%d", d.num, d.den);
	else {
		// print as fractional digits
		int int_part = r_get_integer(d);
		PRINT("%d", int_part);
		RES();

		struct digit_add_data ctx = {fmt, sizeof(fmt), &fmt_len, true, false, settings};

		// printf("DECOMPOSE: %d %d %d\n", d.num, d.den, is_reoccuring(d));
		if (!r_decompose(d, write_digit, &ctx)) goto invalid;
		if (ctx.reoccuring) {
			if (settings.latex)
				PRINT("}"); // close overline
			else
				PRINT(")"); // closing bracket
		}
	}

	RES();
#undef PRINT
#undef RES

	// calculate padding
	if (space && !first) ++fmt_len;

	if (len) *len = fmt_len;
	if (!out) return;

	if (fmt_len > pad) pad = 0;
	else
		pad -= fmt_len;

	// print padding
	for (int i = 0; i < pad; ++i)
		fprintf(out, " ");

	// print sign
	if (is_neg)
		fprintf(out, "-");
	else if (!first) {
		if (!one) fprintf(out, "+");
		else if (!settings.latex)
			fprintf(out, " ");
	}

	if (space && !first)
		fprintf(out, " ");

	fwrite(fmt, 1, fmt_len, out);
	return;

invalid:
	if (out) fprintf(out, "?");
	if (len) *len = 1;
	return;
}

void print_matrix(struct matrix mat, struct display_settings settings) {
	int lens[mat.n];
	for (int col = 0; col < mat.n; ++col) {
		lens[col] = 0;
		if (settings.latex) continue; // no need for padding
		for (int row = 0; row < mat.m; ++row) {
			// calculate padding required for each column
			int cur = 0;
			rational *v = matrix_get(mat, row, col);
			if (!v) raise(SIGABRT);
			print_d(*v, settings, false, true, false, 0, NULL, &cur);
			if (cur > lens[col]) lens[col] = cur;
		}
	}
	if (settings.latex) {
		printf("\\[\\left[\\begin{array}{");
		for (int i = 0; i < mat.m; ++i)
			printf("c");
		// specify where augment line is
		if (mat.n > mat.m)
			printf("|");
		for (int i = mat.m; i < mat.n; ++i)
			printf("c");
		printf("}");
	}
	for (int row = 0; row < mat.m; ++row) {
		bool top = row == 0, bottom = row == mat.m - 1;
#define BORDER(top_, middle_, bottom_) (top ? top_ : bottom ? bottom_ \
	                                                        : middle_)
		if (!settings.latex) printf("%s", BORDER("┌", "│", "└"));
		for (int col = 0; col < mat.n; ++col) {
			if (settings.latex) {
				if (col > 0) printf("&");
			} else {
				if (col == mat.m) printf(" %s", BORDER("╷", "│", "╵"));
				if (col != 0) printf(" ");
			}

			int cur = 0;
			rational *v = matrix_get(mat, row, col);
			if (!v) raise(SIGABRT);
			print_d(*v, settings, false, true, false, lens[col], stdout, &cur);
		}
		if (settings.latex) {
			printf("\\\\");
			continue;
		}
		printf(" %s", BORDER("┐", "│", "┘"));
		printf("\n");
	}
	if (settings.latex) printf("\\end{array}\\right]\\]");
	printf("\n");
}

void print_system(struct matrix mat, struct display_settings settings) {
	if (settings.latex) printf("\\[\\begin{cases}");
	for (int row = 0; row < mat.m; ++row) {
		bool first = true;
		bool term = false;
		for (int col = 0, l = 0; col < mat.n; ++col) {
			rational *v = matrix_get(mat, row, col);
			if (!v) raise(SIGABRT);

			if (col == mat.m) {
				first = true;
				if (!term) {
					printf("0"); // if no terms
					if (!settings.latex) printf(" ");
				}
				printf("=");
				if (!settings.latex) printf(" ");
			}
			if (v->num == 0) goto inc_l;

			// print coefficient
			print_d(*v, settings, first, col == mat.m, true, 0, stdout, NULL);
			first = false;

			if (col != mat.m) {
				if ((size_t) l >= sizeof(letters) / sizeof(letters[0])) {
					// prevent out of bounds
					printf("?");
					continue;
				}
				printf("%s", letters[l]);
				if (!settings.latex) printf(" ");
				if (col < mat.m) term = true;
			inc_l:
				++l;
			} else if (!settings.latex)
				printf(" ");
		}
		if (settings.latex) {
			if (row != mat.m - 1) printf("\\\\");
		} else
			printf("\n");
	}
	if (settings.latex) printf("\\end{cases}\\]");
	printf("\n");
}

void clone_matrix(struct matrix src, struct matrix *dest) {
	dest->m = src.m;
	dest->n = src.n;
	for (int row = 0; row < src.m; ++row)
		for (int col = 0; col < src.n; ++col) {
			rational *v = get(row, col, &src);
			if (!v) raise(SIGABRT);
			set(row, col, v, dest);
		}
}

#define F(x, y) ((rational) {(x), (y)})

int main(int argc, char *argv[]) {
	int opt;
	bool invalid = false;
	struct display_settings settings = {.frac = true, .latex = false};

	// argument handling
	while ((opt = getopt_long(argc, argv, ":hVfl", (struct option[]) {
	                                                       {"help",        no_argument, 0, 'h'},
	                                                       {"version",     no_argument, 0, 'V'},
	                                                       {"no-fraction", no_argument, 0, 'f'},
	                                                       {"latex",       no_argument, 0, 'l'},
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
				printf("-l --latex: Print matrices and systems in LaTeX format\n");
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
					case 'l':
						settings.latex = true;
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
	        (rational[]) {
	                      F(2, 1), F(-3, 1), F(6, -1), F(6, 1), F(2, 1),
	                      F(3, 1), F(4, 1), F(-5, 1), F(7, 1), F(1, 1),
	                      F(4, 1), F(-5, 1), F(6, 1), F(8, 1), F(-9, 1)},
	        .m = 3, .n = 5
    };

	if (mat.n < mat.m) return 1;

	struct matrix old_mat;
	rational old_mat_[mat.m * mat.n];
	old_mat.matrix = old_mat_;
	clone_matrix(mat, &old_mat);

	int r = gauss_jordan(get, set, m_sub, m_mul, m_div, m_nonzero, mat.m, mat.n, &mat);
	if (r < 0) {
		eprintf("Invalid matrix\n");
		free_tmplist();
		return 1;
	}

	print_system(old_mat, settings);
	print_matrix(old_mat, settings);

	print_matrix(mat, settings);
	print_system(mat, settings);

	bool identity = true;

	for (int row = 0; row < mat.m; ++row)
		for (int col = 0; col < mat.m; ++col) {
			rational *g = matrix_get(mat, row, col);
			if (!g) raise(SIGABRT);
			// check if the matrix is the identity matrix
			if (g->num != (row == col ? g->den : 0.0))
				identity = false;
		}

	if (settings.latex) printf("\\[\\text{");
	printf("%solution found", identity ? "S" : "No s");
	if (settings.latex) printf("}\\]");
	printf("\n");

	free_tmplist();

	return 0;
}
