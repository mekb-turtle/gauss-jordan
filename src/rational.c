#include "rational.h"
#include <limits.h>

inline int gcd(int a, int b) {
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	while (b != 0) {
		int temp = b;
		b = a % b;
		a = temp;
	}
	return a;
}

inline int lcm(int a, int b) {
	if (a == 0 || b == 0) return 0;
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	return (a * b) / gcd(a, b);
}

bool r_decompose(rational frac, bool (*write)(int digit, bool is_reoccuring, void *data), void *data) {
	if (frac.den == 0) return false;

	frac = r_abs(frac);

	frac.num %= frac.den; // extract fractional part

	frac = r_simplify(frac);

	struct {
		int quot, rem;
	} cur = {frac.num / frac.den, frac.num % frac.den}, num[512];

	int reoccur_i = -1, end_i = -1, i;
	for (i = 0; cur.rem != 0 && (size_t) i < sizeof(num) / sizeof(num[0]); ++i) {
		// perform long division
		cur.rem *= 10;
		cur.quot = cur.rem / frac.den;
		cur.rem %= frac.den;
		for (int j = 0; j < i; ++j) {
			if (num[j].quot == cur.quot && num[j].rem == cur.rem) {
				reoccur_i = j;
				end_i = i;
				goto break_loop;
			}
		}
		num[i] = cur;
	}
	if (reoccur_i == -1) reoccur_i = i; // no reoccurring part or too long

break_loop:

	// write non-reoccurring part
	for (int k = 0; k < reoccur_i; ++k) {
		if (!write(num[k].quot, false, data)) return false;
	}

	// write reoccurring part
	for (int k = reoccur_i; k < end_i; ++k) {
		if (!write(num[k].quot, true, data)) return false;
	}
	return true;
}

bool r_is_negative(rational frac) {
	return (frac.num < 0) != (frac.den < 0);
}

bool r_is_integer(rational frac) {
	return frac.den != 0 && frac.num % frac.den == 0;
}

int r_get_integer(rational frac) {
	return frac.den != 0 ? frac.num / frac.den : 0;
}

int r_cmp(rational a, rational b) {
	a = r_simplify(a);
	b = r_simplify(b);
	a = r_set_denominator(a, b.den);
	b = r_set_denominator(b, a.den);
	return a.num - b.num;
}

bool is_reoccuring(rational frac) {
	// reduce powers of 2 and 5
	while (frac.den % 2 == 0) frac.den /= 2;
	while (frac.den % 5 == 0) frac.den /= 5;
	// factors of 10 are 2 and 5, so if the denominator is not 1, it's reoccurring
	return frac.den != 1;
}

int r_get_decimal(rational frac, bool (*write)(char digit, void *data), void *data) {
	size_t i = 0;
	int quot = frac.num / frac.den, rem = frac.num % frac.den;

	for (; rem != 0; ++i) {
		// perform long division
		rem *= 10;
		quot = rem / frac.den;
		rem %= frac.den;
		if (quot < 0 || quot > 9) return -1;
		char digit = quot + '0';
		if (!write(digit, data)) return -1;
	}
	return i > INT_MAX ? INT_MAX : i;
}

inline rational r_abs(rational a) {
	if (a.num < 0) a.num = -a.num;
	if (a.den < 0) a.den = -a.den;
	return a;
}

inline rational r_set_denominator(rational a, int den) {
	int l = lcm(a.den, den);
	if (a.den == 0 || l == 0) return a;
	a.num *= l / a.den;
	a.den = l;
	return r_simplify(a);
}

inline rational r_add(rational a, rational b) {
	a = r_set_denominator(a, b.den);
	if (a.den == 0) return a;
	b.num *= a.den / b.den;
	a.num += b.num;
	return r_simplify(a);
}

inline rational r_subtract(rational minuend, rational subtrahend) {
	return r_add(minuend, r_negate(subtrahend));
}

inline rational r_multiply(rational a, rational b) {
	a.num *= b.num;
	a.den *= b.den;
	return r_simplify(a);
}

inline rational r_divide(rational dividend, rational divisor) {
	return r_multiply(dividend, r_inverse(divisor));
}

inline rational r_negate(rational a) {
	if (a.den < 0)
		a.den = -a.den;
	else
		a.num = -a.num;
	return a;
}

inline rational r_inverse(rational a) {
	int t = a.num;
	a.num = a.den;
	a.den = t;
	return a;
}

static inline int pow_(int base, int exp) {
	int power = base;
	for (; exp > 0; --exp) power *= base;
	return power;
}

inline rational r_pow(rational base, int exp) {
	if (exp < 0) {
		base = r_inverse(base);
		exp = -exp;
	}
	base.num = pow_(base.num, exp);
	base.den = pow_(base.den, exp);
	return base;
}

inline rational r_simplify(rational a) {
	if (a.num == 0) {
		a.den = 1;
		return a;
	}
	int g = gcd(a.num, a.den);
	if (g == 0) return a;
	a.num /= g;
	a.den /= g;
	return a;
}
