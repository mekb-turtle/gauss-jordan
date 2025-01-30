#include "rational.h"
#include <stdio.h>

static inline int gcd(int a, int b) {
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	while (b != 0) {
		int temp = b;
		b = a % b;
		a = temp;
	}
	return a;
}
static inline int lcm(int a, int b) {
	if (a < 0) a = -a;
	if (b < 0) b = -b;
	return (a * b) / gcd(a, b);
}

rational r_add(rational a, rational b) {
	int l = lcm(a.den, b.den);
	a.num *= l / a.den;
	b.num *= l / b.den;
	a.den = l;
	a.num += b.num;
	return r_simplify(a);
}
rational r_subtract(rational a, rational b) {
	return r_add(a, r_negate(b));
}
rational r_multiply(rational a, rational b) {
	a.num *= b.num;
	a.den *= b.den;
	return r_simplify(a);
}
rational r_divide(rational a, rational b) {
	return r_multiply(a, r_inverse(b));
}
rational r_negate(rational a) {
	if (a.den < 0)
		a.den = -a.den;
	else
		a.num = -a.num;
	return a;
}
rational r_inverse(rational a) {
	int t = a.num;
	a.num = a.den;
	a.den = t;
	return a;
}
rational r_simplify(rational a) {
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
