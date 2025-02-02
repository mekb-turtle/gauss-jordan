#include <stdbool.h>
#include <stddef.h>

typedef struct rational {
	// numerator and denominator
	int num, den;
} rational;

int gcd(int a, int b);
int lcm(int a, int b);

// set denominator to LCM of both denominators
rational r_set_denominator(rational a, int den);

// writes non-reoccurring part and reoccurring part of a fraction using .write
// .sig_f is the number of significant figures to write
// .write should return false if the function should stop writing and false will be returned
// returns if the function was successful
bool r_decompose(rational frac, bool (*write)(int digit, bool is_reoccuring, void *data), void *data);

// self-explanatory
bool r_is_negative(rational frac);
bool r_is_integer(rational frac);

int r_get_integer(rational frac); // get integer part of fraction

// returns negative if a < b, zero if a == b, positive if a > b
int r_cmp(rational a, rational b);

// checks if fraction is reoccurring, i.e prime factors of the denominator exist other than 2 and 5
bool is_reoccuring(rational frac);

// reduce fraction to lowest numerator and denominator
rational r_simplify(rational a);

rational r_abs(rational a);                                 // absolute value
rational r_add(rational a, rational b);                     // a + b
rational r_subtract(rational minuend, rational subtrahend); // minuend - subtrahend
rational r_multiply(rational a, rational b);                // a * b
rational r_divide(rational dividend, rational divisor);     // dividend / divisor
rational r_negate(rational a);                              // -a
rational r_inverse(rational a);                             // 1 / a
rational r_pow(rational base, int exp);                     // base ^ exp
