typedef struct rational {
	int num, den;
} rational;
rational r_add(rational a, rational b);
rational r_subtract(rational a, rational b);
rational r_multiply(rational a, rational b);
rational r_divide(rational a, rational b);
rational r_negate(rational a);
rational r_inverse(rational a);
rational r_simplify(rational a);
