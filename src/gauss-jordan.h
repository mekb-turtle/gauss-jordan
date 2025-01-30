typedef void *(*get_func)(int x, int y, void *data);
typedef void (*set_func)(int x, int y, void *value, void *data);
typedef void *(*m_func)(void *a, void *b);

// returns positive number on success, -1 if invalid
int gauss_jordan(get_func get, set_func set, m_func sub, m_func mul, m_func div, int (*nonzero)(void *a), int m, int n, void *data);
