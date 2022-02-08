#if !defined(GENERATOR_HEADER_)
#define GENERATOR_HEADER_

void make_edge (int64_t, int64_t * restrict, int64_t * restrict,
		uint8_t * restrict);
void make_edge_endpoints (int64_t, int64_t * restrict, int64_t * restrict);
void edge_list (int64_t * restrict, int64_t * restrict, uint8_t * restrict,
		const int64_t, const int64_t);
void edge_list_64 (int64_t * restrict, int64_t * restrict, uint64_t * restrict,
                   const int64_t, const int64_t);
void edge_list_aos_64 (int64_t * restrict,
                       const int64_t, const int64_t);

#endif /* GENERATOR_HEADER_ */
