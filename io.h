#if !defined(IO_H_)
#define IO_H_ 1

GrB_Info make_mtx_from_file (GrB_Matrix *A_out, GrB_Index * NV_out, GrB_Index * NE_out, FILE *f);
GrB_Info make_mtx_from_binfile (GrB_Matrix *A_out, GrB_Index * NV_out, GrB_Index * NE_out, FILE *f);

void make_file_from_mtx (GrB_Matrix A, const char *name, FILE *f);
void make_binfile_from_mtx (GrB_Matrix A, const char *name, FILE *f);

#endif
