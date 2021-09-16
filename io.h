#if !defined(IO_H_)
#define IO_H_ 1

int open_filename (const char* filename);

GrB_Info make_mtx_from_file (GrB_Matrix *A_out, GrB_Index * NV_out, GrB_Index * NE_out, int fd);
GrB_Info make_mtx_from_binfile (GrB_Matrix *A_out, GrB_Index * NV_out, GrB_Index * NE_out, int fd);

void make_file_from_mtx (GrB_Matrix A, const char *name, int fd);
void make_binfile_from_mtx (GrB_Matrix A, const char *name, int fd);

#endif
