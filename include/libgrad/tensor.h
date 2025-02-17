#ifndef MATRIX_H
#define MATRIX_H

#include <stdint.h>
#include <stdbool.h>
#include <stdio.h>

#define MPI_ENABLED

typedef struct mat2d mat2d;

mat2d* mat2d_create(size_t n, size_t k);
int mat2d_clone(mat2d **out, mat2d *in);
void mat2d_destroy(mat2d *mat);

double mat2d_get(mat2d *mat, size_t indx1, size_t indx2);
size_t mat2d_get_rows(mat2d *mat);
size_t mat2d_get_cols(mat2d *mat);
size_t mat2d_get_size(mat2d *mat);
double* mat2d_get_data(mat2d *mat);
double* mat2d_get_row_ref(mat2d *mat, size_t row);
double* mat2d_get_row_cloned(mat2d *mat, size_t row);
void mat2d_set(mat2d *mat, size_t indx1, size_t indx2, double value);

bool mat2d_eq(struct mat2d* left, struct mat2d* right);

int mat2d_inv(struct mat2d** out, struct mat2d* in);
int mat2d_T(struct mat2d** out, struct mat2d* in);
int mat2d_dot(struct mat2d** out, struct mat2d* left, struct mat2d* right);

void mat2d_fill_random(struct mat2d *mat);
void mat2d_fill_eye(struct mat2d *mat);
void mat2d_fill_zero(struct mat2d *mat);
void mat2d_fill_one(struct mat2d *mat);
void mat2d_fill_value(struct mat2d *mat, double value);

int mat2d_read_from_file(mat2d **out, const char *filename);
int mat2d_read_from_binfile(mat2d **out, const char *filename);
int mat2d_write_to_text_file(const mat2d *mat, const char *filename);
int mat2d_write_to_binfile(const mat2d *mat, const char *filename);

void mat2d_debug(mat2d *mat, FILE *file);
void mat2d_debug_console(mat2d *mat);

#endif