#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "libmatrix/matrix.h"

#define EPS 1e-6
#define MAX_MATRIX_SIZE 1000000

struct mat2d {
    double *data;
    size_t rows, cols;
};

struct mat2d* mat2d_create(size_t rows, size_t cols) {
    struct mat2d* mat = calloc(sizeof(struct mat2d) + sizeof(double) * rows * cols, 1);
    assert(mat != NULL);
    mat->data = (double*)(mat + 1);
    mat->rows = rows;
    mat->cols = cols;
    return mat;
}

int mat2d_clone(mat2d **out, mat2d *in) {
    struct mat2d* result = mat2d_create(in->rows, in->cols);
    assert(result != NULL);
    memcpy(result->data, in->data, sizeof(double) * in->rows * in->cols);
    *out = result;
    return 0;
}

struct range2d {
    size_t start_row, start_col; 
    size_t rows, cols;
};

int mat2d_make_submatrix(struct mat2d** outmat, struct mat2d* inmat, struct range2d range) {
    struct mat2d* out = mat2d_create(range.rows, range.cols);
    assert(out!= NULL);

    assert(range.start_row + range.rows <= inmat->rows);
    assert(range.start_col + range.cols <= inmat->cols);
    for (size_t i = range.start_row, j = 0; j < range.rows; i++, j++) {
        memcpy(
            &out->data[j * range.cols], 
            &inmat->data[i * inmat->cols + range.start_col], 
            sizeof(double) * range.cols
        );
    }

    *outmat = out;
    return 0;
}

bool mat2d_eq(struct mat2d* left, struct mat2d* right) {
    assert(left != NULL && right != NULL);
    assert(left->rows == right->rows && right->cols == right->cols);

    for (size_t i = 0; i < left->rows; i++) {
        for (size_t j = 0; j < left->cols; j++) {
            if (fabs(mat2d_get(left, i, j) - mat2d_get(right, i, j)) > EPS) {
                return false;
            }
        }
    }
    return true;
}

int mat2d_inv(struct mat2d** out, struct mat2d* in) {
    mat2d* tmp = NULL;
    mat2d_clone(&tmp, in);

    mat2d* inverse = mat2d_create(in->rows, in->rows);
    mat2d_fill_eye(inverse);

    size_t n = in->rows;
    for (size_t i = 0; i < n; i++) {
        double diag = mat2d_get(tmp, i, i);
        if (diag < EPS && diag > -EPS) {
            return 0;
        }
        for (size_t j = 0; j < n; j++) {
            mat2d_set(tmp, i, j, mat2d_get(tmp, i, j) / diag);
            mat2d_set(inverse, i, j, mat2d_get(inverse, i, j) / diag);
        }

        for (size_t k = 0; k < n; k++) {
            if (k != i) {
                double factor = mat2d_get(tmp, k, i);
                for (size_t j = 0; j < n; j++) {
                    double tmp1 = mat2d_get(tmp, k, j) - factor * mat2d_get(tmp, i, j);
                    mat2d_set(tmp, k, j, tmp1);
                    double tmp2 = mat2d_get(inverse, k, j) - factor * mat2d_get(inverse, i, j);
                    mat2d_set(inverse, k, j, tmp2);
                }
            }
        }
    }

    mat2d_destroy(tmp);
    *out = inverse;
    return 0;
}

int mat2d_T(struct mat2d** out, struct mat2d* in) {
    struct mat2d* tmp = mat2d_create(in->rows, in->cols);
    if (tmp == NULL) {
        return -1;
    }

    for (size_t i = 0; i < in->rows; i++) {
        for (size_t j = 0; j < in->cols; i++) {
            mat2d_set(tmp, j, i, mat2d_get(in, i, j));
        }
    }

    *out = tmp;
    return 0;
}

int mat2d_dot(struct mat2d **out, struct mat2d* left, struct mat2d* right) {
    assert(left->cols == right->rows);

    mat2d* result = mat2d_create(left->rows, right->cols);
    assert(result != NULL);

    for (size_t i = 0; i < left->rows; i++) {
        for (size_t j = 0; j < right->cols; j++) {
            double acc = 0.0;
            for (size_t k = 0; k < left->cols; ++k) {
                acc += mat2d_get(left, i, k) * mat2d_get(right, k, j);
            }
            mat2d_set(result, i, j, acc);
        }
    }

    *out = result;
    return 0;
}

size_t mat2d_get_size(mat2d *mat) {
    return sizeof(struct mat2d) + sizeof(double) * mat->rows * mat->cols;
}

size_t mat2d_get_cols(mat2d *mat) {
    return mat->cols;
}

size_t mat2d_get_rows(mat2d *mat) {
    return mat->rows;
}

double* mat2d_get_data(mat2d *mat) {
    return mat->data;
}

double mat2d_get(struct mat2d *mat, size_t indx1, size_t indx2) {
    assert(indx1 < mat->rows);
    assert(indx2 < mat->cols);
    return mat->data[indx1 * mat->cols + indx2];
}

double* mat2d_get_row_ref(mat2d *mat, size_t row) {
    assert(row < mat->rows);
    return (double*)&mat->data[row * mat->cols];
}

double* mat2d_get_row_cloned(mat2d *mat, size_t row) {
    assert(row < mat->rows);
    size_t row_size = sizeof(double) * mat->cols;
    double *tmp = malloc(row_size);
    memcpy(tmp, &mat->data[row * mat->cols], row_size);
    return tmp;
}

void mat2d_set(struct mat2d *mat, size_t indx1, size_t indx2, double value) {
    assert(indx1 < mat->rows);
    assert(indx2 < mat->cols);
    mat->data[indx1 * mat->cols + indx2] = value;
}

// mat maybe null
void mat2d_destroy(struct mat2d *mat) {
    free(mat);
}

void mat2d_debug(mat2d *mat, FILE *file) {
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            fprintf(file, "%8.3lf ", mat2d_get(mat, i, j));
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

void mat2d_debug_console(mat2d *mat) {
    mat2d_debug(mat, stdout);
}

void mat2d_fill_random(struct mat2d *mat) {
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            mat2d_set(mat, i, j, (double)random() / RAND_MAX);
        }
    }
}

void mat2d_fill_eye(struct mat2d *mat) {
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            mat2d_set(mat, i, j, i == j ? 1.0 : 0.0);
        }
    }
}

void mat2d_fill_zero(struct mat2d *mat) {
    mat2d_fill_value(mat, 0.0);
}

void mat2d_fill_one(struct mat2d *mat) {
    mat2d_fill_value(mat, 1.0);
}

void mat2d_fill_value(struct mat2d *mat, double value) {
    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            mat2d_set(mat, i, j, value);
        }
    }
}

int mat2d_read_from_file(mat2d **out, const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        return -1;
    }

    size_t rows = 0, cols = 0;
    char line[1024];

    while (fgets(line, sizeof(line), file)) {
        if (rows == 0) {
            char *token = strtok(line, ",");
            while (token) {
                cols++;
                token = strtok(NULL, ",");
            }
        }
        rows++;
    }

    struct mat2d *mat = mat2d_create(rows, cols);
    if (mat == NULL) {
        perror("Memory allocation failed");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    rewind(file);

    size_t row_index = 0;
    double *data = mat->data;
    while (fgets(line, sizeof(line), file)) {
        char *token = strtok(line, ",");
        size_t col_index = 0;
        while (token) {
            data[row_index * cols + col_index] = atof(token);
            token = strtok(NULL, ",");
            col_index++;
        }
        row_index++;
    }

    fclose(file);

    *out = mat;
    return 0;
}

int mat2d_write_to_text_file(const struct mat2d *mat, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Failed to open file");
        return -1;
    }

    for (size_t i = 0; i < mat->rows; i++) {
        for (size_t j = 0; j < mat->cols; j++) {
            fprintf(file, "%8.3f", mat->data[i * mat->cols + j]);
            if (j < mat->cols - 1) {
                fprintf(file, " ");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
    return 0;
}

int mat2d_read_from_binfile(struct mat2d **out, const char *filename) {
    FILE *file = fopen(filename, "rb");
    if (!file) {
        return -1;
    }

    size_t rows, cols;
    if (fread(&rows, sizeof(size_t), 1, file) != 1 || 
        fread(&cols, sizeof(size_t), 1, file) != 1) {
        fclose(file);
        return -1;
    }

    struct mat2d *mat = mat2d_create(rows, cols);
    if (mat == NULL) {
        perror("Memory allocation failed");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    size_t elements = rows * cols;
    if (fread(mat->data, sizeof(double), elements, file) != elements) {
        free(mat);
        fclose(file);
        return -1;
    }

    fclose(file);

    *out = mat;
    return 0;
}

int mat2d_write_to_binfile(const struct mat2d *mat, const char *filename) {
    FILE *file = fopen(filename, "wb");
    if (!file) {
        perror("Failed to open file");
        return -1;
    }

    if (fwrite(&mat->rows, sizeof(size_t), 1, file) != 1 ||
        fwrite(&mat->cols, sizeof(size_t), 1, file) != 1) {
        perror("Error writing matrix dimensions");
        fclose(file);
        return -1;
    }

    size_t elements = mat->rows * mat->cols;
    if (fwrite(mat->data, sizeof(double), elements, file) != elements) {
        perror("Error writing matrix data");
    }

    fclose(file);
    return 0;
}
