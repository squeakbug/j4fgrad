#include <string.h>
#include <stdlib.h>

#include <mpi.h>
#include <omp.h>

#include "libmatrix/app.h"
#include "libmatrix/task.h"
#include "libmatrix/matrix.h"

struct mat2d_inv_task {
    size_t sent_rows;
    struct mat2d *forward_mat;
    struct mat2d *reverse_mat;
};

int mat2d_inv_MPI_v2(struct mat2d_inv_task *task) {
    struct mat2d *mat = task->forward_mat;
    struct mat2d *inv = task->reverse_mat;
    size_t row = 0;
    size_t row_in_shard = 0;
    size_t master_indx = row % mat2d_app_get_size();
    bool is_master = master_indx == mat2d_app_get_rank();

    double *row_data = malloc(sizeof(double) * mat2d_get_cols(mat));
    if (row_data == NULL) {
        return -1;
    }
    double *inv_row_data = malloc(sizeof(double) * mat2d_get_cols(mat));
    if (inv_row_data == NULL) {
        return -1;
    }
    while (row < mat2d_get_cols(mat)) {
        if (is_master) {
            double diag = mat2d_get(mat, row_in_shard, row);
            #pragma omp parallel for schedule(dynamic)
            for (size_t j = 0; j < mat2d_get_cols(mat); j++) {
                mat2d_set(mat, row_in_shard, j, mat2d_get(mat, row_in_shard, j) / diag);
                mat2d_set(inv, row_in_shard, j, mat2d_get(inv, row_in_shard, j) / diag);
            }
            double *row_data = mat2d_get_row_ref(mat, row_in_shard);
            MPI_Bcast(row_data, mat2d_get_cols(mat), MPI_DOUBLE, master_indx, MPI_COMM_WORLD);
            double *inv_row_data = mat2d_get_row_ref(inv, row_in_shard);
            MPI_Bcast(inv_row_data, mat2d_get_cols(mat), MPI_DOUBLE, master_indx, MPI_COMM_WORLD);
            row_in_shard += 1;
        } else {
            MPI_Bcast(
                row_data,  mat2d_get_cols(mat), 
                MPI_DOUBLE, master_indx, MPI_COMM_WORLD
            );
            MPI_Bcast(
                inv_row_data, mat2d_get_cols(mat), 
                MPI_DOUBLE, master_indx, MPI_COMM_WORLD
            );

            #pragma omp parallel for schedule(dynamic)
            for (size_t i = 0; i < mat2d_get_rows(mat); i++) {
                double factor = mat2d_get(mat, i, row) / row_data[row];
                for (size_t j = 0; j < mat2d_get_cols(mat); j++) {
                    double tmp1 = mat2d_get(mat, i, j) - factor * row_data[j];
                    mat2d_set(mat, i, j, tmp1);
                    double tmp2 = mat2d_get(inv, i, j) - factor * inv_row_data[j];
                    mat2d_set(inv, i, j, tmp2);
                }
            }
        }

        row += 1;
        master_indx = row % mat2d_app_get_size();
        is_master = master_indx == mat2d_app_get_rank();
    }

    free(row_data);
    free(inv_row_data);

    return 0;
}