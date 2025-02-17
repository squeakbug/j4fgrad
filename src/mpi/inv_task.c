#include <string.h>
#include <stdlib.h>

#include <mpi.h>

#include "libmatrix/app.h"
#include "libmatrix/task.h"
#include "libmatrix/matrix.h"

#define MAX_NODE_CNT 100
#define MAX_SHARDES_CNT 100

struct mat2d_inv_task {
    size_t sent_rows;
    struct mat2d *forward_mat;
    struct mat2d *reverse_mat;
};

static size_t get_sent_rows(struct mat2d *mat, size_t i) {
    return mat2d_get_cols(mat) / mat2d_app_get_size()
        + (mat2d_app_get_size() - i + (mat2d_get_cols(mat) % mat2d_app_get_size()) - 1)
        / mat2d_app_get_size();
}

// Возвращает через аргументы размер партиции
int mad2d_app_redistribute_matrix_size(
    struct mat2d_inv_task *task,
    size_t *rows,
    size_t *cols
) {
    size_t global_indx = mat2d_app_get_rank();
    size_t global_size = mat2d_app_get_size();
    size_t root_indx = mat2d_app_get_root_indx();
    if (global_indx == root_indx) {
        struct mat2d *mat = task->forward_mat;
        size_t shard_rows = (mat2d_get_rows(mat) + global_size - 1) / global_size;
        size_t shard_cols = mat2d_get_cols(mat);
        for (size_t i = 0; i < global_size; ++i) {
            if (i != root_indx) {
                size_t sent_rows = get_sent_rows(mat, i);

                MPI_Send(&shard_rows, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
                MPI_Send(&shard_cols, 1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
                MPI_Send(&sent_rows,  1, MPI_UNSIGNED_LONG, i, 0, MPI_COMM_WORLD);
            }
        }
        *rows = shard_rows;
        *cols = shard_cols;
        task->sent_rows = get_sent_rows(mat, root_indx);
    } else {
        MPI_Status status;
        size_t tmp_rows, tmp_cols;
        size_t tmp_sent_rows;

        MPI_Recv(
            &tmp_rows, 1, MPI_UNSIGNED_LONG, 
            root_indx, 0, MPI_COMM_WORLD, &status
        );
        MPI_Recv(
            &tmp_cols, 1, MPI_UNSIGNED_LONG, 
            root_indx, 0, MPI_COMM_WORLD, &status
        );
        MPI_Recv(
            &tmp_sent_rows, 1, MPI_UNSIGNED_LONG, 
            root_indx, 0, MPI_COMM_WORLD, &status
        );
        *rows = tmp_rows;
        *cols = tmp_cols;
        task->sent_rows = tmp_sent_rows;
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

static struct mat2d *mad2d_get_shard(
    struct mat2d *mat,
    size_t shard_indx,
    size_t shardes_cnt
) {
    size_t shard_rows = (mat2d_get_rows(mat) + shardes_cnt - 1) / shardes_cnt;
    struct mat2d *shard = mat2d_create(shard_rows, mat2d_get_cols(mat));
    if (shard == NULL) {
        return NULL;
    }
    for (size_t i = shard_indx; i < mat2d_get_rows(mat); i += shardes_cnt) {
        size_t row_in_shard_indx = i / shardes_cnt;

        double *src_row_ref = mat2d_get_row_ref(mat, i);
        double *dest_row_ref = mat2d_get_row_ref(shard, row_in_shard_indx);
        memcpy(dest_row_ref, src_row_ref, sizeof(double) * mat2d_get_cols(mat));
    }
    return shard;
}

static void mad2d_place_shard(
    struct mat2d *mat,
    struct mat2d *shard,
    size_t shard_indx,
    size_t shardes_cnt
) {
    size_t shard_rows = (mat2d_get_rows(mat) + shardes_cnt - 1) / shardes_cnt;
    for (size_t i = 0; i < mat2d_get_rows(shard); i += 1) {
        size_t row_in_mat_indx = shard_indx + shardes_cnt * i;

        double *src_row_ref = mat2d_get_row_ref(shard, i);
        double *dest_row_ref = mat2d_get_row_ref(mat, row_in_mat_indx);
        memcpy(dest_row_ref, src_row_ref, sizeof(double) * mat2d_get_cols(mat));
    }
}

static int mad2d_app_redistribute_matrix_data(struct mat2d** mat_inout) {
    size_t global_size = mat2d_app_get_size();
    size_t global_indx = mat2d_app_get_rank();
    size_t root_indx = mat2d_app_get_root_indx();
    MPI_Request requests[MAX_SHARDES_CNT];
    struct mat2d* mat = *mat_inout;
    struct mat2d **shards = malloc(sizeof(struct mat2d *) * (global_size - 1));

    if (global_indx == root_indx) {
        size_t pending_indx = 0;
        for (size_t i = 0; i < global_size; ++i) {
            if (i != root_indx) {
                struct mat2d *shard = mad2d_get_shard(mat, i, global_size);
                if (shard == NULL) {
                    return -1;
                }
                shards[pending_indx] = shard;
                size_t shard_size = mat2d_get_cols(shard) * mat2d_get_rows(shard);
                MPI_Isend(
                    mat2d_get_data(shard), shard_size, MPI_DOUBLE, 
                    i, 0, MPI_COMM_WORLD, &requests[pending_indx]
                );
                pending_indx++;
            }
        }
        MPI_Waitall(pending_indx, requests, MPI_STATUSES_IGNORE);
        for (size_t i = 0; i < pending_indx; ++i) {
            mat2d_destroy(shards[i]);
        }

        struct mat2d *shard = mad2d_get_shard(mat, root_indx, global_size);
        mat2d_destroy(*mat_inout);
        *mat_inout = shard;
    } else {
        size_t shard_size = mat2d_get_cols(mat) * mat2d_get_rows(mat);
        MPI_Irecv(
            mat2d_get_data(mat), shard_size, MPI_DOUBLE, 
            root_indx, 0, MPI_COMM_WORLD, &requests[0]
        );
        MPI_Waitall(1, requests, MPI_STATUSES_IGNORE);
    }

    free(shards);

    MPI_Barrier(MPI_COMM_WORLD);
}

int mad2d_app_redistribute_matrix(struct mat2d_inv_task *task) {
    size_t rows, cols;
    int rc = mad2d_app_redistribute_matrix_size(task, &rows, &cols);
    if (rc == -1) {
        return -1;
    }
    size_t global_indx = mat2d_app_get_rank();
    size_t root_indx = mat2d_app_get_root_indx();
    if (global_indx != root_indx) {
        struct mat2d* mat = mat2d_create(rows, cols);
        if (mat == NULL) {
            return -1;
        }
        task->forward_mat = mat;
    }
    mad2d_app_redistribute_matrix_data(&task->forward_mat);

    if (global_indx != root_indx) {
        struct mat2d* inv = mat2d_create(rows, cols);
        if (inv == NULL) {
            return -1;
        }
        task->reverse_mat = inv;
    }
    mad2d_app_redistribute_matrix_data(&task->reverse_mat);
}

int mad2d_app_unite_matrix(struct mat2d_inv_task *task) {
    struct mat2d *mat = task->reverse_mat;
    size_t global_size = mat2d_app_get_size();
    size_t global_indx = mat2d_app_get_rank();
    size_t root_indx = mat2d_app_get_root_indx();

    if (global_indx == root_indx) {
        struct mat2d *result = mat2d_create(
            mat2d_get_cols(mat),
            mat2d_get_cols(mat)
        );
        if (result == NULL) {
            return -1;
        }
        for (size_t i = 0; i < global_size; ++i) {
            size_t sent_rows = get_sent_rows(mat, i);
            struct mat2d *tmp = mat2d_create(sent_rows, mat2d_get_cols(mat));
            if (tmp == NULL) {
                return -1;
            }
            if (i != root_indx) {
                MPI_Status status;
                MPI_Recv(
                    mat2d_get_data(tmp), sent_rows * mat2d_get_cols(mat), MPI_DOUBLE,
                    i, 0, MPI_COMM_WORLD, &status
                );
                
            } else {
                memcpy(
                    mat2d_get_data(tmp), 
                    mat2d_get_data(mat), 
                    sizeof(double) * sent_rows * mat2d_get_cols(mat)
                );
            }
            mad2d_place_shard(result, tmp, i, global_size);
            mat2d_destroy(tmp);
        }
        mat2d_destroy(task->reverse_mat);
        task->reverse_mat = result;
    } else {
        MPI_Status status;
        MPI_Send(
            mat2d_get_data(mat), task->sent_rows * mat2d_get_cols(mat), MPI_DOUBLE,
            root_indx, 0, MPI_COMM_WORLD
        );
    }
}

int mat2d_inv_MPI_v1(struct mat2d_inv_task *task) {
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
            MPI_Bcast(row_data, mat2d_get_cols(mat), MPI_DOUBLE, master_indx, MPI_COMM_WORLD);
            MPI_Bcast(inv_row_data, mat2d_get_cols(mat), MPI_DOUBLE, master_indx, MPI_COMM_WORLD);

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

struct mat2d* mat2d_app_get_forward_matrix(struct mat2d_inv_task *task) {
    return task->forward_mat;
}

struct mat2d* mat2d_app_get_reverse_matrix(struct mat2d_inv_task *task) {
    return task->reverse_mat;
}

void mat2d_app_set_forward_matrix(struct mat2d_inv_task *task, struct mat2d* mat) {
    task->forward_mat = mat;
}

void mat2d_app_set_reverse_matrix(struct mat2d_inv_task *task, struct mat2d* mat) {
    task->reverse_mat = mat;
}
