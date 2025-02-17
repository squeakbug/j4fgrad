#include <stdio.h>
#include <stddef.h>

#include <mpi.h>

#include "libmatrix/app.h"

#define MPI_MAX_PROCESSOR_NAME 1000

struct mat2d_app {
    size_t global_size;
    size_t global_indx;
    size_t root_indx;
};

static struct mat2d_app app;

int mat2d_app_init(int argc, char **argv) {
    int global_indx, global_size, nlen;
    char name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &global_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_indx);
    MPI_Get_processor_name(name, &nlen);

    app.global_size = global_size;
    app.global_indx = global_indx;
    app.root_indx = 0;

    if (global_size < 2) {
        printf("Not enough global size\n");
        return -1;
    }

    printf("Hello from host %s[%d] %d of %d\n", name, nlen, global_indx, global_size);

    return 0;
}

int mat2d_app_get_rank() {
    return app.global_indx;
}

int mat2d_app_get_size() {
    return app.global_size;
}

int mat2d_app_get_root_indx() {
    return app.root_indx;
}

void mat2d_app_destroy() {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}
