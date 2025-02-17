#include <stdio.h>
#include <stddef.h>

#include <mpi.h>

#include "libmatrix/app.h"

#define MPI_MAX_PROCESSOR_NAME 1000

struct mat2d_app {
    
};

static struct mat2d_app app;

int mat2d_app_init(int argc, char **argv) {
    return 0;
}

int mat2d_app_get_rank() {
    return 0;
}

int mat2d_app_get_size() {
    return 1;
}

int mat2d_app_get_root_indx() {
    return 0;
}

void mat2d_app_destroy() {
    
}
