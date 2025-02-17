#include <assert.h>

#include "libmatrix/app.h"
#include "libmatrix/matrix.h"
#include "libmatrix/task.h"

int mat2d_dot_mpi(struct mat2d **out, struct mat2d* left, struct mat2d* right) {
    assert(mat2d_get_cols(left) == mat2d_get_rows(right));
    assert((mat2d_app_get_size() - 1) % 4 == 0);

    mat2d* result = mat2d_create(
        mat2d_get_rows(left), 
        mat2d_get_cols(right)
    );
    assert(result != NULL);

    return 0;
}