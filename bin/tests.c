#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <time.h>
#include <unistd.h>

#include "libmatrix/matrix.h"

int test_rev() {
    struct mat2d *mat = mat2d_create(2, 2);
    mat2d_fill_eye(mat);
    mat2d_set(mat, 0, 1, 2.0);
    mat2d_set(mat, 1, 0, 2.0);
    mat2d_debug(mat, stdout);

    int rc = 0;
    struct mat2d *inv = NULL;
    if ((rc = mat2d_inv(&inv, mat)) == 0) {
        printf("mat.inv: \n");
        mat2d_debug(inv, stdout);
        printf("\n");
    } else {
        printf("Error while matrix product: %d\n", rc);
    }

    mat2d_destroy(mat);
    mat2d_destroy(inv);

    return rc;
}

int test_dot() {
    struct mat2d *mat1 = mat2d_create(2, 4);
    mat2d_fill_random(mat1);
    mat2d_debug(mat1, stdout);

    struct mat2d *mat2 = mat2d_create(4, 2);
    mat2d_fill_random(mat2);
    mat2d_debug(mat2, stdout);

    int rc = 0;
    struct mat2d *result = NULL;
    if ((rc = mat2d_dot(&result, mat1, mat2)) == 0) {
        printf("mat1 x mat2: \n");
        mat2d_debug(result, stdout);
        printf("\n");
    } else {
        printf("Error while matrix product: %d\n", rc);
    }

    mat2d_destroy(mat1);
    mat2d_destroy(mat2);
    mat2d_destroy(result);

    return rc;
}

int main(int argc, char **argv)
{
    test_rev();
    test_dot();
    return 0;
}