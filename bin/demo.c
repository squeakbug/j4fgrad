#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include <time.h>
#include <unistd.h>

#include "libmatrix/matrix.h"

int main(int argc, char **argv)
{
    int rc = 0;
    if ((rc = mat2d_app_init(argc, argv)) != 0) {
        printf("main: %s\n", strerror(errno));
        goto out;
    }

    struct mat2d *mat_in;
    const char *filetype = argv[1];
    const char *in_filename = argv[2];
    const char *out_filename = argv[3];
    if (strcmp(filetype, "b") == 0) {
        rc = mat2d_read_from_binfile(&mat_in, in_filename);
        if (rc == -1) {
            printf("main: %s\n", strerror(errno));
            goto out;
        }
    } else if (strcmp(filetype, "t") == 0) {
        rc = mat2d_read_from_file(&mat_in, in_filename);
        if (rc == -1) {
            printf("main: %s\n", strerror(errno));
            goto out;
        }
    } else {
        printf("main: bad args\n");
        goto out;
    }

    struct mat2d *mat_out;
    rc = mat2d_inv(&mat_out, mat_in);
    if (rc == -1) {
        printf("main: %s\n", strerror(errno));
        goto out;
    }

    if (strcmp(filetype, "b") == 0) {
        rc = mat2d_write_to_binfile(mat_out, out_filename);
        if (rc == -1) {
            printf("main: %s\n", strerror(errno));
            goto out;
        }
    } else if (strcmp(filetype, "t") == 0) {
        rc = mat2d_write_to_text_file(mat_out, out_filename);
        if (rc == -1) {
            printf("main: %s\n", strerror(errno));
            goto out;
        }
    } else {
        printf("main: bad args\n");
        goto out;
    }

    mat2d_destroy(mat_in);
    mat2d_app_destroy();
    return rc;

free_inmat:
    mat2d_destroy(mat_in);

out:
    return -1;
}