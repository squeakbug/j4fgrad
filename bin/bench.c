#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include <time.h>
#include <unistd.h>
#include <mpi.h>

#include "libmatrix/app.h"
#include "libmatrix/matrix.h"
#include "libmatrix/task.h"

enum filetype {
    FT_TEXT,
    FT_BIN
};

struct bench_cfg {
    enum filetype ft;
    size_t repeat_cnt;
    const char *in_filename;
    const char *out_filename;
};

uint64_t get_time_ns() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (uint64_t)ts.tv_sec * 1e9 + ts.tv_nsec;
}

int bench_inv_single(struct bench_cfg cfg) {
    struct mat2d *mat_in;
    switch (cfg.ft) {
    case FT_BIN:
        mat2d_read_from_binfile(&mat_in, cfg.in_filename);
        break;
    case FT_TEXT:
        mat2d_read_from_file(&mat_in, cfg.in_filename);
        break;
    };

    uint64_t start = get_time_ns();
    for (size_t i = 0; i < cfg.repeat_cnt; ++i) {
        struct mat2d *inv = NULL;
        mat2d_inv(&inv, mat_in);
    }
    printf("end = %8.3lf\n", (double)(get_time_ns() - start) / cfg.repeat_cnt);

    switch (cfg.ft) {
    case FT_BIN:
        mat2d_write_to_binfile(mat_in, cfg.out_filename);
        break;
    case FT_TEXT:
        mat2d_write_to_text_file(mat_in, cfg.out_filename);
        break;
    };

    return 0;
}

int bench_inv(struct bench_cfg cfg) {
    int my_rank = mat2d_app_get_rank();

    struct mat2d *mat_in;
    switch (cfg.ft) {
    case FT_BIN:
        mat2d_read_from_binfile(&mat_in, cfg.in_filename);
        break;
    case FT_TEXT:
        mat2d_read_from_file(&mat_in, cfg.in_filename);
        break;
    };

    uint64_t start = get_time_ns();
    for (size_t i = 0; i < cfg.repeat_cnt; ++i) {
        struct mat2d *mat_out;
        int rc = mat2d_inv(&mat_out, mat_in);
        if (rc == -1) {
            printf("main: %s\n", strerror(errno));
        }
    }
    printf("end = %8.3lf\n", (double)(get_time_ns() - start) / cfg.repeat_cnt);

    return 0;
}

int parse_filetype(const char *str, enum filetype *ft) {
    if (strcmp(str, "b") == 0) {      
        *ft = FT_BIN;
    } else if (strcmp(str, "t") == 0) {
        *ft = FT_TEXT;
    } else {
        return -1;
    }
    return 0;
}

int parse_repeat_cnt(const char *str, size_t *repeats) {
    char *ep = NULL;
    long tmp = strtol(str, &ep, 10);
    if (ep == str + strlen(str)) {
        if (tmp < 0) {
            return -1;
        } else {
            *repeats = (size_t)tmp;
        }
    } else {
        return -1;
    }
    return 0;
}

int main(int argc, char **argv)
{
    if (argc < 5) {
        printf("usage: bench.elf <repeats> <filetype> <input_filename> <output_filename>");
        return -1;
    }

    int rc = 0;
    if ((rc = mat2d_app_init(argc, argv)) != 0) {
        printf("Goog by dpi\n");
        mat2d_app_destroy();
        return -1;
    }

    size_t repeat_cnt;
    rc = parse_repeat_cnt(argv[1], &repeat_cnt);
    if (rc == -1) {
        printf("Error while parsing repeat count\n");
        return -1;
    }

    enum filetype ft;
    rc = parse_filetype(argv[2], &ft);
    if (rc == -1) {
        printf("Error while parsing filetype\n");
        return -1;
    }

    struct bench_cfg cfg = {
        .ft = ft,
        .repeat_cnt = repeat_cnt,
        .in_filename = argv[3],
        .out_filename = argv[4]  
    };

    bench_inv(cfg);

    mat2d_app_destroy();

out:
    return rc;
}