#ifndef TASK_H
#define TASK_H

typedef struct mat2d_inv_task mat2d_inv_task;

int init();
void destroy();

struct mat2d* mat2d_app_get_forward_matrix(struct mat2d_inv_task *task);
struct mat2d* mat2d_app_get_reverse_matrix(struct mat2d_inv_task *task);
void mat2d_app_set_forward_matrix(
    struct mat2d_inv_task *task, 
    struct mat2d* mat
);
void mat2d_app_set_reverse_matrix(
    struct mat2d_inv_task *task, 
    struct mat2d* mat
);
int mad2d_app_redistribute_matrix();
int mad2d_app_unite_matrix();

#endif