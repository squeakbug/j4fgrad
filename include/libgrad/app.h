#ifndef APP_H
#define APP_H

int mat2d_app_init(int argc, char **argv);
void mat2d_app_destroy();

int mat2d_app_get_rank();
int mat2d_app_get_size();
int mat2d_app_get_root_indx();

#endif