compile with openMP

On suppose pd=3

#include "decomp.h", qui se charge des autres include

Appeller apply_homography_final(float *img_depart, float *img_arrivé, int w_depart, int h_depart, int w_final, int h_final, homography double H[3][3])

dans la version non interactive on donne l'image puis les 9 coeffs, il écrit le résultat dans une nouvelle image (au format .pgm)

gcc-5 -I/usr/X11/include -fopenmp viho.c iio.c ftr.c -I/usr/local/include/libiomp -L/usr/X11/lib -lfftw3 -lX11
