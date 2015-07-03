//gcc-5 -fopenmp -O3 viho.c -I/usr/local/include/libiomp -I/usr/X11/include -I/Users/Hhhh/ENS/Stage_L3_math/homographies/code/jpeg-6b -L/usr/X11/lib -lfftw3 -lX11 -L/usr/local/Cellar/libtiff/4.0.3 -ltiff -ljpeg -lpng


#include "iio.c"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "decomp.h"
#include <time.h>

#define WOUT 512
#define HOUT 512




/**
  * this program warps an image by an homography
  * the warping is such that
  *
  * - maps to -
  * x   -->   a
  * y   -->   b
  * z   -->   c
  * t   -->   d
  */



// compute the inverse homography (inverse of a 3x3 matrix)
double invert_homography(double invH[3][3], double H[3][3]){
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0], *r = invH[0];
	double det = a[0]*a[4]*a[8] + a[2]*a[3]*a[7] + a[1]*a[5]*a[6]
		   - a[2]*a[4]*a[6] - a[1]*a[3]*a[8] - a[0]*a[5]*a[7];
	r[0] = (a[4]*a[8]-a[5]*a[7])/det;
	r[1] = (a[2]*a[7]-a[1]*a[8])/det;
	r[2] = (a[1]*a[5]-a[2]*a[4])/det;
	r[3] = (a[5]*a[6]-a[3]*a[8])/det;
	r[4] = (a[0]*a[8]-a[2]*a[6])/det;
	r[5] = (a[2]*a[3]-a[0]*a[5])/det;
	r[6] = (a[3]*a[7]-a[4]*a[6])/det;
	r[7] = (a[1]*a[6]-a[0]*a[7])/det;
	r[8] = (a[0]*a[4]-a[1]*a[3])/det;
	return det;
}

// C = AoB, composition of two homographies (product of 3x3 matrices)
void compose_homographies(double C[3][3], double A[3][3], double B[3][3]){
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 3; j++)
	{
		C[i][j] = 0;
		for (int k = 0; k < 3; k++)
			C[i][j] += A[i][k] * B[k][j];
	}
}



// Find the homography that changes the canonical projective basis into the given four points (x, y, z, w)
void homography_from_four_points(double H[3][3], double x[2], double y[2], double z[2], double w[2]){
	// fix the degree of freedom (assuming the four points are finite)
	double t = 1;

	// translation coefficients
	double p = x[0];
	double q = x[1];

	// "core" 2x2 system
	double A = w[0] - z[0];
	double B = y[0] - z[0];
	double C = w[1] - z[1];
	double D = y[1] - z[1];
	double P = z[0] - y[0] - w[0] + p;
	double Q = z[1] - y[1] - w[1] + q;
	double DET = A * D - B * C;
	double r = (D * P - B * Q) / DET;
	double s = (A * Q - C * P) / DET;
	if (!isnormal(DET))
		fprintf(stderr, "denormal! DET = %g\n", DET);

	// solve the rest of the diagonal system
	double a = w[0] * ( 1 + r ) - p;
	double b = y[0] * ( 1 + s ) - p;
	double c = w[1] * ( 1 + r ) - q;
	double d = y[1] * ( 1 + s ) - q;

	// fill-in the output
	H[0][0] = a; H[0][1] = b; H[0][2] = p;
	H[1][0] = c; H[1][1] = d; H[1][2] = q;
	H[2][0] = r; H[2][1] = s; H[2][2] = t;
}

// Find the homography that moves the four points (x,y,z,t) to (a,b,c,d)
void homography_from_eight_points(double H[3][3], double x[2], double y[2], double z[2], double w[2], double a[2], double b[2], double c[2], double d[2]){
	double H1[3][3], H2[3][3], iH1[3][3];
	homography_from_four_points(H1, x, y, z, w);
	homography_from_four_points(H2, a, b, c, d);
	invert_homography(iH1, H1);
	compose_homographies(H, H2, iH1);
}



int main(int argc,char *argv[]){
	if (argc != 18) {
		printf("usage :\n\t[image.png] x0 x1 y0 y1 z0 z1 t0 t1 a0 a1 b0 b1 c0 c1 d0 d1\n");
		return 1;
	}



	//image input
	char *filename_in = argv[1];
	float *img;
	int w,h,pd;
	img = iio_read_image_float_vec(filename_in, &w, &h, &pd);
	float *img_f = malloc(3*WOUT*HOUT*sizeof(float));



	//homography input
	double x[2] = {strtod(argv[2],NULL), strtod(argv[3],NULL)},
        y[2] = {strtod(argv[4],NULL), strtod(argv[5],NULL)},
        z[2] = {strtod(argv[6],NULL), strtod(argv[7],NULL)},
        t[2] = {strtod(argv[8],NULL), strtod(argv[9],NULL)},
        a[2] = {strtod(argv[10],NULL), strtod(argv[11],NULL)},
        b[2] = {strtod(argv[12],NULL), strtod(argv[13],NULL)},
        c[2] = {strtod(argv[14],NULL), strtod(argv[15],NULL)},
        d[2] = {strtod(argv[16],NULL), strtod(argv[17],NULL)};
	double H[3][3];
	homography_from_eight_points(H,a,b,c,d,x,y,z,t);
	//this H is such that H(a)=x, H(b)=y,...
	//but apply_homo_final resamples by the inverse homography



    //resample (and time)
	clock_t debutcpu,fincpu;
	double debutreal,finreal;
	debutcpu = clock();
	debutreal = omp_get_wtime();
	if(pd==3){
        apply_homo_final(img,img_f,w,h,WOUT,HOUT,H);
	}else{//suppose pd=1
        float *img3 = malloc(3*w*h*sizeof(float));
        for(int i=0;i<w*h;i++){
            for(int l = 0;l<3;l++){
                img3[3*i+l]=img[i];
            }
        }
        apply_homo_final(img3,img_f,w,h,WOUT,HOUT,H);
	}

	fincpu = clock();
	finreal = omp_get_wtime();
	printf("cputime :%fs\ntime : %fs\n",(double)(fincpu-debutcpu)/CLOCKS_PER_SEC,(double)(finreal-debutreal));



	//output
	iio_save_image_float_vec("img_f.png",img_f,WOUT,HOUT,3);
	free(img);
	free(img_f);



	return 0;
}
