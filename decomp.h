#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "affine.h"
#include "homo_box.h"

/**
  * how are indexed coefficients of matrices :
  *
  * for the homographies
  * 0 1 2
  * 3 4 5
  * 6 7 8
  *
  * for the affinities
  * 0 1 2
  * 3 4 5
  */

//y = H(x)
static void apply_homography(double y[2], double H[3][3], double x[2])
{
	double X = H[0][0] * x[0] + H[0][1] * x[1] + H[0][2];
	double Y = H[1][0] * x[0] + H[1][1] * x[1] + H[1][2];
	double Z = H[2][0] * x[0] + H[2][1] * x[1] + H[2][2];
	y[0] = X / Z;
	y[1] = Y / Z;
}


void smallestRectangle(double A[6], int wIn, int hIn, int *muOut, int *nuOut, int *wOut, int *hOut){
/**
  * @param
  *     A an affine map
  *     [0,wIn]*[0,hIn] a rectangle of the plan
  *     [muOut,muOut+wOut]*[nuOut,nuOut+hOut] an output rectangle
  *
  * @return
  *     [muOut,muOut+wOut]*[nuOut,nuOut+hOut] is the smallest rectangle (with integer-coordinated summits) containing A([0,wIn]*[0,hIn])
  *     (rectangle starting at (muOut,nuOut) with dimensions wOut*hOut)
  */
    //summits of the rectangle [0,wIn]*[0,hIn]
    double xIn[4] = {0,wIn,wIn,0};
    double yIn[4] = {0,0,hIn,hIn};
    //summits of the parallelogram A([0,wIn]*[0,hIn])
    double xOut[4];
    double yOut[4];
    for(int k=0;k<4;k++){
        xOut[k]=xIn[k]*A[0]+yIn[k]*A[1]+A[2];
        yOut[k]=xIn[k]*A[3]+yIn[k]*A[4]+A[5];
    }

    //extrema of the coordinates of the output rectangle
    double xOutMin_double = xOut[0];
    double yOutMin_double = yOut[0];
    double xOutMax_double = xOut[0];
    double yOutMax_double = yOut[0];
    for(int k=1;k<4;k++){
        xOutMin_double = fmin(xOutMin_double,xOut[k]);
        yOutMin_double = fmin(yOutMin_double,yOut[k]);
        xOutMax_double = fmax(xOutMax_double,xOut[k]);
        yOutMax_double = fmax(yOutMax_double,yOut[k]);
    }
    int xOutMin = floor(xOutMin_double);
    int yOutMin = floor(yOutMin_double);
    int xOutMax = ceil(xOutMax_double);
    int yOutMax = ceil(yOutMax_double);

    //output
    *muOut = xOutMin;
    *nuOut = yOutMin;
    *wOut = xOutMax - xOutMin;
    *hOut = yOutMax - yOutMin;
}


// H = A H0 B
void decomp(double H[9],double A[6],double H0[9],double B[6]){
/**
  * @param
  *     H an input homography
  *		H0 an output homography
  *		A,B output affinities
  */
    double a=H[0], b=H[1], p=H[2], c=H[3], d=H[4], q=H[5], r=H[6], s=H[7], t=H[8];
    double t0, t1;
	//assume r,s != 0,0

    //an infinity of (t0,t1) are possible. Here is a simple one
    if(fabs(a*s-b*r)<fabs(c*s-d*r)){
        t0 = 0.;
        t1 = -((a*s-b*r)*(a*r+b*s)+(c*s-d*r)*(c*r+d*s))/(pow(r,2.)+pow(s,2.))/(c*s-d*r);
    }else{
        t0 = -((a*s-b*r)*(a*r+b*s)+(c*s-d*r)*(c*r+d*s))/(pow(r,2.)+pow(s,2.))/(a*s-b*r);;
        t1 = 0.;
    }



    //translation of (t0,t1) (H becomes T_(t0,t1)*H)
    a += t0*r;
    b += t0*s;
    p += t0*t;
    c += t1*r;
    d += t1*s;
    q += t1*t;



    double Nphi = sqrt(pow(r,2.)+pow(s,2.));
    B[0] = r/Nphi;
	B[1] = s/Nphi;
	B[2] = 0.;
	B[3] = -s/Nphi;
	B[4] = r/Nphi;
	B[5] = 0.;

    double Npsi = sqrt(pow(a*s-b*r,2.)+pow(c*s-d*r,2.));
	A[0] = (c*s-d*r)/Npsi;
	A[1] = (a*s-b*r)/Npsi;
	A[2] = -t0;
	A[3] = -(a*s-b*r)/Npsi;
	A[4] = (c*s-d*r)/Npsi;
	A[5] = -t1;

    H0[0] = -(a*d-b*c)*Nphi/Npsi;
    H0[1] = 0.;
    H0[2] = (p*(c*s-d*r)-q*(a*s-b*r))/Npsi;
    H0[3] = 0.;
    H0[4] = -Npsi/Nphi;
    H0[5] = (p*(a*s-b*r)+q*(c*s-d*r))/Npsi;
    H0[6] = Nphi;
    H0[7] = 0.;
    H0[8] = t;

}


void apply_homo_final(float *img,float *img_f,int w,int h,int w_f,int h_f,double H[3][3]){

/**
  * @param
  *     img : input image
  *     img_f : output image
  *     w,h : width and heights of the input image
  *		w_f,h_f : width and heights of the output image
  *     H : inverse homography, such that img_f(x,y)=img(H(x,y))
  */


	double *a = *H;

	if(a[8]!=0 && a[6]/a[8]==0 && a[7]/a[8]==0){    //H is an affinity
		double A[6] = {
            a[0]/a[8], a[1]/a[8], a[2]/a[8],
            a[3]/a[8], a[4]/a[8], a[5]/a[8]
        };
		apply_affinity(img,img_f,w,h,w_f,h_f,A);
	}else{  //H is an homography
		double A[6];
		double H0[9];
		double B[6];

		decomp(a,A,H0,B);   //compute the decomposition H = A H0 B



        //declare the size (w,h) and the position (mu,nu) of every intermediate image
        /*
         * For an image at the position (mu,nu), the coordinates of the (i,j)-th pixel are (i+mu,j+nu)
         * rectangleX is the rectangle of the plan [muX,muX+wX]*[nuX,nuX+hX] containing imgX
         * Thus, rectangle1 is the input image and rectangle_f is the output image
         *
         *
         * img1 is the input image (mu1=nu1=0, w1=w, h1=h)
         * img2 and img3 are the intermediate images
         * img4 is the output image (mu4=nu4=0, w4=w_f, h4=h_f)
         */
        int mu2,nu2,w2,h2,
            mu3,nu3,w3,h3;

		//rectangle2 = invA(rectangle1) where invA = A^(-1)
		double det_A = A[0]*A[4]-A[3]*A[1];
        double invA[6] = {
            A[4]/det_A,
            -A[1]/det_A,
            (A[1]*A[5]-A[4]*A[2])/det_A,
            -A[3]/det_A,
            A[0]/det_A,
            (-A[0]*A[5]+A[3]*A[2])/det_A
        };
        smallestRectangle(invA,w,h,&mu2,&nu2,&w2,&h2);
        if((w2-w)%2!=0){w2++;}  //w2 must have the same parity than w
        if((h2-h)%2!=0){h2++;}  //h2 must have the same parity than h

        //rectangle3 = B(rectangle4)
        smallestRectangle(B,w_f,h_f,&mu3,&nu3,&w3,&h3);
        if((w3-w_f)%2!=0){w3++;}    //w3 must have the same parity than w_f
        if((h3-h_f)%2!=0){h3++;}    //h3 must have the same parity than h_f



        //the affinities must be corrected so that they fit the positions of the rectangles
        /*
         * the exact formulas are
         * A[2] = mu2*A[0] + nu2*A[1] + A[2] - mu1;
         * A[5] = mu2*A[3] + nu2*A[4] + A[5] - nu1;
         * B[2] = mu4*B[0] + nu4*B[1] + B[2] - mu3;
         * B[5] = mu4*B[3] + nu4*B[4] + B[5] - nu3;
         * but a lot of term are zero
         */
        A[2] = mu2*A[0] + nu2*A[1] + A[2];
        A[5] = mu2*A[3] + nu2*A[4] + A[5];
        B[2] = - mu3;
        B[5] = - nu3;



        ///Application of the decomposition

        float *img2 = malloc(3*w2*h2*sizeof(float));
		apply_affinity(img,img2,w,h,w2,h2,A);

		float *img3 = malloc(3*w3*h3*sizeof(float));
		apply_homo(img2,img3,w2,h2,w3,h3,mu2,nu2,mu3,nu3,H0);
		free(img2);

		apply_affinity(img3,img_f,w3,h3,w_f,h_f,B);
		free(img3);
	}


	double p[2];

	//truncate the output image, because it has been symmetrized to prevent ringing
	for(int i=0;i<w_f;i++){
		for(int j=0;j<h_f;j++){
			p[0]=i; p[1]=j;
			apply_homography(p,H,p);
			p[0] = (p[0] - 0.5) * w / (w - 1.0);
			p[1] = (p[1] - 0.5) * h / (h - 1.0);
			if(p[0]<0 || p[0]>w || p[1]<0 || p[1]>h){
				for(int l=0;l<3;l++){img_f[(j*w_f+i)*3+l]=0;}
			}
		}
	}

}
