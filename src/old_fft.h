#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define PI 3.1415926535897932384

#define FFT_FORWARD 1
#define FFT_BACKWARD 2
#define FFT_ESTIMATE 3

// debug
#define MAXIMUM_FFTL 4194304

typedef double old_fft_complex[2];
typedef struct {
	int n;
	int sign;
	unsigned int flags;
	old_fft_complex *c_in;
	double *in;
	old_fft_complex *c_out;
	double *out;
	double *input;
	int *ip;
	double *w;
} old_fft_plan;

old_fft_plan old_fft_plan_dft_1d(int n, old_fft_complex *in, old_fft_complex *out, int sign, unsigned int flags);
old_fft_plan old_fft_plan_dft_c2r_1d(int n, old_fft_complex *in, double *out, unsigned int flags);
old_fft_plan old_fft_plan_dft_r2c_1d(int n, double *in, old_fft_complex *out, unsigned int flags);
void old_fft_execute(const old_fft_plan p);
void old_fft_destroy_plan(old_fft_plan p);

// hidden functions
void old_rdft(int n, int isgn, double *a, int *ip, double *w);
void old_cdft(int n, int isgn, double *a, int *ip, double *w);