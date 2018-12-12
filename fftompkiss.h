/* This file is automatically generated. DO NOT EDIT! */

#ifndef _fftompkiss_h
#define _fftompkiss_h


int fft_init(int pad1 /* padding on the first axis */,
        int nz_,   int nx_,   int ny_ /* input data size */, 
        int *nz2, int *nx2, int *ny2 /* padded data size */,
        bool rtoc /* real to complex fft flag */,
        bool padio /* inp and out are padded*/);
/*< initialize >*/


void fft(void *inp /* [n1*n2*n3] */, 
        sf_complex *out /* [n1*n2*n3] */);
/*< 3-D FFT >*/


void ifft(void *out /* [n1*n2*n3] */, 
        sf_complex *inp /* [n1*n2*n3] */);
/*< 3-D inverse FFT >*/


void fft_finalize();
/*< clean up fft >*/

#endif
