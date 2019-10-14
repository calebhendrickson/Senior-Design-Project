// fftw1.cpp : This file contains the 'main' function. Program execution begins and ends there.
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <tchar.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <errno.h>
#include <fftw3.h>
#include <sndfile.h>
#include "wavfile.h"
using namespace std;

#define REAL 0
#define IMAG 1
#define M_PI 3.14159265358979323846

#define N0 44100

//COMPUTE 1-D FAST FOURIER TRANSFORM
void fft(double *in, fftw_complex *out, int duration)
{
	// create a DFTplan, this plan only computes N0/2 + 1 vals
	fftw_plan plan = fftw_plan_dft_r2c_1d((N0*duration), in, out, FFTW_ESTIMATE);
	//execute
	fftw_execute(plan);
	//do some cleaning
	fftw_destroy_plan(plan);
	fftw_cleanup();
}


// CURRENTLY NOT NEEDED
// COMPUTE INVERSE FFT
void ifft(fftw_complex* in, double* out, int duration) {
	// create a DFTplan, this plan only computes N0/2 + 1 vals
	fftw_plan plan = fftw_plan_dft_c2r_1d(((N0*duration)), in, out, FFTW_ESTIMATE);
	//execute
	fftw_execute(plan);
	//do some cleaning
	fftw_destroy_plan(plan);
	fftw_cleanup();
}

// displays complex numbers in the form a +/- bi
void displayComplex(fftw_complex *y) {
	for (int i = 0; i < (N0/2) + 1; i++) {
		cout << i << "   ";
		if (y[i][IMAG] < 0)
			cout << y[i][REAL] << " - " << abs(y[i][IMAG]) << "i" << endl;
		else
			cout << y[i][REAL] << " + " << y[i][IMAG] << "i" << endl;
	}
}

//displays frequency vs amplitude table
void displayFreq(fftw_complex *y, double samplerate, int duration) {
	for (int i = 0; i < ((N0 * duration) / 2) + 1; i++) {
		if (i % 20 == 0) {
			cout << (samplerate / (N0 * (double)duration)) * i << "Hz   " << sqrt((y[i][REAL] * y[i][REAL]) + (y[i][IMAG] * y[i][IMAG])) << "i" << endl;
	    }
	}
}

// low pass SINC filter array
void lowPassFilter(double* lpf, double* n, double *w, int duration, int bandwidth) {
	int j = 0;
	double max = 0;
	double fc = 0.1;
	for (int i = 0; i < bandwidth; i++) {
		*(lpf + i) = sin(M_PI * (2 * fc * (((*n + i) - ((double)bandwidth - 1)) / 2))) / M_PI * (2 * fc * (((*n + i) - ((double)bandwidth - 1)) / 2));
	}
	//Compute blackman window
	for (int i = 0; i < bandwidth; i++) {
		*(w + i) = 0.42 - (0.5 * cos((2 * M_PI * (*(n + i))) / ((double)bandwidth - 1))) + (0.08 * cos((4 * M_PI * (*(n + i))) / ((double)bandwidth - 1)));
	}
	for (int i = 0; i < bandwidth; i++) {
		*(lpf + i) = *(lpf + i) * (*(w + i));
	}
	double sum = 0;
	for (int i = 0; i < bandwidth; i++) {
		sum = sum + *(lpf + i);
	}
	//NORMALIZE
	cout << "NORMALIZE" << endl;
	for (int i = 0; i < bandwidth; i++) {
		*(lpf + i) = *(lpf + i) / (sum * 10);
		// not sure why multiplying by ten helps here
		// output coefficients to console to see what we will be convolving with our original signal
		cout << *(lpf + i) << endl;
	}
}

// NO LONGER IN USE
// Hann window function
void hannWindow(double *windowlpf, int duration) {
	// generate hann window coefficient array
	for (int i = 0; i < (N0*duration) - 1; i++) {
		// windowed low pass filter array (more effective than truncating lpf filter coefficients
		*(windowlpf + i) = (0.5 * (1 - cos((2 * M_PI * i) / (((double)44100*duration) - 1))));
		// cout << *(windowlpf + i) << endl;
	}
}

// to convole(apply) the windowed low pass impulse response to the original soundwave
void convolve(double *Signal, int SignalLen,
	double *Kernel, int KernelLen,
	double *Result)
{
	int n;
	for (n = 0; n < SignalLen + KernelLen - 1; n++)
	{
		int kmin, kmax, k;
		*(Result+n) = 0;
		kmin = (n >= KernelLen - 1) ? n - (KernelLen - 1) : 0;
		kmax = (n < SignalLen - 1) ? n : SignalLen - 1;

		for (k = kmin; k <= kmax; k++)
		{
			*(Result + n) += *(Signal+k) * (*(Kernel + (n - k)));
		}
	}
}

// reads discrete pcm data from a wav file, windows the data, and calculates the dft and amplitude vs frequency table
int main()
{
	// INIT
	SF_INFO sfinfoREAD;
	char path[] = "./vtfnv-xbvo6.wav";
	SNDFILE* sndfileREAD;
	sfinfoREAD.format = 0;
	cout << "FORMAT1" << endl;
	cout << sfinfoREAD.format << endl;
	// OPEN SOUND FILE
	sndfileREAD = sf_open(path, SFM_READ, &sfinfoREAD);

	// determine samplerate and duration of audio file
	sf_count_t frames = sfinfoREAD.frames;
	double filesamplerate = sfinfoREAD.samplerate;
	double duration = frames / filesamplerate;

	double samplerate = N0;
	// samplerate = 44,100 per second

	//compute fundamental frequency 1 / T
	// T is the time period to capture one cycle
	double fund_frequency = (1 / (int)duration);

	// determine size of output array
	int n_out = (((N0*(int)duration) / 2) + 1);
	// allocate memory for output arrays
	fftw_complex* y;
	y = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n_out);

	fftw_complex* freqresponse;
	freqresponse = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * n_out);


	double transition = 0.08;
	double tran = ceil(4 / transition);
	int bandwidth = (int)tran;
	if (bandwidth % 2 != 0) {
		bandwidth = bandwidth + 1;
	}

	// allocate and init bandwidth arrays
	double* n;
	n = (double*)fftw_malloc(sizeof(double) * bandwidth);
	for (int i = 0; i < bandwidth; i++) {
		*(n + i) = i;
	}

	double* w;
	w = (double*)fftw_malloc(sizeof(double) * bandwidth);
	for (int i = 0; i < bandwidth; i++) {
		*(w + i) = 0;
	}

	// allocate memory for input arrays
	double *samples;
	samples = (double*)fftw_malloc(sizeof(double) * N0 * (int)duration);

	
	double* lpf;
	lpf = (double*)fftw_malloc(sizeof(double) * bandwidth);


	double* sincOutput;
	sincOutput = (double*)fftw_malloc((sizeof(double) * N0 * (int)duration) + (sizeof(double) * bandwidth) -1);

	// initialize pointers for input arrays
	for (int i = 0; i < N0 * (int)duration; i++) {
		*(samples + i) = 0;
	}

	for (int i = 0; i < bandwidth; i++) {
		*(lpf + i) = 0;
	}

	// read data from the sound file
	sf_read_double(sndfileREAD, samples, 44100 * ((sf_count_t)duration));

	// compute the fft of x and store results in y
	fft(samples, y, (int)duration);

	// display time vs frequency plot(spectrogram)
	cout << "\nFrequency Table (Input src)" << endl;
	displayFreq(y, samplerate, (int)duration);

	// low pass filter function
	lowPassFilter(lpf, n, w, (int)duration, bandwidth);

	int sampleSize = N0 * (int)duration;
	int filterSize = bandwidth;
	convolve(samples, sampleSize, lpf, bandwidth, sincOutput);

	fft(sincOutput, freqresponse, (int)duration);
	//technically should be duration + bandwidth


	cout << "\nFrequency Response FFT =" << endl;
	
	// SEEMS TO BE RETURNING CORRECT RESPONSE
	// display frequency response of the fft of the windowed lpf
	displayFreq(freqresponse, samplerate, (int)duration);

	// CHECKING FORMAT
	cout << "FRAMES" << endl;
	cout << frames << endl;
	cout << "DURATION" << endl;
	cout << (sf_count_t)duration << endl;
	cout << "SAMPLERATE" << endl;
	cout << samplerate << endl;
	cout << "FORMAT2" << endl;
	cout << sfinfoREAD.format << endl;
	
	for (int i = 0; i < N0 * (int)duration; i++) {
		// for testing purposes
		if (i % 10000) {
			cout << *(samples + i) << endl;
		}
	}

	SF_INFO sfinfoWRITE;
	SNDFILE* sndfileWRITE;
	sfinfoWRITE.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16;
	cout << "FORMAT3" << endl;
	cout << sfinfoWRITE.format << endl;
	sfinfoWRITE.samplerate = N0;
	sfinfoWRITE.channels = 1;
	sndfileWRITE = sf_open("./outputFile.wav", SFM_WRITE, &sfinfoWRITE);
	
	sf_write_double(sndfileWRITE, sincOutput, N0 * (sf_count_t)duration);
	

	// wrap up
	sf_close(sndfileREAD);
	sf_close(sndfileWRITE);
	fftw_free(samples);
	fftw_free(y);
	fftw_free(freqresponse);
	fftw_free(lpf);
	fftw_free(n);
	fftw_free(w);
	fftw_free(sincOutput);
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
