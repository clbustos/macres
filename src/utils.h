//
// Created by cdx on 10/1/21.
//

#ifndef MACRES_UTILS_H
#define MACRES_UTILS_H

#include "old_fft.h"
#include "world.h"
void f0Lpf(double *f0, int num_frames, int flag_d);
void gFactor(int pCount, int fftl, double **residualSpecgram, int *residualSpecgramLength, double gRatio);
void equalizingPitch(double *f0, int num_frames, char *scaleParam, int modulationParam, int flag_t);
int stretchTime(double *f0, int num_frames, int fftl, int *residualSpecgramIndex,
                double *f02, int num_frames2, int *residualSpecgramIndex2, int os, int st, int ed, int Length2, double vRatio, int mode);
int decipherPitch(char *pitch_input, int *destination, int dest_size);
void consonantAmp2(double *f0, double *volume, int num_frames, int flag_b);
void autoVolume(double *f0, int num_frames, int sample_rate, double *volume, int flag_A);
void Opening(double *f0, int num_frames, int sample_rate, old_fft_complex **waveSpecgram,int equLen, int fftl, int flag_O);
void breath2(double *f0, int num_frames, int sample_rate, double *waveform, int xLen, old_fft_complex **waveSpecgram,int equLen, int fftl, int flag_B);
void createWaveSpec(double *waveform, int xLen, int fftl, int equLen, old_fft_complex **waveSpecgram);
void rebuildWave(double *waveform, int xLen, int fftl, int equLen, old_fft_complex **waveSpecgram);

void platinumToFile(int pCount, double **residualSpecgram, int *residualSpecgramLength, int *residualSpectralIndex, const char* filename);

#endif //MACRES_UTILS_H
