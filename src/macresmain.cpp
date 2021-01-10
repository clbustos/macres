// UTAU vocal synthesis engine『MacRes』 version 0.0.1    3/10/2017
// Originally branched from tn_fnds.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <windows.h>

#include "utils.h"
#include "world.h"
#include "world/dio.h"
#include "wavread.h"

#include <math.h>

inline DWORD timeGetTime() { return (DWORD)time(NULL) * 1000; }




// 13 inputs:
// 1 Input file
// 2 Output file
// 3 Pitch (ex: "C#4")
// 4 Consonant velocity
// 5 Flags
// 6 offset_ms from oto
// 7 Desired note length
// 8 Consonant length
// 9 cutoff from oto
// 10 Volume
// 11 Modulation
// 12 Tempo
// 13 Pitchbends



#pragma comment(lib, "winmm.lib")




/**
 * Main method.
 */
int main(int argc, char *argv[])
{
	enum PARAMS {
		INPUTFILE=1,
		OUTPUTFILE=2,
		NOTENUM=3,
		VELOCITY=4, 
		FLAGS=5,
		OFFSET_MS=6, 
		NOTELENGTH=7,
		FIXEDLENGTH=8,
		END=9,
		INTENSITY=10,
		MODULATION=11,
		TEMPO=12,
		PITCHBENDS=13
	};
	int i;

	double *waveform,*f0,*time_axis,*y;
	double **residualSpecgram;
	int *residualSpecgramLength;
	int *residualSpecgramIndex;
	int pCount;
	int fftl;

	int num_samples;
	int num_frames;

	if(argc < 3)
	{
		printf("Params: inputfile outputfile notenum velocity flags offset_ms notelength");
		printf(" fixedlength end intensity modulation tempo pitchbends\n");
		return 0;
	}

/*
	printf("argc:%d\n", argc);
	for(i = 0;i < argc-1;i++)
		printf("%s\n", argv[i]);
*/

	// Parse flags
	char *string_buf;
	int cur_char_index;
	int flag_B = 50; // Breath

	if(argc > 5 && (string_buf = strchr(argv[FLAGS],'B')) != 0)
	{
		cur_char_index = string_buf - argv[FLAGS];
		if ((cur_char_index == 0) || (argv[FLAGS][cur_char_index - 1] != 'M'))
		{
			sscanf(string_buf+1, "%d", &flag_B);
			flag_B = max(0, min(100, flag_B));
		}
	}

	int flag_b = 0; // Consonant strength
	if(argc > 5 && (string_buf = strchr(argv[FLAGS],'b')) != 0)
	{
		cur_char_index = string_buf - argv[FLAGS];
		if ((cur_char_index == 0) || (argv[FLAGS][cur_char_index - 1] != 'M'))
		{
			sscanf(string_buf+1, "%d", &flag_b);
			flag_b = max(0, min(100, flag_b));
		}
	}

	int flag_t = 0; // 't' flag, to transpose X/120 semitones 
	if(argc > 5 && (string_buf = strchr(argv[FLAGS],'t')) != 0)
	{
		cur_char_index = string_buf - argv[FLAGS];
		if ((cur_char_index == 0) || (argv[FLAGS][cur_char_index - 1] != 'M'))
		{
			sscanf(string_buf+1, "%d", &flag_t);
		}
	}

	double flag_g = 0.0; // Gender flag
	double gRatio;
	if(argc > 5 && (string_buf = strchr(argv[FLAGS],'g')) != 0)
	{
		cur_char_index = string_buf - argv[FLAGS];
		if ((cur_char_index == 0) || (argv[FLAGS][cur_char_index - 1] != 'M'))
		{
			sscanf(string_buf+1, "%lf", &flag_g);
			if (flag_g > 100) flag_g = 100;
			if (flag_g < -100) flag_g= -100;
		}
	}
	gRatio = pow(10, -flag_g / 200);

	// W flag（frequency enforcement settings）F<0 unvoiced  F=0 no effect  50>=F<=1000
	// Set up the designated frequency.
	double flag_W = 0.0;
	double f0Rand = 0;
	if(argc > 5 && (string_buf = strchr(argv[FLAGS], 'W')) != 0)
	{
		cur_char_index = string_buf - argv[FLAGS];
		if ((cur_char_index == 0) || (argv[FLAGS][cur_char_index - 1] != 'M'))
		{
			sscanf(string_buf+1, "%lf", &flag_W);
			if (flag_W > 1000) flag_W = 1000;
			if ((flag_W < 50) && (flag_W > 0)){f0Rand =  flag_W / 50; flag_W = 0;}
			if (flag_W < 0) flag_W = -1;
		}
	}

	// Original flag: Hang/multiply? LPF on DIO's F0 analysis result. 0~20 def 5
	// Can only support the default for now.
	int flag_d = 5;
	//if(argc > 5 && (string_buf = strchr(argv[FLAGS],'d')) != 0)
	//{
	//	sscanf(string_buf+1, "%d", &flag_d);
	//	flag_d = max(0, min(20, flag_d));
	//}

	int flag_A = 0; // Original flag: Correct volume of combined pitch changes.
	if(argc > 5 && (string_buf = strchr(argv[FLAGS],'A')) != 0)
	{
		cur_char_index = string_buf - argv[FLAGS];
		if ((cur_char_index == 0) || (argv[FLAGS][cur_char_index - 1] != 'M'))
		{
			sscanf(string_buf+1, "%d", &flag_A);
			flag_A = max(0, min(100, flag_A));
		}
	}

	int flag_O = 0; // Original flag: voice strength
	if(argc > 5 && (string_buf = strchr(argv[FLAGS],'O')) != 0)
	{
		cur_char_index = string_buf - argv[FLAGS];
		if ((cur_char_index == 0) || (argv[FLAGS][cur_char_index - 1] != 'M'))
		{
			sscanf(string_buf+1, "%d", &flag_O);
			flag_O = max(-100, min(100, flag_O));
		}
	}

	// Original flag: Change the vowel stretching method: default is loops.
	// However you can choose to use UTAU's default method.
	int flag_e = 0;
	if(argc > 5 && (string_buf = strchr(argv[FLAGS],'e')) != 0)
	{
		cur_char_index = string_buf - argv[FLAGS];
		if ((cur_char_index == 0) || (argv[FLAGS][cur_char_index - 1] != 'M'))
		{
			flag_e = 1;
		}
	}

	FILE *file;

	int offset_ms; // This is also the starting point, or stp.
	int cutoff_ms;
	offset_ms = atoi(argv[OFFSET_MS]);
	cutoff_ms = atoi(argv[END]);

	int sample_rate, bits_per_sample;
	waveform = ReadWaveFile(argv[INPUTFILE], &sample_rate, &bits_per_sample, &num_samples, &offset_ms, &cutoff_ms);

	if(waveform == NULL)
	{
		fprintf(stderr, "Error: Your input file does not exist.\n");
		return 0;
	}

	printf("File information\n");
	printf("Sampling : %d Hz %d Bit\n", sample_rate, bits_per_sample);
	printf("Length %d [sample]\n", num_samples);
	printf("Length %f [sec]\n", (double)num_samples/(double)sample_rate);

	// Calculate beforehand the number of samples in F0 (one per FRAMEPERIOD ms).
	num_frames = GetSamplesForDIO(sample_rate, num_samples, FRAMEPERIOD);
	f0 = (double *)malloc(sizeof(double)*num_frames);
	time_axis  = (double *)malloc(sizeof(double)*num_frames);

	// Start to estimate F0 contour (fundamental frequency) using DIO.
	DWORD elapsedTime;
	if(flag_W == 0) // F flag: F0 enforcement settings.
	{
		printf("\nAnalysis\n");
		elapsedTime = timeGetTime();
/**
		old_dio(waveform, num_samples, sample_rate, FRAMEPERIOD, time_axis, f0);

		printf("DIO: %d [msec]\n", timeGetTime() - elapsedTime);

        F0ToFile(f0, num_frames, "test_1.dat");
**/

        DioOption dio_option;
        InitializeDioOption(&dio_option);
        dio_option.frame_period=FRAMEPERIOD;
        dio_option.f0_floor = FLOOR_F0;
        dio_option.f0_ceil = 640;
        dio_option.channels_in_octave = 2.0;
        dio_option.speed= 1;
        Dio(waveform, num_samples, sample_rate, &dio_option, time_axis, f0);


		//F0ToFile(f0, num_frames, "test_2.dat");
		//F0's Low Pass Filter
		if (flag_d !=0)
		{
			f0Lpf(f0, num_frames, flag_d);
		}
		//F0ToFile(f0, num_frames);

	}
	else
	{
		for(i = 0;i < num_frames;i++)
		{
			f0[i] = (flag_W == -1) ? 0.0 : flag_W;
			time_axis[i] = (double)i * FRAMEPERIOD/1000.0;
		}
	}

	fftl = getFFTLengthForStar(sample_rate);

	// Acyclic indicator analysis, using PLATINUM.
	elapsedTime = timeGetTime();
	residualSpecgramIndex = (int *)malloc(sizeof(int) * num_frames);

	pCount = pt101(waveform, num_samples, sample_rate, time_axis, f0, &residualSpecgram, &residualSpecgramLength, residualSpecgramIndex);
	printf("PLATINUM: %d [msec]\n", timeGetTime() - elapsedTime);
    printf("pCount: %d num_frames: %d \n",pCount, num_frames);

   // platinumToFile(pCount, residualSpecgram, residualSpecgramLength, residualSpecgramIndex, "residual.dat");

	// Apply gender flag if necessary.
	if(flag_g != 0)
	{
		 gFactor(pCount, fftl, residualSpecgram, residualSpecgramLength, gRatio);
	}

	// Multiply the window.
	PulseResidualWindow(residualSpecgram, residualSpecgramLength, pCount);

	// Expand and contract the time.
	int lengthMsec, snum_framesgthMsec, inpunum_framesgthMsec;
	double velocity;
	double vRatio;

	inpunum_framesgthMsec = (int)(num_frames*FRAMEPERIOD);// Length of input noise available.
	lengthMsec = atoi(argv[NOTELENGTH]);               // Desired note length
	snum_framesgthMsec = atoi(argv[FIXEDLENGTH]);             // Length of consonant section.
	velocity = (double)atoi(argv[VELOCITY]);         // Consonant velocity.
	vRatio = pow(2.0, (1.0 - (velocity / 100.0))); // Consonant expansion/contraction ratio.

	// Guarantee memory for the control parameters.
	double *fixedF0;
	int *fixedResidualSpecgramIndex;
	double *fixedVolume;         // Volume of frame unit.

	int num_frames2;

    num_frames2 = (int)(0.5+(double)(lengthMsec)/FRAMEPERIOD);

	fixedF0					= (double *) malloc(sizeof(double)   * num_frames2);
	fixedResidualSpecgramIndex	= (int *) malloc(sizeof(int) * num_frames2);
	fixedVolume	= (double *) malloc(sizeof(double) * num_frames2);

	// Guarantee memory for the final waveform.
	int num_samples2;
	num_samples2 = (int)((lengthMsec)/1000.0*(double)sample_rate);
	y  = (double *)malloc(sizeof(double)*num_samples2);
	for(i = 0;i < num_samples2;i++) y[i] = 0.0;
    //printf("length:%d, %f\n",num_samples2, (double)num_samples2/(double)sample_rate*1000.0);
    //printf("%d, %d, %d\n",lengthMsec, offset_ms, sample_rate);

	// Fiddle with F0 before synthesis according to input params.
	equalizingPitch(f0, num_frames, argv[NOTENUM], atoi(argv[MODULATION]), flag_t);

	// Growing and shrinking the time.
	int os, st, ed;
	os = offset_ms;
	st = snum_framesgthMsec + offset_ms;
	ed = inpunum_framesgthMsec - cutoff_ms;

    // number of frames depends on DIO settings, 
	// set by FRAMEPERIOD constant
	num_frames2 = stretchTime(f0, num_frames, fftl, residualSpecgramIndex,
			fixedF0, num_frames2, fixedResidualSpecgramIndex,
			os/(int)FRAMEPERIOD, st/(int)FRAMEPERIOD, min(ed/(int)FRAMEPERIOD, num_frames-1),
			lengthMsec, vRatio, flag_e);
	if (num_frames2 == -1) {
		fprintf(stderr, "Error while stretching sample.\n");
		return 0;
	}
	printf("Num frames (pitch detection by DIO)=%d\n", num_frames2);

	// Let the world4utau library handle the pitchbends.
	int *pitch = NULL;
	double tempo = 120;
	int pLen = num_frames2;
	int pStep = 256;
	if (argc > 13)
	{
		string_buf = argv[TEMPO];
		//printf("Tempo original: %s\n", string_buf);
		sscanf(string_buf + 1, "%lf", &tempo);
		// 96 pitch steps in a beat.
		pStep = (int)(60.0 / 96.0 / tempo * sample_rate + 0.5);
		pLen = num_samples2 / pStep + 1;
		printf("Tempo: %0.3f. Step number: %d. Step length: %d\n", tempo, pStep, pLen);
		pitch = (int*)malloc((pLen+1) * sizeof(int));
		memset(pitch, 0, (pLen+1) * sizeof(int));
		decipherPitch(argv[PITCHBENDS], pitch, pLen);
	}
	else
	{
		pitch = (int*)malloc((pLen+1) * sizeof(int));
		memset(pitch, 0, (pLen+1) * sizeof(int));
	}

	double cur_millis; // Current time in milliseconds.
	double amt_into_cur_step;
	int cur_step;
	// For every frame in the resized F0, adjust pitch.
	for (i = 0; i < num_frames2; i++)
  	{
		cur_millis = FRAMEPERIOD * i;
		amt_into_cur_step = cur_millis * 0.001 * sample_rate / pStep;
		cur_step = (int)floor(amt_into_cur_step);
		amt_into_cur_step -= cur_step;
		if (cur_step >= pLen) cur_step = pLen - 1;
		fixedF0[i] *= pow(2, (pitch[cur_step] * (1.0 - amt_into_cur_step) +
				pitch[cur_step + 1] * amt_into_cur_step) / 1200.0);
	}
	//createFinalPitch(fixedF0, num_frames2, pitchBend, bLen, num_samples2, offset_ms, sample_rate, tempo);

	// W flag's pitch noise: can't envision the death voice change very well. (???)
	//if(f0Rand != 0.0)
	//{
	//f0Noise(fixedF0, num_frames2, f0Rand);

	//}

	// Apply the 'A' flag.
	autoVolume(fixedF0, num_frames2, sample_rate, fixedVolume, flag_A);

	// Apply the 'b' flag if necessary.
	if(flag_b != 0)
	{
		consonantAmp2(fixedF0, fixedVolume, num_frames2, flag_b);
	}

	// Adjust the unvoiced cycles that have been put out of alignment because of the g flag.
	double fixedDefault_f0 = DEFAULT_F0 * gRatio;

	// Synthesis step.
	printf("\nSynthesis\n");
	elapsedTime = timeGetTime();
	synthesisPt101(fixedDefault_f0, fixedF0, num_frames2, residualSpecgram, residualSpecgramLength, fixedResidualSpecgramIndex,
		fixedVolume, fftl, FRAMEPERIOD, sample_rate, y, num_samples2);

	printf("WORLD: %d [msec]\n", timeGetTime() - elapsedTime);

	// Equalizing.
	int equfftL = 1024; // Equalizer's fft length
	int equLen = (num_samples2 / (equfftL/2)) - 1; //繰り返し回数
	old_fft_complex **waveSpecgram;  // spectrogram
	waveSpecgram = (old_fft_complex **)malloc(sizeof(old_fft_complex *) * equLen);
	for(i = 0;i < equLen;i++) waveSpecgram[i] = (old_fft_complex *)malloc(sizeof(old_fft_complex) * (equfftL/2+1));

	// Create wave specgram from y.
	if(flag_B > 50 || flag_O != 0)
	{
		createWaveSpec(y, num_samples2, equfftL, equLen, waveSpecgram);
	}

	// Apply 'O' flag if necessary.
	if(flag_O != 0)
	{
		Opening(fixedF0, num_frames2, sample_rate, waveSpecgram, equLen, equfftL, flag_O);
	}

	// Put the equalizer results (wave_specgram) back into a waveform (y).
	if(flag_O != 0)
	{
		rebuildWave(y, num_samples2, equfftL, equLen, waveSpecgram);
	}

	// Apply noise if B flag over 50.
	if(flag_B > 50)
	{
		 breath2(fixedF0, num_frames2, sample_rate, y, num_samples2, waveSpecgram, equLen, equfftL, flag_B);
	}

	// Offset setup.
	//num_samples2 = (int)((lengthMsec)/1000.0*(double)sample_rate);

	// Write output to file.
	char header[44];
	short *output;
	double maxAmp;
	output = (short *)malloc(sizeof(short) * num_samples2);

	// Amplitude normalization.
	maxAmp = 0.0;
	double volume;
	volume = (double)atoi(argv[INTENSITY]) / 100.0;
	for(i = 0;i < num_samples2;i++) maxAmp = maxAmp < fabs(y[i]) ? fabs(y[i]) : maxAmp;
	for(i = 0;i < num_samples2;i++) output[i] = (short)(32768.0*(y[i]*0.5 * volume/maxAmp));

	file = fopen(argv[1], "rb");
	size_t result = fread(header, sizeof(char), 22, file);
	assert(result == 22);
	fclose(file);

	*((short int*)(&header[22])) = 1;							// channels	 	 2
	*((int*)(&header[24])) = sample_rate;						// samplerate 	 4
	*((int*)(&header[28])) = sample_rate * bits_per_sample / 8;	// bytepersec 	 4
	*((short int*)(&header[32])) = bits_per_sample / 8;			// bytepersample 2
	*((short int*)(&header[34])) = bits_per_sample;				// bitspersample 2

	header[36] = 'd'; header[37] = 'a'; header[38] = 't'; header[39] = 'a';

	file = fopen(argv[2],"wb");
	fwrite(header, sizeof(char), 44, file);
	fwrite(output, sizeof(short), num_samples2, file);
	fseek(file, 40, SEEK_SET);
	num_samples2*=2;
	fwrite(&num_samples2, sizeof(int), 1, file);
	fclose(file);
	free(output);

	free(pitch);
	free(waveform);
	free(time_axis);
	free(f0);
	free(fixedF0);
	free(y);
	for(i = 0;i < pCount;i++)
	{
		free(residualSpecgram[i]);
	}
	free(residualSpecgram);
	free(fixedResidualSpecgramIndex);
	free(fixedVolume);
	free(residualSpecgramIndex);
	free(residualSpecgramLength);

	for(i = 0;i < equLen;i++) free(waveSpecgram[i]);
	free(waveSpecgram);

	fprintf(stderr, "Complete.\n");

	return 0;
}
