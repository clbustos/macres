//
// Created by cdx on 10/1/21.
//

#include "utils.h"
#include "old_fft.h"
void F0ToFile(double* f0, int num_frames, const char *filename)
{
    FILE *file;
    int i;

    file = fopen(filename, "w");
    fprintf(file, "time, value\n");
    for(i = 0; i < num_frames; i++)
    {
        fprintf(file,"%d,%f\n",i ,f0[i]);
    }
    fclose(file);
}

void speqToFile(old_fft_complex * spec, int fftl)
{
    FILE *file;
    int i;

    file = fopen("/tmp/speqlist.txt", "w");
    for(i = 0; i < fftl/2+1; i++)
    {
        fprintf(file,"%d,%f\n",i ,spec[i][0]);
    }
    fclose(file);
}

void createFinalPitch(double *f0, int num_frames, double *pitchBend, int bLen, int num_samples, int offset_ms, int sample_rate, double tempo)
{
    int i;
    double *time1, *time2, *pitch;
    int pStep;
    pStep = (int)(60.0 / 96.0 / tempo * sample_rate + 0.5);
    time1 = (double *)malloc(sizeof(double) * num_frames);
    time2 = (double *)malloc(sizeof(double) * bLen);
    pitch = (double *)malloc(sizeof(double) * num_frames);

    for(i = 0;i < num_frames;i++) time1[i] = (double)i * FRAMEPERIOD;
    for(i = 0;i < bLen;i++) time2[i] = (double)i * pStep / (double)sample_rate * 1000.0 + offset_ms/1000.0;
    time2[0] = 0;
    old_interp1(time2, pitchBend, bLen, time1, num_frames, pitch);

    for(i = (int)(offset_ms*FRAMEPERIOD/1000);i < num_frames;i++) f0[i] *= pitch[i];

    //for(i = 0;i < num_frames;i+=10)
    //{
    //	printf("%f\n", pitch[i]);
    //}

    free(time1); free(time2); free(pitch);
}

int base64decoderForUtau(char x, char y)
{
    int ans1, ans2, ans;

    if(x=='+') ans1 = 62;
    if(x=='/') ans1 = 63;
    if(x>='0' && x <= '9') ans1 = x+4;
    if(x>='A' && x <= 'Z') ans1 = x-65;
    if(x>='a' && x <= 'z') ans1 = x-71;

    if(y=='+') ans2 = 62;
    if(y=='/') ans2 = 63;
    if(y>='0' && y <= '9') ans2 = y+4;
    if(y>='A' && y <= 'Z') ans2 = y-65;
    if(y>='a' && y <= 'z') ans2 = y-71;

    ans = (ans1<<6) | ans2;
    if(ans >= 2048) ans -= 4096;
    return ans;
}

int getF0Contour(char *input, double *output)
{
    int i, j, count, length;
    i = 0;
    count = 0;
    double tmp;

    tmp = 0.0;
    while(input[i] != '\0')
    {
        if(input[i] == '#')
        {
            length = 0;
            for(j = i+1;input[j]!='#';j++)
            {
                length = length*10 + input[j]-'0';
            }
            i = j+1;
            for(j = 0;j < length;j++)
            {
                output[count++] = tmp;
            }
        }
        else
        {
            tmp = pow(2.0, (double)base64decoderForUtau(input[i], input[i+1]) / 1200.0);
            output[count++] = tmp;
            i+=2;
        }
    }

    return count;
}

/**
 * Copied from Ameya's world4utau.cpp
 */
double getFreqAvg(double f0[], int num_frames)
{
    int i, j;
    double value = 0, r;
    double p[6], q;
    double freq_avg = 0;
    double base_value = 0;
    for (i = 0; i < num_frames; i++)
    {
        value = f0[i];
        if (value < 1000.0 && value > 55.0)
        {
            r = 1.0;
            // Continuously weight nearby variables more heavily.=
            for (j = 0; j <= 5; j++)
            {
                if (i > j) {
                    q = f0[i - j - 1] - value;
                    p[j] = value / (value + q * q);
                } else {
                    p[j] = 1/(1 + value);
                }
                r *= p[j];
            }
            freq_avg += value * r;
            base_value += r;
        }
    }
    if (base_value > 0) freq_avg /= base_value;
    return freq_avg;
}

/**
 * Base 64 representation of numbers
 * Copied from Ameya's world4utau.cpp
 * Is very weird: 'A'=0, 'Z'=25 an so on, until '/'=63
 * Order: A-Z a-z 0-9  + /
 */
int get64(int c)
{
    if (c >= '0' && c <='9')
    {
        return c - '0' + 52;
    }
    else if (c >= 'A' && c <='Z')
    {
        return c - 'A';
    }
    else if (c >= 'a' && c <='z')
    {
        return c - 'a' + 26;
    }
    else if (c == '+')
    {
        return 62;
    }
    else if (c == '/')
    {
        return 63;
    }
    else
    {
        return 0;
    }
}

/**
 * Decipher pitch.  Copied from Ameya's world4utau.cpp.
 */
int decipherPitch(char *pitch_input, int *destination, int dest_size)
{
    int input_size = 0;
    int i, j;
    int pitch_unchanged_num;
    int pitch_value = 0;
    int cur_dest_index = 0;
    if (pitch_input != NULL)
    {
        input_size = strlen(pitch_input);
        for (i = 0; i < input_size; i += 2)
        {
            // Decode "#somenumber#" for repeated pitch values.
            if (pitch_input[i] == '#')
            {
                i++;
                sscanf(pitch_input + i, "%d", &pitch_unchanged_num);
                for (j = 0; j < pitch_unchanged_num && cur_dest_index < dest_size; j++) {
                    destination[cur_dest_index] = pitch_value;
                    cur_dest_index++;
                }
                while (pitch_input[i] != '#' && pitch_input[i] != 0) i++;
                i--;
            }
                // Decode a 2-char pitch value between -2048 and 2047.
            else
            {
                pitch_value = get64(pitch_input[i]) * 64 + get64(pitch_input[i + 1]);
                if (pitch_value > 2047) pitch_value -= 4096;
                if (cur_dest_index < dest_size) {
                    destination[cur_dest_index] = pitch_value;
                    cur_dest_index++;
                }
            }
        }
    }
    fprintf(stderr, "Number of pitch steps found: %d\n", cur_dest_index);
    for (i = 0; i < cur_dest_index; i++) {
        fprintf(stderr, "%d ", destination[i]);
    }
    fprintf(stderr, "\n");
    return input_size;
}


/**
 * @param f0 pointer to array of pitchs
 * @param num_frames number of frames to analyze
 * @param scaleParam pitch, in format [A-G][ #b][0-8]
 * @param modulationParam Percent of modulation resect to average pitch
 * @param flag_t transpose in semitones/120
 */
void equalizingPitch(double *f0, int num_frames, char *scaleParam, int modulationParam, int flag_t)
{
    int i;
    // First, find the average value.
    double averageF0;
    double modulation;

    modulation = (double)modulationParam / 100.0;

    averageF0 = getFreqAvg(f0, num_frames);

    int scale;
    int octave;
    double targetF0;
    int bias = 0;


    // check if pitch name is correct
    if(scaleParam==NULL) {
        fprintf(stderr, "No pitch set\n");
        exit(-1);
    }
    if(strlen(scaleParam)<2 || strlen(scaleParam)>3) {
        fprintf(stderr, "Incorrect scaleParam definition [AG][ #b][0-8]\n");
        exit(-1);
    }

    // Identify the desired pitch.
    if(scaleParam[1] == '#')
    {
        bias = 1;
    } else if(scaleParam[1] == 'b'){
        bias = -1;
    }

    // You could have E# in notation. Really.
    switch(scaleParam[0])
    {
        case 'C':
            scale = -9+bias;
            break;
        case 'D':
            scale = -7+bias;
            break;
        case 'E':
            scale = -5+bias;
            break;
        case 'F':
            scale = -4+bias;
            break;
        case 'G':
            scale = -2+bias;
            break;
        case 'A':
            scale = bias;
            break;
        case 'B':
            scale = 2+bias;
            break;
    }
    octave = scaleParam[1+abs(bias)]-'0' - 4;
    targetF0 = 440 * pow(2.0,(double)octave) * pow(2.0, (double)scale/12.0);
    targetF0 *= pow(2, (double)flag_t/120);

    double tmp;

    if(averageF0 != 0.0)
    {
        for(i = 0;i < num_frames;i++)
        {
            if(f0[i] != 0.0)
            {
                tmp = ((f0[i]-averageF0) * modulation) + averageF0;
                f0[i] = tmp * targetF0 / averageF0;
            }
        }
    }
    else
    {
        for(i = 0;i < num_frames;i++)
        {
            if(f0[i] != 0.0)
            {
                f0[i] = targetF0;
            }
        }
    }
}

int stretchTime(double *f0, int num_frames, int fftl, int *residualSpecgramIndex,
                double *f02, int num_frames2, int *residualSpecgramIndex2, int os, int st, int ed, int Length2, double vRatio, int mode)
{
    int i, k;
    int st2, ed2;

    st2 = min(num_frames2, (int)((st-os) * vRatio + 0.5));  // Frame of the consonant section after expanding and contracting.
    ed2 = min(num_frames2, (int)(Length2 + 0.5));     // Number of samples after synthesis.
    // First half
    for(i = 0;i < st2;i++)
    {
        k = max(0, min(num_frames-1, int(i/vRatio) + os));
        f02[i] = f0[k];
        residualSpecgramIndex2[i] = residualSpecgramIndex[k];
    }
    // Second half (Extension of the loop formula).
    if(mode == 0)
    {
        i = st2;
        while(i < ed2)
        {
            bool i_updated = false;
            for(k = st; k < ed - 2; k++)
            {
                if(i > ed2-1) break;
                f02[i] = f0[k];
                residualSpecgramIndex2[i] = residualSpecgramIndex[k];
                i++;
                i_updated = true;
            }
            for(k = ed -1; k > st; k--)
            {
                if(i > ed2-1) break;
                f02[i] = f0[k];
                residualSpecgramIndex2[i] = residualSpecgramIndex[k];
                i++;
                i_updated = true;
            }
            if (!i_updated) {
                fprintf(stderr, "Would loop forever while stretching sample.\n");
                return -1;
            }
        }
    }
    else
    {
        // Second half (Extension of the UTAU formula).
        if(ed2-st2 > ed-st) // Extending.
        {
            double ratio;
            ratio = (double)(ed-st)/(ed2-st2);
            for(i = st2;i < ed2; i++)
            {
                k = max(0, min(num_frames-1, (int)((i - st2) * ratio + 0.5 + st)));
                f02[i] = f0[k];
                residualSpecgramIndex2[i] = residualSpecgramIndex[k];
            }
        }
        else
        {
            for(i = st2;i < ed2; i++)
            {
                k = st + (i - st2);
                f02[i] = f0[k];
                residualSpecgramIndex2[i] = residualSpecgramIndex[k];
            }
        }
    }

    return ed2;
}

// F0 low pass filter
void f0Lpf(double *f0, int num_frames, int flag_d)
{
    int i;
    int addcount = 0;
    double addvalue = 0;
    double* newf0;
    newf0 = (double*)malloc(sizeof(double) * num_frames);
    for(i = 0; i < min(num_frames-1, flag_d); i++)
    {
        if(f0[i] != 0.0)
        {
            addvalue += f0[i];
            addcount += 1;
        }
    }
    for(i = 0; i < num_frames; i++)
    {
        if(i - flag_d -1 >= 0)
        {
            if(f0[i - flag_d -1] != 0.0)
            {
                addvalue -= f0[i - flag_d -1];
                addcount -= 1;
            }
        }
        if(i + flag_d <= num_frames - 1)
        {
            if(f0[i + flag_d] != 0.0)
            {
                addvalue += f0[i + flag_d];
                addcount += 1;
            }
        }
        if(f0[i] != 0)
        {
            newf0[i] = addvalue / addcount;
        }
        else
        {
            newf0[i] = 0.0;
        }
    }
    for(i = 0; i < num_frames; i++) f0[i] = newf0[i];
}





/**
 * Apply the 'B' flag (breath)
 */
void breath2(double *f0, int num_frames, int sample_rate, double *waveform, int xLen, old_fft_complex **waveSpecgram,int equLen, int fftl, int flag_B)
{
    int i, j;

    // Prepare the noise FFT.
    double *noiseData;
    double *noiseBuff;
    double *noise;
    old_fft_plan			noise_f_fft;				// fft set
    old_fft_plan			noise_i_fft;				// fft set
    old_fft_complex		*noiseSpec;	// Spectrogram

    noiseData = (double *)malloc(sizeof(double) * xLen);
    for(i=0;i < xLen; i++) noiseData[i] = (double)rand()/RAND_MAX - 0.5;
    noise = (double *)malloc(sizeof(double) * xLen);
    for(i=0;i < xLen; i++) noise[i] = 0.0;
    //for(i=0;i < xLen; i++) noiseData[i] *= noiseData[i] * (noiseData[i] < 0)? -1 : 1;//Play around with noise distribution
    noiseBuff = (double *)malloc(sizeof(double) * fftl);
    noiseSpec = (old_fft_complex *)malloc(sizeof(old_fft_complex) * fftl);
    noise_f_fft = old_fft_plan_dft_r2c_1d(fftl, noiseBuff, noiseSpec, FFT_ESTIMATE);
    noise_i_fft = old_fft_plan_dft_c2r_1d(fftl, noiseSpec, noiseBuff, FFT_ESTIMATE);

    // Prepare the wave FFT.
    old_fft_complex		*waveSpec;	// Spectrogram
    waveSpec = (old_fft_complex *)malloc(sizeof(old_fft_complex) * fftl);

    int offset_ms;
    double volume;

    int SFreq, MFreq, EFreq;

    SFreq = (int)(fftl * 1500 / sample_rate);//Breath beginning frequency
    MFreq = (int)(fftl * 5000 / sample_rate);//Breath beginning frequency
    EFreq = (int)(fftl * 20000 / sample_rate);//Breath frequency band

    double nowIndex;
    int sIndex, eIndex;
    double nowF0;
    int specs, spece;
    double hs, he;
    int baion;

    for(i = 0; i < equLen; i++)
    {
        offset_ms = i * fftl / 2;
        // Copy the data.
        for(j = 0;j < fftl; j++) noiseBuff[j] = noiseData[offset_ms + j] *
                                                (0.5 - 0.5*cos(2.0*PI*(double)j/(double)fftl));// Multiply the window.

        // Apply fft.
        old_fft_execute(noise_f_fft);

        //Spectrogram wraparound（super slapdash）
        for(j = 0;j < fftl/2+1; j++) waveSpec[j][0] = sqrt(waveSpecgram[i][j][0] * waveSpecgram[i][j][0] + waveSpecgram[i][j][1] * waveSpecgram[i][j][1]);
        for(j = 0;j < fftl/2+1; j++) waveSpec[j][0] = log10(waveSpec[j][0]+0.00000001); // logarithmic change
        for(j = 0;j < fftl/2+1; j++) waveSpec[j][1] = waveSpec[j][0];

        nowIndex = max(0.0, min((double)num_frames-1, (double)(offset_ms + fftl / 2) / sample_rate * 1000 / FRAMEPERIOD));
        sIndex = min(num_frames -2, (int)nowIndex);
        eIndex = sIndex + 1;

        nowF0 = (f0[sIndex] == 0 && f0[eIndex] == 0) ?  DEFAULT_F0 :
                (f0[sIndex] == 0) ? f0[eIndex] :
                (f0[eIndex] == 0) ? f0[sIndex] :
                (f0[eIndex] - f0[sIndex]) * (nowIndex - sIndex) + f0[sIndex];

        specs = 0;
        hs = 0.0;
        j = 0;
        baion = 1;
        spece = 0;
        for(baion = 1;spece != fftl/2+1;baion++)
        {
            spece = min(fftl/2+1, (int)((double)fftl / sample_rate * nowF0 * baion + 0.5));
            he = waveSpec[spece][1];
            for(j = specs;j < spece;j++)
            {
                waveSpec[j][0] = (he-hs)/(spece-specs)*(j-specs)+hs;
            }
            specs = spece;
            hs = he;
        }

        for(j = 0;j < fftl/2+1; j++) waveSpec[j][0] = pow(10, waveSpec[j][0]);//振幅化

        // Transform the noise spectrogram.
        for(j = 0;j < SFreq; j++)
        {
            noiseSpec[j][0] = 0.0;
            noiseSpec[j][1] = 0.0;
        }

        for(;j < MFreq; j++)
        {
            volume = waveSpec[j][0] * (0.5 - 0.5 * cos(PI * (j - SFreq) / (double)(MFreq - SFreq)));
            noiseSpec[j][0] *= volume;
            noiseSpec[j][1] *= volume;
        }
        for(;j < EFreq; j++)
        {
            volume = waveSpec[j][0] * (0.5 - 0.5 * cos(PI + PI * (j - MFreq) / (double)(EFreq - MFreq)));
            noiseSpec[j][0] *= volume;
            noiseSpec[j][1] *= volume;
        }

        for(;j < fftl/2+1; j++)
        {
            noiseSpec[j][0] = 0.0;
            noiseSpec[j][1] = 0.0;
        }

        noiseSpec[0][1] = 0.0;
        noiseSpec[fftl/2][1] = 0.0;

        // Reverse FFT
        old_fft_execute(noise_i_fft);
        for(j = 0;j < fftl; j++) noiseBuff[j] /= fftl;

        // Multiply the window
        //for(j = 0;j < fftl; j++) noiseBuff[j] *= 0.5 - 0.5*cos(2.0*PI*(double)j/(double)fftl);

        // Add on the noise.
        for(j = 0;j < fftl; j++)
        {
            noise[offset_ms + j] += noiseBuff[j] * 0.2;
        }
    }

    // Noise synthesis.
    double noiseRatio = max(0.0, (double)(flag_B - 50) / 50.0);
    double waveRatio = 1 - noiseRatio;
    for(i = 0;i < xLen;i++) waveform[i] = waveform[i] * waveRatio + noise[i] * noiseRatio;

    // Clean up.
    old_fft_destroy_plan(noise_f_fft);
    old_fft_destroy_plan(noise_i_fft);
    free(noise);
    free(noiseData);
    free(noiseBuff);
    free(noiseSpec);
    free(waveSpec);
}

/**
 * Apply the 'O' flag (voice strength)
 */
void Opening(double *f0, int num_frames, int sample_rate, old_fft_complex **waveSpecgram,int equLen, int fftl, int flag_O)
{
    int i, j;
    double opn = (double) flag_O / 100.0;
    int sFreq = (int)(fftl * 500 / sample_rate); // Control frequency 1
    int eFreq = (int)(fftl * 2000 / sample_rate); // Control frequency 2
    double sRatio = -10.0; // Control frequency 1's amplitude factor in decibels.
    double eRatio = 10.0; // Control frequency 2's amplitude factor in decibels.

    // Make a volume map for each frequency.
    double volume;
    double *volumeMap;
    volumeMap = (double *)malloc(sizeof(double) * fftl/2+1);

    volume = pow(10, sRatio * opn / 20);
    for(j = 0;j < sFreq;j++)
    {
        volumeMap[j] = volume;
    }
    for(;j < eFreq;j++)
    {
        volume = pow(10, ((0.5+0.5*cos(PI+PI/(eFreq-sFreq)*(j-sFreq)))*(eRatio-sRatio)+sRatio) * opn / 20);
        volumeMap[j] = volume;
    }
    volume = pow(10, eRatio * opn / 20);
    for(;j < fftl/2+1;j++)
    {
        volumeMap[j] = volume;
    }

    // Change the volume for each frequency.
    int f0Frame;
    for(i = 0;i < equLen;i++)
    {
        f0Frame = max(0, min(num_frames-1, (int)((double)((i+1) * fftl / 2) / sample_rate * 1000 / FRAMEPERIOD + 0.5)));
        if(f0[f0Frame] == 0.0) continue;
        for(j = 0;j < fftl/2+1;j++)
        {
            waveSpecgram[i][j][0] *= volumeMap[j];
            waveSpecgram[i][j][1] *= volumeMap[j];
        }
    }

    free(volumeMap);
}

/**
 * Apply the 'b' flag（emphasizing the unvoiced consonant)
 */
void consonantAmp2(double *f0, double *volume, int num_frames, int flag_b)
{
    int i;
    int frameLen = 5; // Number of smoothing frames (front and back).
    int addCount = 0;
    double addVolume = 0;
    double ratio = (double) flag_b / 20.0; // Scaling factor (5 when b=100)

    for(i = 0;i < min(num_frames, frameLen+1); i++)
    {
        addCount++;
        addVolume += (f0[i] == 0) ? ratio : 0.0;
    }
    for(i = 0;i < num_frames-1; i++)
    {
        volume[i] *= (addCount != 0) ? addVolume / addCount + 1.0 : 1.0;

        if(i >= frameLen)
        {
            addCount--;
            addVolume -= (f0[i-frameLen] == 0) ? ratio : 0.0;
        }
        if(i <= num_frames-1-frameLen-1)
        {
            addCount++;
            addVolume += (f0[i+frameLen+1] == 0) ? ratio : 0.0;
        }
    }
}

/**
 * Apply the 'g' flag (gender where < 0 is feminine, > 0 is masculine)
 */
void gFactor(int pCount, int fftl, double **residualSpecgram, int *residualSpecgramLength, double gRatio)
{
    int i, j;
    double position;
    int sindex, eindex;
    int NewLength;

    for(i = 0; i < pCount-1; i++)
    {
        if(residualSpecgramLength[i] == 0.0) continue;

        NewLength = max(0, min(fftl-1, (int)(residualSpecgramLength[i] / gRatio + 0.5)));
        if (gRatio>1)
        {
            for(j = 0;j < NewLength;j++)
            {
                position = min((double)residualSpecgramLength[i]-1.0001, (double)(j * gRatio));
                sindex = (int)position;
                eindex = sindex + 1;
                residualSpecgram[i][j] = residualSpecgram[i][sindex] +
                                         (double)(residualSpecgram[i][eindex] - residualSpecgram[i][sindex]) *
                                         (double)(position - sindex);
            }
        }
        else
        {
            for(j = NewLength-1;j >= 0;j--)
            {
                position = min((double)residualSpecgramLength[i]-1.0001, (double)(j * gRatio));
                sindex = (int)position;
                eindex = sindex + 1;
                residualSpecgram[i][j] = residualSpecgram[i][sindex] +
                                         (double)(residualSpecgram[i][eindex] - residualSpecgram[i][sindex]) *
                                         (double)(position - sindex);
            }
        }
        residualSpecgramLength[i] = NewLength;
    }
}

// Correct the cycle of noise segments that have changed b/c of the g flag.（Match the F0 in cycles.）
// So that the noise part of F0 can become nonzero, apply this just before synthesisPt101.
void f0FixG(double *f0, int num_frames2, double gRatio)
{
    int i;
    for(i = 0;i < num_frames2;i++)
    {
        if(f0[i] == 0.0)
        {
            f0[i] = DEFAULT_F0 * gRatio;
        }
    }
}

// Add noise to the f0 sequence.
void f0Noise(double *f0, int num_frames, double f0Rand)
{
    int i, j;
    int Pit = 1;//(int)(5 / FRAMEPERIOD + 0.5); // Number of frames in pitch change interval (50ms?)
    double sRand, eRand;
    double NowRand;

    eRand = 0;
    i = 0;
    while(i <= num_frames-1)
    {
        sRand = eRand;
        //eRand = (double)rand()/(RAND_MAX+1) * f0Rand * 2  - f0Rand;
        //eRand = (double)rand()/(RAND_MAX+1) * -f0Rand;
        eRand = (double)(int)(rand()/(RAND_MAX/3)-1) * f0Rand;
        for(j = 0;(j < Pit) && (i+j <= num_frames-1); j++)
        {
            if(f0[i+j] != 0.0)
            {
                NowRand = (eRand - sRand) / Pit * j + sRand;
                f0[i+j] *= pow(2, NowRand);
            }
        }
        i += j;
    }
}

// Convert frequency into pitch.
double FrqToPit(double Frq)
{
    return log(Frq / 220) * 1.44269504088896 * 1200 + 5700;
}

// A flag（correct the volume of all combined pitch changes）.
void autoVolume(double *f0, int num_frames, int sample_rate, double *volume, int flag_A)
{
    int i;

    if(flag_A == 0)
    {
        for(i = 0;i < num_frames; i++) volume[i] = 1.0;
        return;
    }

    double AutoPow;
    for(i = 0;i < num_frames-1; i++)
    {
        if(f0[i] == 0.0)
        {
            volume[i] = 1.0;
            continue;
        }

        if (f0[i+1] != 0.0)
        {
            AutoPow = (FrqToPit(f0[i+1]) - FrqToPit(f0[i])) * (441 / (sample_rate * FRAMEPERIOD)) * flag_A;
            volume[i] = min(1.2, pow(2, AutoPow * 1));

            continue;
        }

        if(i > 0)
        {
            if(f0[i-1] != 0.0)
            {
                volume[i] = volume[i-1];
                continue;
            }
        }
        volume[i] = 1.0;
    }
    if(f0[num_frames-1] != 0.0 && f0[num_frames-2] != 0.0) volume[num_frames-1] = volume[num_frames-2];
}

/**
 * Create wave specgram for equalizing.
 */
void createWaveSpec(double *waveform, int xLen, int fftl, int equLen, old_fft_complex **waveSpecgram)
{
    int i, j;

    double *waveBuff;
    old_fft_plan			wave_f_fft;				// fft set
    old_fft_complex		*waveSpec;					// Specgram
    waveBuff = (double *)malloc(sizeof(double) * fftl);
    waveSpec = (old_fft_complex *)malloc(sizeof(old_fft_complex) * fftl);
    wave_f_fft = old_fft_plan_dft_r2c_1d(fftl, waveBuff, waveSpec, FFT_ESTIMATE);

    int offset_ms;

    for(i = 0;i < equLen;i++)
    {
        offset_ms = i * fftl / 2;
        // Copy the data.
        for (j = 0; j < fftl; j++) {
            waveBuff[j] = waveform[offset_ms + j] *
                          (0.5 - 0.5 * cos(2.0*PI*(double)j/(double)fftl)); // Multiply the window.
        }

        // Apply fft.
        old_fft_execute(wave_f_fft);

        // Load spectrogram into memory.
        for(j = 0;j < fftl/2+1; j++)
        {
            waveSpecgram[i][j][0] = waveSpec[j][0];
            waveSpecgram[i][j][1] = waveSpec[j][1];
        }
    }

    old_fft_destroy_plan(wave_f_fft);
    free(waveBuff);
    free(waveSpec);

}

/**
 * Rebuild wave from spectrogram
 */
void rebuildWave(double *waveform, int xLen, int fftl, int equLen, old_fft_complex **waveSpecgram)
{
    int i, j;
    double *waveBuff;
    old_fft_plan			wave_i_fft;				// fft set
    old_fft_complex		*waveSpec;	// Spectrogram.
    waveBuff = (double *)malloc(sizeof(double) * fftl);
    waveSpec = (old_fft_complex *)malloc(sizeof(old_fft_complex) * fftl);
    wave_i_fft = old_fft_plan_dft_c2r_1d(fftl, waveSpec, waveBuff, FFT_ESTIMATE);

    int offset_ms;
    for(i = 0;i < xLen;i++) waveform[i] = 0;

    for(i = 0;i < equLen;i++)
    {
        offset_ms = i * fftl / 2;

        // Load spectrogram into memory.
        for(j = 0;j < fftl/2+1; j++)
        {
            waveSpec[j][0] = waveSpecgram[i][j][0];
            waveSpec[j][1] = waveSpecgram[i][j][1];
        }


        // Apply fft.
        old_fft_execute(wave_i_fft);

        for(j = 0;j < fftl; j++) waveBuff[j] /= fftl;

        // Copy the data.
        for(j = 0;j < fftl; j++) waveform[offset_ms + j]  += waveBuff[j];

    }

    old_fft_destroy_plan(wave_i_fft);
    free(waveBuff);
    free(waveSpec);

}

void platinumToFile(int pCount, double **residualSpecgram, int *residualSpecgramLength, int *residualSpectralIndex, const char* filename) {
    //printf("pCount: %d\n", pCount);
    FILE *file;

    file = fopen(filename, "w");
    fprintf(file, "i,j, spectral_index, spectral_length, residual\n");

    for(int i=0; i<pCount; i++) {
        //printf("Index: %d -> index: %d, length:%d\n", i, residualSpectralIndex[i], residualSpecgramLength[i]);
        if(residualSpecgramLength[i]>0) {
            for(int j=0; j<residualSpecgramLength[i]; j++) {
                fprintf(file, "%d, %d, %d, %d, %lf\n", i, j, residualSpectralIndex[i], residualSpecgramLength[i],  residualSpecgram[i][j]);
            }
            //printf("\n");
        } else {
            fprintf(file, "%d, %d, %d, %d, %lf\n", i, 0, residualSpectralIndex[i], residualSpecgramLength[i],  0.0);
        }
    }
    fclose(file);
}