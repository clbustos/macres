// MacRes v0.0.1   1/7/18

#include "world.h"

#include <stdio.h> // for debug
#include <stdlib.h>
#include <math.h>

// for debugging
void pulseToFile(int pCount, double *pulseLocations, int *residualSpecgramLength)
{
	FILE *file;
	int i;

	file = fopen("/tmp/pulse.txt", "w");
	for(i = 0; i < pCount; i++)
	{
		fprintf(file,"%d,%f,%d\n",i ,pulseLocations[i], residualSpecgramLength[i]);
	}
	fclose(file);
}

void zeroXToFile(double *x, int xLen, double *f0interpolatedRaw, double *totalPhase, int pCount, double *pulseLocations, int fs, int vuvNum, int *wedgeList)
{
	FILE *file;
	int i,j;

	int *zeroX;
	zeroX = (int *)malloc(sizeof(int) * xLen);
	for(i=0;i<xLen;i++)zeroX[i] = 0;
	for(i=0;i<pCount;i++)zeroX[int(pulseLocations[i] * fs)] = 1;
	file = fopen("/tmp/zerocross.txt", "w");
	for(i = 0; i < xLen; i++)
	{
		fprintf(file,"%d,%f,%f,%f,%d",i ,x[i] ,f0interpolatedRaw[i], totalPhase[i], zeroX[i]);
		for(j=0;j<vuvNum;j++)
		{
			if(wedgeList[j] == i){
				fprintf(file,",2"); 
			}
		}
		fprintf(file,"\n");
	}

	free(zeroX);

	fclose(file);
}


void getOneFrameResidualSignal(double *x, int xLen, int fs, int positionIndex, double framePeriod, double f0, int fftl, double *pulseLocations, int pCount,
							double *residualSpec)
{
	int i;
	double T0;
	int index, tmpIndex, wLen;

	double tmp;
	double tmpValue = 100000.0; // safeGuard
	for(i = 0;i < pCount;i++)
	{
		tmp = fabs(pulseLocations[i] - (double)positionIndex*framePeriod);//フレームに最も近いパルス？を探す
		if(tmp < tmpValue)
		{
			tmpValue = tmp;
			tmpIndex = i;
		}
		index = 1+(int)(0.5+pulseLocations[tmpIndex]*fs);//最も近いパルス？のサンプル位置?
	}

	T0 = (double)fs/f0;//１周期のサンプル数（実数）
	wLen = (int)(0.5 + T0*2.0);//2周期のサンプル数

	if(wLen+index-(int)(0.5+T0) >= xLen)
	{
		for(i = 0;i < fftl;i++) residualSpec[i] = 0.0;
		return;
	}

	for(i = 0;i < wLen;i++)
	{
		tmpIndex = i+index - (int)(0.5+T0);//最も近いパルス？の前後1周期分のインデックス
		residualSpec[i] = x[min( xLen-1, max(0, tmpIndex))] * 
		(0.5 - 0.5*cos(2.0*PI*(double)(i+1)/((double)(wLen+1))));//窓を掛ける
	}
	for(;i < fftl/2;i++)
	{
		residualSpec[i] = 0.0;
	}
}
void getOnePulseResidualSignal(double *x, int xLen, int fs, double framePeriod, double *f0, int fftl, double *pulseLocations, int pCount,
							   double **residualSpecgram, int *residualSpecgramLength)
{
	int i, j;
	double T0;
	int index, tmpIndex, wLen;
	double ff0;
	double f0fi;
	int f0si, f0ei;

	residualSpecgramLength[pCount-1] = 0;
	for(j = 0;j < pCount-1; j++)
	{
		//サンプル位置を計算
		index = 1+(int)(0.5+pulseLocations[j]*fs);

		//f0を計算
		f0fi = pulseLocations[j] / (double)framePeriod;
		f0si = (int)f0fi;
		f0ei = f0si + 1;

		ff0 = (f0[f0si] == 0.0 || f0[f0ei] == 0.0)? DEFAULT_F0:
				(f0[f0si] == 0)? f0[f0ei]:            
				(f0[f0ei] == 0)? f0[f0si]:  
								 f0[f0si] + (f0[f0ei] - f0[f0si]) * (double)(f0fi - f0si);//近くのフレームの前後から直線補完
		


		T0 = (double)fs/ff0;//１周期のサンプル数（実数）
		wLen = (int)(0.5 + T0*2.0);//2周期のサンプル数

		if(wLen+index-(int)(0.5+T0) >= xLen)
		{
			residualSpecgramLength[j] = 0;
			continue;
		}
		residualSpecgramLength[j] = min(fftl-1, wLen);

		for(i = 0;i < residualSpecgramLength[j];i++)
		{
			tmpIndex = i+index - (int)(0.5+T0);//最も近いパルス？の前後1周期分のインデックス
			residualSpecgram[j][i] = x[min( xLen-1, max(0, tmpIndex))];
//			residualSpecgram[j][i] = x[min( xLen-1, max(0, tmpIndex))] * 
//				(0.5 - 0.5*cos(2.0*PI*(double)(i+1)/((double)(wLen+1))));//窓を掛ける
		}
		for(i =  residualSpecgramLength[j]; i < fftl; i++) residualSpecgram[j][i] = 0.0;
	}
}

void PulseResidualWindow(double **residualSpecgram, int *residualSpecgramLength, int pCount)
{
	int i, j;
	for(i = 0;i < pCount-1; i++)
	{
		for(j = 0;j < residualSpecgramLength[i]; j++)
		{
			residualSpecgram[i][j] = residualSpecgram[i][j] * 
				(0.5 - 0.5*cos(2.0*PI*(double)(j+1)/((double)(residualSpecgramLength[i]+1))));//窓を掛ける
		}
	}
}

void getFrameResidualIndex(int tLen, int pCount, double framePeriod, double *pulseLocations, int *residualSpecgramIndex)
{
	int i, j;
	double tmp;
	double tmpValue;
	int tmpIndex;

	for(j = 0;j < tLen;j++)
	{
		tmpValue = 100000.0; // safeGuard
		tmpIndex = pCount-1;
		for(i = 0;i < pCount-1;i++)
		{
			tmp = fabs(pulseLocations[i] +  - (double)j*framePeriod);//フレームに最も近いパルス？を探す
			if(tmp < tmpValue)
			{
				tmpValue = tmp;
				tmpIndex = i;
			}
		}
		residualSpecgramIndex[j] = tmpIndex;//最も近いパルス？のインデックス
	}
}

int getPulseLocations(double *x, int xLen, double *totalPhase, int vuvNum, int *stList, int *edList, int fs, double framePeriod, int *wedgeList, double *pulseLocations)
{
	/*
	int i,j;
	int pCount = 0;
	double wedgePhase;
	for(i = 0;i < vuvNum;i++)
	{
		int stIndex, edIndex; // sample
		// 0点が含まれないように外す
		stIndex = max(0, (int)((double)fs*(stList[i])*framePeriod/1000.0));  //島の先頭のサンプル位置
		edIndex = min(xLen-1, (int)((double)fs*(edList[i]+1)*framePeriod/1000.0+0.5) -1);//島の末尾のサンプル位置
		wedgePhase = totalPhase[wedgeList[i]];
		for(j = stIndex;j < edIndex-1;j++)
			// check Zero-point(find zero cross)
			if(fmod(totalPhase[j+1]-wedgePhase, 2*PI) < fmod(totalPhase[j]-wedgePhase, 2*PI))
			{
				pulseLocations[pCount] =  (double)j/(double)fs;
				pCount++;
			}
	}
	return pCount;
	*/

	int i, j;
	int stIndex, edIndex;

	int pCount = 0;
	int numberOfLocation;
	double *tmpPulseLocations, *basePhase;
	tmpPulseLocations	= (double *)malloc(sizeof(double) * xLen);
	basePhase			= (double *)malloc(sizeof(double) * xLen);


	double tmp;
	for(i = 0;i < vuvNum;i++)
	{
		stIndex = max(0, (int)((double)fs*(stList[i])*framePeriod/1000.0));  //島の先頭のサンプル位置
		edIndex = min(xLen-1, (int)((double)fs*(edList[i]+1)*framePeriod/1000.0+0.5) -1);//島の末尾のサンプル位置

		tmp = totalPhase[wedgeList[i]];

		for(j = stIndex;j < edIndex;j++){
//			basePhase[j] = fmod(totalPhase[j+1]-tmp, 2*PI) - fmod(totalPhase[j]-tmp, 2*PI);//フレーム間の位相差？
			basePhase[j] = fmod(totalPhase[j]-tmp+PI*0.5, 2*PI);  //tn_fnds 各サンプルの位相を補正
		}

//		basePhase[0] = 0; numberOfLocation = 0;
		numberOfLocation = 0;
		for(j = stIndex;j < edIndex-1;j++) 
		{
//			if(fabs(basePhase[j]) > PI/2.0)
			if(basePhase[j+1] < basePhase[j])  //位相が2*PIを超えて0に戻った
			{
				tmpPulseLocations[numberOfLocation++] = (double)j/(double)fs;//ゼロクロス位置の時刻
			}
		}
		for(j = 0;j < numberOfLocation;j++) pulseLocations[pCount++] = tmpPulseLocations[j];//
	}

	free(basePhase);
	free(tmpPulseLocations);
	return pCount;

}

void getWedgeList(double *x, int xLen, int vuvNum, int *stList, int *edList, int fs, double framePeriod, double *f0, int *wedgeList)
{
	int i, j;
	double LowestF0 = 40.0;
	int center, T0;
	double peak;
	int peakIndex = 0;
	double *tmpWav;
	double currentF0;
	tmpWav = (double *)malloc(sizeof(double) * (int)(fs*2/LowestF0));

	for(i = 0;i < vuvNum;i++)
	{
		center		= (int)((stList[i]+edList[i]+1)/2);           //島の中央のフレーム位置
		currentF0	= f0[center] == 0.0 ? DEFAULT_F0 : f0[center];//島の中央のF0 ノイズ領域の場合はデフォルト
		T0			= (int)((fs / currentF0)+0.5);                //島の中央のF0の１周期のサンプル数
//		peakIndex = (int)(((1+center)*framePeriod*fs/1000.0)+0.5);//島の中央のサンプル位置
		peakIndex = (int)(((  center)*framePeriod*fs/1000.0)+0.5);//島の中央のサンプル位置
//		for(j = 0;j < T0*2;j++)
		for(j = 0;j < T0*2+1;j++)
		{
//			tmpWav[j] = x[peakIndex-T0+j-1];
			tmpWav[j] = x[max(0, min(xLen-1, peakIndex-T0+j-1))];//島の中央の２周期分のデータ
		}
		peak = 0.0;
		peakIndex = 0;
		for(j = 0;j < T0*2+1;j++)//波形のピークのサンプル位置を検出
		{
			if(fabs(tmpWav[j]) > peak)
			{
				peak = tmpWav[j];
				peakIndex = j;
			}
		}
//		wedgeList[i] = max(0, min(xLen-1, (int)(0.5 + ((center+1)*framePeriod*fs/1000.0)-T0+peakIndex+1.0) - 1));//島の中央のフレームのピークのサンプル位置
		wedgeList[i] = max(0, min(xLen-1, (int)(0.5 + ((center  )*framePeriod*fs/1000.0)-T0+peakIndex+1.0) - 1));//島の中央のフレームのピークのサンプル位置
	}
	free(tmpWav);
}

// PLATINUM Version 0.0.4.
// Method to extract excitation signals for voice synthesis system based on PLATINUM
// Check: https://www.jstage.jst.go.jp/article/ast/33/2/33_2_123/_pdf
// Check World library before 0.1.4

    void pt100(double *x, int xLen, int fs, double *timeAxis, double *f0,
             double **residualSpecgram)
{
	int i, j, index;
	double framePeriod = (timeAxis[1]-timeAxis[0])*1000.0;

	int	fftl = (int)pow(2.0, 1.0+(int)(log(3.0*fs/FLOOR_F0+1) / log(2.0)));
	int tLen = GetNumDIOSamples(fs, xLen, framePeriod);

	int vuvNum;
	vuvNum = 0;
	for(i = 1;i < tLen;i++)
	{
		if(f0[i]!=0.0 && f0[i-1]==0.0) vuvNum++;
	}
	vuvNum+=vuvNum-1; // 島数の調整 (有声島と無声島)
	if(f0[0] == 0) vuvNum++;
	if(f0[tLen-1] == 0) vuvNum++;

	int stCount, edCount;
	int *stList, *edList;
	stList = (int *)malloc(sizeof(int) * vuvNum);
	edList = (int *)malloc(sizeof(int) * vuvNum);
	edCount = 0;

	stList[0] = 0;
	stCount = 1;
	index = 1;
	if(f0[0] != 0)
	{
		for(i = 1;i < tLen;i++)
		{
			if(f0[i]==0 && f0[i-1]!=0)
			{
				edList[0] = i-1;
				edCount++;
				stList[1] = i;
				stCount++;
				index = i;
			}
		}
	}

	edList[vuvNum-1] = tLen-1;
	for(i = index;i < tLen;i++)
	{
		if(f0[i]!=0.0 && f0[i-1]==0.0) 
		{
			edList[edCount++] = i-1;
			stList[stCount++] = i;
		}
		if(f0[i]==0.0 && f0[i-1]!=0.0) 
		{
			edList[edCount++] = i-1;
			stList[stCount++] = i;
		}
	}

	int *wedgeList;
	wedgeList = (int *)malloc(sizeof(int) * vuvNum);
	getWedgeList(x, xLen, vuvNum, stList, edList, fs, framePeriod, f0, wedgeList);//島中央のピーク位置を取得

	double *signalTime, *f0interpolatedRaw, *totalPhase;
	double *fixedF0;
	fixedF0				= (double *)malloc(sizeof(double) * tLen);
	signalTime			= (double *)malloc(sizeof(double) * xLen);
	f0interpolatedRaw	= (double *)malloc(sizeof(double) * xLen);
	totalPhase			= (double *)malloc(sizeof(double) * xLen);

	for(i = 0;i < tLen;i++) fixedF0[i] = f0[i] == 0 ? DEFAULT_F0 : f0[i]; //F0が0ならデフォルトに補正
	for(i = 0;i < xLen;i++) signalTime[i] = (double)i / (double)fs;       //サンプル位置の時刻
	interp1(timeAxis, fixedF0, tLen, signalTime, xLen, f0interpolatedRaw);//各サンプルのF0
	totalPhase[0] = f0interpolatedRaw[0]*2*PI/(double)fs;                 //各サンプルの位相
	for(i = 1;i < xLen;i++) totalPhase[i] = totalPhase[i-1] + f0interpolatedRaw[i]*2*PI/(double)fs;

	double *pulseLocations;
	pulseLocations		= (double *)malloc(sizeof(double) * xLen);
	int pCount;
	pCount = getPulseLocations(x, xLen, totalPhase, vuvNum, stList, edList, fs, framePeriod, wedgeList, pulseLocations);//位相が大きく動くサンプル位置の時刻

	double *tmpResidualSpec;
	tmpResidualSpec = (double *)malloc(sizeof(double) * fftl);
	double currentF0;

	for(j = 0;j < fftl/2;j++) residualSpecgram[0][j] = 0.0;
	for(i = 1;i < tLen;i++)
	{
		currentF0 = f0[i] <= FLOOR_F0 ? DEFAULT_F0 : f0[i];  //フレームのF0 下限ならデフォルトに補正
		getOneFrameResidualSignal(x, xLen, fs, i, framePeriod/1000.0, currentF0, fftl, pulseLocations, pCount, //最も近いパルス？前後1周期分の波形に窓をかけたものを取得
						tmpResidualSpec);
		for(j = 0;j < fftl/2;j++) residualSpecgram[i][j] = tmpResidualSpec[j];
	}

	free(fixedF0);
	free(tmpResidualSpec);
	free(pulseLocations);
	free(totalPhase); free(f0interpolatedRaw); free(signalTime);
	free(wedgeList);
	free(edList);
	free(stList);
	return;
}

//residualSpecgram 波形のコピーが入る。メモリを確保せずに渡す。この関数により必要なメモリが確保される。
//residualSpecgramLength 波形の長さが入る。メモリを確保せずに渡す。この関数により必要なメモリが確保される。
//residualSpecgramIndex　各フレームの波形を指定するインデックスが入る。
//戻り値波形の数(pCount)
int pt101(double *x, int xLen, int fs, double *timeAxis, double *f0, 
		 double ***residualSpecgram, int **residualSpecgramLength, int *residualSpecgramIndex)
{
	int i, index;
	double framePeriod = (timeAxis[1]-timeAxis[0])*1000.0;

	int	fftl = (int)pow(2.0, 1.0+(int)(log(3.0*fs/FLOOR_F0+1) / log(2.0)));
	int tLen = GetNumDIOSamples(fs, xLen, framePeriod);

	int vuvNum;
//	vuvNum = 0;
	vuvNum = 1;	//tn_fuds
	for(i = 1;i < tLen;i++)
	{
		if(f0[i]!=0.0 && f0[i-1]==0.0) vuvNum++;	//無声→有声
		if(f0[i]==0.0 && f0[i-1]!=0.0) vuvNum++;	//有声→無声  tn_fnds
	}
//	vuvNum+=vuvNum-1; // 島数の調整 (有声島と無声島)  tn_fnds コメントアウト
//	if(f0[0] == 0) vuvNum++;  tn_fnds コメントアウト
//	if(f0[tLen-1] == 0) vuvNum++;  tn_fnds コメントアウト

	int stCount, edCount;
	int *stList, *edList;
	stList = (int *)malloc(sizeof(int) * vuvNum);
	edList = (int *)malloc(sizeof(int) * vuvNum);
	edCount = 0;

	stList[0] = 0;
	stCount = 1;
	index = 1;
	if(f0[0] != 0)	//有声から始まる場合
	{
		for(i = 1;i < tLen;i++)
		{
			if(f0[i]==0 && f0[i-1]!=0)	//有声→無声
			{
				edList[0] = i-1;
				edCount++;
				stList[1] = i;
				stCount++;
				index = i;

				break;	//tn_fnds
			}
		}
	}

	edList[vuvNum-1] = tLen-1;
	for(i = index;i < tLen;i++)
	{
		if(f0[i]!=0.0 && f0[i-1]==0.0) //無声→有声
		{
			edList[edCount++] = i-1;
			stList[stCount++] = i;
		}
		if(f0[i]==0.0 && f0[i-1]!=0.0) //有声→無声
		{
			edList[edCount++] = i-1;
			stList[stCount++] = i;
		}
	}

	int *wedgeList;
	wedgeList = (int *)malloc(sizeof(int) * vuvNum);
	getWedgeList(x, xLen, vuvNum, stList, edList, fs, framePeriod, f0, wedgeList);//島中央のピーク位置を取得

	double *signalTime, *f0interpolatedRaw, *totalPhase;
	double *fixedF0;
	fixedF0				= (double *)malloc(sizeof(double) * tLen);
	signalTime			= (double *)malloc(sizeof(double) * xLen);
	f0interpolatedRaw	= (double *)malloc(sizeof(double) * xLen);
	totalPhase			= (double *)malloc(sizeof(double) * xLen);

	for(i = 0;i < tLen;i++) fixedF0[i] = f0[i] == 0 ? DEFAULT_F0 : f0[i]; //F0が0ならデフォルトに補正
	for(i = 0;i < xLen;i++) signalTime[i] = (double)i / (double)fs;       //サンプル位置の時刻
	interp1(timeAxis, fixedF0, tLen, signalTime, xLen, f0interpolatedRaw);//各サンプルのF0
	totalPhase[0] = f0interpolatedRaw[0]*2*PI/(double)fs;                 //各サンプルの位相
	for(i = 1;i < xLen;i++) totalPhase[i] = totalPhase[i-1] + f0interpolatedRaw[i]*2*PI/(double)fs;

	double *pulseLocations;
	pulseLocations		= (double *)malloc(sizeof(double) * xLen);
	int pCount;
	pCount = getPulseLocations(x, xLen, totalPhase, vuvNum, stList,
				edList, fs, framePeriod, wedgeList, pulseLocations);//位相が大きく動くサンプル位置の時刻

//tn_fndsデバッグ
//zeroXToFile(x, xLen, f0interpolatedRaw, totalPhase, pCount, pulseLocations, fs, vuvNum, wedgeList);

	pCount++;//長さ０のダミーパルスを追加

	*residualSpecgram	= (double **)malloc(sizeof(double *) * pCount);
	for(i = 0;i < pCount;i++) (*residualSpecgram)[i] = (double *)malloc(sizeof(double) * fftl);
	*residualSpecgramLength = (int *)malloc(sizeof(int) * pCount);
	
	getOnePulseResidualSignal(x, xLen, fs, framePeriod/1000.0, f0, fftl, pulseLocations, pCount,
							   *residualSpecgram, *residualSpecgramLength);

	getFrameResidualIndex(tLen, pCount, framePeriod/1000, pulseLocations, residualSpecgramIndex);


//	pulseToFile(pCount, pulseLocations, *residualSpecgramLength);

	free(fixedF0);
	free(pulseLocations);
	free(totalPhase); free(f0interpolatedRaw); free(signalTime);
	free(wedgeList);
	free(edList); free(stList);
	return pCount;
}
