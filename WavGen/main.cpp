#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iterator>
#include <intrin.h>

using namespace std;


/* RIFF HEADER */

// RIFF ASCII
uint32_t ChunkID = _byteswap_ulong(0x52494646);
// File size, must be modified later
uint32_t ChunkSize = 36;
// WAV ASCII
uint32_t Format = _byteswap_ulong(0x57415645);

/* FMT CHUNK */

// fmt ASCII
uint32_t FromatChunkID = _byteswap_ulong(0x666d7420);
// size of chunk after this entry in bytes
uint32_t FormatChunkSize = 16;
// uncompressed audio
uint32_t AudioFormat = 1;
// stereo
uint32_t NumChannels = 2;
uint32_t AudioFormatAndChannels = ((NumChannels << 16) | AudioFormat);
// standard sample rate
uint32_t SampleRate = 44100;
// standard bps
uint32_t BitsPerSample = 16;
uint32_t ByteRate = SampleRate * NumChannels * BitsPerSample / 8;
uint32_t BlockAlign = (NumChannels * BitsPerSample / 8);
uint32_t BlockAlignAndBitsPerSample = ((BitsPerSample << 16) | BlockAlign);

/* DATA CHUNK */

// data ASCII
uint32_t DataChunkID = _byteswap_ulong(0x64617461);
uint32_t DataChunkSize;

void
tone(
	vector<int16_t>& vec,
	double durationSec,
	double frequencyHzBeg,
	double frequencyHzEnd,
	double volumeScaleBeg,
	double volumeScaleEnd,
	void(*interpolationFreqCurve)(vector<double>& out, double fB, double fE, unsigned int numOfSamples),
	void(*interpolationVolumeCurve)(vector<double>& out, double vB, double vE, unsigned int numOfSamples));

void
inverseHamronicTone(
	vector<int16_t>& vec,
	double durationSec,
	double frequencyHz,
	double volumeScale,
	unsigned int numHarmonics);

void
linearHamronicTone(
	vector<int16_t>& vec,
	double durationSec,
	double frequencyHz,
	double volumeScale,
	unsigned int numHarmonics);

void
linearInterpolation(
	vector<double>& out,
	double b,
	double e,
	unsigned int numOfSamples);

int
main()
{
	vector<int16_t> chLeft;
	vector<int16_t> chRight;
	
	tone(chLeft, 4, 220, 880, 0, 0.5, linearInterpolation, linearInterpolation);
	tone(chRight, 4, 220, 880, 0, 0.5, linearInterpolation, linearInterpolation);
	tone(chLeft, 4, 880, 220, 0.5, 1, linearInterpolation, linearInterpolation);
	tone(chRight, 4, 880, 220, 0.5, 1, linearInterpolation, linearInterpolation);
	
	   
	DataChunkSize = NumChannels * sizeof(chLeft[0]) * chLeft.size();
	ChunkSize = DataChunkSize + 36;
	   	  
	ofstream out("./example.wav", ios::out | ios::binary);
	out.write((char*)&ChunkID, 4);
	out.write((char*)&ChunkSize, 4);
	out.write((char*)&Format, 4);
	out.write((char*)&FromatChunkID, 4);
	out.write((char*)&FormatChunkSize, 4);
	out.write((char*)&AudioFormatAndChannels, 4);
	out.write((char*)&SampleRate, 4);
	out.write((char*)&ByteRate, 4);
	out.write((char*)&BlockAlignAndBitsPerSample, 4);
	out.write((char*)&DataChunkID, 4);
	out.write((char*)&DataChunkSize, 4);

	unsigned int i;
	for (i = 0; i < chLeft.size(); ++i)
	{
		out.write((char*)&(chLeft[i]), 2);
		out.write((char*)&(chRight[i]), 2);
	}

	return 0;
}



void
tone(
	vector<int16_t>& vec,
	double durationSec,
	double frequencyHzBeg,
	double frequencyHzEnd,
	double volumeScaleBeg,
	double volumeScaleEnd,
	void(*interpolationFreqCurve)(vector<double>& out, double fB, double fE, unsigned int numOfSamples),
	void(*interpolationVolumeCurve)(vector<double>& out, double vB, double vE, unsigned int numOfSamples)
	)
{
	unsigned int numOfSamples = (unsigned int)(durationSec * SampleRate);

	vector<double> frequencies;
	vector<double> volumes;
	interpolationFreqCurve(frequencies, frequencyHzBeg, frequencyHzEnd, numOfSamples);
	interpolationVolumeCurve(volumes, volumeScaleBeg, volumeScaleEnd, numOfSamples);
	
	double sinPos = 0.0;
	for (unsigned int i = 0; i < numOfSamples; ++i)
	{
		vec.push_back((int16_t)(INT16_MAX * volumes[i] * sin(sinPos)));
		sinPos += 2 * M_PI * frequencies[i] / (double)SampleRate;
	}
}




void
linearInterpolation(
	vector<double>& out,
	double b,
	double e,
	unsigned int numOfSamples)
{
	double delta = (e - b) / numOfSamples;
	for (unsigned int i = 0; i < numOfSamples; ++i)
	{
		out.push_back(b + (i * delta));
	}
}








void
inverseHamronicTone(
	vector<int16_t>& vec,
	double durationSec,
	double frequencyHz,
	double volumeScale,
	unsigned int numHarmonics)
{
	unsigned int i;
	unsigned int j;
	unsigned int numOfSamples = (unsigned int)(durationSec * SampleRate);

	double sinPos = 0.0;
	double sinStep = 2 * M_PI * frequencyHz / (double)SampleRate;

	double factor = 0;
	for (i = 1; i <= numHarmonics; ++i)
	{
		factor += (1.0 / i);
	}
	volumeScale /= factor;

	int16_t val;

	vector<int16_t> tone;

	for (i = 0; i < numOfSamples; ++i)
	{
		val = (int16_t)(INT16_MAX * volumeScale * sin(sinPos));
		tone.push_back(val);
		sinPos += sinStep;
	}

	for (j = 2; j <= numHarmonics; ++j)
	{
		cout << j << endl;
		sinPos = 0.0;
		sinStep = 2 * M_PI * j * frequencyHz / (double)SampleRate;

		for (i = 0; i < numOfSamples; ++i)
		{
			val = (int16_t)(INT16_MAX * (volumeScale / j) * sin(sinPos));
			tone[i] += val;
			sinPos += sinStep;
		}
	}

	for (auto t : tone)
	{
		vec.push_back(t);
	}
}


void
linearHamronicTone(
	vector<int16_t>& vec,
	double durationSec,
	double frequencyHz,
	double volumeScale,
	unsigned int numHarmonics)
{
	unsigned int i;
	unsigned int j;
	unsigned int numOfSamples = (unsigned int)(durationSec * SampleRate);

	double sinPos = 0.0;
	double sinStep = 2 * M_PI * frequencyHz / (double)SampleRate;

	double factor = 0;
	for (i = 0; i < numHarmonics; ++i)
	{
		factor += (numHarmonics - i);
	}
	volumeScale /= factor;

	int16_t val;

	vector<int16_t> tone;

	for (i = 0; i < numOfSamples; ++i)
	{
		val = (int16_t)(INT16_MAX * volumeScale * (numHarmonics) * sin(sinPos));
		tone.push_back(val);
		sinPos += sinStep;
	}

	for (j = 2; j <= numHarmonics; ++j)
	{
		cout << j << endl;
		sinPos = 0.0;
		sinStep = 2 * M_PI * j * frequencyHz / (double)SampleRate;

		for (i = 0; i < numOfSamples; ++i)
		{
			val = (int16_t)(INT16_MAX * volumeScale * (numHarmonics - j + 1) * sin(sinPos));
			tone[i] += val;
			sinPos += sinStep;
		}
	}

	for (auto t : tone)
	{
		vec.push_back(t);
	}
}

