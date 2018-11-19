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
simpleTone(
	vector<int16_t>& vec,
	double durationSec,
	double frequencyHz,
	double volumeScale);

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

int
main()
{
	vector<int16_t> chLeft;
	vector<int16_t> chRight;
	simpleTone(chLeft, 2.0, 440, 1);
	simpleTone(chRight, 2.0, 440, 1);
	inverseHamronicTone(chLeft, 2.0, 440, 1, 2);
	inverseHamronicTone(chRight, 2.0, 440, 1, 2);
	linearHamronicTone(chLeft, 2.0, 440, 1, 2);
	linearHamronicTone(chRight, 2.0, 440, 1, 2);
	   
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
simpleTone(
	vector<int16_t>& vec,
	double durationSec,
	double frequencyHz,
	double volumeScale)
{
	unsigned int i;
	unsigned int numOfSamples = (unsigned int) (durationSec * SampleRate);

	double sinPos = 0.0;
	double sinStep = 2 * M_PI * frequencyHz / (double)SampleRate;

	uint16_t val;

	for (i = 0; i < numOfSamples; ++i)
	{
		val = (int16_t)(INT16_MAX * volumeScale * sin(sinPos));
		vec.push_back(val);
		sinPos += sinStep;
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

