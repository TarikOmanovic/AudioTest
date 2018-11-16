#include <iostream>
#include <fstream>
#include <vector>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iterator>
#include <intrin.h>

using namespace std;

int
main()
{
	/* RIFF HEADER */

	// RIFF ASCII
	uint32_t ChunkID = _byteswap_ulong(0x52494646);
	// File size, must be modified later
	uint32_t ChunkSize = 88236;
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
	uint32_t AudioFormatAndChannels = ((NumChannels << 16)| AudioFormat);
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
	uint32_t DataChunkSize = 88200;


	
	int i;
	double sinPosLeft = 0.0;
	double sinPosRight = 0.0;
	double toneFreqLeft = 100.0;
	double toneFreqRight = 440.0;
	double sinStepLeft = 2 * M_PI * toneFreqLeft / SampleRate;
	double sinStepRight = 2 * M_PI * toneFreqRight / SampleRate;
	vector<uint32_t> vecWav;

	for (i = 0; i < 88200; ++i)
	{
		/* Just fill the stream with sine! */
		uint32_t left = (uint32_t)(((UINT16_MAX / 2) * sin(sinPosLeft)) + (UINT16_MAX / 2)) << 16;
		uint32_t right = (uint32_t)(((UINT16_MAX / 2) * sin(sinPosRight)) + (UINT16_MAX / 2));
		vecWav.push_back(left | right);
		sinPosLeft += sinStepLeft;
		sinPosRight += sinStepRight;
	}

	


	ofstream out("./example.wav", ios::out | ios::binary);
	ofstream outData("./data.txt", ios::out);
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

	for (const auto &e : vecWav) out.write((char*)&e, 4);
	for (const auto &e : vecWav) outData << e;




	return 0;
}

