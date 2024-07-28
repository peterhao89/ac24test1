#pragma once
#include "basics.h"
#include <fstream>
#include <string>
#include<windows.h>
#include<omp.h>
#include <iostream>
#include<vector>
#include<algorithm>
#include<thread>
using namespace std;
//NTL_CLIENT

typedef uint32_t TYPE;
uint64_t p = 0x100000000;
//uint32_t p = 0x7fffffff;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0, p - 1);

const int rConst[8] = { 31,17,0,24,31,17,0,24 };
const int sConst[8] = { 24,17,31,16,24,17,31,16 };

#define ROL32(x, s) (((x) << (s) & 0xffffffff) | ((x) >> (32 - (s))))
#define ROR32(x, s) (((x) >> (s)) | ((x) << (32 - (s)) & 0xffffffff))

void Alzette(uint32_t a_in, uint32_t b_in, uint32_t& a_out, uint32_t& b_out, int start, int end) //end>start, #round=end-start
{
	uint32_t c_0 = 0xB7E15162;
	uint32_t c_1 = 0xBF715880;
	uint32_t c_2 = 0x38B4DA56;
	uint32_t c_3 = 0x324E7738;
	uint32_t c_4 = 0xBB1185EB;
	uint32_t c_5 = 0x4F7C7B57;
	uint32_t c_6 = 0xCFBFA1C8;
	uint32_t c_7 = 0xC2B3293D;
	uint32_t a = a_in, b = b_in;
	uint32_t key = 0;
	for (int i = start; i < end; i++)
	{
		if (i < 4)
		{
			key = c_0;
		}
		else if (i < 8)
		{
			key = c_2;
		}
		else if (i < 12)
		{
			key = c_6;
		}
		else
		{
			key = c_7;
		}
		a = (a + ROR32(b, rConst[i % 4])) & 0xffffffff;
		b = b ^ ROR32(a, sConst[i % 4]);
		a = a ^ key;
	}
	a_out = a;
	b_out = b;
}

void getAlzetteExpDlCorGiven1Diff1MaskIdx(string fileName, uint8_t inDiffIdx[], int diffLen, uint8_t outMaskIdx[], int maskLen, int start, int end, uint64_t datasize)
{
	float startTime = omp_get_wtime();
	double cnt = 0;
	int useCoreNumber = thread::hardware_concurrency() / 2; cout << useCoreNumber << endl;
	omp_set_num_threads(useCoreNumber);
	ofstream fout(fileName, ios::out);
#pragma omp parallel for reduction(+:cnt)
	for (int64_t expr = 0; expr < datasize; expr++)
	{
		TYPE input[2] = { 0 }, output[2] = { 0 };
		input[0] = (dis(gen)) & 0xffffffff;
		input[1] = (dis(gen)) & 0xffffffff;
		Alzette(input[0], input[1], output[0], output[1], start, end);
		TYPE inDiff[2] = { 0 }, input_prime[2] = { 0 }, output_prime[2] = { 0 };
		input_prime[0] = input[0];
		input_prime[1] = input[1];
		for (auto j = 0; j < diffLen; j++)
		{
			int8_t a = inDiffIdx[j] / 32, b = inDiffIdx[j] % 32;
			input_prime[a] ^= ((TYPE)1 << b);
		}
		Alzette(input_prime[0], input_prime[1], output_prime[0], output_prime[1], start, end);
		TYPE outDiff[2] = { 0x0 };
		outDiff[0] = output[0] ^ output_prime[0];
		outDiff[1] = output[1] ^ output_prime[1];
		bool bit_out = 0;
		for (auto i = 0; i < maskLen; i++)
		{
			uint8_t a = outMaskIdx[i] / 32, b = outMaskIdx[i] % 32;
			bit_out ^= (outDiff[a] >> b) & 0x1;
		}
		if (bit_out == 0) cnt++; else cnt--;
	}
	if (cnt)
	{
		double cntt = cnt / datasize;
		double cor = log(abs(cntt)) / log(2);
		std::printf("Alzette: %dDL(R%d to R%d) experimental COR: [", end - start, start, end);
		fout << "Alzette: " << dec << end - start << "DL(R" << start << " to R" << end << ") experimental COR: [";
		for (auto j = 0; j < diffLen; j++)
		{
			int a = inDiffIdx[j];
			cout << a; fout << a;
			if (j < diffLen - 1) { cout << ",";  fout << ","; }
			else { cout << "]->["; fout << "]->["; }
		}
		for (auto i = 0; i < maskLen; i++)
		{
			int a = outMaskIdx[i];
			cout << a; fout << a;
			if (i < maskLen - 1) { cout << ",";  fout << ","; }
			else { cout << "], "; fout << "], "; }

		}
		cout << cntt << ", " << cor << endl; fout << cntt << ", " << cor << endl;
	}
	float endTime = omp_get_wtime();
	std::printf("Multi time: %f\n", endTime - startTime);
	fout << "Multi time:" << endTime - startTime << endl;
	fout.close();
}

struct outMaskConfig
{
	uint8_t idx;
	bool single_or_consecutive = 1; //1, single; 0, consecutive
};
void getAlzetteExpDlCorGiven1Diff2SingleOrConsecutiveMask(string fileName, uint8_t inDiffIdx[], int diffLen, outMaskConfig outMaskIdx[2], int start, int end, uint64_t datasize)
{
	float startTime = omp_get_wtime();
	double cnt0 = 0, cnt1 = 0;
	ofstream fout(fileName, ios::out); //fout << 1 << endl;
	int useCoreNumber = thread::hardware_concurrency() / 2; cout << useCoreNumber << endl;
	omp_set_num_threads(useCoreNumber);
#pragma omp parallel for reduction(+:cnt0,cnt1)
	for (int64_t expr = 0; expr < datasize; expr++)
	{
		TYPE input[2] = { 0 }, output[2] = { 0 };
		input[0] = (dis(gen)) & 0xffffffff;
		input[1] = (dis(gen)) & 0xffffffff;
		Alzette(input[0], input[1], output[0], output[1], start, end);
		TYPE inDiff[2] = { 0 }, input_prime[2] = { 0 }, output_prime[2] = { 0 };
		input_prime[0] = input[0];
		input_prime[1] = input[1];
		for (auto j = 0; j < diffLen; j++)
		{
			int8_t a = inDiffIdx[j] / 32, b = inDiffIdx[j] % 32;
			input_prime[a] ^= ((TYPE)1 << b);
		}
		Alzette(input_prime[0], input_prime[1], output_prime[0], output_prime[1], start, end);
		TYPE outDiff[2] = { 0x0 };
		outDiff[0] = output[0] ^ output_prime[0];
		outDiff[1] = output[1] ^ output_prime[1];
		bool bit_out0 = 0;
		if (outMaskIdx[0].single_or_consecutive)
		{
			uint8_t a = outMaskIdx[0].idx / 32, b = outMaskIdx[0].idx % 32;
			bit_out0 = (outDiff[a] >> b) & 0x1;
		}
		else
		{
			uint8_t a1 = (outMaskIdx[0].idx + 1) / 32, b1 = (outMaskIdx[0].idx + 1) % 32;
			uint8_t a0 = outMaskIdx[0].idx / 32, b0 = outMaskIdx[0].idx % 32;
			bit_out0 = ((outDiff[a1] >> b1) & 0x1) ^ ((outDiff[a0] >> b0) & 0x1);
		}
		if (bit_out0 == 0) cnt0++; else cnt0--;
		bool bit_out1 = 0;
		if (outMaskIdx[1].single_or_consecutive)
		{
			uint8_t a = outMaskIdx[1].idx / 32, b = outMaskIdx[1].idx % 32;
			bit_out1 = (outDiff[a] >> b) & 0x1;
		}
		else
		{
			uint8_t a1 = (outMaskIdx[1].idx + 1) / 32, b1 = (outMaskIdx[1].idx + 1) % 32;
			uint8_t a0 = outMaskIdx[1].idx / 32, b0 = outMaskIdx[1].idx % 32;
			bit_out1 = ((outDiff[a1] >> b1) & 0x1) ^ ((outDiff[a0] >> b0) & 0x1);
		}
		if (bit_out1 == 0) cnt1++; else cnt1--;
	}
	double cntt0 = cnt0 / datasize, cntt1 = cnt1 / datasize;
	double cor0 = log(abs(cntt0)) / log(2), cor1 = log(abs(cntt1)) / log(2);
	for (auto i = 0; i < 2; i++)
	{
		std::printf("Alzette: %dDL(R%d to R%d) experimental COR: [", end - start, start, end);
		fout << "Alzette: " << dec << end - start << "DL(R" << start << " to R" << end << ") experimental COR: [";
		for (auto j = 0; j < diffLen; j++)
		{
			int a = inDiffIdx[j];
			cout << a; fout << a;
			if (j < diffLen - 1) { cout << ",";  fout << ","; }
			else { cout << "]->["; fout << "]->["; }
		}
		int b = outMaskIdx[i].idx;
		if (outMaskIdx[i].single_or_consecutive)
		{
			cout << b << "], "; fout << b << "], ";
		}
		else
		{
			cout << b + 1 << "," << b << "], "; fout << b + 1 << "," << b << "], ";
		}
		if (i == 0) { cout << cntt0 << ", " << cor0 << endl; fout << cntt0 << ", " << cor0 << endl; }
		else { cout << cntt1 << ", " << cor1 << endl; fout << cntt1 << ", " << cor1 << endl; }
	}

	float endTime = omp_get_wtime();
	std::printf("Multi time: %f\n", endTime - startTime);
	fout << "Multi time:" << endTime - startTime << endl;
	fout.close();
}

void getAlzetteSubsSetDiffWt1(string fileName, uint8_t inDiffIdx[], int diffNum, uint8_t outMaskIdx[], int maskNum, int start, int end, double corWt, uint64_t datasize)
{
	float start_time = omp_get_wtime();
	double* cnt = new double[diffNum * maskNum]();
	ofstream fout(fileName, ios::out);
	for (auto i = 0; i < diffNum * maskNum; i++) cnt[i] = 0;
	fout << "Threshold of the absolute correlation weight:  " << -1 * corWt << endl;
	for (int64_t exper = 0; exper < datasize; exper++) //for (int x = 0; x < range; x++)  for (int y = 0; y < range; y++)
	{
		TYPE input[2] = { 0 }, output[2] = { 0 };
		input[0] = (dis(gen)) & 0xffffff;
		input[1] = (dis(gen)) & 0xffffff;
		TYPE output_a = 0x0, output_b = 0x0;
		Alzette(input[0], input[1], output[0], output[1], start, end);
		for (auto j = 0; j < diffNum; j++)
		{
			TYPE inDiff[2] = { 0 };
			uint8_t a = inDiffIdx[j] / 32, b = inDiffIdx[j] % 32;
			inDiff[a] = ROL32(1, b); // ((TYPE)1 << b);
			TYPE input_prime[2] = { 0 }, output_prime[2] = { 0 };
			input_prime[0] = input[0] ^ inDiff[0];
			input_prime[1] = input[1] ^ inDiff[1];
			Alzette(input_prime[0], input_prime[1], output_prime[0], output_prime[1], start, end);
			TYPE outDiff[2] = { 0x0 };
			outDiff[0] = output[0] ^ output_prime[0];
			outDiff[1] = output[1] ^ output_prime[1];
			for (auto i = 0; i < maskNum; i++)
			{
				uint8_t a = outMaskIdx[i] / 32, b = outMaskIdx[i] % 32;
				bool bit_out = (outDiff[a] >> b) & 0x1;
				if (bit_out == 0)
				{
					cnt[j * maskNum + i]++;
				}
				else
				{
					cnt[j * maskNum + i]--;
				}
			}
		}
	}
	for (auto j = 0; j < diffNum; j++)
	{
		int totStrongNum = 0;
		int a = inDiffIdx[j];
		fout << "bit index: ";
		if (a < 32) { fout << a + 32 << endl; }
		else { fout << a - 32 << endl; }
		fout << "strong unbalanced bits are ";
		//int a1 = a / 32, a0 = a % 32;
		//cout << "v" << a1 << "[" << a0 << "]: ";
		//fout << "v" << a1 << "[" << a0 << "]: ";
		for (auto i = 0; i < maskNum; i++)
		{
			double cor = cnt[j * maskNum + i];
			if (cor)
			{
				cor = cor / datasize;
				cor = log(abs(cor)) / log(2);
				if (cor >= corWt)
				{
					totStrongNum++;
					int b = outMaskIdx[i];
					if (b < 32) { fout << b + 32 << "(" << cor << ")  "; }
					else { fout << b - 32 << "(" << cor << ")  "; }
				}
			}
		}
		fout << endl;
		fout << "the number of strong unbalanced bit is " << totStrongNum << endl;
		//cout << totStrongNum << endl;
		//fout << totStrongNum << endl;
	}
	delete[] cnt;
	float end_time = omp_get_wtime();
	std::printf("Multi time: %f\n", end_time - start_time);
	fout << "Multi time:" << end_time - start_time << endl;
	fout.close();
}

void getAlzetteExpLcCorGivenInOutMasks(string fileName, uint32_t inMask[2], uint32_t outMask[2], int start, int end, uint64_t datasize)
{
	float startTime = omp_get_wtime();
	double cnt = 0;
	int useCoreNumber = thread::hardware_concurrency() / 2;
	omp_set_num_threads(useCoreNumber);
	ofstream fout(fileName, ios::out | ios::app);
#pragma omp parallel for reduction(+:cnt)
	for (int64_t expr = 0; expr < datasize; expr++)
	{
		TYPE input[2] = { 0 }, output[2] = { 0 };
		input[0] = (dis(gen)) & 0xffffffff;
		input[1] = (dis(gen)) & 0xffffffff;
		Alzette(input[0], input[1], output[0], output[1], start, end);
		uint8_t bit_out = maskdot(inMask[0], input[0], 32) ^ maskdot(inMask[1], input[1], 32) ^ maskdot(outMask[0], output[0], 32) ^ maskdot(outMask[1], output[1], 32);
		if (bit_out == 0) cnt++; else cnt--;
	}
	if (cnt)
	{
		double cntt = cnt / datasize;
		double cor = log(abs(cntt)) / log(2);
		std::printf("Alzette: %dLC(R%d to R%d) experimental COR: ", end - start, start, end);
		printf("(%x,%x)->(%x,%x), ", inMask[0], inMask[1], outMask[0], outMask[1]);
		fout << "Alzette: " << dec << end - start << "LC(R" << start << " to R" << end << ") experimental COR: ";
		fout << "(" << hex << inMask[0] << "," << inMask[1] << ")->(" << outMask[0] << "," << outMask[1] << "), ";
		cout << cntt << ", " << cor << endl; fout << cntt << ", " << cor << endl;
	}
	float endTime = omp_get_wtime();
	std::printf("Multi time: %f\n", endTime - startTime);
	fout << "Multi time:" << endTime - startTime << endl;
	fout.close();
}


void main()
{
#if 1 //verify 7DL for our new 13- and 14-round DL
	int start = 3, end = 10;
	uint8_t inDiffIdx[1] = { 29 };
	uint8_t outMaskIdx[2] = { 1 + 32 };
	string fileName = "Alzette " + to_string(end - start) + "DL(R" + to_string(start) + " to R" + to_string(end) + "), verify [29][] to [][1].txt";
	getAlzetteExpDlCorGiven1Diff1MaskIdx(fileName, inDiffIdx, 1, outMaskIdx, 1, start, end, pow(2, 34));
#endif


#if 0 //verify 6DL for our new 11-round DL
	int start = 2, end = 8;
	uint8_t inDiffIdx[1] = { 16 + 32 };
	uint8_t outMaskIdx[2] = { 6 + 32, 5 + 32 };
	string fileName = "Alzette " + to_string(end - start) + "DL(R" + to_string(start) + " to R" + to_string(end) + "), verify [][16] to [][6,5].txt";
	getAlzetteExpDlCorGiven1Diff1MaskIdx(fileName, inDiffIdx, 1, outMaskIdx, 2, start, end, pow(2, 34));
#endif

#if 0 //verify 7DL for our new 12-round DL
	int start = 2, end = 9;
	uint8_t inDiffIdx[1] = { 25 + 32 };
	uint8_t outMaskIdx[2] = { 15 + 32, 14 + 32 };
	string fileName = "Alzette " + to_string(end - start) + "DL(R" + to_string(start) + " to R" + to_string(end) + "), verify [][25] to [][15,14].txt";
	getAlzetteExpDlCorGiven1Diff1MaskIdx(fileName, inDiffIdx, 1, outMaskIdx, 2, start, end, pow(2, 38));
#endif


#if 0//prepare SUBs sets for 6-round DL
	int start = 2, end = 8;
	uint8_t inDiffIdx[64] = { 0 }; for (auto i = 0; i < 32; i++) { inDiffIdx[i] = i + 32; inDiffIdx[i + 32] = i; }
	uint8_t outMaskIdx[64] = { 0 }; for (auto i = 0; i < 32; i++) { outMaskIdx[i] = i + 32; outMaskIdx[i + 32] = i; }
	string fileName = "Alzette " + to_string(end - start) + "DL(R" + to_string(start) + " to R" + to_string(end) + ") SUBs sets.txt";
	getAlzetteSubsSetDiffWt1(fileName, inDiffIdx, 64, outMaskIdx, 64, start, end, -9.5, pow(double(2), 31));
#endif

#if 0
	uint32_t endMask[][2] = { {0x8b0a0205,0x81020484},{0x8e0a0205,0x81020684},{0x898a0205,0x81020584},{0x898a0205,0x81020584} };
	int start = 10, end = 14;
	uint32_t inMask[2] = { 0, getInt32(1) };
	string fileName = "Alzette " + to_string(end - start) + "LC(R" + to_string(start) + " to R" + to_string(end) + "), verify.txt";
	for (auto i = 0; i < 4; i++)
	{
		getAlzetteExpLcCorGivenInOutMasks(fileName, inMask, endMask[i], start, end, pow(2, 26));
	}
#endif
}