#pragma once
#include "basics.h"
#include <fstream>
#include <string>
//#include<windows.h>
#include<omp.h>
#include <iostream>
#include<vector>
using namespace std;
//NTL_CLIENT

typedef uint64_t TYPE;
uint64_t p = 0x100000000;
//uint64_t p = 0xffffffff;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0, p - 1);
uint64_t rand_64_bit() {
	TYPE left = dis(gen) & 0xffffffff;
	TYPE right = dis(gen) & 0xffffffff;
	return ((left << 32) | right);
}

#define ROL64(x, s) (((x) << (s) & 0xffffffffffffffff) | ((x) >> (64 - (s))))
#define ROR64(x, s) (((x) >> (s)) | ((x) << (64 - (s)) & 0xffffffffffffffff))
void SipHash(uint64_t a_in, uint64_t b_in, uint64_t c_in, uint64_t d_in, uint64_t& a_out, uint64_t& b_out, uint64_t& c_out, uint64_t& d_out, int start, int end)
{
	int l[2] = { 16,21 }, r[2] = { 13,17 };
	uint64_t a = a_in, b = b_in, c = c_in, d = d_in;
	uint64_t z1, z2, w1, w2;
	for (int i = start; i < end; i++)
	{
		z1 = (a + b) & 0xffffffffffffffff;
		w1 = z1 ^ ROL64(b, r[i % 2]);
		z2 = (c + d) & 0xffffffffffffffff;
		w2 = z2 ^ ROL64(d, l[i % 2]);
		a = z2;
		b = w1;
		c = ROL64(z1, 32);
		d = w2;
	}
	a_out = a;
	b_out = b;
	c_out = c;
	d_out = d;
}

void getSipHashExpDlCorDiffWt1MaskWt1_manyDiffMasks(string fileName, uint8_t inDiffIdx[], int diffNum, uint8_t outMaskIdx[], int maskNum, int start, int end, double corWt, uint64_t datasize)
{
	float startTime = omp_get_wtime();
	double* cnt = new double[diffNum * maskNum]();
	ofstream fout(fileName, ios::out);
	for (auto i = 0; i < diffNum * maskNum; i++) cnt[i] = 0;
	for (int64_t expr = 0; expr < datasize; expr++)
	{
		TYPE input[4] = { 0 }, output[4] = { 0 };
		for (auto i = 0; i < 4; i++) input[i] = rand_64_bit() & 0xffffffffffffffff;
		SipHash(input[0], input[1], input[2], input[3], output[0], output[1], output[2], output[3], start, end);
		for (auto j = 0; j < diffNum; j++)
		{
			TYPE inDiff[4] = { 0 };
			uint8_t a = inDiffIdx[j] / 64, b = inDiffIdx[j] % 64;
			inDiff[a] = (TYPE(1) << b);
			TYPE input_prime[4] = { 0 }, output_prime[4] = { 0 };
			for (auto k = 0; k < 4; k++) input_prime[k] = input[k] ^ inDiff[k];
			SipHash(input_prime[0], input_prime[1], input_prime[2], input_prime[3], output_prime[0], output_prime[1], output_prime[2], output_prime[3], start, end);
			TYPE outDiff[4] = { 0 };
			for (auto k = 0; k < 4; k++) outDiff[k] = output[k] ^ output_prime[k];
			for (auto i = 0; i < maskNum; i++)
			{
				uint8_t a = outMaskIdx[i] / 64, b = outMaskIdx[i] % 64;
				bool bit_out = (outDiff[a] >> b) & 0x1;
				if (bit_out == 0) cnt[j * maskNum + i]++; else cnt[j * maskNum + i]--;
			}
		}
	}
	for (auto j = 0; j < diffNum; j++) for (auto i = 0; i < maskNum; i++)
	{
		double cntt = cnt[j * maskNum + i];
		if (cntt)
		{
			double cor = cntt / datasize;
			double cor_wt = log(abs(cor)) / log(2);
			if (cor_wt >= corWt)
			{
				int a = inDiffIdx[j], b = outMaskIdx[i];
				fout << "SipHash: " << dec << (end - start) / 2.0 << "DL(R" << start << " to R" << end << ") experimental COR: ";
				if (a < 64)
				{
					//std::printf("([%d],[],[],[])->", a);
					fout << "([" << a << "],[],[],[])->";
				}
				else if (a < 128)
				{
					//std::printf("([],[%d],[],[])->", a - 64);
					fout << "([],[" << a - 64 << "],[],[])->";
				}
				else if (a < 192)
				{
					//std::printf("([],[],[%d],[])->", a - 64 * 2);
					fout << "([],[],[" << a - 64 * 2 << "],[])->";
				}
				else if (a < 256)
				{
					//std::printf("([],[],[],[%d])->", a - 64 * 3);
					fout << "([],[],[],[" << a - 64 * 3 << "])->";
				}
				if (b < 64)
				{
					//std::printf("([%d],[],[],[]): ", b);
					fout << "([" << b << "],[],[],[]): ";
				}
				else if (b < 128)
				{
					//std::printf("([],[%d],[],[]): ", b - 64);
					fout << "([],[" << b - 64 << "],[],[]): ";
				}
				else if (b < 192)
				{
					//std::printf("([],[],[%d],[]): ", b - 64 * 2);
					fout << "([],[],[" << b - 64 * 2 << "],[]): ";
				}
				else if (b < 256)
				{
					//std::printf("([],[],[],[%d]): ", b - 64 * 3);
					fout << "([],[],[],[" << b - 64 * 3 << "]): ";
				}
				fout << cor << ":" << cor_wt << endl;
			}
		}
	}
	float endTime = omp_get_wtime();
	std::printf("Multi time: %f\n", endTime - startTime);
	fout << "Time: " << endTime - startTime << endl;
	fout.close();
}
void getSipHashExpDlCorGiven1Diff1MaskIdx(string fileName, uint8_t inDiffIdx[], int diffLen, uint8_t outMaskIdx[], int maskLen, int start, int end, uint64_t datasize)
{
	float startTime = omp_get_wtime();
	double cnt = 0;
	int useCoreNumber = thread::hardware_concurrency() / 2;
	omp_set_num_threads(useCoreNumber);
	ofstream fout(fileName, ios::out); //fout << 1 << endl;
#pragma omp parallel for reduction(+:cnt)
	for (int64_t expr = 0; expr < datasize; expr++)
	{
		TYPE input[4] = { 0 }, output[4] = { 0 };
		for (auto i = 0; i < 4; i++) input[i] = rand_64_bit() & 0xffffffffffffffff;
		SipHash(input[0], input[1], input[2], input[3], output[0], output[1], output[2], output[3], start, end);

		TYPE inDiff[4] = { 0 }, input_prime[4] = { 0 }, output_prime[4] = { 0 };
		for (auto i = 0; i < 4; i++) input_prime[i] = input[i];
		for (auto j = 0; j < diffLen; j++)
		{
			int8_t a = inDiffIdx[j] / 64, b = inDiffIdx[j] % 64;
			input_prime[a] ^= ((TYPE)1 << b);
		}
		SipHash(input_prime[0], input_prime[1], input_prime[2], input_prime[3], output_prime[0], output_prime[1], output_prime[2], output_prime[3], start, end);
		TYPE outDiff[4] = { 0x0 };
		for (auto i = 0; i < 4; i++) outDiff[i] = output[i] ^ output_prime[i];
		bool bit_out = 0;
		for (auto i = 0; i < maskLen; i++)
		{
			uint8_t a = outMaskIdx[i] / 64, b = outMaskIdx[i] % 64;
			bit_out ^= (outDiff[a] >> b) & 0x1;
		}
		if (bit_out == 0) cnt++; else cnt--;
	}
	if (cnt == datasize) { cout << "COR=1" << endl; fout << "COR=1" << endl; }
	if (cnt)
	{
		double cntt = cnt / datasize;
		double cor = log(abs(cntt)) / log(2);
		std::printf("SipHash: %dDL(R%d to R%d) experimental COR: [", (end - start) / 2.0, start / 2.0, end / 2.0);
		fout << "SipHash: " << dec << (end - start) / 2.0 << "DL(R" << start / 2.0 << " to R" << end / 2.0 << ") experimental COR: [";
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
		cout << cnt << ", " << cntt << ", " << cor << endl;
		fout << cnt << ", " << cntt << ", " << cor << endl;
	}
	float endTime = omp_get_wtime();
	std::printf("Multi time: %f\n", endTime - startTime);
	fout << "Time: " << endTime - startTime << endl;
	fout.close();
}


int main()
{
#if 0 //verify our new 2.5DL 
	int start = 0, end = 2.5;
	uint8_t inDiffIdx[1] = { 63 + 64 * 2 };
	uint8_t outMaskIdx[1] = { 1 + 64 };
	string fileName = "SipHash 2.5DL(R0 to R2.5), verify v2[63] to v1[1].txt";
	getSipHashExpDlCorGiven1Diff1MaskIdx(fileName, inDiffIdx, 1, outMaskIdx, 1, 2 * start, 2 * end, pow(2, 25));
#endif

#if 1 //verify our new 2.5DL 
	double start = 0, end = 3;
	uint8_t inDiffIdx[1] = { 6 + 64 * 2 };
	uint8_t outMaskIdx[1] = { 19 + 64 };
	string fileName = "SipHash 3DL(R0 to R3), verify v2[6] to v1[19].txt";
	getSipHashExpDlCorGiven1Diff1MaskIdx(fileName, inDiffIdx, 1, outMaskIdx, 1, 2 * start, 2 * end, pow(2, 25));
#endif

#if 0// verify all 3DL with v1[19] or v3[48] as output masks
	int start = 0, end = 3;
	uint8_t inDiffIdx[256] = { 0 }; for (auto i = 0; i < 256; i++) inDiffIdx[i] = i;
	uint8_t outMaskIdx[2] = { 19 + 64, 48 + 64 * 3 };
	string fileName = "SipHash " + to_string(end - start) + "DL (R" + to_string(start) + " to R" + to_string(end) + "), outMask = v1[19] or v3[48], experimental correlations.txt";
	getSipHashExpDlCorDiffWt1MaskWt1_manyDiffMasks(fileName, inDiffIdx, 256, outMaskIdx, 2, 2 * start, 2 * end, -5, pow(2, 25));
#endif

#if 0 //verify 3DL by Niu in [36]
	int start = 0, end = 6;
	uint8_t inDiffIdx[1] = { 63 };
	uint8_t outMaskIdx[8] = { 14,46,14 + 64,46 + 64,14 + 128,46 + 128, 14 + 192,46 + 192 };
	string fileName = "SipHash 3DL (R0 to R3), verify v0[63] to (v0[46,14],v1[46,14],v2[46,14],v3[46,14]).txt";
	getSipHashExpDlCorGiven1Diff1MaskIdx(fileName, inDiffIdx, 1, outMaskIdx, 8, start, end, pow(2, 30));
#endif
}