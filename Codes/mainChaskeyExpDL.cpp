#pragma once
#include "basics.h"
#include <fstream>
#include <string>
#include<windows.h>
#include<omp.h>
#include <iostream>
#include<vector>
using namespace std;
//NTL_CLIENT

typedef uint32_t TYPE;
uint64_t p = 0x100000000;
//uint32_t p = 0x7fffffff;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0, p - 1);

#define ROL32(x, s) (((x) << (s) & 0xffffffff) | ((x) >> (32 - (s))))
#define ROR32(x, s) (((x) >> (s)) | ((x) << (32 - (s)) & 0xffffffff))

void chaskey_permutation(uint32_t v0, uint32_t v1, uint32_t v2, uint32_t v3, uint32_t& vv0, uint32_t& vv1, uint32_t& vv2, uint32_t& vv3, int start, int end)
{
	uint32_t z0, z1, w0, w1;
	int l[2] = { 5,7 }, r[2] = { 8,13 };
	for (int i = start; i < end; i++)
	{
		z0 = (v0 + v1) & 0xffffffff;
		w0 = z0 ^ ROL32(v1, r[i % 2]);
		z1 = (v2 + v3) & 0xffffffff;
		w1 = z1 ^ ROL32(v3, l[i % 2]);
		v0 = ROL32(z1, 16);
		v1 = w0;
		v2 = z0;
		v3 = w1;
	}
	vv0 = v0;
	vv1 = v1;
	vv2 = v2;
	vv3 = v3;
}

void getChaskeyExpDlCorGiven1Diff1MaskIdx(string fileName, uint8_t oneInDiffIdx[], int diffLen, uint8_t oneOutMaskIdx[], int maskLen, int start, int end, uint64_t datasize)
{
	float start_time = omp_get_wtime();
	double cnt = 0;
	ofstream fout(fileName, ios::out);
	int useCoreNumber = thread::hardware_concurrency() / 2; cout << useCoreNumber << endl;
	omp_set_num_threads(useCoreNumber);
#pragma omp parallel for reduction(+:cnt)
	for (int64_t exper = 0; exper < datasize; exper++)
	{
		TYPE input[4] = { 0 }, output[4] = { 0 };
		for (auto i = 0; i < 4; i++) input[i] = (dis(gen)) & 0xffffffff;
		chaskey_permutation(input[0], input[1], input[2], input[3], output[0], output[1], output[2], output[3], start, end);
		//std::printf("(%x,%x,%x,%x),(%x,%x,%x,%x)\n", input[0], input[1], input[2], input[3], output[0], output[1], output[2], output[3]);
		TYPE input_prime[4] = { 0 }, output_prime[4] = { 0 };
		for (auto i = 0; i < 4; i++) input_prime[i] = input[i];
		for (auto j = 0; j < diffLen; j++)
		{
			int8_t a = oneInDiffIdx[j] / 32, b = oneInDiffIdx[j] % 32;
			input_prime[a] ^= ((TYPE)1 << b);
		}
		chaskey_permutation(input_prime[0], input_prime[1], input_prime[2], input_prime[3], output_prime[0], output_prime[1], output_prime[2], output_prime[3], start, end);
		//std::printf("(%x,%x,%x,%x),(%x,%x,%x,%x)\n", input_prime[0], input_prime[1], input_prime[2], input_prime[3], output_prime[0], output_prime[1], output_prime[2], output_prime[3]);
		TYPE outDiff[4] = { 0x0 };
		for (auto i = 0; i < 4; i++) outDiff[i] = output_prime[i] ^ output[i];
		bool bit_out = 0;
		for (auto i = 0; i < maskLen; i++)
		{
			uint8_t a = oneOutMaskIdx[i] / 32, b = oneOutMaskIdx[i] % 32;
			bit_out ^= (outDiff[a] >> b) & 0x1;
		}
		if (bit_out == 0) cnt++;
		else cnt--;
	}
	double cntD = cnt / datasize;
	double cor = log(abs(cntD)) / log(2);
	std::printf("Chaskey %dDL(r%d->r%d) experimental COR: ", end - start, start, end);
	fout << "Chaskey " << dec << end - start << "DL(r" << start << "->r" << end << ") experimental COR: ";
	for (auto j = 0; j < diffLen; j++)
	{
		int a = oneInDiffIdx[j];
		int a0 = a / 32, a1 = a % 32;
		cout << "v" << a0 << "[" << a1 << "]";
		fout << "v" << a0 << "[" << a1 << "]";
		if (j < diffLen - 1) { cout << ",";  fout << ","; }
		else { cout << "->"; fout << "->"; }
	}
	for (auto i = 0; i < maskLen; i++)
	{
		int a = oneOutMaskIdx[i];
		//cout << a; fout << a;
		int a0 = a / 32, a1 = a % 32;
		cout << "v" << a0 << "[" << a1 << "], ";
		fout << "v" << a0 << "[" << a1 << "], ";
	}
	cout << cntD << ", " << cor << endl; fout << cntD << ", " << cor << endl;
	float end_time = omp_get_wtime();
	std::printf("Multi time: %f\n", end_time - start_time);
	fout << "Multi time: " << end_time - start_time << endl;
	fout.close();
}

void main()
{
#if 1 // Chaskey 7DL=1.5DC + 4DL + 1.5LC, verify our new 4DL, v0[31] to v3[11]
	double start = 1.5, end = 5.5;
	uint8_t oneInDiffIdx[1] = { 31 + 32 * 0 };
	uint8_t oneOutMaskIdx[1] = { 11 + 32 * 3 };
	string fileName = "Chaskey 4DL(R1.5 to R5.5), verify v0[31] to v3[11].txt";
	getChaskeyExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, 2 * start, 2 * end, pow(2, 34));
#endif
}