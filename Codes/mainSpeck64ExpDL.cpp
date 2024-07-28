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

void speck_64(uint32_t a_in, uint32_t b_in, uint32_t& a_out, uint32_t& b_out, int start, int end, uint32_t kk[])
{
	uint32_t a = a_in, b = b_in;
	for (int i = start; i < end; i++)
	{
		a = ROR32(a, 8);
		a = (a + b) & 0xffffffff;
		a = a ^ kk[i];
		b = ROL32(b, 3) ^ a;
	}
	a_out = a;
	b_out = b;
}
void getSpeck64ExpDlCorGiven1Diff1MaskIdx(string fileName, uint8_t oneInDiffIdx[], int diffLen, uint8_t oneOutMaskIdx[], int maskLen, int start, int end, uint64_t datasize, int testNum)
{
	int round = end - start;
	TYPE master_key[4], mid;
	TYPE* kk = new TYPE[end + 1]();
	TYPE* l = new TYPE[end + 3]();
	float start_time = omp_get_wtime();
	double tot = 0;
	ofstream fout(fileName, ios::out);
	int useCoreNumber = thread::hardware_concurrency() / 2; cout << useCoreNumber << endl;
	omp_set_num_threads(useCoreNumber);
	for (int r = 0; r < testNum; r++)
	{
		for (int i = 0; i < 4; i++)
		{
			master_key[i] = (dis(gen)) & 0xffffffff;
		}
		kk[0] = master_key[0];
		for (int i = 0; i < 3; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < end; i++)
		{
			mid = ((ROR32(l[i], 8) + kk[i]) & 0xffffffff) ^ (i & 0xffffffff);
			kk[i + 1] = mid ^ (ROL32(kk[i], 3));
			l[i + 3] = mid;
		}
		double cnt = 0;
#pragma omp parallel for reduction(+:cnt,tot)
		for (int64_t exper = 0; exper < datasize; exper++)
		{
			TYPE input[2] = { 0 };
			input[0] = (dis(gen)) & 0xffffffff;
			input[1] = (dis(gen)) & 0xffffffff;
			TYPE output_a = 0x0, output_b = 0x0;
			speck_64(input[0], input[1], output_a, output_b, start, end, kk);
			TYPE inDiff[2] = { 0 }, input_prime[2] = { 0 };
			input_prime[0] = input[0];
			input_prime[1] = input[1];
			for (auto j = 0; j < diffLen; j++)
			{
				int8_t a = oneInDiffIdx[j] / 32, b = oneInDiffIdx[j] % 32;
				input_prime[a] ^= ((TYPE)1 << b);
			}
			TYPE output_a_prime = 0x0, output_b_prime = 0x0;
			speck_64(input_prime[0], input_prime[1], output_a_prime, output_b_prime, start, end, kk);
			TYPE outDiff[2] = { 0x0 };
			outDiff[0] = output_a ^ output_a_prime;
			outDiff[1] = output_b ^ output_b_prime;
			bool bit_out = 0;
			for (auto i = 0; i < maskLen; i++)
			{
				uint8_t a = oneOutMaskIdx[i] / 32, b = oneOutMaskIdx[i] % 32;
				bit_out ^= (outDiff[a] >> b) & 0x1;
			}
			if (bit_out == 0)
			{
				cnt++; tot++;
			}
			else
			{
				cnt--; tot--;
			}
		}
		double cntD = cnt / datasize;
		double cor = log(abs(cntD)) / log(2);
		std::printf("Speck64 %dDL experimental COR: [", round);
		fout << "Speck64 " << dec << round << "DL experimental COR : [";
		for (auto j = 0; j < diffLen; j++)
		{
			int a = oneInDiffIdx[j];
			cout << a; fout << a;
			if (j < diffLen - 1) { cout << ",";  fout << ","; }
			else { cout << "]->["; fout << "]->["; }
		}
		for (auto i = 0; i < maskLen; i++)
		{
			int a = oneOutMaskIdx[i];
			cout << a; fout << a;
			if (i < maskLen - 1) { cout << ",";  fout << ","; }
			else { cout << "], "; fout << "], "; }

		}
		cout << cntD << ", " << cor << endl; fout << cntD << ", " << cor << endl;
	}
	if (tot)
	{
		double totD = tot / datasize / testNum;
		double cor = log(abs(totD)) / log(2);
		std::printf("Speck64 %dDL AVERAGE experimental COR: [", round);
		fout << "Speck64 " << dec << round << "DL AVERAGE experimental COR : [";
		for (auto j = 0; j < diffLen; j++)
		{
			int a = oneInDiffIdx[j];
			cout << a; fout << a;
			if (j < diffLen - 1) { cout << ",";  fout << ","; }
			else { cout << "]->["; fout << "]->["; }
		}
		for (auto i = 0; i < maskLen; i++)
		{
			int a = oneOutMaskIdx[i];
			cout << a; fout << a;
			if (i < maskLen - 1) { cout << ",";  fout << ","; }
			else { cout << "], "; fout << "], "; }

		}
		cout << totD << ", " << cor << endl; fout << totD << ", " << cor << endl;
	}
	float end_time = omp_get_wtime();
	std::printf("Multi time: %f\n", end_time - start_time);
	fout << "Multi time:" << end_time - start_time << endl;
	fout.close();
}

void getSpeck64SubsSetDiffWt1(string fileName, uint8_t inDiffIdx[], int diffNum, uint8_t outMaskIdx[], int maskNum, int start, int end, double corWt, uint64_t datasize, int testNum)
{
	int round = end - start;
	TYPE master_key[4], mid;
	TYPE* kk = new TYPE[end + 1]();
	TYPE* l = new TYPE[end + 3]();
	float start_time = omp_get_wtime();
	double* tot = new double[diffNum * maskNum]();
	double* cnt = new double[diffNum * maskNum]();
	ofstream fout(fileName, ios::out);
	fout << "Threshold of the absolute correlation weight:  " << -1 * corWt << endl;
	for (int r = 0; r < testNum; r++)
	{
		for (int i = 0; i < 4; i++)
		{
			master_key[i] = (dis(gen)) & 0xffffffff;
		} //cout << r << ": the master key is "; std::printf("(%x,%x,%x,%x)\n", master_key[0], master_key[1], master_key[2], master_key[3]);
		kk[0] = master_key[0];
		for (int i = 0; i < 3; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < end; i++)
		{
			mid = ((ROR32(l[i], 8) + kk[i]) & 0xffffffff) ^ (i & 0xffffffff);
			kk[i + 1] = mid ^ (ROL32(kk[i], 3));
			l[i + 3] = mid;
		}
		for (auto i = 0; i < diffNum * maskNum; i++) cnt[i] = 0;
		TYPE diff_out_a, diff_out_b;
		for (int64_t exper = 0; exper < datasize; exper++) //for (int x = 0; x < range; x++)  for (int y = 0; y < range; y++)
		{
			TYPE input_a = (dis(gen)) & 0xffffffff;
			TYPE input_b = (dis(gen)) & 0xffffffff;
			TYPE output_a = 0x0, output_b = 0x0;
			speck_64(input_a, input_b, output_a, output_b, start, end, kk);
			for (auto j = 0; j < diffNum; j++)
			{
				TYPE inDiff[2] = { 0 };
				uint8_t a = inDiffIdx[j] / 32, b = inDiffIdx[j] % 32;
				inDiff[a] = ROL32(TYPE(1), b); // ((TYPE)1 << b);
				TYPE input_a_prime = (input_a) ^ inDiff[0];
				TYPE input_b_prime = (input_b) ^ inDiff[1];
				TYPE output_a_prime = 0x0, output_b_prime = 0x0;
				speck_64(input_a_prime, input_b_prime, output_a_prime, output_b_prime, start, end, kk);
				TYPE outDiff[2] = { 0x0 };
				outDiff[0] = output_a ^ output_a_prime;
				outDiff[1] = output_b ^ output_b_prime;
				for (auto i = 0; i < maskNum; i++)
				{
					uint8_t a = outMaskIdx[i] / 32, b = outMaskIdx[i] % 32;
					bool bit_out = (outDiff[a] >> b) & 0x1;
					if (bit_out == 0) tot[j * maskNum + i]++; else tot[j * maskNum + i]--;
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
		for (auto i = 0; i < maskNum; i++)
		{
			double cor = tot[j * maskNum + i] / datasize / testNum;
			if (cor)
			{
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
	}
	delete[] tot;
	float end_time = omp_get_wtime();
	std::printf("Multi time: %f\n", end_time - start_time);
	fout << "Multi time:" << end_time - start_time << endl;
	fout.close();


}

void main()
{
#if 0 //verify our new 5DL for 11-round
	int start = 3, end = 8;
	int round = end - start;
	uint8_t oneInDiffIdx[1] = { 29 + 32 };
	uint8_t oneOutMaskIdx[] = { 13 + 32 };
	string fileName = "Speck64 " + to_string(round) + "DL, verify [][29] to [][13].txt";
	getSpeck64ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 32), 100);
#endif

#if 0 //verify our new 6DL for 12-round
	int start = 3, end = 9;
	int round = end - start;
	uint8_t oneInDiffIdx[1] = { 29 + 32 };
	uint8_t oneOutMaskIdx[] = { 5 + 32 };
	string fileName = "Speck64 " + to_string(round) + "DL, verify [][29] to [][5].txt";
	getSpeck64ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 32), 100);
#endif

#if 1 //verify our new 7DL for 13-round
	int start = 3, end = 10;
	int round = end - start;
	uint8_t oneInDiffIdx[1] = { 15 };
	uint8_t oneOutMaskIdx[] = { 5 + 32 };
	string fileName = "Speck64 " + to_string(round) + "DL, [15][]to[][5].txt";
	getSpeck64ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 34), 100);
#endif

#if 0 //verify 5DL for 11-round by Lv 
	int start = 3, end = 8;
	int round = end - start;
	uint8_t oneInDiffIdx[1] = { 7 };
	uint8_t oneOutMaskIdx[] = { 5 + 32 };
	string fileName = "Speck64 " + to_string(round) + "DL, verify [7][] to [][5] by Lv.txt";
	getSpeck64ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 30), 50);
#endif

#if 0 //verify 5DL for 12-round by Lv 
	int start = 4, end = 9;
	int round = end - start;
	uint8_t oneInDiffIdx[1] = { 21 + 32 };
	uint8_t oneOutMaskIdx[] = { 5 + 32 };
	string fileName = "Speck64 " + to_string(round) + "DL, verify [][21] to [][5] by Lv.txt";
	getSpeck64ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 30), 50);
#endif


#if 0 //prepare SUBs sets for 6DL
	int start = 3, end = 9;
	int round = end - start;
	uint8_t inDiffIdx[64] = { 0 }; for (auto i = 0; i < 32; i++) { inDiffIdx[i] = i + 32; inDiffIdx[i + 32] = i; }
	uint8_t outMaskIdx[64] = { 0 }; for (auto i = 0; i < 32; i++) { outMaskIdx[i] = i + 32; outMaskIdx[i + 32] = i; }
	string fileName = "Speck64 " + to_string(round) + "DL SUBs sets.txt";
	getSpeck64SubsSetDiffWt1(fileName, inDiffIdx, 64, outMaskIdx, 64, start, end, -8, pow(2, 29), 10);
#endif

}