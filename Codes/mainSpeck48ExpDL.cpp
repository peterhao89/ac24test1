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

typedef uint32_t TYPE;
uint64_t p = 0x100000000;
//uint32_t p = 0x7fffffff;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0, p - 1);

#define ROL24(x, s) (((x) << (s) & 0xffffff) | ((x) >> (24 - (s))))
#define ROR24(x, s) (((x) >> (s)) | ((x) << (24 - (s)) & 0xffffff))

void speck_48(uint32_t a_in, uint32_t b_in, uint32_t& a_out, uint32_t& b_out, int start, int end, uint32_t kk[])
{
	uint32_t a = a_in, b = b_in;
	for (int i = start; i < end; i++)
	{
		a = ROR24(a, 8);
		a = (a + b) & 0xffffff;
		a = a ^ kk[i];
		b = ROL24(b, 3) ^ a;
	}
	a_out = a;
	b_out = b;
}

void getSpeck48SubsSetDiffWt1(string fileName, uint8_t inDiffIdx[], int diffNum, uint8_t outMaskIdx[], int maskNum, int start, int end, double corWt, uint64_t datasize, int testNum)
{
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
			master_key[i] = (dis(gen)) & 0xffffff;
		}
		kk[0] = master_key[0];
		for (int i = 0; i < 3; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < end; i++)
		{
			mid = ((ROR24(l[i], 8) + kk[i]) & 0xffffff) ^ (i & 0xffffff);
			kk[i + 1] = mid ^ (ROL24(kk[i], 3));
			l[i + 3] = mid;
		}
		for (auto i = 0; i < diffNum * maskNum; i++) cnt[i] = 0;
		for (int64_t expr = 0; expr < datasize; expr++)
		{
			TYPE input_a = (dis(gen)) & 0xffffff;
			TYPE input_b = (dis(gen)) & 0xffffff;
			TYPE output_a = 0x0, output_b = 0x0;
			speck_48(input_a, input_b, output_a, output_b, start, end, kk);
			for (auto j = 0; j < diffNum; j++)
			{
				TYPE inDiff[2] = { 0 };
				uint8_t a = inDiffIdx[j] / 24, b = inDiffIdx[j] % 24;
				inDiff[a] = ROL24(1, b);
				TYPE input_a_prime = (input_a) ^ inDiff[0];
				TYPE input_b_prime = (input_b) ^ inDiff[1];
				TYPE output_a_prime = 0x0, output_b_prime = 0x0;
				speck_48(input_a_prime, input_b_prime, output_a_prime, output_b_prime, start, end, kk);
				TYPE outDiff[2] = { 0x0 };
				outDiff[0] = output_a ^ output_a_prime;
				outDiff[1] = output_b ^ output_b_prime;
				for (auto i = 0; i < maskNum; i++)
				{
					uint8_t a = outMaskIdx[i] / 24, b = outMaskIdx[i] % 24;
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
		if (a < 24) { fout << a + 24 << endl; }
		else { fout << a - 24 << endl; }
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
					if (b < 24) { fout << b + 24 << "(" << cor << ")  "; }
					else { fout << b - 24 << "(" << cor << ")  "; }
				}
			}
		}
		fout << endl;
		fout << "the number of strong unbalanced bit is " << totStrongNum << endl;
	}
	delete[] tot;
	delete[] cnt;
	float end_time = omp_get_wtime();
	std::printf("Multi time: %f\n", end_time - start_time);
	fout << "Multi time:" << end_time - start_time << endl;
	fout.close();
}

void getSpeck48ExpDlCorGiven1Diff1MaskIdx(string fileName, uint8_t oneInDiffIdx[], int diffLen, uint8_t oneOutMaskIdx[], int maskLen, int start, int end, uint64_t datasize, int testNum)
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
			master_key[i] = (dis(gen)) & 0xffffff;
		}
		kk[0] = master_key[0];
		for (int i = 0; i < 3; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < end; i++)
		{
			mid = ((ROR24(l[i], 8) + kk[i]) & 0xffffff) ^ (i & 0xffffff);
			kk[i + 1] = mid ^ (ROL24(kk[i], 3));
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
			speck_48(input[0], input[1], output_a, output_b, start, end, kk);
			TYPE inDiff[2] = { 0 }, input_prime[2] = { 0 };
			input_prime[0] = input[0];
			input_prime[1] = input[1];
			for (auto j = 0; j < diffLen; j++)
			{
				int8_t a = oneInDiffIdx[j] / 24, b = oneInDiffIdx[j] % 24;
				input_prime[a] ^= ((TYPE)1 << b);
			}
			TYPE output_a_prime = 0x0, output_b_prime = 0x0;
			speck_48(input_prime[0], input_prime[1], output_a_prime, output_b_prime, start, end, kk);
			TYPE outDiff[2] = { 0x0 };
			outDiff[0] = output_a ^ output_a_prime;
			outDiff[1] = output_b ^ output_b_prime;
			bool bit_out = 0;
			for (auto i = 0; i < maskLen; i++)
			{
				uint8_t a = oneOutMaskIdx[i] / 24, b = oneOutMaskIdx[i] % 24;
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
		std::printf("Speck48 %dDL experimental COR: [", round);
		fout << "Speck48 " << dec << round << "DL experimental COR : [";
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
		std::printf("Speck48 %dDL AVERAGE experimental COR: [", round);
		fout << "Speck48 " << dec << round << "DL AVERAGE experimental COR : [";
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


int main()
{
#if 0 //prepare SUBs set for 5DL
	int start = 2, end = 7;
	int round = end - start;
	uint8_t inDiffIdx[48] = { 0 }; for (auto i = 0; i < 24; i++) { inDiffIdx[i] = i + 24; inDiffIdx[i + 24] = i; }
	uint8_t outMaskIdx[48] = { 0 }; for (auto i = 0; i < 24; i++) { outMaskIdx[i] = i + 24; outMaskIdx[i + 24] = i; }
	string fileName = "Speck48 " + to_string(round) + "DL SUBs sets.txt";
	getSpeck48SubsSetDiffWt1(fileName, inDiffIdx, 48, outMaskIdx, 48, start, end, -6.5, pow(2, 30), 10);
#endif

#if 1 //verify 6DL for 11-round in both [30] and [41] 
	int start = 2, end = 8;
	int round = end - start;
	uint8_t oneInDiffIdx[1] = { 15 };
	uint8_t oneOutMaskIdx[1] = { 5 + 24 };
	string fileName = "Speck48 " + to_string(round) + "DL, verify [15][] to [][5].txt";
	getSpeck48ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 32), 50);
#endif

}