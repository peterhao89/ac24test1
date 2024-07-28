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

typedef uint64_t TYPE;
TYPE p = 0x100000000;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0, p - 1);
uint64_t rand_64_bit()
{
	TYPE left = dis(gen) & 0xffffffff;
	TYPE right = dis(gen) & 0xffffffff;
	return ((left << 32) | right);
}

#define ROL64(x, s) (((x) << (s) & 0xffffffffffffffff) | ((x) >> (64 - (s))))
#define ROR64(x, s) (((x) >> (s)) | ((x) << (64 - (s)) & 0xffffffffffffffff))

void speck_128(uint64_t a_in, uint64_t b_in, uint64_t& a_out, uint64_t& b_out, int round, uint64_t kk[])
{
	uint64_t a = a_in; uint64_t b = b_in;
	for (int i = 0; i < round; i++)
	{
		a = ROR64(a, 8);
		a = (a + b) & 0xffffffffffffffff;
		a = a ^ kk[i];
		b = ROL64(b, 3) ^ a;
	}
	a_out = a;
	b_out = b;
}

void getSpeck128ExpDlCorGiven1Diff1MaskIdx(string fileName, uint8_t oneInDiffIdx[], int diffLen, uint8_t oneOutMaskIdx[], int maskLen, int round, uint64_t datasize, int testNum)
{
	TYPE master_key[4], mid;
	TYPE* kk = new TYPE[round + 1]();
	TYPE* l = new TYPE[round + 3]();
	float start_time = omp_get_wtime();
	double tot = 0;
	ofstream fout(fileName, ios::out);
	int useCoreNumber = thread::hardware_concurrency() / 2; cout << useCoreNumber << endl;
	omp_set_num_threads(useCoreNumber);
	for (int r = 0; r < testNum; r++)
	{
		for (int i = 0; i < 4; i++)
		{
			master_key[i] = rand_64_bit() & 0xffffffffffffffff;
		}
		kk[0] = master_key[0];
		for (int i = 0; i < 3; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < round; i++)
		{
			mid = ((ROR64(l[i], 8) + kk[i]) & 0xffffffffffffffff) ^ (i & 0xffffffffffffffff);
			kk[i + 1] = mid ^ (ROL64(kk[i], 3));
			l[i + 3] = mid;
		}
		double cnt = 0;
#pragma omp parallel for reduction(+:cnt,tot)
		for (int64_t expr = 0; expr < datasize; expr++)
		{
			TYPE input[2] = { 0 };
			input[0] = rand_64_bit() & 0xffffffffffffffff;
			input[1] = rand_64_bit() & 0xffffffffffffffff;
			TYPE output_a = 0x0, output_b = 0x0;
			speck_128(input[0], input[1], output_a, output_b, round, kk);
			TYPE inDiff[2] = { 0 }, input_prime[2] = { 0 };
			input_prime[0] = input[0];
			input_prime[1] = input[1];
			for (auto j = 0; j < diffLen; j++)
			{
				int8_t a = oneInDiffIdx[j] / 64, b = oneInDiffIdx[j] % 64;
				input_prime[a] ^= ((TYPE)1 << b);
			}
			TYPE output_a_prime = 0x0, output_b_prime = 0x0;
			speck_128(input_prime[0], input_prime[1], output_a_prime, output_b_prime, round, kk);
			TYPE outDiff[2] = { 0x0 };
			outDiff[0] = output_a ^ output_a_prime;
			outDiff[1] = output_b ^ output_b_prime;
			bool bit_out = 0;
			for (auto i = 0; i < maskLen; i++)
			{
				uint8_t a = oneOutMaskIdx[i] / 64, b = oneOutMaskIdx[i] % 64;
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
		if (cnt)
		{
			cnt = cnt / datasize;
			double cor = log(abs(cnt)) / log(2);
			//if (cor >= corWt)
			{
				std::printf("Speck128 %dDL experimental COR: [", round);
				fout << "Speck128 " << dec << round << "DL experimental COR : [";
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
				cout << cnt << ", " << cor << endl; fout << cnt << ", " << cor << endl;
			}
		}
	}
	if (tot)
	{
		tot = tot / datasize / testNum;
		double cor = log(abs(tot)) / log(2);
		//if (cor >= corWt)
		{
			std::printf("Speck128 %dDL AVERAGE experimental COR: [", round);
			fout << "Speck128 " << dec << round << "DL AVERAGE experimental COR : [";
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
			cout << tot << ", " << cor << endl; fout << tot << ", " << cor << endl;
		}
	}
	float end_time = omp_get_wtime();
	std::printf("Multi time: %f\n", end_time - start_time);
	fout << "Multi time:" << end_time - start_time << endl;
	fout.close();
}

void getSpeck128SubsSetDiffWt1(string fileName, uint8_t inDiffIdx[], int diffNum, uint8_t outMaskIdx[], int maskNum, int round, double corWt, uint64_t datasize, int testNum)
{
	TYPE master_key[4], mid;
	TYPE* kk = new TYPE[round + 1]();
	TYPE* l = new TYPE[round + 3]();
	float start_time = omp_get_wtime();
	double* tot = new double[diffNum * maskNum]();
	ofstream fout(fileName, ios::out);
	fout << "Threshold of the absolute correlation weight:  " << -1 * corWt << endl;
	for (int r = 0; r < testNum; r++)
	{
		for (int i = 0; i < 4; i++)
		{
			master_key[i] = rand_64_bit() & 0xffffffffffffffff;
		}  //cout << r << ": the master key is "; std::printf("(%I64x,%I64x,%I64x,%I64x)\n", master_key[0], master_key[1], master_key[2], master_key[3]);
		kk[0] = master_key[0];
		for (int i = 0; i < 3; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < round; i++)
		{
			mid = ((ROR64(l[i], 8) + kk[i]) & 0xffffffffffffffff) ^ (i & 0xffffffffffffffff);
			kk[i + 1] = mid ^ (ROL64(kk[i], 3));
			l[i + 3] = mid;
		}
		for (int64_t exper = 0; exper < datasize; exper++) //for (int x = 0; x < range; x++)  for (int y = 0; y < range; y++)
		{
			TYPE input_a = rand_64_bit() & 0xffffffffffffffff;
			TYPE input_b = rand_64_bit() & 0xffffffffffffffff;
			TYPE output_a = 0x0, output_b = 0x0;
			speck_128(input_a, input_b, output_a, output_b, round, kk);
			for (auto j = 0; j < diffNum; j++)
			{
				TYPE inDiff[2] = { 0 };
				uint8_t a = inDiffIdx[j] / 64, b = inDiffIdx[j] % 64;
				inDiff[a] = ROL64(TYPE(1), b); // ((TYPE)1 << b);
				TYPE input_a_prime = (input_a) ^ inDiff[0];
				TYPE input_b_prime = (input_b) ^ inDiff[1];
				TYPE output_a_prime = 0x0, output_b_prime = 0x0;
				speck_128(input_a_prime, input_b_prime, output_a_prime, output_b_prime, round, kk);
				TYPE outDiff[2] = { 0x0 };
				outDiff[0] = output_a ^ output_a_prime;
				outDiff[1] = output_b ^ output_b_prime;
				for (auto i = 0; i < maskNum; i++)
				{
					uint8_t a = outMaskIdx[i] / 64, b = outMaskIdx[i] % 64;
					bool bit_out = (outDiff[a] >> b) & 0x1;
					if (bit_out == 0)
					{
						tot[j * maskNum + i]++;
					}
					else
					{
						tot[j * maskNum + i]--;
					}
				}
			}
		}
	}
	for (auto j = 0; j < diffNum; j++)
	{
		int totStrongNum = 0;
		int a = inDiffIdx[j];
		fout << "bit index: ";
		if (a < 64) { fout << a + 64 << endl; }
		else { fout << a - 64 << endl; }
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
					if (b < 64) { fout << b + 64 << "(" << cor << ")  "; }
					else { fout << b - 64 << "(" << cor << ")  "; }
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

#if 0 //prepare SUBs set for 8DL
	int round = 8;
	uint8_t inDiffIdx[128] = { 0 }; for (auto i = 0; i < 64; i++) { inDiffIdx[i] = i + 64; inDiffIdx[i + 64] = i; }
	uint8_t outMaskIdx[128] = { 0 };
	for (auto i = 0; i < 64; i++) { outMaskIdx[i] = i + 64; outMaskIdx[i + 64] = i; }
	string fileName = "Speck128 " + to_string(round) + "DL SUBs sets.txt";
	getSpeck128SubsSetDiffWt1(fileName, inDiffIdx, 128, outMaskIdx, 128, round, -8, pow(2, 28), 10);
#endif

#if 0 //verify our new 9DL for 18-round
	int round = 9;
	uint8_t oneInDiffIdx[1] = { 50 };
	uint8_t oneOutMaskIdx[] = { 2 + 64 };
	string fileName = "Speck128 " + to_string(round) + "DL, verify [50][] to [][2].txt";
	getSpeck128ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, round, pow(2, 34), 50);
#endif

#if 1//verify AC23 paper 8DL, convert to 9DL 
	int round = 9;
	uint8_t oneInDiffIdx[1] = { 53 };
	uint8_t oneOutMaskIdx[] = { 5 + 64 };
	string fileName = "Speck128 " + to_string(round) + "DL, verify [53][] to [][5].txt";
	getSpeck128ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, round, pow(2, 34), 50);
#endif
}