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
uint64_t p = 0x100000000;
//uint64_t p = 0xffffffff;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0, p - 1);
uint64_t rand_48_bit()
{
	TYPE left = dis(gen) & 0xffffff;
	TYPE right = dis(gen) & 0xffffff;
	return ((left << 24) | right) & 0xffffffffffff;
}
#define ROL48(x, s) (((x) << (s) & 0xffffffffffff) | ((x) >> (48 - (s))))
#define ROR48(x, s) (((x) >> (s)) | ((x) << (48 - (s)) & 0xffffffffffff))

void speck_96(uint64_t a_in, uint64_t b_in, uint64_t& a_out, uint64_t& b_out, int start, int end, uint64_t kk[])
{
	//int alpha = 8, beta = 3;
	uint64_t a = a_in; uint64_t b = b_in;
	for (int i = start; i < end; i++)
	{
		a = ROR48(a, 8);
		a = (a + b) & 0xffffffffffff;
		a = a ^ kk[i];
		b = ROL48(b, 3) ^ a;
	}
	a_out = a;
	b_out = b;
}

void getSpeck96ExpDlCorGiven1Diff1MaskIdx(string fileName, uint8_t oneInDiffIdx[], int diffLen, uint8_t oneOutMaskIdx[], int maskLen, int start, int end, uint64_t datasize, int testNum)
{
	TYPE master_key[3], mid;
	TYPE* kk = new TYPE[end + 1]();
	TYPE* l = new TYPE[end + 2]();
	float start_time = omp_get_wtime();
	double tot = 0;
	ofstream fout(fileName, ios::out);
	int useCoreNumber = thread::hardware_concurrency() / 2;
	omp_set_num_threads(useCoreNumber);
	for (int r = 0; r < testNum; r++)
	{
		for (int i = 0; i < 3; i++)
		{
			master_key[i] = rand_48_bit() & 0xffffffffffff;
		} //cout << r << ": the master key is "; std::printf("(%I64x,%I64x,%I64x,%I64x)\n", master_key[0], master_key[1], master_key[2], master_key[3]);
		kk[0] = master_key[0];
		for (int i = 0; i < 2; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < end; i++)
		{
			mid = ((ROR48(l[i], 8) + kk[i]) & 0xffffffffffff) ^ (i & 0xffffffffffff);
			kk[i + 1] = mid ^ (ROL48(kk[i], 3));
			l[i + 2] = mid;
		}
		double cnt = 0;
#pragma omp parallel for reduction(+:cnt,tot)
		for (int64_t exper = 0; exper < datasize; exper++)
		{
			TYPE input[2] = { 0 };
			input[0] = rand_48_bit() & 0xffffffffffff;
			input[1] = rand_48_bit() & 0xffffffffffff;
			TYPE output_a = 0x0, output_b = 0x0;
			speck_96(input[0], input[1], output_a, output_b, start, end, kk);
			TYPE inDiff[2] = { 0 }, input_prime[2] = { 0 };
			input_prime[0] = input[0];
			input_prime[1] = input[1];
			for (auto j = 0; j < diffLen; j++)
			{
				int8_t a = oneInDiffIdx[j] / 48, b = oneInDiffIdx[j] % 48;
				input_prime[a] ^= ((TYPE)1 << b);
			}
			TYPE output_a_prime = 0x0, output_b_prime = 0x0;
			speck_96(input_prime[0], input_prime[1], output_a_prime, output_b_prime, start, end, kk);
			TYPE outDiff[2] = { 0x0 };
			outDiff[0] = output_a ^ output_a_prime;
			outDiff[1] = output_b ^ output_b_prime;
			bool bit_out = 0;
			for (auto i = 0; i < maskLen; i++)
			{
				uint8_t a = oneOutMaskIdx[i] / 48, b = oneOutMaskIdx[i] % 48;
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
		int round = end - start;
		double cntD = cnt / datasize;
		double cor = log(abs(cntD)) / log(2);
		std::printf("Speck96 %dDL experimental COR: [", round);
		fout << "Speck96 " << dec << round << "DL experimental COR : [";
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
		cout << cntD << "," << cor << endl; fout << cntD << "," << cor << endl;
	}
	int round = end - start;
	if (tot)
	{
		double totD = tot / datasize / testNum;
		double cor = log(abs(totD)) / log(2);
		std::printf("Speck96 %dDL AVERAGE experimental COR: [", round);
		fout << "Speck96 " << dec << round << "DL AVERAGE experimental COR : [";
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
		cout << totD << "," << cor << endl; fout << totD << "," << cor << endl;
	}
	float end_time = omp_get_wtime();
	std::printf("Multi time: %f\n", end_time - start_time);
	fout << "Multi time:" << end_time - start_time << endl;
	fout.close();
}

void getSpeck96SubsSetDiffWt1(string fileName, uint8_t inDiffIdx[], int diffNum, uint8_t outMaskIdx[], int maskNum, int start, int end, double corWt, uint64_t datasize, int testNum)
{
	TYPE master_key[3], mid;
	TYPE* kk = new TYPE[end + 1]();
	TYPE* l = new TYPE[end + 2]();
	float start_time = omp_get_wtime();
	double* tot = new double[diffNum * maskNum]();
	ofstream fout(fileName, ios::out);
	fout << "Threshold of the absolute correlation weight:  " << -1 * corWt << endl;
	for (int r = 0; r < testNum; r++)
	{
		for (int i = 0; i < 3; i++)
		{
			master_key[i] = rand_48_bit() & 0xffffffffffff;
		}
		kk[0] = master_key[0];
		for (int i = 0; i < 2; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < end; i++)
		{
			mid = ((ROR48(l[i], 8) + kk[i]) & 0xffffffffffff) ^ (i & 0xffffffffffff);
			kk[i + 1] = mid ^ (ROL48(kk[i], 3));
			l[i + 2] = mid;
		}
		for (int64_t exper = 0; exper < datasize; exper++)
		{
			TYPE input_a = rand_48_bit() & 0xffffffffffff;
			TYPE input_b = rand_48_bit() & 0xffffffffffff;
			TYPE output_a = 0x0, output_b = 0x0;
			speck_96(input_a, input_b, output_a, output_b, start, end, kk);
			for (auto j = 0; j < diffNum; j++)
			{
				TYPE inDiff[2] = { 0 };
				uint8_t a = inDiffIdx[j] / 48, b = inDiffIdx[j] % 48;
				inDiff[a] = ROL48(TYPE(1), b); // ((TYPE)1 << b);
				TYPE input_a_prime = (input_a) ^ inDiff[0];
				TYPE input_b_prime = (input_b) ^ inDiff[1];
				TYPE output_a_prime = 0x0, output_b_prime = 0x0;
				speck_96(input_a_prime, input_b_prime, output_a_prime, output_b_prime, start, end, kk);
				TYPE outDiff[2] = { 0x0 };
				outDiff[0] = output_a ^ output_a_prime;
				outDiff[1] = output_b ^ output_b_prime;
				for (auto i = 0; i < maskNum; i++)
				{
					uint8_t a = outMaskIdx[i] / 48, b = outMaskIdx[i] % 48;
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
		if (a < 48) { fout << a + 48 << endl; }
		else { fout << a - 48 << endl; }
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
					if (b < 48) { fout << b + 48 << "(" << cor << ")  "; }
					else { fout << b - 48 << "(" << cor << ")  "; }
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

void getSpeck96ExpLcCorGivenDiffMask(string fileName, TYPE inMask[2], TYPE outMask[2], int start, int end, uint64_t datasize, int testNum)
{
	TYPE master_key[3], mid;
	TYPE* kk = new TYPE[end + 1]();
	TYPE* l = new TYPE[end + 2]();
	float start_time = omp_get_wtime();
	double tot = 0;
	ofstream fout(fileName, ios::out);
	int useCoreNumber = thread::hardware_concurrency() / 2;
	omp_set_num_threads(useCoreNumber);
	for (int r = 0; r < testNum; r++)
	{
		for (int i = 0; i < 3; i++)
		{
			master_key[i] = rand_48_bit() & 0xffffffffffff;
		}
		kk[0] = master_key[0];
		for (int i = 0; i < 2; i++)
		{
			l[i] = master_key[i + 1];
		}
		for (int i = 0; i < end; i++)
		{
			mid = ((ROR48(l[i], 8) + kk[i]) & 0xffffffffffff) ^ (i & 0xffffffffffff);
			kk[i + 1] = mid ^ (ROL48(kk[i], 3));
			l[i + 2] = mid;
		}
		double cnt = 0;
#pragma omp parallel for reduction(+:cnt,tot)
		for (int64_t exper = 0; exper < datasize; exper++) //for (int x = 0; x < range; x++)  for (int y = 0; y < range; y++)
		{
			TYPE input[2] = { 0 };
			input[0] = rand_48_bit() & 0xffffffffffff;
			input[1] = rand_48_bit() & 0xffffffffffff;
			TYPE output[2] = { 0 };
			speck_96(input[0], input[1], output[0], output[1], start, end, kk);
			bool bit0 = maskdot(inMask[0], input[0], 48);
			bool bit1 = maskdot(inMask[1], input[1], 48);
			bool bit2 = maskdot(outMask[0], output[0], 48);
			bool bit3 = maskdot(outMask[1], output[1], 48);
			bool bit_out = bit0 ^ bit1 ^ bit2 ^ bit3;
			if (bit_out == 0) cnt++;
			else cnt--;
		}
		tot += abs(cnt);
		cnt = cnt / datasize;
		double cor = log(abs(cnt)) / log(2);
		std::printf("Speck96 %dLC experimental COR: ", end - start);
		fout << "Speck96 " << dec << (end - start) << "LC experimental COR : ";
		printf("(%I64x, %I64x)-->(%I64x, %I64x), ", inMask[0], inMask[1], outMask[0], outMask[1]);
		fout << "(" << inMask[0] << "," << inMask[1] << "),-->(" << outMask[0] << "," << outMask[1] << "), ";
		cout << cnt << ", " << cor << endl; fout << cnt << ", " << cor << endl;
	}
	tot = tot / datasize / testNum;
	double cor = log(abs(tot)) / log(2);
	std::printf("Speck96 %dLC AVERAGE experimental COR: ", end - start);
	fout << "Speck96 " << dec << (end - start) << "LC AVERAGE experimental COR : ";
	printf("(%I64x, %I64x)-->(%I64x, %I64x), ", inMask[0], inMask[1], outMask[0], outMask[1]);
	fout << "(" << inMask[0] << "," << inMask[1] << "),-->(" << outMask[0] << "," << outMask[1] << "), ";
	cout << tot << ", " << cor << endl; fout << tot << ", " << cor << endl;

	float end_time = omp_get_wtime();
	std::printf("Multi time: %f\n", end_time - start_time);
	fout << "Multi time:" << end_time - start_time << endl;
	fout.close();
}

void main()
{
#if 0 //prepare SUBs sets for 6DL
	int start = 4, end = 10;
	int round = end - start;
	uint8_t inDiffIdx[96] = { 0 }; for (auto i = 0; i < 48; i++) { inDiffIdx[i] = i + 48; inDiffIdx[i + 48] = i; }
	uint8_t outMaskIdx[96] = { 0 }; for (auto i = 0; i < 48; i++) { outMaskIdx[i] = i + 48; outMaskIdx[i + 48] = i; }
	string fileName = "Speck96 " + to_string(round) + "DL SUBs sets.txt";
	getSpeck96SubsSetDiffWt1(fileName, inDiffIdx, 96, outMaskIdx, 96, start, end, -8, pow(2, 29), 10);
#endif

#if 0 //verify our new 7DL for 15-round
	int start = 4, end = 11;
	int round = end - start;
	uint8_t inDiffIdx[1] = { 43 + 48 };
	uint8_t outMaskIdx[2] = { 10 + 48, 9 + 48 };
	string fileName = "Speck96 " + to_string(round) + "DL, verify [][43] to [][10,9].txt";
	getSpeck96ExpDlCorGiven1Diff1MaskIdx(fileName, inDiffIdx, 1, outMaskIdx, 2, start, end, pow(2, 35), 100);
#endif

#if 1 //verify our 8DL for 16-round
	int start = 4, end = 12;
	int round = end - start;
	uint8_t inDiffIdx[1] = { 7 };
	uint8_t outMaskIdx[1] = { 26 + 48 };
	string fileName = "Speck96 " + to_string(round) + "DL, verify [7][] to [][26].txt";
	getSpeck96ExpDlCorGiven1Diff1MaskIdx(fileName, inDiffIdx, 1, outMaskIdx, 1, start, end, pow(2, 36), 50);
#endif



#if 0//verify AC23 paper 6DL, v0[45]->(v0[21],v1[13,12,10])
	int start = 4, end = 10;
	int round = end - start;
	uint8_t oneInDiffIdx[1] = { 45 };
	uint8_t oneOutMaskIdx[4] = { 21, 13 + 48, 12 + 48, 10 + 48 };
	string fileName = "Speck96 " + to_string(round) + "DL, verify AC [45][] to [21][13,12,10].txt";
	getSpeck96ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 4, start, end, pow(2, 29), 10);
#endif

#if 0// verify AC23 paper 6DL, convert to 7DL v0[45]->v1[13]
	int start = 4, end = 11;
	int round = end - start;
	uint8_t oneInDiffIdx[1] = { 45 };
	uint8_t oneOutMaskIdx[1] = { 13 + 48 };
	string fileName = "Speck96 " + to_string(round) + "DL, verify AC [45][] to [][13].txt";
	getSpeck96ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 33), 50);
#endif

#if 0 //verify 4LC cor, (0, 4000000)->(20000000, 20000000)->(101a00000, 101800000)->(801083000, 800090000)->(404c41cd20, 4040404100), Prob: 8
	int start = 12, end = 16;
	int round = end - start;
	TYPE inMask[2] = { 0, 0x4000000 };
	TYPE outMask[2] = { 0x404c41cd20, 0x4040404100 };
	string fileName = "Speck96 " + to_string(end - start) + "LC verify.txt";
	getSpeck96ExpLcCorGivenDiffMask(fileName, inMask, outMask, start, end, pow(2, 33), 50);
#endif

}