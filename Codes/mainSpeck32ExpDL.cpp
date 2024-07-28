#include "basics.h"
#include <fstream>
#include <string>
#include<windows.h>
#include<omp.h>
#include <iostream>
using namespace std;
//NTL_CLIENT
typedef uint16_t TYPE;

uint32_t p = 0x10000;
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dis(0, p - 1);

#define ROL16(x, s) (((x) << (s) & 0xffff) | ((x) >> (16 - (s))))
#define ROR16(x, s) (((x) >> (s)) | ((x) << (16 - (s)) & 0xffff))

void speck_32(TYPE a_in, TYPE b_in, TYPE& a_out, TYPE& b_out, int start, int end, TYPE kk[])
{
    TYPE a = a_in, b = b_in;
    for (int i = start; i < end; i++)
    {
        a = ROR16(a, 7);
        a = (a + b) & 0xffff;
        a = a ^ kk[i];
        b = ROL16(b, 2) ^ a;
    }
    a_out = a;
    b_out = b;
}

void getSpeck32ExpDlCorGiven1Diff1MaskIdx(string fileName, uint8_t oneInDiffIdx[], int diffLen, uint8_t oneOutMaskIdx[], int maskLen, int start, int end, uint64_t datasize, int testNum)
{
    int round = end - start;
    TYPE master_key[4], mid;
    TYPE* kk = new TYPE[end + 1]();
    TYPE* l = new TYPE[end + 3]();
    float start_time = omp_get_wtime();
    double tot = 0;
    ofstream fout(fileName, ios::out | ios::app);
    int useCoreNumber = thread::hardware_concurrency() / 2; cout << useCoreNumber << endl;
    omp_set_num_threads(useCoreNumber);
    for (int r = 0; r < testNum; r++)
    {
        for (int i = 0; i < 4; i++)
        {
            master_key[i] = (dis(gen)) & 0xffff;
        } 
        kk[0] = master_key[0];
        for (int i = 0; i < 3; i++)
        {
            l[i] = master_key[i + 1];
        }
        for (int i = 0; i < end; i++)
        {
            mid = ((ROR16(l[i], 7) + kk[i]) & 0xffff) ^ (i & 0xffff);
            kk[i + 1] = mid ^ (ROL16(kk[i], 2));
            l[i + 3] = mid;
        }
        double cnt = 0; 
#pragma omp parallel for reduction(+:cnt,tot)
        for (int64_t expr = 0; expr < datasize; expr++) 
        {
            TYPE input[2] = { 0 };
            input[0] = (dis(gen)) & 0xffff;
            input[1] = (dis(gen)) & 0xffff;
            TYPE output_a = 0x0, output_b = 0x0;
            speck_32(input[0], input[1], output_a, output_b, start, end, kk);
            TYPE inDiff[2] = { 0 }, input_prime[2] = { 0 };
            input_prime[0] = input[0];
            input_prime[1] = input[1];
            for (auto j = 0; j < diffLen; j++)
            {
                int8_t a = oneInDiffIdx[j] / 16, b = oneInDiffIdx[j] % 16;
                input_prime[a] ^= ((TYPE)1 << b);
            }
            TYPE output_a_prime = 0x0, output_b_prime = 0x0;
            speck_32(input_prime[0], input_prime[1], output_a_prime, output_b_prime, start, end, kk);
            TYPE outDiff[2] = { 0x0 };
            outDiff[0] = output_a ^ output_a_prime;
            outDiff[1] = output_b ^ output_b_prime;
            bool bit_out = 0;
            for (auto i = 0; i < maskLen; i++)
            {
                uint8_t a = oneOutMaskIdx[i] / 16, b = oneOutMaskIdx[i] % 16;
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
		cnt = cnt / datasize;
		double cor = log(abs(cnt)) / log(2);
		std::printf("Speck32 %dDL experimental COR: [", round);
		fout << "Speck32 " << dec << round << "DL experimental COR : [";
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
    if (tot)
    {
        tot = tot / datasize / testNum;
        double cor = log(abs(tot)) / log(2);
        //if (cor >= corWt)
        {
            std::printf("Speck32 %dDL AVERAGE experimental COR: [", round);
            fout << "Speck32 " << dec << round << "DL AVERAGE experimental COR : [";
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
    float end_Time = omp_get_wtime();
    std::printf("Multi time: %f\n", end_Time - start_time);
    fout << "Multi time:" << end_Time - start_time << endl;
    fout.close();
}

void main()
{
#if 1
    int start = 3, end = 7;
    int round = end - start;
    uint8_t oneInDiffIdx[1] = { 3 + 16 };
    uint8_t oneOutMaskIdx[1] = { 3 + 16 };
    string fileName = "Speck32 " + to_string(round) + "DL, verify [][3] to [][3].txt";
    getSpeck32ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 32), 100);
#endif

#if 0
    int start = 3, end = 7;
    int round = end - start;
    uint8_t oneInDiffIdx[1] = { 10 + 16 };
    uint8_t oneOutMaskIdx[1] = { 3 + 16 };
    string fileName = "Speck32 " + to_string(round) + "DL, verify [][10] to [][3].txt";
    getSpeck32ExpDlCorGiven1Diff1MaskIdx(fileName, oneInDiffIdx, 1, oneOutMaskIdx, 1, start, end, pow(2, 32), 100);
#endif

}