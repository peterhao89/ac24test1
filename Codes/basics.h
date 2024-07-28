#pragma once
//#include <NTL/RR.h>
#include <ostream>
#include <vector>
#include <fstream>
#include <random>
#include <string>
#include<cstdlib>
#include<cstdint>
#include<cmath>
#include<array>
//#include "../Util/BoostCommon.hpp"
#include <stdint.h>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>
#include<memory.h>
#include <smmintrin.h>
#include<thread>
//using namespace std;
//NTL_CLIENT

#define ROL16(x, s) (((x) << (s) & 0xffff) | ((x) >> (16 - (s))))
#define ROR16(x, s) (((x) >> (s)) | ((x) << (16 - (s)) & 0xffff))
#define ROL24(x, s) (((x) << (s) & 0xffffff) | ((x) >> (24 - (s))))
#define ROR24(x, s) (((x) >> (s)) | ((x) << (24 - (s)) & 0xffffff))
#define ROL32(x, s) (((x) << (s) & 0xffffffff) | ((x) >> (32 - (s))))
#define ROR32(x, s) (((x) >> (s)) | ((x) << (32 - (s)) & 0xffffffff))
#define ROL48(x, s) (((x) << (s) & 0xffffffffffff) | ((x) >> (48 - (s))))
#define ROR48(x, s) (((x) >> (s)) | ((x) << (48 - (s)) & 0xffffffffffff))
#define ROL64(x, s) (((x) << (s) & 0xffffffffffffffff) | ((x) >> (64 - (s))))
#define ROR64(x, s) (((x) >> (s)) | ((x) << (64 - (s)) & 0xffffffffffffffff))
#define bit64(x, n) (((x) >> (n)) & 1)
#define nc 4
#define wordsize 32
#define sizeAlzette 32
#define sizeSpeck32 16
#define sizeSpeck48 24
#define sizeSpeck64 32
#define sizeSpeck96 48
#define sizeSpeck128 64
extern const int rConst[8];
extern const int sConst[8];
//typedef NTL::RR RR;


uint32_t getInt16(int i);
uint32_t getInt16(int i, bool yorn);
uint32_t getInt24(int i);
uint32_t getInt24(int i, bool yorn);
uint32_t getInt32(int i);
uint32_t getInt32(int i, bool yorn);
uint64_t getInt48(int i);
uint64_t getInt48(int i, bool yorn);
uint64_t getInt64(int i);
uint64_t getInt64(int i, bool yorn);
bool maskdot(uint64_t mask, uint64_t value, int len);

void chachaHalfQR(uint32_t in[4], uint32_t out[4], int r, int l);
void chachaStartEnd(uint32_t in[16], uint32_t out[16], int start, int end);
void chachaHalfQR_avx(__m128i a, __m128i b, __m128i c, __m128i d, __m128i& aa, __m128i& bb, __m128i& cc, __m128i& dd, int r, int l);
void chachaStartR2toR6_avx(__m128i in[4], __m128i out[4]);
void chachaStartR2toR7_avx(__m128i in[4], __m128i out[4]);
void chachaStartR2toR6R7_avx(__m128i in[4], __m128i outR6[4], __m128i outR7[4]);
bool getOneLsbOfInt32FromInt128(__m128i val128, uint8_t b);
bool getLinComOfInt32FromInt128(__m128i val128, uint8_t b, uint32_t mask);
void print128_num(__m128i var);



