#include "basics.h"
#include <cmath>
#include <fstream>
//NTL_CLIENT

using namespace std;

uint32_t getInt16(int i)
{
	return ((uint32_t)1 << (i + 16) % 16);
}
uint32_t getInt16(int i, bool yorn)
{
	if (yorn) return ((uint32_t)1 << (i + 16) % 16);
	else return 0;
}
uint32_t getInt24(int i)
{
	return ((uint32_t)1 << (i + 24) % 24);
}
uint32_t getInt24(int i, bool yorn)
{
	if (yorn) return ((uint32_t)1 << (i + 24) % 24);
	else return 0;
}
uint32_t getInt32(int i)
{
	return ((uint32_t)1 << (i + 32) % 32);
}
uint32_t getInt32(int i, bool yorn)
{
	if (yorn) return ((uint32_t)1 << (i + 32) % 32);
	else return 0;
}
uint64_t getInt48(int i)
{
	return ((uint64_t)1 << (i + 48) % 48);
}
uint64_t getInt48(int i, bool yorn)
{
	if (yorn) return ((uint64_t)1 << (i + 48) % 48);
	else return 0;
}
uint64_t getInt64(int i)
{
	return ((uint64_t)1 << (i + 64) % 64);
}
uint64_t getInt64(int i, bool yorn)
{
	if (yorn) return ((uint64_t)1 << (i + 64) % 64);
	else return 0;
}

bool maskdot(uint64_t mask, uint64_t value, int len)
{
	bool tmp = 0;
	for (int i = 0; i < len; i++)
	{
		if (((mask >> i) & 0x1)) tmp ^= ((value >> i) & 0x1);
	}
	return tmp;
}