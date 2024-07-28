# Supplementary Materials for the Asiacrypt 2024 Submission 124
 



The source codes for verifying our results are given in the `Codes` folder. 
The output results after tests are given in the `Results` folder, where 11 subfolders named ('targetCipher'+results) are included.  

For example, we provide all the results involved in Section 5.1 in folders (‘Speck32 Results’), (‘Speck48 Results’), (‘Speck64 Results’), (‘Speck96 Results’), (‘Speck128 Results’), respectively, where the SUBs sets, the experimental DL correlations, etc. are included. Besides, we also provide the theoretical DC trails ending with 1-bit output differences and the LC trails starting with 1-bit or consecutive 2-bit input masks derived from the MILP-based approach during our constructions of the hourglass(-like) structural DL distinguishers. Moreover, some intermediate results are also given. 

The source codes are accelerated with AVX2 instructions and can be compiled by running the following commands:
```
cd Codes
mkdir build
cd build
cmake ..
make
```
Then, we can get the executable files and verify the results of arbitrary `targetCipher` by running the file named as `test{targetCipher}`.  
For example, we can verify our experimental DL results on Alzette by running the file `testAlzette` as
``
./testAlzette
``

Our verfications are fully implemented on a workstation with two CPUs (Intel(R) Xeon(R) Silver 4210 CPU @2.20GHz 2.19 GHz) and 64.0GB RAM. The operation system is Windows 10 and the codes are compiled with Visual Studio 2022 accelerated with openmp utilizing 20 threads.

Some of the testing information and results on the Speck family are as follows:

|Cipher |DL approximation |Sample Size|   Key Number|Running Time(sec)|Ref.
|----|----|----|----|----|----|
|Speck32 |   4DL: v1[3]->v1[3]|       2^32  |         100  |   22832  |     Section 5.1|
|Speck32 |   4DL: v1[10]->v1[3] |     2^32 |          100   |        23232 |        [8]
|Speck48 |   6DL: v0[15]->v1[5]|     2^32 |           50 |          8576    |    Supp. Mat. F.2|
|Speck64  |  5DL: v1[29]->v1[13]   |  2^32  |         100 |          20864  |     Section 5.1|
|Speck64   | 6DL: v1[29]->v1[5]   |   2^32   |        100    |       15488    |   Section 5.1|
|Speck64  |  7DL: v0[15]->v1[5]      |2^34    |       100     |      57312    |   Section 5.1|
|Speck64 |   5DL: v0[7]->v1[5]      | 2^30     |       50    |       2698     |     [30]|
|Speck64  |  5DL: v1[21]->v1[5]      |2^30   |         50     |      2548     |     [30]|
|Speck96   | 7DL: v0[43]->v1[10,9]   |2^35       |    100          | 91731   |    Section 5.1|
|Speck96    |8DL: v0[7]->v1[26]      |2^36      |      50    |       285116 |     Section 5.1  |
|Speck96  |  7DL: v0[45]->v1[13]  |   2^35   |        100   |        75364    |     [41]  |
|Speck128  | 9DL: v0[50]->v1[2]   |   2^34        |    50        |  89905     |  Section 5.1 |
|Speck128  | 9DL: v0[53]->v1[15]     |2^34     |       50     |      98351   |      [41] |  

 Some of the testing information and results on Alzette, SipHash, Chaskey and ChaCha are
 
|Cipher |DL approximation |Sample Size| Running Time(sec)|Ref.
|----|----|----|----|----|
|ChaCha | 2.5DL(R1->R3.5): (v1[22],v5[29,25,17,5],v9[30,18,10],v13[30,18])->v9[0]|   2^43|  1893940 | Section 6.1|
|Alzette | 6DL(R2->R8): v1[16]->v1[6,5] |     2^34  |         983  |    Section 5.2|
|Alzette | 7DL(R2->R9): v1[25]->v1[15,14] |   2^38  |        16659  |   Section 5.2|
|Alzette | 7DL(R3->R10): v0[29]->v1[1]  |     2^34  |        1074   |  Section 5.2|
|SipHash | 2.5DL(R0->R2.5): v2[63]->v1[1] |   2^25  |         5.5   |   Section 5.3|
|SipHash | 3DL(R0->R3): v2[6]->v1[19]  |      2^25  |          6   |    Section 5.3|
|Chaskey | 4DL(R1.5->R5.5): v0[31]->v3[11] |   2^34  |         1702 |   Section 5.4 |