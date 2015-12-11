#pragma once
#include "Common.h"

class Tool
{
public:
	Tool(void);
	~Tool(void);
	static void Distance(double *cood, double *R, int N);
	static int* RandPerm(int N, int K);
	static double Local(double *cood, int N, PEnergy PE_E, PForce PE_F);
	static double LocalWithR(double *R, int N, PEnergy PE_E, PForce PE_F);
	static void SphereCutSplice(double *fr, double *mr, double *child1, double *child2, int natom);
	static double LJ_E(double *cood);
	static void printCood(double *cood, int N, char *root);
	static void readCood(double *cood, int N, char *root);
	static void printfDiamond(double *cood, int N, double E,char *path);
	static bool readDiamond(double *cood, int N,int differ = 0);
	static void coodToDiamond(char *root, int N, char *path, PEnergy energy);

	static void testR(int N);
	static bool testThree(double *R, int N);
};

