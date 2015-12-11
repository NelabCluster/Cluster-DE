#pragma once
#include "DE.h"
#include "List.h"
class DE_Manage
{
public:
	DE_Manage(void);
	~DE_Manage(void);
	void inputPara();
	void start();
	void Repeat();

	void SinglePop();
	void MultiPop_Line();
	double MultiPop_Master(int N);
	void MultiPop_Pool();
	void AddPop(DE *p);
	void removeAllPop();
	
private:
	int popNum;
	DE *pop;
	ListNode *numList;
	int iterations;
	int repeat;
	int current;
	double bestResult;
	DE_Individual *best;
	int atomNum;

	void readBestResult();

	void printfEnergy(FILE *fp, int index);
	void printfFitnessRatio(FILE *fp, int index);
};

