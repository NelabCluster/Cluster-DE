#pragma once

#include "Common.h"

enum DE_Mut_Mode{
	DERAND1,
	DERAND2,
	DEBEST1,
	DEBEST2,
	DERANDTOBEST1
};

struct DE_Individual
{
	double energy;
	double *cood;

	static DE_Individual* Individual(int N);
	void Copy(DE_Individual *de,int N);
};

class DE
{
public:
	DE(void);
	DE(DE_Mut_Mode mode, int N);
	~DE(void);
	int N;
	int popSize;
	DE_Mut_Mode mutMode;
	DE_Individual *X;
	DE_Individual bestPop;
	DE *next;
	PEnergy P_E;
	PForce P_F;

	void OneGeneration();
	void RecIndividual(DE_Individual* ind);
	double getFintessRatio();

private:
	DE_Individual *V;
	DE_Individual *U;
	void Initialization();
	void Mutation();
	void CrossOver();
	void Greedy();
	void Adjustment();
};

