#pragma once
#include "Common.h"

class BH_Individual
{
public:
	int N;
	double *cood;
	double energy;

	BH_Individual();
	BH_Individual(int N);
	void setN(int N);
	~BH_Individual();
	void operator=(BH_Individual*);
};

class BasinHopping
{
public:
	BasinHopping(int N);
	~BasinHopping(void);

	void start();

private:
	int N;
	BH_Individual current;
	BH_Individual best;
	PEnergy P_E;
	PForce P_F;
	float successRate;

	void initDate();
	
};

