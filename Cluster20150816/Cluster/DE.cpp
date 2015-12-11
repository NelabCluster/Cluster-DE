#include "DE.h"
#include "Tool.h"
#include "PE_Tool.h"
#include "LocalTool.h"
#include <cstring>
extern string EnergyName;
DE_Individual* DE_Individual::Individual(int N)
{
	DE_Individual *obj = (DE_Individual*)malloc(sizeof(DE_Individual));
	obj->cood = (double *)malloc(3 * N * sizeof(double));
	return obj;
}

void DE_Individual::Copy(DE_Individual *de, int N)
{
	int i;
	for(i = 0; i < 3 * N; i++)
	{
		this->cood[i] = de->cood[i];
	}
	this->energy = de->energy;
}

DE::DE(void)
{
	int i;

	cout<<"N="<<N<<endl;
	this->popSize = 50;//(N+10)/10*10;

	if (EnergyName ==  "FS")
	{	
		this->P_E = PE_Tool::FS_E;
		this->P_F = PE_Tool::FS_F;
	} else
	{
		this->P_E = PE_Tool::Johnson_E;
		this->P_F = PE_Tool::Johnson_F;
	}


	this->next = NULL;
	this->mutMode = DERAND1;
	this->X = (DE_Individual *)malloc(this->popSize * sizeof(DE_Individual));
	this->V = (DE_Individual *)malloc(this->popSize * sizeof(DE_Individual));
	this->U = (DE_Individual *)malloc(this->popSize * sizeof(DE_Individual));
	this->bestPop.cood = (double *)malloc(3*this->N * sizeof(double));
	for( i = 0; i < this->popSize; i++)
	{	
		this->X[i].cood = (double *)malloc(3*this->N * sizeof(double));
		this->V[i].cood = (double *)malloc(3*this->N * sizeof(double));
		this->U[i].cood = (double *)malloc(3*this->N * sizeof(double));
	}
	this->Initialization();
}

DE::~DE(void)
{
	int i;
	free(this->bestPop.cood);
	for( i = 0; i < this->popSize; i++)
	{	
		free(this->X[i].cood);
		free(this->V[i].cood);
		free(this->U[i].cood);
	}
	free(this->X);
	free(this->V);
	free(this->U);
}

DE::DE(DE_Mut_Mode mode, int N)
{
	this->mutMode = mode;
	this->N = N;
	
	int i;

	cout<<"N="<<N<<endl;
	this->popSize = 50;//(N+10)/10*10;

	if (EnergyName ==  "FS")
	{	
		this->P_E = PE_Tool::FS_E;
		this->P_F = PE_Tool::FS_F;
	} else
	{
		this->P_E = PE_Tool::Johnson_E;
		this->P_F = PE_Tool::Johnson_F;
	}


	this->next = NULL;
	this->X = (DE_Individual *)malloc(this->popSize * sizeof(DE_Individual));
	this->V = (DE_Individual *)malloc(this->popSize * sizeof(DE_Individual));
	this->U = (DE_Individual *)malloc(this->popSize * sizeof(DE_Individual));
	this->bestPop.cood = (double *)malloc(3*this->N * sizeof(double));
	for( i = 0; i < this->popSize; i++)
	{	
		this->X[i].cood = (double *)malloc(3*this->N * sizeof(double));
		this->V[i].cood = (double *)malloc(3*this->N * sizeof(double));
		this->U[i].cood = (double *)malloc(3*this->N * sizeof(double));
	}
	this->Initialization();
}

void DE::OneGeneration()
{
	this->Mutation();
	this->CrossOver();
	this->Greedy();
}

void DE::RecIndividual(DE_Individual *ind)
{
	int r  = rand() % this->popSize;;
	this->X[r].Copy(ind,this->N);
}

double DE::getFintessRatio()
{
	double fmax = this->X[0].energy;
	double fmin = this->X[0].energy;
	for (int i = 1; i < popSize; i ++)
	{
		if (this->X[i].energy > fmax) fmax = this->X[i].energy;
		if (this->X[i].energy < fmin) fmin = this->X[i].energy;
	}

	if ( fabs(fmax - fmin) < 0.0001 ) return 1.0;

	double up = 0;
	for (int i = 0; i < popSize; i++)
		up += this->X[i].energy - fmin;

	double down = popSize * (fmax - fmin);

	return up/down;
}

void DE::Initialization()
{
	int i,j,differ;
	double r0= 2.75;
	DE_Individual *best;
	
	differ = 3;
	best = this->X;

	for (int i=0;i<=differ;i++)
	{
		if (Tool::readDiamond(X[i].cood,N,i+1))
		{
			for(j = 0; j < i+1; j++)
			{
				X[i].cood[N-j-1] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
				X[i].cood[2*N-j-1] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
				X[i].cood[3*N-j-1] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
			}
		} 
		else
		{
			for(j = 0; j < 3*N; j++)
			{
				this->X[i].cood[j] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
			}
		}
		X[i].energy = LocalTool::shareDefault()->localWithCoodAndPE(this->X[i].cood,this->N,this->P_E,this->P_F);
		if(this->X[i].energy < best->energy)
			best = this->X+i;
	}


	for( i = differ+1; i < this->popSize; i++)
	{	
		for(j = 0; j < 3*N; j++)
		{
			this->X[i].cood[j] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
		}
		this->X[i].energy = LocalTool::shareDefault()->localWithCoodAndPE(this->X[i].cood,this->N,this->P_E,this->P_F);
		if(this->X[i].energy < best->energy)
			best = this->X+i;
	}
	this->bestPop.Copy(best,this->N);
}

void DE::Mutation()
{
	int i,j;
	switch(this->mutMode)
	{
	case DERAND1:
		for(i = 0; i < this->popSize; i++)
		{
			double F = 2 * RAND;
			int *p = Tool::RandPerm(this->popSize,3);
			for(j = 0; j < 3 * this->N; j++)
			{
				this->V[i].cood[j] = this->X[p[2]].cood[j] + F * (this->X[p[0]].cood[j] - this->X[p[1]].cood[j]);
			}
			free(p);
		}
		break;
	case DEBEST1:
		for(i = 0; i < this->popSize; i++)
		{
			double F = 2 * RAND;
			int *p = Tool::RandPerm(this->popSize,2);
			for(j = 0; j < 3 * this->N; j++)
			{
				this->V[i].cood[j] = this->bestPop.cood[j] + F * (this->X[p[0]].cood[j] - this->X[p[1]].cood[j]);
			}
			free(p);
		}
		break;
	case DEBEST2:
		for(i = 0; i < this->popSize; i++)
		{
			double F = 2 * RAND;
			int *p = Tool::RandPerm(this->popSize,4);
			for(j = 0; j < 3 * this->N; j++)
			{
				this->V[i].cood[j] = this->bestPop.cood[j] + F * (this->X[p[0]].cood[j] + this->X[p[1]].cood[j] - this->X[p[2]].cood[j] - this->X[p[3]].cood[j]);
			}
			free(p);
		}
		break;
	case DERAND2:
		for(i = 0; i < this->popSize; i++)
		{
			double F = 2 * RAND;
			int *p = Tool::RandPerm(this->popSize,5);
			for(j = 0; j < 3 * this->N; j++)
			{
				this->V[i].cood[j] = this->X[p[4]].cood[j] + F * (this->X[p[0]].cood[j] + this->X[p[1]].cood[j] - this->X[p[2]].cood[j] - this->X[p[3]].cood[j]);
			}
			free(p);
		}
		break;
	case DERANDTOBEST1:
		for(i = 0; i < this->popSize; i++)
		{
			double F = 2 * RAND;
			int *p = Tool::RandPerm(this->popSize,2);
			for(j = 0; j < 3 * this->N; j++)
			{
				this->V[i].cood[j] = this->X[i].cood[j] + F * (this->X[p[0]].cood[j] + this->bestPop.cood[j] - this->X[i].cood[j] - this->X[p[1]].cood[j]);
			}
			free(p);
		}
		break;
		break;
	}
}

void DE::CrossOver()
{	
	int i;
	DE_Individual childDE1,childDE2;

	childDE1.cood = (double *)malloc(3 * this->N * sizeof(double));
	childDE2.cood = (double *)malloc(3 * this->N * sizeof(double));
	for(i = 0; i< this->popSize; i++)
	{
		Tool::SphereCutSplice(this->X[i].cood,this->V[i].cood,childDE1.cood,childDE2.cood,this->N);

		childDE1.energy = LocalTool::shareDefault()->localWithCoodAndPE(childDE1.cood,this->N,this->P_E,this->P_F);
		childDE2.energy = LocalTool::shareDefault()->localWithCoodAndPE(childDE2.cood,this->N,this->P_E,this->P_F);

		if(childDE1.energy < childDE2.energy)
		{
			this->U[i].Copy(&childDE1,this->N);
		} else
		{
			this->U[i].Copy(&childDE2,this->N);
		}
	}
	free(childDE1.cood);
	free(childDE2.cood);
	
}

void DE::Greedy()
{
	int i;
	for(i = 0; i < this->popSize; i++)
	{
		if(this->U[i].energy < this->X[i].energy){
			this->X[i].Copy(&this->U[i],this->N);
			cout<<"ÐÂ¹¹ÐÍ"<<this->X[i].energy<<endl;
		}
		if(this->X[i].energy < this->bestPop.energy)
			this->bestPop.Copy(this->X+i,this->N);
	}
}

void DE::Adjustment()
{

}