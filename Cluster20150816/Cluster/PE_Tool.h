#pragma once
#include "Common.h"

enum PotentialEnergyType
{
	PotentialEnergyType_FS,
	PotentialEnergyType_Johnson
};

class PE_Tool
{
public:
	PE_Tool(void);
	~PE_Tool(void);
	
	static PEnergy GetPEnergy(PotentialEnergyType type);
	static PForce GetPForce(PotentialEnergyType type);

	static double FS_E(double *R,int N);
	static double FS_F(double *cood,double *R,double *FF, int N);
	static double FS_Fr(double *cood,double *R,double *FF, int N);

	static double Johnson_E(double *R,int N);
	static double Johnson_F(double *cood, double *R, double *FF, int N);
};

