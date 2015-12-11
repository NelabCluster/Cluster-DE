#pragma once
#include "PE_Tool.h"
#include "Common.h"
//采用单例模式

enum LocalOptimizeMethod{
	LocalOptimizeMethod_SD,
	LocalOptimizeMethod_LBFGS,
	LocalOptimizeMethod_ALL
};

class LocalTool
{
public:
	LocalTool(void);
	~LocalTool(void);
	
	PEnergy _PE_E;
	PForce _PE_F;
	double *_cood;
	double *_R;

	static LocalTool* shareDefault();
    double localWithCoodAndPE(double *cood,int N, PEnergy PE_E, PForce PE_F);
	
	void setMethod(LocalOptimizeMethod m);
	static void setPE(PotentialEnergyType type);
private:
	LocalOptimizeMethod method;


};

