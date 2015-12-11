#include "BasinHopping.h"
#include "LocalTool.h"
#include "Tool.h"

BH_Individual::BH_Individual(void){
	cout<<"bh-ind"<<endl;
}

BH_Individual::BH_Individual(int N){
	cout<<"bh-ind"<<endl;
	this->setN(N);
}

BH_Individual::~BH_Individual(void){

	free(this->cood);
}

void BH_Individual::setN(int N)
{
	this->N = N;
	this->cood = (double *)malloc(3 * N * sizeof(double));
}

void BH_Individual::operator=(BH_Individual *ind)
{
	if (N != ind->N && N != 0)
	{
		throw "BH个体原子数不对等";
		return;
	}

	for (int i=0; i < 3*N; i++)
		cood[i] = ind->cood[i];
	energy = ind->energy;
}

BasinHopping::BasinHopping(int N)
{
	cout<<"bh"<<endl;
	this->N = N;
	this->current.setN(N);
	this->best.setN(N);
	initDate();
}

BasinHopping::~BasinHopping(void)
{
}

void BasinHopping::initDate(){
	this->P_E = PE_Tool::FS_E;
	this->P_F = PE_Tool::FS_F;

	double r0= 2.75;
	for(int i=0; i<3 *N; i++)
		current.cood[i] = (RAND-0.5) * r0 * pow((double)N,1.0/3);

	BH_Individual temp(N);
	temp = current;
	current.energy = Tool::Local(temp.cood,N,P_E,P_F);
	best = current;
}

void BasinHopping::start(){
	successRate = 0;
	int success = 0;
	int all = 0;


}

void MC()
{

}