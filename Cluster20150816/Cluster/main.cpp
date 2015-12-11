#include "Common.h"
#include "DE_Manage.h"
#include "BasinHopping.h"
#include "NELbfgs.h"
#include "PE_Tool.h"
#include "Tool.h"
#include "LocalTool.h"


string EnergyName = "FS";

void main()
{
	srand((unsigned)time(NULL));
	//FS  or Johnson
	EnergyName = "FS";
	DE_Manage *manage = new DE_Manage();
//	manage->Repeat();
	manage->inputPara();
	manage->start();
	system("pause");

}

void twoLocalMethod()
{

	double *cood1,*cood2;
	int N = 90;
	double r0= 2.75;
	cood1 = new double[3 * N];
	cood2 = new double[3 * N];
	for (int i = 0;i < 3 * N;i ++) {
		cood1[i] = (RAND-0.5) * r0 * pow((double)N,1.0/3);
		cood2[i] = cood1[i];
	}

	double E1 = LocalTool::shareDefault()->localWithCoodAndPE(cood1,N,PE_Tool::FS_E,PE_Tool::FS_F);
	double E2 = Tool::Local(cood2,N,PE_Tool::FS_E,PE_Tool::FS_F);

	cout<<"E1:"<<E1<<endl;
	cout<<"E2:"<<E2<<endl;

	delete(cood1);
	delete(cood2);
	getchar();
}