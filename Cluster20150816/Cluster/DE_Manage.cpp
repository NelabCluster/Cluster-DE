#include "DE_Manage.h"
#include "Common.h"
 #include <direct.h>
#include "Tool.h"
extern string EnergyName;

DE_Manage::DE_Manage(void)
{
	this->popNum = 0;
	this->pop = NULL;
	this->numList = NULL;
}

DE_Manage::~DE_Manage(void)
{
	cout<<"DE_Manage销毁"<<endl;
}

void DE_Manage::AddPop(DE *p)
{
	int i;
	DE *temp;

	if(this->popNum == 0) 
		this->pop = p;
	else
	{
		temp = this->pop;
		for( i = 1; i < this->popNum; i++)
			temp = temp->next;
		temp->next = p;
	}
	p->next = this->pop;
	this->popNum ++;
}

void DE_Manage::removeAllPop()
{
	int i;
	DE *p,*q;
	p = q = this->pop;
	for(i = 0; i<this->popNum; i++)
	{
		q = p->next;
		delete p;
		p = q;
	}
	this->pop = NULL;
	this->popNum = 0;
}

double DE_Manage::MultiPop_Master(int N)
{
	int i,j,S;
	DE *p;
	char root[20];
	double bestEnergy;
	FILE *energyFile,*fitnessRatioFile;

	p = new DE(DERAND1,N);
	this->AddPop(p);
	p = new DE(DEBEST1,N);
	this->AddPop(p);
	p = new DE(DERANDTOBEST1,N);
	this->AddPop(p);
	
	N = p->N;
	
	sprintf(root,"%s_%d.txt",EnergyName.c_str(),N);
	S = 10;
	best = DE_Individual::Individual(p->N);
	
	atomNum = N;
	readBestResult();

	char energyFileName[50];
	sprintf(energyFileName,"%s_%d",EnergyName.c_str(),N);
	mkdir(energyFileName);
	sprintf(energyFileName,"%s_%d\\%d.txt",EnergyName.c_str(),N,current);
	energyFile = fopen(energyFileName,"w+");

	char fitnessRatioFileName[50];
	sprintf(fitnessRatioFileName,"%s_%d\\%d_diversity.txt",EnergyName.c_str(),N,current);
	fitnessRatioFile = fopen(fitnessRatioFileName,"w+");

	for(i = 0; i < iterations; i++)
	{
		p = this->pop;				//进化一代
		for(j = 0; j < this->popNum; j++)	
		{
			p->OneGeneration();
			p = p->next;
		}

		best->energy = 0;			//得到每个种群最优个体
		p = this->pop;
		for(j = 0; j < this->popNum; j++)	
		{
			if(p->bestPop.energy < best->energy)
				best->Copy(&p->bestPop,p->N);
			p = p->next;
		}
		cout<<i<<" "<<best->energy;

		p = this->pop;
		for(j = 0; j < this->popNum; j++)	
		{
			cout<<" "<<p->bestPop.energy;
			p = p->next;
		}
		cout<<endl;
		
		this->printfEnergy(energyFile, i);
		this->printfFitnessRatio(fitnessRatioFile, i);


		if(i % S != 0)
			continue;

		cout<<"种群交换最优个体"<<endl;
		p = this->pop;				//将最优个体发送给每个种群
		for(j = 0; j < this->popNum; j++)	
		{
			p->RecIndividual(best);
			p = p->next;
		}

		if (best->energy < bestResult)
		{
			Tool::printfDiamond(best->cood,N,best->energy,root);
			bestResult = best->energy;
		}
		
		bestEnergy = best->energy;
	}

	fclose(energyFile);
	fclose(fitnessRatioFile);
	this->removeAllPop();
	free(best->cood);
	free(best);
	
	return bestEnergy;
}

void DE_Manage::inputPara(){
	int i,num;
	char c;
	char numStr[100];

	cout<<"请输入原子个数(如13,14,15-18):"<<endl;
	do 
	{
		gets(numStr);
		i = 0;
		c = numStr[i];
		while (c!=0)
		{
			if(!(c==','||c=='-'||(c>='0'&&c<='9')))
				break;
			c = numStr[++i];
		}
		if (c == 0)
			break;
		cout<<"存在非法字符,请重新输入:"<<endl;
	} while (1);

	i = 0;
	c = numStr[i];
	num = 0;
	while (c!='\0')
	{
		if (c >= '0' && c <= '9')
			num = num * 10 + (c - '0');

		if (c==','){
			AddToTail(&numList,num);
			num = 0;
		}

		if (c=='-')
		{
			int startNum = num;
			int endNum = 0;
			while (c!=',')
			{
				c = numStr[++i];
				if (c >= '0' && c <= '9')
					endNum = endNum * 10 + (c - '0');
			}
			for (int j = startNum; j <= endNum; j++)
				AddToTail(&numList,j);
			num = 0;
		}
		c = numStr[++i];
	}
	cout<<"请输入最大迭代次数:"<<endl;
	cin>>iterations;
	cout<<"请输入重复次数："<<endl;
	cin>>repeat;
}

void DE_Manage::start(){
	if (numList==NULL)
	{
		cout<<"没有设置原子数序列"<<endl;
		return;
	}
	ListNode *p;
	p = numList;
	while (p!=NULL)
	{	
		current = 0;
		for (current=0;current<repeat;current++)
		{
			MultiPop_Master(p->m_nValue);
		}
		
		p = p->m_pNext;
	}
}

void DE_Manage::Repeat(){
	int N,success =0;
	double E,objEnergy,successRate = 0;

	cout<<"请输入原子数："<<endl;
	cin>>N;
	cout<<"请输入最大迭代次数："<<endl;
	cin>>iterations;
	cout<<"请输入目标能量："<<endl;
	cin>>objEnergy;
	cout<<"请输入重复次数："<<endl;
	cin>>repeat;

	current = 0;
	for (current = 0; current < repeat; current++ )
	{
		E = MultiPop_Master(N);
		if (fabs(E - objEnergy) < 0.01)
		{
			success ++;
		}
	}
	successRate = success*100.0 / repeat;
	cout<<"一共进行"<<repeat<<"次，其中成功"<<success<<"次，成功率为"<<successRate<<"%"<<endl;
}

void DE_Manage::readBestResult(){
	char path[100];
	char line[100];
	int atomN;
	double E;
	sprintf(path,"%s_%d.txt",EnergyName.c_str(),atomNum);

	FILE *fp = fopen(path,"r");
	if (fp != NULL)
	{
		fscanf(fp,"%d",&atomN);
		fscanf(fp,"%s  %lf",line,&E);
		fclose(fp);
		bestResult = E;
	} 
	else
	{
		bestResult = 0;
	}
}

void DE_Manage::printfEnergy(FILE *fp, int index)
{
	if (fp == NULL) return;

	fprintf(fp,"%d\t%f",index,best->energy);
	DE *p = this->pop;
	for(int i = 0; i < this->popNum; i ++)	
	{
		fprintf(fp,"\t%f",p->bestPop.energy);
		p = p->next;
	}
	fprintf(fp,"\n");
}

void DE_Manage::printfFitnessRatio(FILE *fp, int index)
{
	if (fp == NULL) return;

	fprintf(fp,"%d",index);
	DE *p = this->pop;
	for(int i = 0; i < this->popNum; i ++)	
	{
		fprintf(fp,"\t%f",p->getFintessRatio());
		p = p->next;
	}
	fprintf(fp,"\n");
}