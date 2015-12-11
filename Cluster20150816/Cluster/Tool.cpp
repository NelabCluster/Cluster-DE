#include "Tool.h"
#include "PE_Tool.h"
#include <iostream>
#include <stdlib.h>
extern string EnergyName;
using namespace std;

Tool::Tool(void)
{
	
}

Tool::~Tool(void)
{
}

void Tool:: printCood(double *cood, int N, char *root)
{
	FILE *fp = fopen(root,"w");
	for(int i=0;i<N;i++)
	{
		fprintf(fp,"%lf\n",cood[i]);
	}
	fclose(fp);
}

void Tool:: readCood(double *cood, int N, char *root)
{
	FILE *fp = fopen(root,"r");
	for(int i=0;i<N;i++)
	{
		fscanf(fp,"%lf",cood+i);
	}
	fclose(fp);
}

bool Tool::readDiamond(double *cood, int N,int differ /* = 0 */)
{
	char path[100];
	char line[100];
	int atomN;
	double E;
	int afterDiffer = N - differ;
	sprintf(path,"%s_%d.txt",EnergyName.c_str(),afterDiffer);

	FILE *fp = fopen(path,"r");
	if (fp == NULL)
	{
		return false;
	}
	fscanf(fp,"%d",&atomN);
	fscanf(fp,"%s  %lf",line,&E);
	for (int i=0;i<afterDiffer;i++)
	{
		fscanf(fp,"%s %lf %lf %lf",line,cood+i,cood+N+i,cood+2*N+i);
	}
	fclose(fp);
	return true;
}

void Tool::printfDiamond(double *cood, int N, double E,char *path)
{
	FILE *fp;
	int i;
	double *x,*y,*z;
	x = cood; y = cood + N; z = cood + 2*N;
	fp = fopen(path,"w");
	fprintf(fp,"%d\n",N);
	fprintf(fp,"Fe,NCS  %lf\n",E);
	for(i=0;i<N;i++){
			fprintf(fp,"Fe\t%lf\t%lf\t%lf\n",x[i],y[i],z[i]);
	}
	fclose(fp);
}

void Tool::coodToDiamond(char *root, int N, char *path, PEnergy energy)
{
	double E;
	double *cood,*R;
	cood = (double *)malloc(3*N*sizeof(double));
	R = (double *)malloc(N*N*sizeof(double));
	readCood(cood,3*N,root);
	Distance(cood,R,N);
	E = energy(R,N);
	printfDiamond(cood,N,E,path);
	free(cood);
	free(R);
}

void Tool::Distance(double *cood, double *R, int N)
{
	int i,j;
	double *x,*y,*z;
	
	x = cood; y = cood + N; z = cood + 2 * N; 
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			*(R+i*N+j)=(x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);;
			*(R+i*N+j)= sqrt(*(R+i*N+j));
			*(R+j*N+i) = *(R+i*N+j);
		}
	}
}

double Tool::Local(double *cood, int N,PEnergy PE_E, PForce PE_F)
{
	int i;
	int s = 1;
	double alpha = 0.01;
	double *R,*tempR,*tempCood;
	double *tempX,*tempY,*tempZ;
	double E0,E1;
	double *FF;
	double *FX,*FY,*FZ,*x,*y,*z;
	double Fmax;
	double e = 1E-6;
	int smax = 3000;
	double alphaE = 1E-10;

	R = (double *)calloc(N*N,sizeof(double));
	tempR = (double *)calloc(N * N,sizeof(double));
	tempCood = (double *)calloc(3 * N,sizeof(double));
	FF = (double *)calloc(3 * N,sizeof(double));

	x = cood; y = cood + N; z = cood + 2 * N;
	FX = FF; FY = FF + N; FZ = FF + 2 * N; 
	tempX = tempCood; tempY = tempCood + N; tempZ = tempCood + 2*N;
//	Tool::printCood(cood,N);
//	Tool::readCood(cood,3* N,"s.txt");
	Distance(cood,R,N);
//	E0 = PE_Tool::FS_E(R,N);
//	Fmax = PE_Tool::FS_F(cood,R,FF,N);
//	E0 = PE_Tool::Johnson_E(R,N);
//	Fmax = PE_Tool::Johnson_F(cood,R,FF,N);
//	Tool::printCood(FF,N,"FF.txt");
	E0 = PE_E(R,N);
	Fmax = PE_F(cood,R,FF,N);
	while(1)
	{
		//cout<<"Fmax"<<Fmax<<endl;
		if(Fmax <= e || s > smax || alpha < e)
		{

			if(s>smax)
				cout<<"Fmax"<<Fmax<<"  "<<s<<"  "<<alpha<<"  "<<E0<<endl;
				
				//getchar();
			//	if(E0>-55)
			//	{
				//	Tool::printfDiamond(cood,N,E0,"badCluster.txt");
				//	printf("badCluster");
				//	getchar();
			//	}
			break;
		}
	/*	for(i = 0; i < 3 * N; i++){
			FF[i] = (FF[i]>0.5)?0.5:((FF[i]<-0.5)?-0.5:0);
		}*/
		for(i=0; i<N ;i++)
		{
			tempX[i] = x[i] + alpha * FX[i]/Fmax;
			tempY[i] = y[i] + alpha * FY[i]/Fmax;
			tempZ[i] = z[i] + alpha * FZ[i]/Fmax;
		//	tempX[i] = x[i] + alpha * FX[i];
		//	tempY[i] = y[i] + alpha * FY[i];
		//	tempZ[i] = z[i] + alpha * FZ[i];
		}

		Distance(tempCood,tempR,N);
	//	E1 = PE_Tool::FS_E(tempR,N);
	//	E1 = PE_Tool::Johnson_E(tempR,N);
		E1 = PE_E(tempR,N);
		if(E1 < E0)
		{
			alpha = alpha * 1.1;
			E0 = E1;
			for(i = 0 ;i < N;i ++)
			{
				x[i] = tempX[i];
				y[i] = tempY[i];
				z[i] = tempZ[i];
			}
			for(i=0;i<N*N;i++)
				R[i] = tempR[i];
		//	Fmax = PE_Tool::FS_F(cood,R,FF,N);
		//	Fmax = PE_Tool::Johnson_F(cood,R,FF,N);
			Fmax = PE_F(cood,R,FF,N);
		}
		else
			alpha = alpha * 0.6;
		s ++;
	}
	free(R);
	free(tempR);
	free(tempCood);
	free(FF);
	return E0;
}

double Tool::LocalWithR(double *R,int N,PEnergy PE_E, PForce PE_F)
{
	int i,j;
	int s = 1;
	double alpha = 0.01;
	double *tempR;
	double E0,E1;
	double *FF;
	double Fmax;
	double e = 1E-6;
	int smax = 3000;
	double alphaE = 1E-10;

	tempR = (double *)calloc(N * N,sizeof(double));
	FF = (double *)calloc(N * N,sizeof(double));

	E0 = PE_E(R,N);
	Fmax = PE_F(NULL,R,FF,N);
	while(1)
	{
		if(Fmax <= e || s > smax || alpha < e)
		{
			if(s>smax)
				cout<<"Fmax"<<Fmax<<"  "<<s<<"  "<<alpha<<"  "<<E0<<endl;
			break;
		}

		for(i=0; i<N-1 ;i++)
		{
			for(j=i+1;j<N;j++)
				*(tempR+i*N+j) = *(R+i*N+j) + alpha * (*(FF+i*N+j))/Fmax;
		}

		E1 = PE_E(tempR,N);
		if(E1 < E0)
		{
			alpha = alpha * 1.1;
			E0 = E1;

			for(i=0; i<N-1 ;i++)
			{
				for(j=i+1;j<N;j++)
					*(R+i*N+j) = *(tempR+i*N+j);
			}

			Fmax = PE_F(NULL,R,FF,N);
		}
		else
			alpha = alpha * 0.6;
		s ++;
	}

	free(tempR);
	free(FF);
	return E0;
}

int *Tool::RandPerm(int N, int K)
{
	int *allP;
	int *p;
	allP = (int *)malloc(N * sizeof(int));
	p = (int *)malloc(K * sizeof(int));
	for(int i = 0; i < N; i++)
		allP[i] = i;
	for(int i=0;i < K;i++)
	{
		int point = rand()%(N-i);
		int temp = allP[i];
		allP[i] = allP[point+i];
		allP[point+i] = allP[i];
		p[i] = allP[i];
	}
	free(allP);
	return p;
}

void Tool::SphereCutSplice(double *fr, double *mr, double *child1, double *child2, int N)
{
	int j,k,kexist,kmid,ncross;
	double x,y,z,mid,frcentre[3],mrcentre[3];
	double *fr3_x,*fr3_y,*fr3_z,*mr3_x,*mr3_y,*mr3_z;
	int *existId;
	double *frd,*mrd,*frdrank,*mrdrank;

//	cout<<1;

	existId = (int *)malloc(N * sizeof(int));
	frd = (double *)malloc(N * sizeof(double));
	mrd = (double *)malloc(N * sizeof(double));
	frdrank = (double *)malloc(N * sizeof(double));
	mrdrank = (double *)malloc(N * sizeof(double));

	for(j=0;j < N;j++){//father
		fr3_x = fr;
		fr3_y = fr + N;
		fr3_z = fr + 2 * N;
		mr3_x = mr;//mother
		mr3_y = mr + N;
		mr3_z = mr + 2 * N;
	}
	x=0;y=0;z=0;
	for(j=0;j < N;j++){
		x = x + fr3_x[j];
		y = y + fr3_y[j];
		z = z + fr3_z[j];	
	}
	frcentre[0] = x/N;
	frcentre[1] = y/N;
	frcentre[2] = z/N;
	x=0;y=0;z=0;
	for(j = 0;j < N;j++){
		x = x + mr3_x[j];
		y = y + mr3_y[j];
		z = z + mr3_z[j];	
	}
	mrcentre[0] = x / N;
	mrcentre[1] = y / N;
	mrcentre[2] = z / N;
	//centre
	for(j = 0; j < N;j++){
		fr3_x[j] = fr3_x[j] - frcentre[0];
		fr3_y[j] = fr3_y[j] - frcentre[1];
		fr3_z[j] = fr3_z[j] - frcentre[2];
		mr3_x[j] = mr3_x[j] - mrcentre[0];
		mr3_y[j] = mr3_y[j] - mrcentre[1];
		mr3_z[j] = mr3_z[j] - mrcentre[2];
	}
	//distance
	for(j = 0;j < N;j++){
		frd[j] = sqrt(fr3_x[j] * fr3_x[j] + fr3_y[j] * fr3_y[j] + fr3_z[j] * fr3_z[j]);
		mrd[j] = sqrt(mr3_x[j] * mr3_x[j] + mr3_y[j] * mr3_y[j] + mr3_z[j] * mr3_z[j]);
		frdrank[j]=frd[j];
		mrdrank[j]=mrd[j];
	}
	//let not two ones are exactly the same
	for( k = 0;k < N-1;k++){
		for(j=k+1;j<N;j++){
			if(frd[k]==frd[j]){
				frd[k]=frd[k]+(2 + RAND)*(1e-16);
				frdrank[k]=frd[k];
			}
			if(mrd[k]==mrd[j]){
				mrd[k]=mrd[k]+(2 + RAND)*(1e-16);
				mrdrank[k]=mrd[k];
			}
		}
	}
	//rank
	for(k=N;k>1;k--){
		for(j=1;j<k;j++){
			if(frdrank[j-1]>frdrank[j]){
				mid=frdrank[j];
				frdrank[j]=frdrank[j-1];
				frdrank[j-1]=mid;
			}
		}
		for(j=1;j<k;j++){
			if(mrdrank[j-1]>mrdrank[j]){
				mid=mrdrank[j];
				mrdrank[j]=mrdrank[j-1];
				mrdrank[j-1]=mid;
			}
		}
	}
	//crossover by sphere
	//check if there exists such a sphere
	kexist=0;
	for(k=1;k<N-1;k++){
		if(!((frdrank[k-1]>mrdrank[k]) || (frdrank[k]<mrdrank[k-1]))){
			existId[kexist]=k;
			kexist=kexist+1;
		}
	}

	//cout<<"k="<<kexist<<";";
	if(kexist!=0){
		//random number of crossover atoms
		
		kmid=(int)ceil(kexist* RAND);
		if(kmid==0)
			kmid = 1;
	//	cout<<"kmid="<<kmid<<";";
		ncross=existId[kmid-1];
		k=0;
		for(j=0;j<N;j++){
			if(frd[j]<frdrank[ncross]-1e-16){
				*(child1 + k) = fr3_x[j];
				*(child1 + N + k) = fr3_y[j];
				*(child1 + 2 * N + k) = fr3_z[j];
				k=k+1;
			}
			if(mrd[j]>mrdrank[ncross-1]+1e-16){
				*(child1 + k) = mr3_x[j];
				*(child1 + N + k) = mr3_y[j];
				*(child1 + 2 * N + k) = mr3_z[j];
				k=k+1;
			}
		}
		//cout<<"child1="<<k<<";";
		k=0;
		for(j=0;j<N;j++){
			if(frd[j]>frdrank[ncross-1]+1e-16){
				*(child2 + k) = fr3_x[j];
				*(child2 + N + k) = fr3_y[j];
				*(child2 + 2 * N + k) = fr3_z[j];
				k=k+1;
			}
			if(mrd[j]<mrdrank[ncross]-1e-16){
				*(child2 + k) = mr3_x[j];
				*(child2 + N + k) = mr3_y[j];
				*(child2 + 2 * N + k) = mr3_z[j];
				k=k+1;
			}
		}
		//cout<<"child2="<<k<<";";
	}
	else{
		for(j=0;j<3 * N;j++){
			*(child1 + j) = fr[j];
			*(child2 + j) = mr[j];
		}
	}
	//cout<<2;
	free(existId);
	free(frd);
	free(mrd);
	free(frdrank);
	free(mrdrank);
	//cout<<endl;
}

void Tool::testR(int N){
	double *cood = new double[3*N];
	double *R = new double[N*N];
	double E0,E1;

	//随机产生一个结构
	for (int i=0; i < 3*N; i++)
	{
		cood[i] = (RAND-0.5) * 2.75 * pow((double)N,1.0/3);
	}
	Distance(cood,R,N);
	E0 = PE_Tool::FS_E(R,N);
	//对R进行最速下降法
	E1 = LocalWithR(R,N,PE_Tool::FS_E,PE_Tool::FS_Fr);
	printf("原子总数为%d;\n",N);
	printf("初始能量:%lf\n",E0);
	printf("优化r后的能量:%lf\n",E1);

	//得到R下，三角形个数和无法围成合法三角形的个数
	testThree(R,N);
}

bool Tool::testThree(double *R,int N){
	int numOfAll = 0;
	int numOfBad = 0;
	bool rational = true;
	printf("无法构成三角形的节点和三条边的长度如下：\n");
	for (int i = 0; i < N-2; i++)
	{
		for(int j=i+1; j < N-1; j++){
			for (int m = j+1; m < N; m++)
			{
				numOfAll ++;
				double a = *(R + i * N + j);
				double b = *(R + i* N + m);
				double c = *(R + j * N + m);
				if (a + b<c || a + c < b || b + c <a)
				{
					rational = false;
					numOfBad ++;
					printf("%d,%d,%d\t%lf %lf %lf\n",i,j,m,a,b,c);
				} 
			}
		}
	}
	printf("一共%d三角形,其中坏的有%d个\n",numOfAll,numOfBad);
	return rational;
}



