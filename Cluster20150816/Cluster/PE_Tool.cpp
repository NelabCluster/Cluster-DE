#include "PE_Tool.h"


PE_Tool::PE_Tool(void)
{
}

PE_Tool::~PE_Tool(void)
{
}




double PE_Tool::FS_E(double *R,int N)
{
	int i,j;
	double d,A,beta,c,c0,c1,c2;
	double r;
	double tempV,tempP;
	double *VEN,*PEN;
	double E = 0;
	double minR = 1.6,tempPmin;

	//Fe
	d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;

	VEN = (double *)calloc(N,sizeof(double));
	PEN = (double *)calloc(N,sizeof(double));
	
//	tempVmin = (minR-c)*(minR-c)*(c0+c1*minR+c2*minR*minR);
	tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
	for(i = 0; i < N-1; i++)
		for(j = i+1; j < N; j++)
		{
			r = *(R + i*N + j);
			tempV = (r>=c)?0 : (r-c)*(r-c)*(c0+c1*r+c2*r*r);
			tempP = (r>=d)?0 : ((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);

			VEN[i] += tempV/2;
			VEN[j] += tempV/2;
			PEN[i] += tempP;
			PEN[j] += tempP;
		}

	for(i=0 ;i < N; i++){
		PEN[i] = sqrt(PEN[i]);
		E += VEN[i] - A * PEN[i];
	//	cout<<VEN[i]<<" "<<A * PEN[i]<<" "<<VEN[i]-A*PEN[i]<<endl;
	}
//	Tool::printCood(PEN,N,"PEN.txt");
	free(VEN);
	free(PEN);
	return E;
}

double PE_Tool::FS_F(double *cood,double *R,double *FF, int N)
{
	int i,j;	
	double d,A,beta,c,c0,c1,c2;
	double r;
	double tempdV,tempP,tempdP;
	double *PEN;
	double maxF = 0,FK;
	double minR = 1.6,tempPmin;
	double *x,*y,*z,*FX,*FY,*FZ;

	PEN = (double *)calloc(N,sizeof(double));
	
	//Fe
	d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;
	x = cood; y = cood + N; z = cood + 2 * N; FX = FF; FY = FF + N; FZ = FF + 2 * N;

	for(i = 0; i < 3 * N; i++)
		FF[i] = 0;

	tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;
	//tempdPmin = 2 * (minR -d) + 3 * beta * (minR - d) * (minR - d) / d;
	//tempdVmin = (2*(minR - c)*(c0 + c1 * minR + c2 * minR * minR) + (c1 + 2 * c2 * minR)* (minR - c) * (minR - c)) * 2;
	for(i=0;i<N-1; i++)
		for(j=i+1;j<N;j++)
		{
			r = *(R + i*N + j);
			tempP = (r>=d)?0:((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);

			PEN[i] += tempP;
			PEN[j] += tempP;
		}

	for(i = 0; i < N ;i ++)
		PEN[i] = (PEN[i] == 0)?0:1 / sqrt(PEN[i]) / 2;

	for(i = 0; i < N-1;i ++)
		for(j=i+1;j<N;j++)
		{
			r = *(R + i*N + j);
			
		//	tempdV = (r>=c)?0: (r<=minR)?tempdVmin : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
		//	tempdP = (r>=d)?0 : (r<=minR)?tempdPmin : 2 * (r -d) + 3 * beta * (r - d) * (r - d) / d;
			tempdV = (r>=c)?0 : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
			tempdP = (r>=d)? 0 : (2 * (r -d) + 3 * beta * (r - d) * (r - d) / d);

			tempdP = (PEN[i] + PEN[j]) * tempdP;
			
			FK = -tempdV / 2 + A * tempdP;
			
			FX[i]=FX[i]+FK*(x[i]-x[j])/r;
		//	cout<<j<<" "<<FX[i]<<" "<<FK<<" "<<r<<endl;
			FX[j]=FX[j]-FK*(x[i]-x[j])/r;
			FY[i]=FY[i]+FK*(y[i]-y[j])/r;
			FY[j]=FY[j]-FK*(y[i]-y[j])/r;
			FZ[i]=FZ[i]+FK*(z[i]-z[j])/r;
			FZ[j]=FZ[j]-FK*(z[i]-z[j])/r;
		}
	for(i = 0;i < 3 * N; i++)
		maxF = (FF[i]>maxF)?FF[i]:maxF;

	free(PEN);
	return maxF;
}

double PE_Tool::FS_Fr(double *cood,double *R,double *FF, int N)
{
	int i,j;	
	double d,A,beta,c,c0,c1,c2;
	double r;
	double tempdV,tempP,tempdP;
	double *PEN;
	double maxF = 0,FK;
	double minR = 1.6,tempPmin;

	PEN = (double *)calloc(N,sizeof(double));

	//Fe
	d = 3.569745;A = 1.828905;beta = 1.8;c = 3.40;c0 = 1.2371147;c1 = -0.3592185; c2 = -0.0385607;

	for(i = 0; i < N * N; i++)
		FF[i] = 0;

	tempPmin = (minR-d)*(minR-d)+beta*(minR-d)*(minR-d)*(minR-d)/d;

	for(i=0;i<N-1; i++)
		for(j=i+1;j<N;j++)
		{
			r = *(R + i*N + j);
			tempP = (r>=d)?0:((r<=minR)?tempPmin:(r-d)*(r-d)+beta*(r-d)*(r-d)*(r-d)/d);

			PEN[i] += tempP;
			PEN[j] += tempP;
		}

	for(i = 0; i < N ;i ++)
		PEN[i] = (PEN[i] == 0)?0:1 / sqrt(PEN[i]) / 2;

	for(i = 0; i < N-1;i ++)
		for(j=i+1;j<N;j++)
		{
			r = *(R + i*N + j);

			tempdV = (r>=c)?0 : (2*(r - c)*(c0 + c1 * r + c2 * r * r) + (c1 + 2 * c2 * r)* (r - c) * (r - c)) * 2;
			tempdP = (r>=d)? 0 : (2 * (r -d) + 3 * beta * (r - d) * (r - d) / d);

			tempdP = (PEN[i] + PEN[j]) * tempdP;

			FK = -tempdV / 2 + A * tempdP;
			*(FF + i*N + j) = FK;
			maxF = (abs(FK)>maxF)?abs(FK):maxF;
		}

		free(PEN);
		return maxF;
}

double PE_Tool::Johnson_E(double *R, int N)
{
	double rc = 4.5;

	double re = 2.481987;
	double fe = 1.885957;
	double roue = 20.041463;
	double alpha = 9.818270;
	double beta = 5.236411;
	double A = 0.392811; 
	double B = 0.646243;
	double kappa = 0.170306;
	double lambda = 0.340613;
	double Fn0 = -2.534992;
	double Fn1 = -0.059605;
	double Fn2 = 0.193065;
	double Fn3 = -2.282322;
	double F0 = -2.54;
	double F1 = 0;
	double F2 = 0.200269;
	double F3 = -0.148770;
	double eta = 0.391750;
	double Fe = -2.539945;

	double roun = 0.85 * roue;
	double rouo = 1.15 * roue;

	int i,j;
	double E,tempRou;
	double tempR,Phi1,Phi2;
	double *Phi,*rou;

	Phi = (double *)calloc(N,sizeof(double));
	rou = (double *)calloc(N,sizeof(double));

	for(i = 0; i < N-1; i++)
	{
		for(j = i+1; j < N; j++)
		{
			tempR = *(R + i*N + j);
			if(tempR > rc)
			{
				continue;
			} else {
			tempR= tempR / re;
			Phi1 = exp(-alpha * (tempR - 1)) / (1 + pow(tempR - kappa,20));
			Phi2 = exp(-beta * (tempR - 1)) / (1 + pow(tempR - lambda,20));
			Phi[i] += (A * Phi1 - B * Phi2) / 2;
			Phi[j] += (A * Phi1 - B * Phi2) / 2;
			rou[i] += fe * Phi2;
			rou[j] += fe * Phi2;
			}
		}
	}
	
	for(i = 0; i < N; i++)
	{
		tempRou = rou[i];
		if(tempRou < roun)
		{
			tempRou = tempRou / roun - 1;
			tempRou = Fn0 + Fn1 * tempRou + Fn2 * tempRou * tempRou + Fn3 * tempRou * tempRou * tempRou;	
		} else if(tempRou < rouo)
		{
			tempRou = tempRou / roue - 1;
			tempRou = F0 + F1 * tempRou + F2 * tempRou * tempRou + F3 * tempRou * tempRou * tempRou;
		} else
		{
			tempRou = pow(tempRou/roue,eta);
			tempRou = Fe * (1-log(tempRou)) * tempRou; 
		}
		rou[i] = tempRou;
	}

	E = 0;
	for(i = 0; i < N; i ++)
	{
		E += Phi[i] + rou[i];
	}

	free(Phi);
	free(rou);
	return E;
}

double PE_Tool::Johnson_F(double *cood, double *R, double *FF, int N)
{
	double rc = 4.5;

	double re = 2.481987;
	double fe = 1.885957;
	double roue = 20.041463;
	double alpha = 9.818270;
	double beta = 5.236411;
	double A = 0.392811; 
	double B = 0.646243;
	double kappa = 0.170306;
	double lambda = 0.340613;
	double Fn0 = -2.534992;
	double Fn1 = -0.059605;
	double Fn2 = 0.193065;
	double Fn3 = -2.282322;
	double F0 = -2.54;
	double F1 = 0;
	double F2 = 0.200269;
	double F3 = -0.148770;
	double eta = 0.391750;
	double Fe = -2.539945;

	double roun = 0.85 * roue;
	double rouo = 1.15 * roue;

	int i,j;
	double r,tempR,Phi2,dPhi1,dPhi2,dPhi,drou,tempRou,FK,maxF;
	double *FX,*FY,*FZ,*x,*y,*z;
	double *rou;
	
	rou = (double*)calloc(N,sizeof(double));
	x = cood; y = cood + N; z = cood + 2 * N; FX = FF; FY = FF + N; FZ = FF + 2 * N;
	for(i = 0; i < 3 * N; i++)
		FF[i] = 0;
	
	for(i = 0; i < N-1; i++)
	{
		for(j = i+1; j < N; j++)
		{
			tempR = *(R + i*N + j);
			tempR = tempR / re;
			Phi2 = exp(-beta * (tempR - 1)) / (1 + pow(tempR - lambda,20));
			
			rou[i] += fe * Phi2;
			rou[j] += fe * Phi2;
		}	
	}

	for(i = 0; i < N; i++)
	{
		tempRou = rou[i];
		if(tempRou < roun)
		{
			tempRou = tempRou / roun - 1;
			tempRou = (Fn1 + 2*Fn2 * tempRou + 3*Fn3 * tempRou * tempRou) / roun;
		} else if(tempRou < rouo)
		{
			tempRou = tempRou / roue - 1;
			tempRou = (F1 + 2*F2 * tempRou + 3*F2 * tempRou * tempRou) / roue;
		} else
		{
			tempRou = tempRou / roue;
			tempRou = -( Fe * eta * log(pow(tempRou,eta)) * pow(tempRou,(eta-1)) ) / roue; 
		}
		rou[i] = tempRou;
	}

	for( i = 0; i < N-1; i++)
	{
		for( j = i+1; j < N; j++)
		{
			r = *(R + i*N +j);
			tempR = r / re;
			dPhi2 = -(exp(-beta * (tempR-1)) * ( (beta + beta * pow((tempR-lambda),20) + 20 * pow((tempR-lambda),19)) / (1 + pow((tempR-lambda),40) + 2*pow((tempR-lambda),20)) ) ) / re;
			if(r<1.6){
				dPhi = -100 + 36.1698 * r;
			}
			else{
			dPhi1 = -(exp(-alpha * (tempR-1)) * ( (alpha + alpha * pow((tempR-kappa),20) + 20 * pow((tempR-kappa),19)) / (1 + pow((tempR-kappa),40) + 2*pow((tempR-kappa),20)) )) / re;
			dPhi = (r > rc)?0:A * dPhi1 - B * dPhi2;
			}
			drou = fe * dPhi2;
			drou = (r>rc)?0:(rou[i] + rou[j]) * drou;
			FK = -dPhi - drou;
			
			FX[i]=FX[i]+FK*(x[i]-x[j])/r;
			FX[j]=FX[j]-FK*(x[i]-x[j])/r;
			FY[i]=FY[i]+FK*(y[i]-y[j])/r;
			FY[j]=FY[j]-FK*(y[i]-y[j])/r;
			FZ[i]=FZ[i]+FK*(z[i]-z[j])/r;
			FZ[j]=FZ[j]-FK*(z[i]-z[j])/r;
		}
	}

	maxF = 0;
	for(i = 0;i < 3 * N; i++)
	{
	/*	if(FF[i]>10)
			FF[i] = 10;
		if(FF[i]<-10)
			FF[i] = -10;*/
		maxF = (FF[i]>maxF)?FF[i]:maxF;
	}
	
	free(rou);
	return maxF;
}

PEnergy PE_Tool::GetPEnergy(PotentialEnergyType type)
{
	switch(type)
	{
	case PotentialEnergyType_FS:
		return FS_E;
		break;
	case PotentialEnergyType_Johnson:
		return Johnson_E;
		break;
	}
	return NULL;
}

PForce PE_Tool::GetPForce(PotentialEnergyType type)
{
	switch(type)
	{
	case PotentialEnergyType_FS:
		return FS_F;
		break;
	case PotentialEnergyType_Johnson:
		return Johnson_F;
		break;
	}
	return NULL;
}