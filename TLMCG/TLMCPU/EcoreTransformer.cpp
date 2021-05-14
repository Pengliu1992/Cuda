#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define SigmaE  0
#define Mu0 1.256637e-6
#define pi 3.1415927
#define Delt (0.05e-6)
#define MaxStep 2001
#define Np 132
#define Ns 39
#define CoilArea 0.5
#define TLM_Tolerance 1e-6
#define Reluct0 (500)
#define UseMemoryAdm 0
#define MaxNodeNumber 1200
#define MaxElemNumber 2400
#define VRMS 220000
#define NR1 (1.3332*0*Np/40)
#define NR2 (0.6666*0*Ns/20)
#define Freq 300000

#define Level 5
#define NoC 4//(level-1)
#define NoM 8//2*(level-1)

//define GET_TIME
#ifndef _TIMER_H_
#define _TIMER_H_
typedef struct
{
	double Y;
	double Y1;
	double Y2;
}interp_t;
#include <time.h>

#define GET_TIME(now){ \
	double t; \
	t=clock(); \
	now =t/CLOCKS_PER_SEC ;\
}
#endif

/*  Solve Diffusion Equation  */
void SM(double dt, double *Vc, double *Vzu2, double *Vau2, double Iu, double Id, double Uref, double *Vp, double *Vp_up, double *Vp_dn, double *vi, double *car);
int saveoutBH(double **B,double **H, int id, int count);
void saveoutL(double* L1, double *L2, double* M12,double* M21,int N);
void saveoutIVIV(double* I1, double *V1, double* I2,double* V2,int N);
void saveoutParaInt(int* Para,int N);
void saveoutPara(double* Para,int N);
int saveoutAxb(double **K,double **L, double **U,int N);
void LoadMeshInfo(int* NumNodes,int* NumElem,double** X,double** Y,int** VertI, int** VertJ, int** VertK);
void LU_De(double A[MaxNodeNumber][MaxNodeNumber], int n);
void LU_Solve(double** L, double ** U,double* B, double** X,int n);
void Gauss_Solve(double** A, double **B,double** X, int n);
void Load_BH(void);
void saveoutV(double V[],int N);
void LoadVe(double Ve[],int N);
void FEM_Data_Prepare();
void Apply_Memory_Adm_Matr();
void FEM_TLM(double Ip1,double Is1,double Ip2,double Is2,double Ip3,double Is3,double *AA);
void Inductance_Calculation(double Ip, double Is, double* Lw1,double*Lw2,double *Lm1,double *Lm2);
int Newton_Raphson_Element(double *X,double *Y,double *Z,double x0, double y0, double z0, double J0,double J1,double J2,double Ke01,double Ke12,double Ke20,int Id,double Ve[],double*RecB,double*RecH);
double Det(double x0,double x1,double x2,double y0,double y1,double y2,double z0,double z1,double z2);
int saveoutput(double *X,double *Y,double** S, int N, int Count);
double* Cubic_Spline(double y[], double x[]);
interp_t Interpo(double x);
//void LuSparse(double *B, double **X,int *ColL,int *RowL,double *ValL,int *ColU,int *RowU,double *ValU,double *Ud);
void LU_SparseSolve(double Ks_temp[MaxNodeNumber],double A[MaxNodeNumber],int N);
//Global parameters for FEM
const int CubicN=15;
double *CubicX1;
double CubicX[15]={1,1.21,1.44,1.69,1.96,2.25,2.56,2.89,3.24,3.61,4,4.41,4.84,5.29,5.76}; 
double CubicY[15]={663.15,970.45,1.4210e3,1.8947e3,2.7441e3,3.6172e3,4.9736e3,7.2343e3,1.1368e4,1.6931e4,3.0607e4,5.3051e4,7.9577e4,1.1368e5,1.3263e5};
double *ValL,*ValU,*Ud;
int *ColL,*RowL,*ColU,*RowU;
int MaxLength,NL,NU;
double G01[MaxElemNumber],G12[MaxElemNumber],G20[MaxElemNumber];

int NumNodes, NumElem;
int NCoil1,NCoil2,NCoil3,NCoil4,NCoil5,NCoil6,NCoil7,NCoil8,NCoil9,NCoil10,NCoil11,NCoil12;
int Coil1Id[80],Coil2Id[80],Coil3Id[80],Coil4Id[80],Coil5Id[80],Coil6Id[80],Coil7Id[80],Coil8Id[80],Coil9Id[80],Coil10Id[80],Coil11Id[80],Coil12Id[80];
int ClassEle[MaxElemNumber];
double Area[MaxElemNumber],a1[MaxElemNumber],b1[MaxElemNumber],c1[MaxElemNumber],a2[MaxElemNumber],b2[MaxElemNumber],c2[MaxElemNumber],a3[MaxElemNumber],b3[MaxElemNumber],c3[MaxElemNumber];
double V12_CIn[MaxElemNumber],V23_CIn[MaxElemNumber],V31_CIn[MaxElemNumber],V10_In[MaxElemNumber],V20_In[MaxElemNumber],V30_In[MaxElemNumber];
int *VertI,*VertJ,*VertK;
int BoundaryNodeId[MaxNodeNumber];
double BoundaryNodeValue[MaxNodeNumber];
double YCij[MaxElemNumber],YCi0[MaxElemNumber], Ve[MaxElemNumber],*sigma;
double RecordB,RecordH;
double *X,*Y;
double A0[MaxNodeNumber], A0p[MaxNodeNumber], A0s[MaxNodeNumber];
//Global parameters for circuit


double R1=5,R2=100000,L1=0,L2=0;


double Lw1,Lw2,Lm1,Lm2,dIp=1,dIs=1;	
double Lw10, Lw20, Lm10, Lm20, Ip20=0, Is20=0;
double Il1[MaxStep], Il2[MaxStep],Vpr[MaxStep],Ipr[MaxStep],Vse[MaxStep],Ise[MaxStep],Ikm[MaxStep],Vkm[MaxStep],SelfL1[MaxStep],SelfL2[MaxStep],MutualL12[MaxStep],MutualL21[MaxStep];
double *Vs1,*Vs2,*Vs3;
double Ihisp=0,Ihis1=0,Ihis2=0;
double T_start,T_end;
int Time_Step=0;
double Kb0, k01, k02, Kb1, k11, k12, Kb2, k21, k22, k00, k10, k20,ka1,kb1,kc1,ka2,kb2,kc2;
double Ip1, Ip2, Ip3, Is1, Is2, Is3;
double ParaC1 = 100e-12, ParaC2 = 100e-12, ParaC = 800e-12;
double Iterative_error = 1;
int main(){
	GET_TIME(T_start);
	FEM_Data_Prepare();
	
	Ip1 = 0; Ip2 = 0; Ip3 = 0;
	Is1 = 0; Is2 = 0; Is3 = 0;
	
	
	Inductance_Calculation(Ip2, Is2, &Lw1, &Lw2, &Lm1, &Lm2);
	//saveoutPara(A0,NumNodes);
	

	for (Time_Step = 1; Time_Step < MaxStep; Time_Step++)
	{
		if (Time_Step > 600)
			R2 = 100;

		if (Time_Step > 1200)
			R2 = 0.01;

		for (int ItCount = 0; ItCount < 1; ItCount++)
		{
			Inductance_Calculation(Ip2, Is2, &Lw1, &Lw2, &Lm1, &Lm2);
			//Solved simultaneously
			ka1 = 2 * Lw1 / Delt + R1 + NR1;
			kb1 = 2 * Lm2 / Delt;
			kc1 = 2 * Lw1 / Delt*Ihis1 + 2 * Lm2 / Delt*Ihis2 + Vs1[Time_Step];
			ka2 = 2 * Lm1 / Delt;                
			kb2 = R2 + NR2 + 2 * Lw2 / Delt;
			kc2 = 2 * Lw2 / Delt*Ihis2 + 2 * Lm1 / Delt*Ihis1;


			Ipr[Time_Step] = (kc2*kb1 - kc1*kb2) / (ka2*kb1 - ka1*kb2);
			Ise[Time_Step] = (kc2*ka1 - kc1*ka2) / (kb2*ka1 - kb1*ka2);

			Vpr[Time_Step] = Vs1[Time_Step] - Ipr[Time_Step] * R1;
			Vse[Time_Step] = -Ise[Time_Step] * R2;

			Ip2 = Ipr[Time_Step];
			Is2 = Ise[Time_Step];

		}


		Ihis1 = Ipr[Time_Step] * 2 - Ihis1;
		Ihis2 = Ise[Time_Step] * 2 - Ihis2;

	/*	if (Time_Step % 10 == 9)
		{

			Apply_Memory_Adm_Matr();

			printf("Memory Admt Matrix Applied on %d Time Step\n", Time_Step);
		}*/
	
	}
	GET_TIME(T_end);
	
	printf("\nElapsed time for FE-TLM Simulation is %.8f\n\n", (T_end - T_start));
	//saveoutV(Ve, 2020);

	saveoutIVIV(Ipr, Vpr, Ise, Vse, MaxStep);
	
	

	/*saveoutput(X,Y,Aout,NumNodes,Time_Step);
	
	saveoutL(SelfL1,SelfL2,MutualL12,MutualL21,MaxStep);*/

	return 0;
		
}
void Apply_Memory_Adm_Matr()
{
	int i,j;
	double x1, x2, x3, y1, y2, y3;
	static double Vr[MaxElemNumber];
	static double Kc[MaxNodeNumber][MaxNodeNumber];

	for (i = 0; i<NumNodes; i++)
		for (j = 0; j<NumNodes; j++)
		{
			Kc[i][j] = 0;
		}

	for (i = 0; i<NumElem; i++)
	{
		sigma[i] = 0;
		Vr[i] = 1 / Mu0;
		if (ClassEle[i] == 3)
		{
			Vr[i] = Ve[i];//material initialized
		}

	}
	//Generate Memory admittance matrix

	//Get Stiffness Matrix and Admittance matrix

	for (i = 0; i<NumElem; i++)
	{
		x1 = X[VertI[i]]; x2 = X[VertJ[i]]; x3 = X[VertK[i]];
		y1 = Y[VertI[i]]; y2 = Y[VertJ[i]]; y3 = Y[VertK[i]];
		Area[i] = 0.5*(x1*(y2 - y3) + x2*(y3 - y1) + x3*(y1 - y2));
		a1[i] = x2*y3 - x3*y2; b1[i] = y2 - y3; c1[i] = x3 - x2;
		a2[i] = x3*y1 - x1*y3; b2[i] = y3 - y1; c2[i] = x1 - x3;
		a3[i] = x1*y2 - x2*y1; b3[i] = y1 - y2; c3[i] = x2 - x1;

		YCij[i] = 2 * (-sigma[i] * Area[i] / 12) / Delt;
		YCi0[i] = 2 * (4 * sigma[i] * Area[i] / 12) / Delt;

		G01[i] = 1.0 / 4 * Vr[i] / Area[i] * (b1[i] * b2[i] + c1[i] * c2[i]);
		G20[i] = 1.0 / 4 * Vr[i] / Area[i] * (b1[i] * b3[i] + c1[i] * c3[i]);
		G12[i] = 1.0 / 4 * Vr[i] / Area[i] * (b3[i] * b2[i] + c3[i] * c2[i]);


		Kc[VertI[i]][VertI[i]] = Kc[VertI[i]][VertI[i]] + (-G01[i] - G20[i]) + YCij[i] + YCij[i] + YCi0[i];
		Kc[VertI[i]][VertJ[i]] = Kc[VertI[i]][VertJ[i]] + G01[i] - YCij[i];
		Kc[VertI[i]][VertK[i]] = Kc[VertI[i]][VertK[i]] + G20[i] - YCij[i];

		Kc[VertJ[i]][VertI[i]] = Kc[VertJ[i]][VertI[i]] + G01[i] - YCij[i];
		Kc[VertJ[i]][VertJ[i]] = Kc[VertJ[i]][VertJ[i]] + (-G01[i] - G12[i]) + YCij[i] + YCij[i] + YCi0[i];
		Kc[VertJ[i]][VertK[i]] = Kc[VertJ[i]][VertK[i]] + G12[i] - YCij[i];

		Kc[VertK[i]][VertI[i]] = Kc[VertK[i]][VertI[i]] + G20[i] - YCij[i];
		Kc[VertK[i]][VertJ[i]] = Kc[VertK[i]][VertJ[i]] + G12[i] - YCij[i];
		Kc[VertK[i]][VertK[i]] = Kc[VertK[i]][VertK[i]] + (-G20[i] - G12[i]) + YCij[i] + YCij[i] + YCi0[i];

	}

	for (i = 0; i<NumNodes; i++)
	{
		if (BoundaryNodeId[i] == 1)
		{
			for (j = 0; j<NumNodes; j++)
			{
				if (j == i)
					Kc[i][j] = 1;
				else
				{
					Kc[i][j] = 0;
					Kc[j][i] = 0;

				}
			}
		}
	}

	LU_De(Kc, NumNodes);

}
void FEM_Data_Prepare()
{
	int i,j,k=0;
	double x1,x2,x3,y1,y2,y3;
	static double Kc[MaxNodeNumber][MaxNodeNumber],Ks[MaxNodeNumber],Vr[MaxElemNumber];
	double Newton_Relative_Error=1;

	double Z0=1;
	double CenterX[MaxElemNumber],CenterY[MaxElemNumber];
	double* Paraout;
	
	VertI=(int*)malloc(MaxElemNumber*sizeof(int));
	VertJ=(int*)malloc(MaxElemNumber *sizeof(int));
	VertK=(int*)malloc(MaxElemNumber *sizeof(int));
	X=(double*)malloc(MaxNodeNumber *sizeof(double));
	Y=(double*)malloc(MaxNodeNumber *sizeof(double));
	//Load Mesh Information from file

	LoadMeshInfo(&NumNodes,&NumElem,&X,&Y,&VertI,&VertJ,&VertK);

	printf("Number of Nodes: %d   Number of Elements: %d\n",NumNodes,NumElem);
	

	Vs1=(double*)malloc(MaxStep*sizeof(double));
	Vs2=(double*)malloc(MaxStep*sizeof(double));
	Vs3=(double*)malloc(MaxStep*sizeof(double));

	MaxLength=(int)(NumNodes*NumNodes);	
	ValL=(double*)malloc(MaxLength*sizeof(double));
	ValU=(double*)malloc(MaxLength*sizeof(double));
	Ud=(double*)malloc(NumNodes*sizeof(double));

	ColL=(int*)malloc(MaxLength*sizeof(int));
	ColU=(int*)malloc(MaxLength*sizeof(int));
	RowL=(int*)malloc(MaxLength*sizeof(int));
	RowU=(int*)malloc(MaxLength*sizeof(int));

	Paraout=(double*)malloc(NumElem*sizeof(double));

	sigma=(double*)malloc(NumElem*sizeof(double));


	for(i=0;i<NumElem;i++)
	{
		V12_CIn[i]=0;
		V23_CIn[i]=0;
		V31_CIn[i]=0;

		V10_In[i]=0;
		V20_In[i]=0;
		V30_In[i]=0;
	}

		for(i=0;i<NumNodes;i++)
			for(j=0;j<NumNodes;j++)	
			{		
				Kc[i][j]=0;	
			}

			for (i=0;i<NumNodes;i++)
			{
				Ks[i]=0;//b vector
				BoundaryNodeId[i]=0;//Boundary Id all 0;
				BoundaryNodeValue[i]=0;//Boundary values initialized
			}


			for (i=0;i<MaxStep;i++)
			{
				Vpr[i]=0;
				Ipr[i]=0;
				Vse[i]=0;
				Ise[i]=0;
				Ikm[i]=0;
				Vkm[i]=0;
				Vs1[i]=VRMS*sin(i*Delt*2*pi*Freq);
			}

			//Boundary node ID and values
			for(i=0;i<NumNodes;i++)
			{
				if(fabs(X[i]-3.5)<1e-5)
				{
					BoundaryNodeId[i]=1;
					BoundaryNodeValue[i]=0;
				}
				if(fabs(X[i]+3.5)<1e-5)
				{
					BoundaryNodeId[i]=1;
					BoundaryNodeValue[i]=0;
				}
				if(fabs(Y[i]-2.5)<1e-5)
				{

					BoundaryNodeId[i]=1;
					BoundaryNodeValue[i]=0;
				}
				if(fabs(Y[i]+2.5)<1e-5)
				{
					BoundaryNodeId[i]=1;
					BoundaryNodeValue[i]=0;
				}
			}
			//Element Class  0: Air,  3: Iron,  x1 x2 x3 x4 Winding for x phase
			for (i=0;i<NumElem;i++)

			{
				ClassEle[i]=0;
				CenterX[i]=(X[VertI[i]]+X[VertJ[i]]+X[VertK[i]])/3.0;
				CenterY[i]=(Y[VertI[i]]+Y[VertJ[i]]+Y[VertK[i]])/3.0;

				if(CenterX[i]>-2.6 && CenterX[i]<2.6 && CenterY[i]>-1.8 && CenterY[i]<1.8)
				{ClassEle[i]=3;}

				if(CenterX[i]>-2.1 && CenterX[i]<-0.5 && CenterY[i]>-1.3 && CenterY[i]<1.3)
				{ClassEle[i]=0;}

				if(CenterX[i]>0.5 && CenterX[i]<2.1 && CenterY[i]>-1.3 && CenterY[i]<1.3)
				{ClassEle[i]=0;}

				if(CenterX[i]>-3.1 && CenterX[i]<-2.85 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=11;}

				if(CenterX[i]>-2.85 && CenterX[i]<-2.6 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=12;}

				if(CenterX[i]>-2.1 && CenterX[i]<-1.85 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=13;}

				if(CenterX[i]>-1.85 && CenterX[i]<-1.6 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=14;}

				if(CenterX[i]>-1 && CenterX[i]<-0.75 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=21;}

				if(CenterX[i]>-0.75 && CenterX[i]<-0.5 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=22;}

				if(CenterX[i]>0.5 && CenterX[i]<0.75 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=23;}

				if(CenterX[i]>0.75 && CenterX[i]<1 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=24;}

				if(CenterX[i]>1.6 && CenterX[i]<1.85 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=31;}

				if(CenterX[i]>1.85 && CenterX[i]<2.1 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=32;}

				if(CenterX[i]>2.6 && CenterX[i]<2.85 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=33;}

				if(CenterX[i]>2.85 && CenterX[i]<3.1 && CenterY[i]>-1 && CenterY[i]<1)
				{ClassEle[i]=34;}
				
			}

			NCoil1=0; NCoil2=0; NCoil3=0; NCoil4=0;
			NCoil5=0; NCoil6=0; NCoil7=0; NCoil8=0;
			NCoil9=0; NCoil10=0; NCoil11=0; NCoil12=0;
			//Find Coil Id for integration
			for(i=0;i<NumElem;i++)
			{
				if(ClassEle[i]==11){ Coil1Id[NCoil1]=i; NCoil1++; }
				if(ClassEle[i]==12){ Coil2Id[NCoil2]=i; NCoil2++; }
				if(ClassEle[i]==13){ Coil3Id[NCoil3]=i; NCoil3++; }
				if(ClassEle[i]==14){ Coil4Id[NCoil4]=i; NCoil4++; }
				if(ClassEle[i]==21){ Coil5Id[NCoil5]=i; NCoil5++; }
				if(ClassEle[i]==22){ Coil6Id[NCoil6]=i; NCoil6++; }
				if(ClassEle[i]==23){ Coil7Id[NCoil7]=i; NCoil7++; }
				if(ClassEle[i]==24){ Coil8Id[NCoil8]=i; NCoil8++; }
				if(ClassEle[i]==31){ Coil9Id[NCoil9]=i; NCoil9++; }
				if(ClassEle[i]==32){ Coil10Id[NCoil10]=i; NCoil10++; }
				if(ClassEle[i]==33){ Coil11Id[NCoil11]=i; NCoil11++; }
				if(ClassEle[i]==34){ Coil12Id[NCoil12]=i; NCoil12++; }
			}

			//LoadVe(Ve,NumElem);
			for (i=0;i<NumElem;i++)
			{
				if(ClassEle[i]==3)
				{
					Vr[i]=Reluct0;
					if(UseMemoryAdm==1){Vr[i]=Ve[i];}//material initialized
					sigma[i]=SigmaE;
				}
				else
				{
					Vr[i]=1/Mu0;
					sigma[i]=0;
				}

			}
			//Generate Memory admittance matrix

			//Get Stiffness Matrix and Admittance matrix

			for(i=0;i<NumElem;i++)
			{
				x1=X[VertI[i]]; x2=X[VertJ[i]];x3=X[VertK[i]];
				y1=Y[VertI[i]]; y2=Y[VertJ[i]];y3=Y[VertK[i]];
				Area[i]=0.5*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
				a1[i]=x2*y3-x3*y2; b1[i]=y2-y3; c1[i]=x3-x2;
				a2[i]=x3*y1-x1*y3; b2[i]=y3-y1; c2[i]=x1-x3;
				a3[i]=x1*y2-x2*y1; b3[i]=y1-y2; c3[i]=x2-x1;

				YCij[i]=2*(-sigma[i]*Area[i]/12)/Delt;
				YCi0[i]=2*(4*sigma[i]*Area[i]/12)/Delt;

				G01[i]=1.0/4*Vr[i]/Area[i]*(b1[i]*b2[i]+c1[i]*c2[i]);
				G20[i]=1.0/4*Vr[i]/Area[i]*(b1[i]*b3[i]+c1[i]*c3[i]);
				G12[i]=1.0/4*Vr[i]/Area[i]*(b3[i]*b2[i]+c3[i]*c2[i]);


				Kc[VertI[i]][VertI[i]]=Kc[VertI[i]][VertI[i]]+(-G01[i]-G20[i])+YCij[i]+YCij[i]+YCi0[i];
				Kc[VertI[i]][VertJ[i]]=Kc[VertI[i]][VertJ[i]]+G01[i]-YCij[i];
				Kc[VertI[i]][VertK[i]]=Kc[VertI[i]][VertK[i]]+G20[i]-YCij[i];

				Kc[VertJ[i]][VertI[i]]=Kc[VertJ[i]][VertI[i]]+G01[i]-YCij[i];
				Kc[VertJ[i]][VertJ[i]]=Kc[VertJ[i]][VertJ[i]]+(-G01[i]-G12[i])+YCij[i]+YCij[i]+YCi0[i];
				Kc[VertJ[i]][VertK[i]]=Kc[VertJ[i]][VertK[i]]+G12[i]-YCij[i];

				Kc[VertK[i]][VertI[i]]=Kc[VertK[i]][VertI[i]]+G20[i]-YCij[i];
				Kc[VertK[i]][VertJ[i]]=Kc[VertK[i]][VertJ[i]]+G12[i]-YCij[i];
				Kc[VertK[i]][VertK[i]]=Kc[VertK[i]][VertK[i]]+(-G20[i]-G12[i])+YCij[i]+YCij[i]+YCi0[i];

			}

			//Modify the stiffness matrix
			for(i=0;i<NumNodes;i++)
			{
				if(BoundaryNodeId[i]==1)
				{
					for(j=0;j<NumNodes;j++)
					{
						if(j==i)
							Kc[i][j]=1;
						else
						{
							Kc[i][j]=0;
							Kc[j][i]=0;

						}
					}
				}
			}

			LU_De(Kc,NumNodes);

}
void Inductance_Calculation(double Ip, double Is, double* Lw1,double*Lw2,double *Lm1,double *Lm2)
{
double sumA5=0,sumA6=0,sumA7=0,sumA8=0;
double Ac5,Ac6,Ac7,Ac8;
double Ep,Es,Ep1,Es1,Ep2,Es2;
int i = 0;

		FEM_TLM(0, 0, Ip, Is, 0, 0, A0);
		FEM_TLM(0, 0, Ip+dIp, Is, 0, 0, A0p);
		FEM_TLM(0, 0, Ip, Is+dIs, 0, 0, A0s);

		sumA5 = 0; sumA6 = 0; sumA7 = 0; sumA8 = 0;
		for (i = 0; i<NCoil5; i++)
			sumA5 = sumA5 + (A0[VertI[Coil5Id[i]]] + A0[VertJ[Coil5Id[i]]] + A0[VertK[Coil5Id[i]]]) / 3 * Area[Coil5Id[i]];
		Ac5 = sumA5 / CoilArea;

		for (i = 0; i<NCoil6; i++)
			sumA6 = sumA6 + (A0[VertI[Coil6Id[i]]] + A0[VertJ[Coil6Id[i]]] + A0[VertK[Coil6Id[i]]]) / 3 * Area[Coil6Id[i]];
		Ac6 = sumA6 / CoilArea;

		for (i = 0; i<NCoil7; i++)
			sumA7 = sumA7 + (A0[VertI[Coil7Id[i]]] + A0[VertJ[Coil7Id[i]]] + A0[VertK[Coil7Id[i]]]) / 3 * Area[Coil7Id[i]];
		Ac7 = sumA7 / CoilArea;

		for (i = 0; i<NCoil8; i++)
			sumA8 = sumA8 + (A0[VertI[Coil8Id[i]]] + A0[VertJ[Coil8Id[i]]] + A0[VertK[Coil8Id[i]]]) / 3 * Area[Coil8Id[i]];
		Ac8 = sumA8 / CoilArea;

		Ep = Np*Ac7 - Np*Ac6;
		Es = Ns*Ac5 - Ns*Ac8;

		sumA5 = 0; sumA6 = 0; sumA7 = 0; sumA8 = 0;
		for (i = 0; i<NCoil5; i++)
			sumA5 = sumA5 + (A0p[VertI[Coil5Id[i]]] + A0p[VertJ[Coil5Id[i]]] + A0p[VertK[Coil5Id[i]]]) / 3 * Area[Coil5Id[i]];
		Ac5 = sumA5 / CoilArea;

		for (i = 0; i<NCoil6; i++)
			sumA6 = sumA6 + (A0p[VertI[Coil6Id[i]]] + A0p[VertJ[Coil6Id[i]]] + A0p[VertK[Coil6Id[i]]]) / 3 * Area[Coil6Id[i]];
		Ac6 = sumA6 / CoilArea;

		for (i = 0; i<NCoil7; i++)
			sumA7 = sumA7 + (A0p[VertI[Coil7Id[i]]] + A0p[VertJ[Coil7Id[i]]] + A0p[VertK[Coil7Id[i]]]) / 3 * Area[Coil7Id[i]];
		Ac7 = sumA7 / CoilArea;

		for (i = 0; i<NCoil8; i++)
			sumA8 = sumA8 + (A0p[VertI[Coil8Id[i]]] + A0p[VertJ[Coil8Id[i]]] + A0p[VertK[Coil8Id[i]]]) / 3 * Area[Coil8Id[i]];
		Ac8 = sumA8 / CoilArea;

		Ep1 = Np*Ac7 - Np*Ac6;
		Es1 = Ns*Ac5 - Ns*Ac8;


		sumA5 = 0; sumA6 = 0; sumA7 = 0; sumA8 = 0;
		for (i = 0; i<NCoil5; i++)
			sumA5 = sumA5 + (A0s[VertI[Coil5Id[i]]] + A0s[VertJ[Coil5Id[i]]] + A0s[VertK[Coil5Id[i]]]) / 3 * Area[Coil5Id[i]];
		Ac5 = sumA5 / CoilArea;

		for (i = 0; i<NCoil6; i++)
			sumA6 = sumA6 + (A0s[VertI[Coil6Id[i]]] + A0s[VertJ[Coil6Id[i]]] + A0s[VertK[Coil6Id[i]]]) / 3 * Area[Coil6Id[i]];
		Ac6 = sumA6 / CoilArea;

		for (i = 0; i<NCoil7; i++)
			sumA7 = sumA7 + (A0s[VertI[Coil7Id[i]]] + A0s[VertJ[Coil7Id[i]]] + A0s[VertK[Coil7Id[i]]]) / 3 * Area[Coil7Id[i]];
		Ac7 = sumA7 / CoilArea;

		for (i = 0; i<NCoil8; i++)
			sumA8 = sumA8 + (A0s[VertI[Coil8Id[i]]] + A0s[VertJ[Coil8Id[i]]] + A0s[VertK[Coil8Id[i]]]) / 3 * Area[Coil8Id[i]];
		Ac8 = sumA8 / CoilArea;

		Ep2 = Np*Ac7 - Np*Ac6;
		Es2 = Ns*Ac5 - Ns*Ac8;

			*Lw1=(Ep1-Ep)/dIp;
			*Lw2=(Es2-Es)/dIs;
			*Lm1=(Es1-Es)/dIp;
			*Lm2=(Ep2-Ep)/dIs;
}
void FEM_TLM(double Ip1,double Is1,double Ip2,double Is2,double Ip3,double Is3,double *AA)
{
double Ke00,Ke01,Ke02,Ke10,Ke11,Ke12,Ke20,Ke21,Ke22,J0,J1,J2;
int i,k=0;
double Error_Temp,TLM_error=1;

static double Js[MaxElemNumber];
static double Ks_temp[MaxNodeNumber];
static double A[MaxNodeNumber];
static double Iteration_temp[MaxNodeNumber];  
static double V12_Re[MaxElemNumber];
static double V23_Re[MaxElemNumber];
static double V31_Re[MaxElemNumber];
static double V12_In[MaxElemNumber];
static double V23_In[MaxElemNumber];
static double V31_In[MaxElemNumber];

for(i=0;i<NumElem;i++)
{
	V12_In[i]=0;
	V23_In[i]=0;
	V31_In[i]=0;
}

	for(i=0;i<NumElem;i++)
		{
			if(ClassEle[i]==11)
				Js[i]=Is1*Ns/CoilArea;
			else if(ClassEle[i]==12)
				Js[i]=-Ip1*Np/CoilArea;
			else if(ClassEle[i]==13)
				Js[i]=Ip1*Np/CoilArea;
			else if(ClassEle[i]==14)
				Js[i]=-Is1*Ns/CoilArea;
			else if(ClassEle[i]==21)
				Js[i]=Is2*Ns/CoilArea;
			else if(ClassEle[i]==22)
				Js[i]=-Ip2*Np/CoilArea;
			else if(ClassEle[i]==23)
				Js[i]=Ip2*Np/CoilArea;
			else if(ClassEle[i]==24)
				Js[i]=-Is2*Ns/CoilArea;
			else if(ClassEle[i]==31)
				Js[i]=Is3*Ns/CoilArea;
			else if(ClassEle[i]==32)
				Js[i]=-Ip3*Np/CoilArea;
			else if(ClassEle[i]==33)
				Js[i]=Ip3*Np/CoilArea;
			else if(ClassEle[i]==34)
				Js[i]=-Is3*Ns/CoilArea;
			else
				Js[i]=0;
		}

		while(TLM_error>TLM_Tolerance)
			{

				for (i=0;i<NumNodes;i++)
					Ks_temp[i]=0;//b vector


				for(i=0;i<NumElem;i++)
				{	
					Ke00=Js[i]*Area[i]/3.0+2*V12_In[i]*(-G01[i])-2*V31_In[i]*(-G20[i])+2*V12_CIn[i]*YCij[i]-2*V31_CIn[i]*YCij[i]+2*V10_In[i]*YCi0[i];
					Ke01=Js[i]*Area[i]/3.0+2*V23_In[i]*(-G12[i])-2*V12_In[i]*(-G01[i])+2*V23_CIn[i]*YCij[i]-2*V12_CIn[i]*YCij[i]+2*V20_In[i]*YCi0[i];
					Ke02=Js[i]*Area[i]/3.0+2*V31_In[i]*(-G20[i])-2*V23_In[i]*(-G12[i])+2*V31_CIn[i]*YCij[i]-2*V23_CIn[i]*YCij[i]+2*V30_In[i]*YCi0[i];

					Ks_temp[VertI[i]]=Ks_temp[VertI[i]]+Ke00;
					Ks_temp[VertJ[i]]=Ks_temp[VertJ[i]]+Ke01;
					Ks_temp[VertK[i]]=Ks_temp[VertK[i]]+Ke02;
					
				}
				//Apply boundary
				for(i=0;i<NumNodes;i++)
				{
					if(BoundaryNodeId[i]==1)
						Ks_temp[i]=BoundaryNodeValue[i];
				}
				
				LU_SparseSolve(Ks_temp,A,NumNodes);	

				for(i=0;i<NumElem;i++)
				{	
					V12_Re[i]=A[VertI[i]]-A[VertJ[i]]-V12_In[i];
					V23_Re[i]=A[VertJ[i]]-A[VertK[i]]-V23_In[i];
					V31_Re[i]=A[VertK[i]]-A[VertI[i]]-V31_In[i];

					Ke01=1.0/4/Area[i]*(b1[i]*b2[i]+c1[i]*c2[i]);

					Ke12=1.0/4/Area[i]*(b3[i]*b2[i]+c3[i]*c2[i]);

					Ke20=1.0/4/Area[i]*(b1[i]*b3[i]+c1[i]*c3[i]);

					
					if(ClassEle[i]==3)
					{
						J0=(b1[i]*b2[i]+c1[i]*c2[i])/4/Area[i]/Area[i];
						J1=(b3[i]*b2[i]+c3[i]*c2[i])/4/Area[i]/Area[i];
						J2=(b1[i]*b3[i]+c1[i]*c3[i])/4/Area[i]/Area[i];

						Newton_Raphson_Element(&V12_In[i],&V23_In[i],&V31_In[i],V12_Re[i],V23_Re[i],V31_Re[i],J0,J1,J2,Ke01,Ke12,Ke20,i,Ve,&RecordB,&RecordH);

						
					}
					else 
					{
						V12_In[i]=0;
						V23_In[i]=0;
						V31_In[i]=0;

						
					
						
					}

				}
			
				for(i=0;i<NumNodes;i++)
				{
					if(BoundaryNodeId[i]!=1)
					{
						Error_Temp=Error_Temp+fabs((Iteration_temp[i]-A[i])/A[i]);
						Iteration_temp[i]=A[i];
					}
				}
				TLM_error=Error_Temp/NumNodes;
				Error_Temp=0;
				k++;
				//printf("TLM error for %d th iteration is %f\n",k,TLM_error);
			}
			printf("TLM iteration for %d th time step is %d\n", Time_Step, k);
			for(i=0;i<NumNodes;i++)
			AA[i]=A[i];
					
}
void LoadMeshInfo(int* NumNodes,int* NumElem, double** X,double** Y,int** VertI, int** VertJ, int** VertK)
{
	FILE* ip;
	int i,j,flag=0;
	char c;
	char filename[50];
	int line_size =50;
	char line[50];
	char temp[50];
	float x,y;
	sprintf(filename, "WidebandMeshInfo.mphtxt");

	if ((ip = fopen(filename,"r")) == NULL){
		printf("error opening the input data.\n");
	}  
	while (fgets(line, sizeof(line), ip)!= NULL)  {

		if(strstr(line,"# number of mesh points") !=NULL)
		{
			sscanf(line,"%d",&(*NumNodes));
		}
	}
	fclose(ip);


	if ((ip = fopen(filename,"r")) == NULL){
		printf("error opening the input data.\n");
	} 

	while (fgets(line, sizeof(line), ip)!= NULL)  {

		if(strstr(line,"# Mesh point coordinates") !=NULL)
		{
			for(i=0;i<*NumNodes;i++)
			{
				fgets(line, sizeof(line), ip);
				sscanf(line,"%lf %lf\n",&(*X)[i],&(*Y)[i]);	
			}
		}
	}
	fclose(ip);



	if ((ip = fopen(filename,"r")) == NULL){
		printf("error opening the input data.\n");
	} 

	while (fgets(line, sizeof(line), ip)!= NULL)  {

		if(strstr(line,"# number of elements") !=NULL)
		{
			flag=flag+1;
			if(flag==3)
			{
				sscanf(line,"%d",&(*NumElem));
				fgets(line,sizeof(line), ip);

				for(i=0;i<*NumElem;i++)
				{
					fgets(line, sizeof(line), ip);
					sscanf(line,"%d %d %d\n",&(*VertI)[i],&(*VertJ)[i],&(*VertK)[i]);
				}
			}
		}	
	}
	fclose(ip);
}
void LoadVe(double Ve[],int N)
{
	FILE* op;
	int i;

	if ((op = fopen("VectorVe.txt","r")) == NULL){
		printf("Error opening the output file.\n");
	}

	for (i = 0; i < N; i++)
		{	
			fscanf(op, "%lf ",&(Ve[i]));
			fscanf(op, "%\n");
	    }

	fclose(op);
}
void LU_De(double A[MaxNodeNumber][MaxNodeNumber], int n)
{
	int i,j,k;
	FILE *ip;
	static double AdmitL[MaxNodeNumber][MaxNodeNumber], AdmitU[MaxNodeNumber][MaxNodeNumber];

	/*AdmitL = (double**)malloc(NumNodes * sizeof(double*));
	for (i = 0; i<NumNodes; i++)
		AdmitL[i] = (double*)malloc(NumNodes * sizeof(double));

	AdmitU = (double**)malloc(NumNodes * sizeof(double*));
	for (i = 0; i<NumNodes; i++)
		AdmitU[i] = (double*)malloc(NumNodes * sizeof(double));*/

	for(j=0; j<n; j++)
	{
		for(i=0; i<n; i++)
		{
			if(i<=j)
			{
				AdmitU[i][j]=A[i][j];
				for(k=0; k<=i-1; k++)
					AdmitU[i][j]-= AdmitL[i][k]* AdmitU[k][j];
				if(i==j)
					AdmitL[i][j]=1;
				else
					AdmitL[i][j]=0;
			}
			else
			{
				AdmitL[i][j]=A[i][j];
				for(k=0; k<=j-1; k++)
					AdmitL[i][j]-= AdmitL[i][k]* AdmitU[k][j];
				AdmitL[i][j]/= AdmitU[j][j];
				AdmitU[i][j]=0;
			}
		}
	}

	for(j=0;j<n;j++)
		Ud[j]= AdmitU[j][j];

	NL=0;NU=0;
	for(j=0;j<n;j++){
		for(i=j+1;i<n;i++)
		{
			if (AdmitL[i][j]!=0)
			{
				ColL[NL]=j;
				RowL[NL]=i;
				ValL[NL]= AdmitL[i][j];
				NL++;
			}
		}
	}

	for(j=0;j<n;j++){
		for(i=0;i<j;i++)
		{
			if (AdmitU[i][j]!=0)
			{
				ColU[NU]=j;
				RowU[NU]=i;
				ValU[NU]= AdmitU[i][j]/ AdmitU[i][i];
				NU++;
			}
		}
	}


	if ((ip = fopen("LSparse.txt","w")) == NULL){
		printf("error opening the input data.\n");
		
	}
	fprintf(ip, "%d\n", NL);

	for (i = 0; i < NL; i++)
		fprintf(ip, "%d %d %.6E\n",RowL[i],ColL[i],ValL[i]);

	fprintf(ip, "\n");
	fclose(ip);


	if ((ip = fopen("USparse.txt","w")) == NULL){
		printf("error opening the input data.\n");
		 
	}
	fprintf(ip, "%d\n", NU);

	for (i = 0; i < NU; i++)
		fprintf(ip, "%d %d %.6E\n",RowU[i],ColU[i],ValU[i]);

	fprintf(ip, "\n");
	fclose(ip);



	if ((ip = fopen("Udiag.txt","w")) == NULL){
		printf("error opening the input data.\n");
		 
	}
	fprintf(ip, "%d\n",n);
	for (i = 0; i <n; i++)
		fprintf(ip, "%.6E\n",Ud[i]);

	fprintf(ip, "\n");
	fclose(ip);

	/*for (i = 0; i<NumNodes; i++)
	{
		free(AdmitL[i]);
		free(AdmitU[i]);
	}*/

}
void LU_Solve(double** L, double ** U,double* B, double** X,int n)
{
	double *Y=(double*)malloc(n*sizeof(double));
	int i,j;
	for(i=0; i<n; i++)
	{
		Y[i]=B[i];
		for(j=0; j<i; j++)
		{
			Y[i]-=L[i][j]*Y[j];
		}
	}

	for(i=n-1; i>=0; i--)
	{
		(*X)[i]= Y[i];
		for(j=i+1; j<n; j++)
		{
			(*X)[i]-=U[i][j]*(*X)[j];
		}
		(*X)[i]/=U[i][i];
	}


}
void Gauss_Solve(double** A, double **B,double** X, int n)
{
	int i,j,k,max_id;
	double temp,max_value;
	double *temp1;
	for(k=0;k<n;k++)
	{

		//Gauss Elimination
		for(i=k+1;i<n;i++)
		{
			temp=A[i][k]/A[k][k];
			for(j=k;j<n;j++)
				A[i][j]=A[i][j]-temp*A[k][j];
			if (j==n)
				(*B)[i]=(*B)[i]-temp*(*B)[k];
		}
	}
	//Jordan Elimination
	for(k=n-1;k>0;k--)
		for (i=0;i<k;i++)
		{
			(*B)[i]=(*B)[i]-(*B)[k]*A[i][k]/A[k][k];
			A[i][k]=0;
		}

		for(i=0;i<n;i++)
			(*B)[i]=(*B)[i]/A[i][i];

		*X=*B;
}
double* Cubic_Spline(double y[], double x[])
{
	int i,j,N;
	double *L,**A,*x1,*b;

	L=(double*) malloc(sizeof(double)*CubicN);
	x1=(double*) malloc(sizeof(double)*CubicN);
	b=(double*) malloc(sizeof(double)*CubicN);

	for(i=0;i<CubicN-1;i++)
	{
		L[i]=x[i+1]-x[i];
	}

	A=(double**)malloc(CubicN*sizeof(double*));//Time term Matrix
	for (i=0;i<CubicN;i++)
	{
		A[i]=(double*)malloc(CubicN*sizeof(double));
		b[i]=0;
	}

	for (i=0;i<CubicN;i++)
		for (j=0;j<CubicN;j++)
		{
			A[i][j]=0;
		}
		for (i=1;i<=CubicN-2;i++)
		{	
			A[i][i]=2/L[i]/L[i]+2/L[i-1]/L[i-1];
			A[i][i-1]=1/L[i-1]/L[i-1];
			A[i][i+1]=1/L[i]/L[i];

			b[i]=-3/L[i-1]/L[i-1]*y[i-1]+(3/L[i-1]/L[i-1]-3/L[i]/L[i])*y[i]+3/L[i]/L[i]*y[i+1];
		}

		A[0][0]=1;A[CubicN-1][CubicN-1]=1;
		b[0]=0;
		b[CubicN-1]=0;
		Gauss_Solve(A,&b,&x1,CubicN);
		return x1;
}
void Load_BH(void)
{
	CubicX[0]=1;CubicX[1]=2;CubicX[2]=3;
	CubicY[0]=2;CubicY[1]=3;CubicY[2]=5;
	CubicX1=Cubic_Spline(CubicY,CubicX);
}
interp_t Interpo(double X)
{
	int i,id;
	interp_t T;
	double x;
	for(i=0;i<CubicN-1;i++)
	{
		if(X>=CubicX[i] && X<=CubicX[i+1])
			id =i;
	}
	if(X<CubicX[0])
	{id=0;X=CubicX[0];}
	if(X>CubicX[CubicN-1])
	{id=CubicN-1;X=CubicX[CubicN-1];}

	x=(X-CubicX[id])/(CubicX[id+1]-CubicX[id]);
	T.Y=(2*x*x*x-3*x*x+1)*CubicY[id]+(-2*x*x*x+3*x*x)*CubicY[id+1]+(x*x*x-2*x*x+x)*CubicX1[id]+(x*x*x-x*x)*CubicX1[id+1];
	T.Y1=1/(CubicX[id+1]-CubicX[id])*((6*x*x-6*x)*CubicY[id]+(-6*x*x+6*x)*CubicY[id+1]+(3*x*x-4*x+1)*CubicX1[id]+(3*x*x-2*x)*CubicX1[id+1]);
	T.Y2=1/(CubicX[id+1]-CubicX[id])*x*((12*x-6)*CubicY[id]+(-12*x+6)*CubicY[id+1]+(6*x-4)*CubicX1[id]+(6*x-2)*CubicX1[id+1]);
	return T;
}
int saveoutput(double *X,double *Y,double** S, int N, int Count){
	FILE* op;
	int i, j;

	if ((op = fopen("FEM_TLM_Solution.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	fprintf(op, "%d\n",N);
	for (i = 0; i < N; i++)
	{
		fprintf(op, "%lf %lf ",X[i],Y[i]);
		for(j=0;j<Count;j++)
		fprintf(op, "%lf ",S[j][i]*1e6);

		fprintf(op, "\n");
	}
	fclose(op);
	return 0;
}
int saveoutBH(double **B,double **H, int id, int count){
	FILE* op;
	int i, j;

	if ((op = fopen("RecordB.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	for (i = 0; i < id; i++)
	{
		for(j=0;j<count;j++)
			fprintf(op, "%lf ",B[j][i]);
		fprintf(op, "\n");
	}
	fclose(op);

	if ((op = fopen("RecordH.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	for (i = 0; i < id; i++)
	{
		for(j=0;j<count;j++)
			fprintf(op, "%lf ",H[j][i]);
		fprintf(op, "\n");
	}
	fclose(op);

	return 0;
}
int Newton_Raphson_Element(double *X,double *Y,double *Z,double x0, double y0, double z0, double J0,double J1,double J2,double Ke01,double Ke12,double Ke20,int Id,double Ve[],double*RecB,double*RecH)
{
	double K00,K01,K02,K10,K11,K12,K20,K21,K22;
	double b0,b1,b2;
	double x=0, y=0, z=0,dx,dy,dz,Error=1,xtemp,ytemp,ztemp;
	double B2;
	double Vb2;
	double B;
	double Vi,dB2x,dB2y,dB2z;
	int count=0;
	if(x0==0 && y0==0 && z0==0)
	{
		x=0;y=0;z=0;
	}
	else
	{
		while(Error>1e-6 && count<50)
		//while(count<10)
	{
		xtemp=x;ytemp=y;ztemp=z;
		B2=-J0*(x+x0)*(x+x0)-J1*(y+y0)*(y+y0)-J2*(z+z0)*(z+z0);
		B=sqrt(B2);

		/*if(B<1.3)
		{
			Vi=15;
			Vb2=0;
			x=x0*(1-1/G01[Id]*Vi*Ke01)/(1+1/G01[Id]*Vi*Ke01);
			y=y0*(1-1/G12[Id]*Vi*Ke12)/(1+1/G12[Id]*Vi*Ke12);
			z=z0*(1-1/G20[Id]*Vi*Ke20)/(1+1/G20[Id]*Vi*Ke20);
		}
		else
		{
			Vi=15+1e5*pow(B-1.3,8)/B;
			Vb2=(B*8e5*pow(B-1.3,7)-1e5*pow(B-1.3,8))/B/B/2/B;
			b0=-(-Ke01*Vi*(x+x0)-(x0-x)*(-G01[Id]));
			b1=-(-Ke12*Vi*(y+y0)-(y0-y)*(-G12[Id]));
			b2=-(-Ke20*Vi*(z+z0)-(z0-z)*(-G20[Id]));
			dB2x=-2*J0*(x+x0);dB2y=-2*J1*(y+y0);dB2z=-2*J2*(z+z0);
			K00=-Ke01*(x+x0)*Vb2*dB2x-Ke01*Vi+(-G01[Id]);
			K01=-Ke01*Vb2*dB2y*(x+x0);
			K02=-Ke01*Vb2*dB2z*(x+x0);

			K10=-Ke12*Vb2*dB2x*(y+y0);
			K11=-Ke12*(y+y0)*Vb2*dB2y-Ke12*Vi+(-G12[Id]);
			K12=-Ke12*Vb2*dB2z*(y+y0);

			K20=-Ke20*Vb2*dB2x*(z+z0);
			K21=-Ke20*Vb2*dB2y*(z+z0);
			K22=-Ke20*(z+z0)*Vb2*dB2z-Ke20*Vi+(-G20[Id]);

			dx=Det(b0,K01,K02,b1,K11,K12,b2,K21,K22)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);
			dy=Det(K00,b0,K02,K10,b1,K12,K20,b2,K22)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);
			dz=Det(K00,K01,b0,K10,K11,b1,K20,K21,b2)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);
			x=x+1.2*dx;
			y=y+1.2*dy;
			z=z+1.2*dz;
			}*/
		/*Vi=2500+1000*B-3000*B*B+2100*B*B*B;
		Vb2=(1000-6000*B+3*2100*B*B)/2/B;*/
		

		if (B<0.6)
		{
			Vi = Reluct0;
			Vb2 = 0;
		}
		else
		{	
			Vi = Reluct0 + 3000*pow(B - 0.6, 3) / B;
			Vb2 = (B*9000*pow(B - 0.6, 2) - 3000*pow(B - 0.6, 3)) / B / B / 2 / B;
		}
	    
		b0=-(-Ke01*Vi*(x+x0)-(x0-x)*(-G01[Id]));
		b1=-(-Ke12*Vi*(y+y0)-(y0-y)*(-G12[Id]));
		b2=-(-Ke20*Vi*(z+z0)-(z0-z)*(-G20[Id]));
		dB2x=-2*J0*(x+x0);dB2y=-2*J1*(y+y0);dB2z=-2*J2*(z+z0);
		K00=-Ke01*(x+x0)*Vb2*dB2x-Ke01*Vi+(-G01[Id]);
		K01=-Ke01*Vb2*dB2y*(x+x0);
		K02=-Ke01*Vb2*dB2z*(x+x0);

		K10=-Ke12*Vb2*dB2x*(y+y0);
		K11=-Ke12*(y+y0)*Vb2*dB2y-Ke12*Vi+(-G12[Id]);
		K12=-Ke12*Vb2*dB2z*(y+y0);

		K20=-Ke20*Vb2*dB2x*(z+z0);
		K21=-Ke20*Vb2*dB2y*(z+z0);
		K22=-Ke20*(z+z0)*Vb2*dB2z-Ke20*Vi+(-G20[Id]);

		dx=Det(b0,K01,K02,b1,K11,K12,b2,K21,K22)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);
		dy=Det(K00,b0,K02,K10,b1,K12,K20,b2,K22)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);
		dz=Det(K00,K01,b0,K10,K11,b1,K20,K21,b2)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);

		x=x+1.2*dx;
		y=y+1.2*dy;
		z=z+1.2*dz;


		Ve[Id]=Vi;
		if(x!=0 && y!=0 && z!=0)
			Error=(fabs((x-xtemp)/x)+fabs((y-ytemp)/y)+fabs((z-ztemp)/z));

		if(fabs(x)<1e-15 && fabs(y)<1e-15 && fabs(z)<1e-15)
		break;

		count++;
	}
	}

	*X=x;
	*Y=y;
	*Z=z;
	*RecB=B;
	*RecH=Vi*B;

	return 0;
}
double Det(double x0,double x1,double x2,double y0,double y1,double y2,double z0,double z1,double z2)
{
	double T=0;
	T=x0*y1*z2+x1*y2*z0+x2*y0*z1-x0*y2*z1-x1*y0*z2-x2*y1*z0;
	return T;
}
int saveoutAxb(double **K,double **L,double **U, int N){
	FILE* op;
	int i, j;

	if ((op = fopen("K.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	fprintf(op, "%d\n",N);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++) 
			fprintf(op, "%lf ",K[i][j]);
		fprintf(op, "\n");
	}
	fprintf(op, "\n");
	fprintf(op, "\n");
	fclose(op);

	if ((op = fopen("L.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	fprintf(op, "%d\n",N);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++) 
			fprintf(op, "%lf ",L[i][j]);
		fprintf(op, "\n");
	}
	fprintf(op, "\n");
	fprintf(op, "\n");
	fclose(op);

	if ((op = fopen("U.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	fprintf(op, "%d\n",N);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++) 
			if(j!=i)
				fprintf(op, "%lf ",U[i][j]);
			else
				fprintf(op, "%lf ",U[i][j]);

		fprintf(op, "\n");
	}
	fprintf(op, "\n");
	fprintf(op, "\n");
	fclose(op);


	return 0;
}
void LU_SparseSolve(double Ks_temp[MaxNodeNumber],double A[MaxNodeNumber],int N)
{
	int i;
	for(i=0;i<N;i++)
		A[i]=Ks_temp[i];

	for(i=0;i<NL;i++)
		A[RowL[i]]=A[RowL[i]]-A[ColL[i]]*ValL[i];


	for(i=0;i<N;i++)
		A[i]=A[i]/Ud[i];

	for(i=NU-1;i>=0;i--)
		A[RowU[i]]=A[RowU[i]]-A[ColU[i]]*ValU[i];
}
void saveoutV(double V[],int N)
{
	FILE* ip;
	int i;

	if ((ip = fopen("VectorVeout.txt", "w")) == NULL) {
		printf("error opening the input data.\n");

	}
	for (i = 0; i < N; i++)
		fprintf(ip, "%f\n",V[i]);

	fclose(ip);
}
void saveoutIVIV(double* I1, double *V1, double* I2,double* V2,int N)
{
	FILE* ip;
	int i;

	if ((ip = fopen("IVIV.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	for (i = 0; i < N; i++)
	{
		
		fprintf(ip, "%f %f %f %f\n",I1[i],V1[i],I2[i],V2[i]);

	}

	fclose(ip);
}
void saveoutL(double* L1, double *L2, double* M12,double* M21,int N)
{
	FILE* op;
	int i;

	if ((op = fopen("SaveL.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	for (i = 0; i < N; i++)
	{

		fprintf(op, "%f %f %f %f\n",L1[i],L2[i],M12[i],M21[i]);

	}

	fclose(op);
}
void saveoutPara(double* Para,int N)
{
	FILE* op;
	int i;
	int k;

	if ((op = fopen("A0.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	for (i = 0; i < N; i++)
	{
		fprintf(op, "%.6E %.6E %.6E\n",X[i],Y[i],A0[i]);
	}


	fclose(op);


	if ((op = fopen("A0p.txt", "w")) == NULL) {
		printf("Error opening the output file.\n");

	}

	for (i = 0; i < N; i++)
	{
		fprintf(op, "%.6E %.6E %.6E\n", X[i], Y[i], A0p[i]);
	}


	fclose(op);


	if ((op = fopen("A0s.txt", "w")) == NULL) {
		printf("Error opening the output file.\n");

	}

	for (i = 0; i < N; i++)
	{
		fprintf(op, "%.6E %.6E %.6E\n", X[i], Y[i], A0s[i]);
	}


	fclose(op);
}
void saveoutParaInt(int* Para,int N)
{
	FILE* op;
	int i;
	int k;

	if ((op = fopen("ParaoutInt.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		 
	}

	fprintf(op, "{");
	for (i = 0; i < N; i++)
	{
		fprintf(op, "%d,",Para[i]);
	}
	fprintf(op, "}");

	fclose(op);
}

void SM(double dt,double *Vc, double *Vzu2, double *Vau2, double Iu, double Id, double Uref, double *Vp, double *Vp_up, double *Vp_dn, double *vi, double *car)
{
	double C = 6e-3, ZC = dt / (2 * C);
	double E = 180000.0, Vcref = E / NoC;
	double G11, G12, G22, f1, f2, Vnew1, Vnew2, Js;
	double K1 = 0.5, K2 = 150, K3 = 1.5, K4 = 150;

	double Vcsum = 0;
	for (int k = 0; k<NoM; k++)
	{
		Vcsum = Vcsum + Vc[k];
	}
	double Vcave = Vcsum / NoM;

	double Vzu1 = (Vcref - Vcave)*K1;
	*Vzu2 = *Vzu2 + (Vcref - Vcave)*dt*K2;
	double Vzu = Vzu1 + *Vzu2;
	double Vzuref = (Iu + Id) / 2;
	double Vau1 = (Vzuref - Vzu)*K3;
	*Vau2 = *Vau2 + (Vzuref - Vzu)*dt*K4;
	double Vau = Vau1 + *Vau2;

	double K5, K6;
	if (Iu >= 0)
		K5 = 0.35;
	else
		K5 = -0.35;

	double Bu[NoM];
	for (int i = 0; i<NoC; i++)
	{
		Bu[i] = K5*(Vcref - Vc[i]);
	}

	if (Id >= 0)
		K6 = 0.35;
	else
		K6 = -0.35;

	for (int i = NoC; i<NoM; i++)
	{
		Bu[i] = K6*(Vcref - Vc[i]);
	}

	double m[NoM];
	for (int i = 0; i<NoC; i++)
	{
		m[i] = (Vau + Bu[i] + 0.5*Vcref - Uref) / Vcref;
	}
	for (int i = NoC; i<NoM; i++)
	{
		m[i] = (Vau + Bu[i] + 0.5*Vcref + Uref) / Vcref;
	}

	double Vg[NoM], Gup[NoM], Gdn[NoM];
	for (int i = 0; i<NoC; i++)
	{
		if (m[i]>car[i])
		{
			Vg[i] = 10;
		}
		else
		{
			Vg[i] = -10;
		}
	}

	for (int i = NoC; i<NoM; i++)
	{
		if (m[i]>car[i - NoC])
		{
			Vg[i] = 10;
		}
		else
		{
			Vg[i] = -10;
		}
	}

	for (int i = 0; i<NoM; i++)
	{
		if (Vg[i]>5)
		{
			Gup[i] = 1e3;
			Gdn[i] = 1e-6;
		}
		else
		{
			Gup[i] = 1e-6;
			Gdn[i] = 1e3;
		}
	}

	for (int i = 0; i<NoC; i++)
	{
		G11 = 1 / ZC + Gup[i];
		G12 = -Gup[i];
		G22 = Gdn[i] + Gup[i];
		f1 = 2 * vi[i] / ZC;
		f2 = Iu;
		Vnew1 = (G12*f2 - G22*f1) / (G12*G12 - G11*G22);
		Vnew2 = (G12*f1 - G11*f2) / (G12*G12 - G11*G22);
		vi[i] = Vnew1 - vi[i];
		Vc[i] = Vnew1;
		Vp[i] = Vnew2;
	}
	for (int i = NoC; i<NoM; i++)
	{
		G11 = 1 / ZC + Gup[i];
		G12 = -Gup[i];
		G22 = Gdn[i] + Gup[i];
		f1 = 2 * vi[i] / ZC;
		f2 = Id;
		Vnew1 = (G12*f2 - G22*f1) / (G12*G12 - G11*G22);
		Vnew2 = (G12*f1 - G11*f2) / (G12*G12 - G11*G22);
		vi[i] = Vnew1 - vi[i];
		Vc[i] = Vnew1;
		Vp[i] = Vnew2;
	}

	*Vp_up = 0;
	*Vp_dn = 0;

	for (int i = 0; i<NoC; i++)
	{
		*Vp_up = *Vp_up + Vp[i];
	}
	for (int i = NoC; i<NoM; i++)
	{
		*Vp_dn = *Vp_dn + Vp[i];
	}

}