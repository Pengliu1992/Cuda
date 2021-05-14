#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define SigmaE  0
#define Mu0 1.256637e-6
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
void LoadMeshInfo(int* NumNodes,int* NumElem,double** X,double** Y,int** VertI, int** VertJ, int** VertK);
void LU_De(double** A, int n, double*** L, double*** U);
void LU_Solve(double** L, double ** U,double* B, double** X,int n);
void Gauss_Solve(double*** A, double **B,double** X, int n);
void Load_BH(void);
int Newton_Raphson_Element(double *X,double *Y,double *Z,double x0, double y0, double z0, double J0,double J1,double J2,double Ke01,double Ke12,double Ke20,double Z0);
double Det(double x0,double x1,double x2,double y0,double y1,double y2,double z0,double z1,double z2);
int saveoutput(double *X,double *Y,double* S, int N);
double* Cubic_Spline(double y[], double x[]);
interp_t Interpo(double x);
//define BH info.
const int CubicN=3;
double CubicX[3], CubicY[3];
double *CubicX1;

int main(){
	const int MaxNodes=5000;
	const int MaxElem=MaxNodes*2;
	
	int i,j,k=0;
	int NumNodes, NumElem;
	int *VertI=(int*)malloc(MaxElem*sizeof(int)),*VertJ=(int*)malloc(MaxElem*sizeof(int)),*VertK=(int*)malloc(MaxElem*sizeof(int));
	double *X=(double*)malloc(MaxNodes*sizeof(double)), *Y=(double*)malloc(MaxNodes*sizeof(double));
	double *Area,*a1,*b1,*c1,*a2,*b2,*c2,*a3,*b3,*c3;
	double x1,x2,x3,y1,y2,y3,Ke00,Ke01,Ke02,Ke10,Ke11,Ke12,Ke20,Ke21,Ke22,J0,J1,J2,J3,J4;
	double **Kc,**Kt,*Ks,*Js,*Vr,*A,*DeltA,**Jaco,*VB2,*JacoB,**AdmitL,**AdmitU,**Kc_temp,*Ks_temp,**Kc_temp1;
	double T_start,T_end,interp,Newton_Relative_Error=1,Error_Temp;
	interp_t Test;
	double *V1,*V2,*V3,*V12,*V23,*V31,*V12_Re,*V12_In,*V23_Re,*V23_In,*V31_Re,*V31_In;
	int *BoundaryNodeId;
	double *BoundaryNodeValue;
	double Z0=10,TLM_error=1,*Iteration_temp;
//Load Mesh Information from file
	for (i=0;i<NumNodes;i++)
	{
		VertI[i]=0;
		VertJ[i]=0;
		VertK[i]=0;
	}
	LoadMeshInfo(&NumNodes,&NumElem,&X,&Y,&VertI,&VertJ,&VertK);
	printf("Number of Nodes: %d   Number of Elements: %d\n",NumNodes,NumElem);
//Malloc and Initialize
																																																																									CubicX1=(double*)malloc(CubicN*sizeof(double));
	Area=(double*)malloc(NumElem*sizeof(double));
	a1=(double*)malloc(NumElem*sizeof(double));
	b1=(double*)malloc(NumElem*sizeof(double));
	c1=(double*)malloc(NumElem*sizeof(double));
	a2=(double*)malloc(NumElem*sizeof(double));
	b2=(double*)malloc(NumElem*sizeof(double));
	c2=(double*)malloc(NumElem*sizeof(double));
	a3=(double*)malloc(NumElem*sizeof(double));
	b3=(double*)malloc(NumElem*sizeof(double));
	c3=(double*)malloc(NumElem*sizeof(double));
	V1=(double*)malloc(NumElem*sizeof(double));
	V2=(double*)malloc(NumElem*sizeof(double));
	V3=(double*)malloc(NumElem*sizeof(double));
	V12=(double*)malloc(NumElem*sizeof(double));
	V23=(double*)malloc(NumElem*sizeof(double));
	V31=(double*)malloc(NumElem*sizeof(double));
	V12_Re=(double*)malloc(NumElem*sizeof(double));
	V23_Re=(double*)malloc(NumElem*sizeof(double));
	V31_Re=(double*)malloc(NumElem*sizeof(double));
	V12_In=(double*)malloc(NumElem*sizeof(double));
	V23_In=(double*)malloc(NumElem*sizeof(double));
	V31_In=(double*)malloc(NumElem*sizeof(double));


	Kc=(double**)malloc(NumNodes*sizeof(double*));//Space term Matrix
	for (i=0;i<NumNodes;i++)
		Kc[i]=(double*)malloc(NumNodes*sizeof(double));
	Kc_temp=(double**)malloc(NumNodes*sizeof(double*));//Space term Matrix
	for (i=0;i<NumNodes;i++)
		Kc_temp[i]=(double*)malloc(NumNodes*sizeof(double));

	Kc_temp1=(double**)malloc(NumNodes*sizeof(double*));//Space term Matrix
	for (i=0;i<NumNodes;i++)
		Kc_temp1[i]=(double*)malloc(NumNodes*sizeof(double));

	Jaco=(double**)malloc(NumNodes*sizeof(double*));//Jacobian Matrix
	for (i=0;i<NumNodes;i++)
		Jaco[i]=(double*)malloc(NumNodes*sizeof(double));

	Kt=(double**)malloc(NumNodes*sizeof(double*));//Time term Matrix
	for (i=0;i<NumNodes;i++)
		Kt[i]=(double*)malloc(NumNodes*sizeof(double));

	AdmitL=(double**)malloc(NumNodes*sizeof(double*));//Time term Matrix
	for (i=0;i<NumNodes;i++)
		AdmitL[i]=(double*)malloc(NumNodes*sizeof(double));

	AdmitU=(double**)malloc(NumNodes*sizeof(double*));//Time term Matrix
	for (i=0;i<NumNodes;i++)
	AdmitU[i]=(double*)malloc(NumNodes*sizeof(double));

	BoundaryNodeValue=(double*)malloc(NumNodes*sizeof(double));
	BoundaryNodeId=(int*)malloc(NumNodes*sizeof(int));
	A=(double*)malloc(NumNodes*sizeof(double));  //Final Node Solution
	Iteration_temp=(double*)malloc(NumNodes*sizeof(double));  
	Ks_temp=(double*)malloc(NumNodes*sizeof(double));  //Final Node Solution
	JacoB=(double*)malloc(NumNodes*sizeof(double));  //NR iteration RHS vector
	DeltA=(double*)malloc(NumNodes*sizeof(double));//Increment Between Iterations
	Ks=(double*)malloc(NumNodes*sizeof(double)); //RHS vector
	Js=(double*)malloc(NumElem*sizeof(double)); //Current Density
	Vr=(double*)malloc(NumElem*sizeof(double)); //Material reluctivity
	VB2=(double*)malloc(NumElem*sizeof(double));//Derivative
	for(i=0;i<NumNodes;i++)
		for(j=0;j<NumNodes;j++)	
		{		
			Kc[i][j]=0;	
			Kc_temp[i][j]=0;	
			Kt[i][j]=0;
			AdmitL[i][j]=0;
			AdmitU[i][j]=0;
		}

	for (i=0;i<NumNodes;i++)
	{
		Ks[i]=0;//b vector
		A[i]=1; //A vector initialized
		BoundaryNodeId[i]=0;//Boundary Id all 0;
		BoundaryNodeValue[i]=0;//Boundary values initialized
	}

	

	for (i=0;i<NumElem;i++)
		
	{
		Vr[i]=1;//material initialized
		VB2[i]=0.1;//Vr=0.1B
	}
	for (i=0;i<NumElem;i++)
		Js[i]=-1;//current	
	//Apply Boundary Conditions
		for(i=0;i<NumNodes;i++)
	{
		if(X[i]==0)
		{
			BoundaryNodeId[i]=1;
			BoundaryNodeValue[i]=X[i]+Y[i];
		}
		if(X[i]==1)
		{
			BoundaryNodeId[i]=1;
			BoundaryNodeValue[i]=X[i]+Y[i];
		}
		if(Y[i]==0)
		{
			BoundaryNodeId[i]=1;
			BoundaryNodeValue[i]=X[i]+Y[i];
		}
		if(Y[i]==1)
		{
			BoundaryNodeId[i]=1;
			BoundaryNodeValue[i]=X[i]+Y[i];
		}
	}
		
//Get Stiffness Matrix and Admittance matrix
	GET_TIME(T_start);

	for(i=0;i<NumElem;i++)
	{
		x1=X[VertI[i]]; x2=X[VertJ[i]];x3=X[VertK[i]];
		y1=Y[VertI[i]]; y2=Y[VertJ[i]];y3=Y[VertK[i]];
		Area[i]=0.5*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
		a1[i]=x2*y3-x3*y2; b1[i]=y2-y3; c1[i]=x3-x2;
		a2[i]=x3*y1-x1*y3; b2[i]=y3-y1; c2[i]=x1-x3;
		a3[i]=x1*y2-x2*y1; b3[i]=y1-y2; c3[i]=x2-x1;

		Ke00=2/Z0;
		Ke01=-1/Z0;
		Ke02=-1/Z0;

		Ke10=-1/Z0;
		Ke11=2/Z0;
		Ke12=-1/Z0;

		Ke20=-1/Z0;
		Ke21=-1/Z0;
		Ke22=2/Z0;

		Kc[VertI[i]][VertI[i]]=Kc[VertI[i]][VertI[i]]+Ke00;
		Kc[VertI[i]][VertJ[i]]=Kc[VertI[i]][VertJ[i]]+Ke01;
		Kc[VertI[i]][VertK[i]]=Kc[VertI[i]][VertK[i]]+Ke02;

		Kc[VertJ[i]][VertI[i]]=Kc[VertJ[i]][VertI[i]]+Ke10;
		Kc[VertJ[i]][VertJ[i]]=Kc[VertJ[i]][VertJ[i]]+Ke11;
		Kc[VertJ[i]][VertK[i]]=Kc[VertJ[i]][VertK[i]]+Ke12;

		Kc[VertK[i]][VertI[i]]=Kc[VertK[i]][VertI[i]]+Ke20;
		Kc[VertK[i]][VertJ[i]]=Kc[VertK[i]][VertJ[i]]+Ke21;
		Kc[VertK[i]][VertK[i]]=Kc[VertK[i]][VertK[i]]+Ke22;

		Ke00=Js[i]*Area[i]/3.0; 
		Ke01=Js[i]*Area[i]/3.0;
		Ke02=Js[i]*Area[i]/3.0;

		Ks[VertI[i]]=Ks[VertI[i]]+Ke00;
		Ks[VertJ[i]]=Ks[VertJ[i]]+Ke01;
		Ks[VertK[i]]=Ks[VertK[i]]+Ke02;
	}
	for(i=0;i<NumNodes;i++)
		for(j=0;j<NumNodes;j++)	
			Kc_temp[i][j]=Kc[i][j];	

	for(i=0;i<NumNodes;i++)
		for(j=0;j<NumNodes;j++)	
			Kc_temp1[i][j]=Kc[i][j];	
//Modify the stiffness matrix
	for(i=0;i<NumNodes;i++)
	{
		if(BoundaryNodeId[i]==1)
		{
			for(j=0;j<NumNodes;j++)
			{
				if(j==i)
					Kc_temp[i][j]=1;
				else
				{
					Kc_temp[i][j]=0;
					Kc_temp[j][i]=0;

				}
			}
		}
	}

	LU_De(Kc_temp,NumNodes,&AdmitL,&AdmitU);

	

	GET_TIME(T_end);
	printf("\nElapsed time for assembling is %.8f\n\n",(T_end-T_start));
	//Newton-Raphson Iteration 
	GET_TIME(T_start);
		
		for(i=0;i<NumElem;i++)
		{
			V12_In[i]=0.1;
			V23_In[i]=0.1;
			V31_In[i]=0.1;
		}

//Initialize TLM 
while(TLM_error>1e-6)
{
	for (i=0;i<NumNodes;i++)
	Ks_temp[i]=0;//b vector

	

	for(i=0;i<NumElem;i++)
	{		
		Ke01=1.0/4/Area[i]*(b1[i]*b2[i]+c1[i]*c2[i]);
		
		Ke12=1.0/4/Area[i]*(b3[i]*b2[i]+c3[i]*c2[i]);

		Ke20=1.0/4/Area[i]*(b1[i]*b3[i]+c1[i]*c3[i]);
		
		
		J0=0.1*(b1[i]*b2[i]+c1[i]*c2[i])/4/Area[i]/Area[i];
		J1=0.1*(b3[i]*b2[i]+c3[i]*c2[i])/4/Area[i]/Area[i];
		J2=0.1*(b1[i]*b3[i]+c1[i]*c3[i])/4/Area[i]/Area[i];

		V12_Re[i]=A[VertI[i]]-A[VertJ[i]]-V12_In[i];
		V23_Re[i]=A[VertJ[i]]-A[VertK[i]]-V23_In[i];
		V31_Re[i]=A[VertK[i]]-A[VertI[i]]-V31_In[i];

		Newton_Raphson_Element(&V12_In[i],&V23_In[i],&V31_In[i],V12_Re[i],V23_Re[i],V31_Re[i],J0,J1,J2,Ke01,Ke12,Ke20,Z0);
		
		/*	Vr[i]=0.1*(J3*J3+J4*J4)/4/Area[i]/Area[i];

		V12_In[i]=V12_Re[i]*(1+Z0*Vr[i]*Ke01)/(1-Z0*Vr[i]*Ke01);
		V23_In[i]=V23_Re[i]*(1+Z0*Vr[i]*Ke12)/(1-Z0*Vr[i]*Ke12);
		V31_In[i]=V31_Re[i]*(1+Z0*Vr[i]*Ke20)/(1-Z0*Vr[i]*Ke20);*/

		Ke00=Js[i]*Area[i]/3.0+2*V12_In[i]/Z0-2*V31_In[i]/Z0; 
		Ke01=Js[i]*Area[i]/3.0+2*V23_In[i]/Z0-2*V12_In[i]/Z0;
		Ke02=Js[i]*Area[i]/3.0+2*V31_In[i]/Z0-2*V23_In[i]/Z0;
		
		
		Ks_temp[VertI[i]]=Ks_temp[VertI[i]]+Ke00;
		Ks_temp[VertJ[i]]=Ks_temp[VertJ[i]]+Ke01;
		Ks_temp[VertK[i]]=Ks_temp[VertK[i]]+Ke02;
		
		
	}

	
	//Apply boundary


	for(i=0;i<NumNodes;i++)
		if(BoundaryNodeId[i]==1)
			Ks_temp[i]=BoundaryNodeValue[i];

	for(i=0;i<NumNodes;i++)
		if(BoundaryNodeId[i]==0)
			for(j=0;j<NumNodes;j++)
				if(BoundaryNodeId[j]==1)
					Ks_temp[i]=Ks_temp[i]-Ks_temp[j]*Kc_temp1[j][i];

	
	LU_Solve(AdmitL,AdmitU,Ks_temp,&A,NumNodes);
	
	
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

	
}
	
	GET_TIME(T_end);
	

	printf("\nNode Solution is: \n");
	for(j=0;j<NumNodes;j++)
	{printf("X: %f   Y: %f   S= %f ",X[j],Y[j],A[j]);
	printf("\n");}

	printf("\nElapsed time for Newton_Raphson Iteration is %.8f\n\n",(T_end-T_start));
	//Load_BH();
	//saveoutput(X,Y,A,NumNodes);
	//Test=Interpo(1.5);
	//printf("%f\n",Test.Y);
	
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
		sprintf(filename, "MeshInfo.mphtxt");

	if ((ip = fopen(filename,"r")) == NULL){
		printf("error opening the input data.\n");
		return;
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
		return;
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
		return;
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
void LU_De(double** A, int n, double*** L, double*** U)
{
	int i,j,k;
	for(j=0; j<n; j++)
	{
		for(i=0; i<n; i++)
		{
			if(i<=j)
			{
				(*U)[i][j]=A[i][j];
				for(k=0; k<=i-1; k++)
					(*U)[i][j]-=(*L)[i][k]*(*U)[k][j];
				if(i==j)
					(*L)[i][j]=1;
				else
					(*L)[i][j]=0;
			}
			else
			{
				(*L)[i][j]=A[i][j];
				for(k=0; k<=j-1; k++)
					(*L)[i][j]-=(*L)[i][k]*(*U)[k][j];
				(*L)[i][j]/=(*U)[j][j];
				(*U)[i][j]=0;
			}
		}
	}
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
void Gauss_Solve(double*** A, double **B,double** X, int n)
{
	int i,j,k,max_id;
	double temp,max_value;
	double *temp1;
		for(k=0;k<n;k++)
		{
		//pivot
		max_id=k;
		max_value=fabs((*A)[k][k]);
		for(i=k+1;i<n;i++)
		if(fabs((*A)[i][k])>max_value)
		{max_value=fabs((*A)[i][k]); max_id=i;}
		temp1=(*A)[k];
		(*A)[k]=(*A)[max_id];
		(*A)[max_id]=temp1;
		//Gauss Elimination
		for(i=k+1;i<n;i++)
		{
			temp=(*A)[i][k]/(*A)[k][k];
			for(j=k;j<n;j++)
				(*A)[i][j]=(*A)[i][j]-temp*(*A)[k][j];
			if (j==n)
				(*B)[i]=(*B)[i]-temp*(*B)[k];
		}
		}
		//Jordan Elimination
		for(k=n-1;k>0;k--)
		for (i=0;i<k;i++)
		{
			(*B)[i]=(*B)[i]-(*B)[k]*(*A)[i][k]/(*A)[k][k];
			(*A)[i][k]=0;
		}
	
		for(i=0;i<n;i++)
			(*B)[i]=(*B)[i]/(*A)[i][i];

		*X=*B;
		return;
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

	A[0][0]=2;A[0][1]=1;A[CubicN-1][CubicN-2]=1;A[CubicN-1][CubicN-1]=2;
	b[0]=-3*y[0]+3*y[1];
	b[CubicN-1]=-3*y[CubicN-2]+3*y[CubicN-1];

	Gauss_Solve(&A,&b,&x1,CubicN);
	return x1;
}
void Load_BH(void)
{
	//int i;
	//for (i=0;i<CubicN;i++)
	//{
	//	CubicX[i]=i+1;
	//	CubicY[i]=(i+1);
	//}
	CubicX[0]=1;CubicX[1]=2;CubicX[2]=3;
	CubicY[0]=2;CubicY[1]=3;CubicY[2]=5;
	CubicX1=Cubic_Spline(CubicY,CubicX);
}
interp_t Interpo(double X)
{
	int i,id;
	interp_t T;
	double x;
	for(i=0;i<CubicN;i++)
	{
		if(X>=CubicX[i] && X<=CubicX[i+1])
		id =i;
	}
	x=(X-CubicX[id])/(CubicX[id+1]-CubicX[id]);
	T.Y=(2*x*x*x-3*x*x+1)*CubicY[id]+(-2*x*x*x+3*x*x)*CubicY[id+1]+(x*x*x-2*x*x+x)*CubicX1[id]+(x*x*x-x*x)*CubicX1[id+1];
	T.Y1=1/(CubicX[id+1]-CubicX[id])*((6*x*x-6*x)*CubicY[id]+(-6*x*x+6*x)*CubicY[id+1]+(3*x*x-4*x+1)*CubicX1[id]+(3*x*x-2*x)*CubicX1[id+1]);
	T.Y2=1/(CubicX[id+1]-CubicX[id])*x*((12*x-6)*CubicY[id]+(-12*x+6)*CubicY[id+1]+(6*x-4)*CubicX1[id]+(6*x-2)*CubicX1[id+1]);
	return T;
}
int saveoutput(double *X,double *Y,double* S, int N){
	FILE* op;
	int i, j;

	if ((op = fopen("C:\\Users\\Peng Liu\\Desktop\\FEMout.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		return 1;
	}

	fprintf(op, "%d\n",N);
	for (i = 0; i < N; i++) 
			fprintf(op, "%lf %lf %lf\n",X[i],Y[i],S[i]);
	fclose(op);
	return 0;
}
int Newton_Raphson_Element(double *X,double *Y,double *Z,double x0, double y0, double z0, double J0,double J1,double J2,double Ke01,double Ke12,double Ke20,double Z0)
{
	double K00,K01,K02,K10,K11,K12,K20,K21,K22;
	double b0,b1,b2;
	double x=0, y=0, z=0,dx,dy,dz,Error=1;
	

	while(Error>1e-3)
	{
		b0=-(Ke01*J0*pow(x+x0,3)+Ke01*J1*(x+x0)*pow(y+y0,2)+Ke01*J2*(x+x0)*pow(z+z0,2)-(x0-x)/Z0);
		b1=-(Ke12*J0*(y+y0)*pow(x+x0,2)+Ke12*J1*pow(y+y0,3)+Ke12*J2*(y+y0)*pow(z+z0,2)-(y0-y)/Z0);
		b2=-(Ke20*J0*(z+z0)*pow(x+x0,2)+Ke20*J1*(z+z0)*pow(y+y0,2)+Ke20*J2*pow(z+z0,3)-(z0-z)/Z0);

		K00=Ke01*J0*3*pow(x+x0,2)+Ke01*J1*pow(y+y0,2)+Ke01*J2*pow(z+z0,2)+1/Z0;
		K01=Ke01*J1*2*(x+x0)*(y+y0);
		K02=Ke01*J2*2*(x+x0)*(z+z0);

		K10=Ke12*J1*2*(x+x0)*(y+y0);
		K11=Ke12*J0*pow(x+x0,2)+Ke12*J1*3*pow(y+y0,2)+Ke12*J2*pow(z+z0,2)+1/Z0;
		K12=Ke12*J2*2*(x+x0)*(z+z0);

		K20=Ke20*J0*2*(z+z0)*(x+x0);
		K21=Ke20*J1*2*(z+z0)*(y+y0);
		K22=Ke20*J0*pow(x+x0,2)+Ke20*J1*pow(y+y0,2)+Ke20*J2*3*pow(z+z0,2)+1/Z0;

		dx=Det(b0,K01,K02,b1,K11,K12,b2,K21,K22)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);
		dy=Det(K00,b0,K02,K10,b1,K12,K20,b2,K22)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);
		dz=Det(K00,K01,b0,K10,K11,b1,K20,K21,b2)/Det(K00,K01,K02,K10,K11,K12,K20,K21,K22);

		if(x!=0 && y!=0 && z!=0)
		Error=(fabs(dx/x)+fabs(dy/y)+fabs(dz/z));

		x=x+dx;
		y=y+dy;
		z=z+dz;
	}
	
	*X=x;
	*Y=y;
	*Z=z;
	return 0;
}

double Det(double x0,double x1,double x2,double y0,double y1,double y2,double z0,double z1,double z2)
{
	double T=0;
	T=x0*y1*z2+x1*y2*z0+x2*y0*z1-x0*y2*z1-x1*y0*z2-x2*y1*z0;
	return T;
}