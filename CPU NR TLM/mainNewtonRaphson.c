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
int saveoutput(double *X,double *Y,double* S, int N);
double* Cubic_Spline(double y[], double x[]);
int saveoutAxb(double **A,double *b, int N);
interp_t Interpo(double x);
//define BH info.
const int CubicN=3;
double CubicX[3], CubicY[3];
double *CubicX1;

void main(){
	const int MaxNodes=12000;
	const int MaxElem=MaxNodes*2;
	
	int i,j,k=0;
	int NumNodes, NumElem;
	int *VertI=(int*)malloc(MaxElem*sizeof(int)),*VertJ=(int*)malloc(MaxElem*sizeof(int)),*VertK=(int*)malloc(MaxElem*sizeof(int));
	double *X=(double*)malloc(MaxNodes*sizeof(double)), *Y=(double*)malloc(MaxNodes*sizeof(double));
	double *Area,*a1,*b1,*c1,*a2,*b2,*c2,*a3,*b3,*c3;
	double x1,x2,x3,y1,y2,y3,Ke00,Ke01,Ke02,Ke10,Ke11,Ke12,Ke20,Ke21,Ke22,J0,J1,J2,J3,J4;
	double **Kc,**Kt,*Ks,*Js,*Vr,*A,*DeltA,**Jaco,*VB2,*JacoB,**AdmitL,**AdmitU;
	double T_start,T_end,interp,Newton_Relative_Error=1,Error_Temp;
	interp_t Test;
	int *BoundaryNodeId;
	double *BoundaryNodeValue;
//Load Mesh Information from file
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

	Kc=(double**)malloc(NumNodes*sizeof(double*));//Space term Matrix
	for (i=0;i<NumNodes;i++)
		Kc[i]=(double*)malloc(NumNodes*sizeof(double));

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
		Js[i]=-1;//current	
	for (i=0;i<NumElem;i++)
		
	{
		Vr[i]=1;//material initialized
		VB2[i]=0.1;//Vr=0.1B^2
	}
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

//Get Stiffness Matrix and Assemble
	GET_TIME(T_start);
	for(i=0;i<NumElem;i++)
	{
		x1=X[VertI[i]]; x2=X[VertJ[i]];x3=X[VertK[i]];
		y1=Y[VertI[i]]; y2=Y[VertJ[i]];y3=Y[VertK[i]];
		Area[i]=0.5*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2));
		a1[i]=x2*y3-x3*y2; b1[i]=y2-y3; c1[i]=x3-x2;
		a2[i]=x3*y1-x1*y3; b2[i]=y3-y1; c2[i]=x1-x3;
		a3[i]=x1*y2-x2*y1; b3[i]=y1-y2; c3[i]=x2-x1;

		Ke00=1.0/4*Vr[i]/Area[i]*(b1[i]*b1[i]+c1[i]*c1[i]);
		Ke01=1.0/4*Vr[i]/Area[i]*(b1[i]*b2[i]+c1[i]*c2[i]);
		Ke02=1.0/4*Vr[i]/Area[i]*(b1[i]*b3[i]+c1[i]*c3[i]);

		Ke10=1.0/4*Vr[i]/Area[i]*(b1[i]*b2[i]+c1[i]*c2[i]);
		Ke11=1.0/4*Vr[i]/Area[i]*(b2[i]*b2[i]+c2[i]*c2[i]);
		Ke12=1.0/4*Vr[i]/Area[i]*(b3[i]*b2[i]+c3[i]*c2[i]);

		Ke20=1.0/4*Vr[i]/Area[i]*(b1[i]*b3[i]+c1[i]*c3[i]);
		Ke21=1.0/4*Vr[i]/Area[i]*(b2[i]*b3[i]+c2[i]*c3[i]);
		Ke22=1.0/4*Vr[i]/Area[i]*(b3[i]*b3[i]+c3[i]*c3[i]);

		Kc[VertI[i]][VertI[i]]=Kc[VertI[i]][VertI[i]]+Ke00;
		Kc[VertI[i]][VertJ[i]]=Kc[VertI[i]][VertJ[i]]+Ke01;
		Kc[VertI[i]][VertK[i]]=Kc[VertI[i]][VertK[i]]+Ke02;

		Kc[VertJ[i]][VertI[i]]=Kc[VertJ[i]][VertI[i]]+Ke10;
		Kc[VertJ[i]][VertJ[i]]=Kc[VertJ[i]][VertJ[i]]+Ke11;
		Kc[VertJ[i]][VertK[i]]=Kc[VertJ[i]][VertK[i]]+Ke12;

		Kc[VertK[i]][VertI[i]]=Kc[VertK[i]][VertI[i]]+Ke20;
		Kc[VertK[i]][VertJ[i]]=Kc[VertK[i]][VertJ[i]]+Ke21;
		Kc[VertK[i]][VertK[i]]=Kc[VertK[i]][VertK[i]]+Ke22;

		Ke00=SigmaE*Area[i]/6; Ke01=SigmaE*Area[i]/12;Ke02=SigmaE*Area[i]/12;
		Ke10=SigmaE*Area[i]/12;Ke11=SigmaE*Area[i]/6; Ke12=SigmaE*Area[i]/12;
		Ke20=SigmaE*Area[i]/12;Ke21=SigmaE*Area[i]/12;Ke22=SigmaE*Area[i]/6;

		Kt[VertI[i]][VertI[i]]=Kt[VertI[i]][VertI[i]]+Ke00;
		Kt[VertI[i]][VertJ[i]]=Kt[VertI[i]][VertJ[i]]+Ke01;
		Kt[VertI[i]][VertK[i]]=Kt[VertI[i]][VertK[i]]+Ke02;

		Kt[VertJ[i]][VertI[i]]=Kt[VertJ[i]][VertI[i]]+Ke10;
		Kt[VertJ[i]][VertJ[i]]=Kt[VertJ[i]][VertJ[i]]+Ke11;
		Kt[VertJ[i]][VertK[i]]=Kt[VertJ[i]][VertK[i]]+Ke12;

		Kt[VertK[i]][VertI[i]]=Kt[VertK[i]][VertI[i]]+Ke20;
		Kt[VertK[i]][VertJ[i]]=Kt[VertK[i]][VertJ[i]]+Ke21;
		Kt[VertK[i]][VertK[i]]=Kt[VertK[i]][VertK[i]]+Ke22;

		Ke00=Js[i]*Area[i]/3.0; Ke01=Js[i]*Area[i]/3.0; Ke02=Js[i]*Area[i]/3.0;
		Ks[VertI[i]]=Ks[VertI[i]]+Ke00;
		Ks[VertJ[i]]=Ks[VertJ[i]]+Ke01;
		Ks[VertK[i]]=Ks[VertK[i]]+Ke02;
	}
	GET_TIME(T_end);
	printf("\nElapsed time for assembling is %.8f\n\n",(T_end-T_start));
	//Newton-Raphson Iteration 
	GET_TIME(T_start);
	while(Newton_Relative_Error>1e-6){
//Apply A Boundary 
		Error_Temp=0;
		for(i=0;i<NumNodes;i++)
			if(BoundaryNodeId[i]==1)
				A[i]=BoundaryNodeValue[i];

		for(i=0;i<NumNodes;i++)
		{
			JacoB[i]=0;
			for(j=0;j<NumNodes;j++)
				Jaco[i][j]=0;
		}

//Get Jacobian Matrix and RHS vector with given A and Vr
	for(i=0;i<NumElem;i++)
	{	
		J0=(b1[i]*b1[i]+c1[i]*c1[i])*A[VertI[i]]+(b1[i]*b2[i]+c1[i]*c2[i])*A[VertJ[i]]+(b1[i]*b3[i]+c1[i]*c3[i])*A[VertK[i]];
		J1=(b1[i]*b2[i]+c1[i]*c2[i])*A[VertI[i]]+(b2[i]*b2[i]+c2[i]*c2[i])*A[VertJ[i]]+(b3[i]*b2[i]+c3[i]*c2[i])*A[VertK[i]];
		J2=(b1[i]*b3[i]+c1[i]*c3[i])*A[VertI[i]]+(b2[i]*b3[i]+c2[i]*c3[i])*A[VertJ[i]]+(b3[i]*b3[i]+c3[i]*c3[i])*A[VertK[i]];
		J3=b1[i]*A[VertI[i]]+b2[i]*A[VertJ[i]]+b3[i]*A[VertK[i]];
		J4=c1[i]*A[VertI[i]]+c2[i]*A[VertJ[i]]+c3[i]*A[VertK[i]];

		

		Ke00=1.0/4*Vr[i]/Area[i]*(b1[i]*b1[i]+c1[i]*c1[i])+VB2[i]/4.0/Area[i]*J0*(2*b1[i]*J3+2*c1[i]*J4)/4/Area[i]/Area[i];
		Ke01=1.0/4*Vr[i]/Area[i]*(b1[i]*b2[i]+c1[i]*c2[i])+VB2[i]/4.0/Area[i]*J0*(2*b2[i]*J3+2*c2[i]*J4)/4/Area[i]/Area[i];
		Ke02=1.0/4*Vr[i]/Area[i]*(b1[i]*b3[i]+c1[i]*c3[i])+VB2[i]/4.0/Area[i]*J0*(2*b3[i]*J3+2*c3[i]*J4)/4/Area[i]/Area[i];

		Ke10=1.0/4*Vr[i]/Area[i]*(b1[i]*b2[i]+c1[i]*c2[i])+VB2[i]/4.0/Area[i]*J1*(2*b1[i]*J3+2*c1[i]*J4)/4/Area[i]/Area[i];
		Ke11=1.0/4*Vr[i]/Area[i]*(b2[i]*b2[i]+c2[i]*c2[i])+VB2[i]/4.0/Area[i]*J1*(2*b2[i]*J3+2*c2[i]*J4)/4/Area[i]/Area[i];
		Ke12=1.0/4*Vr[i]/Area[i]*(b3[i]*b2[i]+c3[i]*c2[i])+VB2[i]/4.0/Area[i]*J1*(2*b3[i]*J3+2*c3[i]*J4)/4/Area[i]/Area[i];

		Ke20=1.0/4*Vr[i]/Area[i]*(b1[i]*b3[i]+c1[i]*c3[i])+VB2[i]/4.0/Area[i]*J2*(2*b1[i]*J3+2*c1[i]*J4)/4/Area[i]/Area[i];
		Ke21=1.0/4*Vr[i]/Area[i]*(b2[i]*b3[i]+c2[i]*c3[i])+VB2[i]/4.0/Area[i]*J2*(2*b2[i]*J3+2*c2[i]*J4)/4/Area[i]/Area[i];
		Ke22=1.0/4*Vr[i]/Area[i]*(b3[i]*b3[i]+c3[i]*c3[i])+VB2[i]/4.0/Area[i]*J2*(2*b3[i]*J3+2*c3[i]*J4)/4/Area[i]/Area[i];
	
		Jaco[VertI[i]][VertI[i]]=Jaco[VertI[i]][VertI[i]]+Ke00;
		Jaco[VertI[i]][VertJ[i]]=Jaco[VertI[i]][VertJ[i]]+Ke01;
		Jaco[VertI[i]][VertK[i]]=Jaco[VertI[i]][VertK[i]]+Ke02;

		Jaco[VertJ[i]][VertI[i]]=Jaco[VertJ[i]][VertI[i]]+Ke10;
		Jaco[VertJ[i]][VertJ[i]]=Jaco[VertJ[i]][VertJ[i]]+Ke11;
		Jaco[VertJ[i]][VertK[i]]=Jaco[VertJ[i]][VertK[i]]+Ke12;

		Jaco[VertK[i]][VertI[i]]=Jaco[VertK[i]][VertI[i]]+Ke20;
		Jaco[VertK[i]][VertJ[i]]=Jaco[VertK[i]][VertJ[i]]+Ke21;
		Jaco[VertK[i]][VertK[i]]=Jaco[VertK[i]][VertK[i]]+Ke22;

		JacoB[VertI[i]]=JacoB[VertI[i]]-1.0/4*Vr[i]/Area[i]*J0+Js[i]*Area[i]/3.0;
		JacoB[VertJ[i]]=JacoB[VertJ[i]]-1.0/4*Vr[i]/Area[i]*J1+Js[i]*Area[i]/3.0;
		JacoB[VertK[i]]=JacoB[VertK[i]]-1.0/4*Vr[i]/Area[i]*J2+Js[i]*Area[i]/3.0;
		//Update Vr
		Vr[i]=0.1*(J3*J3+J4*J4)/4/Area[i]/Area[i];
	}
//Apply Boundary on DeltA
		for(i=0;i<NumNodes;i++)
		{
			if(BoundaryNodeId[i]==1)
			{
				JacoB[i]=0;
				for(j=0;j<NumNodes;j++)
					if(j==i)
						Jaco[i][j]=1;
					else
					{
						Jaco[i][j]=0;
						Jaco[j][i]=0;
					}
				
			}
		}

	
		//saveoutAxb(Jaco,JacoB,NumNodes);
		Gauss_Solve(&Jaco,&JacoB,&DeltA,NumNodes);
	/*	for(i=0;i<NumNodes;i++)
			printf("%f ",DeltA[i]);
			printf("\n");*/
//Update A
		for(i=0;i<NumNodes;i++)
		{
			A[i]=A[i]+DeltA[i];
			if(BoundaryNodeId[i]!=1)
			Error_Temp=Error_Temp+fabs(1.0*DeltA[i]/A[i]);
		}

		Newton_Relative_Error=Error_Temp/NumNodes;
		//printf("Value on sample point for %d th iteration is %f\n",k++,A[1]);

		
	}
	GET_TIME(T_end);
	
	//Print before boundary conditions
	//Apply Boundary Conditions
	//Apply boundary
	//Change the stiffness matrix 
	//for(i=0;i<NumNodes;i++)
	//{
	//	if(BoundaryNodeId[i]==1)
	//	{
	//		Ks[i]=BoundaryNodeValue[i];
	//		for(j=0;j<NumNodes;j++)
	//		{
	//			if(j==i)
	//				Kc[i][j]=1;
	//			else
	//			{
	//				Ks[j]=Ks[j]-Ks[i]*Kc[j][i];
	//				Kc[i][j]=0;
	//				Kc[j][i]=0;
	//				
	//			}
	//		}
	//	}
	//}
	//Solution phase
	/*GET_TIME(T_start);
	LU_De(Kc,NumNodes,&AdmitL,&AdmitU);
	GET_TIME(T_end);
	printf("\nElapsed time for LU decomposition is %.8f\n\n",(T_end-T_start));

	GET_TIME(T_start);
	LU_Solve(AdmitL,AdmitU,Ks,&A,NumNodes);
	GET_TIME(T_end);
	printf("\nElapsed time for back substitution is %.8f\n\n",(T_end-T_start));*/
	//Gauss_Solve(&Kc,&Ks,&A,NumNodes);
//Print Stiffness Matrix and RHS Vector 
/*	printf("The Stiffness Matrix is: \n");
	for(i=0;i<NumNodes;i++)
	{for(j=0;j<NumNodes;j++)
	printf(" %lf ",Kc[i][j]);
	printf("\n");}

	for(i=0;i<NumNodes;i++)
	{for(j=0;j<NumNodes;j++)
	printf(" %lf ",Kt[i][j]);
	printf("\n");}
	printf("\nThe RHS vector is: \n");
	for(j=0;j<NumNodes;j++)
	{printf(" %lf ",Ks[j]);
	printf("\n");}*/
	/*printf("\nL Matrix is: \n");
	for(i=0;i<NumNodes;i++)
	{for(j=0;j<NumNodes;j++)
	printf(" %lf ",AdmitL[i][j]);
	printf("\n");}
	printf("\nU Matrix is: \n");
	for(i=0;i<NumNodes;i++)
	{for(j=0;j<NumNodes;j++)
	printf(" %lf ",AdmitU[i][j]);
	printf("\n");}*/
	
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
		return 1;
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
		return 1;
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
		return 1;
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
		return 0;
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

int saveoutAxb(double **A,double *b, int N){
	FILE* op;
	int i, j;

	if ((op = fopen("C:\\Users\\Peng Liu\\Desktop\\Axbout.txt","w")) == NULL){
		printf("Error opening the output file.\n");
		return 1;
	}

	fprintf(op, "%d\n",N);
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++) 
		fprintf(op, "%lf ",A[i][j]);
		fprintf(op, "\n");
	}
	fprintf(op, "\n");
	fprintf(op, "\n");
	for (j = 0; j < N; j++) 
		fprintf(op, "%lf\n",b[j]);
	fclose(op);
	return 0;
}