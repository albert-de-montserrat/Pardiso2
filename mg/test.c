#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <hdf5.h>
/*
*/
#include"head3mg.c"

/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
/* Counters */
int cycle,nt,n0,n1,f0,f1,fln3,tid;
long int pos0cur0,i,m0,m1,m2,m3,m4,mm1,nodenum=50000000;
/**/
double start1,start2,start3,timesum0;
for(cycle=0;cycle<2;cycle++)
{
#pragma omp parallel
	{
	start1=omp_get_wtime();
	}
#pragma omp parallel
{nt=omp_get_num_threads();}
printf("Number of threads %d\n",nt);
if(cycle==0)
{
nu=malloc(sizeof(float *) * nodenum);
nxx=malloc(sizeof(float *) * nodenum);
nxy=malloc(sizeof(float *) * nodenum);
nxz=malloc(sizeof(float *) * nodenum);
nyz=malloc(sizeof(float *) * nodenum);
ro=malloc(sizeof(float *) * nodenum);
sol0=malloc(sizeof(double *) * nodenum);
sol1=malloc(sizeof(double *) * nodenum);
Nu=malloc(sizeof(float *) * nt);
Nxx=malloc(sizeof(float *) * nt);
Nxy=malloc(sizeof(float *) * nt);
Nxz=malloc(sizeof(float *) * nt);
Nyz=malloc(sizeof(float *) * nt);
Ro=malloc(sizeof(float *) * nt);
Sol0=malloc(sizeof(double *) * nt);
Sol1=malloc(sizeof(double *) * nt);
for(i=0;i<nt;i++)
	{
	Nu[i]=malloc(sizeof(float) * nodenum);
	Nxx[i]=malloc(sizeof(float) * nodenum);
	Nxy[i]=malloc(sizeof(float) * nodenum);
	Nxz[i]=malloc(sizeof(float) * nodenum);
	Nyz[i]=malloc(sizeof(float) * nodenum);
	Ro[i]=malloc(sizeof(float) * nodenum);
	Sol0[i]=malloc(sizeof(double) * nodenum);
	Sol1[i]=malloc(sizeof(double) * nodenum);
	}
}
start2=omp_get_wtime();
//#pragma omp parallel for shared(Sol0,Sol1,Ro,Nu,Nxx,Nxy,Nxz,Nyz) private(m1,tid) firstprivate(nodenum) schedule(static)
/* Obtain cur thread number */
for (tid=0;tid<nt;tid++)
	{
#pragma omp parallel for
for (m1=0;m1<nodenum;m1++)
	{
	Sol0[tid][m1]=0;
	Sol1[tid][m1]=0;
	Ro[tid][m1]=0;
	Nu[tid][m1]=0;
	Nxx[tid][m1]=0;
	Nxy[tid][m1]=0;
	Nxz[tid][m1]=0;
	Nyz[tid][m1]=0;
	}
	}
#pragma omp parallel for
for (m1=0;m1<nodenum;m1++)
	{
	sol0[m1]=0;
	sol1[m1]=0;
	ro[m1]=0;
	nu[m1]=0;
	nxx[m1]=0;
	nxy[m1]=0;
	nxz[m1]=0;
	nyz[m1]=0;
	}
//double rho;
//#pragma omp parallel for shared(Sol0,Sol1,Ro,Nu,Nxx,Nxy,Nxz,Nyz) private(tid,m1) firstprivate(nodenum) schedule(static) 
//#pragma omp parallel for private(m1) schedule(static) //collapse(2)
/* Obtain cur thread number */
for (tid=0;tid<nt;tid++)
	{
#pragma omp parallel for
for (m1=0;m1<nodenum;m1++)
	{
	//rho=Ro[tid][m1];
	Sol0[tid][m1]+=0.1;
	Sol1[tid][m1]+=0.2;
	Ro[tid][m1]+=0.3;
	Nu[tid][m1]+=0.4;
	Nxx[tid][m1]+=0.5;
	Nxy[tid][m1]+=0.6;
	Nxz[tid][m1]+=0.7;
	Nyz[tid][m1]+=0.8;
//if(Ro[tid][m1]-rho>0.300001) {printf("\n %d %ld %e %e %e\n",tid,m1,rho,Ro[tid][m1],Ro[tid][m1]-rho); exit(0);}
	}
	}
for (tid=0;tid<nt;tid++)
	{
#pragma omp parallel for
for (m1=0;m1<nodenum;m1++ )
	{
	/*
        Sol0[0][1]+=Sol0[tid][1];
	Sol1[0][1]+=Sol1[tid][1];
	Ro[0][1]+=Ro[tid][1];
	Nu[0][1]+=Nu[tid][1];
	Nxx[0][1]+=Nxx[tid][1];
	Nxy[0][1]+=Nxy[tid][1];
	Nxz[0][1]+=Nxz[tid][1];
	Nyz[0][1]+=Nyz[tid][1];
	*/
        sol0[m1]+=Sol0[tid][m1];
	sol1[m1]+=Sol1[tid][m1];
	ro[m1]+=Ro[tid][m1];
	nu[m1]+=Nu[tid][m1];
	nxx[m1]+=Nxx[tid][m1];
	nxy[m1]+=Nxy[tid][m1];
	nxz[m1]+=Nxz[tid][m1];
	nyz[m1]+=Nyz[tid][m1];
	}
}
printf("\n %d %e %e %e %e %e %e %e %e\n",tid,sol0[1],sol1[1],ro[1],nu[1],nxx[1],nxy[1],nxz[1],nyz[1]);
//printf("\n %d %e %e %e %e %e %e %e %e\n",tid,Sol0[0][1],Sol1[0][1],Ro[0][1],Nu[0][1],Nxx[0][1],Nxy[0][1],Nxz[0][1],Nyz[0][1]);
start3=omp_get_wtime();
printf("\n Cycle = %d (clock=%e difference=%e)\n",cycle,start3,start3-start2);
}
}
