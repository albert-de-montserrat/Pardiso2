#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <time.h>
#include <hdf5.h>
/*
*/

/* --------------------- UNCLUDE PARTS --------------------- */
#include"head3mg.c"
#include"load3mg.c"
#include"move3mg.c"
#include"mark3mg.c"
#include"gaus3mg_MG.c"
#include"heat3mg.c"
/* --------------------- UNCLUDE PARTS --------------------- */




/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
/* Counters */
int n0,n1,f0,f1,fln3;
long int pos0cur0,m0,m1,m2,m3,m4,mm1;
char command[600];
/**/
double start1,start2,start3,timesum0;
#pragma omp parallel
	{
	start1=omp_get_wtime();
	}
/**/
/**/
/* Load data from input file */
fln3=loadconf()+1;
/**/
/**/
/**/
/* Load from data file */
loader();
/**/
/**/
/**/
/* Reset iterations */
if(1==0 && maxcyc<0) 
	{
	/* Clear vx,vy,vz,P */
	for(m2=0;m2<nodenum;m2++) 
		{
		vx[m2]=vy[m2]=vz[m2]=pr[m2]=0;
		exx[m2]=eyy[m2]=ezz[m2]=exy[m2]=exz[m2]=eyz[m2]=0;
		sxx[m2]=syy[m2]=szz[m2]=sxy[m2]=sxz[m2]=syz[m2]=0;
		}
	maxcyc=-maxcyc;
	for (m3=1;m3<znumz;m3++)
	for (m1=1;m1<xnumx;m1++)
		{
		m4=xynumxy*m3+m1*ynumy+1;
		pr[m4]=pinit;
		for (m2=2;m2<ynumy;m2++)
			{
			m4++;
			pr[m4]=pr[m4-1]+(gy[m2]-gy[m2-2])/2.0*GYKOEF*(ro[m4-1]+ro[m4-ynumy-1]+ro[m4-xynumxy-1]+ro[m4-xynumxy-ynumy-1])/4.0;
			}
		}
	ronurecalc();
/*
*/
	}
/**/
/**/
/**/
/* Output File Cycle */
f1=fln3+filesjob; if (f1>fl0num) f1=fl0num;
if (printmod) printf("\n PROCESSING UP TO %d .prn FILES (from %d to %d) FOR THIS JOB \n",f1-fln3,fln3,f1-1);
for (f0=fln3;f0<filesjob;f0++)
{
/* Reload Cur Output File Name */
for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
fl1otp=fl0otp[f0];
/**/
/* Reload cyc0max_maxxystep_maxtkstep_maxtmstep */
cyc0max=fl0cyc[f0][0];
maxxystep=fl0stp[f0][0];
maxtkstep=fl0stp[f0][1];
maxtmstep=fl0stp[f0][2];
nubeg=fl0stp[f0][3];
nuend=fl0stp[f0][4];
p0koef=fl0stp[f0][5];
p1koef=fl0stp[f0][6];
p2koef=fl0stp[f0][7];
multinum=fl0cyc[f0][1];
/**/
/* Reset velocity field */
if(drexmod){
timesum0=timesum;
for(m1=0;m1<nodenum;m1++)
	{
	vx0[m1]=0;
	vy0[m1]=0;
	vz0[m1]=0;
	/*
	fd0[m1]=0;
	tk0[m1]=0;
	pr0[m1]=0;
   	*/
   	}
	}
/* General Cycle */
for (n0=0;n0<cyc0max;n0++)
	{
	if (printmod) printf("\n! FILE %s  KRUG ! %d\n",fl1out,n0+1);
	/**/
	/**/
	/**/
	/* Set initial time step */
	timestep=maxtmstep;
	if (printmod) printf("\n !!! MAX VALID TIME STEP FOR CYCLE %e YEARS !!! \n",timestep/3.15576e+7);
	/**/
	/**/
	/**/
	/* vX,vY, vZ recalc after Stokes+Contin equation */
	if(movemod)
		{
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		start2=omp_get_wtime();
		if (printmod) printf("\n G, VX, VY, VZ, P ITERATION CYCLE BEGINING (clock=%e)\n",start2);
		vpiterate(n0,f0);
		start3=omp_get_wtime();
		if (printmod) printf("\n G, VX, VY, VZ, P ITERATION CYCLE END (clock=%e difference=%e)\n",start3,start3-start2);
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		}
	/**/
	/**/
	/**/
	/* Tk recalc after Heat conservation equation */
	if(tempmod && timestep)
		{
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		start2=omp_get_wtime();
		if (printmod) printf("\n T RECALC...(clock=%e)\n",start2);
		titerate(n0);
		start3=omp_get_wtime();
		if (printmod) printf("OK!(clock=%e difference=%e)\n",start3,start3-start2);
		if (printmod) printf("\n CURRENT VALID TIME STEP %e YEARS IN CYCLE\n",timestep/3.15576e+7);
		}
	/**/
	/**/
	/**/
	/* Calculate inclusion/matrix average strain rate, deviatoric stress tensor components */
	if(effbulkviscmod)
		{
		if (printmod) printf("\n CALCULATE EFFECTIVE VISCOSITY\n");     
                effbulkvisccalc(1);
		}
	/**/
	/**/
	/**/
	/* Move marker */
	if(markmod && timestep)
		{
		start2=omp_get_wtime();
		if (printmod) printf("\n MARKERS MOVE...(clock=%e)\n",start2);
		movemark();
		start3=omp_get_wtime();
		if (printmod) printf("OK!(clock=%e difference=%e)\n",start3,start3-start2);
		}
	/**/
	/**/
	/**/
	/* ro[],nu[] Recalc */
	if(gridmod)
		{
		start2=omp_get_wtime();
		if (printmod) printf("\n RO, NU, CP etc  RECALC AFTER NEW MARKERS POSITIONS...(clock=%e)\n",start2);
		ronurecalc();
		start3=omp_get_wtime();
		if (printmod) printf("OK!(clock=%e difference=%e)\n",start3,start3-start2);
		}
	/**/
	/**/
	/**/
	/* Increse Timesum */
	timesum+=timestep;
	/* 1year=3.15576*10^7sek */
	if (printmod) printf("\n %e YEARS IN CYCLE     %e YEARS FROM START\n",timestep/3.15576e+7,timesum/3.15576e+7);
	/**/
	/**/
	/* Add velocity field */
	if(drexmod){
	for(m1=0;m1<nodenum;m1++)
		{
		vx0[m1]+=vx[m1];
		vy0[m1]+=vy[m1];
		vz0[m1]+=vz[m1];
		/*
   		tk0[m1]+=tk[m1];
		pr0[m1]+=pr[m1];
		fd0[m1]+=fd[m1];
		*/
		}
	}
	/**/
	/* Calculate inclusion/matrix average strain rate, deviatoric stress tensor components */
	if(effbulkviscmod)
		{
		if (printmod) printf("\n SAVE EFFECTIVE VISCOSITY\n");     
                effbulkvisccalc(2);
		}
	}
/**/
/**/
/* Print Results */
clockbeg=clock();
if (printmod) printf("\n SAVING RESULTS...(clock=%ld)\n",clockbeg/CLOCKS_PER_SEC);
saver(f0+1,n0-1);
clockend=clock();
if (printmod) printf("OK!(clock=%ld difference=%ld)\n",clockend/CLOCKS_PER_SEC,(clockend-clockbeg)/CLOCKS_PER_SEC);
/* Print Results */
/* Save timestep, timesum and velocity interpolated on nodes */
if(f0+1>10 && drexmod)
{
char c[3];
int mgi=0;
long int wn1[100],mcmax1;
double vxyz1[3],eps1[100];
hid_t       file_id, group_id, attr_id, dataset_id, dataspace_id, dcpl,memspace,plist_id;
herr_t      status;
hsize_t     dims,dims2,ndims[2],count[2],offset[2],chunk_dims[2],stride[2],block[2];
int n1;
/**/
char sa[]="vtp";
if(f0+1<10) strcat(sa,"00");
if(f0+1>=10 && f0+1<100) strcat(sa,"0");
sprintf(c,"%d",f0+1);
strcat(sa,c);
strcat(sa,".h5");
file_id = H5Fcreate(sa,H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
printf("\n %s\n",sa);
dims=2;
double datadb[dims];
datadb[0]=timesum-timesum0;
datadb[1]=timesum;
dataspace_id=H5Screate_simple(1,&dims,NULL);
attr_id = H5Acreate(file_id,"General",H5T_NATIVE_DOUBLE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
H5Awrite(attr_id,H5T_NATIVE_DOUBLE,&datadb);
H5Aclose(attr_id);
H5Sclose(dataspace_id);
//
//
/* Eulerian nodes infos */
group_id = H5Gcreate(file_id, "/Nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
/* Number of nodes along X,Y,Z axes */
dims=3;
double datali[dims];
datali[0]=(double)(xnumx);
datali[1]=(double)(ynumy);
datali[2]=(double)(znumz);
dataspace_id=H5Screate_simple(1,&dims,NULL);
attr_id = H5Acreate(group_id,"nxyz",H5T_NATIVE_DOUBLE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
H5Awrite(attr_id,H5T_NATIVE_DOUBLE,&datali);
H5Aclose(attr_id);
H5Sclose(dataspace_id);
/**/
/* Gridlines positions */
dims=xnumx;
loadsave_dataset_double(group_id,dims,gx, "gx",0);
dims=ynumy;
loadsave_dataset_double(group_id,dims,gy, "gy",0);
dims=znumz;
loadsave_dataset_double(group_id,dims,gz, "gz",0);
/**/
/* Nodes parameters */
/* Average velocity field */
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
        {
        /* Pos in vx[], vy[], vz[], pr[], etc */
        mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
        //vxyzcalc1(gx[m1],gy[m2],gz[m3],vxyz1,wn1);
        //prcalc1(gx[m1],gy[m2],gz[m3],eps1,wn1);
        //pr0[mcmax1]=eps1[10];
	vx0[mcmax1]/=cyc0max;
	vy0[mcmax1]/=cyc0max;
	vz0[mcmax1]/=cyc0max;
	//fd0[mcmax1]/=cyc0max;
	}
dims=nodenum;
//loadsave_dataset_double(group_id,dims,pr, "pr0",0);
loadsave_dataset_float(group_id,dims,vx0, "vx0",0);
loadsave_dataset_float(group_id,dims,vy0, "vy0",0);
loadsave_dataset_float(group_id,dims,vz0, "vz0",0);
//loadsave_dataset_float(group_id,dims,tk, "tk",0);
//loadsave_dataset_float(group_id,dims,fd, "fd0",0);
//
H5Gclose(group_id);
H5Fclose(file_id);
/* Move vtp file to output directory */	 	
char d[3];
char sc[]="vtp";
if(f0+1<10) strcat(sc,"00");
if(f0+1>=10 && f0+1<100) strcat(sc,"0");
sprintf(d,"%d",f0+1);
strcat(sc,d);
strcat(sc,".h5");
sprintf(command,"mv %s %s",sc,outdir);
system(command);
}
/**/
/**/
if(omp_get_wtime()-start1>3600*23.0) break;
/* End Output file Names Cycle */
}
/**/
sprintf(command,"mv effbulkvisc.txt %s",outdir);
system(command);
sprintf(command,"mv strainbudget.txt %s",outdir);
system(command);
/* Free memory */
dynmemall(3);
/* End Program */
return 0;
}
/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */



