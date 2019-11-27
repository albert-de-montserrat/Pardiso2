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
#include"head3viz.c"
#include"load3viz.c"
/* --------------------- UNCLUDE PARTS --------------------- */




/* Solve differencial equations by MARKER-IN-CELL + FINITE-DIFFERENCES method */
int main()
{
/* Counters */
int n0,n1,f0,f1;
long int pos0cur0,m0,m1,m2,m3,m4;
char c[3],sb[50],sc[50],sd[50];
char fl0in[MAXFLN][60];
/**/
/**/
/**/
/* Open File amir.t3c */
fl = fopen("viz3D.t3c","rt");
/**/
/* Input directory */
ffscanf();
for (n1=0;n1<500;n1++) indir0[n1]=sa[n1];
sprintf(indir,"/%s",indir0);
/* Output directory */
ffscanf();
for (n1=0;n1<500;n1++) outdir0[n1]=sa[n1];
sprintf(outdir,"/%s",outdir0);
// Make output directory
char command[600];
sprintf(command,"mkdir -p %s",outdir);
system(command);
/* Read Input file name */
ffscanf();
for (m1=0;m1<50;m1++) sd[m1]=sc[m1]=sb[m1]=sa[m1];
/* Read nodes/marker infos to be saved */
ffscanf(); fl0stp2[0]=atoi(sa); /* Eulerian Mesh */
ffscanf(); fl0stp2[1]=atoi(sa); /* Density */
ffscanf(); fl0stp2[2]=atoi(sa); /* Viscosity */
ffscanf(); fl0stp2[3]=atoi(sa); /* Temperature */
ffscanf(); fl0stp2[4]=atoi(sa); /* Velocity */
ffscanf(); fl0stp2[5]=atoi(sa); /* Strain rate invariant */
ffscanf(); fl0stp2[6]=atoi(sa); /* Stress invariant */
ffscanf(); fl0stp2[7]=atoi(sa); /* Pressure */
ffscanf(); fl0stp2[8]=atof(sa); /* Heat dissipated below the specified T(°K) */
ffscanf(); fl0stp2[9]=atoi(sa); /* Diffusion/Dislocaiton creep */
ffscanf(); fl0stp2[10]=atoi(sa); /* Cumulative strain */
ffscanf(); fl0stp2[11]=atoi(sa); /* Marker type */
ffscanf(); fl0stp2[12]=atoi(sa); /* Vorticity */
/**/
/* Read Input File number */
ffscanf();
fl0num=0;
while(sa[0]!='~')
        {
        /* Check file Counter */
        if(fl0num>=MAXFLN) {printf("Space out in fl0out[]"); exit(0);}
        /* Save Input file number */
	for (m1=0;m1<50;m1++) sc[m1]=sb[m1]=sd[m1];
	num[fl0num]=atoi(sa);
	if(num[fl0num]<10) strcat(sb,"00");
	if(num[fl0num]>=10 && num[fl0num]<100) strcat(sb,"0");
	sprintf(c,"%d",num[fl0num]);
	strcat(sb,c);
	strcat(sb,".h5");
	for (n1=0;n1<50;n1++) fl0in[fl0num][n1]=sb[n1];
        /**/
        /* Save Output file name */
	if(num[fl0num]<10) strcat(sc,"00");
	if(num[fl0num]>=10 && num[fl0num]<100) strcat(sc,"0");
	strcat(sc,c);
	strcat(sc,"nodes.h5");
        for (n1=0;n1<50;n1++) fl0out[fl0num][n1]=sc[n1];
	/**/
	fl0num++;
	/**/
        ffscanf();
	/**/
	}
/**/  
/* Output File Cycle */
for (f0=0;f0<fl0num;f0++)
{       
/* Reload Cur Input File Name, Type */
for (n1=0;n1<50;n1++) fl1in[n1]=fl0in[f0][n1];
fl1itp=2;
/* Add directory to file name */
sprintf(fl1in1,"%s%s",outdir,fl1in);
/**/ 
/* Reload Cur Output File Name, Type */
for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
fl1otp=fl0otp[f0];
/* Load fr m data file */
/*
loader();
*/
makeXMFviz(f0);
printf("Result printed from %s to %s \n",fl1in,fl1out);
}
/**/
}


/* Make XMF file to visualize .h5 fields in Paraview -------------*/
void makeXMFviz(int f0)
{
long int m1,m2,m3,m4;
/**/
hid_t       file_id, file_id_0, file_id_1, group_id, attr_id, dataset_id, dataspace_id, dcpl, memspace;
herr_t      status;
hsize_t     dims,ndims[2],count[2],offset[2];
double a[100];
char command[600];
/**/
/* Open an existing file */
file_id_0 = H5Fopen(fl1in1, H5F_ACC_RDONLY, H5P_DEFAULT);
/* Create output hdf5 file */
file_id_1 = H5Fcreate(fl1out,H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
/* General */
attr_id=H5Aopen_by_name(file_id_0,".","General",H5P_DEFAULT,H5P_DEFAULT);
H5Aread (attr_id, H5T_NATIVE_DOUBLE, &a);
H5Aclose(attr_id);
/**/
xsize=a[0];
ysize=a[1];
zsize=a[2];
timesum=a[10]; timesum*=3.15576e+7;
/* periodic */
int Periodic[3];
attr_id=H5Aopen_by_name(file_id_0,".","XYZ-periodic",H5P_DEFAULT,H5P_DEFAULT);
H5Aread (attr_id, H5T_NATIVE_INT, &Periodic);
H5Aclose(attr_id);
xperiodic=Periodic[0];
yperiodic=Periodic[1];
zperiodic=Periodic[2];
/* Form xmf file name */
char fname[50];
char c[3];
if(num[f0]<10)
{
strcpy(fname,"XDMF.00");
}
if(num[f0]>=10 && num[f0]<100)
{
strcpy(fname,"XDMF.0");
}
if(num[f0]>=100 && num[f0]<1000)
{
strcpy(fname,"XDMF.");
}
sprintf(c,"%d",num[f0]);
strcat(fname,c); 
strcat(fname,".xmf");
printf("\nGenerate file %s\n",fname);
/* Open xmf file */  
fl = fopen(fname,"wt");
/* Write xmf file */
fprintf(fl,"<?xml version=\"1.0\" ?> \n");
fprintf(fl,"<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n\n");
fprintf(fl," <Domain>\n\n");
/* Eulerian grid */
if(fl0stp2[0] || fl0stp2[1] || fl0stp2[2] || fl0stp2[3] || fl0stp2[4] || fl0stp2[5] || fl0stp2[6] || fl0stp2[7] || fl0stp2[8] || fl0stp2[9] || fl0stp2[12])
{
/* Read Nodes */
group_id = H5Gopen(file_id_0, "/Nodes", H5P_DEFAULT);
dims=3;
long int b[dims];
attr_id=H5Aopen_by_name(group_id,".","/Nodes/Number",H5P_DEFAULT,H5P_DEFAULT);
H5Aread (attr_id, H5T_NATIVE_LONG, &b);
H5Aclose(attr_id);
xnumx=b[0];
ynumy=b[1];
znumz=b[2];
dynmemall(0);
/* Gridlines positions */
dims=xnumx;
loadsave_dataset_double(group_id,dims,gx, "/Nodes/gx",1);
dims=ynumy;
loadsave_dataset_double(group_id,dims,gy, "/Nodes/gy",1);
dims=znumz;
loadsave_dataset_double(group_id,dims,gz, "/Nodes/gz",1);
gridcheck();
/* Eulerian nodes informations */
fprintf(fl,"    <Grid Name=\"Eulerian Grid\">\n\n");
fprintf(fl,"       <Time Value=\"%f\" />\n\n",timesum/0.9);
fprintf(fl,"         <Topology Type=\"Hexahedron\" NumberOfElements=\"%ld\"> \n",cellnum);
fprintf(fl,"            <DataItem Format=\"HDF\" DataType=\"Int\"  Dimensions=\"%ld 8\">Mesh.h5:/connectivity</DataItem> \n",cellnum);
fprintf(fl,"         </Topology> \n\n");
fprintf(fl,"         <Geometry Type=\"XYZ\"> \n");
fprintf(fl,"            <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">Mesh.h5:/vertices</DataItem> \n",nodenum);
fprintf(fl,"         </Geometry>\n\n");
/* Make mesh strcture file when run in3mg */
if(fl0stp2[0])
	{
	/* General */
	file_id = H5Fcreate("Mesh.h5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	/* Cell vertexes */
	long int xy,y;
	y=ynumy;
	xy=xnumx*ynumy;
	for (m3=0;m3<znumz;m3++)
	for (m1=0;m1<xnumx;m1++)
	for (m2=0;m2<ynumy;m2++)
		{
		m4=m2+y*m1+xy*m3;
		vertices[m4][0]=(float)(gx[m1]);
		vertices[m4][1]=(float)(ysize-gy[m2]);
		/*
		vertices[m4][1]=(float)(gy[m2]);
		*/
		vertices[m4][2]=(float)(gz[m3]);
		}
	ndims[1]=3;
	ndims[0]=nodenum;
	dataspace_id=H5Screate_simple(2,ndims,NULL);
	dataset_id = H5Dcreate(file_id, "/vertices", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vertices);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
	/* Connectivity */
	long int xyc,yc;
	yc=ynumy1;
	xyc=xnumx1*ynumy1;
	for (m3=0;m3<znumz1;m3++)
	for (m1=0;m1<xnumx1;m1++)
	for (m2=0;m2<ynumy1;m2++)
		{
		m4=m2+yc*m1+xyc*m3;
		connectivity[m4][0]=m2+m1*y+m3*xy;
		connectivity[m4][1]=connectivity[m4][0]+1;
		connectivity[m4][2]=connectivity[m4][0]+y+1;
		connectivity[m4][3]=connectivity[m4][0]+y;
		connectivity[m4][4]=connectivity[m4][0]+xy;
		connectivity[m4][5]=connectivity[m4][0]+xy+1;
		connectivity[m4][6]=connectivity[m4][0]+xy+y+1;
		connectivity[m4][7]=connectivity[m4][0]+xy+y;
		}
	ndims[1]=8;
	ndims[0]=cellnum;
	dataspace_id=H5Screate_simple(2,ndims,NULL);
	dataset_id = H5Dcreate(file_id, "/connectivity", H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, connectivity);
	H5Dclose(dataset_id);
	H5Sclose(dataspace_id);
        /**/
        /* Close hdf5 file */
        H5Fclose(file_id);
        /**/
	sprintf(command,"mv Mesh.h5 %s",outdir);
	system(command);
	}
/**/
dims=nodenum;
if(fl0stp2[1]) 
	{
	ro=malloc(sizeof(float) * nodenum);
	loadsave_dataset_float(group_id,dims,ro, "/Nodes/ro",1);
	loadsave_dataset_float(file_id_1,dims,ro, "/ro",0);
	add_eulerian_to_XDMF("Density",fl1out,"/ro",0,nodenum,fl);
	}
/**/
if(fl0stp2[2])
	{
	nu=malloc(sizeof(float) * nodenum);
	loadsave_dataset_float(group_id,dims,nu, "/Nodes/nu",1);
	loadsave_dataset_float(file_id_1,dims,nu, "/nu",0);
	add_eulerian_to_XDMF("Viscosity",fl1out,"/nu",0,nodenum,fl);
	}
/**/
if(fl0stp2[3])
	{
	tk=malloc(sizeof(float) * nodenum);
	loadsave_dataset_float(group_id,dims,tk, "/Nodes/tk",1);
	loadsave_dataset_float(file_id_1,dims,tk, "/tk",0);
	add_eulerian_to_XDMF("Temperature",fl1out,"/tk",0,nodenum,fl);
	}
if(fl0stp2[9]) 
	{
	fd=malloc(sizeof(float) * nodenum);
	loadsave_dataset_float(group_id,dims,fd, "/Nodes/fd",1);
	loadsave_dataset_float(file_id_1,dims,fd, "/fd",0);
	add_eulerian_to_XDMF("Rheology",fl1out,"/fd",0,nodenum,fl);
	}
/**/
if(fl0stp2[4] || fl0stp2[5] || fl0stp2[6] || fl0stp2[8] || fl0stp2[12])
	{
	gridcheck();
	vx=malloc(sizeof(float) * nodenum);
	vy=malloc(sizeof(float) * nodenum);
	vz=malloc(sizeof(float) * nodenum);
	/* Read marker coordinates from 3 * marknum matrix */
        ndims[1]=3;
        ndims[0]=dims;
        dataset_id = H5Dopen(group_id, "/Nodes/Velocity", H5P_DEFAULT);
        count[1] =1;
        count[0] =dims;
        memspace = H5Screate_simple(2, count, NULL);
        offset[0] = 0;
        offset[1] = 0;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vx);
        H5Sclose(dataspace_id);
        offset[1] = 1;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vy);
        H5Sclose(dataspace_id);
        offset[1] = 2;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vz);
        H5Sclose(dataspace_id);
        H5Sclose(memspace);
        H5Dclose(dataset_id);
	/**/
if(fl0stp2[4]) 
	{
	int mgi=0;
	long int mcmax1,wn1[100];
	double vxyz1[3];
	vx1=malloc(sizeof(float) * nodenum);
	vy1=malloc(sizeof(float) * nodenum);
	vz1=malloc(sizeof(float) * nodenum);
	for (m3=0;m3<mgz[mgi];m3++)
	for (m1=0;m1<mgx[mgi];m1++)
	for (m2=0;m2<mgy[mgi];m2++)
        	{
        	/* Pos in vx[], vy[], vz[], pr[], etc */
        	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
		vxyzcalc1(gx[m1],gy[m2],gz[m3],vxyz1,wn1);
        	vx1[mcmax1]=vxyz1[0];	
        	vy1[mcmax1]=vxyz1[1];	
        	vz1[mcmax1]=vxyz1[2];	
		}
	ndims[1]=3;
        ndims[0]=dims;
        dataspace_id=H5Screate_simple(2,ndims,NULL);
        dataset_id = H5Dcreate(file_id_1, "/Velocity", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Sclose(dataspace_id);
        count[1] =1;
        count[0] =dims;
        memspace = H5Screate_simple(2, count, NULL);
        offset[0] = 0;
        offset[1] = 0;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vx1);
        H5Sclose(dataspace_id);
	for(m1=0;m1<nodenum;m1++){
	vy1[m1]=-vy1[m1];
        	}
        offset[1] = 1;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vy1);
        H5Sclose(dataspace_id);
        offset[1] = 2;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vz1);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Sclose(memspace);
	/**/
	add_eulerian_to_XDMF("Velocity",fl1out,"/Velocity",1,nodenum,fl);
	}
	/**/
if(fl0stp2[5] || fl0stp2[6] || fl0stp2[8] || fl0stp2[12]) 
	{
	int m1,m2,m3,mgi=0;
	long int mcmax1,wn1[500];
	long int un1[200];
	double ui1[200],eps1[100];
	/* Recalculate stress-straini rate components */
	nu=malloc(sizeof(float) * nodenum);
	nxx=malloc(sizeof(float) * nodenum);
	nxy=malloc(sizeof(float) * nodenum);
	nxz=malloc(sizeof(float) * nodenum);
	nyz=malloc(sizeof(float) * nodenum);
	exx=malloc(sizeof(float) * nodenum);
	eyy=malloc(sizeof(float) * nodenum);
	ezz=malloc(sizeof(float) * nodenum);
	exy=malloc(sizeof(float) * nodenum);
	exz=malloc(sizeof(float) * nodenum);
	eyz=malloc(sizeof(float) * nodenum);
	sxx=malloc(sizeof(float) * nodenum);
	syy=malloc(sizeof(float) * nodenum);
	szz=malloc(sizeof(float) * nodenum);
	sxy=malloc(sizeof(float) * nodenum);
	sxz=malloc(sizeof(float) * nodenum);
	syz=malloc(sizeof(float) * nodenum);
	wxy=malloc(sizeof(float) * nodenum);
	wxz=malloc(sizeof(float) * nodenum);
	wyz=malloc(sizeof(float) * nodenum);
	loadsave_dataset_float(group_id,dims,nu, "/Nodes/nu",1);
	loadsave_dataset_float(group_id,dims,nxx, "/Nodes/nxx",1);
	loadsave_dataset_float(group_id,dims,nxy, "/Nodes/nxy",1);
	loadsave_dataset_float(group_id,dims,nxz, "/Nodes/nxz",1);
	loadsave_dataset_float(group_id,dims,nyz, "/Nodes/nyz",1);
	nubeg=1e+0;	
        nuend=1e+40;
	/* Load pressure */
	pr=malloc(sizeof(double) * nodenum);
	loadsave_dataset_double(group_id,dims,pr, "/Nodes/pr",1);
	for (m3=0;m3<mgz[mgi];m3++)
	for (m1=0;m1<mgx[mgi];m1++)
	for (m2=0;m2<mgy[mgi];m2++)
        	{
        	/* Pos in vx[], vy[], vz[], pr[], etc */
        	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
		/**/
        	if(m1>0 && m2>0 && m3>0)
                	{
                	/* Sxx,Exx */
                	sxx[mcmax1]=sxxcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exx[mcmax1]=eps1[0];
                	/**/
                	/* Syy,Eyy */
                	syy[mcmax1]=syycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyy[mcmax1]=eps1[0];
                	/**/
                	/* Szz,Ezz */
                	szz[mcmax1]=szzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); ezz[mcmax1]=eps1[0];
                	}
	        /**/
		/* Sxy,Exy */
	        if(xperiodic && !yperiodic)
			{
			if(m2>0 && m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
	                	{       
	                	if(m1==mgx[mgi]-1) {sxy[mcmax1]=sxycalc(0,m2,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];}
	                	else {sxy[mcmax1]=sxycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];}
	                	}       
	        	}
	        else if(yperiodic && !xperiodic)
			{
			if(m1>0 && m1<mgx[mgi]-1 && m3<mgz[mgi]-1)
	        	        {       
	        	        if(m2==mgy[mgi]-1) {sxy[mcmax1]=sxycalc(m1,0,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];}
	        	        else {sxy[mcmax1]=sxycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];}
	        	        }       
	        	}
	        else if(xperiodic && yperiodic)
			{
			if(m3<mgz[mgi]-1)
	        	        {       
	        	        if(m1==mgx[mgi]-1 && m2==mgy[mgi]-1) {sxy[mcmax1]=sxycalc(0,0,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];}
	        	        else if(m1==mgx[mgi]-1) {sxy[mcmax1]=sxycalc(0,m2,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];}
	        	        else if(m2==mgy[mgi]-1) {sxy[mcmax1]=sxycalc(m1,0,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];}
	        	        else {sxy[mcmax1]=sxycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];}
	        	        }       
	        	}
		else
			{
			if(m1>0 && m1<mgx[mgi]-1 && m2>0 && m2<mgy[mgi]-1 && m3<mgz[mgi]-1) 
				{
				sxy[mcmax1]=sxycalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exy[mcmax1]=eps1[0]; wxy[mcmax1]=eps1[1];
				}
			}
		/* Sxz,Exz */
	        if(xperiodic && !zperiodic)
			{
			if(m3>0 && m3<mgz[mgi]-1 && m2<mgy[mgi]-1)
	                	{       
	                	if(m1==mgx[mgi]-1) {sxz[mcmax1]=sxzcalc(0,m2,m3,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1];}
	                	else {sxz[mcmax1]=sxzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1];}
	                	}       
	        	}
	        else if(zperiodic && !xperiodic)
			{
			if(m1>0 && m1<mgx[mgi]-1 && m2<mgy[mgi]-1)
	        	        {       
	        	        if(m3==mgz[mgi]-1) {sxz[mcmax1]=sxzcalc(m1,m2,0,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1];}
	        	        else {sxz[mcmax1]=sxzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1];}
	        	        }       
	        	}
	        else if(xperiodic && zperiodic)
			{
			if(m2<mgy[mgi]-1)
	        	        {       
	        	        if(m1==mgx[mgi]-1 && m3==mgz[mgi]-1) {sxz[mcmax1]=sxzcalc(0,m2,0,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1];}
	        	        else if(m1==mgx[mgi]-1) {sxz[mcmax1]=sxzcalc(0,m2,m3,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1];}
	        	        else if(m3==mgz[mgi]-1) {sxz[mcmax1]=sxzcalc(m1,m2,0,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1];}
	        	        else {sxz[mcmax1]=sxzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1];}
	        	        }       
	        	}
		else
			{
			if(m1>0 && m1<mgx[mgi]-1 && m3>0 && m3<mgz[mgi]-1 && m2<mgy[mgi]-1) 
				{
				sxz[mcmax1]=sxzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); exz[mcmax1]=eps1[0]; wxz[mcmax1]=eps1[1]; 
				}
			}
		/* Syz,Eyz */	
	        if(yperiodic && !zperiodic)
			{
			if(m3>0 && m3<mgz[mgi]-1 && m1<mgx[mgi]-1)
	                	{       
	                	if(m2==mgy[mgi]-1) {syz[mcmax1]=syzcalc(m1,0,m3,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];}
	                	else {syz[mcmax1]=syzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];}
	                	}       
	        	}
	        else if(zperiodic && !yperiodic)
			{
			if(m2>0 && m2<mgy[mgi]-1 && m1<mgx[mgi]-1)
	        	        {       
	        	        if(m3==mgz[mgi]-1) {syz[mcmax1]=syzcalc(m1,m2,0,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];}
	        	        else {syz[mcmax1]=syzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];}
		//if(m3==0){printf("A %d %d %d %e %e %e %e %e %e\n",m1,m2,m3,gy[m2],nu[mcmax1],sxx[mcmax1],sxy[mcmax1],sxz[mcmax1],syz[mcmax1]); getchar();}
	        	        }       
	        	}
	        else if(yperiodic && zperiodic)
			{
			if(m1<mgx[mgi]-1)
	        	        {       
	        	        if(m2==mgy[mgi]-1 && m3==mgz[mgi]-1) {syz[mcmax1]=syzcalc(m1,0,0,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];}
	        	        else if(m2==mgy[mgi]-1) {syz[mcmax1]=syzcalc(m1,0,m3,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];}
	        	        else if(m3==mgz[mgi]-1) {syz[mcmax1]=syzcalc(m1,m2,0,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];}
	        	        else {syz[mcmax1]=syzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];}
	        	        }       
	        	}
		else
			{
			if(m2>0 && m2<mgy[mgi]-1 && m3>0 && m3<mgz[mgi]-1 && m1<mgx[mgi]-1) 
				{
				syz[mcmax1]=syzcalc(m1,m2,m3,0,mgi,eps1,un1,ui1); eyz[mcmax1]=eps1[0]; wyz[mcmax1]=eps1[1];
				}
			}
	/**/
		//if(m3==0){printf("%d %d %d %e %e %e %e %e %e\n",m1,m2,m3,gy[m2],nu[mcmax1],sxx[mcmax1],sxy[mcmax1],sxz[mcmax1],syz[mcmax1]); getchar();}
		}
	/**/
	if(fl0stp2[5] || fl0stp2[6] || fl0stp2[12])
	{
	int initx=1,inity=1,initz=1,a;
	double c;
	sii=malloc(sizeof(float) * nodenum);
	eii=malloc(sizeof(float) * nodenum);
	wii=malloc(sizeof(float) * nodenum);
	if(xperiodic) initx=0; 
	if(yperiodic) inity=0; 
	if(zperiodic) initz=0; 
	/* Strain rate invariant */ 
	for (m3=initz;m3<mgz[mgi];m3++)
	for (m1=initx;m1<mgx[mgi];m1++)
	for (m2=inity;m2<mgy[mgi];m2++)
        	{
        	/* Pos in vx[], vy[], vz[], pr[], etc */
        	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
		/**/
		epscalc1(gx[m1],gy[m2],gz[m3],0,eps1,wn1);
        	eii[mcmax1]=pow(0.5*(eps1[4]*eps1[4]+eps1[5]*eps1[5]+eps1[6]*eps1[6])+eps1[7]*eps1[7]+eps1[8]*eps1[8]+eps1[9]*eps1[9],0.5);
        	sii[mcmax1]=pow(0.5*(eps1[54]*eps1[54]+eps1[55]*eps1[55]+eps1[56]*eps1[56])+eps1[57]*eps1[57]+eps1[58]*eps1[58]+eps1[59]*eps1[59],0.5);
		wii[mcmax1]=pow(eps1[37]*eps1[37]+eps1[38]*eps1[38]+eps1[39]*eps1[39],0.5);
		/*
        	        if(m3==0 && m1==0){printf("%ld %ld %e %e %e %e %e %e %e \n",m2,mcmax1,exx[mcmax1],eyy[mcmax1],ezz[mcmax1],exy[mcmax1],exz[mcmax1],eyz[mcmax1],eii[mcmax1]); getchar();}
		if(wii[mcmax1]){printf("%d %d %d %e %e %e %e \n",m1,m2,m3,wii[mcmax1],eps1[37],eps1[38],eps1[39]); getchar();}
		*/
		//if(m2==0){printf("%d %d %d %e %e %e %e %e %e %e %e %e\n",m1,m2,m3,gy[m2],nu[mcmax1],sii[mcmax1],eii[mcmax1],sii[mcmax1]*eii[mcmax1],wii[mcmax1],eps1[37],eps1[38],eps1[39]); getchar();}
		/*
		a=-1; if(eps1[4]<0) a=1 ; c=ABSV(eps1[4]); exx[mcmax1]=a*log10(c);
		printf("%e %d %e %e\n",eps1[4],a,c,exx[mcmax1]); getchar();
		a=-1; if(eps1[5]<0) a=1 ; c=ABSV(eps1[5]); eyy[mcmax1]=a*log10(c);
		a=-1; if(eps1[6]<0) a=1 ; c=ABSV(eps1[6]); ezz[mcmax1]=a*log10(c);
		a=-1; if(eps1[7]<0) a=1 ; c=ABSV(eps1[7]); exy[mcmax1]=a*log10(c);
		a=-1; if(eps1[8]<0) a=1 ; c=ABSV(eps1[8]); exz[mcmax1]=a*log10(c);
		a=-1; if(eps1[9]<0) a=1 ; c=ABSV(eps1[9]); eyz[mcmax1]=a*log10(c);
		exx[mcmax1]=eps1[4];
		eyy[mcmax1]=eps1[5];
		ezz[mcmax1]=eps1[6];
		exy[mcmax1]=eps1[7];
		exz[mcmax1]=eps1[8];
		eyz[mcmax1]=eps1[9];
		*/
		}
	/* Make boundaries layers equal to adjacent ones */
	for (m3=0;m3<mgz[mgi];m3++)
	for (m1=0;m1<mgx[mgi];m1++)
	for (m2=0;m2<mgy[mgi];m2++)
        	{
		/**/
		/* Left boundary */
		if(!xperiodic && m1==0 ) 
			{
        		mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
			eii[mcmax1]=eii[mcmax1+mgy[mgi]];
			sii[mcmax1]=sii[mcmax1+mgy[mgi]];
			wii[mcmax1]=wii[mcmax1+mgy[mgi]];
			}
		/* Right boundary */
		if(!xperiodic && m1==xnumx1 ) 
			{
        		mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
			eii[mcmax1]=eii[mcmax1-mgy[mgi]];
			sii[mcmax1]=sii[mcmax1-mgy[mgi]];
			wii[mcmax1]=wii[mcmax1-mgy[mgi]];
			}
        	}
	for (m3=0;m3<mgz[mgi];m3++)
        for (m1=0;m1<mgx[mgi];m1++)
        for (m2=0;m2<mgy[mgi];m2++)
                {
                /**/
		/* Top boundary */
		if(!yperiodic && m2==0 ) 
			{
        		mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
			eii[mcmax1]=eii[mcmax1+1];
			sii[mcmax1]=sii[mcmax1+1];
			wii[mcmax1]=wii[mcmax1+1];
			}
		if(1==1){
		/* Bottom boundary */
		if(!yperiodic && m2==ynumy1) 
			{
        		mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
			eii[mcmax1]=eii[mcmax1-1];
			sii[mcmax1]=sii[mcmax1-1];
			wii[mcmax1]=wii[mcmax1-1];
			}
			}
		}
        for (m3=0;m3<mgz[mgi];m3++)
        for (m1=0;m1<mgx[mgi];m1++)
        for (m2=0;m2<mgy[mgi];m2++)
                {
                /**/
		/* Front boundary */
		if(!zperiodic && m3==0) 
			{
        		mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
			eii[mcmax1]=eii[mcmax1+mgxy[mgi]];
			sii[mcmax1]=sii[mcmax1+mgxy[mgi]];
			wii[mcmax1]=wii[mcmax1+mgxy[mgi]];
			}
		if(1==1){
		/* Back boundary */
		if(!zperiodic && m3==znumz1) 
			{
        		mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
			eii[mcmax1]=eii[mcmax1-mgxy[mgi]];
			sii[mcmax1]=sii[mcmax1-mgxy[mgi]];
			wii[mcmax1]=wii[mcmax1-mgxy[mgi]];
			}
			}
		}
	if(fl0stp2[5])
		{
		loadsave_dataset_float(file_id_1,dims,eii, "/eii",0);
		add_eulerian_to_XDMF("Eii",fl1out,"/eii",0,nodenum,fl);
		/*
		loadsave_dataset_float(file_id_1,dims,exx, "/exx",0);
		add_eulerian_to_XDMF("Exx",fl1out,"/exx",0,nodenum,fl);
		loadsave_dataset_float(file_id_1,dims,eyy, "/eyy",0);
		add_eulerian_to_XDMF("Eyy",fl1out,"/eyy",0,nodenum,fl);
		loadsave_dataset_float(file_id_1,dims,ezz, "/ezz",0);
		add_eulerian_to_XDMF("Ezz",fl1out,"/ezz",0,nodenum,fl);
		loadsave_dataset_float(file_id_1,dims,exy, "/exy",0);
		add_eulerian_to_XDMF("Exy",fl1out,"/exy",0,nodenum,fl);
		loadsave_dataset_float(file_id_1,dims,exz, "/exz",0);
		add_eulerian_to_XDMF("Exz",fl1out,"/exz",0,nodenum,fl);
		loadsave_dataset_float(file_id_1,dims,eyz, "/eyz",0);
		add_eulerian_to_XDMF("Eyz",fl1out,"/eyz",0,nodenum,fl);
		*/
		}
	if(fl0stp2[6])
		{
		loadsave_dataset_float(file_id_1,dims,sii, "/sii",0);
		add_eulerian_to_XDMF("Sii",fl1out,"/sii",0,nodenum,fl);
		}
	if(fl0stp2[12])
                {
                loadsave_dataset_float(file_id_1,dims,wii, "/wii",0);
                add_eulerian_to_XDMF("Wii",fl1out,"/wii",0,nodenum,fl);
                }
        }
	/**/
	if(fl0stp2[8])
	{
	double epsval,x,y,z,xkf,ykf,zkf;
	long int v[7],p[8],mcmax1;
	hs=malloc(sizeof(float) * nodenum);
	tk=malloc(sizeof(float) * nodenum);
	loadsave_dataset_float(group_id,dims,tk, "/Nodes/tk",1);
	for (m3=0;m3<mgz[mgi];m3++)
	for (m1=0;m1<mgx[mgi];m1++)
	for (m2=0;m2<mgy[mgi];m2++)
        	{
        	/* Pos in vx[], vy[], vz[], pr[], etc */
        	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
		/**/
		hs[mcmax1]=-20;
        	/*     2 6   */
		/*     |/    */
		/* 1---3---5 */
		/*    /|     */
		/*   0 4     */
		v[3]=mcmax1;
		v[0]=v[3]-mgxy[mgi];
		if(xperiodic && m1==0) v[1]=mgp[mgi]+m3*mgxy[mgi]+(mgx[mgi]-2)*mgy[mgi]+m2;
		else v[1]=v[3]-mgy[mgi];
		v[2]=v[3]-1;
		v[4]=v[3]+1;
		v[5]=v[3]+mgy[mgi];
		v[6]=v[3]+mgxy[mgi];
		if(tk[mcmax1]<=fl0stp2[8] && m1>0 && m2>0 && m3>0 && m1<mgx[mgi]-1 && m2<mgy[mgi]-1 && m3<mgz[mgi]-1)
			{
			if(xperiodic && m1==0) xkf=2.0/(mggx[mgi][1]+mggx[mgi][mgx[mgi]-1]-mggx[mgi][mgx[mgi]-2]);
			else xkf=2.0/(mggx[mgi][m1+1]-mggx[mgi][m1-1]);
			ykf=2.0/(mggy[mgi][m2+1]-mggy[mgi][m2-1]);
			zkf=2.0/(mggz[mgi][m3+1]-mggz[mgi][m3-1]);
			/* SIGxx*EPSxx+SIGyy*EPSyy+SIGzz*EPSzz */
	                /* SIGxx-SIGyy-SIGzz-Nodes num */
	                /*    4----6 */
	                /*   /|   /| */
	                /*  / 5--/-7 */
	                /* 0----2    */
	                /* |    |    */
	                /* 1----3    */
	                p[0]=v[3];
	                p[1]=p[0]+1;
	                p[2]=p[0]+mgy[mgi];
	                p[3]=p[2]+1;
	                p[4]=p[0]+mgxy[mgi];
	                p[5]=p[4]+1;
	                p[6]=p[4]+mgy[mgi];
	                p[7]=p[6]+1;
	                /* Relative distances */
	                if(xperiodic && m1==0) x=(mggx[mgi][mgx[mgi]-1]-mggx[mgi][mgx[mgi]-2])/2.0*xkf;
	                else x=(mggx[mgi][m1]-mggx[mgi][m1-1])/2.0*xkf;
	                y=(mggy[mgi][m2]-mggy[mgi][m2-1])/2.0*ykf;
	                z=(mggz[mgi][m3]-mggz[mgi][m3-1])/2.0*zkf;
			epsval =(1.0-x)*(1.0-y)*(1.0-z)*(sxx[p[0]]*exx[p[0]]+syy[p[0]]*eyy[p[0]]+szz[p[0]]*ezz[p[0]]);
	                epsval+=(1.0-x)*(    y)*(1.0-z)*(sxx[p[1]]*exx[p[1]]+syy[p[1]]*eyy[p[1]]+szz[p[1]]*ezz[p[1]]);
	                epsval+=(    x)*(1.0-y)*(1.0-z)*(sxx[p[2]]*exx[p[2]]+syy[p[2]]*eyy[p[2]]+szz[p[2]]*ezz[p[2]]);
	                epsval+=(    x)*(    y)*(1.0-z)*(sxx[p[3]]*exx[p[3]]+syy[p[3]]*eyy[p[3]]+szz[p[3]]*ezz[p[3]]);
	                epsval+=(1.0-x)*(1.0-y)*(    z)*(sxx[p[4]]*exx[p[4]]+syy[p[4]]*eyy[p[4]]+szz[p[4]]*ezz[p[4]]);
	                epsval+=(1.0-x)*(    y)*(    z)*(sxx[p[5]]*exx[p[5]]+syy[p[5]]*eyy[p[5]]+szz[p[5]]*ezz[p[5]]);
	                epsval+=(    x)*(1.0-y)*(    z)*(sxx[p[6]]*exx[p[6]]+syy[p[6]]*eyy[p[6]]+szz[p[6]]*ezz[p[6]]);
	                epsval+=(    x)*(    y)*(    z)*(sxx[p[7]]*exx[p[7]]+syy[p[7]]*eyy[p[7]]+szz[p[7]]*ezz[p[7]]);
	                /**/
	                /* 2*(SIGxy*EPSxy+SIGxz*EPSxz+SIGyz*EPSyz)*/
	                /* SIGxy-SIGxz-SIGyz-T-Nodes num */
	                /*     SIGxz2 SIGxy3  */
	                /*          |/        */
	                /* SIGyz1--T3--SIGyz3 */
	                /*         /|         */
	                /*   SIGxy0 SIGxz3    */
	                /* SIGxy*EPSxy */
                	epsval+=2.0*((1.0-z)*(sxy[v[0]]*exy[v[0]])+z*(sxy[v[3]]*exy[v[3]]));
                	/* SIGxz*EPSxz */
                	epsval+=2.0*((1.0-y)*(sxz[v[2]]*exz[v[2]])+y*(sxz[v[3]]*exz[v[3]]));
        	        /* SIGyz*EPSyz */
	                epsval+=2.0*((1.0-x)*(syz[v[1]]*eyz[v[1]])+x*(syz[v[3]]*eyz[v[3]]));
			/* Dissipated heat in log10(milliWatt/m^3) */
			hs[mcmax1]=log10(epsval)+6;
			//printf("Hs: %ld %e %e \n",mcmax1,tk[mcmax1],hs[mcmax1]); getchar();
			}
		}
	dims=nodenum;
	loadsave_dataset_float(file_id_1,dims,hs, "/hs",0);
	add_eulerian_to_XDMF("Hs",fl1out,"/hs",0,nodenum,fl);
	}
	}
	}
/**/
if(fl0stp2[7])
	{
	int m1,m2,m3,mgi=0;
	long int mcmax1,wn1[500];
	long int un1[MAXPOS],pos0cur1;
	double ui1[MAXPOS],eps1[100],GYKOEF=9.80665;
	pr=malloc(sizeof(double) * nodenum);
	pr0=malloc(sizeof(double) * nodenum);
	ro=malloc(sizeof(double) * nodenum);
	double *prr=malloc(sizeof(double) * ynumy);
	/* Load pressure */
	loadsave_dataset_double(group_id,dims,pr, "/Nodes/pr",1);
	loadsave_dataset_float(group_id,dims,ro, "/Nodes/ro",1);
	/* Calculate reference lithostatic pressure profile */
	prr[0]=0.0;
	for(m2=1;m2<mgy[mgi];m2++)
		{
		mcmax1=mgp[mgi]+(znumz-2)*mgxy[mgi]+(xnumx-2)*mgy[mgi]+(m2-1);
        	if(m2==1) prr[m2]=prr[m2-1]+GYKOEF*(ro[mcmax1]+ro[mcmax1-mgy[mgi]])/2*(mggy[mgi][m2]-mggy[mgi][m2-1])/2;
        	else prr[m2]=prr[m2-1]+GYKOEF*(ro[mcmax1]+ro[mcmax1-mgy[mgi]])/2*(mggy[mgi][m2]-mggy[mgi][m2-2])/2;
		/*
		printf("%d %ld %e\n",m2,mcmax1,prr[m2]/1e+6); getchar();
		*/
		}
	for (m3=0;m3<mgz[mgi];m3++)
	for (m1=0;m1<mgx[mgi];m1++)
	for (m2=0;m2<mgy[mgi];m2++)
        	{
        	/* Pos in vx[], vy[], vz[], pr[], etc */
        	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
		/**/
		//pr[mcmax1]-=prr[m2];
		/*
		printf("%d %ld %e\n",m2,mcmax1,pr0[m2]/1e+6); getchar();
		*/
		}
	/**/
	if(1==1){
	for (m3=1;m3<mgz[mgi];m3++)
	for (m1=1;m1<mgx[mgi];m1++)
	for (m2=1;m2<mgy[mgi];m2++)
        	{
        	/* Pos in vx[], vy[], vz[], pr[], etc */
        	mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
		/**/
		epscalc1(gx[m1],gy[m2],gz[m3],0,eps1,wn1);
		pr0[mcmax1]=eps1[10];
		}
		}
	/* Make boundaries layers equal to adjacent ones */
        for (m3=0;m3<mgz[mgi];m3++)
        for (m1=0;m1<mgx[mgi];m1++)
        for (m2=0;m2<mgy[mgi];m2++)
                {
                /* Pos in vx[], vy[], vz[], pr[], etc */
                mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
                /**/
               	/* Left boundary */
		if(!xperiodic && m1==0 && m2>0 && m2<ynumy1 && m3>0 && m3<znumz1) pr0[mcmax1]=pr0[mcmax1+mgy[mgi]];
               	/* Right boundary */
               	if(!xperiodic && m1==xnumx1 && m2>0 && m2<ynumy1 && m3>0 && m3<znumz1) pr0[mcmax1]=pr0[mcmax1-mgy[mgi]];
		}
        for (m3=0;m3<mgz[mgi];m3++)
        for (m1=0;m1<mgx[mgi];m1++)
        for (m2=0;m2<mgy[mgi];m2++)
                {
                /* Pos in vx[], vy[], vz[], pr[], etc */
                mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
                /**/
		/* Top boundary */
                if(m2==0 && m3>0 && m3<znumz1) pr0[mcmax1]=pr0[mcmax1+1];
                /* Bottom boundary */
                if(1==0 && m2==ynumy1 && m3>0 && m3<znumz1) pr0[mcmax1]=pr0[mcmax1-1];
        	}
	for (m3=0;m3<mgz[mgi];m3++)
        for (m1=0;m1<mgx[mgi];m1++)
        for (m2=0;m2<mgy[mgi];m2++)
                {
                /* Pos in vx[], vy[], vz[], pr[], etc */
                mcmax1=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
                /**/
                /* Front boundary */
                if(m3==0) pr0[mcmax1]=pr0[mcmax1+mgxy[mgi]];
                /* Back boundary */
                if(1==0 && m3==znumz1) pr0[mcmax1]=pr0[mcmax1-mgxy[mgi]];
		}
	loadsave_dataset_double(file_id_1,dims,pr0, "/pr",0);
	add_eulerian_to_XDMF("Pressure",fl1out,"/pr",0,nodenum,fl);
	}
if(fl0stp2[10]==0)
	{
	fprintf(fl,"   </Grid>\n\n");
	}
H5Gclose(group_id);
}
/**/
/* Lagrangian particles */
if(fl0stp2[10]||fl0stp2[11])
{
/* Lagrangian particles */
group_id = H5Gopen(file_id_0, "/Particles", H5P_DEFAULT);
/* Creating attribute */
/* Read /Nodes attribute */
long int b[4];
attr_id=H5Aopen_by_name(group_id,".","/Particles/Number",H5P_DEFAULT,H5P_DEFAULT);
H5Aread (attr_id, H5T_NATIVE_LONG, &b);
H5Aclose(attr_id);
marknum=b[0];
/*
printf("%ld %ld %ld %ld\n",mnumx,mnumy,mnumz,marknum);
 */
/**/
dynmemall(2);
/* Particle parameters */
dims=marknum;
/* Read marker coordinates from 3 * marknum matrix */
ndims[1]=3;
ndims[0]=dims;
dataset_id = H5Dopen(group_id, "/Particles/Position", H5P_DEFAULT);
count[1] =1;
count[0] =dims;
memspace = H5Screate_simple(2, count, NULL);
offset[0] = 0;
offset[1] = 0;
dataspace_id = H5Dget_space(dataset_id);
H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, markx);
H5Sclose(dataspace_id);
offset[1] = 1;
dataspace_id = H5Dget_space(dataset_id);
H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, marky);
H5Sclose(dataspace_id);
offset[1] = 2;
dataspace_id = H5Dget_space(dataset_id);
H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, markz);
H5Sclose(dataspace_id);
H5Sclose(memspace);
H5Dclose(dataset_id);
/* Interpolate cumulative strain to nodes */
if(fl0stp2[10])
	{
	long int mm1,m10,m20,m30;
	double dx,dy,dz,ddx,ddy,ddz,ival,swt;
	ee=malloc(sizeof(float) * nodenum);
	sol=malloc(sizeof(double) * nodenum);
	for(m1=0;m1<nodenum;m1++)
		{
		ee[m1]=0.0;
		sol[m1]=0.0;
		}
	/* Load cumulative strain on markers */
	loadsave_dataset_float(group_id,dims,marke,"/Particles/marke",1);
	for(mm1=0;mm1<marknum;mm1++)	
	if(marke[mm1]>0)
	{
	/* Cell num calc for A-D Diagonal Corners of cur marker */
	m10=m1serch(markx[mm1]);
	m20=m2serch(marky[mm1]);
	m30=m3serch(markz[mm1]);
        dx=gx[m10+1]-gx[m10];
	dy=gy[m20+1]-gy[m20];
	dz=gz[m30+1]-gz[m30];
	/**/
	/* Add nodes around current cell */
	ival=0;
	for (m3=m30;m3<=m30+1;m3++)
		{
		/* Normalized Z-distance calc */
		ddz=(gz[m3]-markz[mm1])/dz;
		ddz=1.0-ABSV(ddz);
		for (m1=m10;m1<=m10+1;m1++)
			{
			/* Normalized X-distance calc */
			ddx=(gx[m1]-markx[mm1])/dx;
			ddx=1.0-ABSV(ddx);
			for (m2=m20;m2<=m20+1;m2++)
				{
				/* Normalized Y-distance calc */
				ddy=(gy[m2]-marky[mm1])/dy;
				ddy=1.0-ABSV(ddy);
				/* Wt calc */
				swt=ddx*ddy*ddz;
/*
ival+=swt;
{printf("%ld  %ld %ld %ld  %e %e %e   %e %e %e %e  %e %e %e   %e %e",mm1,m1,m2,m3,markx[mm1],marky[mm1],markz[mm1],marke[mm1],dx,dy,dz,ddx,ddy,ddz,swt,ival);getchar();}
*/
				/**/
				/* Node num */
				m4=m3*xynumxy+m1*ynumy+m2;
				/**/
				/* Add viscosity */
			        if(swt>0)
					{
					ival=marke[mm1]*swt;
					ee[m4]+=ival;
					sol[m4]+=swt;
					}
				}
			}
		}
	}
	/* Normalization */
	for(m1=0;m1<nodenum;m1++)
	if(sol[m1]>0)
		{
		if(isinf(ee[m1]/sol[m1])){printf("Strain %ld %e %e\n",m1,ee[m1],sol[m1]); getchar();}
		ee[m1]/=sol[m1];
		}
	/* Save */
	dims=nodenum;
	loadsave_dataset_float(file_id_1,dims,ee, "/ee",0);
	add_eulerian_to_XDMF("Strain",fl1out,"/ee",0,nodenum,fl);
	fprintf(fl,"   </Grid>\n\n");
	}
/**/
if(fl0stp2[11])
	{
/* Marker type */
loadsave_dataset_char(file_id_0,dims, markt,"/Particles/markt",1);
/**/
dims=marknum;
/* Save marker coordinates as a 3 * marknum matrix */
ndims[1]=3;
ndims[0]=dims;
dataspace_id=H5Screate_simple(2,ndims,NULL);
dataset_id = H5Dcreate(file_id_1, "/Position", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
H5Sclose(dataspace_id);
count[1] =1;
count[0] =dims;
memspace = H5Screate_simple(2, count, NULL);
offset[0] = 0;
offset[1] = 0;
dataspace_id = H5Dget_space(dataset_id);
H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, markx);
H5Sclose(dataspace_id);
/*
*/
for(m1=0;m1<marknum;m1++){
marky[m1]=ysize-marky[m1];
        }
offset[1] = 1;
dataspace_id = H5Dget_space(dataset_id);
H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, marky);
H5Sclose(dataspace_id);
/*
*/
offset[1] = 2;
dataspace_id = H5Dget_space(dataset_id);
H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, markz);
H5Sclose(dataspace_id);
H5Dclose(dataset_id);
H5Sclose(memspace);
H5Gclose(group_id);
/* Lagrangian particles informations */
fprintf(fl,"   <Grid Name=\"Lagrangian particles\">\n\n");
fprintf(fl,"      <Time Value=\"%f\" />\n\n",timesum/365.25/86400);
fprintf(fl,"         <Topology Type=\"POLYVERTEX\" NodesPerElement=\"%ld\"> </Topology>\n",marknum);
fprintf(fl,"         <Geometry Type=\"XYZ\">\n");
fprintf(fl,"            <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">%s:/Position</DataItem>\n",marknum,fl1out);
fprintf(fl,"         </Geometry>\n\n");
/* Add Lagrangian fields  */
loadsave_dataset_char(file_id_1,dims, markt,"/markt",0);
add_lagrangian_to_XDMF("Marker Type",fl1out,"/markt",0,marknum,1,fl);
/**/
fprintf(fl,"    </Grid>\n\n");
	}
	}
fprintf(fl," </Domain>\n");
fprintf(fl,"</Xdmf>\n");
/* Close xmf file */
H5Fclose(file_id_0);
H5Fclose(file_id_1);
fclose(fl);
//Move file to output directory 
sprintf(command,"mv %s %s",fl1out,outdir);
system(command);
sprintf(command,"mv %s %s",fname,outdir);
system(command);
//printf("\n Move file %s to %s\n",fl1out,outdir);
}
/* End Make XMF file to visualize .h5 fields in Paraview -------------*/


/* Calculation of P, SIG and  EPSxx, EPSyy, EPSzz, EPSxy, EPSxz, EPSyz for current location by Linear Interpolation */
void epscalc1(double x, double y, double z, int yn, double *eps1, long int *wn1)
/* x,y,z - XYZ location of point for Vx,Vy,Vz calc */
/* yn - cell address available in wn[]  Y(1)/N(0) */
{
/* Vx,Vy-nodes num */
long int m1,m2,m3,m10,m20,m30,m40,m10min[2],m20min[2],m30min[2];
/* Relativ coord */
double X,Y,Z,dx,dy,dz,ddx,ddy,ddz,swt,ival;
/**/
/**/
/**/
/* Clear All EPS,SIG */
eps1[4]=eps1[5]=eps1[6]=eps1[7]=eps1[8]=eps1[9]=eps1[10]=0;
eps1[54]=eps1[55]=eps1[56]=eps1[57]=eps1[58]=eps1[59]=0;
eps1[37]=eps1[38]=eps1[39]=0;
/* Out of grid */
if(x<0) x=0; if(x>xsize) x=xsize;
if(y<0) y=0; if(y>ysize) y=ysize;
if(z<0) z=0; if(z>zsize) z=zsize;
/* Fing Cell indexes */
if(yn==0)
	{
	wn1[0]=m1=m1serch(x);
	wn1[1]=m2=m2serch(y);
	wn1[2]=m3=m3serch(z);
	}
else
	{
	m1=wn1[0];
	m2=wn1[1];
	m3=wn1[2];
	}
/**/
/**/
/**/
/* EPSxx, EPSyy, EPSzz Cell indexes */
m10min[0]=m1;
X=x;
if(X>(gx[m10min[0]]+gx[m10min[0]+1])/2.0) m10min[0]+=1;
if(xperiodic && (m10min[0]<1 || m10min[0]>xnumx-2))
	{
	if(m10min[0]<1) X=xsize+x;
        m10min[0]=xnumx-1;
	m10min[1]=1;
	dx=(gx[1]-gx[0]+gx[xnumx-1]-gx[xnumx-2])/2.0;
	}
else
	{
	if(m10min[0]<1) {m10min[0]=1; X=(gx[m10min[0]]+gx[m10min[0]-1])/2.0; }
        if(m10min[0]>xnumx-2) {m10min[0]=xnumx-2; X=(gx[m10min[0]]+gx[m10min[0]+1])/2.0;}
        m10min[1]=m10min[0]+1;
        dx=(gx[m10min[0]+1]-gx[m10min[0]-1])/2.0;
	}
m20min[0]=m2;
Y=y;
if(Y>(gy[m20min[0]]+gy[m20min[0]+1])/2.0) m20min[0]+=1;
if(yperiodic && (m20min[0]<1 || m20min[0]>ynumy-2))
        {
        if(m20min[0]<1) Y=ysize+y;
        m20min[0]=ynumy-1;
        m20min[1]=1;
        dy=(gy[1]-gy[0]+gy[ynumy-1]-gy[ynumy-2])/2.0;
	}
else
        {
        if(m20min[0]<1) {m20min[0]=1; Y=(gy[m20min[0]]+gy[m20min[0]-1])/2.0; }
        if(m20min[0]>ynumy-2) {m20min[0]=ynumy-2; Y=(gy[m20min[0]]+gy[m20min[0]+1])/2.0;}
        m20min[1]=m20min[0]+1;
	dy=(gy[m20min[0]+1]-gy[m20min[0]-1])/2.0;
        }
m30min[0]=m3;
Z=z;
if(Z>(gz[m30min[0]]+gz[m30min[0]+1])/2.0) m30min[0]+=1;
if(zperiodic && (m30min[0]<1 || m30min[0]>znumz-2))
        {
        if(m30min[0]<1) Z=zsize+z;
        m30min[0]=znumz-1;
        m30min[1]=1;
        dz=(gz[1]-gz[0]+gz[znumz-1]-gz[znumz-2])/2.0;
        }
else
        {
        if(m30min[0]<1) {m30min[0]=1; Z=(gz[m30min[0]]+gz[m30min[0]-1])/2.0; }
        if(m30min[0]>znumz-2) {m30min[0]=znumz-2; Z=(gz[m30min[0]]+gz[m30min[0]+1])/2.0;}
        m30min[1]=m30min[0]+1;
	dz=(gz[m30min[0]+1]-gz[m30min[0]-1])/2.0;
	}
/* EPS Cell dimensions */
/* Interpolate from 8 nodes */
ival=0;
for (m30=0;m30<=1;m30++)
	{
	/* Normalized Z-distance calc */
	if(zperiodic && m30min[0]==znumz-1 && m30==1) ddz=((gz[m30min[m30]]+gz[m30min[m30]-1])/2.0-(Z-zsize))/dz;
	else ddz=((gz[m30min[m30]]+gz[m30min[m30]-1])/2.0-Z)/dz;
	if (m30==0) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=0;m10<=1;m10++)
		{
		/* Normalized X-distance calc */
		if(xperiodic && m10min[0]==xnumx-1 && m10==1) ddx=((gx[m10min[m10]]+gx[m10min[m10]-1])/2.0-(X-xsize))/dx;
		else ddx=((gx[m10min[m10]]+gx[m10min[m10]-1])/2.0-X)/dx;
		if (m10==0) ddx=-ddx;
		ddx=1.0-ddx;
		//printf("%ld %ld %ld %ld %e %e %e %e %e\n",m1,m10,m10min[0],m10min[1],(gx[m10]+gx[m10-1])/2.0,x,X,dx,ddx);
		for (m20=0;m20<=1;m20++)
			{
			/* Normalized Y-distance calc */
			if(yperiodic && m20min[0]==ynumy-1 && m20==1) ddy=((gy[m20min[m20]]+gy[m20min[m20]-1])/2.0-(Y-ysize))/dy;
			else ddy=((gy[m20min[m20]]+gy[m20min[m20]-1])/2.0-Y)/dy;
			if (m20==0) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vx[] */
			m40=m30min[m30]*xynumxy+m10min[m10]*ynumy+m20min[m20];
			/* Calc wt */
			swt=ddx*ddy*ddz;
			if(ABSV(swt)<1e-10) swt=0.0;
//if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Exx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPS for current node */
			eps1[4]+=exx[m40]*swt;
			eps1[5]+=eyy[m40]*swt;
			eps1[6]+=ezz[m40]*swt;
			/* Add SIG for current node */
			eps1[54]+=sxx[m40]*swt;
			eps1[55]+=syy[m40]*swt;
			eps1[56]+=szz[m40]*swt;
			/* Add P   for current node */
			eps1[10]+=pr[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/*
if(ABSV(ival-1.0)>1e-6){printf("EPS  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[4]);getchar();}
*/
/**/
/* EPSxy Cell indexes */
m10min[0]=m1;
X=x;
if(xperiodic)
	{
	if(m10min[0]<0) {m10min[0]=0; X=0;}
	if(m10min[0]>xnumx-2) {m10min[0]=xnumx-2; X=gx[m10min[0]+1];}
	}
else
	{
	if(m10min[0]<1) {m10min[0]=1; X=gx[m10min[0]];}
	if(m10min[0]>xnumx-3) {m10min[0]=xnumx-3; X=gx[m10min[0]+1];}
	}
m20min[0]=m2;
Y=y;
if(yperiodic)
	{
	if(m20min[0]<0) {m20min[0]=0; Y=gy[m20min[0]];}
	if(m20min[0]>ynumy-2) {m20min[0]=ynumy-2; Y=gy[m20min[0]+1];}
	}
else
	{
	if(m20min[0]<1) {m20min[0]=1; Y=gy[m20min[0]];}
	if(m20min[0]>ynumy-3) {m20min[0]=ynumy-3; Y=gy[m20min[0]+1];}
	}
m30min[0]=m3;
Z=z;
if(Z<(gz[m30min[0]]+gz[m30min[0]+1])/2.0) m30min[0]-=1;
if(zperiodic && (m30min[0]<0 || m30min[0]>znumz-2))
        {
        if(m30min[0]<0) Z=zsize+z;
        m30min[0]=znumz-2;
        m30min[1]=0;
        dz=(gz[1]-gz[0]+gz[znumz-1]-gz[znumz-2])/2.0;
	}
else
        {
        if(m30min[0]<1) {m30min[0]=0; Z=(gz[m30min[0]]+gz[m30min[0]+1])/2.0; }
        if(m30min[0]>znumz-3) {m30min[0]=znumz-3; Z=(gz[m30min[0]+1]+gz[m30min[0]+2])/2.0;}
	dz=(gz[m30min[0]+2]-gz[m30min[0]])/2.0;
	m30min[1]=m30min[0]+1;
        }
/* EPSxy Cell dimensions */
dx= gx[m10min[0]+1]-gx[m10min[0]];
dy= gy[m20min[0]+1]-gy[m20min[0]];
/* Interpolate from 8 nodes */
ival=0;
for (m30=0;m30<=1;m30++)
	{
	/* Normalized Z-distance calc */
	if(zperiodic && m30min[0]==znumz-2 && m30==1) ddz=((gz[m30min[m30]]+gz[m30min[m30]+1])/2.0-(Z-zsize))/dz;
	else ddz=((gz[m30min[m30]]+gz[m30min[m30]+1])/2.0-Z)/dz;
	if (m30==0) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=m10min[0];m10<=m10min[0]+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=(gx[m10]-X)/dx;
		if (m10==m10min[0]) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=m20min[0];m20<=m20min[0]+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=(gy[m20]-Y)/dy;
			if (m20==m20min[0]) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vx[] */
			m40=m30min[m30]*xynumxy+m10*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
//if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Exy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPSxy for current node */
			eps1[7]+=exy[m40]*swt;
			/* Add SIGxy for current node */
			eps1[57]+=sxy[m40]*swt;
			/* Add Vorxy for current node */
			eps1[37]+=wxy[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/*
if(ABSV(ival-1.0)>1e-6){printf("EPSxy %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[7]);getchar();}
*/
/**/
/* EPSxz Cell indexes */
m10min[0]=m1;
X=x;
if(xperiodic)
	{
	if(m10min[0]<0) {m10min[0]=0; X=0;}
	if(m10min[0]>xnumx-2) {m10min[0]=xnumx-2; X=gx[m10min[0]+1];}
	}
else
	{
	if(m10min[0]<1) {m10min[0]=1; X=gx[m10min[0]];}
	if(m10min[0]>xnumx-3) {m10min[0]=xnumx-3; X=gx[m10min[0]+1];}
	}
m20min[0]=m2;
Y=y;
if(Y<(gy[m20min[0]]+gy[m20min[0]+1])/2.0) m20min[0]-=1;
if(yperiodic && (m20min[0]<0 || m20min[0]>ynumy-2))
        {
        if(m20min[0]<0) Y=ysize+y;
        m20min[0]=ynumy-2;
        m20min[1]=0;
        dy=(gy[1]-gy[0]+gy[ynumy-1]-gy[ynumy-2])/2.0;
	}
else
        {
        if(m20min[0]<1) {m20min[0]=0; Y=(gy[m20min[0]]+gy[m20min[0]+1])/2.0; }
        if(m20min[0]>ynumy-3) {m20min[0]=ynumy-3; Y=(gy[m20min[0]+1]+gy[m20min[0]+2])/2.0;}
	dy=(gy[m20min[0]+2]-gy[m20min[0]])/2.0;
	m20min[1]=m20min[0]+1;
        }
m30min[0]=m3;
Z=z;
if(zperiodic)
	{
	if(m30min[0]<0) {m30min[0]=0; Z=0;}
	if(m30min[0]>znumz-2) {m30min[0]=znumz-2; Z=gz[m30min[0]+1];}
	}
else
	{
	if(m30min[0]<1) {m30min[0]=1; Z=gz[m30min[0]];}
	if(m30min[0]>znumz-3) {m30min[0]=znumz-3; Z=gz[m30min[0]+1];}
	}
/* EPSxz Cell dimensions */
dx= gx[m10min[0]+1]-gx[m10min[0]];
dz= gz[m30min[0]+1]-gz[m30min[0]];
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min[0];m30<=m30min[0]+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=(gz[m30]-Z)/dz;
	if (m30==m30min[0]) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=m10min[0];m10<=m10min[0]+1;m10++)
		{
		/* Normalized X-distance calc */
		ddx=(gx[m10]-X)/dx;
		if (m10==m10min[0]) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=0;m20<=1;m20++)
			{
			/* Normalized Y-distance calc */
			if(yperiodic && m20min[0]==ynumy-2 && m20==1) ddy=((gy[m20min[m20]]+gy[m20min[m20]+1])/2.0-(Y-ysize))/dy;
			else ddy=((gy[m20min[m20]]+gy[m20min[m20]+1])/2.0-Y)/dy;
			if (m20==0) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vz[] */
			m40=m30*xynumxy+m10*ynumy+m20min[m20];
			/* Calc wt */
			swt=ddx*ddy*ddz;
//if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Exz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
/*
*/
			/* Add EPSxz for current node */
			eps1[8]+=exz[m40]*swt;
			/* Add SIGxz for current node */
			eps1[58]+=sxz[m40]*swt;
			/* Add Vorxz for current node */
			eps1[38]+=wxz[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/*
if(ABSV(ival-1.0)>1e-6){printf("EPSxz %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[8]);getchar();}
*/
/**/
/* EPSyz Cell indexes */
m10min[0]=m1;
X=x;
if(X<(gx[m10min[0]]+gx[m10min[0]+1])/2.0) m10min[0]-=1;
if(xperiodic && (m10min[0]<0 || m10min[0]>xnumx-2))
        {
        if(m10min[0]<0) X=xsize+x;
        m10min[0]=xnumx-2;
        m10min[1]=0;
        dx=(gx[1]-gx[0]+gx[m10min[0]+1]-gx[m10min[0]])/2.0;
	}
else
        {
        if(m10min[0]<1) {m10min[0]=0; X=(gx[m10min[0]]+gx[m10min[0]+1])/2.0; }
        if(m10min[0]>xnumx-3) {m10min[0]=xnumx-3; X=(gx[m10min[0]+1]+gx[m10min[0]+2])/2.0;}
	dx=(gx[m10min[0]+2]-gx[m10min[0]])/2.0;
        m10min[1]=m10min[0]+1;
        }
m20min[0]=m2;
Y=y;
if(yperiodic)
        {
        if(m20min[0]<0) {m20min[0]=0; Y=gy[m20min[0]];}
        if(m20min[0]>ynumy-2) {m20min[0]=ynumy-2; Y=gy[m20min[0]+1];}
        }
else
        {
        if(m20min[0]<1) {m20min[0]=1; Y=gy[m20min[0]];}
        if(m20min[0]>ynumy-3) {m20min[0]=ynumy-3; Y=gy[m20min[0]+1];}
        }
m30min[0]=m3;
Z=z;
if(zperiodic)
        {
        if(m30min[0]<0) {m30min[0]=0; Z=0;}
        if(m30min[0]>znumz-2) {m30min[0]=znumz-2; Z=gz[m30min[0]+1];}
        }
else
        {
        if(m30min[0]<1) {m30min[0]=1; Z=gz[m30min[0]];}
        if(m30min[0]>znumz-3) {m30min[0]=znumz-3; Z=gz[m30min[0]+1];}
        }
/* Vy Cell dimensions */
dy= gy[m20min[0]+1]-gy[m20min[0]];
dz= gz[m30min[0]+1]-gz[m30min[0]];
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min[0];m30<=m30min[0]+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=(gz[m30]-Z)/dz;
	if (m30==m30min[0]) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=0;m10<=1;m10++)
		{
		/* Normalized X-distance calc */
		if(xperiodic && m10min[0]==xnumx-2 && m10==1) ddx=((gx[m10min[m10]]+gx[m10min[m10]+1])/2.0-(X-xsize))/dx;
		else ddx=((gx[m10min[m10]]+gx[m10min[m10]+1])/2.0-X)/dx;
		if (m10==0) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=m20min[0];m20<=m20min[0]+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=(gy[m20]-Y)/dy;
			if (m20==m20min[0]) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vy[] */
			m40=m30*xynumxy+m10min[m10]*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
//if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Eyz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPSyz for current node */
			eps1[9]+=eyz[m40]*swt;
			/* Add SIGyz for current node */
			eps1[59]+=syz[m40]*swt;
			/* Add Voryz for current node */
			eps1[39]+=wyz[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/*
if(ABSV(ival-1.0)>1e-6){printf("EPSyz %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,eps1[9]);getchar();}
*/
/**/
}
/* Calculation of EPSxx, EPSyy, EPSzz, EPSxy, EPSxz, EPSyz for current location by Linear Interpolation */


/* Left side or Value for Sxx  Equation */ 
/* Sxx=2NU(dVx/dX-1/3(dVx/dX+dVy/dY+dVz/dZ)) */
double sxxcalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
/* mgi - cur multigrid level */
{
/* Val buf */
double ival,leftsxx=0,nueff,dvxdx,dvydy,dvzdz;
long int v[8];
/* Distances */
double xkf,ykf,zkf;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
if(xperiodic && m1==mgx[mgi]-1) v[2]=mgp[mgi]+(m3-1)*mgxy[mgi]+m2-1;
if(yperiodic && m2==mgy[mgi]-1) v[1]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi];
if(zperiodic && m3==mgz[mgi]-1) v[4]=mgp[mgi]+(m1-1)*mgy[mgi]+m2-1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0)// && m1>3 && m2>3 && m3>3 && m1<mgx[mgi]-3 && m2<mgy[mgi]-3 && m3<mgz[mgi]-3)
        {
        nueff=nxx[v[7]];
        }
/* Higher levels */
else
        {
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]]+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/8.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/8.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/8.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSxx7, SIGxx7 */
/**/
/* Return Sxx,Exx val ----------------------------*/
if(ynval==0)
	{
	/* Exx=dVx/dX-1/3(dVx/dX+dVy/dY+dVz/dZ) */
	dvxdx=(vx[v[2]]-vx[v[0]])*xkf;
	dvydy=(vy[v[1]]-vy[v[0]])*ykf;
	dvzdz=(vz[v[4]]-vz[v[0]])*zkf;
	leftsxx=dvxdx-(dvxdx+dvydy+dvzdz)/3.0;
	/**/
	/* Save Exx */
	eps1[0]=leftsxx;
	/**/
	/* Calc Sxx=2Nu*Exx */
	leftsxx=2.0*nueff*leftsxx;
	/**/
	return leftsxx;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxx ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSxx7, SIGxx7 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Sxx=2NU(dVx/dX-1/3(dVx/dX+dVy/dY+dVz/dZ)) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[0]*4+1;
ui1[un1[0]+1]=-ynval*nueff*4.0/3.0*xkf;
un1[un1[0]+2]=v[2]*4+1;
ui1[un1[0]+2]=+ynval*nueff*4.0/3.0*xkf;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[0]*4+2;
ui1[un1[0]+3]=+ynval*nueff*2.0/3.0*ykf;
un1[un1[0]+4]=v[1]*4+2;
ui1[un1[0]+4]=-ynval*nueff*2.0/3.0*ykf;
/* Add Vz with koefficients */
un1[un1[0]+5]=v[0]*4+3;
ui1[un1[0]+5]=+ynval*nueff*2.0/3.0*zkf;
un1[un1[0]+6]=v[4]*4+3;
ui1[un1[0]+6]=-ynval*nueff*2.0/3.0*zkf;
/* Add total Num of lines */
un1[0]+=6;
/**/
/**/
/**/
/*
int n1;
if(m3==mgz[mgi]-1){for(n1=0;n1<=un1[0];n1++)printf("Exx %ld %ld %ld %d %e %ld \n",m1,m2,m3,n1,ui1[n1],un1[n1]);getchar();}
*/
return 0;
}
/* Left side or Value for Sxx  Equation */ 



/* Left side or Value for Syy  Equation */ 
/* Syy=2NU(dVy/dY-1/3(dVx/dX+dVy/dY+dVz/dZ)) */
double syycalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftsyy=0,nueff,dvxdx,dvydy,dvzdz;
long int v[8];
/* Distances */
double xkf,ykf,zkf;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
if(xperiodic && m1==mgx[mgi]-1) v[2]=mgp[mgi]+(m3-1)*mgxy[mgi]+m2-1;
if(yperiodic && m2==mgy[mgi]-1) v[1]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi];
if(zperiodic && m3==mgz[mgi]-1) v[4]=mgp[mgi]+(m1-1)*mgy[mgi]+m2-1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0)// && m1>3 && m2>3 && m3>3 && m1<mgx[mgi]-3 && m2<mgy[mgi]-3 && m3<mgz[mgi]-3)
        {
        nueff=nxx[v[7]];
        }
/* Higher levels */
else
        {
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]]+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/8.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/8.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/8.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSyy7, SIGyy7 */
/**/
/* Return Syy,Eyy val ----------------------------*/
if(ynval==0)
	{
	/* Eyy=dVy/dY-1/3(dVx/dX+dVy/dY+dVz/dZ) */
	dvxdx=(vx[v[2]]-vx[v[0]])*xkf;
	dvydy=(vy[v[1]]-vy[v[0]])*ykf;
	dvzdz=(vz[v[4]]-vz[v[0]])*zkf;
	leftsyy=dvydy-(dvxdx+dvydy+dvzdz)/3.0;
	/**/
	/* Save Eyy */
	eps1[0]=leftsyy;
	/**/
	/* Calc Syy=2Nu*Eyy */
	leftsyy=2.0*nueff*leftsyy;
	/**/
	return leftsyy;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Syy ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSzz7, SIGzz7 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Syy=2NU(dVy/dY-1/3(dVx/dX+dVy/dY+dVz/dZ)) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[0]*4+1;
ui1[un1[0]+1]=+ynval*nueff*2.0/3.0*xkf;
un1[un1[0]+2]=v[2]*4+1;
ui1[un1[0]+2]=-ynval*nueff*2.0/3.0*xkf;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[0]*4+2;
ui1[un1[0]+3]=-ynval*nueff*4.0/3.0*ykf;
un1[un1[0]+4]=v[1]*4+2;
ui1[un1[0]+4]=+ynval*nueff*4.0/3.0*ykf;
/* Add Vz with koefficients */
un1[un1[0]+5]=v[0]*4+3;
ui1[un1[0]+5]=+ynval*nueff*2.0/3.0*zkf;
un1[un1[0]+6]=v[4]*4+3;
ui1[un1[0]+6]=-ynval*nueff*2.0/3.0*zkf;
/* Add total Num of lines */
un1[0]+=6;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Ezz %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Syy  Equation */ 



/* Left side or Value for Szz  Equation */ 
/* Szz=2NU(dVz/dZ-1/3(dVx/dX+dVy/dY+dVz/dZ)) */
double szzcalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftszz=0,nueff,dvxdx,dvydy,dvzdz;
long int v[8];
/* Distances */
double xkf,ykf,zkf;
/**/
/**/
/**/
/* Distances Calc */
xkf=1.0/(mggx[mgi][m1]-mggx[mgi][m1-1]);
ykf=1.0/(mggy[mgi][m2]-mggy[mgi][m2-1]);
zkf=1.0/(mggz[mgi][m3]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
if(xperiodic && m1==mgx[mgi]-1) v[2]=mgp[mgi]+(m3-1)*mgxy[mgi]+m2-1;
if(yperiodic && m2==mgy[mgi]-1) v[1]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi];
if(zperiodic && m3==mgz[mgi]-1) v[4]=mgp[mgi]+(m1-1)*mgy[mgi]+m2-1;
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0)// && m1>3 && m2>3 && m3>3 && m1<mgx[mgi]-3 && m2<mgy[mgi]-3 && m3<mgz[mgi]-3)
        {
        nueff=nxx[v[7]];
        }
/* Higher levels */
else
        {
	if(viscmod==0) nueff=(nu[v[0]]+nu[v[1]]+nu[v[2]]+nu[v[3]]+nu[v[4]]+nu[v[5]]+nu[v[6]]+nu[v[7]])/8.0;
	if(viscmod==1) nueff=exp((log(nu[v[0]])+log(nu[v[1]])+log(nu[v[2]])+log(nu[v[3]])+log(nu[v[4]])+log(nu[v[5]])+log(nu[v[6]])+log(nu[v[7]]))/8.0);
	if(viscmod==2) nueff=1.0/((1.0/nu[v[0]]+1.0/nu[v[1]]+1.0/nu[v[2]]+1.0/nu[v[3]]+1.0/nu[v[4]]+1.0/nu[v[5]]+1.0/nu[v[6]]+1.0/nu[v[7]])/8.0);
	}
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSzz7, SIGzz7 */
/**/
/* Return Szz,Ezz val ----------------------------*/
if(ynval==0)
	{
	/* Ezz=dVz/dZ-1/3(dVx/dX+dVy/dY+dVz/dZ) */
	dvxdx=(vx[v[2]]-vx[v[0]])*xkf;
	dvydy=(vy[v[1]]-vy[v[0]])*ykf;
	dvzdz=(vz[v[4]]-vz[v[0]])*zkf;
	leftszz=dvzdz-(dvxdx+dvydy+dvzdz)/3.0;
	/**/
	/* Save Ezz */
	eps1[0]=leftszz;
	/**/
	/* Calc Szz=2Nu*Ezz */
	leftszz=2.0*nueff*leftszz;
	/**/
	return leftszz;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Szz ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx0 - Vx2 */
/* Vy0 | Vy1 */
/* Vz0 / Vz4 */
/* EPSzz7, SIGzz7 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Szz=2NU(dVz/dZ-1/3(dVx/dX+dVy/dY+dVz/dZ)) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[0]*4+1;
ui1[un1[0]+1]=+ynval*nueff*2.0/3.0*xkf;
un1[un1[0]+2]=v[2]*4+1;
ui1[un1[0]+2]=-ynval*nueff*2.0/3.0*xkf;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[0]*4+2;
ui1[un1[0]+3]=+ynval*nueff*2.0/3.0*ykf;
un1[un1[0]+4]=v[1]*4+2;
ui1[un1[0]+4]=-ynval*nueff*2.0/3.0*ykf;
/* Add Vz with koefficients */
un1[un1[0]+5]=v[0]*4+3;
ui1[un1[0]+5]=-ynval*nueff*4.0/3.0*zkf;
un1[un1[0]+6]=v[4]*4+3;
ui1[un1[0]+6]=+ynval*nueff*4.0/3.0*zkf;
/* Add total Num of lines */
un1[0]+=6;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Ezz %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Szz  Equation */ 





/* Left side or Value for Sxy  Equation */ 
/* Sxy=NU(dVx/dY+dVy/dX) */
double sxycalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftsxy=0,nueff,dvxdy,dvydx;
long int v[8];
/* Distances */
double xkf,ykf;
/**/
/**/
/**/
/* Distances Calc */
if(xperiodic && m1==0) xkf=1.0/(mggx[mgi][1]-mggx[mgi][0]+mggx[mgi][mgx[mgi]-1]-mggx[mgi][mgx[mgi]-2]);
else xkf=1.0/(mggx[mgi][m1+1]-mggx[mgi][m1-1]);
if(yperiodic && m2==0) ykf=1.0/(mggy[mgi][1]-mggy[mgi][0]+mggy[mgi][mgy[mgi]-1]-mggy[mgi][mgy[mgi]-2]);
else ykf=1.0/(mggy[mgi][m2+1]-mggy[mgi][m2-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+m3*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
if(xperiodic && m1==0)
	{
	v[1]=mgp[mgi]+m3*mgxy[mgi]+(mgx[mgi]-2)*mgy[mgi]+(m2);
	v[2]=mgp[mgi]+m3*mgxy[mgi]+(m2-1); v[3]=v[2]+1; v[7]=v[3]+mgxy[mgi];
	}
if(yperiodic && m2==0)
        {
        v[2]=mgp[mgi]+m3*mgxy[mgi]+(m1)*mgy[mgi]+mgy[mgi]-2;
        v[1]=mgp[mgi]+m3*mgxy[mgi]+(m1-1)*mgy[mgi]; v[3]=v[1]+mgy[mgi]; v[7]=v[3]+mgxy[mgi];
        }
if(xperiodic && m1==0 && yperiodic && m2==0)
        {
        v[1]=mgp[mgi]+m3*mgxy[mgi]+(mgx[mgi]-2)*mgy[mgi];
        v[2]=mgp[mgi]+m3*mgxy[mgi]+mgy[mgi]-2; v[3]=mgp[mgi]+m3*mgxy[mgi]; v[7]=v[3]+mgxy[mgi];
        }
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0)// && m1>3 && m2>3 && m3>2 && m1<mgx[mgi]-4 && m2<mgy[mgi]-4 && m3<mgz[mgi]-4)
        {
        nueff=nxy[v[3]];
        }
/* Higher levels */
else
        {
        if(viscmod==0) nueff=(nu[v[3]]+nu[v[7]])/2.0;
        if(viscmod==1) nueff=exp((log(nu[v[3]])+log(nu[v[7]]))/2.0);
        if(viscmod==2) nueff=1.0/((1.0/nu[v[3]]+1.0/nu[v[7]])/2.0);
        }
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx2 - Vx3 */
/* Vy1 | Vy3 */
/* EPSxy3, SIGxy3 */
/**/
/* Return Sxy,Exy val ----------------------------*/
if(ynval==0)
	{
	/* Exy=1/2(dVx/dY+dVy/dX) */
	dvxdy=(vx[v[3]]-vx[v[2]])*ykf;
	dvydx=(vy[v[3]]-vy[v[1]])*xkf;
	leftsxy=dvxdy+dvydx;
	/**/
	/* Save Exy */
	eps1[0]=leftsxy;
	eps1[1]=dvxdy-dvydx;
	/**/
	/* Calc Sxy=2Nu*Exy */
	leftsxy=2.0*nueff*leftsxy;
	/**/
	return leftsxy;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxx ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx2 - Vx3 */
/* Vy1 | Vy3 */
/* EPSxy3, SIGxy3 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Sxy=NU(dVx/dY+dVy/dX) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[2]*4+1;
ui1[un1[0]+1]=-ynval*nueff*2.0*ykf;
un1[un1[0]+2]=v[3]*4+1;
ui1[un1[0]+2]=+ynval*nueff*2.0*ykf;
/* Add Vy with koefficients */
un1[un1[0]+3]=v[1]*4+2;
ui1[un1[0]+3]=-ynval*nueff*2.0*xkf;
un1[un1[0]+4]=v[3]*4+2;
ui1[un1[0]+4]=+ynval*nueff*2.0*xkf;
/* Add total Num of lines */
un1[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exy %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxy  Equation */ 



/* Left side or Value for Sxz  Equation */ 
/* Sxz=NU(dVx/dZ+dVz/dX) */
double sxzcalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftsxz=0,nueff,dvxdz,dvzdx;
long int v[8];
/* Distances */
double xkf,zkf;
/**/
/**/
/**/
/* Distances Calc */
if(xperiodic && m1==0) xkf=1.0/(mggx[mgi][1]-mggx[mgi][0]+mggx[mgi][mgx[mgi]-1]-mggx[mgi][mgx[mgi]-2]);
else xkf=1.0/(mggx[mgi][m1+1]-mggx[mgi][m1-1]);
if(zperiodic && m3==0) zkf=1.0/(mggz[mgi][1]-mggz[mgi][0]+mggz[mgi][mgz[mgi]-1]-mggz[mgi][mgz[mgi]-2]);
else zkf=1.0/(mggz[mgi][m3+1]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+m2;v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
if(xperiodic && m1==0)
	{
	v[4]=mgp[mgi]+m3*mgxy[mgi]+(mgx[mgi]-2)*mgy[mgi]+m2;
	v[2]=mgp[mgi]+(m3-1)*mgxy[mgi]+m2; v[6]=v[2]+mgxy[mgi]; v[7]=v[6]+1;
	}
if(zperiodic && m3==0)
	{
	v[2]=mgp[mgi]+(mgz[mgi]-2)*mgxy[mgi]+(m1)*mgy[mgi]+m2;
	v[4]=mgp[mgi]+(m1-1)*mgy[mgi]+m2; v[6]=v[4]+mgy[mgi]; v[7]=v[6]+1;
	}
if(xperiodic && m1==0 && zperiodic && m3==0)
	{
	v[2]=mgp[mgi]+(mgz[mgi]-2)*mgxy[mgi]+m2;
	v[4]=mgp[mgi]+(mgx[mgi]-2)*mgy[mgi]+m2; v[6]=mgp[mgi]+m2; v[7]=v[6]+1;
	}
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0)// && m1>3 && m2>2 && m3>3 && m1<mgx[mgi]-4 && m2<mgy[mgi]-4 && m3<mgz[mgi]-4)
        {
        nueff=nxz[v[6]];
        }
/* Higher levels */
else
        {
        if(viscmod==0) nueff=(nu[v[6]]+nu[v[7]])/2.0;
        if(viscmod==1) nueff=exp((log(nu[v[6]])+log(nu[v[7]]))/2.0);
        if(viscmod==2) nueff=1.0/((1.0/nu[v[6]]+1.0/nu[v[7]])/2.0);
        }
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx2 - Vx6 */
/* Vz4 / Vz6 */
/* EPSxz6, SIGxz6 */
/**/
/* Return Sxz,Exz val ----------------------------*/
if(ynval==0)
	{
	/* Exz=1/2(dVx/dZ+dVz/dX) */
	dvxdz=(vx[v[6]]-vx[v[2]])*zkf;
	dvzdx=(vz[v[6]]-vz[v[4]])*xkf;
	leftsxz=dvxdz+dvzdx;
	/**/
	/* Save Exy */
	eps1[0]=leftsxz;
	eps1[1]=dvxdz-dvzdx;
	/**/
	/* Calc Sxy=2Nu*Exy */
	leftsxz=2.0*nueff*leftsxz;
	/**/
	return leftsxz;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Sxz ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vx2 - Vx6 */
/* Vz4 / Vz6 */
/* EPSxz6, SIGxz6 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Sxz=NU(dVx/dZ+dVz/dX) */
/* Add Vx with koefficients */
un1[un1[0]+1]=v[2]*4+1;
ui1[un1[0]+1]=-ynval*nueff*2.0*zkf;
un1[un1[0]+2]=v[6]*4+1;
ui1[un1[0]+2]=+ynval*nueff*2.0*zkf;
/* Add Vz with koefficients */
un1[un1[0]+3]=v[4]*4+3;
ui1[un1[0]+3]=-ynval*nueff*2.0*xkf;
un1[un1[0]+4]=v[6]*4+3;
ui1[un1[0]+4]=+ynval*nueff*2.0*xkf;
/* Add total Num of lines */
un1[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Exz %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Sxz  Equation */ 



/* Left side or Value for Syz  Equation */ 
/* Syz=NU(dVy/dZ+dVz/dY) */
double syzcalc(long int m1, long int m2, long int m3, double ynval, int mgi, double *eps1, long int *un1, double *ui1)
/* m1,m2,m3 - node X,Y,Z number */
/* ynval - Val Syz Calc Y(0)/N(koefficient) */
{
/* Val buf */
double ival,leftsyz=0,nueff,dvydz,dvzdy;
long int v[8];
/* Distances */
double ykf,zkf;
/**/
/**/
/**/
/* Distances Calc */
if(yperiodic && m2==0) ykf=1.0/(mggy[mgi][1]-mggy[mgi][0]+mggy[mgi][mgy[mgi]-1]-mggy[mgi][mgy[mgi]-2]);
else ykf=1.0/(mggy[mgi][m2+1]-mggy[mgi][m2-1]);
if(zperiodic && m3==0) zkf=1.0/(mggz[mgi][1]-mggz[mgi][0]+mggz[mgi][mgz[mgi]-1]-mggz[mgi][mgz[mgi]-2]);
else zkf=1.0/(mggz[mgi][m3+1]-mggz[mgi][m3-1]);
/**/
/**/
/**/
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
v[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+m1*mgy[mgi]+(m2-1);v[1]=v[0]+1;
v[2]=v[0]+mgy[mgi];v[3]=v[2]+1;
v[4]=v[0]+mgxy[mgi];v[5]=v[4]+1;
v[6]=v[4]+mgy[mgi];v[7]=v[6]+1;
if(yperiodic && m2==0)
	{
	v[1]=mgp[mgi]+(m3-1)*mgxy[mgi]+m1*mgy[mgi]; v[5]=v[1]+mgxy[mgi]; v[7]=v[5]+mgy[mgi];
	v[4]=v[5]+mgy[mgi]-2;
	//v[4]=mgp[mgi]+(m3-1)*mgxy[mgi]+m1*mgy[mgi]+mgy[mgi]-1;
	}
if(zperiodic && m3==0)
	{
	v[4]=mgp[mgi]+m1*mgy[mgi]+(m2-1); v[5]=v[4]+1; v[7]=v[5]+mgy[mgi];
	v[1]=v[5]+(mgz[mgi]-2)*mgxy[mgi];
	//v[1]=mgp[mgi]+(mgz[mgi]-2)*mgxy[mgi]+m1*mgy[mgi]+(m2);
	}
if(yperiodic && m2==0 && zperiodic && m3==0)
	{
	v[5]=mgp[mgi]+m1*mgy[mgi]; v[7]=v[5]+mgy[mgi];
	v[1]=v[5]+(mgz[mgi]-2)*mgxy[mgi];
	//v[1]=mgp[mgi]+(mgz[mgi]-2)*mgxy[mgi]+m1*mgy[mgi];
	v[4]=v[5]+mgy[mgi]-2;
	//v[4]=mgp[mgi]+m1*mgy[mgi]+mgy[mgi]-1; 
	}
/**/
/**/
/**/
/* Effective viscosity calc */
/* Basic level */
if(mgi==0)// && m1>2 && m2>3 && m3>3 && m1<mgx[mgi]-4 && m2<mgy[mgi]-4 && m3<mgz[mgi]-4)
        {
        nueff=nyz[v[5]];
        }
/* Higher levels */
else
        {
        if(viscmod==0) nueff=(nu[v[5]]+nu[v[7]])/2.0;
        if(viscmod==1) nueff=exp((log(nu[v[5]])+log(nu[v[7]]))/2.0);
        if(viscmod==2) nueff=1.0/((1.0/nu[v[5]]+1.0/nu[v[7]])/2.0);
        }
if(nueff>nuend) nueff=nuend;
if(nueff<nubeg) nueff=nubeg;
/**/
/**/
/**/
/* Nu Save */
eps1[2]=nueff;
/**/
/**/
/**/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vy1 | Vy5 */
/* Vz4 / Vz5 */
/* EPSyz5, SIGyz5 */
/**/
/* Return Syz,Eyz val ----------------------------*/
if(ynval==0)
	{
	/* Eyz=1/2(dVy/dZ+dVz/dY) */
	dvydz=(vy[v[5]]-vy[v[1]])*zkf;
	dvzdy=(vz[v[5]]-vz[v[4]])*ykf;
	leftsyz=dvydz+dvzdy;
	/**/
	/* Save Exy */
	eps1[0]=leftsyz;
	eps1[1]=dvydz-dvzdy;
	/**/
	/* Calc Sxy=2Nu*Exy */
	leftsyz=2.0*nueff*leftsyz;
	/**/
	return leftsyz;
	}
/**/
/**/
/**/
/* Add Coefficients for left parts of Syz ----------------*/
/* Staggered Nodes num */
/*   Z=0      Z=1 */
/* 0  2      4  6 */
/* 1  3      5  7 */
/* Vy1 | Vy5 */
/* Vz4 / Vz5 */
/* EPSyz5, SIGyz5 */
/**/
/*  0(P) 1(Vx)  2(Vy)  3(Vz) */
/* Syz=NU(dVy/dZ+dVz/dY) */
/* Add Vy with koefficients */
un1[un1[0]+1]=v[1]*4+2;
ui1[un1[0]+1]=-ynval*nueff*2.0*zkf;
un1[un1[0]+2]=v[5]*4+2;
ui1[un1[0]+2]=+ynval*nueff*2.0*zkf;
/* Add Vz with koefficients */
un1[un1[0]+3]=v[4]*4+3;
ui1[un1[0]+3]=-ynval*nueff*2.0*ykf;
un1[un1[0]+4]=v[5]*4+3;
ui1[un1[0]+4]=+ynval*nueff*2.0*ykf;
/* Add total Num of lines */
un1[0]+=4;
/**/
/**/
/**/
/*
for(n1=0;n1<=wn[0];n1++)printf("Eyz %e %ld \n",ui[n1],wn[n1]);getchar();
*/
return 0;
}
/* Left side or Value for Syz  Equation */ 



/* Number of nearest left X-surface find */
long int m1serch(double x)
/* x - X coordinate */
{
/* Variables */
long int m1,m10=0,m11=xnumx-1;
/**/
/* Serch cycle */
do
	{
	m1=(m10+m11)/2;
	if (gx[m1]>x) m11=m1; else m10=m1;
	}
while((m11-m10)>1);
if(m10>xnumx-2) m10=xnumx-2;
/*
if(x<gx[m10] || x>gx[m10+1]) {printf("XXX %ld %ld %ld  %e %e  %e ",m10,m11,m1,gx[m10],gx[m11],x); getchar();}
*/
return m10;
}
/* Number of nearest left X-surface find */



/* Number of nearest upper Y-surface find */
long int m2serch(double y)
/* y - Y coordinate */
{
/* Variables */
long int m2,m20=0,m21=ynumy-1;
/**/
/* Serch cycle */
do
	{
	m2=(m20+m21)/2;
	if (gy[m2]>y) m21=m2; else m20=m2;
	}
while((m21-m20)>1);
if(m20>ynumy-2) m20=ynumy-2;
/*
if(y<gy[m20] || y>gy[m20+1]) {printf("YYY %ld %ld %ld  %e %e  %e ",m20,m21,m2,gy[m20],gy[m21],y); getchar();}
*/
return m20;
}
/* Number of nearest upper Y-surface find */



/* Number of nearest frontal Z-surface find */
long int m3serch(double z)
/* z - Z coordinate */
{
/* Variables */
long int m3,m30=0,m31=znumz-1;
/**/
/* Serch cycle */
do
	{
	m3=(m30+m31)/2;
	if (gz[m3]>z) m31=m3; else m30=m3;
	}
while((m31-m30)>1);
if(m30>znumz-2) m30=znumz-2;
/*
if(z<gz[m30] || z>gz[m30+1]) {printf("ZZZ %ld %ld %ld  %e %e  %e ",m30,m31,m3,gz[m30],gz[m31],z); getchar();}
*/
return m30;
}
/* Number of nearest frontal Z-surface find */


/* Calculation of Vx,Vy,Vz for current location by Linear Interpolation */
void vxyzcalc1(double x, double y, double z, double *vxyz1, long int *wn1)
/* x,y,z - XYZ location of point for Vx,Vy,Vz calc */
{
/* Vx,Vy-nodes num */
long int m1,m2,m3,m10,m20,m30,m40,m10min[2],m20min[2],m30min[3];
/* Relativ coord */
double X,Y,Z,dx,dy,dz,ddx,ddy,ddz,swt,ival;
/**/
/**/
/**/
/* Clear Vx,Vy,Vz */
vxyz1[0]=vxyz1[1]=vxyz1[2]=0;
/* Out of grid */
if(x<0) x=0; if(x>xsize) x=xsize;
if(y<0) y=0; if(y>ysize) y=ysize;
if(z<0) z=0; if(z>zsize) z=zsize;
/* Fing Cell indexes */
wn1[0]=m1=m1serch(x);
wn1[1]=m2=m2serch(y);
wn1[2]=m3=m3serch(z);
/**/
/* Vx Cell indexes */
m10min[0]=m1;
if(m10min[0]<0) m10min[0]=0;
if(m10min[0]>xnumx-2) m10min[0]=xnumx-2;
m10min[1]=m10min[0]+1;
m20min[0]=m2;
Y=y;
if(Y<(gy[m20min[0]]+gy[m20min[0]+1])/2.0) m20min[0]-=1;
if(yperiodic && (m20min[0]<0 || m20min[0]>ynumy-3))
        {
        if(m20min[0]<0) Y=ysize+y;
        m20min[0]=ynumy-2;
        m20min[1]=0;
        dy=(gy[1]-gy[0]+gy[m20min[0]+1]-gy[m20min[0]])/2.0;
        }
else
        {
        if(m20min[0]<0) {m20min[0]=0; Y=(gy[m20min[0]]+gy[m20min[0]+1])/2.0;}
        if(m20min[0]>ynumy-3) {m20min[0]=ynumy-3; Y=(gy[m20min[0]+1]+gy[m20min[0]+2])/2.0;}
        dy=(gy[m20min[0]+2]-gy[m20min[0]])/2.0;
        m20min[1]=m20min[0]+1;
        }
m30min[0]=m3;
Z=z;
if(Z<(gz[m30min[0]]+gz[m30min[0]+1])/2.0) m30min[0]-=1;
if(zperiodic && (m30min[0]<0 || m30min[0]>znumz-3))
        {
        if(m30min[0]<0) Z=zsize+z;
        m30min[0]=znumz-2;
        m30min[1]=0;
        dz=(gz[1]-gz[0]+gz[znumz-1]-gz[znumz-2])/2.0;
	}
else
        {
        if(m30min[0]<0) {m30min[0]=0; Z=(gz[m30min[0]]+gz[m30min[0]+1])/2.0; }
        if(m30min[0]>znumz-3) {m30min[0]=znumz-3; Z=(gz[m30min[0]+1]+gz[m30min[0]+2])/2.0;}
	dz=(gz[m30min[0]+2]-gz[m30min[0]])/2.0;
	m30min[1]=m30min[0]+1;
        }
/* EPSxy Cell dimensions */
/* Vx Cell dimensions */
dx= gx[m10min[0]+1]-gx[m10min[0]];
/* Interpolate from 8 nodes */
ival=0;
for (m30=0;m30<=1;m30++)
	{
	/* Normalized Z-distance calc */
	if(zperiodic && m30min[0]==znumz-2 && m30==1) ddz=((gz[m30min[m30]]+gz[m30min[m30]+1])/2.0-(Z-zsize))/dz;
	else ddz=((gz[m30min[m30]]+gz[m30min[m30]+1])/2.0-Z)/dz;
	if (m30==0) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=m10min[0];m10<=m10min[1];m10++)
		{
		/* Normalized X-distance calc */
		ddx=(gx[m10]-x)/dx;
		if (m10==m10min[0]) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=0;m20<=1;m20++)
			{
			/* Normalized Y-distance calc */
			if(yperiodic && m20min[0]==ynumy-2 && m20==1) ddy=((gy[m20min[m20]]+gy[m20min[m20]+1])/2.0-(Y-ysize))/dy;
                	else ddy=((gy[m20min[m20]]+gy[m20min[m20]+1])/2.0-Y)/dy;
			if (m20==0) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vx[] */
			m40=m30min[m30]*xynumxy+m10*ynumy+m20min[m20];
			/* Calc wt */
			swt=ddx*ddy*ddz;
//if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
if(m40>nodenum){printf("Vx %ld %ld %ld   %ld %ld %ld %ld   %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt);getchar();}
*/
			/* Add Vx for current node */
			vxyz1[0]+=vx[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
//if(ABSV(ival-1.0)>1e-6){printf("Vx %e %e %e  %ld %ld %ld %ld  %e %e\n",x,y,z,m10min[0],m10min[1],m20min[0],m30min[0],ival,vxyz[0]);getchar();}
/*
*/
/**/
/* Vy Cell indexes */
m10min[0]=m1;
X=x;
if(X<(gx[m10min[0]]+gx[m10min[0]+1])/2.0) m10min[0]-=1;
if(xperiodic && (m10min[0]<0 || m10min[0]>xnumx-3))
	{
	if(m10min[0]<0) X=xsize+x;
	m10min[0]=xnumx-2;
	m10min[1]=0;
	dx=(gx[1]-gx[0]+gx[m10min[0]+1]-gx[m10min[0]])/2.0;
	}
else
	{
	if(m10min[0]<0) {m10min[0]=0; X=(gx[m10min[0]]+gx[m10min[0]+1])/2.0;}
	if(m10min[0]>xnumx-3) {m10min[0]=xnumx-3; X=(gx[m10min[0]+1]+gx[m10min[0]+2])/2.0;}
	dx=(gx[m10min[0]+2]-gx[m10min[0]])/2.0;
        m10min[1]=m10min[0]+1;
	}
/*
printf("%ld %ld %ld %ld %ld %e\n",m1,m2,m3,m10min[0],m10min[1],X); getchar();
*/
m20min[0]=m2;
if(m20min[0]<0) m20min[0]=0;
if(m20min[0]>ynumy-2) m20min[0]=ynumy-2;
m30min[0]=m3;
Z=z;
if(Z<(gz[m30min[0]]+gz[m30min[0]+1])/2.0) m30min[0]-=1;
if(zperiodic && (m30min[0]<0 || m30min[0]>znumz-3))
        {
        if(m30min[0]<0) Z=zsize+z;
        m30min[0]=znumz-2;
        m30min[1]=0;
        dz=(gz[1]-gz[0]+gz[znumz-1]-gz[znumz-2])/2.0;
	}
else
        {
        if(m30min[0]<0) {m30min[0]=0; Z=(gz[m30min[0]]+gz[m30min[0]+1])/2.0; }
        if(m30min[0]>znumz-3) {m30min[0]=znumz-3; Z=(gz[m30min[0]+1]+gz[m30min[0]+2])/2.0;}
	dz=(gz[m30min[0]+2]-gz[m30min[0]])/2.0;
	m30min[1]=m30min[0]+1;
        }
/* Vy Cell dimensions */
dy= gy[m20min[0]+1]-gy[m20min[0]];
/* Interpolate from 8 nodes */
ival=0;
for (m30=0;m30<=1;m30++)
	{
	/* Normalized Z-distance calc */
	if(zperiodic && m30min[0]==znumz-2 && m30==1) ddz=((gz[m30min[m30]]+gz[m30min[m30]+1])/2.0-(Z-zsize))/dz;
	else ddz=((gz[m30min[m30]]+gz[m30min[m30]+1])/2.0-Z)/dz;
	if (m30==0) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=0;m10<=1;m10++)
		{
		/* Normalized X-distance calc */
		if(xperiodic && m10min[0]==xnumx-2 && m10==1) ddx=((gx[m10min[m10]]+gx[m10min[m10]+1])/2.0-(X-xsize))/dx;
		else ddx=((gx[m10min[m10]]+gx[m10min[m10]+1])/2.0-X)/dx;
		if (m10==0) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=m20min[0];m20<=m20min[0]+1;m20++)
			{
			/* Normalized Y-distance calc */
			ddy=(gy[m20]-y)/dy;
			if (m20==m20min[0]) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vy[] */
			m40=m30min[m30]*xynumxy+m10min[m10]*ynumy+m20;
			/* Calc wt */
			swt=ddx*ddy*ddz;
//if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Vy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
/*
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
if(m40>nodenum){printf("Vy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e \n",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt);getchar();}
*/
			/* Add Vy for current node */
			vxyz1[1]+=vy[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
//if(ABSV(ival-1.0)>1e-6){printf("Vy %e %e %e   %ld %ld %ld  %e %e %e",x,y,z,m10min[0],m20min[0],m30min[0],ABSV(ival-1.0),ival,vxyz[1]);getchar();}
/*
*/
/**/
/* Vz Cell indexes */
m10min[0]=m1;
X=x;
if(X<(gx[m10min[0]]+gx[m10min[0]+1])/2.0) m10min[0]-=1;
if(xperiodic && (m10min[0]<0 || m10min[0]>xnumx-3))
	{
	if(m10min[0]<0) X=xsize+x;
	m10min[0]=xnumx-2;
	m10min[1]=0;
	dx=(gx[1]-gx[0]+gx[m10min[0]+1]-gx[m10min[0]])/2.0;
	}
else
	{
	if(m10min[0]<0) {m10min[0]=0; X=(gx[m10min[0]]+gx[m10min[0]+1])/2.0;}
	if(m10min[0]>xnumx-3) {m10min[0]=xnumx-3; X=(gx[m10min[0]+1]+gx[m10min[0]+2])/2.0;}
	dz=(gz[m30min[0]+2]-gz[m30min[0]])/2.0;
        m10min[1]=m10min[0]+1;
	}
m20min[0]=m2;
Y=y;
if(Y<(gy[m20min[0]]+gy[m20min[0]+1])/2.0) m20min[0]-=1;
if(yperiodic && (m20min[0]<0 || m20min[0]>ynumy-3))
        {
        if(m20min[0]<0) Y=ysize+y;
        m20min[0]=ynumy-2;
        m20min[1]=0;
        dy=(gy[1]-gy[0]+gy[m20min[0]+1]-gy[m20min[0]])/2.0;
	}
else
        {
        if(m20min[0]<0) {m20min[0]=0; Y=(gy[m20min[0]]+gy[m20min[0]+1])/2.0;}
        if(m20min[0]>ynumy-3) {m20min[0]=ynumy-3; Y=(gy[m20min[0]+1]+gy[m20min[0]+2])/2.0;}
        dy=(gy[m20min[0]+2]-gy[m20min[0]])/2.0;
        m20min[1]=m20min[0]+1;
        }
m30min[0]=m3;
if(m30min[0]<0) m30min[0]=0;
if(m30min[0]>znumz-2) m30min[0]=znumz-2;
/* Vz Cell dimensions */
if(xperiodic && m10min[0]==xnumx-2) dx=(gx[1]+gx[m10min[0]+1]-gx[m10min[0]])/2.0;
else dx=(gx[m10min[0]+2]-gx[m10min[0]])/2.0;
dz= gz[m30min[0]+1]-gz[m30min[0]];
/* Interpolate from 8 nodes */
ival=0;
for (m30=m30min[0];m30<=m30min[0]+1;m30++)
	{
	/* Normalized Z-distance calc */
	ddz=(gz[m30]-z)/dz;
	if (m30==m30min[0]) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=0;m10<=1;m10++)
		{
		/* Normalized X-distance calc */
		if(xperiodic && m10min[0]==xnumx-2 && m10==1) ddx=((gx[m10min[m10]]+gx[m10min[m10]+1])/2.0-(X-xsize))/dx;
		else ddx=((gx[m10min[m10]]+gx[m10min[m10]+1])/2.0-X)/dx;
		if (m10==0) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=0;m20<=1;m20++)
			{
			/* Normalized Y-distance calc */
			if(yperiodic && m20min[0]==ynumy-2 && m20==1) ddy=((gy[m20min[m20]]+gy[m20min[m20]+1])/2.0-(Y-ysize))/dy;
                	else ddy=((gy[m20min[m20]]+gy[m20min[m20]+1])/2.0-Y)/dy;
			if (m20==0) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vz[] */
			m40=m30*xynumxy+m10min[m10]*ynumy+m20min[m20];
			/* Calc wt */
			swt=ddx*ddy*ddz;
//if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
/*
if(m40>nodenum){printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e  %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt);getchar();}
*/
			/* Add Vx for current node */
			vxyz1[2]+=vz[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
//if(ABSV(ival-1.0)>1e-6){printf("Vz %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min[0],m20min[0],m30min[0],ival,vxyz[2]);getchar();}
/*
*/
/**/
}
/* Calculation of Vx,Vy,Vz for current location by Linear Interpolation */
