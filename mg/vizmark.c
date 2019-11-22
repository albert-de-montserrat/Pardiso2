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
int n0,n1,f0,f1,mtype;
long int pos0cur0,m0,m1,m2,m3,m4;
char c[3],sb[50],sc[50],sd[50],marktypetemp[100];
char fl0in[MAXFLN][60];
/**/
/**/
/* Open File amir.t3c */
fl = fopen("vizmark.t3c","rt");
/**/
/* Save Input file name */
ffscanf();
for (m1=0;m1<50;m1++) sd[m1]=sc[m1]=sb[m1]=sa[m1];
/**/
/* Load marker type to visualize */
ffscanf();
mtype=0;
while(sa[0]!='~')
        {
	marktypetemp[mtype]=atoi(sa);
        /**/
	mtype++;
        /**/
	ffscanf();
        /**/
	}
/**/
/* Load, Save nodes/marker infos */
ffscanf(); fl0stp2[0]=atoi(sa); /* Data Compression */
ffscanf(); fl0stp2[1]=atoi(sa); /* Resolution */
ffscanf(); fl0stp2[2]=atof(sa);/* X init */
ffscanf(); fl0stp2[3]=atof(sa);/* X end */
ffscanf(); fl0stp2[4]=atof(sa);/* Y init */
ffscanf(); fl0stp2[5]=atof(sa);/* Y end */
ffscanf(); fl0stp2[6]=atof(sa);/* Z init */
ffscanf(); fl0stp2[7]=atof(sa);/* Z end */
/**/
/* Load Input File number */
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
	strcat(sc,"composition.h5");
        for (n1=0;n1<50;n1++) fl0out[fl0num][n1]=sc[n1];
	/**/
	fl0num++;
	/**/
        ffscanf();
	/**/
	}
/**/  
/* Close file */
fclose(fl);
printf("\nVisualize lithologies:\n");
for(m1=1;m1<mtype;m1++)
	{
	/* All lithologies */
	if(marktypetemp[0]) marktypetemp[m1]=1;
	marktype[m1-1]=marktypetemp[m1];
	if(marktype[m1-1]) printf(" %ld",m1-1);
	}
/* Output File Cycle */
for (f0=0;f0<fl0num;f0++)
{       
/* Reload Cur Input File Name, Type */
for (n1=0;n1<50;n1++) fl1in[n1]=fl0in[f0][n1];
fl1itp=2;
/**/ 
/* Reload Cur Output File Name, Type */
for (n1=0;n1<50;n1++) fl1out[n1]=fl0out[f0][n1];
fl1otp=fl0otp[f0];
/* Load fr m data file */
/*
loader();
*/
printf("\nComposition of file %s: \n",fl1in);
makeXMFviz(f0);
printf("Result printed to file %s \n",fl1out);
printf("\nFile size is approximately %f MB \n",(float)(ncell*sizeof(long int)*8+nvrtx*(sizeof(char)+sizeof(float)*3))/1.0e+6);
printf("File size on Paraview is approximately %f MB \n",(float)(ncell*sizeof(int)*8+nvrtx*(sizeof(int)+sizeof(float)*2*3))/1.0e+6);
}
/**/
}


/* Make XMF file to visualize .h5 fields in Paraview -------------*/
void makeXMFviz(int f0)
{
int m1,m2,m3,m4;
long int i,j,k,ijk,ijk2;
char fname[50];
char c[3];
/**/
hid_t       file_id, file_id_0, file_id_1, group_id, attr_id, dataset_id, dataspace_id, dcpl, memspace;
herr_t      status;
hsize_t     dims,ndims[2],count[2],offset[2];
/**/
/* Load infos from input file */
loadinfos();
/**/
/* Create output hdf5 file */
file_id_1 = H5Fcreate(fl1out,H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
/**/
/**/
/* Form xmf file name */
if(num[f0]<10)
{
strcpy(fname,"XDMF.composition.00");
}
if(num[f0]>=10 && num[f0]<100)
{
strcpy(fname,"XDMF.composition.0");
}
if(num[f0]>=100 && num[f0]<1000)
{
strcpy(fname,"XDMF.composition.");
}
sprintf(c,"%d",num[f0]);
strcat(fname,c); 
strcat(fname,".xmf");
/* Open xmf file */  
fl = fopen(fname,"wt");
/**/
/**/
/* Check voronoi grid */
gridcell(f0);
/**/
/**/
/* Dynamically allocate memory */
dynmemall(4);
/**/
/**/
/* Calculate Voronoi tessellation */
voronoi();
/**/
/**/
/* Compress data */
if(compress)
{
merge_voronoicells(&ncell,&nvrtx);
}
/**/
/**/
/* 3D plot */
/**/
/**/
/* Cell vertexes */
if(!compress){
for(k=0;k<znum;k++)
for(i=0;i<xnum;i++)
for(j=0;j<ynum;j++)
	{
	ijk=j+i*ynum+k*xy;
	vrtx[ijk][0]=(float)(gxc[i]-xstp/2);
	vrtx[ijk][1]=(float)(ysize-(gyc[j]-ystp/2));
	vrtx[ijk][2]=(float)(gzc[k]-zstp/2);
	/*
	printf("%ld %ld %ld %ld %e %e %e\n",i,j,k,ijk,vrtx[ijk][0],vrtx[ijk][1],vrtx[ijk][2]); getchar();
	*/
	}
}
ndims[0]=nvrtx;
ndims[1]=3;
dataspace_id=H5Screate_simple(2,ndims,NULL);
dataset_id = H5Dcreate(file_id_1, "/vertices", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, vrtx);
H5Dclose(dataset_id);
H5Sclose(dataspace_id);
/**/
/**/
/**/
/* Cell connectivity */
if(!compress){
for(k=0;k<zcell;k++)
for(i=0;i<xcell;i++)
for(j=0;j<ycell;j++)
	{
	ijk=j+i*ycell+k*xyc;
	cntv[ijk][0]=j+i*ynum+k*xy;
	cntv[ijk][1]=cntv[ijk][0]+1;
	cntv[ijk][2]=cntv[ijk][0]+ynum+1;
	cntv[ijk][3]=cntv[ijk][0]+ynum;
	cntv[ijk][4]=cntv[ijk][0]+xy;
	cntv[ijk][5]=cntv[ijk][0]+xy+1;
	cntv[ijk][6]=cntv[ijk][0]+xy+ynum+1;
	cntv[ijk][7]=cntv[ijk][0]+xy+ynum;
	/*
	printf("%ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",i,j,k,ijk,cntv[ijk][0],cntv[ijk][1],cntv[ijk][2],cntv[ijk][3],cntv[ijk][4],cntv[ijk][5],cntv[ijk][6],cntv[ijk][7]); getchar();
	*/
	}
}
ndims[0]=ncell;
ndims[1]=8;
dataspace_id=H5Screate_simple(2,ndims,NULL);
dataset_id = H5Dcreate(file_id_1, "/connectivity", H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, cntv);
H5Dclose(dataset_id);
H5Sclose(dataspace_id);
/**/
/**/
/**/
/* Cell index */
if(!compress){
for(k=0;k<zcell;k++)
for(i=0;i<xcell;i++)
for(j=0;j<ycell;j++)
	{
	ijk=j+i*ycell+k*xyc;
	ijk2=j+i*ynum+k*xy;
	Idxv2[ijk2]=idx[ijk];
	}
/* Assign index to external vertices */
for(k=0;k<znum;k++)
for(i=0;i<xnum;i++)
for(j=0;j<ynum;j++)
	{
	if(k==znum-1 || j==ynum-1 || i==xnum-1)
		{
		ijk2=j+i*ynum+k*xy;
		if(k==znum-1) Idxv2[ijk2]=Idxv2[ijk2-xy];
		if(j==ynum-1) Idxv2[ijk2]=Idxv2[ijk2-1];
		if(i==xnum-1) Idxv2[ijk2]=Idxv2[ijk2-ynum];
		}
	/*
	printf("%ld %d\n",ijk2,Idxv2[ijk2]);
	*/
	}	
dims=nvrtx;
loadsave_dataset_char(file_id_1,dims,Idxv2, "/index",0);
}
/* Compressed data */
else 
{
dims=nvrtx;
loadsave_dataset_char(file_id_1,dims,Idxv2, "/index",0);
}
/**/
/**/
/* Close hdf5 file */
H5Fclose(file_id_1);
/**/
/**/
/* Free memory */
dynmemall(3);
dynmemall(5);
/**/
/**/
/* Make XDMF file */
printf("\nGenerate file %s\n",fname);
printf("\n");
/**/
fprintf(fl,"<?xml version=\"1.0\" ?> \n");
fprintf(fl,"<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n\n");
fprintf(fl," <Domain>\n\n");
/**/
fprintf(fl,"    <Grid Name=\"Eulerian Grid\">\n\n");
fprintf(fl,"       <Time Value=\"%f\" />\n\n",timesum);
fprintf(fl,"         <Topology Type=\"Hexahedron\" NumberOfElements=\"%ld\"> \n",ncell);
fprintf(fl,"            <DataItem Format=\"HDF\" DataType=\"Int\"  Dimensions=\"%ld 8\">%s:/connectivity</DataItem> \n",ncell,fl1out);
fprintf(fl,"         </Topology> \n\n");
fprintf(fl,"         <Geometry Type=\"XYZ\"> \n");
fprintf(fl,"            <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">%s:/vertices</DataItem> \n",nvrtx,fl1out);
fprintf(fl,"         </Geometry>\n\n");
			/**/
			/**/
			/* Add index array */
			printf("%ld\n",nvrtx);
			add_eulerian_to_XDMF("Index",fl1out,"/index",0,nvrtx,fl);
			/**/
			/**/
fprintf(fl,"    </Grid>\n\n");
fprintf(fl," </Domain>\n");
fprintf(fl,"</Xdmf>\n");
/**/
/**/
/* Close xmf file */
fclose(fl);
/**/
}
/* End Make XMF file to visualize .h5 fields in Paraview -------------*/


/* Make grid for voronoi cells */
void gridcell(int f0)
{
int i,resol;
double xinit,xend,yinit,yend,zinit,zend;
/* Make Eulerian grid finer than Lagrangian grid */
compress=fl0stp2[0];
resol=fl0stp2[1];
xinit=fl0stp2[2];
xend=fl0stp2[3];
yinit=fl0stp2[4];
yend=fl0stp2[5];
zinit=fl0stp2[6];
zend=fl0stp2[7];
m10x=m1serch(xsize*xinit);
m20x=m1serch(xsize*xend)+1;
m10y=m2serch(ysize*yinit);
m20y=m2serch(ysize*yend)+1;
m10z=m3serch(zsize*zinit);
m20z=m3serch(zsize*zend)+1;
printf("\n");
printf("from x0=%e 	to x1=%e \n",gx[m10x],gx[m20x]);
printf("from y0=%e 	to y1=%e \n",gy[m10y],gy[m20y]);
printf("from z0=%e 	to z1=%e \n",gz[m10z],gz[m20z]);
printf("\n");
xcell=resol*(m20x-m10x)*mnumx; 
ycell=resol*(m20y-m10y)*mnumy;
zcell=resol*(m20z-m10z)*mnumz; 
xnum=xcell+1;
ynum=ycell+1;
znum=zcell+1;
ncell=xcell*ycell*zcell;
nvrtx=xnum*ynum*znum;
xstp=xsize*(xend-xinit)/xcell;
ystp=ysize*(yend-yinit)/ycell;
zstp=zsize*(zend-zinit)/zcell;
xyc=xcell*ycell;
xy=xnum*ynum;
printf("Resolution is %d\n",resol);
printf("\n");
printf("Number of markers=%ld Number of cells=%ld Number of vertices=%ld\n",marknum,ncell,nvrtx);
printf("\n");
/**/
}

/* Load infos */
void loadinfos()
{
hid_t       file_id, file_id_0, file_id_1, group_id, attr_id, dataset_id, dataspace_id, dcpl, memspace;
herr_t      status;
hsize_t     dims,ndims[2],count[2],offset[2];
double a[100];
long int b[4];
/**/
/* Open an existing file */
file_id_0 = H5Fopen(fl1in, H5F_ACC_RDONLY, H5P_DEFAULT);
/* General */
attr_id=H5Aopen_by_name(file_id_0,".","General",H5P_DEFAULT,H5P_DEFAULT);
H5Aread (attr_id, H5T_NATIVE_DOUBLE, &a);
H5Aclose(attr_id);
/**/
xsize=a[0];
ysize=a[1];
zsize=a[2];
timesum=a[10]; timesum*=3.15576e+7;
/* Read Nodes */
group_id = H5Gopen(file_id_0, "/Nodes", H5P_DEFAULT);
/* Read /Nodes attribute */
attr_id=H5Aopen_by_name(group_id,".","/Nodes/Number",H5P_DEFAULT,H5P_DEFAULT);
H5Aread (attr_id, H5T_NATIVE_LONG, &b);
H5Aclose(attr_id);
xnumx=b[0];
ynumy=b[1];
znumz=b[2];
/* Gridlines positions */
gx=malloc(sizeof(double) * xnumx);
gy=malloc(sizeof(double) * ynumy);
gz=malloc(sizeof(double) * znumz);
/* Load gridline positions */
group_id = H5Gopen(file_id_0, "/Nodes", H5P_DEFAULT);
dims=xnumx;
loadsave_dataset_double(group_id,dims,gx, "/Nodes/gx",1);
dims=ynumy;
loadsave_dataset_double(group_id,dims,gy, "/Nodes/gy",1);
dims=znumz;
loadsave_dataset_double(group_id,dims,gz, "/Nodes/gz",1);
H5Gclose(group_id);
/* Lagrangian particles */
group_id = H5Gopen(file_id_0, "/Particles", H5P_DEFAULT);
/* Read /Nodes attribute */
attr_id=H5Aopen_by_name(group_id,".","/Particles/Number",H5P_DEFAULT,H5P_DEFAULT);
H5Aread (attr_id, H5T_NATIVE_LONG, &b);
H5Aclose(attr_id);
marknum=b[0];
mnumx=b[1];
mnumy=b[2];
mnumz=b[3];
/*
printf("%ld %ld %ld %ld\n",mnumx,mnumy,mnumz,marknum);
 */
/**/
dynmemall(2);
/* Particle parameters */
dims=marknum;
/* Load marker coordinates from 3 * marknum matrix */
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
/* Load Marker type */
loadsave_dataset_char(file_id_0,dims, markt,"/Particles/markt",1);
/**/
H5Gclose(group_id);
}

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


/* Voronoi tessellation */
void voronoi()
{
int n,ln,tid;
int printmod=0;
long int i,j,k,l,m,N,xn,yn,zn;
double dx,dy,dz,dxyz1,dxyz2;
/**/
printf("Calculate voronoi tessellation:\n");
/* Cell centre coordinates */
for(i=0;i<xnum;i++)
	{
	gxc[i]=gx[m10x]+((double)(i)+0.5)*xstp;
	}
for(i=0;i<ynum;i++)
	{
	gyc[i]=gy[m10y]+((double)(i)+0.5)*ystp;
	}
for(i=0;i<znum;i++)
	{
	gzc[i]=gz[m10z]+((double)(i)+0.5)*zstp;
	}
/**/
/* Set initial cell index to -1 */
for(i=0;i<ncell;i++) idx[i]=-1;
/**/
l=m=0;
/* Assign a cell to each marker Each cell will have as index the marker index */
/* If more than a marker is inside the cell, choose the one closer to the centre */
for(i=0;i<marknum;i++)
	{
	if(markx[i]>gx[m10x] && markx[i]<gx[m20x] && marky[i]>gy[m10y] && marky[i]<gy[m20y] && markz[i]>gz[m10z] && markz[i]<gz[m20z])
	{
	/*
	printf("%e %e %e %e %e %e %e %e %e\n",markx[i],gx[m10x],gx[m20x],marky[i],gy[m10y],gy[m20y],markz[i],gz[m10z],gz[m20z]);getchar();
	printf("%ld %d\n",Lc[i][0],idx[Lc[i][0]]); getchar();
	*/
	/* Find nearest cell */
	xn=floor(((double)(markx[i])-gx[m10x])/xstp);
	yn=floor(((double)(marky[i])-gy[m10y])/ystp);
	zn=floor(((double)(markz[i])-gz[m10z])/zstp);
	/* List of claimed cells */
	Lc[m]=yn+xn*ycell+zn*xyc;
	/* Cell index */
	if(idx[Lc[m]]==-1) 
		{
		idx[Lc[m]]=i;
		/* Number of claimed cells */
		m++;
		}
	/* More than 1 marker in the same cell, choose only the one closer to the centroid */
	else
		{
	        /* Check for distance */
	        dx=markx[i]-gxc[xn];
	        dy=marky[i]-gyc[yn];
	        dz=markz[i]-gzc[zn];
	        dxyz1=pow(dx*dx+dy*dy+dz*dz,0.5);
	        k=idx[Lc[m]];
	        dx=markx[k]-gxc[xn];
	        dy=marky[k]-gyc[yn];
	        dz=markz[k]-gzc[zn];
	        dxyz2=pow(dx*dx+dy*dy+dz*dz,0.5);
		if(dxyz1<dxyz2) 
			{
			idx[Lc[m]]=i;
			}
		}
	/*
	printf("%ld %e %e %ld %e %e %ld %e %e %ld %ld\n",i,markx[i],xstp,xn,marky[i],ystp,yn,markz[i],zstp,zn,Lc[i][0]); 
	getchar();
	*/
	}
	}
/*
*/
printf("Initial seeds= %ld \n",m);
printf("Claimed cells= %ld Unclaimed cells= %ld\n",m,ncell-m);
/**/
l=m;
N=0;
while (m<ncell)
	{
	/* All nulls for lt */
	for(i=0;i<N;i++) lt[i] = 0;
	/* N = number of newly claimed cells */
        N=0;
	for(i=0;i<l;i++)
		{
		/* Find coordinates of cell */
		cellsite(Lc[i]);
		/* Find neighbours cells and their coordinates, output is # of neighobours */
		ln=cellneighbours(Lc[i]);
		for(n=0;n<ln;n++) 
			{
			/*
			if(i==494) printf("%ld %ld %ld %ld %ld %ld\n",i,Lc[i][j],Ln[n],xlayers,ylayers,zlayers);
			*/
			if(idx[Ln[n]]==-1)
				{
				/* Ln is the list of neighnours cells */
				/* lt is the list of new claimed cells */
				/* Insert new claimed cell into lt */
				lt[N]=Ln[n];
				/* Assign marker index to the newly claimed cells */
				idx[Ln[n]]=idx[Lc[i]];
				/* Increase counter of newly claimed cells */
				N++;
				/* Increase total number of claimed cells */
				m++;
				/*
				if(idx[Ln[n]]==-1) printf("%ld %ld \n",m,Ln[n]);
				*/
				} 
			/* 
			else if(idx[Ln[n]]!=idx[Lc[i][j]])
			*/
			else if(idx[Ln[n]]!=idx[Lc[i]])
			        {
			        /* Battle !!! */
			        /* Check for distance */
			        dx=markx[i]-gxcn[n];
			        dy=marky[i]-gycn[n];
			        dz=markz[i]-gzcn[n];
			        dxyz1=pow(dx*dx+dy*dy+dz*dz,0.5);
			        k=idx[Ln[n]];
			        dx=markx[k]-gxcn[n];
			        dy=marky[k]-gycn[n];
			        dz=markz[k]-gzcn[n];
			        dxyz2=pow(dx*dx+dy*dy+dz*dz,0.5);
			        /* Claim cell */
			        if(dxyz1<dxyz2)
					{
					/* Insert new claimed cell into lt */
					lt[N]=Ln[n];
					/* Assign marker index to the newly claimed cells */
					idx[Ln[n]]=idx[Lc[i]];
					/* Increase counter of newly claimed cells */
					N++;
					}
				}
			/*
                        if(Ln[n]) printf("%ld %d\n",i,idx[Ln[n]]); getchar();
			*/
			}
		}
       	/**/
	if(N>marknum*6){printf("Space out in Lc: %ld is bigger than %ld\n",N,marknum*6); exit(0);}
	/* Loop over number of claimed cells in the previous step */
	l=N;
	/* Assign boundary cells */
       	for(i=0;i<N;i++)
       		{
		Lc[i]=lt[i];
       		/*
		printf("%ld %ld %ld\n",i,j,Lc[i]); getchar();
		*/
		}
       	if(printmod) printf("\nOK %ld %ld\n",i,N);
	/**/
	/**/
	printf("Claimed cells= %ld Unclaimed cells= %ld\n",m,ncell-m);
	}
/**/
/**/
/* Assign marker type */
for(i=0;i<ncell;i++)
	{
	idx[i]=markt[idx[i]];
	/*
	if(idx[i]==) {printf("%ld \n",i); getchar();}
	*/
	}
}
/* Voronoi tessellation */

/* Compress data = merge voronoi cells */
void merge_voronoicells(long int *Numcell,long int *Numvrtx)
{
long int i,j,k,m,n,p,ijk,ijk1,ix;
int x,y,z,xf,yf,zf,xi,yi,zi;
int x0,y0,z0,x1,y1,z1;
/**/
printf("\nCompress data: \n");
for(ijk=0;ijk<ncell;ijk++) Idxc[ijk]=idx[ijk];
/**/
long int numcell=0,numvrtx=0,cntvtemp[8];
for(ijk=0;ijk<ncell;ijk++)
	{
	/* Find next cell not merged yet and with lithology flag != 0*/
	while(Idxc[ijk]==-1 || marktype[idx[ijk]]==0) {ijk++; if(ijk>=ncell) break;}
	if(ijk>=ncell) break;
	/*
	printf("%ld %d\n",ijk,idx[ijk]);
	*/
	/* Find cell position */
	cellsite(ijk);
	x0=xlayers;
	y0=ylayers;
	z0=zlayers;
	x1=1;
	y1=1;
	z1=1;
	/* Initialize x,y,z flags */
	xf=yf=zf=1;
	/* Merge vornoi cells for cube ijk */
	while(xf || yf || zf)
		{
		/* Increase one cell in the x direction */
		if(xf)
			{
			/* Check all the cells on the right */
			for(m=0;m<y1;m++)
			for(n=0;n<z1;n++)
				{
				/* idx instead of Idx : check for the largest cell */
				if(idx[ijk+x1*ycell+m*1+n*xyc]!=idx[ijk]) {xf=0; break;}
				/*
				if(Idxc[ijk+x1*ycell+m*1+n*xyc]!=Idxc[ijk]) {xf=0; break;}
				*/
				}
			if(xf) 
				{
				if(x0+x1+1==xnum) xf=0;
				else x1++;
				}
			}
		/* Increase one cell in the y direction */
		if(yf)
			{
			/* Check all the cells below */
			for(m=0;m<x1;m++)
			for(n=0;n<z1;n++)
				{
				/* idx instead of Idx : check for the largest cell */
				if(idx[ijk+y1*1+m*ycell+n*xyc]!=idx[ijk]) {yf=0; break;}
				/*
				if(Idxc[ijk+y1*1+m*ycell+n*xyc]!=Idxc[ijk]) {yf=0; break;}
				*/
				}
			if(yf) 
				{
				if(y0+y1+1==ynum) yf=0;
				else y1++;
				}
			}
		/* Increase one cell in the z direction */
		if(zf)
			{
			/* Check all the cells in the front */
			for(m=0;m<x1;m++)
			for(n=0;n<y1;n++)
				{
				/* idx instead of Idx : check for the largest cell */
				if(idx[ijk+n*1+m*ycell+z1*xyc]!=idx[ijk]) {zf=0; break;}
				/*
				if(Idxc[ijk+n*1+m*ycell+z1*xyc]!=Idxc[ijk]) {zf=0; break;}
				*/
				}
			if(zf)
				{
				if(z0+z1+1==znum) zf=0;
				else z1++;
				}
			}
		}
	/* Reset all null cells */
	for(m=0;m<x1;m++)
	for(n=0;n<y1;n++)
	for(p=0;p<z1;p++)
		{
		Idxc[ijk+m*ycell+n*1+p*xyc]=-1;
		}
	/*
	printf("ijk=%ld %d %d %d %d %d %d %d\n",ijk,x0,x1,y0,y1,z0,z1,idx[ijk]); getchar();
	*/
	/* Save cell verteces connectivity, position and index */ 
	for(m=0;m<8;m++) 
		{
		if(m==0) {x=y=z=0;}
		if(m==1) {x=x1; y=z=0;}
		if(m==2) {x=z=0; y=y1;}
		if(m==3) {x=x1; y=y1; z=0;}   
		if(m==4) {x=y=0; z=z1;}
		if(m==5) {x=x1; y=0; z=z1;}    
		if(m==6) {x=0; y=y1; z=z1;}
		if(m==7) {x=x1; y=y1; z=z1;}   
		/* Vertex index */
		ijk1=(y0+y)+(x0+x)*ynum+(z+z0)*xy;
		/* New vertex if not exisiting or it exists but has a different marker type */
		if(Idxv[ijk1]==0 || Idxv[ijk1]>0 && (Idxv2[Idxv[ijk1]]!=idx[ijk]))
			{
			/* Save connectivity */
			cntvtemp[m]=numvrtx;
			/* Save vertices */
			vrtx[numvrtx][0]=(float)(gxc[x0+x]-xstp/2);
			vrtx[numvrtx][1]=(float)(ysize-(gyc[y0+y]-ystp/2));
			vrtx[numvrtx][2]=(float)(gzc[z0+z]-zstp/2);
			/* Save index */
			Idxv[ijk1]=numvrtx;
			/* Assign marker type to the vertex = marker type of the cell */
			Idxv2[numvrtx]=idx[ijk];
			/*
			printf("%d %ld %ld %ld %d %d %d %d %e %e %e\n",Idxv2[numvrtx],ijk1,m,numvrtx,x0,x,y,z,vrtx[numvrtx][0],vrtx[numvrtx][1],vrtx[numvrtx][2]);
			*/
			/* Increase number of verteces */
			numvrtx++;
			}
		/* Existing vertex with same marker type */
		else
			{
			cntvtemp[m]=Idxv[ijk1];
			/*
			Idxv2[cntvtemp[m]]=idx[ijk];
			printf("%d %ld %ld %ld %d %d %d %e %e %e\n",Idxv2[cntvtemp[m]],ijk1,m,cntvtemp[m],x,y,z,vrtx[cntvtemp[m]][0],vrtx[cntvtemp[m]][1],vrtx[cntvtemp[m]][2]);
			*/
			}
		}
		cntv[numcell][0]=cntvtemp[0];
		cntv[numcell][1]=cntvtemp[1];
		cntv[numcell][2]=cntvtemp[3];
		cntv[numcell][3]=cntvtemp[2];
		cntv[numcell][4]=cntvtemp[4];
		cntv[numcell][5]=cntvtemp[5];
		cntv[numcell][6]=cntvtemp[7];
		cntv[numcell][7]=cntvtemp[6];
	/*
	printf(" %ld %ld %ld %ld %ld %ld %ld %ld %ld\n",ijk,cntv[numcell][0],cntv[numcell][1],cntv[numcell][2],cntv[numcell][3],cntv[numcell][4],cntv[numcell][5],cntv[numcell][6],cntv[numcell][7]); getchar();
	*/
	/* Increase cell number */
	numcell++;
	}
printf("Number of final cells=%ld instead of=%ld\n",numcell,ijk);
printf("Number of final vertexes=%ld instead of=%ld\n",numvrtx,nvrtx);
/**/
*Numcell=numcell;
*Numvrtx=numvrtx;
}


/* Find cell site */
void cellsite(long int lc)
{
int leftcell;
lc+=1;
/* Determine position of cell */
zlayers=floor(lc/xyc);
/* Number of left cells (add 1 cell) */
leftcell=lc-zlayers*xyc;
/* Last cell of layer */
if(leftcell==0)
       	{
	ylayers=ycell;
       	xlayers=xcell;
       	}
else
       	{
	zlayers++;
       	xlayers=floor(leftcell/ycell);
        ylayers=leftcell-xlayers*ycell;
        /* Last cell of row */
        if(ylayers==0) ylayers=ycell;
	else xlayers++;
        }
xlayers--;
ylayers--;
zlayers--;
/*
if(zlayers>0){printf("%ld %ld %ld %ld %ld %ld\n",i,j,Lc[i][j],xlayers,ylayers,zlayers); getchar();}
printf("%e  %e  %e  \n",gxc[xlayers],gyc[ylayers],gzc[zlayers]);getchar();
*/
}

/* Cell neighbours */
int cellneighbours(long int lc)
{
int a=0;
/* Left cell */
if(xlayers>0 && idx[lc-ycell]!=idx[lc])
        {
	Ln[a]=lc-ycell;
	gxcn[a]=gxc[xlayers-1];
	gycn[a]=gyc[ylayers];
	gzcn[a]=gzc[zlayers];
	a++;
       	}
/* Right cell */
if(xlayers<xcell-1 && idx[lc+ycell]!=idx[lc])
        {
        Ln[a]=lc+ycell;
        gxcn[a]=gxc[xlayers+1];
        gycn[a]=gyc[ylayers];
        gzcn[a]=gzc[zlayers];
        a++;
        }
/* Back cell */
if(zlayers>0 && idx[lc-xyc]!=idx[lc])
        {
        Ln[a]=lc-xyc;
        gxcn[a]=gxc[xlayers];
        gycn[a]=gyc[ylayers];
        gzcn[a]=gzc[zlayers-1];
        a++;
        }
/* Front cell */
if(zlayers<zcell-1 && idx[lc+xyc]!=idx[lc])
        {
        Ln[a]=lc+xyc;
        gxcn[a]=gxc[xlayers];
        gycn[a]=gyc[ylayers];
        gzcn[a]=gzc[zlayers+1];
        a++;
        }
/* Up cell */
if(ylayers>0 && idx[lc-1]!=idx[lc])
        {
        Ln[a]=lc-1;
        gxcn[a]=gxc[xlayers];
        gycn[a]=gyc[ylayers-1];
        gzcn[a]=gzc[zlayers];
        a++;
        }
/* Down cell */
if(ylayers<ycell-1 && idx[lc+1]!=idx[lc])
        {               
        Ln[a]=lc+1;
        gxcn[a]=gxc[xlayers];
        gycn[a]=gyc[ylayers+1];
        gzcn[a]=gzc[zlayers];
        a++;
        }
return a;
}
