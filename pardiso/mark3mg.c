/* Calculate inclusion/matrix average strain rate, deviatoric stress tensor components */
void effbulkvisccalc(int nodemark)
{
if(nodemark==1)
{
long int mm1,m1,m2,m3,m4,v1[8];
long int numtot=0,numinc=0,nummat=0;
double tau2=0,eiinloc,siinloc,hsnloc,eiivploc;
double Exy,Exz,Eyz,Sxy,Sxz,Syz,Nuxy,Nuxz,Nuyz;
double nubulk1,nubulk2,nubulk3,nubulk4,nubulkn;
double siin,eiin,eiivp,hsn;
double exxn,eyyn,ezzn,exyn,exzn,eyzn;
double sxxn,syyn,szzn,sxyn,sxzn,syzn;
double exxi,eyyi,ezzi,exyi,exzi,eyzi;
double sxxi,syyi,szzi,sxyi,sxzi,syzi;
//Set averages to zeros 
siin=eiin=hsn=eiivp=0;
exxn=eyyn=ezzn=exyn=exzn=eyzn=0;
sxxn=syyn=szzn=sxyn=sxzn=syzn=0;
exxi=eyyi=ezzi=exyi=exzi=eyzi=0;
sxxi=syyi=szzi=sxyi=sxzi=syzi=0;
int mgi=0; //finest grid level
nubulkn = 0;
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	m4=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
	//m1=mgp[mgi]+10*mgxy[mgi]+10*mgy[mgi]+m2;
	//Only nodes with viscosity different from that of the rigid plates
	if(ro[m4]>1500 && m1>0 && m2>0 && m3>0) 
	{
	numtot += 1;
	//nxx = (nu[m4]+nu[m4-1]+nu[m4-mgy[mgi]]+nu[m4-mgy[mgi]-1])/4;
	nubulkn += 1/nxx[m4];
	/*   Z=0      Z=1 */
	/* 0  2      4  6 */
	/* 1  3      5  7 */
        v1[0]=mgp[mgi]+(m3-1)*mgxy[mgi]+(m1-1)*mgy[mgi]+(m2-1);v1[1]=v1[0]+1;
	v1[2]=v1[0]+mgy[mgi];v1[3]=v1[2]+1;
	v1[4]=v1[0]+mgxy[mgi];v1[5]=v1[4]+1;
	v1[6]=v1[4]+mgy[mgi];v1[7]=v1[6]+1;
	/* Stress and Strain rate second invariant */
	Exy = (exy[v1[0]]*exy[v1[0]] + exy[v1[1]]*exy[v1[1]] + exy[v1[2]]*exy[v1[2]] + exy[v1[3]]*exy[v1[3]])/4;
	Exz = (exz[v1[0]]*exz[v1[0]] + exz[v1[2]]*exz[v1[2]] + exz[v1[4]]*exz[v1[4]] + exz[v1[6]]*exz[v1[6]])/4;
	Eyz = (eyz[v1[0]]*eyz[v1[0]] + eyz[v1[1]]*eyz[v1[1]] + eyz[v1[4]]*eyz[v1[4]] + eyz[v1[5]]*eyz[v1[5]])/4;
	Sxy = (sxy[v1[0]]*sxy[v1[0]] + sxy[v1[1]]*sxy[v1[1]] + sxy[v1[2]]*sxy[v1[2]] + sxy[v1[3]]*sxy[v1[3]])/4;
	Sxz = (sxz[v1[0]]*sxz[v1[0]] + sxz[v1[2]]*sxz[v1[2]] + sxz[v1[4]]*sxz[v1[4]] + sxz[v1[6]]*sxz[v1[6]])/4;
	Syz = (syz[v1[0]]*syz[v1[0]] + syz[v1[1]]*syz[v1[1]] + syz[v1[4]]*syz[v1[4]] + syz[v1[5]]*syz[v1[5]])/4;
	eiinloc = pow(0.5*(exx[m4]*exx[m4]+eyy[m4]*eyy[m4]+ezz[m4]*ezz[m4]) + Exy + Exz + Eyz,0.5);
	siinloc = pow(0.5*(sxx[m4]*sxx[m4]+syy[m4]*syy[m4]+szz[m4]*szz[m4]) + Sxy + Sxz + Syz,0.5);	
//if(m1==30 && m3==30){printf("%ld %e %e %e %e %e %e %e %e %e\n",m4,gy[m2],nxx[m4],exx[m4],eyy[m4],ezz[m4],Exy,Exz,Eyz,eiinloc); getchar();}	
//if(m2==30 && m3==30){printf("%ld %e %e %e %e %e %e %e %e %e\n",m4,gx[m1],nxx[m4],sxx[m4],syy[m4],szz[m4],Sxy,Sxz,Syz,siinloc); getchar();}	
	/* VP Strain rate second invariant */
	Exy = (pow(sxy[v1[0]]/2/nxy[v1[0]],2) + pow(sxy[v1[1]]/2/nxy[v1[1]],2) + pow(sxy[v1[2]]/2/nxy[v1[2]],2) + pow(sxy[v1[3]]/2/nxy[v1[3]],2))/4;
	Exz = (pow(sxz[v1[0]]/2/nxz[v1[0]],2) + pow(sxz[v1[2]]/2/nxz[v1[2]],2) + pow(sxz[v1[4]]/2/nxz[v1[4]],2) + pow(sxz[v1[6]]/2/nxz[v1[6]],2))/4;
	Eyz = (pow(syz[v1[0]]/2/nyz[v1[0]],2) + pow(syz[v1[1]]/2/nyz[v1[1]],2) + pow(syz[v1[4]]/2/nyz[v1[4]],2) + pow(syz[v1[5]]/2/nyz[v1[5]],2))/4;
	eiivploc = pow(0.5*(pow(sxx[m4]/2/nxx[m4],2)+pow(syy[m4]/2/nxx[m4],2)+pow(szz[m4]/2/nxx[m4],2)) + Exy + Exz + Eyz,0.5);
	/* Heat dissipation */
	Sxy = (sxy[v1[0]]*sxy[v1[0]]/2/nxy[v1[0]]+sxy[v1[1]]*sxy[v1[1]]/2/nxy[v1[1]]+sxy[v1[2]]*sxy[v1[2]]/2/nxy[v1[2]]+sxy[v1[3]]*sxy[v1[3]]/2/nxy[v1[3]])/4;
	Sxz = (sxz[v1[0]]*sxz[v1[0]]/2/nxz[v1[0]]+sxz[v1[2]]*sxz[v1[2]]/2/nxz[v1[2]]+sxz[v1[4]]*sxz[v1[4]]/2/nxz[v1[4]]+sxz[v1[6]]*sxz[v1[6]]/2/nxz[v1[6]])/4;
	Syz = (syz[v1[0]]*syz[v1[0]]/2/nyz[v1[0]]+syz[v1[1]]*syz[v1[1]]/2/nyz[v1[1]]+syz[v1[4]]*syz[v1[4]]/2/nyz[v1[4]]+syz[v1[5]]*syz[v1[5]]/2/nyz[v1[5]])/4;
	hsnloc = (sxx[m4]*sxx[m4] + syy[m4]*syy[m4] + szz[m4]*szz[m4])/2/nxx[m4] + 2.0*(Sxy + Sxz + Syz);	
//if(m1==30 && m3==30){printf("%ld %e %e %e %e %e %e %e %e %e\n",m4,gy[m2],nxx[m4],sxx[m4],syy[m4],szz[m4],Sxy,Sxz,Syz,hsnloc); getchar();}	
        /* Sum up local values */
	eiin +=eiinloc;
        tau2 +=nxx[m4]*eiinloc;
	siin+= siinloc;
	hsn += hsnloc;
	eiivp += eiivploc;
	//Hard inclusions
	if(ro[m4]>3850){
	numinc += 1;
	exxi += exx[m4];
	eyyi += eyy[m4];
	ezzi += ezz[m4];
	exyi += (exy[v1[0]] + exy[v1[1]] + exy[v1[2]] + exy[v1[3]])/4;
	exzi += (exz[v1[0]] + exz[v1[2]] + exz[v1[4]] + exz[v1[6]])/4;
	eyzi += (eyz[v1[0]] + eyz[v1[1]] + eyz[v1[4]] + eyz[v1[5]])/4;
	sxxi += sxx[m4];
	syyi += syy[m4];
	szzi += szz[m4];
	sxyi += (sxy[v1[0]] + sxy[v1[1]] + sxy[v1[2]] + sxy[v1[3]])/4;
	sxzi += (sxz[v1[0]] + sxz[v1[2]] + sxz[v1[4]] + sxz[v1[6]])/4;
	syzi += (syz[v1[0]] + syz[v1[1]] + syz[v1[4]] + syz[v1[5]])/4;
	}
	else{
	nummat += 1;
	exxn += exx[m4];
	eyyn += eyy[m4];
	ezzn += ezz[m4];
	exyn += (exy[v1[0]] + exy[v1[1]] + exy[v1[2]] + exy[v1[3]])/4;
	exzn += (exz[v1[0]] + exz[v1[2]] + exz[v1[4]] + exz[v1[6]])/4;
	eyzn += (eyz[v1[0]] + eyz[v1[1]] + eyz[v1[4]] + eyz[v1[5]])/4;
	sxxn += sxx[m4];
	syyn += syy[m4];
	szzn += szz[m4];
	sxyn += (sxy[v1[0]] + sxy[v1[1]] + sxy[v1[2]] + sxy[v1[3]])/4;
	sxzn += (sxz[v1[0]] + sxz[v1[2]] + sxz[v1[4]] + sxz[v1[6]])/4;
	syzn += (syz[v1[0]] + syz[v1[1]] + syz[v1[4]] + syz[v1[5]])/4;
	}	
	//if(m1==10 && m3==10){printf("%d %e %e %e %e %e %e\n",m2,gy[m2],nu[m4],siinloc,eiinloc,siinloc*eiinloc,exy[m4]); getchar();}
	}
	}
eiin /= numtot; siin /= numtot; hsn /=numtot; eiivp /=numtot; nubulkn = numtot/nubulkn; tau2 /=numtot; 
if(nummat>0)
{
exxn /= nummat; eyyn /= nummat; ezzn /= nummat; exyn /= nummat; exzn /= nummat; eyzn /= nummat;
sxxn /= nummat; syyn /= nummat; szzn /= nummat; sxyn /= nummat; sxzn /= nummat; syzn /= nummat;
}
if(numinc>0)
{
exxi /= numinc; eyyi /= numinc; ezzi /= numinc; exyi /= numinc; exzi /= numinc; eyzi /= numinc;
sxxi /= numinc; syyi /= numinc; szzi /= numinc; sxyi /= numinc; sxzi /= numinc; syzi /= numinc;
}
nubulk1 = tau2/eiin;
nubulk2 = hsn/4/eiin/eiin;
nubulk3 = hsn/4/eiivp/eiivp;
nubulk4 = siin/2/eiin;
printf("\n Siin = %e Eiin = %e Eiivp = %e Tau2n = %e Hsn = %e Nubulkn = %e Tot Num = %ld\n",siin,eiin,eiivp,tau2,hsn,nubulkn,numtot);
printf("\n Nubulk1 = Tau2n/Eiin %e \n",nubulk1);
printf("\n Nubulk2 = Hs/4/Eiin^2 %e \n",nubulk2);
printf("\n Nubulk3 = Hs/4/Eiivp^2 %e \n",nubulk3);
printf("\n Nubulk4 = Siin/2/Eiin %e \n",nubulk4);
printf("\n STRAIN RATE TENSOR \n");
printf("\n Matrix     <exx>= %e <eyy>= %e <ezz>= %e <exy>= %e <exz>= %e <eyz>= %e Num Mat = %ld\n",exxn,eyyn,ezzn,exyn,exzn,eyzn,nummat);
printf("\n Inclusions <exx>= %e <eyy>= %e <ezz>= %e <exy>= %e <exz>= %e <eyz>= %e Num Inc = %ld\n",exxi,eyyi,ezzi,exyi,exzi,eyzi,numinc);
printf("\n DEVIATORIC STRESS TENSOR \n");
printf("\n Matrix     <sxx>= %e <syy>= %e <szz>= %e <sxy>= %e <sxz>= %e <syz>= %e Num Mat = %ld\n",sxxn,syyn,szzn,sxyn,sxzn,syzn,nummat);
printf("\n Inclusions <sxx>= %e <syy>= %e <szz>= %e <sxy>= %e <sxz>= %e <syz>= %e Num Inc = %ld\n",sxxi,syyi,szzi,sxyi,sxzi,syzi,numinc);
/* Print amount of markers which left the model */
//    char sc[]="";
//    strcat(sc,outdir);
//    strcat(sc,"effbulkvisc.txt");
//    fl = fopen(sc,"at");
//fl = fopen("/storage1/unipd/navarro/Newtonian/nucontr1000/single_sphere_LowRes_pureshear_r03/effbulkvisc.txt","at");
//fl = fopen("/storage1/unipd/navarro/Newtonian/nucontr1000/single_sphere_HiRes_pureshear_r005/effbulkvisc.txt","at");
fl = fopen("effbulkvisc.txt","at");
fprintf(fl,"%e %e %ld %ld %ld %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",timesum,timestep,numinc,nummat,numtot,exxn,eyyn,ezzn,exyn,exzn,eyzn,sxxn,syyn,szzn,sxyn,sxzn,syzn,exxi,eyyi,ezzi,exyi,exzi,eyzi,sxxi,syyi,szzi,sxyi,sxzi,syzi,siin,eiin,eiivp,hsn,nubulk2,nubulk3,nubulk4);
fclose(fl);
}
/**/
/**/
if(nodemark==2)
{
// Check strain
double eeincl=0,eematrix=0,eeplates=0,eetotal,gammatot;
long int mm1,numincl=0,nummatrix=0,numplates=0,numtotal;
for (mm1=0;mm1<marknum;mm1++)
	{
	//Inclusions
	if(markt[mm1]==7 || markt[mm1]==8) {eeincl+=marke[mm1]; numincl+=1;}
	//Matrix    
	if(markt[mm1]==10) {eematrix+=marke[mm1];nummatrix+=1;}
	//Plates    
	if(markt[mm1]==9) {eeplates+=marke[mm1]; numplates+=1;}
	}
//Gammatot=timesum*dVx/dy
gammatot=timesum*2.0/platesdistance;
printf("Gamma = %e --- Inclusion strain = %e (%e per cent) --- Matrix strain = %e (%e per cent) --- Plates strain = %e (%e per cent) --- Total strain = %e\n",gammatot,eeincl,eeincl/eetotal*100,eematrix,eematrix/eetotal*100,eeplates,eeplates/eetotal*100,eetotal);
/* Print amount of markers which left the model */
    char sc1[]="";
    strcat(sc1,outdir);
    strcat(sc1,"strainbudget.txt");
    fl = fopen(sc1,"at");
    //fl = fopen("strainbudget.txt","at");
fprintf(fl,"%e %e %e %e %e %ld %ld %ld %e %e\n",timesum,timestep,eeincl,eematrix,eeplates,numincl,nummatrix,numplates,eeincl/eetotal*100,eematrix/eetotal*100);
fclose(fl);
}
/**/
}

/* Move markers by using Runge-Kutta method */
void movemark()
{
/* Counters */
int n1,n0,ynreset=0,ynres[1000],ynres1,nt,tid;
long int mm1,marknum1,m1,m2,m3,m4,p[8],wn1[500];
/* Val buffer */
double X,Y,Z,v[3][5],xold,yold,zold,vxyz1[10],eps1[500];
double vxwater,vywater,vzwater,dx,dy,dz,ival,dxgrid=0,dpdx,dpdy,dpdz,x,y,z,xkf,ykf,zkf,vxkoef,vykoef,vzkoef;
double ee0,ee1,ee2,ee3,ee4;
double ss0,ss1,ss2,ss3,ss4;
double nubulk,nueff,mnu,eii;
float mfd[2];
/**/
/**/
/**/
/* Hydration front progress, melt extraction  */
if(timesum>1e+11)
	{
/*
	hydration2();
	meltextract();
	dxgrid=gridchange();
*/
	}
/**/
/**/
/* Save number of markers */
marknum1=marknum;
/**/
/**/
/**/
/* Move markers */
#pragma omp parallel shared(ynres,n0,bufv,ynreset,markx,marky,markz,markk,markt,markw,markd,marke,val1,vx,vy,vz,pr,exx,eyy,ezz,exy,exz,eyz,sxx,syy,szz,sxy,sxz,syz) private(ee0,ee1,ee2,ee3,ee4,ynres1,n1,vxyz,xold,yold,zold,mm1,p,v,vxwater,vywater,vzwater,dpdx,dpdy,dpdz,x,y,z,xkf,ykf,zkf,vxkoef,vykoef,vzkoef,X,Y,Z,m1,m2,m3,m4,ival,wn1,eps1,vxyz1) firstprivate(xperiodic,yperiodic,zperiodic,timestep,timesum,xnumx,ynumy,znumz,xynumxy,vyfluid,vymelt,nubeg,nuend,stp100,dxgrid,vdeep,zdeep,tdeep,xsize,ysize,zsize,outgrid,marknum,marknum1,gx,gy,gz,markll,marka0,marka1,markb0,markb1,marke0,marke1,marknu,markdh,markdv,markss,marks0,marks1,markmm,markn0,markn1)
{
/* Obtain cur thread number */
n1=omp_get_thread_num();
/* Obtain total number of threads */
if (n1==0) n0=omp_get_num_threads();
/* Reset marker reset identifier */
ynres1=0;
/*
*/
#pragma omp for schedule(static)// reduction(+:sii,siiplates,nueff,nubulkm,eii,eiiplates,nueffplates,hsplates,hs)
for (mm1=0;mm1<marknum;mm1++)
if((markx[mm1]>=0 && marky[mm1]>=0 && markz[mm1]>=0 && (double)(markx[mm1])<=xsize && (double)(marky[mm1])<=ysize && (double)(markz[mm1])<=zsize) || (outgrid!=1 && (markt[mm1]<50 || markt[mm1]>=100)))
	{
	n1=omp_get_thread_num();
	/* Save initial coordinates */
	xold=markx[mm1];
	yold=marky[mm1];
	zold=markz[mm1];
	/**/
	/* Velocity and second strain rate invariant calc */
	vxyzcalc1(markx[mm1],marky[mm1],markz[mm1],vxyz1,wn1);
	epscalc1(markx[mm1],marky[mm1],markz[mm1],1,eps1,wn1);
	val1[mm1]=pow(0.5*(eps1[4]*eps1[4]+eps1[5]*eps1[5]+eps1[6]*eps1[6])+eps1[7]*eps1[7]+eps1[8]*eps1[8]+eps1[9]*eps1[9],0.5);
	//printf("A %ld %d %e %e %e %e\n",mm1,markt[mm1],xold,yold,zold,eps1[10]);
	if(marke[mm1]>0) ee1=val1[mm1];
	//if(markt[mm1]!=9)
         //       {
                /* Compute effective viscosity */
                //mnu=viscalc(markk[mm1],eps1[10],mm1,markt[mm1],mfd);
        //        mnu=viscalc(273,eps1[10],mm1,markt[mm1],mfd);
        //        nubulkm += 1/mnu;
	//	}
	/**/
	//printf("B %ld %d %e %e %e\n",mm1,markt[mm1],xold,yold,zold);
	if(1==0){
	/* Water marker move */
	vxwater=vywater=vzwater=0;
	if(markt[mm1]>=50 && markt[mm1]<100) 
		{
		/* Water velocity */
		vywater=vyfluid; if(markd[mm1]>1100.0) vywater=vymelt;
		/* Fluid in rock */
		if(vyfluid>0 && (markk[mm1]==0 || markk[mm1]>298.0)) 
			{
			/* Horizontal,Vertical P-cell index */
			m1=wn1[0]; if(markx[mm1]>(gx[m1]+gx[m1+1])/2.0) m1+=1;
			if(m1<1) m1=1; if(m1>xnumx-2) m1=xnumx-2;
			m2=wn1[1]; if(marky[mm1]>(gy[m2]+gy[m2+1])/2.0) m2+=1;
			if(m2<1) m2=1; if(m2>ynumy-2) m2=ynumy-2;
			m3=wn1[2]; if(markz[mm1]>(gz[m3]+gz[m3+1])/2.0) m3+=1;
			if(m3<1) m3=1; if(m3>znumz-2) m3=znumz-2;
			/* Distances Calc */
			xkf=2.0/(gx[m1+1]-gx[m1-1]);
			ykf=2.0/(gy[m2+1]-gy[m2-1]);
			zkf=2.0/(gz[m3+1]-gz[m3-1]);
			/* Pressure gradients */
			x=xkf*(markx[mm1]-(gx[m1-1]+gx[m1])/2.0);
			y=xkf*(marky[mm1]-(gy[m2-1]+gy[m2])/2.0);
			z=xkf*(markz[mm1]-(gz[m3-1]+gz[m3])/2.0);
			/* P-SIGXX_EPSXX-Nodes num */
			/*    4----6 */
			/*   /|   /| */
			/*  / 5--/-7 */
			/* 0----2    */
			/* |    |    */
			/* 1----3    */
			p[0]=m3*xynumxy+m1*ynumy+m2;
			p[1]=p[0]+1;
			p[2]=p[0]+ynumy;
			p[3]=p[2]+1;
			p[4]=p[0]+xynumxy;
			p[5]=p[4]+1;
			p[6]=p[4]+ynumy;
			p[7]=p[6]+1;
			/* P derivatives */
			dpdx=xkf*((1.0-y)*(1.0-z)*(pr[p[2]]-pr[p[0]])+y*(1.0-z)*(pr[p[3]]-pr[p[1]])+(1.0-y)*z*(pr[p[6]]-pr[p[4]])+y*z*(pr[p[7]]-pr[p[5]]));
			dpdy=ykf*((1.0-x)*(1.0-z)*(pr[p[1]]-pr[p[0]])+x*(1.0-z)*(pr[p[3]]-pr[p[2]])+(1.0-x)*z*(pr[p[5]]-pr[p[4]])+x*z*(pr[p[7]]-pr[p[6]]));
			dpdz=zkf*((1.0-x)*(1.0-y)*(pr[p[4]]-pr[p[0]])+x*(1.0-y)*(pr[p[6]]-pr[p[2]])+(1.0-x)*y*(pr[p[5]]-pr[p[1]])+x*y*(pr[p[7]]-pr[p[3]]));
			/* Recalc velocity koefficients */
			vxkoef=(1000.0*GXKOEF-dpdx)/(2300.0*9.81);
			vykoef=(1000.0*GYKOEF-dpdy)/(2300.0*9.81);
			vzkoef=(1000.0*GZKOEF-dpdz)/(2300.0*9.81);
/*
printf("%ld %ld %ld    %ld %ld %ld   %e %e %e  %ld %d %e %e %e    %e %e %e   %e %e %e   %e %e %e",wn[0],wn[1],wn[2],m1,m2,m3,gx[m1],gy[m2],gz[m3],mm1,markt[mm1],markx[mm1],marky[mm1],markz[mm1],x,y,z,dpdx,dpdy,dpdz,vxkoef,vykoef,vzkoef);getchar();
*/
			if(vxkoef>2.0) vxkoef=2.0; if(vxkoef<-2.0) vxkoef=-2.0;
			if(vykoef>2.0) vykoef=2.0; if(vykoef<-2.0) vykoef=-2.0;
			if(vzkoef>2.0) vzkoef=2.0; if(vzkoef<-2.0) vzkoef=-2.0;
			/* Recalc velocity */
			vxwater=vywater*vxkoef;
			vzwater=vywater*vzkoef;
			vywater*=vykoef;
/*
printf("%ld %ld %ld %e %e %e  %ld %d %e %e %e    %e %e %e   %e %e %e   %e %e %e   %e %e %e ",m1,m2,m3,gx[m1],gy[m2],gz[m3],mm1,markt[mm1],markx[mm1],marky[mm1],markz[mm1],x,y,z,dpdx,dpdy,dpdz,vxkoef,vykoef,vzkoef,vxwater,vywater,vzwater);getchar();
*/
			}
		else
		/* Fluid in water */
			{
			vxwater=0;
			vywater=-ABSV(vywater);
			vzwater=0;
			}
		/**/
		}
/*
if(vywater>0) {printf("%e %e %e %e %e %e",tkpor,zmpor,vyfluid,vymelt,dmwamin,vywater);getchar();}
*/
	}
	/**/
	/* Motion Calc ///////////////////////////////// */
	if (markmod==1)
		{
		/* Simple Velocity calc */
/*
		vxyzcalc1(markx[mm1],marky[mm1],markz[mm1],vxyz1,wn1);
*/
		v[0][0]=vxyz1[0]+vxwater; v[1][0]=vxyz1[1]+vywater; v[2][0]=vxyz1[2]+vzwater;
		}
	else
		{
		/* 4 Runge-Kutta koef calc */
/*
		vxyzcalc1(markx[mm1],marky[mm1],markz[mm1],vxyz1,wn1);
*/
		v[0][1]=vxyz1[0]+vxwater; v[1][1]=vxyz1[1]+vywater; v[2][1]=vxyz1[2]+vzwater;
		/**/
		X=(double)(markx[mm1]+v[0][1]*timestep/2.0);
		Y=(double)(marky[mm1]+v[1][1]*timestep/2.0);
		Z=(double)(markz[mm1]+v[2][1]*timestep/2.0);
		if(xperiodic)
			{
                        if(X<0) {X=xsize+X;}
                        if(X>xsize) {X=X-xsize;}
                        if(X<0 || X>xsize) {printf("RK1 MARKER OUT OF X DOMAIN %ld %e %e %e %e\n",mm1,X,markx[mm1],v[0][1],timestep); exit(0);}
                        m1=m1serch(X);
			}
		if(yperiodic)
			{
                        if(Y<0) {Y=ysize+Y;}
                        if(Y>ysize) {Y=Y-ysize;}
                        if(Y<0 || Y>ysize) {printf("RK MARKER OUT OF Y DOMAIN\n"); exit(0);}
                        m2=m2serch(Y);
			}
		if(zperiodic)
			{
                        if(Z<0) {Z=zsize+Z;}
                        if(Z>xsize) {Z=Z-zsize;}
                        if(Z<0 || Z>zsize) {printf("RK MARKER OUT OF Z DOMAIN\n"); exit(0);}
                        m3=m3serch(Z);
			}
		vxyzcalc1(X,Y,Z,vxyz1,wn1);
		if(marke[mm1]>0) epscalc1(X,Y,Z,1,eps1,wn1);
		/**/
		v[0][2]=vxyz1[0]+vxwater; v[1][2]=vxyz1[1]+vywater; v[2][2]=vxyz1[2]+vzwater;
		if(marke[mm1]>0) ee2=pow(0.5*(eps1[4]*eps1[4]+eps1[5]*eps1[5]+eps1[6]*eps1[6])+eps1[7]*eps1[7]+eps1[8]*eps1[8]+eps1[9]*eps1[9],0.5);
		/**/
		X=(double)(markx[mm1]+v[0][2]*timestep/2.0);
		Y=(double)(marky[mm1]+v[1][2]*timestep/2.0);
		Z=(double)(markz[mm1]+v[2][2]*timestep/2.0);
		if(xperiodic)
			{
                        if(X<0) {X=xsize+X;}
                        if(X>xsize) {X=X-xsize;}
                        if(X<0 || X>xsize) {printf("RK2 MARKER OUT OF X DOMAIN %ld %e %e %e %e\n",mm1,X,markx[mm1],v[0][2],timestep); exit(0);}
                        m1=m1serch(X);
			}
		if(yperiodic)
			{
                        if(Y<0) {Y=ysize+Y;}
                        if(Y>ysize) {Y=Y-ysize;}
                        if(Y<0 || Y>ysize) {printf("RK MARKER OUT OF Y DOMAIN\n"); exit(0);}
                        m2=m2serch(Y);
			}
		if(zperiodic)
			{
                        if(Z<0) {Z=zsize+Z;}
                        if(Z>xsize) {Z=Z-zsize;}
                        if(Z<0 || Z>zsize) {printf("RK MARKER OUT OF Z DOMAIN\n"); exit(0);}
                        m3=m3serch(Z);
			}
		vxyzcalc1(X,Y,Z,vxyz1,wn1);
		if(marke[mm1]>0) epscalc1(X,Y,Z,1,eps1,wn1);
		/**/
		v[0][3]=vxyz1[0]+vxwater; v[1][3]=vxyz1[1]+vywater; v[2][3]=vxyz1[2]+vzwater;
		if(marke[mm1]>0) ee3=pow(0.5*(eps1[4]*eps1[4]+eps1[5]*eps1[5]+eps1[6]*eps1[6])+eps1[7]*eps1[7]+eps1[8]*eps1[8]+eps1[9]*eps1[9],0.5);
		/**/
		X=(double)(markx[mm1]+v[0][3]*timestep);
		Y=(double)(marky[mm1]+v[1][3]*timestep);
		Z=(double)(markz[mm1]+v[2][3]*timestep);
		if(xperiodic)
			{
                        if(X<0) {X=xsize+X;}
                        if(X>xsize) {X=X-xsize;}
                        if(X<0 || X>xsize) {printf("RK3 MARKER OUT OF X DOMAIN %ld %e %e %e %e\n",mm1,X,markx[mm1],v[0][3],timestep); exit(0);}
                        m1=m1serch(X);
			}
		if(yperiodic)
			{
                        if(Y<0) {Y=ysize+Y;}
                        if(Y>ysize) {Y=Y-ysize;}
                        if(Y<0 || Y>ysize) {printf("RK MARKER OUT OF Y DOMAIN\n"); exit(0);}
                        m2=m2serch(Y);
			}
		if(zperiodic)
			{
                        if(Z<0) {Z=zsize+Z;}
                        if(Z>xsize) {Z=Z-zsize;}
                        if(Z<0 || Z>zsize) {printf("RK MARKER OUT OF Z DOMAIN\n"); exit(0);}
                        m3=m3serch(Z);
			}
		vxyzcalc1(X,Y,Z,vxyz1,wn1);
		if(marke[mm1]>0) epscalc1(X,Y,Z,1,eps1,wn1);
		/**/
		v[0][4]=vxyz1[0]+vxwater; v[1][4]=vxyz1[1]+vywater; v[2][4]=vxyz1[2]+vzwater;
		if(marke[mm1]>0) ee4=pow(0.5*(eps1[4]*eps1[4]+eps1[5]*eps1[5]+eps1[6]*eps1[6])+eps1[7]*eps1[7]+eps1[8]*eps1[8]+eps1[9]*eps1[9],0.5);
		/**/
		/* Vx,Vy,Vz calc after Runge-Kutta */
		v[0][0]=(v[0][1]+2.0*v[0][2]+2.0*v[0][3]+v[0][4])/6.0;
		v[1][0]=(v[1][1]+2.0*v[1][2]+2.0*v[1][3]+v[1][4])/6.0;
		v[2][0]=(v[2][1]+2.0*v[2][2]+2.0*v[2][3]+v[2][4])/6.0;
		if(marke[mm1]>0) ee0=(ee1+2.0*ee2+2.0*ee3+ee4)/6.0;
		}
	/**/
	/**/
	/**/
	/* Orthogonal motion only */
	if (outgrid==2)
		{
		if(xperiodic)
			{
			if(marky[mm1]<0 || (double)(marky[mm1])>ysize) v[0][0]=v[2][0]=0;		
			if(markz[mm1]<0 || (double)(markz[mm1])>zsize) v[0][0]=v[1][0]=0;		
			}
		else if(yperiodic)
			{
			if(markx[mm1]<0 || (double)(markx[mm1])>xsize) v[1][0]=v[2][0]=0;		
			if(markz[mm1]<0 || (double)(markz[mm1])>zsize) v[0][0]=v[1][0]=0;		
			}
		else if(zperiodic)
			{
			if(markx[mm1]<0 || (double)(markx[mm1])>xsize) v[1][0]=v[2][0]=0;		
			if(marky[mm1]<0 || (double)(marky[mm1])>zsize) v[0][0]=v[2][0]=0;		
			}
		else        
			{
			if(markx[mm1]<0 || (double)(markx[mm1])>xsize) v[1][0]=v[2][0]=0;		
			if(marky[mm1]<0 || (double)(marky[mm1])>ysize) v[0][0]=v[2][0]=0;		
			if(markz[mm1]<0 || (double)(markz[mm1])>zsize) v[0][0]=v[1][0]=0;		
			}
		}
	/**/
	/**/
	/**/
	/* Normal/Immobile markers */
	if(markt[mm1]<100)
		{
		/* X,Y,Z calc after Runge-Kutta */
		markx[mm1]+=timestep*v[0][0]-dxgrid;
		marky[mm1]+=timestep*v[1][0];
		markz[mm1]+=timestep*v[2][0];
		/**/
		/* Periodic boundary condition for markers */
                if(xperiodic)
                        {
                        if(markx[mm1]<0) markx[mm1]=xsize+markx[mm1]; 
                        if(markx[mm1]>xsize) markx[mm1]=markx[mm1]-xsize;
                        if(markx[mm1]<0 || markx[mm1]>xsize) {printf("MARKER OUT OF X DOMAIN\n"); exit(0);}
                        }
                if(yperiodic)
                        {
                        if(marky[mm1]<0) marky[mm1]=ysize+marky[mm1]; 
                        if(marky[mm1]>ysize) marky[mm1]=marky[mm1]-ysize;
                        if(marky[mm1]<0 || marky[mm1]>ysize) {printf("MARKER OUT OF Y DOMAIN\n"); exit(0);}
                        }
                if(zperiodic)
                        {
                        if(markz[mm1]<0) markz[mm1]=zsize+markz[mm1]; 
                        if(markz[mm1]>zsize) markz[mm1]=markz[mm1]-zsize;
                        if(markz[mm1]<0 || markz[mm1]>zsize) {printf("MARKER OUT OF Z DOMAIN\n"); exit(0);}
                        }
		/* Update brittle strain */
		if(marke[mm1]>0) marke[mm1]+=(float)(ee0*timestep);
		/**/
		/* Out of grid marker reset */
		if(markx[mm1]<0 || marky[mm1]<0 || markz[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize || (double)(markz[mm1])>zsize) 
			{
			markk[mm1]=0;
			}
		}
	else
		{
		/* Immobile markers */
		/* X,Y,Z calc after Runge-Kutta */
		markx[mm1]+=timestep*v[0][0]-dxgrid;
		marky[mm1]+=timestep*v[1][0];
		markz[mm1]+=timestep*v[2][0];
                if(xperiodic)
                        {
                        if(markx[mm1]<0) markx[mm1]=xsize+markx[mm1]; 
                        if(markx[mm1]>xsize) markx[mm1]=markx[mm1]-xsize; 
                        if(markx[mm1]<0 || markx[mm1]>xsize) {printf("MARKER OUT OF X DOMAIN\n"); exit(0);}
                        }
                if(yperiodic)
                        {
                        if(marky[mm1]<0) marky[mm1]=ysize+marky[mm1]; 
                        if(marky[mm1]>ysize) marky[mm1]=marky[mm1]-ysize;
                        if(marky[mm1]<0 || marky[mm1]>ysize) {printf("MARKER OUT OF Y DOMAIN\n"); exit(0);}
                        }
                if(zperiodic)
                        {
                        if(markz[mm1]<0) markz[mm1]=zsize+markz[mm1]; 
                        if(markz[mm1]>zsize) markz[mm1]=markz[mm1]-zsize;
                        if(markz[mm1]<0 || markz[mm1]>zsize) {printf("MARKER OUT OF Z DOMAIN\n"); exit(0);}
                        }
		/* Displacement calc */
		if(!ynres1)
			{
			dx=markx[mm1]-markk[mm1];
			dy=marky[mm1]-markd[mm1];
			dz=markz[mm1]-markw[mm1];
			ival=pow(dx*dx+dy*dy+dz*dz,0.5);
			/* Generate new markers, reset immobile marker positions y/n  */
			if(ival>stp100) ynres1=1;
			}
		}
	/**/
	/**/
	/**/
	/* Motion Calc ///////////////////////////////// */
	}
ynres[n1]=ynres1;
}
/* Check if there is a need for reset */
for(n1=0;n1<n0;n1++)
	{
	if(ynres[n1]) ynreset=1;
	}
/**/
/* Reset immobile markers */
if(ynreset && 1==0)
	{
	if(printmod) printf("\n RESET immobile markers, Generate new markers OLD = %ld NEW = %ld \n",marknum,marknum1);
	for (mm1=0;mm1<marknum;mm1++)
	if(markt[mm1]>=100)
		{
		/* Create new normal marker */
		if(markx[mm1]>0 && marky[mm1]>0 && markz[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (double)(markz[mm1])<zsize)
			{
/*
printf("%e %e %e   %e %e %e   %e %e %e   %e %e %e",xold,yold,zold,markx[mm1],marky[mm1],markz[mm1],xnew,ynew,znew,v[0][0],v[1][0],v[2][0]);getchar();
*/
			markt[marknum1]=markt[mm1]-100;
			markx[marknum1]=markx[mm1];
			marky[marknum1]=marky[mm1];
			markz[marknum1]=markz[mm1];
			markk[marknum1]=0;
			markd[marknum1]=-1.0;
			markw[marknum1]=-1.0;
			val1[marknum1]=val1[mm1];
			/*
			markex[marknum1]=0;
			marktm[marknum1]=0;
			markc1[marknum1]=0;
			markc2[marknum1]=0;
			*/
			/* Add aditional markers counter */
			marknum1++;
			}
		/* X,Y reset for immobile marker */
		markx[mm1]=markk[mm1];
		marky[mm1]=markd[mm1];
		markz[mm1]=markw[mm1];
		}
	}
/**/
/* New Mark num check */
if(printmod) printf("\n NEW Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
if(marknum1>marknum0*MRKNEW) {printf("Space out in markx[]"); exit(0);}
/**/
/**/
/**/
/* Reset aditional markers */
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if(1==0 && markt[mm1]<100 && (markx[mm1]<0 || marky[mm1]<0 || markz[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize || (double)(markz[mm1])>zsize))
		{
		/* Decrease aditional markers counter */
		marknum1--;
		/* Type save */
		markt[mm1]=markt[marknum1];
		/* Temperature, Density, Water, Reset */
		markk[mm1]=0;
		markd[mm1]=-1.0;
		markw[mm1]=-1.0;
		/*
		markex[mm1]=0;
		marktm[mm1]=0;
		markc1[mm1]=0;
		markc2[mm1]=0;
		*/
		/* X,Y reload */
		markx[mm1]=markx[marknum1];
		marky[mm1]=marky[marknum1];
		markz[mm1]=markz[marknum1];
		val1[mm1]=val1[marknum1];
		}
	/* Increase markers counter */
	mm1++;
	}
if(printmod) printf("\n Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
/* Set new marker number */
marknum=marknum1;
/**/
}
/* End Move markers by using Runge-Kutta method */




/* ro[],nu[] recalc after marker positions */
void ronurecalc()
{
/* Counters */
long int m1,m2,m3,m4;
long int m10,m20,m30;
int mm2,yn,i,tid,nt;
long int mm1,wn1[100];
double dx,dy,dz,ddx,ddy,ddz;
double mnu,mpb,mtk,mro,mcp,mht,mbb,mkt,mwa,dmwa,swt,swt1,ival;
double mro0,mcp0,mbb0;
double anu,aro,acp,akt,aht,abb,markwlev;
double wnu,wro,wcp,wkt,wht,wbb,dywa,eps1[100],wi1[100],ival2;
float mfd[2];
double wtime;
/**/
/**/
/**/
/* Layering on sediments */
sedimnum++;
m1=(long int)(sedimnum/sedimcyc);
m2=((long int)(m1/2))*2;
if(m2==m1) yn=3; else yn=4;
/**/
/**/
#pragma omp parallel
{nt=omp_get_num_threads();}
/**/
/**/
/* ADD MARKERS TO THE v-CELLS ========================== */
wtime=omp_get_wtime();
/* Clear wt */
/*
#pragma omp parallel for shared(Sol0,Sol1,Ro,Nu,Nxx,Nxy,Nxz,Nyz) private(m1) firstprivate(nodenum,nodenum2) schedule(static) 
*/
/* Obtain cur thread number */
for (tid=0;tid<nt;tid++)
{
#pragma omp parallel for
for (m1=0;m1<nodenum;m1++)
	{
	Sol0[tid][m1]=0;
	Sol1[tid][m1]=0;
	Sol0[tid][nodenum+m1]=0;
	Sol1[tid][nodenum+m1]=0;
	Sol0[tid][nodenum2+m1]=0;
	Sol1[tid][nodenum2+m1]=0;
	Ro[tid][m1]=0;
	Nu[tid][m1]=0;
	Nxx[tid][m1]=0;
	Nxy[tid][m1]=0;
	Nxz[tid][m1]=0;
	Nyz[tid][m1]=0;
	}
}
/**/
/*
printf("%ld %e %e %e    %e %e %e",marknum,xsize,ysize,zsize,markx[0],marky[0],markz[0]); getchar();
*/
/* Add ro[] nu[] */
/*
#pragma omp parallel for shared(sol0,sol1,ro,nu,cp,pr,et,kt,tk,ht,Sol0,Sol1,Ro,Nu,Nxx,Nxy,Nxz,Nyz,markx,marky,markz,markt,markk,markd,markw,marke) private(i,tid,mfd,wn1,wi1,eps1,m1,m2,m3,m4,m10,m20,m30,mm2,mm1,dx,dy,dz,ddx,ddy,ddz,mnu,mpb,mtk,mro,mcp,mht,mbb,mkt,mwa,dmwa,swt,ival,ival2,mro0,mcp0,mbb0,anu,aro,acp,akt,aht,abb,markwlev,wnu,wro,wcp,wkt,wht,wbb,dywa) firstprivate(yn,marknum,nodenum,xsize,ysize,zsize,eroslev,sedilev,waterlev,gx,gy,gz,markbb,markht,markkt,markkf,markcp,markkp,nubeg,nuend,zdeep,vdeep,nudeep,tdeep,dtdeep,drdeep,mnumy,mnumx,xynumxy,ynumy,Mu,mu,mucnst,densimod,viscmod,nukoef,timesum,timeend) schedule(static)
*/
#pragma omp parallel for shared(Sol0,Sol1,Ro,Nu,Nxx,Nxy,Nxz,Nyz,markx,marky,markz,markt) private(tid,mm1,mm2,wn1,eps1,m1,m2,m3,m4,m10,m20,m30,dx,dy,dz,ddx,ddy,ddz,mnu,mpb,mtk,mro,mro0,mfd,swt,swt1,ival,ival2) firstprivate(marknum,nodenum,nodenum2,xsize,ysize,zsize,gx,gy,gz,markbb,markaa,nubeg,nuend,xynumxy,ynumy,densimod,viscmod,nukoef) schedule(static)
for (mm1=0;mm1<marknum;mm1++)
if(markx[mm1]>0 && marky[mm1]>0 && markz[mm1]>0 && markx[mm1]<xsize && marky[mm1]<ysize && markz[mm1]<zsize && markt[mm1]<50)// && markk[mm1]>0)
	{
	/* Thread id number */
	tid=omp_get_thread_num();
	/* Marker type */
	mm2=markt[mm1];
	/**/
	prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
	mpb=eps1[10]*1e-5;
	mtk=273;//(double)(markk[mm1]);
	/**/
	/* Marker Properties */
	mnu=viscalc(mtk,mpb,mm1,mm2,mfd);
	mro0=mro=dencalc1(mtk,mpb,mm2,eps1);
	/* Limit viscosity of the grid */
	if(mnu<nubeg) mnu=nubeg;
	if(mnu>nuend) mnu=nuend;
	/* Limit viscosity */
	/**/
	/* Cell num calc for A-D Diagonal Corners of cur marker */
	m10=wn1[0];
        m20=wn1[1];
        m30=wn1[2];
        dx=gx[m10+1]-gx[m10];
	dy=gy[m20+1]-gy[m20];
	dz=gz[m30+1]-gz[m30];
	/**/
	//printf("%ld %d %e %e %e %e\n",mm1,mm2,markx[mm1],marky[mm1],markz[mm1],mpb); getchar();
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
if(m10==0 && m30==0) {printf("%ld  %ld %ld %ld  %e %e %e   %e %e %e   %e %e %e   %e %e",mm1,m1,m2,m3,markx[m1],marky[m2],markz[m3],dx,dy,dz,ddx,ddy,ddz,swt,ival);getchar();}
*/
				/**/
				/* Node num */
				m4=m3*xynumxy+m1*ynumy+m2;
				/**/
				/* Add viscosity for basic nodes */
				if(ddx>nukoef && ddy>nukoef && ddz>nukoef)
					{
					/* Add Viscosity */
					if(viscmod==0) ival2=mnu*swt;
					if(viscmod==1) ival2=log(mnu)*swt;
					if(viscmod==2) ival2=1.0/mnu*swt;
					Nu[tid][m4]+=ival2;
					Sol1[tid][m4]+=swt;
					}
				/**/
				if(0==0)
				{
				/* Add viscosity for SIGxx, SIGyy, SIGzz (cell center) */
				if(m1>m10 && m2>m20 && m3>m30)
					{
					/* Add Viscosity */
					swt1=8.0*(0.5-ABSV(ddx-0.5))*(0.5-ABSV(ddy-0.5))*(0.5-ABSV(ddz-0.5));
					if(viscmod==0) ival2=mnu*swt1;
					if(viscmod==1) ival2=log(mnu)*swt1;
					if(viscmod==2) ival2=1.0/mnu*swt1;
					Nxx[tid][m4]+=ival2;
					Sol0[tid][nodenum+m4]+=swt1;
					}
				/**/
				/* Add viscosity for SIGxy (z-edges) */
				if(m3==m30 && ddx>0.5 && ddy>0.5)
					{
					/* Add Viscosity */
					swt1=8.0*(ddx-0.5)*(ddy-0.5)*(0.5-ABSV(ddz-0.5));
					if(viscmod==0) ival2=mnu*swt1;
					if(viscmod==1) ival2=log(mnu)*swt1;
					if(viscmod==2) ival2=1.0/mnu*swt1;
					Nxy[tid][m4]+=ival2;
					Sol1[tid][nodenum+m4]+=swt1;
					}
				/**/
				/* Add viscosity for SIGxz (y-edges) */
				if(m2==m20 && ddx>0.5 && ddz>0.5)
					{
					/* Add Viscosity */
					swt1=8.0*(ddx-0.5)*(0.5-ABSV(ddy-0.5))*(ddz-0.5);
					if(viscmod==0) ival2=mnu*swt1;
					if(viscmod==1) ival2=log(mnu)*swt1;
					if(viscmod==2) ival2=1.0/mnu*swt1;
					Nxz[tid][m4]+=ival2;
					Sol0[tid][nodenum2+m4]+=swt1;
					}
				/**/
				/* Add viscosity for SIGyz (x-edges) */
				if(m1==m10 && ddy>0.5 && ddz>0.5)
					{
					/* Add Viscosity */
					swt1=8.0*(0.5-ABSV(ddx-0.5))*(ddy-0.5)*(ddz-0.5);
					if(viscmod==0) ival2=mnu*swt1;
					if(viscmod==1) ival2=log(mnu)*swt1;
					if(viscmod==2) ival2=1.0/mnu*swt1;
					Nyz[tid][m4]+=ival2;
					Sol1[tid][nodenum2+m4]+=swt1;
					}
}
				/**/
				/* Add other properties */
				if(ddx>0.5 && ddy>0.5 && ddz>0.5)
					{
					/* Add Properties */
					swt1=8.0*(ddx-0.5)*(ddy-0.5)*(ddz-0.5);
					ival2=mro*swt1;
					Ro[tid][m4]+=ival2;
					Sol0[tid][m4]+=swt1;
					}
/*
if(swt<0){printf("%ld  %e %e %e   %ld %ld %ld %ld  %e %e %e swt=%e ",mm1,dx,dy,dz,m1,m2,m3,m4,ddx,ddy,ddz,swt);getchar();}
if(swt<0){printf("%e %e  %e %e  %e %e ",gx[m1],gx[m1+1],gy[m2],gy[m2+1],gz[m3],gz[m3+1]);getchar();}
*/
				}
			}
		/* Use Zmax of marker corner */
		}
	}
/**/
/*
*/
#pragma omp parallel for shared(sol0,sol1,ro,nu,nxx,nxy,nxz,nyz) private(m1) firstprivate(nodenum) schedule(static)
for (m1=0;m1<nodenum;m1++)
	{
	sol0[m1]=0;
	sol1[m1]=0;
	sol0[nodenum+m1]=0;
	sol1[nodenum+m1]=0;
	sol0[nodenum2+m1]=0;
	sol1[nodenum2+m1]=0;
	nu[nodenum+m1]=0;
	ro[nodenum+m1]=0;
	nxx[nodenum+m1]=0;
	nxy[nodenum+m1]=0;
	nxz[nodenum+m1]=0;
	nyz[nodenum+m1]=0;
	}
/* Sum up thread informations */
for(tid=0;tid<nt;tid++)
{
//#pragma omp parallel for shared(Sol0,Sol1,Ro,Nu,Nxx,Nxy,Nxz,Nyz,sol0,sol1,nu,ro,nxx,nxy,nxz,nyz) private(tid,m1) firstprivate(nodenum,nodenum2) schedule(static) 
#pragma omp parallel for
for(m1=0;m1<nodenum;m1++)
	{
	if(Sol0[tid][m1]>0)
		{
		sol0[m1]+=Sol0[tid][m1];
		ro[nodenum+m1]+=Ro[tid][m1];
		}
	if(Sol1[tid][m1]>0)
		{
		sol1[m1]+=Sol1[tid][m1];
		nu[nodenum+m1]+=Nu[tid][m1];
		}
	if(Sol0[tid][nodenum+m1]>0)
		{
		sol0[nodenum+m1]+=Sol0[tid][nodenum+m1];
		nxx[nodenum+m1]+=Nxx[tid][m1];
		}
	if(Sol1[tid][nodenum+m1]>0)
		{
		sol1[nodenum+m1]+=Sol1[tid][nodenum+m1];
		nxy[nodenum+m1]+=Nxy[tid][m1];
		}
	if(Sol0[tid][nodenum2+m1]>0)
		{
		sol0[nodenum2+m1]+=Sol0[tid][nodenum2+m1];
		nxz[nodenum+m1]+=Nxz[tid][m1];
		}
	if(Sol1[tid][nodenum2+m1]>0)
		{
		sol1[nodenum2+m1]+=Sol1[tid][nodenum2+m1];
		nyz[nodenum+m1]+=Nyz[tid][m1];
		}
	/*
	if(m1==896292){printf("%d %ld %e %e %e %e\n",tid,m1,sol0[m1],tk[nodenum+m1]/sol0[m1],sol1[m1],nu[nodenum+m1]/sol1[m1]);printf("OK: %ld %e %e %e %e \n",m1,Sol0[tid][m1],Tk[tid][m1],Sol1[tid][m1],Nu[tid][m1]); getchar();}
	*/
	}
}
printf("\n Time for adding marker properties = %e \n",omp_get_wtime()-wtime);
/**/
/* Add last column infos to first column */
if(xperiodic)
	{
	for(m3=0;m3<znumz;m3++)
	for(m2=0;m2<ynumy;m2++)
		{
		m1=m3*xynumxy+m2;
		m4=m1+(xnumx-1)*ynumy;
		sol0[m1]+=sol0[m4];
		sol1[m1]+=sol1[m4];
		sol1[nodenum+m1]+=sol1[nodenum+m4];
		sol0[nodenum2+m1]+=sol0[nodenum2+m4];
		ro[nodenum+m1]+=ro[nodenum+m4];
		nu[nodenum+m1]+=nu[nodenum+m4];
		nxy[nodenum+m1]+=nxy[nodenum+m4];
		nxz[nodenum+m1]+=nxz[nodenum+m4];
		/*
		printf("%ld %ld %e %e %e\n",m1,m4,nu[nodenum+m1],sol1[m1],nu[nodenum+m1]/sol1[m1]); getchar();
		*/
		}
	}
/**/
if(yperiodic)
	{
	for(m3=0;m3<znumz;m3++)
	for(m1=0;m1<xnumx;m1++)
		{
		m2=m3*xynumxy+m1*ynumy;
		m4=m2+ynumy-1;
		sol0[m2]+=sol0[m4];
		sol1[m2]+=sol1[m4];
		sol1[nodenum+m2]+=sol1[nodenum+m4];
		sol1[nodenum2+m2]+=sol1[nodenum2+m4];
		ro[nodenum+m2]+=ro[nodenum+m4];
		nu[nodenum+m2]+=nu[nodenum+m4];
		nxy[nodenum+m2]+=nxy[nodenum+m4];
		nyz[nodenum+m2]+=nyz[nodenum+m4];
		/*
		printf("%ld %ld %e %e %e\n",m1,m4,nu[nodenum+m1],sol1[m1],nu[nodenum+m1]/sol1[m1]); getchar();
		*/
		}
	}
/**/
if(zperiodic)
	{
	for(m1=0;m1<xnumx;m1++)
	for(m2=0;m2<ynumy;m2++)
		{
		m3=m1*ynumy+m2;
		m4=m3+(znumz-1)*xynumxy;
		sol0[m3]+=sol0[m4];
		sol1[m3]+=sol1[m4];
		sol0[nodenum2+m3]+=sol0[nodenum2+m4];
		sol1[nodenum2+m3]+=sol1[nodenum2+m4];
		ro[nodenum+m3]+=ro[nodenum+m4];
		nu[nodenum+m3]+=nu[nodenum+m4];
		nxz[nodenum+m3]+=nxz[nodenum+m4];
		nyz[nodenum+m3]+=nyz[nodenum+m4];
		/*
		printf("%ld %ld %e %e\n",m3,m4,nyz[nodenum+m3],nyz[nodenum+m3]/sol1[nodenum2+m3]); getchar();
		printf("%ld %ld %e %e %e\n",m1,m4,nu[nodenum+m1],sol1[m1],nu[nodenum+m1]/sol1[m1]); getchar();
		*/
		}
	}
/**/
/* Recalc ro[] nu[] */
#pragma omp parallel for shared(sol0,sol1,ro,nu,nxx,nxy,nxz,nyz,cp,kt,et,ht,tk,fd) private(m3) firstprivate(nodenum,nubeg,nuend,viscmod) schedule(static)
/*
printf("sizeof(struct flexarray) = %zu\n", sizeof(float));
*/
for (m3=0;m3<nodenum;m3++)
	{
	/* Recompute viscosity */
	if(sol1[m3]>0)
		{
		nu[m3]=nu[nodenum+m3]/sol1[m3];
		if(viscmod==1) nu[m3]=exp(nu[m3]);
		if(viscmod==2) nu[m3]=1.0/nu[m3];
/*
		if(nukoef>0 && nu[nodenum+m3]>0) nu[m3]=exp(log(nu[m3])*(1.0-nukoef)+log(nu[nodenum+m3])*nukoef);
		if(nukoef>0 && ro[nodenum+m3]>0) ro[m3]=ro[m3]*(1.0-nukoef)+ro[nodenum+m3]*nukoef;
{printf("%ld %e",m3,nu[m3]);getchar();}
if(nu[m3]>3e+22) {printf("%ld %e",m3,nu[m3]);getchar();}
*/
		}
	if(sol0[nodenum+m3]>0)
		{
		nxx[m3]=nxx[nodenum+m3]/sol0[nodenum+m3];
		if(viscmod==1) nxx[m3]=exp(nxx[m3]);
		if(viscmod==2) nxx[m3]=1.0/nxx[m3];
		}
	if(sol1[nodenum+m3]>0)
		{
		nxy[m3]=nxy[nodenum+m3]/sol1[nodenum+m3];
		if(viscmod==1) nxy[m3]=exp(nxy[m3]);
		if(viscmod==2) nxy[m3]=1.0/nxy[m3];
		}
	if(sol0[nodenum2+m3]>0)
		{
		nxz[m3]=nxz[nodenum+m3]/sol0[nodenum2+m3];
		if(viscmod==1) nxz[m3]=exp(nxz[m3]);
		if(viscmod==2) nxz[m3]=1.0/nxz[m3];
		}
	if(sol1[nodenum2+m3]>0)
		{
		nyz[m3]=nyz[nodenum+m3]/sol1[nodenum2+m3];
		if(viscmod==1) nyz[m3]=exp(nyz[m3]);
		if(viscmod==2) nyz[m3]=1.0/nyz[m3];
		}
	/* Recompute other properties */
	if(sol0[m3]>0)
		{
		ro[m3]=ro[nodenum+m3]/sol0[m3];
/*
{printf("%ld %e",m3,ro[m3]);getchar();}
*/
		}
/*
	else
		{
		printf("%ld",m3);getchar();
		}
*/
	if(nu[m3]<nubeg) nu[m3]=nubeg;
	if(nu[m3]>nuend) nu[m3]=nuend;
	if(nxx[m3]<nubeg) nxx[m3]=nubeg;
	if(nxx[m3]>nuend) nxx[m3]=nuend;
	if(nxy[m3]<nubeg) nxy[m3]=nubeg;
	if(nxy[m3]>nuend) nxy[m3]=nuend;
	if(nxz[m3]<nubeg) nxz[m3]=nubeg;
	if(nxz[m3]>nuend) nxz[m3]=nuend;
	if(nyz[m3]<nubeg) nyz[m3]=nubeg;
	if(nyz[m3]>nuend) nyz[m3]=nuend;
	sol0[m3]=0;
	sol1[m3]=0;
	sol0[nodenum+m3]=0;
	sol1[nodenum+m3]=0;
	sol0[nodenum2+m3]=0;
	sol1[nodenum2+m3]=0;
	}
if(xperiodic)
	{
	/* Set last column nodes parameters equale to those of first column */
	for(m3=0;m3<znumz;m3++)
	for(m2=0;m2<ynumy;m2++)
		{
		m1=m3*xynumxy+m2;
		m4=m1+(xnumx-1)*ynumy;
		ro[m4]=ro[m1];
		nu[m4]=nu[m1];
		nxy[m4]=nxy[m1];
		nxz[m4]=nxz[m1];
		/*
		printf("%ld %ld %e %e %e %e %e %e\n",m1,m4,nu[m4],nu[m4-ynumy],ro[m4],ro[m4-ynumy],nxy[m4],nxy[m4-ynumy]); getchar();
		*/
		}
	}
/**/
if(yperiodic)
	{
	/* Set last column nodes parameters equale to those of first column */
	for(m3=0;m3<znumz;m3++)
	for(m1=0;m1<xnumx;m1++)
		{
		m2=m3*xynumxy+m1*ynumy;
		m4=m2+ynumy-1;
		ro[m4]=ro[m2];
		nu[m4]=nu[m2];
		nxy[m4]=nxy[m2];
		nyz[m4]=nyz[m2];
		/*
		printf("%ld %ld %e %e %e %e %e %e\n",m1,m4,nu[m4],nu[m4-ynumy],ro[m4],ro[m4-ynumy],tk[m4],tk[m4-ynumy]); getchar();
		*/
		}
	}
/**/
if(zperiodic)
	{
	for(m1=0;m1<xnumx;m1++)
	for(m2=0;m2<ynumy;m2++)
		{
		m3=m1*ynumy+m2;
		m4=m3+(znumz-1)*xynumxy;
		ro[m4]=ro[m3];
		nu[m4]=nu[m3];
		nxz[m4]=nxz[m3];
		nyz[m4]=nyz[m3];
		/*
		printf("%ld %ld %e\n",m3,m4,nyz[m4]); getchar();
		*/
		}
	}
/**/
if(xperiodic && zperiodic)
	{
	/* Set last column nodes parameters equale to those of first column */
        for(m2=0;m2<ynumy;m2++)
                {
		/* 0  4 */
		/* 2  6 */
		/* 2 */
                m4=m2+(xnumx-1)*ynumy;
                ro[m4]=ro[m2];
                nu[m4]=nu[m2];
                nxz[m4]=nxz[m2];
		/* 4 */
                m4=m2+(znumz-1)*xynumxy;
                ro[m4]=ro[m2];
                nu[m4]=nu[m2];
                nxz[m4]=nxz[m2];
		/* 6 */
                m4=m2+(znumz-1)*xynumxy+(xnumx-1)*ynumy;
                ro[m4]=ro[m2];
                nu[m4]=nu[m2];
                nxz[m4]=nxz[m2];
                /*
                printf("%ld %ld %e %e %e \n",m1,m4,nu[m4],ro[m4],nxz[m4]); getchar();
                */
                }
        }

/* Set Boundary iconditions for T */
/*
if (printmod) printf("\n AVERAGE TEMPERATURE CORRECTION FOR BOUNDARY CONDITIONS ...\n");
tkrecalc();
if (printmod) printf("AVERAGE TEMPERATURE OK!\n");
*/
/**/
/**/
/*
int mgi=0;
for (m3=0;m3<mgz[mgi];m3++)
for (m1=0;m1<mgx[mgi];m1++)
for (m2=0;m2<mgy[mgi];m2++)
	{
	m4=mgp[mgi]+m3*mgxy[mgi]+m1*mgy[mgi]+m2;
	if(m1==10 && m3==10){printf("%ld %e %e %e %e %e\n",m4,nu[m4],nxx[m4],nxy[m4],nxz[m4],nyz[m4]); getchar();}
	}
*/
/* ADD MARKERS TO THE v-CELLS ========================== */
}
/* End ro[],nu[] recalc after marker positions */



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





/* Calculation of TK current location by Linear Interpolation */
void tkcalc1(double x, double y, double z, double *eps1, long int *wn1, double *wi1)
/* x,y,z - XYZ location of point for Vx,Vy,Vz calc */
{
/* Vx,Vy-nodes num */
long int m1,m2,m3,m4,m10,m20,m30,m40;
int n1;
/* Relativ coord */
double dx,dy,dz,ddx,ddy,ddz,swt,ival;
/**/
/**/
/**/
/* Cell num calc for A-D Diagonal Corners of cur marker */
wn1[0]=m10=m1serch(x);
wn1[1]=m20=m2serch(y);
wn1[2]=m30=m3serch(z);
wi1[0]=dx=gx[m10+1]-gx[m10];
wi1[1]=dy=gy[m20+1]-gy[m20];
wi1[2]=dz=gz[m30+1]-gz[m30];
/**/
/* Clear TK */
eps1[0]=eps1[1]=eps1[2]=eps1[3]=0;
/* Add nodes around current cell */
n1=0;
ival=0;
for (m3=m30;m3<=m30+1;m3++)
	{
	/* Normalized Z-distance calc */
	ddz=(gz[m3]-z)/dz;
	ddz=1.0-ABSV(ddz);
	for (m1=m10;m1<=m10+1;m1++)
		{
		/* Normalized X-distance calc */
		ddx=(gx[m1]-x)/dx;
		ddx=1.0-ABSV(ddx);
		for (m2=m20;m2<=m20+1;m2++)
			{
			/* Normalized Y-distance calc */
			ddy=(gy[m2]-y)/dy;
			ddy=1.0-ABSV(ddy);
			/* Wt calc */
			wi1[3+n1]=swt=ddx*ddy*ddz;
			/**/
			/* Node num */
			wn1[3+n1]=m4=m3*xynumxy+m1*ynumy+m2;
/*
printf("%ld  %ld %ld %ld  %e %e %e   %e %e %e   %e %e %e   %e %d %e",m4,m1,m2,m3,x,y,z,dx,dy,dz,ddx,ddy,ddz,swt,n1,tk[m4]);getchar();
*/
			/**/
			/* Add TK for current node */
			eps1[0]+=swt*tk0[m4];
			eps1[1]+=swt*tk1[m4];
			eps1[2]+=swt*tk[m4];
			eps1[3]+=swt*tk2[m4];
			/* Add Wt */
			ival+=swt;
			/* Add counter */
			n1++;
			}
		}
	}
/*
if(eps1[2]<900.0){printf("T  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10,m20,m30,ival,eps1[2]);getchar();}
{printf("T  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10,m20,m30,ival,eps1[2]);getchar();}
if(ABSV(ival-1.0)>1e-6){printf("T  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10,m20,m30,ival,eps1[2]);getchar();}
*/
/**/
}
/* Calculation of TK for current location by Linear Interpolation */




/* Calculation of P for current location by Linear Interpolation */
void  prcalc1(double x, double y, double z, double *eps1, long int *wn1)
/* x,y,z - XYZ location of point for Vx,Vy,Vz calc */
{
/* Vx,Vy-nodes num */
long int m1,m2,m3,m10,m20,m30,m40,m10min[2],m20min[2],m30min[2];
/* Relativ coord */
double X,Y,Z,dx,dy,dz,ddx,ddy,ddz,swt,ival,mpb1;
/**/
/**/
/**/
/* Clear P */
eps1[10]=0;
/* Out of grid */
if(x<0) x=0; if(x>xsize) x=xsize;
if(y<0) y=0; if(y>ysize) y=ysize;
if(z<0) z=0; if(z>zsize) z=zsize;
/* Fing Cell indexes */
wn1[0]=m1=m1serch(x);
wn1[1]=m2=m2serch(y);
wn1[2]=m3=m3serch(z);
/**/
/* P Cell indexes */
m10min[0]=m1;
if(x>(gx[m10min[0]]+gx[m10min[0]+1])/2.0) m10min[0]+=1;
if(xperiodic && (m10min[0]<1 || m10min[0]>xnumx-2))
	{
	if(m10min[0]<1) X=xsize+x; else X=x;
        m10min[0]=xnumx-1;
	m10min[1]=1;
	}
else
	{
	if(m10min[0]<1) m10min[0]=1;
	if(m10min[0]>xnumx-2) m10min[0]=xnumx-2;
	m10min[1]=m10min[0]+1;
	X=x;
	}
/**/
m20min[0]=m2;
if(y>(gy[m20min[0]]+gy[m20min[0]+1])/2.0) m20min[0]+=1;
if(yperiodic && (m20min[0]<1 || m20min[0]>ynumy-2))
	{
	if(m20min[0]<1) Y=ysize+y; else Y=y;
        m20min[0]=ynumy-1;
	m20min[1]=1;
	}
else
	{
	if(m20min[0]<1) m20min[0]=1;
	if(m20min[0]>ynumy-2) m20min[0]=ynumy-2;
	m20min[1]=m20min[0]+1;
	Y=y;
	}
/**/
m30min[0]=m3;
if(z>(gz[m30min[0]]+gz[m30min[0]+1])/2.0) m30min[0]+=1;
if(zperiodic && (m30min[0]<1 || m30min[0]>znumz-2))
        {
        if(m30min[0]<1) Z=zsize+z; else Z=z;
        m30min[0]=znumz-1;
        m30min[1]=1;
        }
else
        {
	if(m30min[0]<1) m30min[0]=1;
	if(m30min[0]>znumz-2) m30min[0]=znumz-2;
	m30min[1]=m30min[0]+1;
        Z=z;
        }
/**/
/* P Cell dimensions */
if(xperiodic && m10min[0]==xnumx-1) dx=(gx[1]+gx[xnumx-1]-gx[xnumx-2])/2.0;
else dx=(gx[m10min[0]+1]-gx[m10min[0]-1])/2.0;
if(yperiodic && m20min[0]==ynumy-1) dy=(gy[1]+gy[ynumy-1]-gy[ynumy-2])/2.0;
else dy=(gy[m20min[0]+1]-gy[m20min[0]-1])/2.0;
if(zperiodic && m30min[0]==znumz-1) dz=(gz[1]+gz[znumz-1]-gz[znumz-2])/2.0;
else dz=(gz[m30min[0]+1]-gz[m30min[0]-1])/2.0;
/* Interpolate from 8 nodes */
ival=0;
for (m30=0;m30<=1;m30++)
	{
	/* Normalized Z-distance calc */
	if(zperiodic && m30min[0]==znumz-1 && m30==1) ddz=((gz[m30min[m30]]+gz[m30min[m30]-1])/2.0-(Z-zsize))/dz;
        else ddz=((gz[m30min[m30]]+gz[m30min[m30]-1])/2.0-Z)/dz;
//ddz=((gz[m30-1]+gz[m30])/2.0-z)/dz;
	if (m30==m30min[0]) ddz=-ddz;
	ddz=1.0-ddz;
	for (m10=0;m10<=1;m10++)
		{
		/* Normalized X-distance calc */
		if(xperiodic && m10min[0]==xnumx-1 && m10==1) ddx=((gx[m10min[m10]]+gx[m10min[m10]-1])/2.0-(X-xsize))/dx;
		else ddx=((gx[m10min[m10]]+gx[m10min[m10]-1])/2.0-X)/dx;
		if (m10min[m10]==m10min[0]) ddx=-ddx;
		ddx=1.0-ddx;
		for (m20=0;m20<=1;m20++)
			{
			/* Normalized Y-distance calc */
			if(yperiodic && m20min[0]==ynumy-1 && m20==1) ddy=((gy[m20min[m20]]+gy[m20min[m20]-1])/2.0-(Y-ysize))/dy;
                	else ddy=((gy[m20min[m20]]+gy[m20min[m20]-1])/2.0-Y)/dy;
			//ddy=((gy[m20-1]+gy[m20])/2.0-y)/dy;
			if (m20==m20min[0]) ddy=-ddy;
			ddy=1.0-ddy;
			/* Pos in vx[] */
			m40=m30min[m30]*xynumxy+m10min[m10]*ynumy+m20min[m20];
			/* Calc wt */
			swt=ddx*ddy*ddz;
/*
if(m1==0 && m2>80){printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e\n",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,ival,pr[m40]*1e-5,eps1[10]);getchar();}
*/
			/* Add P for current node */
			eps1[10]+=pr[m40]*swt;
			/* Add Wt for current node */
			ival+=swt;
			}
		}
	}
/* Lithostatic pressure */
/*
mpb1=(y-1e+4)*1e+5/3.0;
if(mpb1<0) mpb1=0;
if(eps1[10]<0) eps1[10]=0;
if(eps1[10]>2.0*mpb1) eps1[10]=2.0*mpb1;
*/
/**/
/*
if(ABSV(ival-1.0)>1e-6){printf("P  %e %e %e   %ld %ld %ld  %e %e",x,y,z,m10min,m20min,m30min,ival,vxyz[0]);getchar();}
*/
/**/
}
/* Calculation of P for current location by Linear Interpolation */


/* Calc ro for given P,T */
double dencalc1(double mtk, double mpb, int mm2, double *eps1)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm2 - Rock number */
{
/* Val buffer */
double ival;
eps1[20]=0;
/**/
/* Constant density */
if (densimod==0) return markro[mm2];
/**/
/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
/* Adiabatic term: al=bro/(1-bro*(Tk-298.15)) */
ival=markro[mm2]*(1.0-markbb[mm2]*(mtk-298.15))*(1.0+markaa[mm2]*(mpb-1.0)*1e-3);
eps1[20]=markbb[mm2]/(1.0-markbb[mm2]*(mtk-298.15));
return  ival;
}
/* End Calc ro for given P,T */




/* Nu calc after reological equation */
/* P-T-stress dependent rheology without/with brittle/ductile transition */
/* Reological equations */
/* Stress>SScr */
/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
/* Stress<SScr */
/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
/* NU1=NU0/SScr^(n-1) */
/* SScr - dislocation, diffusion transition stress */
/* SSii - second invariant of deviatoric stress tensor */
/* EEii - second invariant of strain rate tensor */
/* E - activation energy, J */
/* V - activation volume, J/bar */
/* R - gase constant 8.314 J/K */
/* Viscosity NU  calc after reological equations */
/* NU=SSii/(2*EEii) */
/* Brittle - Ductile transition */
/* sbrit=MINV(0.85e+5*pb,60e+6+0.6e+5*pb)*lambda;  (Schott & Scmeling, 1998) */
/* sbrit=MINV(0.667e+5*pb,51.2e+6+0.512e+5*pb)*lambda; (Brace & Kohlsstedt, 1980) */
double viscalc(double mtk, double mpb, long int mm1, int mm2, float *mfd)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - Marker number */
/* mm2 - rock type */
{
/* Val buffer */
double e,n,rt=8.314*mtk,k1,e1,epsin,sduct,sbrit,nueff,smin,smax,nmin,nmax,strain,abrit,bbrit,nubrit,nunewt,nupowl;
double k1min,A,sig0,p,q,Pepsin,nupeierls,sigmin,sigmax,siginnew,sigin,mpf;
double y,y1,ylowvisc,ymax,Ea,Eadepth,A0disl,A0diff;
/* Reological Eq par */
double lamb,mpb1;
/* Counters */
long int m1;
int n1,n2;
/* Melted rocks */
if(marke[mm1]==0.0) marke[mm1]=1e-20;
if (mm2>20) return markn0[mm2];
if (markn0[mm2] == markn1[mm2]) return markn0[mm2];
/**/
/**/
/*
printf("%e %e",marky[mm1],mpb);getchar();
*/
/* Calc effective strain rate after second strain rate Tenzor invariant EEii=(1/2SUM(EPSik^2))^(1/2) */
epsin=val1[mm1];
/**/
/* Mechanical non-Newtonian viscosity */
double eii0=1.0;
n=markmm[mm2];
nueff=marknu[mm2]*pow(epsin/eii0,(1-n)/n);
if(nueff<markn0[mm2]) nueff=markn0[mm2]; 
if(nueff>markn1[mm2]) nueff=markn1[mm2];
//printf("%d %e %e\n",mm2,epsin,nueff); getchar();
return nueff;
/********************************************************/
/********************************************************/
/********************************************************/
/********************************************************/
/********************************************************/
/********************************************************/
/* Brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/* Brittle/ductile transition stress calc */
lamb=markll[mm2]; 
abrit=marka0[mm2]; 
bbrit=markb0[mm2]; 
/* A,B coefficients calc depending on integral strain */
strain=marke[mm1];
if(strain>marke1[mm2])
        {
        abrit=marka1[mm2];
        bbrit=markb1[mm2];
        }
else
        {
        if(strain>marke0[mm2] && marke1[mm2]>marke0[mm2])
                {
                abrit=marka0[mm2]+(marka1[mm2]-marka0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
                bbrit=markb0[mm2]+(markb1[mm2]-markb0[mm2])*(strain-marke0[mm2])/(marke1[mm2]-marke0[mm2]);
                }
        }
/*
if(markx[mm1]<1e+5 || markx[mm1]>xsize-1e+5) 
	{
	if (markx[mm1]<1e+5) e=markx[mm1]/1e+5; else e=(xsize-markx[mm1])/1e+5;
	n=bbrit*e+0.1*(1.0-e);
	if(bbrit>n) bbrit=n; 
	}
*/
sbrit=abrit+bbrit*mpb*1e+5*lamb;
if (sbrit>marks1[mm2]) sbrit=marks1[mm2];
/* Yield stress limit for upper plate:  5, 17 rocktypes */
if((mm2==16 || mm2==17) && sbrit>yieldstress)
        {
        /*
        printf("%ld %d %e %e\n",mm1,mm2,sbrit,yieldstress); getchar();
        */
        sbrit=yieldstress;
        }
/* Limit stress */
/*
if(marky[mm1]>vdeep && sbrit>nudeep) sbrit=nudeep; 
if((markx[mm1]>xsize-5e+4 || markx[mm1]<5e+4) && sbrit>1e+8) sbrit=1e+8; 
*/
/* Inverted value of Brittle "Mohr-Coulomb" viscosity calc */
nubrit=0; if(epsin>0 && sbrit>0) nubrit=1.0/(0.5*sbrit/epsin);
if(epsin>0 && sbrit<=0 && (bbrit>0 || abrit*lamb>0)) nubrit=1.0/markn0[mm2];
/* End Brittle "viscosity" calc: sbrit=A(eps)+LAM*B(eps)*P ------------------- */
/**/
/**/
/* Peierls plasticity-creep mechanism, data from Katayama and Karato, 2008 */
nupeierls=0;
if(0==0 && mm2>1 && mm2<20 && mtk<1373.0 && epsin>0) 
	{
	n1=9;
	A=pow(10.0,7.8)*1e-12;/* Constant f(p,q), 1/s/MPa^2 */ 
	/* Dry Peierls stress at 0 K f(p,q), MPa Evans & Goetze, 1979 */
	sig0=9.1e+9; /* Dry Peierls stress at 0 K f(p,q), MPa Evans & Goetze, 1979 */
	if(0==1 && mm2!=9 && mm2!=10 && mm2!=14)
		{
		sig0=2.9e+9; /* Wet Peierls stress at 0 K f(p,q), MPa Katayama & Karato, 2008 */
		n1=11;
		}
	/* Using bisection */
	sigmin=1e+6;
	sigmax=sig0*0.9999;
	k1min=A*pow(sigmin,2.0)*exp(MAXV(-100.0,-(markdh[n1]+markdv[n1]*mpb)/rt*pow(1.0-sigmin/sig0,2.0)));
	for(n2=0;n2<15;n2++)
		{
		siginnew=(sigmax+sigmin)/2.0;
/*
		siginnew=pow(sigmax*sigmin,0.5);
*/
		k1=A*pow(siginnew,2.0)*exp(MAXV(-100,-(markdh[n1]+markdv[n1]*mpb)/rt*pow(1.0-siginnew/sig0,2.0)));
		if((k1<epsin && k1min<epsin) || (k1>epsin && k1min>epsin)) 
			{
			sigmin=siginnew;
			}
		else
			{
			sigmax=siginnew;
			}
/*
	p=1.0;
	q=2.0;
	k1=A*pow(siginnew,2.0)*exp(-(markdh[n1]+markdv[n1]*mpb)/rt*pow(1.0-pow(siginnew/sig0,p),q));
printf("%ld %d   %e %e   %e %e   %e %e",mm1,mm2,mtk,mpb,sigin,epsin,siginnew,1.0/nupeierls);getchar();
printf("%d   %ld %d   %e %e   %e %e   %e %e  %e %e     %e %e",n2,mm1,mm2,mtk,mpb,sigin,epsin,siginnew,k1,k1min,k1max,sigmin,sigmax);getchar();
*/
		}
		
	nupeierls=1.0/(0.5*siginnew/epsin);
	if(siginnew<1e+8) nupeierls=0;
/*
	if(nupeierls>1e-23 && (markx[mm1]<5e+4 || markx[mm1]>xsize-5e+4)) nupeierls=1e-23;
*/
	}
/*
if(siginnew>1e+8){printf("%ld %d   %e %e  %e %e %e %e",mm1,mm2,mtk,mpb,epsin,siginnew,1/nubrit,1/nupeierls);getchar();}
*/
/**/
/**/
/**/
/* Ductile viscosity calc -------------------------------------------*/
/* Inverted value of newtonian NU set */
nunewt=0;
/**/
/* Inverted value of power-low NU set */
nupowl=0;
/**/
/* Check for the presence of ductile rheology */
if (marknu[mm2])
	{
	/* A)  Simple Newtonian rheology */
	/* Newtonian creep: SSii=NU0*2.0*EEii */
	/* Effective viscosity: NU=NU0 */
	/* Effective viscosity member in Stoks: NUs=NU */
	if(markdh[mm2]==0 && markdv[mm2]==0 && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/marknu[mm2];
		}
	/**/
	/**/
	/**/
	/* B)  P-T dependent, stress independent Newtonian rheology */
	/* Newtonian diffusion creep: SSii=NU0*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0*exp[(E+PV)/RT] */
	if((markdh[mm2]!=0 || markdv[mm2]!=0) && (markss[mm2]==0 || markmm[mm2]==1.0))
		{
		/* Inverted value of newtonian NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		/* Test creep Moresi & Solomatov (1995): SSii=NU0*exp[-a*(T-T0)] */
		if(markdh[mm2]<0 && markdv[mm2]>0) 
			{
			e1=markdh[mm2]*(mtk-markdv[mm2]);
			if(e1<-150.0) e1=-150.0;
			}
		/* Test creep Turkotte & Schubert(1982): SSii=NU0*exp[E/RTo(1-(T-T0)/T0)] */
		if(markdh[mm2]<0 && markdv[mm2]<0) 
			{
			e1=(-markdh[mm2])*(1.0-(mtk-(-markdv[mm2]))/(-markdv[mm2]))/8.314/(-markdv[mm2]);
			if(e1>150.0) e1=150.0;
			}
		nunewt=1.0/(marknu[mm2]*exp(e1));
		}
	/**/
	/**/
	/**/
	/* C)  P-T independent, stress dependent rheology without/with brittle/ductile transition */
	/* Stress>SScr */
	/* Power law creep: SSii={NU0*EEii}^(1/n) */
	/* Effective viscosity: NU=1/2*NU0^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian creep: SSii=NU1*EEii */
	/* Effective viscosity: NU=NU1/2 */
	/* Effective viscosity member in Stoks: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if((markdh[mm2]==0 && markdv[mm2]==0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* Koef for stress independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/(0.5*k1);
		/**/
		/* Ductile power-low stress calc */
		sduct=pow(marknu[mm2]*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power-low NU calc */
		if (epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		}
	/**/
	/**/
	/**/
	/* D)  P-T-stress dependent rheology without/with brittle/ductile transition */
	/* Reological equations */
	/* Stress>SScr */
	/* Power law dislocation creep: SSii={NU0*EEii*exp[(E+PV)/RT]}^(1/n) */
	/* Effective viscosity: NU=1/2*{NU0*exp[(E+PV)/RT]}^(1/n)*EEii^[(1-n)/n] */
	/* Effective viscosity member in Stoks: NUs=NU/n */
	/* Stress<SScr */
	/* Newtonian diffusion   creep: SSii=NU1*EEii*exp[(E+PV)/RT] */
	/* Effective viscosity: NU=NU0/2*exp[(E+PV)/RT] */
	/* Effective viscosity member in Stoks: NUs=NU */
	/* NU1=NU0/SScr^(n-1) */
	if((markdh[mm2]!=0 || markdv[mm2]!=0) && markss[mm2]!=0 && markmm[mm2]!=1.0)
		{
		/* T-P exponent for effective NU calc */
		e1=(markdh[mm2]+markdv[mm2]*mpb)/rt;
		if(e1>150.0) e1=150.0;
		e1=exp(e1);
		/**/
		/* Koef for stress independent creep NU1 calc */
		k1=marknu[mm2]/pow(markss[mm2],markmm[mm2]-1.0);
		/**/
		/* Inverted value of newtonian NU calc */
		nunewt=1.0/(0.5*k1*e1);
		/**/
		/* Ductile power-low stress calc */
		sduct=pow(marknu[mm2]*e1*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power-low NU calc */
		if(epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		}
	// E) Depth dependent activation enthalpy
	if(markdv[mm2]>1e+5)
		{
		y=marky[mm1];
		//y=990e+3;
		// Calculate depth dependent activation enthalpy
		ylowvisc=000e+3; // Depth of low viscosity channel
		ymax= 9e+6; // increase this to increase viscosity in the lowermost mantle
		y1 = pow(y-ylowvisc,2.0)*(1/y-1/ymax);
		Ea=530e+3; // Activation energy
		Eadepth= 0.30; // coefficient of scaling for activation volume; increase this to increase viscosity proportionally to depth
		// Dislocation creep
		/* T-P exponent */
		//rt=8.313*1974.7911338;
		//epsin=1e-15;
		e1 = (Ea + Eadepth * y1)/rt;
		e1 = exp(e1);
		//if(mm2==9||mm2==10||mm2==15){printf("1 %ld %d %e\n",mm1,mm2,e1);}
		/* Ductile power-low stress calc */
		A0disl = 4e+17; //upper mantle pre-exponential factor
		if(mm2==15|| mm2==18 || mm2==19) A0disl = A0disl * 1000; //lower mantle pre-exponential factor
		sduct=pow(A0disl*e1*epsin,1.0/markmm[mm2]);
		/**/
		/* Inverted value of power-low NU calc */
		if(epsin>0) nupowl=1.0/(0.5*sduct/epsin);
		/**/
		// Diffusion creep
		/* T-P exponent */
		ylowvisc=150e+3; // Depth of low viscosity channel
		y1 = pow(y-ylowvisc,2.0)*(1/y-1/ymax);
		e1 = (Ea + Eadepth * y1)/rt/2.0;
		e1 = exp(e1);
		A0diff = 4e+16; //upper mantle pre-exponential factor
		if(mm2==15|| mm2==18 || mm2==19) A0diff = A0diff / 1000; //lower mantle pre-exponential factor
		nunewt = 1.0/(0.5*A0diff*e1);
		//if(mm2==9||mm2==10||mm2==15){printf("1 %ld %d %e %e %e\n",mm1,mm2,e1,A0diff,0.5*A0diff*e1);}
		}
	}
/* End Ductile viscosity calc -------------------------------------------*/
/**/
/**/
/**/
/* Inverted value of effective viscosity calc, check */
nueff=nunewt+nupowl;
/* Fraction of creep accomodated by power law creep */
mfd[0]=nupowl/nueff;
/*
if(mm2>1) {printf("A %ld %d %e %e %e %e %e %e\n",mm1,mm2,mtk,mpb,mfd[0],nupowl,nunewt,nueff); getchar();}
*/
if(mfd[0] > 1) {printf("Too high fractdisl: nupowl %e nueff %e \n",1/nupowl,1/nueff); exit(0);}
if(nupeierls>nueff) {nueff=nupeierls; mfd[0]=0;}
if(nubrit>nueff) {nueff=nubrit; mfd[0]=0;}
/*
if(mm2>1) {printf("B %ld %d %e %e %e %e %e %e\n",mm1,mm2,mtk,mpb,mfd[0],nupowl,nunewt,nueff); getchar();}
*/
/*
if(mm2==9||mm2==10||mm2==15){printf("2 %ld %d  %e %e %e %e %e %e %e %e\n",mm1,mm2,mtk,y1,epsin,1/nubrit,1/nupeierls,1/nunewt,1/nupowl,1/nueff);getchar();}
*/
/**/
/* Inverted Viscosity check */
if(nueff<=0) {nueff=1.0/markn1[mm2]; mfd[0]=0;}
/* Additional brittle strain */
if(marke[mm1]==0.0) marke[mm1]=1e-20;
//if(nueff==nubrit) marke[mm1]=ABSV(marke[mm1]);
//if(mm2==7||mm2==9||mm2==14) marke[mm1]=ABSV(marke[mm1]);
//else marke[mm1]=-ABSV(marke[mm1]);
/* Add all the strain for the slab */
/* Viscosity calc */
nueff=1.0/nueff;
/*
printf("A %ld %d  %e %e %e   %e %e  %e %e",mm1,mm2,markx[mm1],marky[mm1],markz[mm1],mtk,mpb,epsin,nueff);getchar();
printf("A %ld %d %e %e %e",mm1,mm2,mtk,mpb,nueff);getchar();
*/
/**/
/* Viscosity check */
if(nueff<markn0[mm2]) nueff=markn0[mm2]; 
if(nueff>markn1[mm2]) nueff=markn1[mm2];
/* Lower cutoff value in the lowermost mantle */
if(1==0 && marky[mm1]>1500e+3)
	{
	double cutoff;
        if(marky[mm1]>2500e+3) cutoff= 1e+22;
        else cutoff = markn1[mm2] - (markn1[mm2]-1e+22)*(marky[mm1]-1500e+3)/1000e+3;
        if(nueff>cutoff) nueff=cutoff;
	}
/*
if(nueff<nubeg) nueff=nubeg; 
if(nueff>nuend) nueff=nuend;
*/
/**/
/*
if(strain>1.5 && nubrit>0 && x>1000000.0 && y<400000.0) {printf("%ld %d  %e %e  %e %e  %e %e %e %e  %e %e %e",mm1,mm2,x,y,mtk,mpb,strain,abrit,bbrit,lamb,1.0/nubrit,epsin,nueff); getchar();}
*/
/* Return calculated viscosity */
/*
if(mm2==9 || mm2==14) {printf("B   %ld %d %e %e %e %e %e",mm1,mm2,mtk,mpb,nueff,markn0[mm2],markn1[mm2]);getchar();}
*/
return nueff;
}
/* Nu calc after reological equation */




/* Thermodynamic database use for ro, Cp */
void tdbasecalc1(double mtk, double mpb, int mm2, long int mm1, double *eps1)
{
/* TD Database variables,  dTK,dPB - TK, PB step for tabulation in TD database */
double H0,H1,H2,H3,R0,R1,R2,R3,W0,W1,W2,W3,n,e;
/* Val Buffers */
int n1,n2,mm3,ynout=0;
double mhh0,mhh1,mdhh,maa,mwa,dmwa,wro,mro,mcp,mbb,xh2o;
double sy1,e1;
/**/
/* Reset TD variables */
eps1[40]=eps1[41]=eps1[42]=eps1[43]=eps1[44]=eps1[45]=0;
/**/
/* TD base type */
switch (mm2)
	{
	/* Dry Upper crust */
	case 5:
	mm3=0 ; break;
	/* Wet Upper crust */
	case 17:
	mm3=0 ; break;
	/* Dry Lower crust */
	case 6:
	mm3=0 ; break;
	/* Wet Lower crust */
	case 18:
	mm3=0 ; break;
	/* Sediments */
	case 2:
	case 3:
	case 4: 
	mm3=0; break;
	/* Molten Sediments */
	case 22:
	case 23:
	case 24: 
	case 25: 
	case 35: 
	case 37:
	mm3=0; break;
	/* Basalts */
	case 7: 
	case 16:
	mm3=0; break;
	/* Molten Basalt */
	case 27: 
	case 36:
	mm3=0; break;
	/* Gabbro */
	case 8: 
	mm3=0; break;
	/* Molten Gabbro */
	case 26: 
	case 28: 
	case 38:
	mm3=0; break;
	/* Dry peridotite */
	case 9:
	case 10: 
	case 12:
	case 14:
	case 15: 
	case 19: 
	mm3=0; break;
	/* Wet peridotite */
	case 13:
	case 11: 
	mm3=0; break;
	/* Molten peridotite */
	case 29: 
	case 30: 
	case 31: 
	case 32: 
	case 34: 
	mm3=0; break;
	/* Unknown type */
	default: {printf("Unknown rock type for TD database %d",mm2); exit(0);}
	}
/* ABCD-4Cell Number */
e=(mtk-tkmin)/tkstp;
if(e<0) e=0;
if(e>(double)(tknum-1)) {ynout=1;e=(double)(tknum-1);}
n=(mpb-pbmin)/pbstp;
if(n<0) n=0;
if(n>(double)(pbnum-1)) {ynout=1;n=(double)(pbnum-1);}
n1=(int)(e);
if(n1>tknum-2) n1=tknum-2;
n2=(int)(n);
if(n2>pbnum-2) n2=pbnum-2;
/* e,n Calc */
e=(e-(double)(n1));
n=(n-(double)(n2));
/* Ro H values */
/* 0 2 */
/* 1 3 */
mm3=0;
R0=td[n1  ][n2  ][mm3][0]*1000.0;
R1=td[n1  ][n2+1][mm3][0]*1000.0;
R2=td[n1+1][n2  ][mm3][0]*1000.0;
R3=td[n1+1][n2+1][mm3][0]*1000.0;
H0=td[n1  ][n2  ][mm3][1]*1000.0*4.1837;
H1=td[n1  ][n2+1][mm3][1]*1000.0*4.1837;
H2=td[n1+1][n2  ][mm3][1]*1000.0*4.1837;
H3=td[n1+1][n2+1][mm3][1]*1000.0*4.1837;
W0=W1=W2=W3=0.0;
/*
W0=td[n1  ][n2  ][mm3][2];
W1=td[n1  ][n2+1][mm3][2];
W2=td[n1+1][n2  ][mm3][2];
W3=td[n1+1][n2+1][mm3][2];
*/
/* Ro calc by interpolation */
mro=((R0*(1.0-n)+R1*n)*(1.0-e)+(R2*(1.0-n)+R3*n)*e);
/* Water wt% calc by interpolation */
mwa=((W0*(1.0-n)+W1*n)*(1.0-e)+(W2*(1.0-n)+W3*n)*e);
if(1==0){
/* Add porocity fluid */
/* Erosion surface */
if(marks0[mm2]>0 && marky[mm1]<zmpor && mtk<tkpor) 
	{
	dmwa=marks0[mm2]*(tkpor-mtk)/(tkpor-273.15)*(zmpor-MAXV(marky[mm1],sedilev))/(zmpor-sedilev);
	mwa+=dmwa;
	wro=1050.0;
	mro=mro/(1.0+dmwa*1e-2*(mro/wro-1.0));
/*
if(sy1>10000.0){printf("TD1 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
if(sy1>10000.0){printf("TD2 %d %d %e %e   %e %e   %e %e",mm2,mm3,mtk-273.15,sy1,mwa,mro,dmwa,zmpor);getchar();}
*/
	}
/* Limit amount of fluid in partially molten mantle */
if(mm2>=29 && mm2<=34 && mwa>maxwater) mwa=maxwater;
/* Add porous water to the hydrated mantle in the region of fluid-fluxed melting */
if(mm2==11 && mwa<maxwater)
	{
	xh2o=markw[mm1];
	markw[mm1]=maxwater;
	meltpart11(mtk,mpb,mm1,eps1);
	if(eps1[21]>0) mwa=maxwater;
	markw[mm1]=xh2o;
	}
/* No dehydration below database */
if(ynout) {mwa=markw[mm1]; if(mwa<0) mwa=0;}
}
mwa=0;
/* Cp calc by interpolation */
mcp=((H2-H0)*(1.0-n)+(H3-H1)*n)/tkstp;
if(mcp<1e+2) mcp=1e+2; else if(mcp>5e+4) mcp=5e+4;
/* Effective adiabatic betta=1/V*dV/dT=ro/T*[-dH/dP+V] calc by interpolation */
mbb=(2.0/(R1+R0)-(H1-H0)/pbstp/1e+5)*(1.0-e)+(2.0/(R3+R2)-(H3-H2)/pbstp/1e+5)*e;
mbb*=mro/mtk;
if(mbb<-1e-2) mbb=-1e-2; else if(mbb>1e-2) mbb=1e-2;
/* Effective compressibility term alpha=1/ro*d(ro)/dP calc by interpolation */
maa=(2.0/(R1+R0)*(R1-R0)*(1.0-e)+2.0/(R3+R2)*(R3-R2)*e)/pbstp/1e+5;
if(maa<0) maa=0;
/* Activation enthalpy recalc using enthalpy changes */
/* Current Enthalpy */
mhh1=((H0*(1.0-n)+H1*n)*(1.0-e)+(H2*(1.0-n)+H3*n)*e);
/* Pmin Enthalpy */
mhh0=(td[n1][0 ][mm3][1]*(1.0-e) + td[n1+1][0 ][mm3][1]*e)*1000.0*4.1837;
/* Enthalpy Difference calc */
mdhh=(mhh1-mhh0);
/*
{printf("TD1 %d %d %e %e   %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb);getchar();}
{printf("TD1 %d %d %e %e   %e %e %e %e %e %e   %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,mgg,mro,mwa,mcp,mbb,maa,mdhh,mhh1,mhh0);getchar();}
eps1[47]=mhh1;
eps1[48]=mhh0;
*/
/* Save TD variables */
eps1[41]=mro;
eps1[42]=mwa;
eps1[43]=mcp;
eps1[44]=mbb;
eps1[45]=maa;
eps1[46]=mdhh;
}
/* Thermodynamic database use for ro, Cp */




/* Rock to rock+melt transformation */
void melting1(double mtk, double mpb, long int mm1, double *eps1)
/* mtk - T, K */
/* mpb - P, bar */
/* mm1 - mark number */
{
/* Marker type */
int mm2=(int)markt[mm1];
/**/
/* Melting related cahnge of the marker type */
/* Check marker type */
if (mm2>1)
if (mpb>0)
switch (mm2)
	{
	/* All rocks except hydrated mantle */
	case 2:
	case 3:
	case 4:
	case 5:
	case 6:
	case 7:
	case 8:
	case 9:
	case 10:
	case 12:
	case 15:
	case 14:
	case 16:
	case 17:
	case 18:
	case 22:
	case 23:
	case 24:
	case 25:
	case 26:
	case 27:
	case 28:
	case 29:
	case 30:
	case 32:
	case 34:
	case 35:
	case 36:
	case 37:
	case 38:
	meltpart11(mtk,mpb,mm1,eps1);
	if(eps1[21]>0 && mm2<20) markt[mm1]+=20;
	if(eps1[21]<=0 && mm2>20) markt[mm1]-=20;
/*
if (mm2!=markt[mm1]){printf("Granite/Basalt %ld %d %d  %e %e  %e",mm1,mm2,markt[mm1],mtk,mpb,eps1[21]);getchar();}
*/
 	return;
	/**/
	/* Hydrated Peridotite */
	case 11:
	case 31:
	meltpart11(mtk,mpb,mm1,eps1);
	if(eps1[21]>0 && mm2==11) markt[mm1]=34;
	if(eps1[21]<=0 && mm2==31) markt[mm1]=14;
/*
if (mm2!=markt[mm1]){printf("Peridotite %ld %d %d  %e %e  %e",mm1,mm2,markt[mm1],mtk,mpb,eps1[21]);getchar();}
*/
 	return;
	/* Others */
	default: return;
	}
}
/* Rock to rock+melt transformation */





/* Melt fraction, density, viscosity, heat capacity calculation */
double meltpart0(double mtk, double mpb, long int mm1, int mm2, double *eps1)
/* mtk - T, K */
/* mpb - P, bar */
/* mm1 - mark number */
/* mm2 - mark type */
{
/* Val buffer */
double xmelt=0,ival,dmpb,dmtk,epsin,sduct,nueff,smin,smax,nmin,nmax,cpadd=0;
long int m1;
int mm3;
float mfd[2];
/**/
/* Check marker type */
if (mm2>=20)
	{
	/* Calculate melt fraction */
	meltpart11(mtk,mpb,mm1,eps1);
	/*
	xmelt=eps1[21]-markex[mm1];if(xmelt<0) xmelt=0;
	*/
	/**/
	/* Solid rock type */
	mm3=mm2-20;
	/**/
	/* Standard adiabatic term: al=bro/(1+bro*(Tk-298.15)) */
	eps1[20]=(markbb[mm2]*xmelt+markbb[mm3]*(1.0-xmelt))/(1.0-(markbb[mm2]*xmelt+markbb[mm3]*(1.0-xmelt))*(mtk-298.15));
	/**/
	/* Density */
	/* Ro=ro0 */
	if (densimod==0) 
		{
		eps1[23]=markro[mm2]*xmelt+markro[mm3]*(1.0-xmelt);
		}
	/* Ro=ro0*(1-bro*(TK-298.15))*(1+aro*(Pkbar-0.001)) */
	else
		{
		eps1[23]=(markro[mm2]*xmelt+markro[mm3]*(1.0-xmelt))*(1.0-(markbb[mm2]*xmelt+markbb[mm3]*(1.0-xmelt))*(mtk-298.15))*(1.0+(markaa[mm2]*xmelt+markaa[mm3]*(1.0-xmelt))*(mpb-1.0)*1e-3);
/* Eclogitization */
if((mm2==27 || mm2==28) && mpb>1.5e+4) eps1[23]+=300.0*(1.0-xmelt);
		}
	/**/
	/* Viscosity */
	/* Effective NU calc check */
	/* Little melt */
	if(xmelt<0.1)
		{
		eps1[24]=viscalc(mtk,mpb,mm1,mm3,mfd);
		}
	else
		{
		/* Significant melt */
		/* Effective NU calc check */
		nueff=marknu[mm2]*exp(2.5+pow((1.0-xmelt)/xmelt,0.48)*(1.0-xmelt));
		/* Viscosity check */
/*
		if(nueff<nubeg) nueff=nubeg; 
		if(nueff>nuend) nueff=nuend;
*/
		if(nueff<markn0[mm2]) nueff=markn0[mm2]; 
		if(nueff>markn1[mm2]) nueff=markn1[mm2];
		eps1[24]=nueff;
		}
	/**/
	/* Heat capacity */
	eps1[25]=markcp[mm2]*xmelt+markcp[mm3]*(1.0-xmelt);
	/**/
	/* heat conductivity */
	eps1[26]=((markkt[mm2]+markkf[mm2]/(mtk+77.0))*exp(markkp[mm2]*mpb))*xmelt+((markkt[mm3]+markkf[mm3]/(mtk+77.0))*exp(markkp[mm3]*mpb))*(1.0-xmelt);
	/**/
	/* Additional melting adiabatic term, heat capacity */
	if(xmelt>0 && xmelt<1.0)
		{
		/* Melting adiabatic term: alm=-ro*(dHlat/dP)/T */
		/* Numerical differentiation */
		dmpb=mpb*0.001;
		meltpart11(mtk,mpb-dmpb,mm1,eps1);
		ival=eps1[22];
		meltpart11(mtk,mpb+dmpb,mm1,eps1);
		ival-=eps1[22];
		ival*=eps1[23]/(mtk*2.0*dmpb*1e+5);
		eps1[20]+=ival;
		/**/
		/* Melting heat capacity term: cpm=dHlat/dT */
		/* Numerical differentiation */
		dmtk=1.0;
		meltpart11(mtk+dmtk,mpb,mm1,eps1);
		ival=eps1[22];
		meltpart11(mtk-dmtk,mpb,mm1,eps1);
		ival-=eps1[22];
		ival/=2.0*dmtk;
		eps1[25]+=ival;
/*
printf("Meltpart: %ld %d %e %e %e  %e %e %e %e %e %e",mm1,mm2,mtk,mpb,xmelt,eps1[20],eps1[22],eps1[23],eps1[24],eps1[25],eps1[26]);getchar();
*/
		}
	return 1.0;
	}
eps1[20]=eps1[21]=eps1[22]=eps1[23]=eps1[24]=eps1[25]=eps1[26]=0;
return 0;
}
/* Rock to rock+melt transformation */




/* Melt fraction, latent heat calculation */
void meltpart11(double mtk, double mpb, long int mm1, double *eps1)
/* mtk - T, K */
/* mpb - P, bar */
/* x,y - XY location of point for Vx,Vy calc */
/* mm1 - mark number */
/* yn  - type of calculation: 0 - Ro, 1 - Nu, 2 - Cp, 3 - kt */
{
/* Val buffer */
double xmelt=0,hlatent=0,ival,mpg,Fcpxout,Tcpxout;
long int m1;
double ykm=mpb*3e-3,ts=0,tl=0,tll=0,mwa=markw[mm1],xh2osat,xh2o,dt,dt0,dtmin,dtmax;
/* Marker type */
int n1,n2,mm2=(int)markt[mm1];
/**/
/**/
/**/
/* Calculate melt fraction using marker type */
if (ykm>0)
switch(mm2)
	{
	/* Sediments: latent heat 300 kJ/kg (Bittner & Schmeling, 1995) */
	case 15:
	case 17:
	case 35:
	case 3:
	case 4:
	case 5:
	case 23:
	case 24:
	case 25:
	case 37:
	/* Wet Solidus Temperature, Johannes, 1985, Poli & Schmidt, 2002 */
	if (ykm<36.0) 
		{
		ts=889.0+536.6/(ykm+1.609)+18.21/(ykm+1.609)/(ykm+1.609);
		}
	else
		{
		ts=831.3+2.0*ykm;
		}
	/* Dry Granite Liquidus, Johannes, 1985 */
	tl=1262.0+3.0*ykm;
	hlatent=300000.0;
	break;
	/**/
	/* Basalt, Gabbro: latent heat 380 kJ/kg (Bittner & Schmeling, 1995) */
	case 7:
	case 16:
	case 18:
	case 27:
	case 36:
	case 38:
	/* Wet solidus, Schmidt & Poli, 1998  */
	if (ykm<48.0) 
		{
		ts=972.6-2111.0/(ykm+10.63)+70033.0/(ykm+10.63)/(ykm+10.63);
		}
	else
		{
		ts=935.4+0.1162*ykm+0.006937*ykm*ykm;
		}
	/* Dry Toleitic Basalt Liquidus, Hess, 1989 */
	tl=1423.15+3.5*ykm;
	hlatent=380000.0;
	break;
	/**/
	/* Gabbro latent heat 380 kJ/kg (Bittner & Schmeling, 1995) */
	case 6:
	case 26:
	case 8:
	case 28: 
	/* Dry Toleitic Basalt solidus, Hess, 1989 */
	ts=1327.15+3.02*ykm;
	/* Dry Toleitic Basalt Liquidus, Hess, 1989 */
	tl=1423.15+3.5*ykm;
	hlatent=380000.0;
	break;
	/**/
	/* Dry/Wet Peridotite: latent heat 400 kJ/kg Turcotte & Schubert, 1982, p.171 */
	case  9:
	case 10:
	case 11:
	case 12:
	case 14:
	case 29:
	case 30:
	case 31:
	case 32:
	case 34:
	/* Dry/wet mantle melting Katz et al., 2003 */
	mpg=mpb*1e-4;
	if(mwa<=0)
		{
		dtmin=dtmax=0;
		n2=1;
		}
	else
		{
		xh2osat=12.0*pow(mpg,0.6)+1.0*mpg;
		xh2o=mwa/(0.01+0.0*(1.0-0.01));
		if(xh2o>xh2osat) xh2o=xh2osat;
		dt=43.0*pow(xh2o,0.75);
		ts=273.15+1085.7+132.9*mpg-5.1*mpg*mpg-dt;
		n2=1;
		dtmin=dtmax=dt;
		if(mtk>=ts)
			{
			xh2osat=12.0*pow(mpg,0.6)+1.0*mpg;
			xh2o=mwa/(0.01+1.0*(1.0-0.01));
			if(xh2o>xh2osat) xh2o=xh2osat;
			dt=43.0*pow(xh2o,0.75);
			dtmin=dt;
			n2=100;
			}
		}
	for(n1=0;n1<n2;n1++)
		{
		dt=(dtmin+dtmax)/2.0;
		ts=273.15+1085.7+132.9*mpg-5.1*mpg*mpg-dt;
		tll=273.15+1475.0+80.0*mpg-3.2*mpg*mpg-dt;
		tl=273.15+1780.0+45.0*mpg-2.0*mpg*mpg-dt;
		xmelt=0;
		if(mtk>=ts)
			{
			if(mtk>=tl)
				{
				xmelt=1.0;
				}
			else
				{
				/* 0.15 = 15 wt% of Cpx in the mantle */
				Fcpxout=0.15/(0.5+0.08*mpg);
				Tcpxout=pow(Fcpxout,1.0/1.5)*(tll-ts)+ts;
				if(mtk>Tcpxout)
					{
					/* Opx melting */
					xmelt=Fcpxout+(1.0-Fcpxout)*pow((mtk-Tcpxout)/(tl-Tcpxout),1.5);
					}
				else
					{
					/* Cpx melting */
					xmelt=pow((mtk-ts)/(tll-ts),1.5);
					}
				/* Compute artificial Tliquidus */
				tl=(mtk-ts)/xmelt+ts;
				}
			}
		/* Water content in melt */
		if(mwa>0)
			{
			xh2osat=12.0*pow(mpg,0.6)+1.0*mpg;
			xh2o=mwa/(0.01+xmelt*(1.0-0.01));
			if(xh2o>xh2osat) xh2o=xh2osat;
			dt0=43.0*pow(xh2o,0.75);
			if(dt0>dt) 
				{
				dtmin=dt;
				}
			else
				{
				dtmax=dt;
				}
			}
		if(ABSV(dtmax-dtmin)<0.001) break;
		}
/*
if(mwa>0){printf("%d %d %e %e %e   %e %e %e   %e %e %e   %e %e %e %e",n1,mm2,mtk-273.15,mpg,mwa,ts,tll,tl,xh2o,xh2osat,xmelt,dt,dt0,dtmin,dtmax);getchar();}
*/
	hlatent=400000.0;
	break;
	/**/
	/* Other rocks - No melting */
	default:
	break;
	}
/**/
/* Melt fraction, latent heat calculation */
eps1[21]=eps1[22]=0;
if(tl)
	{
	/* Melt fraction calc, check */
	xmelt=(mtk-ts)/(tl-ts);
	if(xmelt<0) xmelt=0;
	if(xmelt>1.0) xmelt=1.0;
	eps1[21]=xmelt;
	/* Latent heat calc */
	hlatent*=xmelt;
	eps1[22]=hlatent;
/*
if(mm2<20 && xmelt) {printf("Meltpart1: %d %e %e  %e %e",mm2,mtk,mpb,xmelt,hlatent);getchar();}
*/
	}
/**/
}
/* Melt fraction, latent heat calculation */





/* Hydration front progress after H2O budget */
double hydration2()
{
/* Val buffer */
double ysurf,vfiltr,yfiltr,dydx,dydx1,sy1,sy2,sy3,sy4,sy5,e1,mwamin,x0,y0,x1,y1,vx1,vy1;
double hytimesum,hytimesum0;
/* TD Database variables */
double n,e,dx,dy,dz,eps1[100],wi1[100];
double mtk,mpb,mwa,mro,dmwa,wro;
long int m1,m2,m3,m4,mm1,marknum1=marknum,wn1[100],markfluid=0;
int mm2,mm3,n1,n2;
/**/
printf("\n WATER Transport BEG \n");
/* Marker steps */
dx=dxwater;
dy=dywater;
dz=dzwater;
/**/
/**/
/* Min water contents in the hydraten mantle wt% */
mwamin=0.1;
/* Min Distance from erosion surface for water release */
ysurf=8000.0;
/**/
/**/
/* Clear wa[] wt */
#pragma omp parallel for shared(val1,sol0,sol1,fre0,fre1) private(m1) firstprivate(nodenum,nodenum2,nodenum3) schedule(static)
/*
*/
for (m1=0;m1<nodenum;m1++)
	{
	val1[m1]=0;
	val1[nodenum+m1]=0;
	sol0[m1]=0;
	sol1[m1]=0;
	sol0[nodenum+m1]=1e+30;
	sol1[nodenum+m1]=-1e+30;
	sol0[nodenum2+m1]=1e+30;
	sol1[nodenum2+m1]=-1e+30;
	sol0[nodenum3+m1]=1e+30;
	sol1[nodenum3+m1]=-1e+30;
	fre0[         m1]=1e+30;
	fre0[nodenum +m1]=-1e+30;
	fre0[nodenum2+m1]=1e+30;
	fre0[nodenum3+m1]=-1e+30;
	fre1[         m1]=1e+30;
	fre1[nodenum +m1]=-1e+30;
	}
/**/
/**/
/**/
/* Fluid marker generation cycle */
for (mm1=0;mm1<marknum;mm1++)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/* Count fluid markers */
if(markt[mm1]>=50 && markt[mm1]<100) markfluid++;
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m4=m3serch((double)markz[mm1]);
/* Node num */
m3=m4*xynumxy+m1*ynumy+m2;
/**/
/*Fluid disappearance surface */
sy1=waterlev;
/**/
/*
e1=(markx[mm1]-gx[m1])/(gx[m1+1]-gx[m1]);
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Check markers out of grid and within hydration range */
if(markx[mm1]>0 && marky[mm1]>0 && markz[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (double)(markz[mm1])<zsize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
if((double)(markd[mm1])>=0 && (double)(markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10)
	{
	if(mm2<50)
		{
		/* P, T parameters calc */
		prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
		mpb=eps1[10]*1e-5;
		mtk=(double)(markk[mm1]);
		/**/
		/* Mantle to Antigorite transformation */
		antigor(mtk,mpb,mm1);
/*
*/
		/**/
		/* Rocks to rock+melt transformation */
		melting1(mtk,mpb,mm1,eps1);
/*
*/
		if (markt[mm1]>=20)
			{
			/* Check melting extent */
			if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
			if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
			if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
			if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
			if(fre1[        +m3]>markz[mm1]-dz) fre1[         m3]=markz[mm1]-dz;
			if(fre1[nodenum +m3]<markz[mm1]+dz) fre1[nodenum +m3]=markz[mm1]+dz;
			}
		/* Compute TD variables */
		tdbasecalc1(mtk,mpb,mm2,mm1,eps1);
		mro=eps1[41];
		mwa=eps1[42];
		/**/
		/* Water changes in kg/m3 calc */
		dmwa=mro*(mwa-markw[mm1])*1e-2;
/*
{printf("H2O MARKER %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
{printf("H2O RELEASE %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
*/
		/**/
		/* Add water changes to the current cell, kg/m3 */
		/* Water release */
		if ((markw[mm1]-mwa)>dmwamin)
			{
			/* Save new water content */
			markw[mm1]=mwa;
			/* Generation of fluid marker (NO FLUID From melts */
			if (markt[mm1]<20 && marky[mm1]>sy1)
				{
				markt[marknum1]=markt[mm1]+50;
				markx[marknum1]=markx[mm1];
				marky[marknum1]=marky[mm1];
				markz[marknum1]=markz[mm1];
				markk[marknum1]=markk[mm1];
				markd[marknum1]=1050.0;
				markw[marknum1]=-dmwa;
				/*
				markex[marknum1]=0;
				marktm[marknum1]=0;
				markc1[marknum1]=0;
				markc2[marknum1]=0;
				*/
				/* Add aditional markers counter */
				marknum1++;
				/* Check hydration extent */
				if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
				if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
				if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
				if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
				if(sol0[nodenum3+m3]>markz[mm1]-dz) sol0[nodenum3+m3]=markz[mm1]-dz;
				if(sol1[nodenum3+m3]<markz[mm1]+dz) sol1[nodenum3+m3]=markz[mm1]+dz;
				}
			}
		else
		/* Water consuming */
			{
			if(dmwa>0)
				{
				val1[nodenum+m3]+=dmwa;
				sol1[m3]+=1.0;
				}
			}
		}
	else
	/* Fluid marker count */
		{
		/* Check position */
		if(marky[mm1]>sy1)
			{
			/* Check hydration extent */
			if(sol0[nodenum+m3]>markx[mm1]-dx) sol0[nodenum+m3]=markx[mm1]-dx;
			if(sol1[nodenum+m3]<markx[mm1]+dx) sol1[nodenum+m3]=markx[mm1]+dx;
			if(sol0[nodenum2+m3]>marky[mm1]-dy) sol0[nodenum2+m3]=marky[mm1]-dy;
			if(sol1[nodenum2+m3]<marky[mm1]+dy) sol1[nodenum2+m3]=marky[mm1]+dy;
			if(sol0[nodenum3+m3]>markz[mm1]-dz) sol0[nodenum3+m3]=markz[mm1]-dz;
			if(sol1[nodenum3+m3]<markz[mm1]+dz) sol1[nodenum3+m3]=markz[mm1]+dz;
			}
		else
		/* Erase fluid marker */
			{
			if(!xperiodic) markx[mm1]=-2.0*stp100;
			else if (!yperiodic) marky[mm1]=-2.0*stp100;
			else markz[mm1]=-2.0*stp100;
			markk[mm1]=0;
			}
		}
	}
}
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],val1[m3],val1[nodenum+m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();
}
*/
/**/
/* Rock hydration cycle */
/*
#pragma omp parallel for shared(val1,sol0,sol1,fre0,fre1,markt,markx,marky,markz,markk,markw,markex,marktm,markc1,markc2) private(m1,m2,m3,m4,mm1,mm2,mpb,mtk,mwa,mro,eps1,wn1,dmwa) firstprivate(marknum,nodenum,nodenum2,nodenum3,xynumxy,ynumy,dx,dy,dz,xsize,ysize,zsize) schedule(static)
*/
#pragma omp parallel for shared(val1,sol0,sol1,fre0,fre1,markt,markx,marky,markz,markk,markw) private(m1,m2,m3,m4,mm1,mm2,mpb,mtk,mwa,mro,eps1,wn1,dmwa) firstprivate(marknum,nodenum,nodenum2,nodenum3,xynumxy,ynumy,dx,dy,dz,xsize,ysize,zsize) schedule(static)
/*
*/
for (mm1=0;mm1<marknum;mm1++)
if(markx[mm1]>0 && marky[mm1]>0 && markz[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (double)(markz[mm1])<zsize && markt[mm1]<50)
{
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m4=m3serch((double)markz[mm1]);
/* Node num */
m3=m4*xynumxy+m1*ynumy+m2;
/* Check markers within hydration range */
if(markx[mm1]>sol0[nodenum+m3] && marky[mm1]>sol0[nodenum2+m3] && markz[mm1]>sol0[nodenum3+m3] && (double)(markx[mm1])<sol1[nodenum+m3] && (double)(marky[mm1])<sol1[nodenum2+m3] && markz[mm1]<sol1[nodenum3+m3])
	{
	/* Fluid presence mark */
/*
	markv[mm1]=-1.0;
*/
if(markt[mm1]==9 || markt[mm1]==10 || markt[mm1]==12 || markt[mm1]==14 || markt[mm1]==5 || markt[mm1]==6)
	{
	/* Mantle Hydration */
	if (markt[mm1]!=5 && markt[mm1]!=6)
		{
		mm2=markt[mm1]=11;
		}
	/* Crust  Hydration */
	else
		{
		mm2=markt[mm1]=markt[mm1]+12;
		}
	/**/
	/* P, T parameters calc */
	prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
	mpb=eps1[10]*1e-5;
	mtk=(double)(markk[mm1]);
	/**/
	/* Mantle to Antigorite transformation */
	antigor(mtk,mpb,mm1);
/*
*/
	/**/
	/* Rocks to rock+melt transformation */
	melting1(mtk,mpb,mm1,eps1);
/*
*/
	if (markt[mm1]>=20)
		{
		/* Check melting extent */
		if(fre0[        +m3]>markx[mm1]-dx) fre0[         m3]=markx[mm1]-dx;
		if(fre0[nodenum +m3]<markx[mm1]+dx) fre0[nodenum +m3]=markx[mm1]+dx;
		if(fre0[nodenum2+m3]>marky[mm1]-dy) fre0[nodenum2+m3]=marky[mm1]-dy;
		if(fre0[nodenum3+m3]<marky[mm1]+dy) fre0[nodenum3+m3]=marky[mm1]+dy;
		if(fre1[        +m3]>markz[mm1]-dz) fre1[         m3]=markz[mm1]-dz;
		if(fre1[nodenum +m3]<markz[mm1]+dz) fre1[nodenum +m3]=markz[mm1]+dz;
		}
	/**/
	/* Thermodynamic database use for Ro, Water */
	/* Compute TD variables */
	tdbasecalc1(mtk,mpb,mm2,mm1,eps1);
	mro=eps1[41];
	mwa=eps1[42];
	/**/
	/* Water changes in kg/m3 calc */
	dmwa=mro*(mwa-markw[mm1])*1e-2;
/*
{printf("H2O HYDRATE %ld %d %d %e %e   %e %e  %e %e   %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa);getchar();}
*/
	/**/
	/* Add water changes to the current cell, kg/m3 */
	/* Water consuming */
	if (dmwa>0)
		{
		val1[nodenum+m3]+=dmwa;
		sol1[m3]+=1.0;
		}
	}
	}
}
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],val1[m3],val1[nodenum+m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();
}
*/
/**/
/**/
/**/
/* Fluid marker computing cycle */
/*
#pragma omp parallel for shared(val1,sol0,sol1,fre0,fre1,markt,markx,marky,markz,markk,markw,markd,markex,marktm,markc1,markc2) private(m1,m2,m3,m4,mm1,mm2,sy1) firstprivate(marknum,nodenum,nodenum2,nodenum3,xynumxy,ynumy,xsize,ysize,zsize,stp100,waterlev,vdeep,zdeep) schedule(static)
*/
#pragma omp parallel for shared(val1,sol0,sol1,fre0,fre1,markt,markx,marky,markz,markk,markw,markd) private(m1,m2,m3,m4,mm1,mm2,sy1) firstprivate(marknum,nodenum,nodenum2,nodenum3,xynumxy,ynumy,xsize,ysize,zsize,stp100,waterlev,vdeep,zdeep) schedule(static)
/*
*/
for (mm1=0;mm1<marknum1;mm1++)
{
/**/
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
*/
/**/
/* Check markers out of grid and within hydration range */
if(markt[mm1]>=50 && markt[mm1]<100 && markx[mm1]>0 && marky[mm1]>0 && markz[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (double)(markz[mm1])<zsize)
	{
	/* Marker cell number */
	m1=m1serch((double)markx[mm1]);
	m2=m2serch((double)marky[mm1]);
	m4=m3serch((double)markz[mm1]);
	/* Node num */
	m3=m4*xynumxy+m1*ynumy+m2;
	/**/
	/*Fluid disappearance surface */
	sy1=waterlev;
	/* Water in melt region conversion */
	if(markd[mm1]<1100.0 && markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markz[mm1]>fre1[m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3] && markz[mm1]<fre1[nodenum+m3]) markd[mm1]=1150.0;
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE1 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,val1[m3],val1[nodenum+m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
	/* Check position, no fluid above erosion/sedimentation level */
/*
	if(marky[mm1]>sy1 && marky[mm1]<zdeep-(vdeep-zdeep) && (markd[mm1]<1100.0 || (markx[mm1]>fre0[m3] && marky[mm1]>fre0[nodenum2+m3] && markz[mm1]>fre1[m3] && markx[mm1]<fre0[nodenum+m3] && marky[mm1]<fre0[nodenum3+m3] && markz[mm1]<fre1[nodenum+m3])))
*/
	if(marky[mm1]>sy1 && marky[mm1]<zdeep-(vdeep-zdeep))
		{
		val1[m3]+=markw[mm1];
		sol0[m3]+=1.0;
		}
	else
	/* Erase fluid marker */
		{
		if(!xperiodic) markx[mm1]=-2.0*stp100;
		else if (!yperiodic) marky[mm1]=-2.0*stp100;
		else markz[mm1]=-2.0*stp100;
		markk[mm1]=0;
		}
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE2 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,val1[m3],val1[nodenum+m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
	}
}
/**/
/**/
/**/
/*
for(m1=0;m1<xnumx-1;m1++)
for(m2=0;m2<ynumy-1;m2++)
{
m3=m1*ynumy+m2;
printf("%ld %ld  %e %e  %e %e  %e %e   %e %e %e %e",m1,m2,gx[m1],gx[m1+1],gy[m2],gy[m2+1],val1[m3],val1[nodenum+m3],sol0[nodenum+m3],sol1[nodenum+m3],sol0[nodenum*2+m3],sol1[nodenum*2+m3]);getchar();
}
*/
/* Fluid marker consuming cycle */
#pragma omp parallel for shared(val1,sol0,sol1,fre0,fre1,markt,markx,marky,markz,markk,markw) private(m1,m2,m3,m4,mm1,mm2,mpb,mtk,mwa,mro,eps1,wn1,dmwa) firstprivate(marknum,nodenum,nodenum2,nodenum3,xynumxy,ynumy,dx,dy,dz,xsize,ysize,zsize,stp100) schedule(static)
/*
*/
for (mm1=0;mm1<marknum1;mm1++)
{
/* Marker type */
mm2=(int)markt[mm1]; if (mm2>=100) mm2-=100;
/**/
/* Marker cell number */
m1=m1serch((double)markx[mm1]);
m2=m2serch((double)marky[mm1]);
m4=m3serch((double)markz[mm1]);
/* Node num */
m3=m4*xynumxy+m1*ynumy+m2;
/**/
/*
printf("%ld %e %e %e %e %e %e",mm1,markx[mm1],marky[mm1],e1,sy1,sy2,sy3);getchar();
if(mm1>marknum){printf("%ld %d  %e %e  %e %e ",mm1,mm2,markx[mm1],marky[mm1],e1,sy1);getchar();}
*/
/**/
/* Change water consuming rocks  and fluid makers */
if(markx[mm1]>0 && marky[mm1]>0 && markz[mm1]>0 && (double)(markx[mm1])<xsize && (double)(marky[mm1])<ysize && (double)(markz[mm1])<zsize && (markk[mm1]>0 || markt[mm1]>=50) && markt[mm1]<100)
if((double)(markd[mm1])>=0 && (double)(markw[mm1])>=0 && mm2>1 && mm2!=9 && mm2!=10 && mm2!=12 && mm2!=14 && mm2!=5 && mm2!=6)
	{
	if(mm2<50)
		{
		/* P, T parameters calc */
		prcalc1((double)(markx[mm1]),(double)(marky[mm1]),(double)(markz[mm1]),eps1,wn1);
		mpb=eps1[10]*1e-5;
		mtk=(double)(markk[mm1]);
		/**/
		/* Thermodynamic database use for Ro, Water */
		/* Compute TD variables */
		tdbasecalc1(mtk,mpb,mm2,mm1,eps1);
		mwa=eps1[42];
		/**/
		/* Water changes in kg/m3 calc */
		dmwa=mwa-markw[mm1];
/*
{printf("TD! %ld %d %d %e %e %e %e ",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro);getchar();}
{printf("TDa %d %d %e %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,hpor,ppor,tpor,mwa);getchar();}
{printf("TDb %d %d %e %e %e %e %e %e",mm2,mm3,mtk-273.15,mpb/1000.0,hpor,ppor,tpor,mwa);getchar();}
*/
		/**/
		/* Add water changes to the current cell, kg/m3 */
		/* Water consuming */
		if(dmwa>0)
			{
			if (val1[nodenum+m3]<=val1[m3])
				{
				/* Save complete new water content */
				markw[mm1]=mwa;
				}
			else
				{
				/* COmpute, Save partial new water content */
				markw[mm1]=markw[mm1]+dmwa*val1[m3]/val1[nodenum+m3];
				}
/*
{printf("H2O CONSUME %ld %d %d %e %e   %e %e  %e %e   %e    %e %e",mm1,mm2,mm3,mtk-273.15,mpb/1000.0,mwa,mro,markw[mm1],markd[mm1],dmwa,val1[m3],val1[nodenum+m3]);getchar();}
*/
			}
		}
	else
	/* Fluid marker change */
		{
		if(val1[nodenum+m3]<val1[m3])
			{
			/* Count water changes for fluid marker */
			markw[mm1]*=1.0-val1[nodenum+m3]/val1[m3];
			}
		else
		/* Erase fluid marker */
			{
			if(!xperiodic) markx[mm1]=-2.0*stp100;
			else if (!yperiodic) marky[mm1]=-2.0*stp100;
			else markz[mm1]=-2.0*stp100;
			markk[mm1]=0;
			}
/*
if(markd[mm1]>1100.0) {printf("H2O REDUCE2 %ld %d   %ld %ld %ld   %e %e    %e %e %e",mm1,mm2,m1,m2,m3,val1[m3],val1[nodenum+m3],markx[mm1],markw[mm1],markd[mm1]);getchar();}
*/
		}
	}
}
/*
marknum=marknum1;
return 0;
*/
/**/
/**/
/**/
/* Reset aditional markers */
if(printmod) printf("\n WATER BEG Number of markers: OLD = %ld NEW = %ld \n",marknum,marknum1);
mm1=0;
while(marknum1>marknum && mm1<marknum)
	{
	/* Reload marker */
	if(markt[mm1]<100 && (markx[mm1]<0 || marky[mm1]<0 || markz[mm1]<0 || (double)(markx[mm1])>xsize || (double)(marky[mm1])>ysize || (double)(markz[mm1])>zsize))
		{
		/* Decrease aditional markers counter */
		marknum1--;
		if(markx[marknum1]>=0);
			{
			/* Type save */
			markt[mm1]=markt[marknum1];
			/* X,Y, water reload */
			markx[mm1]=markx[marknum1];
			marky[mm1]=marky[marknum1];
			markz[mm1]=markz[marknum1];
			markw[mm1]=markw[marknum1];
			markd[mm1]=markd[marknum1];
			markk[mm1]=markk[marknum1];
			/*
			markex[mm1]=markex[marknum1];
			marktm[mm1]=marktm[marknum1];
			markc1[mm1]=markc1[marknum1];
			markc2[mm1]=markc2[marknum1];
			*/
			}
		}
	/* Increase markers counter */
	mm1++;
	}
if(printmod) printf("\n WATER END Number of markers: OLD = %ld NEW = %ld FLUID %ld \n",marknum,marknum1,markfluid);
/* Set new marker number */
marknum=marknum1;
/**/
/**/
/**/
return 0;
}
/* Hydration front progress after H2O budget */



/* Antigorite weakening of mantle */
void antigor(double mtk, double mpb, long int mm1)
/* mtk - T, K */
/* mpb - P, bar */
/* mm1 - mark number */
{
/* Val buffer */
double k1;
int mm2=(int)markt[mm1];
/* Y below surface in meters */
double y=mpb*3.0;
/**/
/* Check marker type */
if(mm2<9 || (mm2>13 && mm2!=34)) return;
/**/
/**/
/**/
/*
printf("%d %e %e %e %e %e",m1,x,y,e,hydry,hydryl);getchar();
printf("%ld %e %e %e %e %e %e",mm1,x,y,xsubd,ysubd,vxs,vys,);getchar();
*/
/* Antigorite weakening of mantle above oceanic crust */
/*
printf("%d %e %e %e %e %e",m1,x,y,xsubd,e,hydry);getchar();
printf("%ld %e %e %e %e %e %e ",mm1,x,y,xsubd,ysubd,vxs,vys);getchar();
*/
/* Atg stability field after Schmidt and Poli, 1998 */
if(y>63000.0)
	{
	k1=1013.17699-0.060387633e-3*y-0.004289442e-6*y*y;
	}
else
	{
	k1=751.490422+6.00773668e-3*y-0.034690759e-6*y*y;
	}
/* Change marker Type */
/* Serpentinized (13) - to hydrated (11) */
if(k1<=mtk && markt[mm1]==13) markt[mm1]=11;
/* Hydrated(11) -  to serpentinized (13) */
if(k1>mtk && markt[mm1]==11) markt[mm1]=13;
}
/* Antigorite weakening of mantle */




/* Grid position change */
double gridchange()
{
/* Counters */
long int m3,mm1;
double ival=0;
/**/
/**/
/**/
/* Check mantle wedge displacement */
for (m3=0;m3<=marknum;m3++) 
/* Check markers out of grid */
if ((markt[m3]==10 || markt[m3]==12 || markt[m3]==14 || markt[m3]==34) && markx[m3]>0 && marky[m3]>0 && markz[m3]>0 && (double)(markx[m3])<xsize && (double)(marky[m3])<ysize && (double)(markz[m3])<zsize)
	{
	/* Check Asthenosphere position */
	if (markx[m3]>gx[100] && markx[m3]<gx[300] && marky[m3]<5e+4 && marky[m3]>3e+4 && markx[m3]>ival) ival=markx[m3];
	}
if (printmod) printf("\n Wedge position = %e km \n",ival/1000.0);
/**/
/* Calc leftward displacement of markers relative to the grid */
if (ival>gx[200])
	{
	/* Add current/total grid displacement */
	ival-=gx[200];
	if(ival>maxxystep) ival=maxxystep;
	gridcur+=ival;
	gridtot+=ival;
	}
else
	{
	ival=0;
	}
/**/
/* Print total grid displacement */
if (printmod) printf("\n Total grid displacement due to subduction = %e km \n",gridtot/1000.0);
/* Return grid displacement */
return ival;
}
/* End Grid position change */




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
/*
if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
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
/*
if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Vy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
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
/*
if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Vz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
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
/*
if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Exx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
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
/*
if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Exy %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPSxy for current node */
			eps1[7]+=exy[m40]*swt;
			/* Add SIGxy for current node */
			eps1[57]+=sxy[m40]*swt;
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
/*
if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Exz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
*/
			/* Add EPSxz for current node */
			eps1[8]+=exz[m40]*swt;
			/* Add SIGxz for current node */
			eps1[58]+=sxz[m40]*swt;
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
/*
if(ddx<0||ddx>1.0||ddy<0||ddy>1.0||ddz<0||ddz>1.0||swt<0||swt>1.0){printf("Eyz %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,x,X,y,Y,z,Z);getchar();}
printf("Vx %ld %ld %ld   %ld %ld %ld %ld  %e %e %e %e %e",m1,m2,m3,m10,m20,m30,m40,ddx,ddy,ddz,swt,vx[m40]);getchar();
*/
			/* Add EPSyz for current node */
			eps1[9]+=eyz[m40]*swt;
			/* Add SIGyz for current node */
			eps1[59]+=syz[m40]*swt;
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

