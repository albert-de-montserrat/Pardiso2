/* ADD MATRIX */
int gausmat2(int am, long int mcmax, long int mcmin)
/* wn[] - line koef numbers */
/* wi[] - line koef values */
/* am - mode of action rm=1,2,3 - add, rm=0 - solve */
/* val0[] - matrix contents */
/* lin0[] - line numbers for matrix contents */
/* pos0[] - first pos numbers for line in  val0[], lin0[] */
/* num0[] - pos numbers for line in  val0[], lin0[] */
/* fre0[] - free member for lines */
/* pos0cur - first free position number in val0[], lin0[] */
/* mcmax - current line in num0[] */
/* mcmin - first line in num0[] */
{
/* Counters */
long int m1,m2,m3;
/* Val Buffer */
double ival,ival1;
/**/
/* Space Check */
if (mcmax>=MAXPAR)
	{
	printf("EXIT PROGRAM: Space out in fre0[] %ld",mcmax);
	exit(0);
	}
/**/
/**/
/**/
/* STEP 0: RELOAD KOEF FROM wi[] to val0[] */
/*
printf("\n GAUS %ld %ld   %ld   %ld %e ",am,mcmax,pos0cur,wn[0],wi[0]); getchar();
*/
	/**/
	/* Free member reload from buffer */
	fre0[mcmax]=wi[0];
	/**/
	/* Line koef reload to val1[] from buffer */
	pos0[mcmax]=pos0cur;
	num0[mcmax]=0;
	for (m2=1;m2<=wn[0];m2++)
		{
		/* wi[]>0 */
		if (wi[m2])
			{
			/* Check left/right band width */
			m3=mcmax-wn[m2];
			leftnum=MAXV(leftnum,m3);
			rightnum=MAXV(rightnum,-m3);
			/* Save Cur Koef */
			lin0[pos0cur]=wn[m2];
			val1[pos0cur]=wi[m2];
			pos0cur++;
			/* Check Space */
			if (pos0cur>=MAXBUF) 
				{
				printf("EXIT PROGRAM: Space out in val1[] %d %ld",mcmax,pos0cur);
				exit(0);
				}
			num0[mcmax]++;
			}
		}
	return 0;
}
/* End STEP 0: RELOAD KOEF FROM wi[] to val0[] */




/* SOLVE MATRIX BY ECONOMICAL FRONTAL GAUSS METHOD */
int gausmat3(int am, long int mcmax, long int mcmin)
/* wn[] - line koef numbers */
/* wi[] - line koef values */
/* am - mode of action rm=1,2,3 - add, rm=0 - solve */
/* val0[] - matrix contents */
/* lin0[] - line numbers for matrix contents */
/* pos0[] - first pos numbers for line in  val0[], lin0[] */
/* num0[] - pos numbers for line in  val0[], lin0[] */
/* fre0[] - free member for lines */
/* pos0cur - first free position number in val0[], lin0[] */
/* mcmax - current line in num0[] */
/* mcmin - first line in num0[] */
{
/* Counters */
int n0,n1;
long int m1,m2,m3,m4,m6,m8,leftn,rightn,mcmax1;
/* Val Buffer */
double ival,ival1,*buf1;
/**/
/**/
/**/
/* STEP 1: EXEED KOEF ELIMINATION MATRIX BODY FORMATION am>0 */
/* Band Matrix size check */
m1=(mcmax+1)*rightnum;
printf("\n REQUIRED BAND MATRIX SIZE: val0[%ld] (LINES=%ld BAND=%ld)\n",m1,mcmax+1,rightnum);
/*
*/
/* Check Space */
if (m1>MAXMAT) 
	{
	printf("\n EXIT PROGRAM: Space out in val0[] %ld > %ld \n",m1,MAXMAT);
	exit(0);
	}
/**/
/**/
/**/
/* Parallel region ----------------------------------------------------- */
#pragma omp parallel private(n0,n1,m1,m2,m3,m4,m6,m8,leftn,rightn,mcmax1,ival,buf1) shared(num0,lin0,pos0,val0,val1,fre0,sol0,buf0)
{
/* Obtain total number of threads */
n0=omp_get_num_threads();
/* Obtain cur thread number */
n1=omp_get_thread_num();
/* Master thread print total number of threads */
if (n1==0)
	{
	printf("Number of threads = %d\n",n0);
	}
/**/
mcmax1=mcmax;
leftn=leftnum;
rightn=rightnum;
#pragma omp flush(buf0)
buf1=buf0+(mcmax*n1);
/* Clear buffer */
for (m1=0;m1<=mcmax1;m1++) 
	{
	*(buf1+m1)=0; 
	}
/* Gaussian Elimination */
/* Major line cycle */
for (m1=n1;m1<=mcmax1;m1+=n0)
	{
	#pragma omp flush(pos0,num0)
if (m1<=mcmax1 && (*(num0+m1)) && pos0[m1]!=-1) 
	{
	/* Reload koef from val1[] */
	m6=m1;
	#pragma omp flush(num0,pos0,lin0,val1)
	for (m2=0;m2<(*(num0+m1));m2++) 
		{
		/* Cur line number */
		m4=*(lin0+(*(pos0+m1))+m2);
		/* Check presence of the line */
		if (!(*(num0+m4))) 
			{
			printf(" \n THREAD %d: EXIT PROGRAM: Line  %ld absent in matrix num0=%ld when processing Line %ld \n",n1,m2,num0[m2],m1);
			exit(0);
			}
		/* Reload line */
		*(buf1+m4)+=*(val1+(*(pos0+m1))+m2);
		/* Define smallest index */
		m6=MINV(m4,m6);
		}
	/**/
	/* Gaussian elimination of coefficients in the current line */
	for (m2=m6;m2<m1;m2++)
		{
		if(*(buf1+m2))
			{
			ival=*(buf1+m2);
			*(buf1+m2)=0;
			/* Current Line koef recalc after cur koef and upper line */
			/* Wait for the line untill ready */
			#pragma omp flush(pos0)
			while(pos0[m2] != -1) 
				{
				#pragma omp flush(pos0)
				}
			/* Calc first position for curent line in val0[] */
			m8=m2*rightn;
			/* 1-st coef of any upper line = 1 */
			#pragma omp flush(val0)
			for (m3=m2+1;m3<=m2+rightn;m3++)
				{
				*(buf1+m3)-=ival*val0[m8];
				m8++;
				}
			/**/
			/* Free member recalc after cur koef upper line */
			#pragma omp flush(fre0)
			fre0[m1]-=fre0[m2]*ival;
			}
		}
	/**/
	/* Cur line normalize in val0[] */
	/* Check Singularity */
	/* Check Presence of line */
	ival=*(buf1+m1);
	*(buf1+m1)=0;
	if (!ival)
		{
		printf("THREAD %d EXIT PROGRAM: Matrix is singular at Line %ld",n1,m1);
		exit(0);
		}
	/* Calc first position for curent line in val0[] */
	m4=m1*rightn;
	#pragma omp flush(val0)
	for (m2=m1+1;m2<=m1+rightn;m2++)
		{
		/* Recalc Val reset buf  */	
		val0[m4]=*(buf1+m2)/ival;
		*(buf1+m2)=0;
		m4++;
		}
	/* Free member recalc */
	#pragma omp flush(fre0)
	fre0[m1]/=ival;
	/* Ready line signal for other threads */
	}
	#pragma omp flush(pos0)
	*(pos0+m1)=-1;
	}
/* End STEP 1: EXEED KOEF ELIMINATION MATRIX BODY FORMATION am>0 */
/* End Parallel region ----------------------------------------------------- */
/**/
/**/
/**/
/* STEP 3: SOLUTION CALC CHECK am=0 */
#pragma omp barrier
#pragma omp master 
{
#pragma omp flush(val0,fre0,sol0)
/*
if(mcmax1>40000)
for (m1=0;m1<=mcmax1;m1++) 
{printf("BEG Thread %d    %ld  %e %e \n",n1,m1,ival,*(fre0+m1));}
*/
for (m1=mcmax1;m1>=0;m1--)
	{
	/* Calc sol0[] */
	if (num0[m1])
		{
		/* Recalc koef after Sol Koef at current line is allways = 1 */
		/* Calc first position for curent line in val0[] */
		m4=m1*rightn;
		ival=fre0[m1];
		for (m2=m1+1;m2<m1+1+rightn;m2++)
			{
			ival-=val0[m4]*sol0[m2];
			m4++;
			}
		/* Calc Sol */
		sol0[m1]=ival;
		}
	else
		{
		sol0[m1]=0;
		}
	}
}
}
/* End STEP 3: SOLUTION CALC CHECK am=0 */
return 0;
}
/* End SOLVE MATRIX BY ECONOMICAL GAUSS METHOD */

