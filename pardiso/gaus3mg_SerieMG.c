/* CLEAR mat(rm) */
void matclear(int rm, double *adr)
/* rm - rank of matrix */
/* adr - adres of first element of matrix */
{
  /* String length */
  /* int rm1=rm+1; */
  /* Counters */
  int m1,m2;
  /**/
  for (m1=0;m1<rm;m1++)
  for (m2=0;m2<=rm;m2++)
  {
    *adr=0;
    adr++;
  }
}
/* End CLEAR mat0(rm) */

/* SOLVE mat(rm) BY USING GAUSS METHOD */
void gausmat(int rm, double *adr, double *adr1)
/* rm - rank of matrix */
/* adr - adres of first element of matrix */
/* adr1 - adres of first element of solutions */
{
  /* Counters */
  int m1,m2,m3;
  /* Buffer for val */
  double val1, val2;
  /* String lenth */
  int rm1=rm+1;
  /**/
  /* STEP 2: KOEF OPERATIONS */
  for (m1=0;m1<rm-1;m1++)
  if ((val1=*(adr+m1*rm1+m1)))
  for (m2=m1+1;m2<rm;m2++)
  if ((val2=*(adr+m2*rm1+m1)))
  for (m3=m1+1;m3<=rm;m3++)
  /* F[]  G[] */
  *(adr+m2*rm1+m3)=(*(adr+m2*rm1+m3))/val2-(*(adr+m1*rm1+m3))/val1;
  /* End STEP 2: KOEF OPERATIONS */
  /**/
  /* STEP 3: SOLUTION CALC CHECK */
  for (m1=rm-1;m1>=0;m1--)
  {
    /* Clear Sol */
    *(adr1+m1)=0;
    if ((val1=*(adr+m1*rm1+m1))!=0)
    {
      /* Calc Sol */
      val2=*(adr+m1*rm1+rm)/val1;
      *(adr1+m1)=val2;
      /* Use Sol */
      for (m2=0;m2<=m1-1;m2++)
      if ((val1=*(adr+m2*rm1+m1))!=0)
      *(adr+m2*rm1+rm)-=val1*val2;
    }
  }
  /* End STEP 3: SOLUTION CALC CHECK */
}
/* End SOLVE m0() BY USING GAUSS METHOD */

/* ADD/SOLVE PARDISO MATRIX */
int pardisomat(int am, long int mcmax,long int mcmax2, long int pos0cur1, int tid, long int *un1, double *ui1)
/* wn[] - line koef numbers */
/* wi[] - line koef values */
/* am - mode of action rm=1,2,3 - add, rm=0 - solve */
/* ia[] - first coloumn of each row */
/* ja[] - coloumn position of each koef in the respective row from left to right */
/* a[]  - value of each koef */
/* b[] - Right hand side */
/* bufv[] - koefficient values buffer */
/* bufv1[] - koefficient use y/n */
{
  /* Counters */
  long int m1,m2,m3,leftnum,rightnum;
  int i;
  /* Val Buffer */
  double ival;
  /**/
  /**/
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
  /* Solve Matrix ===================================================== */
  if (am==0)
  {
    /* Solve Matrix with PARDISO */
    /* Matrix data. */
    /*Calling function making arrays */
    MKL_INT  na=mcmax;
    /**/
    MKL_INT mtype = 11; /* Real unsymmetric matrix */
    /* RHS and solution vectors. */
    MKL_INT nrhs = 1; /* Number of right hand sides. */
    /* Internal solver memory pointer pt, */
    /* 32-bit: int pt[64]; 64-bit: long int pt[64] */
    /* or void *pt[64] should be OK on both architectures */
    void *pt[64];
    /* Pardiso control parameters. */
    MKL_INT iparm[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    /* Auxiliary variables. */
    MKL_INT i;
    double ddum; /* Double dummy */
    MKL_INT idum; /* Integer dummy. */
    /* -------------------------------------------------------------------- */
    /* .. Setup Pardiso control parameters. */
    /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++) {
      iparm[i] = 0;
    }
    iparm[0]  = 1;  /* No solver default */
    iparm[1]  = 2;  /* 2:Fill-in reordering from METIS  || 3: OpenMP version of the nested dissection algorithm */
    /* Numbers of processors, value of OMP_NUM_THREADS */
    int nt;
    #pragma omp parallel
    {
      nt=omp_get_num_threads();
    }
    iparm[2]  = nt;
    iparm[3]  = 0;  /* No iterative-direct algorithm */
    iparm[4]  = 0;  /* No user fill-in reducing permutation */
    iparm[5]  = 0;  /* Write solution into x */
    iparm[6]  = 0;  /* Not in use */
    iparm[7]  = 20;  /* Max numbers of iterative refinement steps */
    iparm[8]  = 0;  /* Not in use */
    iparm[9]  = 13; /* Perturb the pivot elements with 1E-13 */
    iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
    iparm[11] = 0;  /* Not in use */
    iparm[12] = 0;  /* Not in use */
    iparm[13] = 0;  /* Output: Number of perturbed pivots */
    iparm[14] = 0;  /* Not in use */
    iparm[15] = 0;  /* Not in use */
    iparm[16] = 0;  /* Not in use */
    iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
    iparm[18] = -1; /* Output: Mflops for LU factorization */
    iparm[19] = 0;  /* Output: Numbers of CG Iterations */
    //iparm[20] = 0;
    //iparm[23] = 0;  // uses a two-level factorization algorithm. This algorithm generally improves scalability in case of parallel factorization on many threads (more than eight).
    //iparm[26] = 1;  // matrix check
    //iparm[59] = 2; // 1 = in-core (if possible); 2 = out of core
    maxfct    = 1; /* Maximum number of numerical factorizations. */
    mnum      = 1; /* Which factorization to use. */
    msglvl    = 0; /* Print statistical information in file */
    error     = 0; /* Initialize error flag */
    /* -------------------------------------------------------------------- */
    /* .. Initialize the internal solver memory pointer. This is only */
    /* necessary for the FIRST call of the PARDISO solver. */
    /* -------------------------------------------------------------------- */
    for (i = 0; i < 64; i++) {
      pt[i] = 0;
    }

    clock_t t;
    t = clock();
    /* -------------------------------------------------------------------- */
    /* .. Reordering and Symbolic Factorization. This step also allocates */
    /* all memory that is necessary for the factorization. */
    /* -------------------------------------------------------------------- */
    phase = 11;
    PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
      &na, a, ia, ja, &idum, &nrhs,
      iparm, &msglvl, &ddum, &ddum, &error);
      if (error != 0) {
        printf("\nERROR during symbolic factorization: %d\n", error);
        exit(1);
      }
      printf("\nReordering completed ... ");
      printf("\nNumber of nonzeros in factors = %d", iparm[17]);
      printf("\nNumber of factorization MFLOPS = %d", iparm[18]);

      /* -------------------------------------------------------------------- */
      /* .. Numerical factorization. */
      /* -------------------------------------------------------------------- */
      phase = 22;
      PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
        &na, a, ia, ja, &idum, &nrhs,
        iparm, &msglvl, &ddum, &ddum, &error);
        if (error != 0) {
          printf("\nERROR during numerical factorization: %d", error);
          exit(2);
        }
        printf("\nFactorization completed ... ");
        /* -------------------------------------------------------------------- */
        /* .. Back substitution and iterative refinement. */
        /* -------------------------------------------------------------------- */
        phase = 33;
        iparm[7] = 2; /* Max numbers of iterative refinement steps. */
        PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
          &na, a, ia, ja, &idum, &nrhs,
          iparm, &msglvl, b, x, &error);
          if (error != 0) {
            printf("\nERROR during solution: %d", error);
            exit(3);
          }
          printf("\nSolve completed ...\n ");

          t = clock() - t;
          double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds

          printf("#### PARDISO TOOK %f seconds to execute \n", time_taken);
          printf ("\n");
          /* -------------------------------------------------------------------- */
          /* .. Termination and release of memory. */
          /* -------------------------------------------------------------------- */
          phase = -1; /* Release internal memory. */
          PARDISO (pt, &maxfct, &mnum, &mtype, &phase,
            &na, &ddum, ia, ja, &idum, &nrhs,
            iparm, &msglvl, &ddum, &ddum, &error);
            /**/
            /**/
            if (printmod) printf("PARDISO SOLVER OK!\n");
            /**/
            /**/
            /* Load the residuals */
            int mgi=multinum;
            long int mcmax1;
            system("rm sol0.txt");
            fl = fopen("sol0.txt","a");
            for (m3=0;m3<mgz[mgi];m3++)
            for (m1=0;m1<mgx[mgi];m1++)
            for (m2=0;m2<mgy[mgi];m2++)
            {
              /* Pos in vx[], vy[], vz[], pr[], etc */
              mcmax1                  = (m3*mgxy[mgi]+m1*mgy[mgi]+m2)*4;
              sol0[mgp[mgi]+mcmax1+0] = x[mcmax1+0];
              sol0[mgp[mgi]+mcmax1+1] = x[mcmax1+1];
              sol0[mgp[mgi]+mcmax1+2] = x[mcmax1+2];
              sol0[mgp[mgi]+mcmax1+3] = x[mcmax1+3];
              fprintf(fl,"%e \n", sol0[mgp[mgi]+mcmax1+0]);
              fprintf(fl,"%e \n", sol0[mgp[mgi]+mcmax1+1]);
              fprintf(fl,"%e \n", sol0[mgp[mgi]+mcmax1+2]);
              fprintf(fl,"%e \n", sol0[mgp[mgi]+mcmax1+3]);
            }
            fclose(fl);

            printf("#### PARDISO SOLUTION PRINTED \n"); getchar();
            /**/

            return 0;
          }
          /* End Solve Matrix ===================================================== */
          /**/
          /**/
          // clock_t t;
          // t = clock();

          /* R.h.s. */
          // b[mcmax]     = ui1[0];
          // fre1[mcmax2] = ui1[0];
          // if (fre1[mcmax2]!=0) {
          //   printf("%e %i",fre1[mcmax2],mcmax2); getchar();
          // }
          /* Reset coefficients */
          for (m2=1;m2<=un1[0];m2++)
          {
            /* wi[]>0 */
            if (ui1[m2])
            {
              bufvs[un1[m2]] = 0.0;
              curs[un1[m2]]  = 1;
            }
          }

          /* COMPUTE BANDWITH  */
          leftnum  = nodenum4+1;
          rightnum = -1;
          for (m2=1;m2<=un1[0];m2++){
            if (ui1[m2])
            {/* Check left/right band width || do this to build ja & a in left2right order */
              m3                   = un1[m2];
              leftnum              = MINV(leftnum,m3);
              rightnum             = MAXV(rightnum,m3);
              bufvs[un1[m2]] += ui1[m2];
              if(un1[m2]<0){printf("negative un1[m2] = %i\n",un1[m2]);}
            }
          }

          if(leftnum<0){printf("negative leftnum = %i\n",leftnum);getchar();}
          ia[mcmax+1]=0;
          for (m2=leftnum;m2<=rightnum;m2++)
          {
            /* Marked Line */
            if (curs[m2]==1) {
              if (bufvs[m2])
              {
                ja[pos0cur1] = m2+1;
                a[pos0cur1]  = bufvs[m2];
                pos0cur1++;
                ia[mcmax+1]++;
                if(ja[pos0cur1]<0) {printf("ciao %i %i\n",ja[pos0cur1], m2+1); getchar();}
              }
            }

            /* Clean marked line */
            bufvs[m2]=0.0;
            curs[m2]=0;
          }
          return 0;
        } /** END OF PARDIDO **/

        /* ADD MATRIX */
        int gausmat2(int am, long int mcmax, long int pos0cur1, long int *un1, double *ui1)
        /* un[] - line koef numbers */
        /* ui[] - line koef values */
        /* am - mode of action rm=1,2,3 - add, rm=0 - solve */
        /* val1[] - matrix contents */
        /* lin0[] - line numbers for matrix contents */
        /* pos0[] - first pos numbers for line in  val0[], lin0[] */
        /* num0[] - pos numbers for line in  val0[], lin0[] */
        /* fre0[] - free member for lines */
        /* pos0cur1 - position number in val0[], lin0[] */
        /* mcmax - current line in num0[] */
        /* mcmin - first line in num0[] */
        {
          /* Counters */
          long int m1,m2,m3;
          /* Val Buffer */
          double ival,ival1;
          /**/
          /**/
          /* STEP 0: RELOAD KOEF FROM ui[] to val1[] */
          /*
          printf("GAUS %ld   %e ",mcmax,fre1[mcmax]); getchar();
          */
          /**/
          /* Free member reload from buffer */
          fre1[mcmax]=ui1[0];
          /*
          if(un[0]>1) {printf(" %ld %ld   %ld %e \n",pos0cur,mcmax,un[0],ui[0]); getchar();}
          */
          /**/
          /* Line koef reload to val1[] from buffer */
          pos0[mcmax]=pos0cur1;
          num0[mcmax]=0;

          for (m2=1;m2<=un1[0];m2++)
          {
              /* Save Cur Koef */
              lin0[pos0cur1]=un1[m2];
              val1[pos0cur1]=ui1[m2];
              pos0cur1++;
              num0[mcmax]++;
          }

          // fl = fopen("coef_par_finest_old.txt","a");
          // for (m2=0;m2<=un1[0];m2++)
          // {
          //     fprintf(fl,"%e %i \n",ui1[m2],un1[m2]);
          // }
          // fclose(fl);
          /*
          if (mcmax>4000) getchar();
          */
          return 0;
        }
        /* End STEP 0: RELOAD KOEF FROM ui[] to val1[] */
