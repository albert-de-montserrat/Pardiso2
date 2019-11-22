/* Load information from configuration file mode.t3c ============== */
int loadconf()
{
    /* Counter */
    int m1,m2,n1,n2,fln3;
    /**/
    /**/
    /**/
    /* Open File file.t3c */
    fl = fopen("file.t3c","rt");
    ffscanf(); fln3=atoi(sa)-1;
    fclose(fl);
    /**/
    /**/
    /**/
    /* Open File mode.t3c */
    fl = fopen("mode.t3c","rt");
    /**/
    /* Output directory */
    ffscanf();
    for (n1=0;n1<500;n1++) outdir0[n1]=sa[n1];
    sprintf(outdir,"/%s",outdir0);
    /* Data File name */
    char sb[50];
    ffscanf();
    for (m1=0;m1<50;m1++) sb[m1]=sa[m1];
    ffscanf();
    if(sa[0] == 'b')
    {
        fl1itp=1;
        strcat(sb,".prn");
    }
    else if(sa[0] == 'h')
    {
        fl1itp=2;
        strcat(sb,".h5");
    }
    for (n1=0;n1<50;n1++) fl1in[n1]=sb[n1];
    /**/
    /* Load first Results File names */
    ffscanf();
    fl0num=0;
    while(sa[0]!='~')
    {
        /* Check file Counter */
        if(fl0num>=MAXFLN) {printf("Space out in fl0out[]"); exit(0);}
        /**/
        /* Save results file name */
        for (n1=0;n1<50;n1++) sb[n1]=sa[n1];
        /**/
        /* Load TYPE_cyc0max_maxxystep_maxtkstep_maxtmstep */
        ffscanf();
        if(sa[0] == 'b')
        {
            fl0otp[fl0num]=1;
            strcat(sb,".prn");
        }
        else if(sa[0] == 'h')
        {
            fl0otp[fl0num]=2;
            strcat(sb,".h5");
        }
        /* Save results file name */
        for (m1=0;m1<50;m1++) fl0out[fl0num][m1]=sb[m1];
        /**/
        /* Load TYPE_cyc0max_maxxystep_maxtkstep_maxtmstep */
        ffscanf();fl0cyc[fl0num][0]=atoi(sa);
        ffscanf();fl0stp[fl0num][0]=atof(sa);
        ffscanf();fl0stp[fl0num][1]=atof(sa);
        ffscanf();fl0stp[fl0num][2]=atof(sa)*3.15576e+7;
        ffscanf();fl0stp[fl0num][3]=atof(sa);
        ffscanf();fl0stp[fl0num][4]=atof(sa);
        ffscanf();fl0stp[fl0num][5]=atof(sa);
        ffscanf();fl0stp[fl0num][6]=atof(sa);
        ffscanf();fl0stp[fl0num][7]=atof(sa);
        ffscanf();fl0stp[fl0num][8]=atof(sa);
        ffscanf();fl0cyc[fl0num][1]=atoi(sa);
        ffscanf();fl0cyc[fl0num][2]=atoi(sa);
        ffscanf();fl0cyc[fl0num][3]=atoi(sa);
        /**/
        /* Incr File Counters */
        fl0num++;
        /**/
        /* Load Next Results File names */
        ffscanf();
    }
    /**/
    /* Service */
    ffscanf();printmod=atoi(sa);
    ffscanf();movemod=atoi(sa);
    ffscanf();tempmod=atoi(sa);
    ffscanf();markmod=atoi(sa);
    ffscanf();gridmod=atoi(sa);
    ffscanf();outgrid=atoi(sa);
    ffscanf();densimod=atoi(sa);
    ffscanf();stp100=atof(sa);
    ffscanf();filesjob=atoi(sa);
    if(stp100<=0) {printf("Invalid stp100 %e",stp100);exit(0);}
    printf("\n Current stp100 %e \n",stp100);
    /**/
    /* Effective viscosity calc */
    ffscanf();effbulkviscmod=atoi(sa);
    ffscanf();platesdistance=atof(sa);
    /**/
    /* Errosion/Sedimentation */
    ffscanf();eroslev=atof(sa);
    ffscanf();sedilev=atof(sa);
    ffscanf();waterlev=atof(sa);
    /**/
    /* V Iteration */
    ffscanf();cyc1max=atoi(sa);
    ffscanf();DIVVMIN=atof(sa);
    ffscanf();STOKSMIN=atof(sa);
    ffscanf();DIVVMAX=atof(sa);
    ffscanf();STOKSMAX=atof(sa);
    ffscanf();multimax=multinum=atoi(sa);
    ffscanf();multicyc=atoi(sa);
    for(m1=0;m1<=1;m1++)
        for(m2=0;m2<=multinum;m2++)
        {
            ffscanf();multinnn[m1][m2]=atoi(sa);
        }
    ffscanf();p0koef=atof(sa);
    ffscanf();p1koef=atof(sa);
    ffscanf();p2koef=atof(sa);
    ffscanf();v0koef=atof(sa);
    ffscanf();v1koef=atof(sa);
    ffscanf();nubeg=atof(sa);
    ffscanf();nuend=atof(sa);
    ffscanf();nukoef=atof(sa);
    ffscanf();viscmod=atoi(sa);
    ffscanf();spheryn=atoi(sa);
    ffscanf();drexmod=atoi(sa);
    /**/
    /* T Iteration */
    ffscanf();cyc2max=atoi(sa);
    ffscanf();HEATMIN=atof(sa);
    ffscanf();multinumt=atoi(sa);
    ffscanf();multicyct=atoi(sa);
    for(m1=0;m1<=1;m1++)
        for(m2=0;m2<=multinumt;m2++)
        {
            ffscanf();multittt[m1][m2]=atoi(sa);
        }
    ffscanf();t0koef=atof(sa);
    ffscanf();t1koef=atof(sa);
    ffscanf();heatdif=atof(sa);
    ffscanf();frictyn=atof(sa);
    ffscanf();adiabyn=atof(sa);
    /**/
    /* Water */
    ffscanf();tkpor=atof(sa);
    ffscanf();zmpor=atof(sa);
    ffscanf();vyfluid=atof(sa);
    ffscanf();vymelt=atof(sa);
    ffscanf();dmwamin=atof(sa);
    ffscanf();tdeep=atof(sa);
    ffscanf();dtdeep=atof(sa);
    ffscanf();drdeep=atof(sa);
    ffscanf();zdeep=atof(sa);
    ffscanf();vdeep=atof(sa);
    ffscanf();nudeep=atof(sa);
    ffscanf();dxwater=atof(sa);
    ffscanf();dywater=atof(sa);
    ffscanf();dzwater=atof(sa);
    ffscanf();maxwater=atof(sa);
    ffscanf();minmelt=atof(sa);
    ffscanf();maxmelt=atof(sa);
    /**/
    /* Velocity change */
    ffscanf();timebeg=atof(sa)*365.25*24.0*3600.0;
    ffscanf();timeend=atof(sa)*365.25*24.0*3600.0;
    ffscanf();velocitykf=atof(sa);
    /*
     printf("%e %e %e %e %e",tkpor,zmpor,vyfluid,vymelt,dmwamin);getchar();
     */
    /* Data File name change after number */
    if(fln3>=0 && fln3<fl0num)
    {
        for (n1=0;n1<50;n1++) fl1in[n1]=fl0out[fln3][n1];
        fl1itp=fl0otp[fln3];
    }
    else
    {
        fln3=-1;
    }
    /**/
    fclose(fl);
    /* End Load information from configuration file mode.t3c */
    /* Add directory to file name */
    sprintf(fl1in1,"%s%s",outdir,fl1in);
    /**/
    /* stop.yn file creation */
    fl = fopen("stop.yn","wt");
    fprintf(fl,"n \n");
    fclose(fl);
    /**/
    /*
     */
    /* Load thermodynamic database */
    /* Dry peridotite */
    /* RO - density */
    fl = fopen("m895_ro","rt");
    ffscanf();
    ffscanf(); tknum=atoi(sa);
    ffscanf(); pbnum=atoi(sa);
    ffscanf(); tkmin=atof(sa);
    ffscanf(); pbmin=atof(sa);
    ffscanf(); tkstp=atof(sa);
    ffscanf(); pbstp=atof(sa);
    ffscanf();
    ffscanf();
    /*
     printf("%d %d %e %e %e %e",tknum,pbnum,tkmin,pbmin,tkstp,pbstp);getchar();
     */
    for (n1=0;n1<pbnum;n1++)
        for (n2=0;n2<tknum;n2++)
        {
            ffscanf(); td[n2][n1][0][0]=atof(sa)/1000.0;
            /*
             printf("%d %d %e",n1,n2,td[n1][n2][0][0]);getchar();
             */
        }
    fclose(fl);
    /**/
    /* H  - enthalpy */
    fl = fopen("m895_hh","rt");
    for (n1=0;n1<9;n1++) ffscanf();
    for (n1=0;n1<pbnum;n1++)
        for (n2=0;n2<tknum;n2++)
        {
            ffscanf(); td[n2][n1][0][1]=atof(sa)/1000.0/4.1837;
            /*
             if (n2>0)
             {
             ival=1000.0*4.1837*(td[n2][n1][0][1]-td[n2-1][n1][0][1])/tkstp;
             if (ival<=1e+2 || ival>=5e+4) {printf("PDRY %d %d %e <%s>",n1,n2,ival,sa);getchar();}
             }
             */
        }
    fclose(fl);
    if (printmod) printf("m895_hh OK \n");
    /**/
    if(1==0){
        /* Wet peridotite */
        /* RO - density */
        fl = fopen("pwet_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][1][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("pwet_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][1][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][1][1]-td[n2-1][n1][1][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("PWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("pwet_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][1][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        if (printmod) printf("pwet_h2o OK \n");
        /**/
        /**/
        /* Molten peridotite */
        /* RO - density */
        fl = fopen("pwetmelt_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][2][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("pwetmelt_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][2][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][2][1]-td[n2-1][n1][2][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("PMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("pwetmelt_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][2][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        if (printmod) printf("pwetmelt_h2o OK \n");
        /**/
        /**/
        /* Wet Gabbro */
        /* RO - density */
        fl = fopen("gabwet_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][3][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("gabwet_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][3][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][3][1]-td[n2-1][n1][3][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("GWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("gabwet_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][3][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        if (printmod) printf("gabwet_h2o OK \n");
        /**/
        /**/
        /* Molten Gabbro */
        /* RO - density */
        fl = fopen("gabwetmelt_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][4][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("gabwetmelt_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][4][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][4][1]-td[n2-1][n1][4][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("GMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("gabwetmelt_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][4][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        if (printmod) printf("gabwetmelt_h2o OK \n");
        /**/
        /**/
        /* Wet sediments */
        /* RO - density */
        fl = fopen("swet_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][5][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("swet_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][5][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][5][1]-td[n2-1][n1][5][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("SWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("swet_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][5][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        if (printmod) printf("swet_h2o OK \n");
        /**/
        /**/
        /* Molten sediments */
        /* RO - density */
        fl = fopen("swetmelt_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][6][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("swetmelt_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][6][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][6][1]-td[n2-1][n1][6][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("SMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("swetmelt_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][6][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        if (printmod) printf("swetmelt_h2o OK \n");
        /**/
        /**/
        /* Wet Basalt */
        /* RO - density */
        fl = fopen("bwet_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][7][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("bwet_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][7][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][7][1]-td[n2-1][n1][7][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("BWET %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("bwet_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][7][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        if (printmod) printf("bwet_h2o OK \n");
        /**/
        /**/
        /* Molten Basalt */
        /* RO - density */
        fl = fopen("bwetmelt_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][8][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("bwetmelt_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][8][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][8][1]-td[n2-1][n1][8][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("bwetmelt_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][8][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        if (printmod) printf("bwetmelt_h2o OK \n");
        /**/
        /**/
        /* Dry Upper crust */
        /* RO - density */
        fl = fopen("ucdry_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][11][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("ucdry_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][11][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][11][1]-td[n2-1][n1][11][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        printf("ucdry_hh OK \n");
        /**/
        /* Wet Upper crust */
        /* RO - density */
        fl = fopen("ucwet_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][12][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("ucwet_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][12][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][12][1]-td[n2-1][n1][12][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("ucwet_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][12][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        printf("ucwet_h2o OK \n");
        /**/
        /**/
        /* Dry Lower crust */
        /* RO - density */
        fl = fopen("lcdry_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][13][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("lcdry_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][13][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][13][1]-td[n2-1][n1][13][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        printf("lcdry_hh OK \n");
        /**/
        /* Wet Lower crust */
        /* RO - density */
        fl = fopen("lcwet_rho","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][14][0]=atof(sa)/1000.0;
            }
        fclose(fl);
        /**/
        /* H  - enthalpy */
        fl = fopen("lcwet_hh","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][14][1]=atof(sa)/1000.0/4.1837;
                /*
                 if (n2>0)
                 {
                 ival=1000.0*4.1837*(td[n2][n1][14][1]-td[n2-1][n1][14][1])/tkstp;
                 if (ival<=1e+2 || ival>=5e+4) {printf("BMLT %d %d %e <%s>",n1,n2,ival,sa);getchar();}
                 }
                 */
            }
        fclose(fl);
        /**/
        /* Wa - water contents, wt% */
        fl = fopen("lcwet_h2o","rt");
        for (n1=0;n1<9;n1++) ffscanf();
        for (n1=0;n1<pbnum;n1++)
            for (n2=0;n2<tknum;n2++)
            {
                ffscanf(); td[n2][n1][14][2]=MAXV(0,atof(sa));
            }
        fclose(fl);
        printf("lcwet_h2o OK \n");
        /**/
    }
    /**/
    /*
     printf("%d %d %e %e %e %e ",n1,n2,td[n1-1][n2-1][0][0],td[n1-1][n2-1][0][1],td[n1-1][n2-1][0][2],td[n1-1][n2-1][0][3]);getchar();
     */
    /**/
    /* Load global thermodynamic database */
    /*
     printf("%d %d %e %e %e %e ",n1,n2,td[n1-1][n2-1][0][0],td[n1-1][n2-1][0][1],td[n1-1][n2-1][0][2],td[n1-1][n2-1][0][3]);getchar();
     */
    /**/
    /**/
    /**/
    /* Return counter */
    return fln3;
}
/* Load information from configuration file mode.t3c ============== */



/* Load Information from data file ------------------------------- */
void loader()
/* bondv[] - bondary value */
/* bondm[] - bondary mode 0=Not, -1=Value, 1,2...=LinNum+1 */
{
    /* Counter */
    int n1;
    char nn1;
    long int m1,m2,m3,m4,mv[10];
    long int mm1;
    char szint,szlong,szfloat,szdouble,szcur;
    float ival0,iv[1000];
    double ival1;
    /**/
    /**/
    /**/
    /* Load Past Results from data file-------------------------------- */
    if (printmod) printf("Load Past results from %s ...",fl1in1);
    /**/
    /**/
    /**/
    /* Load in Binary Format ---------------------------- */
    if (fl1itp == 1)
    {
        fl = fopen(fl1in1,"rb");
        /**/
        /* Sizes of var definition */
        szint=sizeof(n1);
        szlong=sizeof(m1);
        szfloat=sizeof(ival0);
        szdouble=sizeof(ival1);
        /* Check sizes of variables */
        fread(&szcur,1,1,fl);
        if (szcur!=szint) {printf("Current INT size <%d> is different from given in file <%d> \n",szint,szcur); exit(0);}
        fread(&szcur,1,1,fl);
        if (szcur!=szlong) {printf("Current LONG INT size <%d> is different from given in file <%d> \n",szlong,szcur); exit(0);}
        fread(&szcur,1,1,fl);
        if (szcur!=szfloat) {printf("Current FLOAT size <%d> is different from given in file <%d> \n",szfloat,szcur); exit(0);}
        fread(&szcur,1,1,fl);
        if (szcur!=szdouble) {printf("Current DOUBLE size <%d> is different from given in file <%d> \n",szdouble,szcur); exit(0);}
        /**/
        /* Grid Parameters */
        fread(&xnumx,szlong,1,fl);
        fread(&ynumy,szlong,1,fl);
        fread(&znumz,szlong,1,fl);
        fread(&mnumx,szlong,1,fl);
        fread(&mnumy,szlong,1,fl);
        fread(&mnumz,szlong,1,fl);
        fread(&xsize,szdouble,1,fl);
        fread(&ysize,szdouble,1,fl);
        fread(&zsize,szdouble,1,fl);
        fread(&pxinit,szlong,1,fl);
        fread(&pyinit,szlong,1,fl);
        fread(&pzinit,szlong,1,fl);
        fread(&pinit,szdouble,1,fl);
        fread(&GXKOEF,szdouble,1,fl);
        fread(&GYKOEF,szdouble,1,fl);
        fread(&GZKOEF,szdouble,1,fl);
        fread(&rocknum,szint,1,fl);
        fread(&bondnum,szlong,1,fl);
        fread(&marknum,szlong,1,fl);
        fread(&n1,szint,1,fl);
        fread(&timesum,szdouble,1,fl);timesum*=3.15576e+7;
        fread(&gridcur,szdouble,1,fl);
        fread(&gridtot,szdouble,1,fl);
        /**/
        /* Dynamically allocate memory */
        dynmemall(0);
        /* Calc,Check Grid parameters */
        gridcheck();
        /**/
        /* Dynamically allocate memory */
        dynmemall(1);
        /* Rock Types information */
        fread(markn0,szdouble,rocknum,fl);
        fread(markn1,szdouble,rocknum,fl);
        fread(marks0,szdouble,rocknum,fl);
        fread(marks1,szdouble,rocknum,fl);
        fread(marknu,szdouble,rocknum,fl);
        fread(markdh,szdouble,rocknum,fl);
        fread(markdv,szdouble,rocknum,fl);
        fread(markss,szdouble,rocknum,fl);
        fread(markmm,szdouble,rocknum,fl);
        fread(markll,szdouble,rocknum,fl);
        fread(marka0,szdouble,rocknum,fl);
        fread(marka1,szdouble,rocknum,fl);
        fread(markb0,szdouble,rocknum,fl);
        fread(markb1,szdouble,rocknum,fl);
        fread(marke0,szdouble,rocknum,fl);
        fread(marke1,szdouble,rocknum,fl);
        fread(markro,szdouble,rocknum,fl);
        fread(markbb,szdouble,rocknum,fl);
        fread(markaa,szdouble,rocknum,fl);
        fread(markcp,szdouble,rocknum,fl);
        fread(markkt,szdouble,rocknum,fl);
        fread(markkf,szdouble,rocknum,fl);
        fread(markkp,szdouble,rocknum,fl);
        fread(markht,szdouble,rocknum,fl);
        /**/
        /* Nodes information */
        /* Vx,Vy,ro[],nu[],tk[],cp[],kt[],ht[] */
        for (m1=0;m1<nodenum;m1++)
        {
            fread(iv,szfloat,11,fl);
            pr[m1]=(double)(iv[0]);
            vx[m1]=(double)(iv[1]);
            vy[m1]=(double)(iv[2]);
            vz[m1]=(double)(iv[3]);
            ro[m1]=(double)(iv[4]);
            nu[m1]=(double)(iv[5]);
            tk[m1]=(double)(iv[6]);
            cp[m1]=(double)(iv[7]);
            et[m1]=(double)(iv[8]);
            kt[m1]=(double)(iv[9]);
            ht[m1]=(double)(iv[10]);
        }
        /**/
        /* Gridlines positions */
        fread(iv,szfloat,xnumx,fl);
        for (m1=0;m1<xnumx;m1++)
        {
            gx[m1]=(double)(iv[m1]);
        }
        fread(iv,szfloat,ynumy,fl);
        for (m2=0;m2<ynumy;m2++)
        {
            gy[m2]=(double)(iv[m2]);
        }
        fread(iv,szfloat,znumz,fl);
        for (m3=0;m3<znumz;m3++)
        {
            gz[m3]=(double)(iv[m3]);
        }
        /**/
        /* Bondary Conditions Equations */
        for (m1=0;m1<bondnum;m1++)
        {
            fread(mv,szlong,2,fl);
            m2=mv[0];m3=mv[1];
            bondm[m2]=m3;
            /* Check boundary array */
            if(m3>MAXBON) {printf("Space out in bondv[]"); exit(0);}
            fread(iv,szfloat,2,fl);
            bondv1[m3][0]=(double)(iv[0]);
            bondv1[m3][1]=(double)(iv[1]);
            fread(&m2,szlong,1,fl);bondn1[m3]=m2;
            /*
             printf("Bond %ld %ld %ld %e %ld %e",bondnum,m1,m3,bondv1[m3][0],bondn1[m3],bondv1[m3][1]);getchar();
             */
        }
        /*
         printf("Bond %ld %ld %ld ",bondnum,m3,m2);getchar();
         */
        /**/
        /* Markers X,Y,types */
        for (mm1=0;mm1<marknum;mm1++)
        {
            /* General information load */
            fread(iv,szfloat,7,fl);
            /*
             fread(iv,szfloat,10,fl);
             */
            markx[mm1]=iv[0];
            marky[mm1]=iv[1];
            markz[mm1]=iv[2];
            markk[mm1]=iv[3];
            markw[mm1]=iv[4];
            markd[mm1]=iv[5];
            marke[mm1]=iv[6];
            /*
             markex[mm1]=iv[6];
             marktm[mm1]=iv[7];
             markc1[mm1]=iv[8];
             markc2[mm1]=iv[9];
             */
            fread(&nn1,1,1,fl);markt[mm1]=(int)nn1;
            /*
             printf("MARK %ld   %ld %e %e %e %d ",marknum,mm1,markx[mm1],marky[mm1],markz[mm1],markt[mm1]);getchar();
             */
        }
        fclose(fl);
    }
    /* Load in Binary Format ---------------------------- */
    else if (fl1itp == 2)
    {
        hid_t       file_id, group_id, attr_id, dataset_id, dataspace_id, memtype, memspace, filetype, space,npoints;  /* identifiers */
        herr_t      status;
        hsize_t     dims,Dims,ndims[2],count[2],offset[2];
        int i,n0;
        double a[100];
        /* Open an existing file */
        file_id = H5Fopen(fl1in1, H5F_ACC_RDONLY, H5P_DEFAULT);
        /* General */
        attr_id=H5Aopen_by_name(file_id,".","General",H5P_DEFAULT,H5P_DEFAULT);
        H5Aread (attr_id, H5T_NATIVE_DOUBLE, &a);
        H5Aclose(attr_id);
        /**/
        xsize=a[0];
        ysize=a[1];
        zsize=a[2];
        pxinit=a[3];
        pyinit=a[4];
        pzinit=a[5];
        pinit=a[6];
        GXKOEF=a[7];
        GYKOEF=a[8];
        GZKOEF=a[9];
        timesum=a[10]; timesum*=3.15576e+7;
        gridcur=a[11];
        gridtot=a[12];
        n1=(int)(a[13]);
        /* Read material */
        group_id = H5Gopen(file_id, "/Material", H5P_DEFAULT);
        /**/
        /* Reading attribute */
        attr_id=H5Aopen_by_name(group_id,".","/Material/Number",H5P_DEFAULT,H5P_DEFAULT);
        H5Aread (attr_id, H5T_NATIVE_INT, &rocknum);
        H5Aclose(attr_id);
        /**/
        dims=3;
        int Periodic[3];
        attr_id=H5Aopen_by_name(file_id,".","XYZ-periodic",H5P_DEFAULT,H5P_DEFAULT);
        H5Aread (attr_id, H5T_NATIVE_INT, &Periodic);
        H5Aclose(attr_id);
        xperiodic=Periodic[0];
        yperiodic=Periodic[1];
        zperiodic=Periodic[2];
        /**/
        /* Material properties */
        dataset_id = H5Dopen(group_id, "/Material/Properties", H5P_DEFAULT);
        space=H5Dget_space(dataset_id);
        npoints=H5Sget_simple_extent_npoints(space);
        double rdata[rocknum][npoints/rocknum];
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,rdata);
        for (i=0;i<rocknum;i++)
        {
            markn0[i]=rdata[i][0];
            markn1[i]=rdata[i][1];
            marks0[i]=rdata[i][2];
            marks1[i]=rdata[i][3];
            marknu[i]=rdata[i][4];
            markdh[i]=rdata[i][5];
            markdv[i]=rdata[i][6];
            markss[i]=rdata[i][7];
            markmm[i]=rdata[i][8];
            markll[i]=rdata[i][9];
            marka0[i]=rdata[i][10];
            marka1[i]=rdata[i][11];
            markb0[i]=rdata[i][12];
            markb1[i]=rdata[i][13];
            marke0[i]=rdata[i][14];
            marke1[i]=rdata[i][15];
            markro[i]=rdata[i][16];
            markbb[i]=rdata[i][17];
            markaa[i]=rdata[i][18];
            markcp[i]=rdata[i][19];
            markkt[i]=rdata[i][20];
            markkf[i]=rdata[i][21];
            markkp[i]=rdata[i][22];
            markht[i]=rdata[i][23];
        }
        H5Dclose(dataset_id);
        H5Sclose(space);
        H5Gclose(group_id);
        /**/
        /**/
        /* Read Nodes */
        group_id = H5Gopen(file_id, "/Nodes", H5P_DEFAULT);
        /* Read /Nodes attribute */
        attr_id=H5Aopen_by_name(group_id,".","/Nodes/Number",H5P_DEFAULT,H5P_DEFAULT);
        dims=3;
        long int b[dims];
        H5Aread (attr_id, H5T_NATIVE_LONG, &b);
        H5Aclose(attr_id);
        xnumx=b[0];
        ynumy=b[1];
        znumz=b[2];
        /* Dynamically allocate memory */
        dynmemall(0);
        /* Gridlines positions */
        dims=xnumx;
        loadsave_dataset_double(group_id,dims,gx, "/Nodes/gx",1);
        dims=ynumy;
        loadsave_dataset_double(group_id,dims,gy, "/Nodes/gy",1);
        dims=znumz;
        loadsave_dataset_double(group_id,dims,gz, "/Nodes/gz",1);
        /**/
        gridcheck();
        /**/
        /* Dynamically allocate memory */
        dynmemall(1);
        /* Nodes parameters */
        dims=nodenum;
        loadsave_dataset_double(group_id,dims,pr, "/Nodes/pr",1);
        loadsave_dataset_float(group_id,dims,ro, "/Nodes/ro",1);
        loadsave_dataset_float(group_id,dims,nu, "/Nodes/nu",1);
        loadsave_dataset_float(group_id,dims,nxx, "/Nodes/nxx",1);
        loadsave_dataset_float(group_id,dims,nxy, "/Nodes/nxy",1);
        loadsave_dataset_float(group_id,dims,nxz, "/Nodes/nxz",1);
        loadsave_dataset_float(group_id,dims,nyz, "/Nodes/nyz",1);
        //loadsave_dataset_float(group_id,dims,tk, "/Nodes/tk",1);
        //loadsave_dataset_float(group_id,dims,cp, "/Nodes/cp",1);
        //loadsave_dataset_float(group_id,dims,et, "/Nodes/et",1);
        //loadsave_dataset_float(group_id,dims,kt, "/Nodes/kt",1);
        //loadsave_dataset_float(group_id,dims,ht, "/Nodes/ht",1);
        //loadsave_dataset_float(group_id,dims,fd, "/Nodes/fd",1);
        /**/
        //for(m3=0;m3<znumz;m3++)
        //for(m1=0;m1<xnumx;m1++)
        //for(m2=0;m2<ynumy;m2++)
        //{m4=m3*mgxy[0]+m1*mgy[0]+m2;
        //if(tk[m4]==0){printf(" %ld %e %e %e %e\n",m4,gx[m1],gy[m2],gz[m3],tk[m4]);getchar();}}
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
        /*
         for(m1=0;m1<nodenum;m1++){
         printf("%ld %e %e %e\n",m1,vx[m1],vy[m1],vz[m1]); getchar();
         }
         */
        /**/
        /**/
        /* Boundary conditions */
        /* Reading attribute */
        long int c[2];
        attr_id=H5Aopen_by_name(group_id,".","/Nodes/BC number",H5P_DEFAULT,H5P_DEFAULT);
        H5Aread (attr_id, H5T_NATIVE_LONG, &c);
        H5Aclose(attr_id);
        bondnum=c[0];
        dims=c[1];
        /* BC value */
        dataset_id = H5Dopen(group_id, "/Nodes/bondv", H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datad);
        H5Dclose(dataset_id);
        /* BC index */
        dataset_id = H5Dopen(group_id, "/Nodes/bondn", H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datal);
        H5Dclose(dataset_id);
        /**/
        for (m1=0;m1<dims;m1++)
        {
            m2=datal[m1][0];
            m3=datal[m1][1];
            bondm[m2]=m3;
            bondn1[m3]=datal[m1][2];
            /* Check boundary array */
            if(m3>MAXBON) {printf("Space out in bondv[]"); exit(0);}
            bondv1[m3][0]=(double)(datad[m1][0]);
            bondv1[m3][1]=(double)(datad[m1][1]);
            /*
             if(m2==4892680){printf("%ld %d %e %e\n",m1,bondm[m2],datad[m1][0],datad[m1][1]); getchar();}
             */
        }
        /**/
        H5Gclose(group_id);
        /**/
        /**/
        /* Lagrangian particles */
        group_id = H5Gopen(file_id, "/Particles", H5P_DEFAULT);
        /* Creating attribute */
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
        /**/
        //loadsave_dataset_float(group_id,dims,markk,"/Particles/markk",1);
        //loadsave_dataset_float(group_id,dims,markw,"/Particles/markw",1);
        //loadsave_dataset_float(group_id,dims,markd,"/Particles/markd",1);
        loadsave_dataset_float(group_id,dims,marke,"/Particles/marke",1);
        /* Marker type */
        loadsave_dataset_char(group_id,dims, markt,"/Particles/markt",1);
        /**/
        H5Gclose(group_id);
        /**/
        /**/
        /* Close the file */
        H5Fclose(file_id);
        /**/
    }
    /**/
    /**/
    /**/
    /* bondv1, bondn1 */
    memtot[1]+=bondnum*(2*4+2);
    if(printmod)printf("\n OUTPUT FILE ABOUT %f GB \n",(float)(memtot[1]/1e+9));
    /**/
    if(printmod)printf("\n OK!\n");
    /**/
    /**/
    /* Calc,Check Grid parameters */
    gridcheck();
    /**/
}
/* Load Information from data file ------------------------------- */

/* Load/Save datasets for in h5 format */
void loadsave_dataset_char(hid_t group_id, hsize_t ndims, char buf[],const char *name, int mode)
{
    hid_t dataspace_id, dataset_id;
    if (mode==0)
    {
        dataspace_id=H5Screate_simple(1,&ndims,NULL);
        dataset_id = H5Dcreate(group_id,name, H5T_NATIVE_CHAR, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    else
    {
        dataset_id = H5Dopen(group_id,name, H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
    }
}

void loadsave_dataset_int(hid_t group_id, hsize_t ndims, int *buf,const char *name, int mode)
{
    hid_t dataspace_id, dataset_id;
    if (mode==0)
    {
        dataspace_id=H5Screate_simple(1,&ndims,NULL);
        dataset_id = H5Dcreate(group_id,name, H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    else
    {
        dataset_id = H5Dopen(group_id,name, H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
    }
}

void loadsave_dataset_long(hid_t group_id, hsize_t ndims, long int *buf,const char *name, int mode)
{
    hid_t dataspace_id, dataset_id;
    if (mode==0)
    {
        dataspace_id=H5Screate_simple(1,&ndims,NULL);
        dataset_id = H5Dcreate(group_id,name, H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    else
    {
        dataset_id = H5Dopen(group_id,name, H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
    }
}

void loadsave_dataset_float(hid_t group_id, hsize_t ndims, float *buf,const char *name, int mode)
{
    hid_t dataspace_id, dataset_id;
    if (mode==0)
    {
        dataspace_id=H5Screate_simple(1,&ndims,NULL);
        dataset_id = H5Dcreate(group_id,name, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    else
    {
        dataset_id = H5Dopen(group_id,name, H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
    }
}

void loadsave_dataset_double(hid_t group_id, hsize_t ndims, double *buf,const char *name, int mode)
{
    hid_t dataspace_id, dataset_id;
    if (mode==0)
    {
        dataspace_id=H5Screate_simple(1,&ndims,NULL);
        dataset_id = H5Dcreate(group_id,name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    else
    {
        dataset_id = H5Dopen(group_id,name, H5P_DEFAULT);
        H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,buf);
        H5Dclose(dataset_id);
    }
}
/* Load/Save datasets for in h5 format */


/* Print Results to data file ----------------------------------- */
void saver(int f0, int n0)
/* n0 - circle number */
/* f0 - file number */
{
    /* Counters */
    int n1;
    char nn1;
    long int m1,m2,m3,m4;
    long int mm1;
    /* Buffers for XY */
    double x,y,z;
    char szint,szlong,szfloat,szdouble,szcur;
    float ival0;
    double ival1;
    float iv[1000];
    /**/
    /**/
    /**/
    if (printmod) printf("Print %d circle results to %s...",n0+1,fl1out);
    /**/
    /**/
    /**/
    /* Save data in binary format ---------------------------- */
    if (fl1otp == 1)
    {
        fl = fopen(fl1out,"wb");
        /**/
        /* Sizes of var definition */
        szint=sizeof(n1);
        szlong=sizeof(m1);
        szfloat=sizeof(ival0);
        szdouble=sizeof(ival1);
        fwrite(&szint,1,1,fl);
        fwrite(&szlong,1,1,fl);
        fwrite(&szfloat,1,1,fl);
        fwrite(&szdouble,1,1,fl);
        /**/
        /* Grid Parameters */
        fwrite(&xnumx,szlong,1,fl);
        fwrite(&ynumy,szlong,1,fl);
        fwrite(&znumz,szlong,1,fl);
        fwrite(&mnumx,szlong,1,fl);
        fwrite(&mnumy,szlong,1,fl);
        fwrite(&mnumz,szlong,1,fl);
        fwrite(&xsize,szdouble,1,fl);
        fwrite(&ysize,szdouble,1,fl);
        fwrite(&zsize,szdouble,1,fl);
        fwrite(&pxinit,szlong,1,fl);
        fwrite(&pyinit,szlong,1,fl);
        fwrite(&pzinit,szlong,1,fl);
        fwrite(&pinit,szdouble,1,fl);
        fwrite(&GXKOEF,szdouble,1,fl);
        fwrite(&GYKOEF,szdouble,1,fl);
        fwrite(&GZKOEF,szdouble,1,fl);
        fwrite(&rocknum,szint,1,fl);
        fwrite(&bondnum,szlong,1,fl);
        fwrite(&marknum,szlong,1,fl);
        fwrite(&n0,szint,1,fl);
        ival1=timesum/3.15576e+7;fwrite(&ival1,szdouble,1,fl);
        fwrite(&gridcur,szdouble,1,fl);
        fwrite(&gridtot,szdouble,1,fl);
        /**/
        /* Rock Types information */
        fwrite(markn0,szdouble,rocknum,fl);
        fwrite(markn1,szdouble,rocknum,fl);
        fwrite(marks0,szdouble,rocknum,fl);
        fwrite(marks1,szdouble,rocknum,fl);
        fwrite(marknu,szdouble,rocknum,fl);
        fwrite(markdh,szdouble,rocknum,fl);
        fwrite(markdv,szdouble,rocknum,fl);
        fwrite(markss,szdouble,rocknum,fl);
        fwrite(markmm,szdouble,rocknum,fl);
        fwrite(markll,szdouble,rocknum,fl);
        fwrite(marka0,szdouble,rocknum,fl);
        fwrite(marka1,szdouble,rocknum,fl);
        fwrite(markb0,szdouble,rocknum,fl);
        fwrite(markb1,szdouble,rocknum,fl);
        fwrite(marke0,szdouble,rocknum,fl);
        fwrite(marke1,szdouble,rocknum,fl);
        fwrite(markro,szdouble,rocknum,fl);
        fwrite(markbb,szdouble,rocknum,fl);
        fwrite(markaa,szdouble,rocknum,fl);
        fwrite(markcp,szdouble,rocknum,fl);
        fwrite(markkt,szdouble,rocknum,fl);
        fwrite(markkf,szdouble,rocknum,fl);
        fwrite(markkp,szdouble,rocknum,fl);
        fwrite(markht,szdouble,rocknum,fl);
        /**/
        /* Nodes information */
        /* Vx,Vy,ro[],nu[],tk[],cp[],kt[],ht[] */
        for (m1=0;m1<nodenum;m1++)
        {
            iv[0]=(float)(pr[m1]);
            iv[1]=(float)(vx[m1]);
            iv[2]=(float)(vy[m1]);
            iv[3]=(float)(vz[m1]);
            iv[4]=(float)(ro[m1]);
            iv[5]=(float)(nu[m1]);
            iv[6]=(float)(tk[m1]);
            iv[7]=(float)(cp[m1]);
            iv[8]=(float)(et[m1]);
            iv[9]=(float)(kt[m1]);
            iv[10]=(float)(ht[m1]);
            fwrite(iv,szfloat,11,fl);
        }
        /**/
        /* Gridlines positions */
        for (m1=0;m1<xnumx;m1++)
        {
            iv[m1]=(float)(gx[m1]);
        }
        fwrite(iv,szfloat,xnumx,fl);
        for (m2=0;m2<ynumy;m2++)
        {
            iv[m2]=(float)(gy[m2]);
        }
        fwrite(iv,szfloat,ynumy,fl);
        for (m3=0;m3<znumz;m3++)
        {
            iv[m3]=(float)(gz[m3]);
        }
        fwrite(iv,szfloat,znumz,fl);
        /**/
        /* Bondary Conditions Equations */
        m3=0;
        for (m1=0;m1<nodenum*10;m1++)
            if(bondm[m1])
            {
                fwrite(&m1,szlong,1,fl);
                m2=bondm[m1];fwrite(&m2,szlong,1,fl);
                iv[0]=(float)(bondv1[m2][0]);
                iv[1]=(float)(bondv1[m2][1]);
                fwrite(iv,szfloat,2,fl);
                m2=bondn1[m2];fwrite(&m2,szlong,1,fl);
                m3++;
            }
        printf("BOND %ld %ld ",bondnum,m3);
        /**/
        /* Markers X,Y,types */
        for (mm1=0;mm1<marknum;mm1++)
        {
            /* General information Save */
            iv[0]=markx[mm1];
            iv[1]=marky[mm1];
            iv[2]=markz[mm1];
            iv[3]=markk[mm1];
            iv[4]=markw[mm1];
            iv[5]=markd[mm1];
            iv[6]=marke[mm1];
            /*
             iv[7]=marktm[mm1];
             iv[8]=markc1[mm1];
             iv[9]=markc2[mm1];
             fwrite(iv,szfloat,10,fl);
             */
            fwrite(iv,szfloat,7,fl);
            nn1=(char)markt[mm1];fwrite(&nn1,1,1,fl);
            /*
             printf("MARK %ld   %ld %e %e %e %d ",marknum,mm1,markx[mm1],marky[mm1],markz[mm1],markt[mm1]);getchar();
             */
        }
        printf("MARK %ld ",marknum);
        fclose(fl);
    }
    /* Save data in binary format ---------------------------- */
    else if (fl1otp == 2)
    {
        hid_t       file_id, group_id, attr_id, dataset_id, dataspace_id, dcpl,memspace,plist_id;
        herr_t      status;
        hsize_t     dims,dims2,ndims[2],count[2],offset[2],chunk_dims[2],stride[2],block[2];
        int n1;
        /**/
        /* General */
        file_id = H5Fcreate(fl1out,H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        dims=14;
        double datadb[dims];
        datadb[0]=xsize;
        datadb[1]=ysize;
        datadb[2]=zsize;
        datadb[3]=pxinit;
        datadb[4]=pyinit;
        datadb[5]=pzinit;
        datadb[6]=pinit;
        datadb[7]=GXKOEF;
        datadb[8]=GYKOEF;
        datadb[9]=GZKOEF;
        datadb[10]=timesum/3.15576e+7;
        datadb[11]=gridcur;
        datadb[12]=gridtot;
        datadb[13]=(double)(n0+1);
        dataspace_id=H5Screate_simple(1,&dims,NULL);
        attr_id = H5Acreate(file_id,"General",H5T_NATIVE_DOUBLE,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(attr_id,H5T_NATIVE_DOUBLE,&datadb);
        H5Aclose(attr_id);
        H5Sclose(dataspace_id);
        /**/
        dims=3;
        int Periodic[dims];
        Periodic[0]=xperiodic;
        Periodic[1]=yperiodic;
        Periodic[2]=zperiodic;
        dataspace_id=H5Screate_simple(1,&dims,NULL);
        attr_id = H5Acreate(file_id,"XYZ-periodic",H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(attr_id,H5T_NATIVE_INT,&Periodic);
        H5Aclose(attr_id);
        H5Sclose(dataspace_id);
        /**/
        /* Material informations */
        group_id = H5Gcreate(file_id, "/Material", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        /* Creating attribute */
        dims=1;
        int datai[dims];
        datai[0]=rocknum;
        dataspace_id=H5Screate_simple(1,&dims,NULL);
        attr_id = H5Acreate(group_id,"/Material/Number",H5T_NATIVE_INT,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(attr_id,H5T_NATIVE_INT,&datai);
        H5Aclose(attr_id);
        H5Sclose(dataspace_id);
        /* Rock properties */
        ndims[0]=rocknum;
        ndims[1]=24;
        double datart[ndims[0]][ndims[1]];
        for (n1=0;n1<rocknum;n1++)
        {
            datart[n1][0]=markn0[n1];
            datart[n1][1]=markn1[n1];
            datart[n1][2]=marks0[n1];
            datart[n1][3]=marks1[n1];
            datart[n1][4]=marknu[n1];
            datart[n1][5]=markdh[n1];
            datart[n1][6]=markdv[n1];
            datart[n1][7]=markss[n1];
            datart[n1][8]=markmm[n1];
            datart[n1][9]=markll[n1];
            datart[n1][10]=marka0[n1];
            datart[n1][11]=marka1[n1];
            datart[n1][12]=markb0[n1];
            datart[n1][13]=markb1[n1];
            datart[n1][14]=marke0[n1];
            datart[n1][15]=marke1[n1];
            datart[n1][16]=markro[n1];
            datart[n1][17]=markbb[n1];
            datart[n1][18]=markaa[n1];
            datart[n1][19]=markcp[n1];
            datart[n1][20]=markkt[n1];
            datart[n1][21]=markkf[n1];
            datart[n1][22]=markkp[n1];
            datart[n1][23]=markht[n1];
        }
        dataspace_id=H5Screate_simple(2,ndims,NULL);
        dataset_id = H5Dcreate(group_id, "/Material/Properties", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datart);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        H5Gclose(group_id);
        /**/
        /**/
        /* Eulerian nodes */
        group_id = H5Gcreate(file_id, "/Nodes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        /* Creating attribute */
        dims=3;
        long int datali[dims];
        datali[0]=xnumx;
        datali[1]=ynumy;
        datali[2]=znumz;
        dataspace_id=H5Screate_simple(1,&dims,NULL);
        attr_id = H5Acreate(group_id,"/Nodes/Number",H5T_NATIVE_LONG,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(attr_id,H5T_NATIVE_LONG,&datali);
        H5Aclose(attr_id);
        H5Sclose(dataspace_id);
        /**/
        /* Nodes parameters */
        dims=nodenum;
        loadsave_dataset_double(group_id,dims,pr, "/Nodes/pr",0);
        /*
         loadsave_dataset_float(group_id,dims,vx, "/Nodes/vx",0);
         loadsave_dataset_float(group_id,dims,vy, "/Nodes/vy",0);
         loadsave_dataset_float(group_id,dims,vz, "/Nodes/vz",0);
         */
        loadsave_dataset_float(group_id,dims,ro, "/Nodes/ro",0);
        loadsave_dataset_float(group_id,dims,nu, "/Nodes/nu",0);
        loadsave_dataset_float(group_id,dims,nxx, "/Nodes/nxx",0);
        loadsave_dataset_float(group_id,dims,nxy, "/Nodes/nxy",0);
        loadsave_dataset_float(group_id,dims,nxz, "/Nodes/nxz",0);
        loadsave_dataset_float(group_id,dims,nyz, "/Nodes/nyz",0);
        //loadsave_dataset_float(group_id,dims,tk, "/Nodes/tk",0);
        //loadsave_dataset_float(group_id,dims,cp, "/Nodes/cp",0);
        //loadsave_dataset_float(group_id,dims,et, "/Nodes/et",0);
        //loadsave_dataset_float(group_id,dims,kt, "/Nodes/kt",0);
        //loadsave_dataset_float(group_id,dims,ht, "/Nodes/ht",0);
        //loadsave_dataset_float(group_id,dims,fd, "/Nodes/fd",0);
        /**/
        //for(m3=0;m3<znumz;m3++)
        //for(m1=0;m1<xnumx;m1++)
        //for(m2=0;m2<ynumy;m2++)
        //{m4=m3*mgxy[0]+m1*mgy[0]+m2;
        //if(tk[m4]<=1){printf(" %ld %e %e %e %e\n",m4,gx[m1],gy[m2],gz[m3],tk[m4]);getchar();}}
        /* Save marker coordinates as a 3 * marknum matrix */
        ndims[1]=3;
        ndims[0]=dims;
        dataspace_id=H5Screate_simple(2,ndims,NULL);
        dataset_id = H5Dcreate(group_id, "/Nodes/Velocity", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Sclose(dataspace_id);
        count[1] =1;
        count[0] =dims;
        memspace = H5Screate_simple(2, count, NULL);
        offset[0] = 0;
        offset[1] = 0;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vx);
        H5Sclose(dataspace_id);
        offset[1] = 1;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vy);
        H5Sclose(dataspace_id);
        offset[1] = 2;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        /*
         */
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, vz);
        H5Sclose(dataspace_id);
        H5Dclose(dataset_id);
        H5Sclose(memspace);
        /**/
        /* Gridlines positions */
        dims=xnumx;
        loadsave_dataset_double(group_id,dims,gx, "/Nodes/gx",0);
        dims=ynumy;
        loadsave_dataset_double(group_id,dims,gy, "/Nodes/gy",0);
        dims=znumz;
        loadsave_dataset_double(group_id,dims,gz, "/Nodes/gz",0);
        /**/
        /**/
        /* Boundary conditions */
        dims=0;
        for (m1=0;m1<nodenum8+mgp[multinum]+mgn[multinum];m1++)
            if(bondm[m1])
            {
                dims++;
            }
        m3=0;
        for (m1=0;m1<nodenum8+mgp[multinum]+mgn[multinum];m1++)
            if(bondm[m1])
            {
                datal[m3][0]=m1;
                m2=bondm[m1];
                datal[m3][1]=m2;
                datal[m3][2]=bondn1[m2];
                datad[m3][0]=bondv1[m2][0];
                datad[m3][1]=bondv1[m2][1];
                m3++;
            }
        /* BC value */
        dims2=2;
        long int c[dims2];
        c[0]=bondnum;
        c[1]=dims;
        dataspace_id=H5Screate_simple(1,&dims2,NULL);
        attr_id = H5Acreate(group_id,"/Nodes/BC number",H5T_NATIVE_LONG,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(attr_id,H5T_NATIVE_LONG,&c);
        H5Aclose(attr_id);
        H5Sclose(dataspace_id);
        /**/
        ndims[0]=dims;
        ndims[1]=3;
        dataspace_id=H5Screate_simple(2,ndims,NULL);
        dataset_id = H5Dcreate(group_id, "/Nodes/bondn", H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datal);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        /**/
        ndims[0]=dims;
        ndims[1]=2;
        dataspace_id=H5Screate_simple(2,ndims,NULL);
        dataset_id = H5Dcreate(group_id, "/Nodes/bondv", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &datad);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
        /**/
        /**/
        H5Gclose(group_id);
        /**/
        /**/
        /* Lagrangian particles */
        group_id = H5Gcreate(file_id, "/Particles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        /* Creating attribute */
        dims=4;
        long int mark[dims];
        mark[0]=marknum;
        mark[1]=mnumx;
        mark[2]=mnumy;
        mark[3]=mnumz;
        dataspace_id=H5Screate_simple(1,&dims,NULL);
        attr_id = H5Acreate(group_id,"/Particles/Number",H5T_NATIVE_LONG,dataspace_id,H5P_DEFAULT,H5P_DEFAULT);
        H5Awrite(attr_id,H5T_NATIVE_LONG,&mark);
        H5Aclose(attr_id);
        H5Sclose(dataspace_id);
        /**/
        /* Particle parameters */
        dims=marknum;
        /* Save marker coordinates as a 3 * marknum matrix */
        ndims[1]=3;
        ndims[0]=dims;
        dataspace_id=H5Screate_simple(2,ndims,NULL);
        dataset_id = H5Dcreate(group_id, "/Particles/Position", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
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
        offset[1] = 1;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, marky);
        H5Sclose(dataspace_id);
        offset[1] = 2;
        dataspace_id = H5Dget_space(dataset_id);
        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, offset, NULL, count, NULL);
        H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, memspace, dataspace_id, H5P_DEFAULT, markz);
        H5Sclose(dataspace_id);
        H5Sclose(memspace);
        H5Dclose(dataset_id);
        /**/
        //loadsave_dataset_float(group_id,dims,markk,"/Particles/markk",0);
        //loadsave_dataset_float(group_id,dims,markw,"/Particles/markw",0);
        //loadsave_dataset_float(group_id,dims,markd,"/Particles/markd",0);
        loadsave_dataset_float(group_id,dims,marke,"/Particles/marke",0);
        /* Marker type */
        loadsave_dataset_char(group_id,dims,markt,"/Particles/markt",0);
        /**/
        H5Gclose(group_id);
        /**/
        /**/
        /* Close hdf5 file */
        H5Fclose(file_id);
        /**/
        /*
         makeXMF(f0);
         */
        char command[600];
        sprintf(command,"mv %s %s",fl1out,outdir);
        system(command);
        printf("\n Move file %s to %s\n",fl1out,outdir);
    }
    /**/
    /**/
    if (printmod) printf("OK!\n");
    /**/
    /* file.t3c file creation */
    fl = fopen("file.t3c","wt");
    fprintf(fl,"%d \n",f0);
    fclose(fl);
    /**/
    /**/
    /*
     for (mm1=0;mm1<marknum;mm1++)
     {
     printf("SAVE MARK %ld   %ld %e %e %e %d ",marknum,mm1,markx[mm1],marky[mm1],markz[mm1],markt[mm1]);getchar();
     }
     */
    /* stop.yn file information read */
    fl = fopen("stop.yn","rt");
    /* Read String */
    ffscanf();
    /**/
    /* Stop Y/N */
    if (sa[0]=='y' || sa[0]=='Y')
    {
        fclose(fl);
        printf("PROGRAM TERMINATED FROM stop.yn \n");
        exit(0);
    }
    /**/
    /* Change printmod */
    if (sa[0]>='0' && sa[0]<='9')
    {
        printmod=atoi(sa);
    }
    fclose(fl);
}
/* Print Results to data file ----------------------------------- */

/* Make XMF file to visualize .h5 fields in Paraview -------------*/
void makeXMF(int f0)
{
    int m1,m2,m3,m4;
    /**/
    /* Make mesh strcture file when run in3mg */
    if(f0==0)
    {
        hid_t       file_id, group_id, attr_id, dataset_id, dataspace_id, dcpl;
        herr_t      status;
        hsize_t     dims,ndims[2];
        /**/
        /* General */
        file_id = H5Fcreate("Mesh.h5",H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
        /* Cell vertexes */
        long int xy,y;
        y=ynumy;
        xy=xnumx*ynumy;
        /*
         float vertices[nodenum][3];
         for (m3=0;m3<znumz;m3++)
         for (m1=0;m1<xnumx;m1++)
         for (m2=0;m2<ynumy;m2++)
         {
         m4=m2+y*m1+xy*m3;
         vertices[m4][0]=(float)(gx[m1]);
         vertices[m4][1]=(float)(ysize-gy[m2]);
         vertices[m4][1]=(float)(gy[m2]);
         vertices[m4][2]=(float)(gz[m3]);
         }
         ndims[1]=3;
         ndims[0]=nodenum;
         dataspace_id=H5Screate_simple(2,ndims,NULL);
         dataset_id = H5Dcreate(file_id, "vertices", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &vertices);
         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);
         */
        /* Connectivity */
        long int xyc,yc;
        yc=ynumy1;
        xyc=xnumx1*ynumy1;
        /*
         long int connectivity[cellnum][8];
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
         dataset_id = H5Dcreate(file_id, "connectivity", H5T_NATIVE_LONG, dataspace_id, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
         H5Dwrite(dataset_id, H5T_NATIVE_LONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &connectivity);
         H5Dclose(dataset_id);
         H5Sclose(dataspace_id);
         */
        /**/
        /* Close hdf5 file */
        H5Fclose(file_id);
        /**/
    }
    /* Form xmf file name */
    char fname[50];
    char c[3];
    if(f0<10)
    {
        strcpy(fname,"XDMF.00");
    }
    if(f0>=10 && f0<100)
    {
        strcpy(fname,"XDMF.0");
    }
    if(f0>=100 && f0<1000)
    {
        strcpy(fname,"XDMF.");
    }
    sprintf(c,"%d",f0);
    strcat(fname,c);
    strcat(fname,".xmf");
    printf("\nGenerate file %s\n",fname);
    /* Open xmf file */
    fl = fopen(fname,"wt");
    /* Write file */
    fprintf(fl,"<?xml version=\"1.0\" ?> \n");
    fprintf(fl,"<Xdmf xmlns:xi=\"http://www.w3.org/2001/XInclude\" Version=\"2.0\">\n\n");
    fprintf(fl," <Domain>\n\n");
    /* Eulerian nodes informations */
    fprintf(fl,"    <Grid Name=\"Eulerian Grid\">\n\n");
    fprintf(fl,"       <Time Value=\"%f\" />\n\n",timesum);
    fprintf(fl,"         <Topology Type=\"Hexahedron\" NumberOfElements=\"%ld\"> \n",cellnum);
    fprintf(fl,"            <DataItem Format=\"HDF\" DataType=\"Int\"  Dimensions=\"%ld 8\">Mesh.h5:/connectivity</DataItem> \n",cellnum);
    fprintf(fl,"         </Topology> \n\n");
    fprintf(fl,"         <Geometry Type=\"XYZ\"> \n");
    fprintf(fl,"            <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">Mesh.h5:/vertices</DataItem> \n",nodenum);
    fprintf(fl,"         </Geometry>\n\n");
    /* Add Eulerian fields  */
    add_eulerian_to_XDMF("Density",fl1out,"/Nodes/ro",0,nodenum,fl);
    add_eulerian_to_XDMF("Temperature",fl1out,"/Nodes/tk",0,nodenum,fl);
    add_eulerian_to_XDMF("Viscosity",fl1out,"/Nodes/nu",0,nodenum,fl);
    add_eulerian_to_XDMF("Velocity",fl1out,"/Nodes/Velocity",1,nodenum,fl);
    /**/
    fprintf(fl,"   </Grid>\n\n");
    /**/
    if(1==1){
        /* Lagrangian particles informations */
        fprintf(fl,"   <Grid Name=\"Lagrangian particles\">\n\n");
        fprintf(fl,"      <Time Value=\"%f\" />\n\n",timesum);
        fprintf(fl,"         <Topology Type=\"POLYVERTEX\" NodesPerElement=\"%ld\"> </Topology>\n",marknum);
        fprintf(fl,"         <Geometry Type=\"XYZ\">\n");
        fprintf(fl,"            <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%ld 3\">%s:/Particles/Position</DataItem>\n",marknum,fl1out);
        fprintf(fl,"         </Geometry>\n\n");
        /* Add Lagrangian fields  */
        add_lagrangian_to_XDMF("Marker Type",fl1out,"/Particles/markt",0,marknum,1,fl);
        /**/
        fprintf(fl,"    </Grid>\n\n");
    }
    fprintf(fl," </Domain>\n");
    fprintf(fl,"</Xdmf>\n");
    /* Close xmf file */
    fclose(fl);
}
/* End Make XMF file to visualize .h5 fields in Paraview -------------*/

/* Add Eulerian field ------------------------------------------------*/
void add_eulerian_to_XDMF(char *field_name,char *file_name,char *dataset_id,int is_vector,long int dimensions,FILE *fl)
{
    int d;
    char *type_vec_scal;
    /**/
    if(is_vector)
    {
        d=3;
        type_vec_scal="Vector";
    }
    else
    {
        d=1;
        type_vec_scal="Scalar";
    }
    /**/
    fprintf(fl,"        <Attribute Type=\"%s\" Center=\"Node\" Name=\"%s\">\n",type_vec_scal,field_name);
    fprintf(fl,"            <DataItem ItemType=\"HyperSlab\" Dimensions=\"%ld %d\">\n",nodenum,d);
    fprintf(fl,"               <DataItem Dimensions=\"3 2\" Format=\"XML\"> 0 0 1 1 %ld %d </DataItem>\n",nodenum,d);
    fprintf(fl,"               <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" Dimensions=\"%ld %d\"> %s:%s</DataItem>\n",nodenum,d,file_name,dataset_id);
    fprintf(fl,"             </DataItem>\n");
    fprintf(fl,"   	    </Attribute>\n\n");
}
/* End Add Eulerian field ------------------------------------------------*/

/* Add Lagrangian field ------------------------------------------------*/
void add_lagrangian_to_XDMF(char *field_name,char *file_name, char *dataset_id,int is_vector,long int dimensions,int is_int,FILE *fl)
{
    char *type,*type_vec_scal;
    int d;
    /**/
    if(is_vector)
    {
        d=3;
        type_vec_scal="Vector";
    }
    else
    {
        d=1;
        type_vec_scal="Scalar";
    }
    if(is_int) type="\"Int\"";
    else type="\"Float\" Precision=\"8\"";
    /**/
    fprintf(fl,"         <Attribute Type=\"%s\" Center=\"Node\" Name=\"%s\">\n",type_vec_scal,field_name);
    fprintf(fl,"           <DataItem Format=\"HDF\" NumberType=%s Dimensions=\"%ld %d\">%s:%s</DataItem>\n",type,dimensions,d,file_name,dataset_id);
    fprintf(fl,"         </Attribute>\n\n");
}
/* End Add Lagrangian field ------------------------------------------------*/


/* LOAD WITHOUT EMPTY LINES =================================== */
void ffscanf()
{
    /* Counter */
    int m1;
    /**/
    /* Read cycle */
    do
    {
        /* Input string */
        m1=fscanf(fl,"%s",sa);
        /* Check end of file */
        if (m1<1)
        {
            printf("\n Unexpected end of file\n");
            fclose(fl);
            exit(0);
        }
        /* Delete last symbol <32 */
        for(m1=strlen(sa)-1;m1>=0;m1--)
            if (*(sa+m1)<=32)
                *(sa+m1)=0;
            else
                break;
    }
    while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES =================================== */


/* LOAD WITHOUT EMPTY LINES from fl1 =================================== */
void ffscanf1()
{
    /* Counter */
    int n1;
    /**/
    /* Read cycle */
    do
    {
        /* Input string */
        n1=fscanf(fl1,"%s",sa);
        /* Check end of file */
        if (n1<1)
        {
            printf("\n Unexpected end of file\n");
            fclose(fl1);
            exit(0);
        }
        /* Delete last symbol <32 */
        for(n1=strlen(sa)-1;n1>=0;n1--)
            if (*(sa+n1)<=32)
                *(sa+n1)=0;
            else
                break;
    }
    while (strlen(sa)<1 || sa[0]=='/');
}
/* End LOAD WITHOUT EMPTY LINES from fl1 =================================== */


/* Calc,Check Parameters of Grid */
void gridcheck()
{
    /* Counter */
    int n1,n2;
    /* Nodes Num */
    nodenum=xnumx*ynumy*znumz;
    /**/
    /* Cells Num */
    cellnum=(xnumx-1)*(ynumy-1)*(znumz-1);
    /**/
    /* Rock types Num */
    if (rocknum>MAXTMR){printf("Space out in marknu[]"); exit(0);}
    /**/
    /* Koef for processing */
    xstpx=xsize/(double)(xnumx-1);
    ystpy=ysize/(double)(ynumy-1);
    zstpz=zsize/(double)(znumz-1);
    kfx=1.0/xstpx;
    kfy=1.0/ystpy;
    kfz=1.0/zstpz;
    kfxx=kfx*kfx;
    kfyy=kfy*kfy;
    kfxy=kfx*kfy;
    kfxz=kfx*kfz;
    kfyz=kfy*kfz;
    kfzz=kfz*kfz;
    /**/
    /* Spec counters */
    /* Node num in one xy slide */
    xynumxy=xnumx*ynumy;
    /* Cell num in one xy slide */
    xynumxy1=(xnumx-1)*(ynumy-1);
    /* Cell num in x,y,z directions */
    xnumx1=xnumx-1;
    ynumy1=ynumy-1;
    znumz1=znumz-1;
    /* Counters for sol[] */
    nodenum2=nodenum*2;
    nodenum3=nodenum*3;
    nodenum4=nodenum*4;
    nodenum8=nodenum*8;
    /* Multigrid variables */
    mgn[0]=nodenum;
    mgp[0]=0;
    mgs[0]=1;
    mgx[0]=xnumx;
    mgy[0]=ynumy;
    mgz[0]=znumz;
    mgxy[0]=xnumx*ynumy;
    /*
     printf("%ld %ld %ld %d",xnumx,ynumy,znumz,multinum);getchar();
     */
  //   printf("\n %ld %d %i\n",mgx[0],mgp[0],mgn[0] );
    for(n1=1;n1<=multinum;n1++)
    {
        mgp[n1]=mgp[n1-1]+mgn[n1-1];
        mgs[n1]=mgs[n1-1]*2;
        mgx[n1]=((mgx[n1-1]-5)/2)+5;
        mgy[n1]=((mgy[n1-1]-5)/2)+5;
        mgz[n1]=((mgz[n1-1]-5)/2)+5;
        mgxy[n1]=mgx[n1]*mgy[n1];
        mgn[n1]=mgx[n1]*mgy[n1]*mgz[n1];
        //printf("%ld %ld %ld %ld\n",mgx[n1],mgy[n1],mgz[n1],mgn[n1]); getchar();
        //printf("%ld %d %ld\n",mgp[n1],mgp[n1-1],mgn[n1-1] );
    }
    //printf("%i %i %i %i \n", mgn[multinum], mgx[multinum], mgy[multinum], mgz[multinum]);exit(0);
    /* Multigrid nodes positions */
    double mgdx,mgdy,mgdz;
    for(n1=0;n1<=multinum;n1++)
    {
        /* Set external X grid nodes */
        mggx[n1][0]=gx[0];
        mggx[n1][1]=gx[1];
        mggx[n1][2]=gx[2];
        mggx[n1][mgx[n1]-1]=gx[xnumx-1];
        mggx[n1][mgx[n1]-2]=gx[xnumx-2];
        mggx[n1][mgx[n1]-3]=gx[xnumx-3];
        /* Set internal X grid nodes */
        // for(n2=3;n2<mgx[n1]-3;n2++)
        // {
        //     mggx[n1][n2]=gx[(n2-2)*mgs[n1]+2];
        // }

        // if (n1!=multinum){
          for(n2=3;n2<mgx[n1]-3;n2++)
          {
              mggx[n1][n2]=gx[(n2-2)*mgs[n1]+2];
          }
        // }else{
        //   mggx[n1][0]=0;
        //   for(n2=1;n2<mgx[n1];n2++)
        //   {
        //       mggx[n1][n2]=mggx[n1][n2-1] + xsize/(double)(mgx[n1]-1);
        //   }
        // }
        // if (n1==multinum){
        //   system("rm GX.txt");
        //   fl = fopen("GX.txt","at");
        //   for(n2=0;n2<mgx[n1];n2++){
        //     fprintf(fl,"%e\n",mggx[n1][n2]);
        //   }
        //   fclose(fl);
        // }

        /* Set external Y grid nodes */
        mggy[n1][0]=gy[0];
        mggy[n1][1]=gy[1];
        mggy[n1][2]=gy[2];
        mggy[n1][mgy[n1]-1]=gy[ynumy-1];
        mggy[n1][mgy[n1]-2]=gy[ynumy-2];
        mggy[n1][mgy[n1]-3]=gy[ynumy-3];
        /* Set internal Y grid nodes */
        // for(n2=3;n2<mgy[n1]-3;n2++)
        // {
        //     mggy[n1][n2]=gy[(n2-2)*mgs[n1]+2];
        // }

        // if(n1!=multinum){
          for(n2=3;n2<mgy[n1]-3;n2++)
          {
              mggy[n1][n2]=gy[(n2-2)*mgs[n1]+2];
          }
        // }else{
        //   mggy[n1][0]=0;
        //   for(n2=1;n2<mgx[n1];n2++)
        //   {
        //       mggy[n1][n2]=mggy[n1][n2-1]+ysize/(double)(mgy[n1]-1);
        //   }
        // }

        // if (n1==multinum){
        //   system("rm GY.txt");
        //   fl = fopen("GY.txt","at");
        //   for(n2=0;n2<mgy[n1];n2++){
        //     // printf("GX %d %d %e %e\n",n1,n2,mggx[n1][n2],gx[n2]);getchar();
        //     fprintf(fl,"%e %i\n",mggy[n1][n2],n1);
        //   }
        //   fclose(fl);
        // }

         // for(n2=0;n2<mgy[n1];n2++)printf("GY %d %d %e \n",n1,n2,mggy[n1][n2]);getchar();

        /* Set external Z grid nodes */
        mggz[n1][0]=gz[0];
        mggz[n1][1]=gz[1];
        mggz[n1][2]=gz[2];
        mggz[n1][mgz[n1]-1]=gz[znumz-1];
        mggz[n1][mgz[n1]-2]=gz[znumz-2];
        mggz[n1][mgz[n1]-3]=gz[znumz-3];
        /* Set internal Z grid nodes */
        // for(n2=3;n2<mgz[n1]-3;n2++)
        // {
        //     mggz[n1][n2]=gz[(n2-2)*mgs[n1]+2];
        // }

        // if(n1!=multinum){
          for(n2=3;n2<mgz[n1]-3;n2++)
          {
              mggz[n1][n2]=gz[(n2-2)*mgs[n1]+2];
          }
        // }else{
        //   mggz[n1][0]=0;
        //   for(n2=1;n2<mgz[n1];n2++)
        //   {
        //       mggz[n1][n2]=mggz[n1][n2-1]+zsize/(double)(mgz[n1]-1);
        //   }
        // }


        // if (n1==multinum){
        //   system("rm GZ.txt");
        //   fl = fopen("GZ.txt","at");
        //   for(n2=0;n2<mgy[n1];n2++){
        //     // printf("GX %d %d %e %e\n",n1,n2,mggx[n1][n2],gx[n2]);getchar();
        //     fprintf(fl,"%e %i\n",mggz[n1][n2],n1);
        //   }
        //   fclose(fl);
        // }

        /*
         for(n2=0;n2<mgz[n1];n2++)printf("GZ %d %d %e \n",n1,n2,mggz[n1][n2]);getchar();
         */
    }
}
/* Calc,Check Parameters of Grid */

/* Dynamically allocate/free memory */
void dynmemall(int mode)
/*
 mode = 0,1,2: dynamically allocate memory
 mode = 3: free memory
 */
{
    int nt,i;
    long int mcmax,mcmax1,pos0cur1;
#pragma omp parallel
    {nt=omp_get_num_threads();}
    /**/
    /**/
    memtot[0]=0;
    if(mode==0)
    {
        gx=malloc(sizeof(double) * xnumx);
        gy=malloc(sizeof(double) * ynumy);
        gz=malloc(sizeof(double) * znumz);
        /* 2D arrays */
        mggx=malloc(sizeof(double *) * 10);
        mggy=malloc(sizeof(double *) * 10);
        mggz=malloc(sizeof(double *) * 10);
        for(i=0;i<10;i++)
        {
            mggx[i]=malloc(sizeof(double) * xnumx);
            mggy[i]=malloc(sizeof(double) * ynumy);
            mggz[i]=malloc(sizeof(double) * znumz);
        }
        memtot[0]+=(xnumx+ynumy+znumz)*11*8;
    }
    /**/
    /**/
    if(mode==1)
    {
        memtot[1]=memtot[0];
        /* Buffers */
        mcmax=mgp[multinum]+mgn[multinum];
        pos0cur1=mcmax*57;
        mcmax1=mcmax*4;
        printf("\n mcmax=%ld mcmax1=%ld Pos0cur1=%ld\n",mcmax,mcmax1,pos0cur1) ;
        sol0=malloc(sizeof(double) * mcmax1);
        sol1=malloc(sizeof(double) * mcmax1);
        val1=malloc(sizeof(double) * pos0cur1);
        fre0=malloc(sizeof(double) * mcmax1);
        fre1=malloc(sizeof(double) * mcmax1);
        bufv=malloc(sizeof(double) * mcmax1);
        memtot[0]+=8*(5*mcmax1+pos0cur1);
        bufn=malloc(sizeof(int) * mcmax1);
        lin0=malloc(sizeof(int) * pos0cur1);
        num0=malloc(sizeof(int) * mcmax1);
        pos0=malloc(sizeof(int) * mcmax1);
        bondm=malloc(sizeof(int) * (mcmax+nodenum8));
        memtot[0]+=2*(3*mcmax1+pos0cur1+mcmax+nodenum8);
        /*
         */
        bufv0 = malloc(sizeof(double) * mgn[multinum]*4 );
        cur0  = malloc(sizeof(int) * mgn[multinum]*4 );
        ia=malloc(sizeof(int) * (mgn[multinum]*4+1));
        x=malloc(sizeof(double) * mgn[multinum]*4);
        b=malloc(sizeof(double) * mgn[multinum]*4);
        ja=malloc(sizeof(int) * mgn[multinum]*57);
        a=malloc(sizeof(double) * mgn[multinum]*57);
        /*
         */
        /* Nodes */
        pr =malloc(sizeof(double) * nodenum);
        pr0=malloc(sizeof(double) * nodenum);
        nu =malloc(sizeof(float) * nodenum*2);
        nxx=malloc(sizeof(float) * nodenum*2);
        nxy=malloc(sizeof(float) * nodenum*2);
        nxz=malloc(sizeof(float) * nodenum*2);
        nyz=malloc(sizeof(float) * nodenum*2);
        ro =malloc(sizeof(float) * nodenum*2);
        nd =malloc(sizeof(float) * nodenum*2);
        vx =malloc(sizeof(float) * nodenum);
        vy =malloc(sizeof(float) * nodenum);
        vz =malloc(sizeof(float) * nodenum);
        vx0=calloc(nodenum,sizeof(float));
        vy0=calloc(nodenum,sizeof(float));
        vz0=calloc(nodenum,sizeof(float));
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
        memtot[0]+=(4*32+2*8)*nodenum;
        memtot[1]+=4*13*nodenum;
        /* Threads restriction arrays */
        Nu    = malloc(sizeof(float *) * nt);
        Nxx   = malloc(sizeof(float *) * nt);
        Nxy   = malloc(sizeof(float *) * nt);
        Nxz   = malloc(sizeof(float *) * nt);
        Nyz   = malloc(sizeof(float *) * nt);
        Ro    = malloc(sizeof(float *) * nt);
        Sol0  = malloc(sizeof(double *) * nt);
        Sol1  = malloc(sizeof(double *) * nt);
        bufv1 = malloc(sizeof(double *) * nt );
        cur1  = malloc(sizeof(int *) * nt );
        for(i=0;i<nt;i++)
        {
            Nu[i]   = malloc(sizeof(float) * mcmax1);
            Nxx[i]  = malloc(sizeof(float) * nodenum*2);
            Nxy[i]  = malloc(sizeof(float) * nodenum*2);
            Nxz[i]  = malloc(sizeof(float) * nodenum*2);
            Nyz[i]  = malloc(sizeof(float) * nodenum*2);
            Ro[i]   = malloc(sizeof(float) * mcmax1);
            Sol0[i] = malloc(sizeof(double) * mcmax1);
            Sol1[i] = malloc(sizeof(double) * mcmax1);
            bufv1[i]= malloc(sizeof(double) * mgn[multinum]*4 );
            cur1[i] = malloc(sizeof(int) * mgn[multinum]*4 );
        }
        memtot[0]+=nt*(nodenum*4*8+mcmax1*(4*2+8*2));
    }
    printf(" %i \n",mgn[multinum]*4);

    if(mode==2)
    {
        /* Markers */
        marknum0 =(long int)(marknum*MRKNEW);
        printf("\n Number of markers = %ld, Number of markers with extra ones = %ld\n",marknum,marknum0);
        markx=malloc(sizeof(float) * marknum0);
        marky=malloc(sizeof(float) * marknum0);
        markz=malloc(sizeof(float) * marknum0);
        marke=malloc(sizeof(float) * marknum0);
        markt=malloc(sizeof(char) * marknum0);
        memtot[0]+=marknum0*(4*4+1);
        memtot[1]+=marknum0*(4*4+1);
        /* Rock properties*/
        memtot[0]+=rocknum*24*8;
        memtot[1]+=rocknum*24*8;
        /* bondv1, bondn1, datal, datad */
        memtot[0]+=MAXBON*(2*4+2+3*4+2*4);
        /* TD database */
        memtot[0]+=4*351*351*15*3;
        /**/
        if(printmod)printf("\n ALLOCATE ABOUT %f GB of memory with %d procs\n",(float)(memtot[0]*1.05/1e+9),nt);
    }
    /**/
    /* Free memory */
    if(mode==3)
    {

        free(sol0);free(sol1);free(val1);free(fre0);free(fre1);
        free(bufv);free(lin0);free(num0);free(pos0);free(bufn);
        free(bondm);
        free(gx);free(gy);free(gz);
        free(pr);free(pr0);
        free(nu);free(ro);free(nd);
        free(nxx);free(nxy);free(nxz);free(nyz);
        free(vx);free(vy);free(vz);
        free(exx);free(eyy);free(ezz);free(exy);free(exz);free(eyz);
        free(sxx);free(syy);free(szz);free(sxy);free(sxz);free(syz);
        /**/
        free(ia);free(b);free(x);free(ja);free(a);
        /**/
        /* 2D arrays */
        for(i=0;i<nt;i++)
        {
            free(Nu[i]);
            free(Nxx[i]);
            free(Nxy[i]);
            free(Nxz[i]);
            free(Nyz[i]);
            free(Ro[i]);
            free(Sol0[i]);
            free(Sol1[i]);
            free(bufv1[i]);
            free(cur1[i]);
        }

        free(Nxx);free(Nxy);free(Nxz);free(Nyz);
        free(Nu);free(Ro);
        free(Sol0);free(Sol1);
        free(bufv1);free(cur1);
        for(i=0;i<10;i++)
        {
            free(mggx[i]);
            free(mggy[i]);
            free(mggz[i]);
        }
        free(mggx);free(mggy);free(mggz);
        /* Markers */
        free(markx);free(marky);free(markz);
        free(marke);free(markt);
    }
    /**/
}
