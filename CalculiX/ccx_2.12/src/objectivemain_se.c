/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                    */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include "CalculiX.h"

static char *lakon1,*matname1;

static ITG *kon1,*ipkon1,*ne1,*nelcon1,*nrhcon1,*nalcon1,*ielmat1,*ielorien1,
    *norien1,*ntmat1_,*ithermal1,*iprestr1,*iperturb1,*iout1,*nmethod1,
    *nplicon1,*nplkcon1,*npmat1_,*mi1,*ielas1,*icmd1,*ncmat1_,*nstate1_,
    *istep1,*iinc1,calcul_qa1,nener1,ikin1,*istartdesi1,*ialdesi1,
    num_cpus,*ne01,*mortar1,*ielprop1,*ndesi1,*nodedesi1,idesvar1,
    *nobject1,iobject1;
    
static double *co1,*v1,*stx1,*elcon1,*rhcon1,*alcon1,*alzero1,*orab1,*t01,*t11,
    *prestr1,*vold1,*veold1,*dtime1,*time1,*xdesi1,
    *ttime1,*plicon1,*plkcon1,*xstateini1,*xstiff1,*xstate1,*stiini1,
    *vini1,*ener1,*eei1,*enerini1,*springarea1,*reltime1,*thicke1,*emeini1,
    *prop1,*pslavsurf1,*pmastsurf1,*clearini1,*distmin1,*g01,*dgdx1,
    *sti1,*xmass1=NULL,*xener1=NULL;
    

void objectivemain_se(double *co,ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
       double *v,double *stn,ITG *inum,double *stx,double *elcon,ITG *nelcon,
       double *rhcon,ITG *nrhcon,double *alcon,ITG *nalcon,double *alzero,
       ITG *ielmat,ITG *ielorien,ITG *norien,double *orab,ITG *ntmat_,
       double *t0,
       double *t1,ITG *ithermal,double *prestr,ITG *iprestr,char *filab,
       double *eme,double *emn,
       double *een,ITG *iperturb,double *f,double *fn,ITG *nactdof,ITG *iout,
       double *qa,double *vold,double *b,ITG *nodeboun,ITG *ndirboun,
       double *xboun,ITG *nboun,ITG *ipompc,ITG *nodempc,double *coefmpc,
       char *labmpc,ITG *nmpc,ITG *nmethod,double *cam,ITG *neq,double *veold,
       double *accold,double *bet,double *gam,double *dtime,double *time,
       double *ttime,double *plicon,ITG *nplicon,double *plkcon,
       ITG *nplkcon,double *xstateini,double *xstiff,double *xstate,ITG *npmat_,
       double *epn,char *matname,ITG *mi,ITG *ielas,ITG *icmd,ITG *ncmat_,
       ITG *nstate_,
       double *stiini,double *vini,ITG *ikboun,ITG *ilboun,double *ener,
       double *enern,double *emeini,double *xstaten,double *eei,double *enerini,
       double *cocon,ITG *ncocon,char *set,ITG *nset,ITG *istartset,
       ITG *iendset,
       ITG *ialset,ITG *nprint,char *prlab,char *prset,double *qfx,double *qfn,
       double *trab,
       ITG *inotr,ITG *ntrans,double *fmpc,ITG *nelemload,ITG *nload,
       ITG *ikmpc,ITG *ilmpc,
       ITG *istep,ITG *iinc,double *springarea,double *reltime, ITG *ne0,
       double *xforc, ITG *nforc, double *thicke,
       double *shcon,ITG *nshcon,char *sideload,double *xload,
       double *xloadold,ITG *icfd,ITG *inomat,double *pslavsurf,
       double *pmastsurf,ITG *mortar,ITG *islavact,double *cdn,
       ITG *islavnode,ITG *nslavnode,ITG *ntie,double *clearini,
       ITG *islavsurf,ITG *ielprop,double *prop,double *energyini,
       double *energy,double *distmin,ITG *ndesi,ITG *nodedesi,
       ITG *nobject,char *objectset,double *g0,double *dgdx,
       double *sti,double *df,ITG *nactdofinv,ITG *jqs,
       ITG *irows,ITG *idisplacement,ITG *nzs,char *jobnamec,ITG *isolver,
       ITG *icol,ITG *irow,ITG *jq,ITG *kode,double *cs,char *output,
       ITG *istartdesi,ITG *ialdesi,double *xdesi,char *orname,
       ITG *icoordinate,ITG *iev,double *d,double *z,double *au,double *ad,
       double *aub,double*adb,ITG *cyclicsymmetry,ITG *nzss,ITG *nev,
       ITG *ishapeenergy,double *fint,ITG *nlabel,ITG *igreen,ITG *nasym){

    char description[13]="            ",cflag[1]=" ",*filabl=NULL;

    ITG calcul_qa,nener,ikin,i,j,k,m,iobject,im,symmetryflag=0,inputformat=0,
	mt=mi[1]+1,mode=-1,noddiam=-1,ngraph=1,idesvar,nea,neb,nodeset,
        kscale=1,*iponoel=NULL,*inoel=NULL,idir,iorien,network=0,
        inorm=0,irand=0;

    double sigma=0.,ptime=0.,*temp=NULL,*bfix=NULL,*vnew=NULL,*dstn=NULL,
        freq,*c=NULL,orabsav[7],rotvec[3],a[9],pgauss[3];
    
    if(*nasym!=0){symmetryflag=2;inputformat=1;}

    /* variables for multithreading procedure */
    
    ITG sys_cpus,*ithread=NULL;
    char *env,*envloc,*envsys;
    
    num_cpus = 0;
    sys_cpus=0;
    
    /* explicit user declaration prevails */
    
    envsys=getenv("NUMBER_OF_CPUS");
    if(envsys){
        sys_cpus=atoi(envsys);
        if(sys_cpus<0) sys_cpus=0;
    }
    
    /* automatic detection of available number of processors */
    
    if(sys_cpus==0){
        sys_cpus = getSystemCPUs();
        if(sys_cpus<1) sys_cpus=1;
    }
    
    /* local declaration prevails, if strictly positive */
    
    envloc = getenv("CCX_NPROC_RESULTS");
    if(envloc){
        num_cpus=atoi(envloc);
        if(num_cpus<0){
            num_cpus=0;
        }else if(num_cpus>sys_cpus){
            num_cpus=sys_cpus;
        }
        
    }
    
    /* else global declaration, if any, applies */
    
    env = getenv("OMP_NUM_THREADS");
    if(num_cpus==0){
        if (env)
            num_cpus = atoi(env);
        if (num_cpus < 1) {
            num_cpus=1;
        }else if(num_cpus>sys_cpus){
            num_cpus=sys_cpus;
        }
    }
    
// next line is to be inserted in a similar way for all other paralell parts
    
    if(*ne<num_cpus) num_cpus=*ne;
    
    pthread_t tid[num_cpus];

    if((*idisplacement==1)||((*ishapeenergy==1)&&(iperturb[1]==1))){

	/* factor the system */

	if(*isolver==0){
#ifdef SPOOLES
	    spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
			   &symmetryflag,&inputformat,&nzs[2]);
#else
	    printf("*ERROR in objectivemain_se: the SPOOLES library is not linked\n\n");
	    FORTRAN(stop,());
#endif
	}
	else if(*isolver==4){
#ifdef SGI
	    token=1;
	    sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
	    printf("*ERROR in objectivemain_se: the SGI library is not linked\n\n");
	    FORTRAN(stop,());
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	    tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
	    printf("*ERROR in objectivemain_se: the TAUCS library is not linked\n\n");
	    FORTRAN(stop,());
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	    pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
			   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
	    printf("*ERROR in objectivemain_se: the PARDISO library is not linked\n\n");
	    FORTRAN(stop,());
#endif
	}
	
    }
    
    /* loop over all objective functions */

    for(m=0;m<*nobject;m++){
	if(strcmp1(&objectset[m*243],"MASS")==0){
	    iobject=m+1;
	    iobject1=iobject;
	    
	    /* OBJECTIVE: MASS */
	    
	    NNEW(xmass1,double,*ne);

            /* deactivating the elements which are not part of the 
               target function */

	    FORTRAN(actideacti,(set,nset,istartset,iendset,ialset,objectset,
                                ipkon,&iobject,ne));
 
            /* call without perturbation */
   
	    idesvar=0;

	    /* calculating the objective function and the derivatives */
	    
	    NNEW(g01,double,num_cpus**nobject);
	    
	    co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;v1=v;nelcon1=nelcon;rhcon1=rhcon;
	    ielmat1=ielmat;ielorien1=ielorien;norien1=norien;ntmat1_=ntmat_;vold1=vold;
            matname1=matname;mi1=mi;thicke1=thicke;mortar1=mortar;ielprop1=ielprop;
            prop1=prop;distmin1=distmin;ndesi1=ndesi;nodedesi1=nodedesi;
	    nobject1=nobject;iobject1=iobject;ne1=ne;istartdesi1=istartdesi;
            ialdesi1=ialdesi;xdesi1=xdesi;idesvar1=idesvar;
	    
	    if(((*nmethod!=4)&&(*nmethod!=5))||(iperturb[0]>1)){
		printf(" Using up to %" ITGFORMAT " cpu(s) for the mass sensitivity.\n\n", num_cpus);
	    }
	    
	    NNEW(ithread,ITG,num_cpus);
	    
	    /* Total difference of the mass */
	    /* create threads and wait */
	    
	    for(i=0; i<num_cpus; i++)  {
		ithread[i]=i;
		pthread_create(&tid[i], NULL, (void *)objectivemt_mass_dx, (void *)&ithread[i]);
	    }
	    
	    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
	    /* Assembling g0 */
	    
	    g0[m]=g01[m];
	    for(j=1;j<num_cpus;j++){
		g0[m]+=g01[m+j**nobject];
	    }
	    SFREE(g01);SFREE(ithread);
    
	    /* loop over the design variables (perturbation) */
    
	    for(idesvar=1;idesvar<=*ndesi;idesvar++){

		nea=istartdesi[idesvar-1];
		neb=istartdesi[idesvar]-1;

		FORTRAN(objective_mass_dx,(co,kon,ipkon,lakon,nelcon,rhcon,
			ielmat,ielorien,norien,ntmat1_,matname,mi,
			thicke,mortar,&nea,&neb,ielprop,prop,distmin,ndesi,nodedesi,
			nobject,g0,dgdx,&iobject,xmass1,
			istartdesi,ialdesi,xdesi,&idesvar));
	    }

	    SFREE(xmass1);

            /* reactivating all elements */

	    for(i=0;i<*ne;i++){
		if(ipkon[i]<-1) ipkon[i]=-2-ipkon[i];
	    }
	    
	}else if(strcmp1(&objectset[m*243],"SHAPEENERGY")==0){
	    iobject=m+1;
	    iobject1=iobject;
	    
	    /* OBJECTIVE: SHAPE ENERGY */
	    
	    NNEW(xener1,double,*ne);

            /* deactivating the elements which are not part of the 
               target function */

	    FORTRAN(actideacti,(set,nset,istartset,iendset,ialset,objectset,
                                ipkon,&iobject,ne));
 
            /* call without perturbation */
   
	    idesvar=0;
	    
	    /* calculating the objective function and the derivatives */
	    
	    NNEW(g01,double,num_cpus**nobject);
	    
	    co1=co;kon1=kon;ipkon1=ipkon;lakon1=lakon;ne1=ne;
	    stx1=stx;elcon1=elcon;nelcon1=nelcon;rhcon1=rhcon;
	    nrhcon1=nrhcon;alcon1=alcon;nalcon1=nalcon;alzero1=alzero;
	    ielmat1=ielmat;ielorien1=ielorien;norien1=norien;orab1=orab;
	    ntmat1_=ntmat_;t01=t0;t11=t1;ithermal1=ithermal;prestr1=prestr;
	    iprestr1=iprestr;iperturb1=iperturb;iout1=iout;
	    vold1=vold;nmethod1=nmethod;veold1=veold;dtime1=dtime;
	    time1=time;ttime1=ttime;plicon1=plicon;nplicon1=nplicon;
	    plkcon1=plkcon;nplkcon1=nplkcon;xstateini1=xstateini;
	    xstiff1=xstiff;xstate1=xstate;npmat1_=npmat_;matname1=matname;
	    mi1=mi;ielas1=ielas;icmd1=icmd;ncmat1_=ncmat_;nstate1_=nstate_;
	    stiini1=stiini;vini1=vini;ener1=ener;eei1=eei;enerini1=enerini;
	    istep1=istep;iinc1=iinc;springarea1=springarea;reltime1=reltime;
	    calcul_qa1=calcul_qa;nener1=nener;ikin1=ikin;ne01=ne0;thicke1=thicke;
	    emeini1=emeini;pslavsurf1=pslavsurf;pmastsurf1=pmastsurf;mortar1=mortar;
            clearini1=clearini;ielprop1=ielprop;prop1=prop;
	    distmin1=distmin;ndesi1=ndesi;nodedesi1=nodedesi;
	    nobject1=nobject;iobject1=iobject;sti1=sti;istartdesi1=istartdesi;
            ialdesi1=ialdesi;xdesi1=xdesi;idesvar1=idesvar;
	    
	    if(((*nmethod!=4)&&(*nmethod!=5))||(iperturb[0]>1)){
		printf(" Using up to %" ITGFORMAT " cpu(s) for the shape energy sensitivity.\n\n", num_cpus);
	    }
	    
	    NNEW(ithread,ITG,num_cpus);
	    
	    /* Total difference of the internal shape energy */
	    /* create threads and wait */
	    
	    for(i=0; i<num_cpus; i++)  {
		ithread[i]=i;
		pthread_create(&tid[i], NULL, (void *)objectivemt_shapeener_dx, (void *)&ithread[i]);
	    }
	    
	    for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
    
	    /* Assembling g0 */
	    
	    g0[m]=g01[m];
	    for(j=1;j<num_cpus;j++){
		g0[m]+=g01[m+j**nobject];
	    }
	    SFREE(g01);SFREE(ithread);
    
	    /* loop over the design variables (perturbation) */
    
	    for(idesvar=1;idesvar<=*ndesi;idesvar++){

		nea=istartdesi[idesvar-1];
		neb=istartdesi[idesvar]-1;
    
		FORTRAN(objective_shapeener_dx,(co,kon,ipkon,lakon,ne,
                       stx,elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,
	               ielmat,ielorien,norien,orab,ntmat1_,t0,t1,ithermal,prestr,
	               iprestr,iperturb,iout,vold,
	               nmethod,veold,dtime,time,ttime,plicon,nplicon,plkcon,
	               nplkcon,xstateini,xstiff,xstate,npmat1_,matname,mi,ielas,
	               icmd,ncmat1_,nstate1_,stiini,vini,ener,enerini,istep,iinc,
                       springarea,reltime,&calcul_qa,&nener,&ikin,          
                       ne0,thicke,emeini,pslavsurf,pmastsurf,mortar,clearini,
                       &nea,&neb,ielprop,prop,distmin,ndesi,nodedesi,
	               nobject,g0,dgdx,&iobject,sti,xener1,
		       istartdesi,ialdesi,xdesi,&idesvar));

	    }

            if(iperturb[1]==1){
	    	       
	       /* solve the system */
		    
		    if(*isolver==0){
#ifdef SPOOLES
			spooles_solve(fint,&neq[1]);
#endif
		    }
		    else if(*isolver==4){
#ifdef SGI
			sgi_solve(fint,token);
#endif
		    }
		    else if(*isolver==5){
#ifdef TAUCS
			tau_solve(fint,&neq[1]);
#endif
		    }
		    else if(*isolver==7){
#ifdef PARDISO
			pardiso_solve(fint,&neq[1],&symmetryflag);
#endif
}	      
	      
	       /* solve the system */	      
	       	    
	    }
	    
	    SFREE(xener1);

            /* reactivating all elements */

	    for(i=0;i<*ne;i++){
		if(ipkon[i]<-1) ipkon[i]=-2-ipkon[i];
	    }

            /* composing the total derivative */

	    FORTRAN(objective_shapeener_tot,(dgdx,df,vold,ndesi,&iobject,
					     mi,nactdofinv,jqs,irows,iperturb,fint));
	    
	}else if((strcmp1(&objectset[m*243],"EIGENFREQUENCY")==0)||
                 (strcmp1(&objectset[m*243],"GREEN")==0)){
	    iobject=m+1;
	    
	    /* OBJECTIVE: EIGENFREQUENCY */

	    if(*igreen!=1){

                /* determination of the sensitivity of the eigenvalues */

		if(!*cyclicsymmetry){
		    
		    FORTRAN(objective_freq,(dgdx,df,v,ndesi,&iobject,
					    mi,nactdofinv,jqs,irows));
		    
		    /* change sign since df contains -(dK/dX-lambda*dM/DX).U */
		    
		    for(idesvar=0;idesvar<*ndesi;idesvar++){dgdx[idesvar]=-dgdx[idesvar];}
		}else{
		    
		    FORTRAN(objective_freq_cs,(dgdx,df,v,ndesi,&iobject,
					       mi,nactdofinv,jqs,irows,nk,nzss));
		}
	    }

	    g0[m]=d[*iev];

            /* in case the design variables are the orientations
               the sensitivity of the eigenvectors is also
               determined */

	    if(*icoordinate!=1){
		if(*igreen!=1) FORTRAN(writedeigdx,(iev,d,ndesi,orname,dgdx));

  	        /* createinum is called in order to determine the nodes belonging
	           to elements; this information is needed in frd_se */
	    
		NNEW(inum,ITG,*nk);
		FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],
                    nelemload,nload,nodeboun,nboun,ndirboun,ithermal,co,
                    vold,mi,ielmat));
		
                /* the frequency is also needed for frd_se */

		if(d[*iev]>=0.){
		    freq=sqrt(d[*iev])/6.283185308;
		}else{
		    freq=0.;
		}

                /* determine the derivative of the eigenvectors */

		NNEW(bfix,double,*neq);
		NNEW(b,double,*neq);
		NNEW(temp,double,mt**nk);

		if(*igreen!=1){
		    
		    /* bfix = M * eigenvector */
		    
		    FORTRAN(op,(neq,&z[*iev**neq],bfix,adb,aub,jq,irow));

		}else{		

		    sigma=d[*iev];
		    
		    /* factor the system */
		    
		    if(*isolver==0){
#ifdef SPOOLES
			spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
				   &symmetryflag,&inputformat,&nzs[2]);
#else
			printf("*ERROR in objectivemain_se: the SPOOLES library is not linked\n\n");
			FORTRAN(stop,());
#endif
		    }
		    else if(*isolver==4){
#ifdef SGI
			token=1;
			sgi_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],token);
#else
			printf("*ERROR in objectivemain_se: the SGI library is not linked\n\n");
			FORTRAN(stop,());
#endif
		    }
		    else if(*isolver==5){
#ifdef TAUCS
			tau_factor(ad,&au,adb,aub,&sigma,icol,&irow,&neq[1],&nzs[1]);
#else
			printf("*ERROR in objectivemain_se: the TAUCS library is not linked\n\n");
			FORTRAN(stop,());
#endif
		    }
		    else if(*isolver==7){
#ifdef PARDISO
			pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq[1],&nzs[1],
				       &symmetryflag,&inputformat,jq,&nzs[2]);
#else
			printf("*ERROR in objectivemain_se: the PARDISO library is not linked\n\n");
			FORTRAN(stop,());
#endif
		    }
		}	
		
                /* loop over all design variables */
	
		for(idesvar=0;idesvar<*ndesi;idesvar++){

                    /* setting up the RHS of the system */

		    if(*igreen!=1){
			for(j=0;j<*neq;j++){
			    b[j]=dgdx[idesvar]*bfix[j];
			}
		    }else{
			DMEMSET(b,0,*neq,0.);
		    }

		    for(j=jqs[idesvar]-1;j<jqs[idesvar+1]-1;j++){
			b[irows[j]-1]+=df[j];
		    }

		    if(*igreen==1){

			/* solve the system */
			
			if(*isolver==0){
#ifdef SPOOLES
			    spooles_solve(b,&neq[1]);
#endif
			}
			else if(*isolver==4){
#ifdef SGI
			    sgi_solve(b,token);
#endif
			}
			else if(*isolver==5){
#ifdef TAUCS
			    tau_solve(b,&neq[1]);
#endif
			}
			else if(*isolver==7){
#ifdef PARDISO
			    pardiso_solve(b,&neq[1],&symmetryflag);
#endif
			}
		    }else{
		    
			NNEW(c,double,*nev);
			for(j=0;j<*nev;j++){
			    if(j==*iev) continue;
			    for(k=0;k<*neq;k++){
				c[j]+=z[j**neq+k]*b[k];
			    }
			    c[j]/=(d[j]-d[*iev]);
			}
			DMEMSET(b,0,*neq,0.);
			for(j=0;j<*nev;j++){
			    if(j==*iev) continue;
			    for(k=0;k<*neq;k++){
				b[k]+=c[j]*z[j**neq+k];
			    }
			}
			SFREE(c);
		    }

		    /* store the answer in temp w.r.t. node and direction
		       instead of w.r.t. dof */
		    
		    DMEMSET(temp,0,mt**nk,0.);
		    FORTRAN(resultsnoddir,(nk,temp,nactdof,b,ipompc,nodempc,
				      coefmpc,nmpc,mi));
		    
		    /* storing the sensitivity of the eigenmodes to file */
		    
		    ++*kode;
		    frd_sen(co,nk,stn,inum,nmethod,kode,filab,
                       &freq,nstate_,
		       istep,iinc,&mode,&noddiam,description,mi,&ngraph,
                       ne,cs,set,nset,istartset,iendset,ialset,
		       jobnamec,output,temp,&iobject,objectset,ntrans,
		       inotr,trab,&idesvar,orname,icoordinate,&inorm,
                       &irand); 

		}  // enddo loop idesvar

		if(*igreen==1){
		    
		    /* clean the system */
		    
		    if(*isolver==0){
#ifdef SPOOLES
			spooles_cleanup();
#endif
		    }
		    else if(*isolver==4){
#ifdef SGI
			sgi_cleanup(token);
#endif
		    }
		    else if(*isolver==5){
#ifdef TAUCS
			tau_cleanup();
#endif
		    }
		    else if(*isolver==7){
#ifdef PARDISO
			pardiso_cleanup(&neq[1],&symmetryflag);
#endif
		    }
		}

		SFREE(temp);SFREE(bfix);SFREE(b);SFREE(inum);

	    }

	}else if(strcmp1(&objectset[m*243],"DISPLACEMENT")==0){
	    iobject=m+1;

            /* OBJECTIVE: DISPLACEMENT */
	    
	    /* createinum is called in order to determine the nodes belonging
	       to elements; this information is needed in frd_se */
	    
	    NNEW(inum,ITG,*nk);
	    FORTRAN(createinum,(ipkon,inum,kon,lakon,nk,ne,&cflag[0],nelemload,
		    nload,nodeboun,nboun,ndirboun,ithermal,co,vold,mi,ielmat));
	    
	    NNEW(b,double,*neq);
	    NNEW(temp,double,mt**nk);

            /* if the design variables are the coordinates:
               check for the existence of a target node set */

            /* calculating the objective function */

	    if(*icoordinate==1){
		nodeset=0;
		for(i=0;i<*nset;i++){
		    if(strcmp1(&objectset[m*243+162]," ")==0) continue;
		    if(strcmp2(&objectset[m*243+162],&set[i*81],81)==0){
			nodeset=i+1;
			break;
		    }
		}
		FORTRAN(objective_disp,(&nodeset,istartset,iendset,
			ialset,nk,&idesvar,&iobject,mi,g0,
			nobject,vold));
	    }
	    
	    for(idesvar=0;idesvar<*ndesi;idesvar++){
		
                /* copying the RHS from field df */

		DMEMSET(b,0,*neq,0.);
		for(j=jqs[idesvar]-1;j<jqs[idesvar+1]-1;j++){
		    b[irows[j]-1]=df[j];
		}

                /* solve the system */

		if(*isolver==0){
#ifdef SPOOLES
		    spooles_solve(b,&neq[1]);
#endif
		}
		else if(*isolver==4){
#ifdef SGI
		    sgi_solve(b,token);
#endif
		}
		else if(*isolver==5){
#ifdef TAUCS
		    tau_solve(b,&neq[1]);
#endif
		}
		else if(*isolver==7){
#ifdef PARDISO
		    pardiso_solve(b,&neq[1],&symmetryflag);
#endif
		}

		if(*icoordinate!=1){
		    
		    /* store the answer in temp w.r.t. node and direction
		       instead of w.r.t. dof */
		    
		    DMEMSET(temp,0,mt**nk,0.);
		    FORTRAN(resultsnoddir,(nk,temp,nactdof,b,ipompc,nodempc,
				      coefmpc,nmpc,mi));
		    
		    /* storing the results to file */
		    
		    ++*kode;
		    frd_sen(co,nk,stn,inum,nmethod,kode,filab,
                       &ptime,nstate_,
		       istep,iinc,&mode,&noddiam,description,mi,&ngraph,
                       ne,cs,set,nset,istartset,iendset,ialset,
		       jobnamec,output,temp,&iobject,objectset,ntrans,
		       inotr,trab,&idesvar,orname,icoordinate,&inorm,
                       &irand); 

		}else{
		    FORTRAN(objective_disp_dx,(&nodeset,istartset,iendset,
			    ialset,nk,&idesvar,&iobject,mi,nactdof,dgdx,
                            ndesi,nobject,vold,b));
		}
	    }

	    SFREE(b);SFREE(temp);SFREE(inum);
	    
	}else if(strcmp1(&objectset[m*243],"STRESS")==0){
	    iobject=m+1;

	    NNEW(filabl,char,87**nlabel);
	    for(i=0;i<87**nlabel;i++){strcpy1(&filabl[i]," ",1);}
	    strcpy1(&filabl[174],"S   ",4);
	    
	    /* OBJECTIVE: STRESS */

            /* calculating the stress in the unperturbed state */
  
	    NNEW(v,double,mt**nk);
	    NNEW(fn,double,mt**nk);
	    NNEW(stn,double,6**nk);
	    NNEW(inum,ITG,*nk);
	    NNEW(stx,double,6*mi[0]**ne);
	    NNEW(eei,double,6*mi[0]**ne);
	    
	    memcpy(&v[0],&vold[0],sizeof(double)*mt**nk);
	    *iout=2;
	    *icmd=3;
	    
	    results(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
		    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		    ielorien,norien,orab,ntmat_,t0,t1,ithermal,
		    prestr,iprestr,filabl,eme,emn,een,iperturb,
		    f,fn,nactdof,iout,qa,vold,b,nodeboun,
		    ndirboun,xboun,nboun,ipompc,
		    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		    bet,gam,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
		    xstateini,xstiff,xstate,npmat_,epn,matname,mi,ielas,icmd,
		    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		    xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		    ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		    nelemload,nload,ikmpc,ilmpc,istep,iinc,springarea,
		    reltime,ne0,xforc,nforc,thicke,shcon,nshcon,
		    sideload,xload,xloadold,icfd,inomat,pslavsurf,pmastsurf,
		    mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		    islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
                    inoel,&nener,orname,&network);
	    
	    *icmd=0;
	    
	    SFREE(v);SFREE(fn);SFREE(stx);SFREE(eei);
	    
	    NNEW(b,double,*neq);
	    NNEW(temp,double,mt**nk);

            /* if the design variables are the coordinates:
               check for the existence of a target node set */

            /* calculating the objective function */

	    if(*icoordinate==1){
		nodeset=0;
		for(i=0;i<*nset;i++){
		    if(strcmp1(&objectset[m*243+162]," ")==0) continue;
		    if(strcmp2(&objectset[m*243+162],&set[i*81],81)==0){
			nodeset=i+1;
			break;
		    }
		}
		FORTRAN(objective_stress,(&nodeset,istartset,iendset,
			ialset,nk,&idesvar,&iobject,mi,g0,
			nobject,stn,objectset));
	    }
	    
	    for(idesvar=0;idesvar<*ndesi;idesvar++){
		
                /* copying the RHS from field df */

		DMEMSET(b,0,*neq,0.);
		for(j=jqs[idesvar]-1;j<jqs[idesvar+1]-1;j++){
		    b[irows[j]-1]=df[j];
		}

                /* solve the system */

		if(*isolver==0){
#ifdef SPOOLES
		    spooles_solve(b,&neq[1]);
#endif
		}
		else if(*isolver==4){
#ifdef SGI
		    sgi_solve(b,token);
#endif
		}
		else if(*isolver==5){
#ifdef TAUCS
		    tau_solve(b,&neq[1]);
#endif
		}
		else if(*isolver==7){
#ifdef PARDISO
		    pardiso_solve(b,&neq[1],&symmetryflag);
#endif
		}

                /* calculating the perturbed displacements */

		NNEW(vnew,double,mt**nk);
		    
		FORTRAN(resultsnoddir,(nk,vnew,nactdof,b,ipompc,nodempc,
				       coefmpc,nmpc,mi));

		for(i=0;i<mt**nk;i++){vnew[i]=vold[i]+(*distmin)*vnew[i];}

                /* calculating the stress in the perturbed state */
  
		NNEW(v,double,mt**nk);
		NNEW(fn,double,mt**nk);
		NNEW(stx,double,6*mi[0]**ne);
		NNEW(eei,double,6*mi[0]**ne);
		NNEW(dstn,double,6**nk);
		
		memcpy(&v[0],&vnew[0],sizeof(double)*mt**nk);
		*iout=2;
		*icmd=3;
	   
		/* calculate a delta in the orientation
		   in case the material orientation is the design variable */
	
		if(*icoordinate!=1){
		    iorien=idesvar/3;
		    
		    /* save nominal orientation */
		    
		    memcpy(&orabsav[0],&orab[7*iorien],sizeof(double)*7);
		    
		    /* calculate the transformation matrix */
		    
		    FORTRAN(transformatrix,(&orab[7*iorien],pgauss,a));
		    
		    /* calculate the rotation vector from the transformation matrix */
		    
		    FORTRAN(rotationvector,(a,rotvec));
		    idir=idesvar-iorien*3;
		    
		    /* add a small variation to the rotation vector component */
		    
		    rotvec[idir]+=*distmin;
		    
		    /* determine the new transformation matrix */
		    
		    FORTRAN(rotationvectorinv,(a,rotvec));
		    
		    /* determine two new points in the x-y plane */
		    
		    for(i=0;i<6;i++){orab[7*iorien+i]=a[i];}
		}
		
		results(co,nk,kon,ipkon,lakon,ne,v,dstn,inum,stx,
		    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
		    ielorien,norien,orab,ntmat_,t0,t1,ithermal,
		    prestr,iprestr,filabl,eme,emn,een,iperturb,
		    f,fn,nactdof,iout,qa,vold,b,nodeboun,
		    ndirboun,xboun,nboun,ipompc,
		    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,&neq[1],veold,accold,
		    bet,gam,dtime,time,ttime,plicon,nplicon,plkcon,nplkcon,
		    xstateini,xstiff,xstate,npmat_,epn,matname,mi,ielas,icmd,
		    ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,emeini,
		    xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,iendset,
		    ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,fmpc,
		    nelemload,nload,ikmpc,ilmpc,istep,iinc,springarea,
		    reltime,ne0,xforc,nforc,thicke,shcon,nshcon,
		    sideload,xload,xloadold,icfd,inomat,pslavsurf,pmastsurf,
		    mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
		    islavsurf,ielprop,prop,energyini,energy,&kscale,iponoel,
		    inoel,&nener,orname,&network);
	    
		*icmd=0;
		
		SFREE(v);SFREE(fn);SFREE(stx);SFREE(eei);

                /* calculate the stress sensitivity */

		for(i=0;i<6**nk;i++){dstn[i]=(dstn[i]-stn[i])/(*distmin);}
		SFREE(vnew);

		if(*icoordinate!=1){
	   
		    /* restoring the nominal orientation  */
	
		    if(idesvar>0){
			memcpy(&orab[7*iorien],&orabsav[0],sizeof(double)*7);
		    }
		    
		    /* storing the results to file */
		    
		    ++*kode;
		    frd_sen(co,nk,dstn,inum,nmethod,kode,filab,
                       &ptime,nstate_,
		       istep,iinc,&mode,&noddiam,description,mi,&ngraph,
                       ne,cs,set,nset,istartset,iendset,ialset,
		       jobnamec,output,temp,&iobject,objectset,ntrans,
		       inotr,trab,&idesvar,orname,icoordinate,&inorm,
                       &irand); 

		}else{
		    FORTRAN(objective_stress_dx,(&nodeset,istartset,iendset,
			    ialset,nk,&idesvar,&iobject,dgdx,
			    ndesi,nobject,stn,dstn,objectset,g0));
		}

		SFREE(dstn);

	    }

	    SFREE(b);SFREE(temp);SFREE(inum);SFREE(stn);SFREE(filabl);
	    
	}

    }
    
    if(*idisplacement==1){

	/* clean the system */

	if(*isolver==0){
#ifdef SPOOLES
	    spooles_cleanup();
#endif
	}
	else if(*isolver==4){
#ifdef SGI
	    sgi_cleanup(token);
#endif
	}
	else if(*isolver==5){
#ifdef TAUCS
	    tau_cleanup();
#endif
	}
	else if(*isolver==7){
#ifdef PARDISO
	    pardiso_cleanup(&neq[1],&symmetryflag);
#endif
	}
    }
    
    return;
    
}

/* ----------------------------------------------------------------*/
/* subroutine for multithreading: Differentiation of shape energy  */
/* ----------------------------------------------------------------*/

void *objectivemt_shapeener_dx(ITG *i){
    
    ITG nea,neb,nedelta,indexg0,indexdgdx;
    
    indexg0=*i**nobject1;
    indexdgdx=*i**nobject1**ndesi1;

    nedelta=(ITG)floor(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if((*i==num_cpus-1)&&(neb<*ne1)) neb=*ne1;
    
    FORTRAN(objective_shapeener_dx,(co1,kon1,ipkon1,lakon1,ne1,
	  stx1,elcon1,nelcon1,rhcon1,nrhcon1,alcon1,nalcon1,alzero1,
	  ielmat1,ielorien1,norien1,orab1,ntmat1_,t01,t11,ithermal1,prestr1,
	  iprestr1,iperturb1,iout1,vold1,
	  nmethod1,veold1,dtime1,time1,ttime1,plicon1,nplicon1,plkcon1,
	  nplkcon1,xstateini1,xstiff1,xstate1,npmat1_,matname1,mi1,ielas1,
	  icmd1,ncmat1_,nstate1_,stiini1,vini1,ener1,enerini1,istep1,iinc1,
          springarea1,reltime1,&calcul_qa1,&nener1,&ikin1,          
          ne01,thicke1,emeini1,pslavsurf1,pmastsurf1,mortar1,clearini1,
          &nea,&neb,ielprop1,prop1,distmin1,ndesi1,nodedesi1,
	  nobject1,&g01[indexg0],&dgdx1[indexdgdx],&iobject1,sti1,xener1,
	  istartdesi1,ialdesi1,xdesi1,&idesvar1));

    return NULL;
}

/* ---------------------------------------------------*/
/* subroutine for multithreading of objective_mass    */
/* ---------------------------------------------------*/

void *objectivemt_mass_dx(ITG *i){

    ITG nea,neb,nedelta,indexg0,indexdgdx;
    
    indexg0=*i**nobject1;
    indexdgdx=*i**nobject1**ndesi1;

    nedelta=(ITG)floor(*ne1/(double)num_cpus);
    nea=*i*nedelta+1;
    neb=(*i+1)*nedelta;
    if((*i==num_cpus-1)&&(neb<*ne1)) neb=*ne1;

    FORTRAN(objective_mass_dx,(co1,kon1,ipkon1,lakon1,nelcon1,rhcon1,          
          ielmat1,ielorien1,norien1,ntmat1_,matname1,mi1,
          thicke1,mortar1,&nea,&neb,ielprop1,prop1,distmin1,ndesi1,nodedesi1,
	  nobject1,&g01[indexg0],&dgdx1[indexdgdx],&iobject1,xmass1,
	  istartdesi1,ialdesi1,xdesi1,&idesvar1));
          
    return NULL;
}
