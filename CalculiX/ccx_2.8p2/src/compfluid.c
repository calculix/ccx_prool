/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                     */

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
#include <pthread.h>
#include "CalculiX.h"
#ifdef SPOOLES
#include "spooles.h"
#endif
#ifdef SGI
#include "sgi.h"
#endif
#ifdef TAUCS
#include "tau.h"
#endif
#ifdef PARDISO
#include "pardiso.h"
#endif

static ITG num_cpus;

void compfluid(double **cop, ITG *nk, ITG **ipkonfp, ITG **konp, char **lakonfp,
    char **sidefacep, ITG *ifreestream, 
    ITG *nfreestream, ITG *isolidsurf, ITG *neighsolidsurf,
    ITG *nsolidsurf, ITG **iponoelp, ITG **inoelp, ITG *nshcon, double *shcon,
    ITG *nrhcon, double *rhcon, double **voldp, ITG *ntmat_,ITG *nodeboun, 
    ITG *ndirboun, ITG *nboun, ITG **ipompcp,ITG **nodempcp, ITG *nmpc,
    ITG **ikmpcp, ITG **ilmpcp, ITG *ithermal, ITG *ikboun, ITG *ilboun,
    ITG *iturbulent, ITG *isolver, ITG *iexpl, double *vcontu, double *ttime,
    double *time, double *dtime, ITG *nodeforc,ITG *ndirforc,double *xforc,
    ITG *nforc, ITG *nelemload, char *sideload, double *xload,ITG *nload,
    double *xbody,ITG *ipobody,ITG *nbody, ITG *ielmatf, char *matname,
    ITG *mi, ITG *ncmat_, double *physcon, ITG *istep, ITG *iinc,
    ITG *ibody, double *xloadold, double *xboun,
    double **coefmpcp, ITG *nmethod, double *xforcold, double *xforcact,
    ITG *iamforc,ITG *iamload, double *xbodyold, double *xbodyact,
    double *t1old, double *t1, double *t1act, ITG *iamt1, double *amta,
    ITG *namta, ITG *nam, double *ampli, double *xbounold, double *xbounact,
    ITG *iamboun, ITG *itg, ITG *ntg, char *amname, double *t0, 
    ITG **nelemfacep,
    ITG *nface, double *cocon, ITG *ncocon, double *xloadact, double *tper,
    ITG *jmax, ITG *jout, char *set, ITG *nset, ITG *istartset,
    ITG *iendset, ITG *ialset, char *prset, char *prlab, ITG *nprint,
    double *trab, ITG *inotr, ITG *ntrans, char *filab, char **labmpcp, 
    double *sti, ITG *norien, double *orab, char *jobnamef,char *tieset,
    ITG *ntie, ITG *mcs, ITG *ics, double *cs, ITG *nkon, ITG *mpcfree,
    ITG *memmpc_,double **fmpcp,ITG *nef,ITG **inomatp,double *qfx,
    ITG *neifa,ITG *neiel,ITG *ielfa,ITG *ifaext,double *vfa,double *vel,
    ITG *ipnei,ITG *nflnei,ITG *nfaext,char *typeboun,ITG *neij,
    double *tincf,ITG *nactdoh,ITG *nactdohinv,ITG *ielorienf){

    /* main computational fluid dynamics routine */
  
  char cflag[1],*labmpc=NULL,*lakonf=NULL,*sideface=NULL;

  char matvec[7]="MATVEC",msolve[7]="MSOLVE";

  ITG *ipointer=NULL,*mast1=NULL,*irow=NULL,*icol=NULL,*jq=NULL,
      nzs=20000000,neq,kode,compressible,*ifabou=NULL,*ja=NULL,
      *nodempc=NULL,*ipompc=NULL,*ikmpc=NULL,*ilmpc=NULL,nfabou,im,
      *ipkonf=NULL,*kon=NULL,*nelemface=NULL,*inoel=NULL,last=0,
      *iponoel=NULL,*inomat=NULL,ithermalref,*integerglob=NULL,iit,
      iconvergence=0,symmetryflag,inputformat,i,*inum=NULL,iitf,ifreefa,
      *iponofa=NULL,*inofa=NULL,is,ie,*ia=NULL,nstate_,*ielpropf=NULL,
      icent=0,isti=0,iqfx=0,nfield,ndim,iorienglob,force=0,icfdout=1;

  ITG nelt,isym,itol,itmax,iunit,lrgw,*igwk=NULL,ligw,ierr,*iwork=NULL,iter,
      nsave,lenw,leniw;

  double *coefmpc=NULL,*fmpc=NULL,*umfa=NULL,reltime,*doubleglob=NULL,
      *co=NULL,*vold=NULL,*coel=NULL,*cosa=NULL,*gradvel=NULL,*gradvfa=NULL,
      *xxn=NULL,*xxi=NULL,*xle=NULL,*xlen=NULL,*xlet=NULL,timef,dtimef,
      *cofa=NULL,*area=NULL,*xrlfa=NULL,reltimef,ttimef,*hcfa=NULL,*cpel=NULL,
      *au=NULL,*ad=NULL,*b=NULL,*volume=NULL,*body=NULL,sigma=0.,betam,
      *adb=NULL,*aub=NULL,*advfa=NULL,*ap=NULL,*bp=NULL,*xxj=NULL,
      *v=NULL,*velo=NULL,*veloo=NULL,*gammat=NULL,*cosb=NULL,dmin,tincfguess,
      *hel=NULL,*hfa=NULL,*auv=NULL,*adv=NULL,*bv=NULL,*sel=NULL,*gamma=NULL,
      *gradtfa=NULL,*gradtel=NULL,*umel=NULL,*cpfa=NULL,*gradpel=NULL,
      *fn=NULL,*eei=NULL,*xstate=NULL,*ener=NULL,*thicke=NULL,*eme=NULL,
      ptimef,*stn=NULL,*qfn=NULL,*hcel=NULL,*aua=NULL,a1,a2,a3,beta,
      *prop=NULL;

  double tol,*rgwk=NULL,err,*sb=NULL,*sx=NULL,*rwork=NULL;

  nodempc=*nodempcp;ipompc=*ipompcp;ikmpc=*ikmpcp;ilmpc=*ilmpcp;
  coefmpc=*coefmpcp;labmpc=*labmpcp;fmpc=*fmpcp;co=*cop;
  ipkonf=*ipkonfp;lakonf=*lakonfp;kon=*konp;
  nelemface=*nelemfacep;sideface=*sidefacep;inoel=*inoelp;
  iponoel=*iponoelp;vold=*voldp;inomat=*inomatp;

#ifdef SGI
  ITG token;
#endif

  /* relative time at the end of the mechanical increment */

  reltime=(*time)/(*tper);

  /* open frd-file for fluids */

  FORTRAN(openfilefluid,(jobnamef));

  /* variables for multithreading procedure */

  ITG sys_cpus;
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
  
  envloc = getenv("CCX_NPROC_CFD");
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
  
  if(*nef<num_cpus) num_cpus=*nef;
  
  printf(" Using up to %" ITGFORMAT " cpu(s) for CFD.\n", num_cpus);
  
  pthread_t tid[num_cpus];

  
  kode=0;
  
  /*  *iexpl==0:  structure:implicit, fluid:incompressible
      *iexpl==1:  structure:implicit, fluid:compressible
      *iexpl==2:  structure:explicit, fluid:incompressible
      *iexpl==3:  structure:explicit, fluid:compressible */

  if((*iexpl==1)||(*iexpl==3)){
      compressible=1;
  }else{
      compressible=0;
  }

  /* if initial conditions are specified for the temperature, 
     it is assumed that the temperature is an unknown */

  ithermalref=*ithermal;
  if(*ithermal==1){
    *ithermal=2;
  }

  /* determining the matrix structure */
  
  NNEW(ipointer,ITG,3**nk);
  NNEW(mast1,ITG,nzs);
  NNEW(irow,ITG,nzs);
  NNEW(ia,ITG,nzs);
  NNEW(icol,ITG,*nef);
  NNEW(jq,ITG,*nef+1);
  NNEW(ja,ITG,*nef+1);
//  NNEW(nactdoh,ITG,*nef);

  mastructf(nk,kon,ipkonf,lakonf,nef,icol,jq,&mast1,&irow,
	    isolver,&neq,ipointer,&nzs,ipnei,neiel,mi);

//  printf("Unterschied start\n");
//  for(i=0;i<*ne;i++){if(i+1!=nactdoh[i]){printf("Unterschied i=%d,nactdoh[i]=%d\n",i+1,nactdoh[i]);}}
//  printf("Unterschied end\n");

  SFREE(ipointer);SFREE(mast1);
 
  /* calculation geometric data */

  NNEW(coel,double,3**nef);
  NNEW(volume,double,*nef);
  NNEW(cosa,double,*nflnei);
  NNEW(cosb,double,*nflnei);
  NNEW(xxn,double,3**nflnei);
  NNEW(xxi,double,3**nflnei);
  NNEW(xxj,double,3**nflnei);
  NNEW(xle,double,*nflnei);
  NNEW(xlen,double,*nflnei);
  NNEW(xlet,double,*nflnei);
  NNEW(cofa,double,3**nface);
  NNEW(area,double,*nface);
  NNEW(xrlfa,double,3**nface);

  FORTRAN(initialcfd,(nef,ipkonf,kon,lakonf,co,coel,cofa,nface,
	  ielfa,area,ipnei,neiel,xxn,xxi,xle,xlen,xlet,xrlfa,cosa,
	  volume,neifa,xxj,cosb,vel,&dmin));

  /* storing pointers to the boundary conditions in ielfa */

  NNEW(ifabou,ITG,7**nfaext);
  FORTRAN(applyboun,(ifaext,nfaext,ielfa,ikboun,ilboun,
       nboun,typeboun,nelemload,nload,sideload,isolidsurf,nsolidsurf,
       ifabou,&nfabou,nface,nodeboun,ndirboun,ikmpc,ilmpc,labmpc,nmpc,
       nactdohinv));
  RENEW(ifabou,ITG,nfabou);

  /* catalogueing the nodes for output purposes (interpolation at
     the nodes */
  
  NNEW(iponofa,ITG,*nk);
  NNEW(inofa,ITG,2**nface*4);

  FORTRAN(cataloguenodes,(iponofa,inofa,&ifreefa,ielfa,ifabou,ipkonf,
			  kon,lakonf,nface,nk));

  RENEW(inofa,ITG,2*ifreefa);

  /* material properties for athermal calculations 
     = calculation for which no initial thermal conditions
     were defined */

  NNEW(umfa,double,*nface);
  NNEW(umel,double,*nef);
      
  /* calculating the density at the element centers */
  
  FORTRAN(calcrhoel,(nef,vel,rhcon,nrhcon,ielmatf,ntmat_,
		     ithermal,mi));
  
  /* calculating the density at the face centers */
  
  FORTRAN(calcrhofa,(nface,vfa,rhcon,nrhcon,ielmatf,ntmat_,
		     ithermal,mi,ielfa));
  
  /* calculating the dynamic viscosity at the face centers */
  
  FORTRAN(calcumfa,(nface,vfa,shcon,nshcon,ielmatf,ntmat_,
		    ithermal,mi,ielfa,umfa));
  
  /* calculating the dynamic viscosity at the element centers */
  
  FORTRAN(calcumel,(nef,vel,shcon,nshcon,ielmatf,ntmat_,
		    ithermal,mi,umel));
  
  if(*ithermal!=0){
      NNEW(hcfa,double,*nface);
      NNEW(cpel,double,*nef);
      NNEW(cpfa,double,*nface);
  }
      

  if(*nbody>0) NNEW(body,double,4**nef);

  /* v is a auxiliary field: set to zero for the calls to
     tempload */

  NNEW(v,double,5**nk);

  /* next section is for stationary calculations */
  
  if(*nmethod==1){
      
      /* boundary conditions at the end of the mechanical
	 increment */
      
      FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
	     xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,time,&reltime,ttime,dtime,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,v,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
	     ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
             iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));

      /* body forces (gravity, centrifugal and Coriolis forces */

      if(*nbody>0){
	  FORTRAN(inicalcbody,(nef,body,ipobody,ibody,xbody,coel,vel,lakonf,
                            nactdohinv,&icent));
      }
  }

  /* extrapolating the velocity from the elements centers to the face
     centers, thereby taking the boundary conditions into account */

  FORTRAN(extrapolate_vel,(nface,ielfa,xrlfa,vel,vfa,
			   ifabou,xbounact,ipnei,nef));

  /* applying MPC's to the faces */

  if(*nmpc>0){
      is=1;ie=3;
      FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfa,coefmpc,
			nodempc,ipnei,neifa,labmpc,xbounact,nactdoh));
  }

  /* extrapolation of the pressure at the element centers
     to the face centers */
  
  FORTRAN(extrapolate_pel,(nface,ielfa,xrlfa,vel,vfa,ifabou,
			   xbounact,nef));

  /* applying MPC's to the faces */

  if(*nmpc>0){
      is=4;ie=4;
      FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfa,coefmpc,
			nodempc,ipnei,neifa,labmpc,xbounact,nactdoh));
  }

  /* extrapolation of the temperature at the element centers
     to the face centers */

  if(*ithermal>0){
      FORTRAN(extrapolate_tel,(nface,ielfa,xrlfa,vel,vfa,
			       ifabou,xbounact,ipnei,nef));

      /* applying MPC's to the faces */
      
      if(*nmpc>0){
	  is=0;ie=0;
	  FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfa,coefmpc,
			    nodempc,ipnei,neifa,labmpc,xbounact,nactdoh));
      }
  }

  /* calculating the maximum velocity (for the determination of the
     time increment) */

  FORTRAN(calcguesstincf,(nface,&dmin,vfa,umfa,&tincfguess));

  /* start of the major loop */

  NNEW(gradvel,double,9**nef);
  NNEW(gradvfa,double,9**nface);

  NNEW(advfa,double,*nface);
  NNEW(hfa,double,3**nface);

  NNEW(ap,double,*nface);
  NNEW(bp,double,*nface);

  NNEW(au,double,2*nzs+neq);
  NNEW(aua,double,nzs+neq);
//  NNEW(au,double,2*nzs);
  NNEW(ad,double,neq);
  NNEW(b,double,neq);

  NNEW(auv,double,2*nzs+neq);
//  NNEW(auv,double,2*nzs);
  NNEW(adv,double,neq);
  NNEW(bv,double,3*neq);
  NNEW(hel,double,3*neq);
  NNEW(sel,double,3*neq);

  NNEW(rwork,double,neq);

  NNEW(gradpel,double,3**nef);

  if(*ithermal>0){
      NNEW(gradtel,double,3**nef);
      NNEW(gradtfa,double,3**nface);
  }

  NNEW(inum,ITG,*nk);

  NNEW(velo,double,6**nef);
  NNEW(veloo,double,6**nef);

  /* initializing velo and veloo */

  memcpy(&veloo[0],&vel[0],sizeof(double)*6**nef);
  memcpy(&velo[0],&vel[0],sizeof(double)*6**nef);

  iit=0;

  a1=1.5;
  a2=-2.;
  a3=0.5;

  if(*tincf<=0.) *tincf=tincfguess;
  printf("time increment for the CFD-calculations = %e\n\n",*tincf);

  ttimef=*ttime;
  timef=*time-*dtime;
  dtimef=*tincf;

  do{

      iit++;

      printf("iteration = %d\n",iit);

      timef+=dtimef;
      if((*time<timef)&&(*nmethod==4)){
	  dtimef-=(timef-*time);
	  timef=*time;
	  last=1;
	  beta=dtimef/(*tincf);
	  a1=(2.+beta)/(1.+beta);
	  a2=-(1.+beta)/beta;
	  a3=1./(beta*(1.+beta));
      }

      /* conditions for transient calculations */
     
      if(*nmethod==4){

          /* boundary conditions at end of fluid increment */

	  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,
             xloadold,xload,xloadact,iamload,nload,ibody,xbody,nbody,
             xbodyold,xbodyact,t1old,t1,t1act,iamt1,nk,amta,
             namta,nam,ampli,&timef,&reltimef,&ttimef,&dtimef,ithermal,nmethod,
             xbounold,xboun,xbounact,iamboun,nboun,
             nodeboun,ndirboun,nodeforc,ndirforc,istep,iinc,
	     co,v,itg,ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
             ntrans,trab,inotr,vold,integerglob,doubleglob,tieset,istartset,
             iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));

	  /* body forces (gravity, centrifugal and Coriolis forces) */
      
	  if(*nbody>0){
	      FORTRAN(calcbody,(nef,body,ipobody,ibody,xbody,coel,vel,lakonf,
				nactdohinv));
	  }

      }else if(icent==1){
	  
	  /* body forces (gravity, centrifugal and Coriolis forces;
             only if centrifugal forces are active => the ensuing
             Coriolis forces depend on the actual velocity) */
	  
	  FORTRAN(calcbody,(nef,body,ipobody,ibody,xbody,coel,vel,lakonf,
				nactdohinv));
      }

      /* updating of the material properties */

      if(*ithermal>0){
	  
	  /* calculating the density at the element centers */
	  
	  FORTRAN(calcrhoel,(nef,vel,rhcon,nrhcon,ielmatf,ntmat_,
			     ithermal,mi));
	  
	  /* calculating the density at the face centers */
	  
	  FORTRAN(calcrhofa,(nface,vfa,rhcon,nrhcon,ielmatf,ntmat_,
			     ithermal,mi,ielfa));
	  
	  /* calculating the dynamic viscosity at the face centers */
	  
	  FORTRAN(calcumfa,(nface,vfa,shcon,nshcon,ielmatf,ntmat_,
			     ithermal,mi,ielfa,umfa));

          /* calculating the dynamic viscosity at the element centers */
	  
	  FORTRAN(calcumel,(nef,vel,shcon,nshcon,ielmatf,ntmat_,
			    ithermal,mi,umel));
	  
	  /* calculating the heat conduction at the face centers */
	  
	  FORTRAN(calchcfa,(nface,vfa,cocon,ncocon,ielmatf,ntmat_,
			    mi,ielfa,hcfa));

	  /* calculating the specific heat at constant pressure at the 
             element centers (secant value) */
	  
	  FORTRAN(calccpel,(nef,vel,shcon,nshcon,ielmatf,ntmat_,
			    mi,cpel,physcon));

	  /* calculating the specific heat at constant pressure at the 
             face centers (secant value) */
	  
	  FORTRAN(calccpfa,(nface,vfa,shcon,nshcon,ielmatf,ntmat_,
			    mi,ielfa,cpfa,physcon));

      }
	  
      /* calculating the gradient of the velocity at the element
	 centers */
      
      DMEMSET(gradvel,0,9**nef,0.);
      FORTRAN(calcgradvel,(nef,lakonf,ipnei,vfa,area,xxn,gradvel,neifa,
                         volume));
      
      /* extrapolating the gradient of the velocity from the element
	 centers to the face centers */
      
      FORTRAN(extrapolate_gradvel,(nface,ielfa,xrlfa,gradvel,gradvfa));
      
      /* calculate gamma (Ph.D. Thesis Jasak) */

      /*       betam=0.1;
      NNEW(gamma,double,*nface);
      FORTRAN(calcgamma,(nface,ielfa,vel,gradvfa,gamma,xlet,xxn,xxj,
      ipnei,&betam,nef));*/

      /* filling the lhs and rhs's for the balance of momentum
	 equations */
      
      DMEMSET(adv,0,neq,0.);
      DMEMSET(auv,0,2*nzs,0.);
      DMEMSET(bv,0,3*neq,0.);
      
      FORTRAN(mafillv,(nef,ipnei,neifa,neiel,vfa,xxn,area,
		       auv,adv,jq,irow,&nzs,bv,vel,cosa,umfa,xlet,xle,gradvfa,xxi,
		       body,volume,&compressible,ielfa,lakonf,ifabou,nbody,&neq,
		       &dtimef,velo,veloo,sel,xrlfa,gamma,xxj,nactdohinv,&a1,
                       &a2,&a3));
//      SFREE(gamma);
      
      /* LU decomposition (asymmetric system) */
      
/*      inputformat=1;
      symmetryflag=2;
      
      if(*isolver==0){
#ifdef SPOOLES
	  spooles_factor(adv,auv,adb,aub,&sigma,icol,irow,&neq,&nzs,
			 &symmetryflag,&inputformat,&nzs);
#else
	  printf("*ERROR in compfluid: the SPOOLES library is not linked\n\n");
	  FORTRAN(stop,());
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_factor(adv,auv,adb,aub,&sigma,icol,irow,&neq,&nzs,
			 &symmetryflag,&inputformat,jq,&nzs);
#else
	  printf("*ERROR in compfluid: the PARDISO library is not linked\n\n");
	  FORTRAN(stop,());
#endif
}*/
      
      /* solving the system of equations (for x, y and z
	 separately */
      
/*      for(i=0;i<3;i++){
	  
	  if(*isolver==0){
#ifdef SPOOLES
	      spooles_solve(&bv[i*neq],&neq);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	      pardiso_solve(&bv[i*neq],&neq,&symmetryflag);
#endif
	  }
	  }*/
      
      /* free memory */
      
/*      if(*isolver==0){
#ifdef SPOOLES
	  spooles_cleanup();
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_cleanup(&neq,&symmetryflag);
#endif
}*/

      memcpy(&auv[2*nzs],adv,sizeof(double)*neq);
      nelt=2*nzs+neq;
      isym=0;
      itol=0;
      tol=0.;
      itmax=0;
      iunit=0;
      lrgw=131+16*neq;
      NNEW(rgwk,double,lrgw);
      NNEW(igwk,ITG,20);
      igwk[0]=10;
      igwk[1]=10;
      igwk[2]=0;
      igwk[3]=1;
      igwk[4]=10;
      ligw=20;
      for(i=0;i<neq;i++){rwork[i]=1./adv[i];}
      FORTRAN(predgmres,(&neq,&bv[0],&vel[neq],&nelt,irow,jq,auv,
		      &isym,&itol,&tol,&itmax,&iter,
                      &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
                      &ligw,rwork,iwork));
//	      printf(" err1=%e\n",err);
//      memcpy(&bv[0],&vel[neq],sizeof(double)*neq);
      FORTRAN(predgmres,(&neq,&bv[neq],&vel[2*neq],&nelt,irow,jq,auv,
		      &isym,&itol,&tol,&itmax,&iter,
                      &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
                      &ligw,rwork,iwork));
//	      printf(" err2=%e\n",err);
//      memcpy(&bv[neq],&vel[2*neq],sizeof(double)*neq);
      FORTRAN(predgmres,(&neq,&bv[2*neq],&vel[3*neq],&nelt,irow,jq,auv,
		      &isym,&itol,&tol,&itmax,&iter,
                      &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
                      &ligw,rwork,iwork));
//	      printf(" err3=%e\n",err);
//      memcpy(&bv[2*neq],&vel[3*neq],sizeof(double)*neq);
      SFREE(rgwk);SFREE(igwk);
      
      /* storing the solution into vel */
      
//      FORTRAN(calcvel,(ne,nactdoh,vel,bv,&neq,ne));

      /* calculating the pressure gradient at the element
         centers */

      DMEMSET(gradpel,0,3**nef,0.);
      FORTRAN(calcgradpel,(nef,lakonf,ipnei,vfa,area,xxn,gradpel,neifa,
			       volume));
      
      for(iitf=0;iitf<2;iitf++){

	  memcpy(&hel[0],&sel[0],sizeof(double)*(3*neq));
	  
          /* completing hel with the neighboring velocity contributions */
	  
          FORTRAN(complete_hel,(&neq,&vel[neq],hel,adv,auv,jq,irow,&nzs));
//          FORTRAN(complete_hel,(&neq,bv,hel,adv,auv,jq,irow,&nzs));
	  
	  /* generating ad and h at the face centers (advfa and hfa) */
	  
	  FORTRAN(extrapolate_ad_h,(nface,ielfa,xrlfa,adv,advfa,hel,hfa));
	  
	  /* calculating the lhs and rhs of the equation system to determine
	     p (balance of mass) */
	  
	  DMEMSET(b,0,neq,0.);

	  if(iitf==0){

              /* first iteration: calculating both lhs and rhs */

	      DMEMSET(ad,0,neq,0.);
	      DMEMSET(au,0,nzs,0.);
	      
	      FORTRAN(mafillp,(nef,lakonf,ipnei,neifa,neiel,vfa,area,
			       advfa,xlet,cosa,volume,au,ad,jq,irow,ap,
                               ielfa,ifabou,xle,b,xxn,&compressible,&neq,
                               &nzs,hfa,gradpel,bp,xxi,neij,xlen,cosb));

	      FORTRAN(convert2slapcol,(au,ad,irow,ia,jq,ja,&nzs,&neq,aua));
	      
	      /* LU decomposition of the p system (symmetric system) */
	      
/*	      inputformat=0;
	      symmetryflag=0;
	      
	      if(*isolver==0){
#ifdef SPOOLES
		  spooles_factor(ad,au,adb,aub,&sigma,icol,irow,&neq,&nzs,
				 &symmetryflag,&inputformat,&nzs);
#else
		  printf("*ERROR in compfluid: the SPOOLES library is not linked\n\n");
		  FORTRAN(stop,());
#endif
	      }
	      else if(*isolver==7){
#ifdef PARDISO
		  pardiso_factor(ad,au,adb,aub,&sigma,icol,irow,&neq,&nzs,
				 &symmetryflag,&inputformat,jq,&nzs);
#else
		  printf("*ERROR in compfluid: the PARDISO library is not linked\n\n");
		  FORTRAN(stop,());
#endif
}*/
	  }else{

	      /* calculating the pressure gradient at the element
		 centers */

	      DMEMSET(gradpel,0,3**nef,0.);
	      FORTRAN(calcgradpel,(nef,lakonf,ipnei,vfa,area,xxn,gradpel,neifa,
	      volume));

	      /* second, third.. iteration: calculate the rhs only */
	      
	      FORTRAN(rhsp,(nef,lakonf,nactdoh,ipnei,neifa,neiel,vfa,area,
		      advfa,xlet,cosa,volume,au,ad,jq,irow,ap,ielfa,ifabou,xle,
		      b,xxn,&compressible,&neq,&nzs,hfa,bp,neij,xxi,
                      gradpel,xlen));

	  }
	  
	  /* solving the system of equations for p  */
	  
/*	  if(*isolver==0){
#ifdef SPOOLES
	      spooles_solve(b,&neq);
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	      pardiso_solve(b,&neq,&symmetryflag);
#endif
}*/

          /* dslugm; attention: convert2slapcol changes au! */

//	  FORTRAN(convert2slapcol,(au,ad,irow,ia,jq,ja,&nzs,&neq,aua));

	  nelt=nzs+neq;
	  isym=1;
	  nsave=10;
	  itol=0;
	  tol=0.;
	  itmax=110;
	  iunit=0;
	  lenw=131+17*neq+2*nelt;
	  NNEW(rgwk,double,lenw);
	  leniw=32+4*neq+2*nelt;
	  NNEW(igwk,ITG,leniw);
	  
	  FORTRAN(dslugm,(&neq,&b[0],&vel[4*neq],&nelt,ia,ja,aua,
			  &isym,&nsave,&itol,&tol,&itmax,&iter,
			  &err,&ierr,&iunit,rgwk,&lenw,igwk,&leniw));
	  SFREE(rgwk);SFREE(igwk);

	  if(ierr!=0){
	      printf("       error message from dslugm: ierr = %d\n",ierr);
	      printf(" err=%e\n",err);
	      }
	  
	  /* storing the solution p into vel(4,*) */
	  
//	  FORTRAN(calcpel,(nef,nactdoh,vel,b,nef));
	  
	  /* extrapolation of the pressure at the element centers
	     to the face centers */
	  
	  FORTRAN(extrapolate_pel,(nface,ielfa,xrlfa,vel,vfa,ifabou,
				   xbounact,nef));

	  /* applying MPC's to the faces */
	  
	  if(*nmpc>0){
	      is=4;ie=4;
	      FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfa,coefmpc,
				nodempc,ipnei,neifa,labmpc,xbounact,nactdoh));
	  }
	  
	  /* correction of the velocity at the element centers due
             to the pressure change */

	  FORTRAN(correctvel,(hel,adv,vfa,ipnei,area,&vel[neq],xxn,neifa,
			      lakonf,nef,&neq));

      }
      
      /* free memory */
      
      /*   if(*isolver==0){
#ifdef SPOOLES
	  spooles_cleanup();
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_cleanup(&neq,&symmetryflag);
#endif
}*/
      
      if(*ithermal>0){

	  /* adding the velocity correction at the face centers
	     due to the balance of mass =>
	     the resulting mass flux is correct,
	     the face velocity vectors do not have to be correct
	     only needed if extra equations (temperature,
	     turbulence have to be solved)  */

	  FORTRAN(correctvfa,(nface,ielfa,area,vfa,ap,bp,xxn,
			      ifabou,ipnei,nef,neifa,hfa,vel,xbounact));

          /* calculating the temperature gradient at the element
             centers */

	  DMEMSET(gradtel,0,3**nef,0.);
	  FORTRAN(calcgradtel,(nef,lakonf,ipnei,vfa,area,xxn,gradtel,neifa,
			       volume));

	  /* extrapolating the temperature gradient from the element
             centers to the face centers */

	  FORTRAN(extrapolate_gradtel,(nface,ielfa,xrlfa,gradtel,gradtfa));
      
      /* calculate gammat (Ph.D. Thesis Jasak) */

//	  betam=0.1;
//	  NNEW(gammat,double,*nface);
//	  FORTRAN(calcgammat,(nface,ielfa,vel,gradtfa,gammat,xlet,xxn,xxj,
//			      ipnei,&betam,ne));

          /* calculating the lhs and rhs of the energy equation */

	  DMEMSET(ad,0,neq,0.);
	  DMEMSET(au,0,2*nzs,0.);
	  DMEMSET(b,0,neq,0.);

	  FORTRAN(mafillt,(nef,ipnei,neifa,neiel,vfa,xxn,area,
	       au,ad,jq,irow,&nzs,b,vel,umel,xlet,xle,gradtfa,xxi,
	       body,volume,&compressible,ielfa,lakonf,ifabou,nbody,&neq,
	       &dtimef,velo,veloo,cpfa,hcfa,cpel,gradvel,xload,gammat,
	       xrlfa,xxj,nactdohinv,&a1,&a2,&a3));

//	  SFREE(gammat);
	  
          /* solving the asymmetric system of equations */

	  /* inputformat=1;
	  symmetryflag=2;
	  
	  if(*isolver==0){
#ifdef SPOOLES
	      spooles(ad,au,adb,aub,&sigma,b,icol,irow,&neq,&nzs,
		      &symmetryflag,&inputformat,&nzs);
#else
	      printf("*ERROR in compfluid: the SPOOLES library is not linked\n\n");
	      FORTRAN(stop,());
#endif
	  }
	  else if(*isolver==7){
#ifdef PARDISO
	      pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,&neq,&nzs,
			   &symmetryflag,&inputformat,jq,&nzs);
#else
	      printf("*ERROR in compfluid: the PARDISO library is not linked\n\n");
	      FORTRAN(stop,());
#endif
}*/
      
      /* free memory */
      
      /*    if(*isolver==0){
#ifdef SPOOLES
	  spooles_cleanup();
#endif
      }
      else if(*isolver==7){
#ifdef PARDISO
	  pardiso_cleanup(&neq,&symmetryflag);
#endif
}*/

	  memcpy(&au[2*nzs],ad,sizeof(double)*neq);
	  nelt=2*nzs+neq;
	  isym=0;
	  itol=0;
	  tol=0.;
	  itmax=0;
	  iunit=0;
	  lrgw=131+16*neq;
	  NNEW(rgwk,double,lrgw);
	  NNEW(igwk,ITG,20);
	  igwk[0]=10;
	  igwk[1]=10;
	  igwk[2]=0;
	  igwk[3]=1;
	  igwk[4]=10;
	  ligw=20;
	  for(i=0;i<neq;i++){rwork[i]=1./ad[i];}
	  FORTRAN(predgmres,(&neq,&b[0],&vel[0],&nelt,irow,jq,au,
			     &isym,&itol,&tol,&itmax,&iter,
			     &err,&ierr,&iunit,sb,sx,rgwk,&lrgw,igwk,
			     &ligw,rwork,iwork));
	  SFREE(rgwk);SFREE(igwk);
	  if(ierr>0){
	      printf("*WARNING in compfluid: error message from predgmres=%e\n",err);
	  }
	  
	  /* storing the solution t into vel(0,*) */
	  
//	  FORTRAN(calctel,(ne,nactdoh,vel,b,ne));

	  /* extrapolation of the temperature at the element centers
	     to the face centers */

	  FORTRAN(extrapolate_tel,(nface,ielfa,xrlfa,vel,vfa,
				   ifabou,xbounact,ipnei,nef));

	  /* applying MPC's to the faces */
	  
	  if(*nmpc>0){
	      is=0;ie=0;
	      FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfa,coefmpc,
				nodempc,ipnei,neifa,labmpc,xbounact,nactdoh));
	  }

      }
      
      /* storing the solution into vel */
      
//      FORTRAN(calcvel,(ne,nactdoh,vel,bv,&neq,ne));

      /* extrapolating the velocity from the elements centers to the face
	 centers, thereby taking the boundary conditions into account */
      
      FORTRAN(extrapolate_vel,(nface,ielfa,xrlfa,vel,vfa,
			       ifabou,xbounact,ipnei,nef));

      /* applying MPC's to the faces */
      
      if(*nmpc>0){
	  is=1;ie=3;
	  FORTRAN(applympc,(nface,ielfa,&is,&ie,ifabou,ipompc,vfa,coefmpc,
			    nodempc,ipnei,neifa,labmpc,xbounact,nactdoh));
      }

//      FORTRAN(writevfa,(vfa,nface,nactdohinv,ielfa));
      
      
      if(((iit/jout[1])*jout[1]==iit)||(iconvergence==1)||
	 (iit==jmax[1])){

	  /* calculating the stress and the heat flow at the
             integration points, if requested */

	  if((strcmp1(&filab[3306],"SF  ")==0)||
             (strcmp1(&filab[3480],"SVF ")==0))isti=1;
          if(strcmp1(&filab[3393],"HFLF")==0)iqfx=1;
	  for(i=0;i<*nprint;i++){
	      if(strcmp1(&prlab[6*i],"SVF")==0) isti=1;
	      if(strcmp1(&prlab[6*i],"HFLF")==0)iqfx=1;
	  }

          /* calculating the heat conduction at the element centers */

	  if(iqfx==1){
	      NNEW(hcel,double,*nef);
	      FORTRAN(calchcel,(vel,cocon,ncocon,ielmatf,ntmat_,mi,
			       hcel,nef));
	  }

	  /* calculating the stress and/or the heat flux at the
             element centers */

	  if((isti==1)||(iqfx==1)){
	      FORTRAN(calcstressheatflux,(sti,umel,gradvel,qfx,hcel,
					  gradtel,nef,&isti,&iqfx,mi));
	      if(iqfx==1)SFREE(hcel);
	  }
 
          /* extrapolating the stresses */

	  if((strcmp1(&filab[3306],"SF  ")==0)||
             (strcmp1(&filab[3480],"SVF ")==0)){
	      nfield=6;
	      ndim=6;
	      if((*norien>0)&&
                 ((strcmp1(&filab[3311],"L")==0)||(strcmp1(&filab[3485],"L")==0))){
		  iorienglob=1;
	      }else{
		  iorienglob=0;
	      }
	      strcpy1(&cflag[0],&filab[2962],1);
	      NNEW(stn,double,6**nk);
	      FORTRAN(extrapolate,(sti,stn,ipkonf,inum,kon,lakonf,
		      &nfield,nk,nef,mi,&ndim,orab,ielorienf,co,&iorienglob,
		      cflag,nelemload,nload,nodeboun,nboun,ndirboun,
		      vold,ithermal,&force,&icfdout,ielmatf,thicke,filab));
	  }

	  /* extrapolating the heat flow */

	  
	  if(strcmp1(&filab[3393],"HFLF")==0){
	      nfield=3;
	      ndim=3;
	      if((*norien>0)&&(strcmp1(&filab[3398],"L")==0)){
		  iorienglob=1;
	      }else{
		  iorienglob=0;
	      }
	      strcpy1(&cflag[0],&filab[3049],1);
	      NNEW(qfn,double,3**nk);
	      FORTRAN(extrapolate,(qfx,qfn,ipkonf,inum,kon,lakonf,
		      &nfield,nk,nef,mi,&ndim,orab,ielorienf,co,&iorienglob,
		      cflag,nelemload,nload,nodeboun,nboun,ndirboun,
		      vold,ithermal,&force,&icfdout,ielmatf,thicke,filab));
	  }
	  
	  /* extrapolating the facial values of the static temperature 
             and/or the velocity and/or the static pressure to the nodes */
	  
	  FORTRAN(extrapolatefluid,(nk,iponofa,inofa,inum,vfa,vold,ielfa,
                      ithermal));

          /* storing the results in dat-format */

	  ptimef=ttimef+timef;
	  FORTRAN(printoutfluid,(set,nset,istartset,iendset,ialset,nprint,
				 prlab,prset,v,t1,fn,ipkonf,lakonf,sti,eei,
                                 xstate,ener,mi,&nstate_,ithermal,co,kon,qfx,
                                 &ptimef,trab,inotr,ntrans,orab,ielorienf,
                                 norien,nk,nef,inum,filab,vold,ielmatf,
                                 thicke,eme,vcontu,physcon,nactdoh,
                                 ielpropf,prop));
	  
	  /* storing the results in frd-format */
	  
	  FORTRAN(frdfluid,(co,nk,kon,ipkonf,lakonf,nef,vold,&kode,&timef,ielmatf,
			    matname,filab,inum,ntrans,inotr,trab,mi,istep,
                            stn,qfn,nactdohinv));

	  if((strcmp1(&filab[3306],"SF  ")==0)||
             (strcmp1(&filab[3480],"SVF ")==0)){SFREE(stn);}
	  if(strcmp1(&filab[3393],"HFLF")==0){SFREE(qfn);}

      }
      
      if(iit==jmax[1]){
	  printf("*INFO: maximum number of fluid increments reached\n\n");
	  FORTRAN(stop,());
      }
      if(last==1){FORTRAN(stop,());}
      
      memcpy(&veloo[0],&velo[0],sizeof(double)*6**nef);
      memcpy(&velo[0],&vel[0],sizeof(double)*6**nef);
      
  }while(1);
  
  FORTRAN(closefilefluid,());

  SFREE(irow);SFREE(ia);SFREE(icol);SFREE(jq);SFREE(ja);
//SFREE(nactdoh);
  
  SFREE(coel);SFREE(cosa);SFREE(xxn);SFREE(xxi);SFREE(xle);SFREE(xlen);
  SFREE(xlet);SFREE(cofa);SFREE(area);SFREE(xrlfa);SFREE(volume);
  SFREE(cosb);SFREE(xxj);

  SFREE(ifabou);SFREE(umfa);SFREE(umel);

  SFREE(gradvel);SFREE(gradvfa);SFREE(au);SFREE(ad);SFREE(b);SFREE(advfa);
  SFREE(ap);SFREE(bp);SFREE(gradpel);SFREE(rwork);SFREE(aua);
  SFREE(hfa);SFREE(hel);SFREE(auv);SFREE(adv);SFREE(bv);SFREE(sel);

  if(*ithermal>0){
      SFREE(gradtel);SFREE(gradtfa);SFREE(hcfa);SFREE(cpel);SFREE(cpfa);
  }

  SFREE(inum);SFREE(v);SFREE(velo);SFREE(veloo);

  SFREE(iponofa);SFREE(inofa);

  if(*nbody>0) SFREE(body);

  *ithermal=ithermalref;

  return;
  
} 
