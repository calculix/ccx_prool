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

static char *sideload1;

static ITG *kontri1,*nloadtr1,*idist=NULL,*ntrit1,*mi1,*jqrad1,
    *irowrad1,*nzsrad1,num_cpus,*ntri1,*ntr1;

static double *vold1,*co1,*pmid1,*e11,*e21,*e31,*adview=NULL,*auview=NULL,*dist=NULL,
    *area1,sidemean1;

void radflowload(ITG *itg,ITG *ieg,ITG *ntg,ITG *ntr,double *adrad,
                 double *aurad,double *bcr,ITG *ipivr,
                 double *ac,double *bc,ITG *nload,char *sideload,
                 ITG *nelemload,double *xloadact,char *lakon,ITG *ipiv,
                 ITG *ntmat_,double *vold,double *shcon,
                 ITG *nshcon,ITG *ipkon,ITG *kon,double *co,
                 ITG *kontri,
                 ITG *ntri,ITG *nloadtr,double *tarea,double *tenv,
                 double *physcon,double *erad,double **adviewp, 
                 double **auviewp,
                 ITG *nflow,ITG *ikboun,
                 double *xbounact,ITG *nboun,ITG *ithermal,
                 ITG *iinc,ITG *iit,double *cs, ITG *mcs, ITG *inocs, 
                 ITG *ntrit,ITG *nk, double *fenv,ITG *istep,double *dtime,
                 double *ttime,double *time,ITG *ilboun,ITG *ikforc,
                 ITG *ilforc,double *xforcact,ITG *nforc,double *cam,
                 ITG *ielmat,ITG *nteq,double *prop,ITG *ielprop,ITG *nactdog,
                 ITG *nacteq,ITG *nodeboun,ITG *ndirboun,
                 ITG *network, double *rhcon, ITG *nrhcon, ITG *ipobody,
                 ITG *ibody, double *xbodyact, ITG *nbody,ITG *iviewfile,
                 char *jobnamef, double *ctrl, double *xloadold,
                 double *reltime, ITG *nmethod, char *set, ITG *mi,
		 ITG * istartset,ITG* iendset,ITG *ialset,ITG *nset,
                 ITG *ineighe, ITG *nmpc, ITG *nodempc,ITG *ipompc,
                 double *coefmpc,char *labmpc, ITG *iemchange,ITG *nam, 
                 ITG *iamload,ITG *jqrad,ITG *irowrad,ITG *nzsrad,
                 ITG *icolrad,ITG *ne,ITG *iaxial,double *qa){
  
  /* network=0: purely thermal
     network=1: general case (temperatures, fluxes and pressures unknown)
     network=2: purely aerodynamic, i.e. only fluxes and pressures unknown */

  ITG nhrs=1,info=0,i,j,iin=0,icntrl,icutb=0,iin_abs=0,mt=mi[1]+1,im,
      symmetryflag=2,inputformat=1,node,channel,*ithread=NULL;

  static ITG ifactorization=0;

  double uamt=0,uamf=0,uamp=0,camt[2],camf[2],camp[2],
    cam1t=0.,cam1f=0.,cam1p=0.,sidemean,
    cam2t=0.,cam2f=0.,cam2p=0.,dtheta=1.,*v=NULL,cama[2],cam1a=0.,
    cam2a=0.,uama=0.,vamt=0.,vamf=0.,vamp=0.,vama=0.,cam0t=0.,cam0f=0.,
    cam0p=0.,cam0a=0.,sigma=0.,*adbrad=NULL,*aubrad=NULL,*q=NULL,
    *area=NULL,*pmid=NULL,*e1=NULL,*e2=NULL,*e3=NULL;
  
  adview=*adviewp;auview=*auviewp;

  /* check whether there are any gas temperature nodes; this check should
     NOT be done on nteq, since also for zero equations the temperature
     of the gas nodes with boundary conditions must be stored in v
     (in initialgas) */ 

  NNEW(v,double,mt**nk);

  /* gas networks */

  if(*ntg!=0) {
      icntrl=0;
      while(icntrl==0) {
	  
	  if(iin==0){
	      
	      for(i=0;i<mt**nk;i++) v[i]=vold[i];

              /* initialization pressurized flow 
                 (no free surface: gas networks or
                  water networks with fully wetted perimeter*/

	      FORTRAN(initialnet,(itg,ieg,ntg,ac,bc,lakon,v,
                           ipkon,kon,nflow,
			   ikboun,nboun,prop,ielprop,nactdog,ndirboun,
			   nodeboun,xbounact,ielmat,ntmat_,shcon,nshcon,
			   physcon,ipiv,nteq,rhcon,nrhcon,ipobody,ibody,
			   xbodyact,co,nbody,network,&iin_abs,vold,set,
			   istep,iit,mi,ineighe,ilboun,&channel,iaxial,
			   nmpc,labmpc,ipompc,nodempc,coefmpc));
      
              /* initialization for channels with free surface */

	      if(channel==1){
		  FORTRAN(initialchannel,(itg,ieg,ntg,ac,bc,lakon,v,
                           ipkon,kon,nflow,
			   ikboun,nboun,prop,ielprop,nactdog,ndirboun,
			   nodeboun,xbounact,ielmat,ntmat_,shcon,nshcon,
			   physcon,ipiv,nteq,rhcon,nrhcon,ipobody,ibody,
			   xbodyact,co,nbody,network,&iin_abs,vold,set,
			   istep,iit,mi,ineighe,ilboun));
	      }

              /* storing the residual in the rhs vector */

	      FORTRAN(resultnet,(itg,ieg,ntg,bc,nload,sideload,
			  nelemload,xloadact,
			  lakon,ntmat_,v,shcon,nshcon,ipkon,kon,co,nflow,
			  iinc,istep,dtime,ttime,time,
			  ikforc,ilforc,xforcact,
                          nforc,ielmat,nteq,prop,ielprop,nactdog,nacteq,&iin,
			  physcon,camt,camf,camp,rhcon,nrhcon,ipobody,
			  ibody,xbodyact,nbody,&dtheta,vold,xloadold,
			  reltime,nmethod,set,mi,ineighe,cama,&vamt,
			  &vamf,&vamp,&vama,nmpc,nodempc,ipompc,coefmpc,
			  labmpc,iaxial));
	  }
	  
	  iin++;
	  iin_abs++;
	  printf("      gas iteration %" ITGFORMAT " \n \n",iin);
	  
          /* filling the lhs matrix */
         
	  FORTRAN(mafillnet,(itg,ieg,ntg,ac,nload,sideload,
			     nelemload,xloadact,lakon,ntmat_,v,
			     shcon,nshcon,ipkon,kon,co,nflow,iinc,
			     istep,dtime,ttime,time,
			     ielmat,nteq,prop,ielprop,nactdog,nacteq,
			     physcon,rhcon,nrhcon,ipobody,ibody,xbodyact,
			     nbody,vold,xloadold,reltime,nmethod,set,mi,
                             nmpc,nodempc,ipompc,coefmpc,labmpc,iaxial));
	  
          /* solving the system of equations */

	  if(*nteq>0){
	      FORTRAN(dgesv,(nteq,&nhrs,ac,nteq,ipiv,bc,nteq,&info)); 
	  }

	    /*spooles(ac,au,adb,aub,&sigma,bc,icol,irow,nteq,nteq,
	      &symmetryflag,&inputformat);*/
	  
	  if (info!=0) {
	      printf(" *WARNING in radflowload: singular matrix\n");
	    
	      FORTRAN(mafillnet,(itg,ieg,ntg,ac,nload,sideload,
				 nelemload,xloadact,lakon,ntmat_,v,
				 shcon,nshcon,ipkon,kon,co,nflow,iinc,
				 istep,dtime,ttime,time,
				 ielmat,nteq,prop,ielprop,nactdog,nacteq,
				 physcon,rhcon,nrhcon,ipobody,ibody,xbodyact,
				 nbody,vold,xloadold,reltime,nmethod,set,mi,
                                 nmpc,nodempc,ipompc,coefmpc,labmpc,iaxial));
	    
	      FORTRAN(equationcheck,(ac,nteq,nactdog,itg,ntg,nacteq,network));
	    
	      iin=0;

	  }
	  else {

              /* storing the residual in the rhs vector */

	      FORTRAN(resultnet,(itg,ieg,ntg,bc,nload,sideload,nelemload,
	       xloadact,lakon,ntmat_,v,shcon,nshcon,ipkon,kon,co,
	       nflow,iinc,istep,dtime,ttime,time,ikforc,ilforc,xforcact,
	       nforc,ielmat,nteq,prop,ielprop,nactdog,nacteq,
	       &iin,physcon,camt,camf,camp,rhcon,nrhcon,ipobody,
	       ibody,xbodyact,nbody,&dtheta,vold,xloadold,
	       reltime,nmethod,set,mi,ineighe,cama,&vamt,
	       &vamf,&vamp,&vama,nmpc,nodempc,ipompc,coefmpc,labmpc,
               iaxial));

              /* printing the largest corrections */
	    
	      if(*network!=2){ 
		  cam2t=cam1t;
		  cam1t=cam0t;
		  cam0t=camt[0];
		  if (camt[0]>uamt) {uamt=camt[0];}
		  printf
                    ("      largest increment of gas temperature= %e\n",uamt);
		  if((ITG)camt[1]==0){
		      printf
		      ("      largest correction to gas temperature= %e\n",
                       camt[0]);
		  }else{
		      printf
		      ("      largest correction to gas temperature= %e in node %" ITGFORMAT "\n",
                       camt[0],(ITG)camt[1]);
		  }
	      }
	      
	      if(*network!=0){
		  cam2f=cam1f;
		  cam1f=cam0f;
		  cam0f=camf[0];
		  if (camf[0]>uamf) {uamf=camf[0];}
		  printf("      largest increment of gas massflow= %e\n",uamf);
		  if((ITG)camf[1]==0){
		      printf("      largest correction to gas massflow= %e\n",
			 camf[0]);
		  }else{
		      printf("      largest correction to gas massflow= %e in node %" ITGFORMAT "\n",
			 camf[0],(ITG)camf[1]);
		  }
		  
		  cam2p=cam1p;
		  cam1p=cam0p;
		  cam0p=camp[0];
		  if (camp[0]>uamp) {uamp=camp[0];}
		  printf("      largest increment of gas pressure= %e\n",uamp);
		  if((ITG)camp[1]==0){
		      printf("      largest correction to gas pressure= %e\n",
                         camp[0]);
		  }else{
		      printf("      largest correction to gas pressure= %e in node %" ITGFORMAT "\n",
                         camp[0],(ITG)camp[1]);
		  }
		  
		  cam2a=cam1a;
		  cam1a=cam0a;
		  cam0a=cama[0];
		  if (cama[0]>uama) {uama=cama[0];}
		  printf("      largest increment of geometry= %e\n",uama);
		  if((ITG)cama[1]==0){
		      printf("      largest correction to geometry= %e\n",
                         cama[0]);
		  }else{
		      printf("      largest correction to geometry= %e in node %" ITGFORMAT "\n",
                         cama[0],(ITG)cama[1]);
		  }
	      }	      
	  }
	  
	  printf("\n");
	  
	  /* for purely thermal calculations no iterations are
	     deemed necessary */
	  
	  if(*network==0) {icntrl=1;}
	  else {

              /* check the convergence */

	      checkconvnet(&icutb,&iin,&uamt,&uamf,&uamp,
		 &cam1t,&cam1f,&cam1p,&cam2t,&cam2f,&cam2p,&cam0t,&cam0f,
		 &cam0p,&icntrl,&dtheta,ctrl,&uama,&cam1a,&cam2a,&cam0a,
		 &vamt,&vamf,&vamp,&vama,qa);
	  }
      }

      /* storing network output as boundary conditions for
         the structure */

      FORTRAN(flowresult,(ntg,itg,cam,vold,v,nload,sideload,
	      nelemload,xloadact,nactdog,network,mi,ne,ipkon,lakon,kon));

      /* extra output for hydraulic jump (fluid channels) */

#ifdef NETWORKOUT
      if(*network!=0){
	FORTRAN(flowoutput,(itg,ieg,ntg,nteq,bc,lakon,ntmat_,
			    v,shcon,nshcon,ipkon,kon,co,nflow, dtime,ttime,time,
			    ielmat,prop,ielprop,nactdog,nacteq,&iin,physcon,
			    camt,camf,camp,&uamt,&uamf,&uamp,rhcon,nrhcon,
			    vold,jobnamef,set,istartset,iendset,ialset,nset,
                            mi,iaxial));
      }
#endif
  }
      
  /* radiation */

  if(*ntr>0){
      
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
      
      envloc = getenv("CCX_NPROC_VIEWFACTOR");
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
      
      if(*ntr<num_cpus) num_cpus=*ntr;
      
      pthread_t tid[num_cpus];
      
      /*the default sink temperature is updated at the start of each
	increment */
     
      for(i=0;i<*ntr;i++){
	  node=nelemload[2*nloadtr[i]-1];
	  if(node!=0){
	      tenv[i]=vold[mt*(node-1)]-physcon[0];
	  }else if(*iit<=0){
	      tenv[i]=xloadact[2*nloadtr[i]-1]-physcon[0];
	  }
      }
     
/*     for pure thermal steps the viewfactors have to be
       calculated only once, for thermo-mechanical steps
       (ithermal=3) they are recalculated in each iteration
       unless they are read from file */

      if(((*ithermal==3)&&(*iviewfile>=0))||(*iit==-1)){
	  if(*iviewfile<0){

              /* reading viewfactors from file */

	      FORTRAN(readview,(ntr,adview,auview,fenv,nzsrad,ithermal,
				jobnamef));

	  }else{

	      /* determining geometric data to calculate the viewfactors */

	      NNEW(area,double,*ntrit);
	      NNEW(pmid,double,3**ntrit);
	      NNEW(e1,double,3**ntrit);
	      NNEW(e2,double,3**ntrit);
	      NNEW(e3,double,4**ntrit);

	      FORTRAN(geomview,(vold,co,pmid,e1,e2,e3,kontri,area,
				cs,mcs,inocs,ntrit,nk,mi,&sidemean));

	      RENEW(adview,double,num_cpus**ntr);
	      RENEW(auview,double,num_cpus*2**nzsrad);
	      
	      NNEW(dist,double,num_cpus**ntrit);
	      NNEW(idist,ITG,num_cpus**ntrit);

	      DMEMSET(adview,0,num_cpus**ntr,0.);
	      DMEMSET(auview,0,num_cpus*2**nzsrad,0.);

	      sideload1=sideload;vold1=vold;co1=co;pmid1=pmid;
	      e11=e1;e21=e2;e31=e3;kontri1=kontri;ntr1=ntr;
              nloadtr1=nloadtr;area1=area;ntri1=ntri;
              ntrit1=ntrit;mi1=mi;jqrad1=jqrad;irowrad1=irowrad;
              nzsrad1=nzsrad;sidemean1=sidemean;

	      /* calculating the viewfactors */
	      
	      printf(" Using up to %" ITGFORMAT " cpu(s) for the viewfactor calculation.\n\n", num_cpus);
  
	      /* create threads and wait */
	      
	      NNEW(ithread,ITG,num_cpus);
	      for(i=0; i<num_cpus; i++)  {
		  ithread[i]=i;
		  pthread_create(&tid[i], NULL, (void *)calcviewmt, (void *)&ithread[i]);
	      }
	      for(i=0; i<num_cpus; i++)  pthread_join(tid[i], NULL);
	      
	      for(i=0;i<*ntr;i++){
		  for(j=1;j<num_cpus;j++){
		      adview[i]+=adview[i+j**ntr];
		  }
	      }
	      RENEW(adview,double,*ntr);
	      
	      for(i=0;i<2**nzsrad;i++){
		  for(j=1;j<num_cpus;j++){
		      auview[i]+=auview[i+j*2**nzsrad];
		  }
	      }
	      RENEW(auview,double,2**nzsrad);

/*	      for(i=0;i<*ntr;i++){
		  printf("radflowload adview = %" ITGFORMAT " %e\n",i,adview[i]);
	      }
	      for(i=0;i<2**nzsrad;i++){
		  printf("radflowload auview = %" ITGFORMAT " %e\n",i,auview[i]);
		  }*/

	      SFREE(dist);SFREE(idist);SFREE(e1);SFREE(e2);SFREE(e3);
              SFREE(pmid);SFREE(ithread);

	      /* postprocessing the viewfactors */

	      FORTRAN(postview,(ntr,sideload,nelemload,kontri,ntri,nloadtr,
				tenv,adview,auview,area,fenv,jqrad,irowrad,
                                nzsrad));

	      SFREE(area);

	      if(*iviewfile>=2){
		  
		  /* writing viewfactors to file */
		  
		  FORTRAN(writeview,(ntr,adview,auview,fenv,nzsrad,
				     jobnamef));
	      }
	      
	      if(*iviewfile==3){
		  
		  /* calculation of viewfactors only */
		  
		  FORTRAN(stop,());
	      }

	  }
      }

      /* assembling the radiation matrix */
      
      FORTRAN(radmatrix,(ntr,adrad,aurad,bcr,sideload,nelemload,
          xloadact,lakon,vold,ipkon,kon,co,nloadtr,tarea,tenv,physcon,
          erad,adview,auview,ithermal,iinc,iit,fenv,istep,dtime,ttime,
          time,iviewfile,xloadold,reltime,nmethod,mi,iemchange,nam,
          iamload,jqrad,irowrad,nzsrad));

      /* factoring the system of equations */

      /* the left hand side of the radiation matrix has probably
         changed if
         - the viewfactors were updated
         - a new step was started and NO CHANGE is not active
         - the emissivity coefficients were changed
         - a new increment was started in a stationary calculation
           (since the emissivity coefficients are ramped)
	   in that case the LU decomposition has to be repeated
           (i.e. call of dgesv) */

      if(((*ithermal==3)&&(*iviewfile>=0))||
         ((*iit==-1)&&(*iviewfile!=-2))||(*iemchange==1)||((*iit==0)&&(abs(*nmethod)==1))){

#if defined(PARDISO)
	if(ifactorization==1) pardiso_cleanup_as(ntr,&symmetryflag);
	pardiso_factor_as(adrad,aurad,adbrad,aubrad,&sigma,icolrad,
			  irowrad,ntr,nzsrad,jqrad);
	ifactorization=1;
#elif defined(SPOOLES)
	if(ifactorization==1) spooles_cleanup_rad();
	spooles_factor_rad(adrad,aurad,adbrad,aubrad,&sigma,
			   icolrad,irowrad,ntr,nzsrad,
			   &symmetryflag,&inputformat);
	ifactorization=1;
#else
	printf("*ERROR in radflowload: the SPOOLES library is not linked\n\n");
	FORTRAN(stop,());
#endif

      }

      /* solving the system of equations */

#if defined(PARDISO)
          pardiso_solve_as(bcr,ntr);

#elif defined(SPOOLES)
          spooles_solve_rad(bcr,ntr);
#endif
	
      if (info!=0){
	  printf("*ERROR IN RADFLOWLOAD: SINGULAR MATRIX*\n");}   
      
      else{ 
	  NNEW(q,double,*ntr);
	  FORTRAN(radresult,(ntr,xloadact,bcr,nloadtr,tarea,
				tenv,physcon,erad,auview,fenv,
			        irowrad,jqrad,nzsrad,q));
	  SFREE(q);
      }
      
  }

  SFREE(v);

  *adviewp=adview;*auviewp=auview;

  return;

} 

/* subroutine for multithreading of calcview */

void *calcviewmt(ITG *i){

    ITG indexad,indexau,indexdi,ntria,ntrib,nedelta;

    indexad=*i**ntr1;
    indexau=*i*2**nzsrad1;
    indexdi=*i**ntrit1;
    
    nedelta=(ITG)ceil(*ntri1/(double)num_cpus);
    ntria=*i*nedelta+1;
    ntrib=(*i+1)*nedelta;
    if(ntrib>*ntri1) ntrib=*ntri1;

//    printf("i=%" ITGFORMAT ",ntria=%" ITGFORMAT ",ntrib=%" ITGFORMAT "\n",i,ntria,ntrib);
//    printf("indexad=%" ITGFORMAT ",indexau=%" ITGFORMAT ",indexdi=%" ITGFORMAT "\n",indexad,indexau,indexdi);

    FORTRAN(calcview,(sideload1,vold1,co1,pmid1,e11,e21,e31,
                      kontri1,nloadtr1,&adview[indexad],
                      &auview[indexau],&dist[indexdi],&idist[indexdi],area1,
		      ntrit1,mi1,jqrad1,irowrad1,nzsrad1,&sidemean1,
                      &ntria,&ntrib));

    return NULL;
}
    
