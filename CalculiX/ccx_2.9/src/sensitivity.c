/*     CalculiX - A 3-dimensional finite element program                   */
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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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

void sensitivity(double *co, ITG *nk, ITG **konp, ITG **ipkonp, char **lakonp,
	     ITG *ne, 
	     ITG *nodeboun, ITG *ndirboun, double *xboun, ITG *nboun, 
	     ITG *ipompc, ITG *nodempc, double *coefmpc, char *labmpc,
             ITG *nmpc, 
	     ITG *nodeforc, ITG *ndirforc,double *xforc, ITG *nforc, 
	     ITG *nelemload, char *sideload, double *xload,
	     ITG *nload, ITG *nactdof, 
	     ITG **icolp, ITG *jq, ITG **irowp, ITG *neq, ITG *nzl, 
	     ITG *nmethod, ITG *ikmpc, ITG *ilmpc, ITG *ikboun, 
	     ITG *ilboun,
	     double *elcon, ITG *nelcon, double *rhcon, ITG *nrhcon,
	     double *alcon, ITG *nalcon, double *alzero, ITG **ielmatp,
	     ITG *ielorien, ITG *norien, double *orab, ITG *ntmat_,
	     double *t0, double *t1, double *t1old,
	     ITG *ithermal,double *prestr, ITG *iprestr, 
	     double *vold,ITG *iperturb, double *sti, ITG *nzs,  
	     ITG *kode, char *filab, double *eme,
             ITG *iexpl, double *plicon, ITG *nplicon, double *plkcon,
             ITG *nplkcon,
             double **xstatep, ITG *npmat_, char *matname, ITG *isolver,
             ITG *mi, ITG *ncmat_, ITG *nstate_, double *cs, ITG *mcs,
             ITG *nkon, double **enerp, double *xbounold,
	     double *xforcold, double *xloadold,
             char *amname, double *amta, ITG *namta,
	     ITG *nam, ITG *iamforc, ITG *iamload,
             ITG *iamt1, ITG *iamboun, double *ttime, char *output, 
             char *set, ITG *nset, ITG *istartset,
             ITG *iendset, ITG *ialset, ITG *nprint, char *prlab,
             char *prset, ITG *nener, double *trab, 
             ITG *inotr, ITG *ntrans, double *fmpc, char *cbody, ITG *ibody,
	     double *xbody, ITG *nbody, double *xbodyold, double *timepar,
	     double *thicke, char *jobnamec,char *tieset,ITG *ntie,
	     ITG *istep,ITG *nmat,ITG *ielprop,double *prop,char *typeboun,
	     ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,ITG *icontact){
  
  char description[13]="            ",*lakon=NULL;

  ITG *inum=NULL,k,*icol=NULL,*irow=NULL,ielas=0,icmd=0,iinc=1,nasym=0,i,j,ic,ir,
      mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,*ncocon=NULL,
      *nshcon=NULL,mode=-1,noddiam=-1,*ipobody=NULL,inewton=0,coriolis=0,iout,
      ifreebody,*itg=NULL,ntg=0,symmetryflag=0,inputformat=0,ngraph=1,im,
      mt=mi[1]+1,ne0,*integerglob=NULL,iglob=0,*ipneigh=NULL,*neigh=NULL,
      icfd=0,*inomat=NULL,*islavact=NULL,*islavnode=NULL,*nslavnode=NULL,
      *islavsurf=NULL,nretain,*iretain=NULL,*noderetain=NULL,*ndirretain=NULL,
      nmethodl,nintpoint,ifacecount,memmpc_,mpcfree,icascade,maxlenmpc,
      ncont=0,*itietri=NULL,*koncont=NULL,nslavs=0,ismallsliding=0,
      *itiefac=NULL,*imastnode=NULL,*nmastnode=NULL,*imastop=NULL,
      *iponoels=NULL,*inoels=NULL,*ipe=NULL,*ime=NULL,iit=-1,iflagact=0,
      icutb=0,*kon=NULL,*ipkon=NULL,*ielmat=NULL,ialeatoric=0;

  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*xstiff=NULL,*stiini=NULL,*tper,
         *f=NULL,*fn=NULL,qa[3],*fext=NULL,*epn=NULL,*xstateini=NULL,
         *vini=NULL,*stx=NULL,*enern=NULL,*xbounact=NULL,*xforcact=NULL,
         *xloadact=NULL,*t1act=NULL,*ampli=NULL,*xstaten=NULL,*eei=NULL,
         *enerini=NULL,*cocon=NULL,*shcon=NULL,*physcon=NULL,*qfx=NULL,
         *qfn=NULL,sigma=0.,*cgr=NULL,*xbodyact=NULL,*vr=NULL,*vi=NULL,
         *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*springarea=NULL,
         *eenmax=NULL,*fnr=NULL,*fni=NULL,*emn=NULL,*clearini=NULL,ptime,
         *emeini=NULL,*doubleglob=NULL,*au=NULL,*ad=NULL,*b=NULL,*aub=NULL,
         *adb=NULL,*pslavsurf=NULL,*pmastsurf=NULL,*cdn=NULL,*cdnr=NULL,
         *cdni=NULL,*submatrix=NULL,*xnoels=NULL,*cg=NULL,*straight=NULL,
         *areaslav=NULL,*xmastnor=NULL,theta=0.,*ener=NULL,*xstate=NULL,
         *fnext=NULL,*energyini=NULL,*energy=NULL;
   
   /* variables introduced for sensitivity calculation */
   
  ITG ndesi,*ndirdesi=NULL,*nodedesi=NULL;
   
  double distmin,*dfextminds=NULL,*df=NULL;
   
#ifdef SGI
  ITG token;
#endif
  
  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time,reltime=1.;

  icol=*icolp;irow=*irowp;

  kon=*konp;ipkon=*ipkonp;lakon=*lakonp;ielmat=*ielmatp;ener=*enerp;
  xstate=*xstatep;

  tper=&timepar[1];

  time=*tper;
  dtime=*tper;

  ne0=*ne;
  
  /* allocating a field for the definition of the designvariables */
  
  NNEW(ndirdesi,ITG,*nk);
  NNEW(nodedesi,ITG,*nk);
  
   /* determining the information of the designvariable set */
   
  FORTRAN(getdesiinfo,(set,istartset,iendset,ialset,nset,istep,
     	    mi,nactdof,&ndesi,ndirdesi,nodedesi,ntie,tieset));  
     	    
  /* calculation of the smallest distance between nodes */
  
  FORTRAN(smalldist,(co,&distmin,lakon,ipkon,kon,ne));

  /* allocating fields for the actual external loading */

  NNEW(xbounact,double,*nboun);
  for(k=0;k<*nboun;++k){xbounact[k]=xbounold[k];}
  NNEW(xforcact,double,*nforc);
  NNEW(xloadact,double,2**nload);
  NNEW(xbodyact,double,7**nbody);
  /* copying the rotation axis and/or acceleration vector */
  for(k=0;k<7**nbody;k++){xbodyact[k]=xbody[k];}
  if(*ithermal==1){
    NNEW(t1act,double,*nk);
    for(k=0;k<*nk;++k){t1act[k]=t1old[k];}
  }
  
  /* assigning the body forces to the elements */ 

  if(*nbody>0){
      ifreebody=*ne+1;
      NNEW(ipobody,ITG,2*ifreebody**nbody);
      for(k=1;k<=*nbody;k++){
	  FORTRAN(bodyforce,(cbody,ibody,ipobody,nbody,set,istartset,
			     iendset,ialset,&inewton,nset,&ifreebody,&k));
	  RENEW(ipobody,ITG,2*(*ne+ifreebody));
      }
      RENEW(ipobody,ITG,2*(ifreebody-1));
  }

  /* allocating a field for the instantaneous amplitude */

  NNEW(ampli,double,*nam);

  FORTRAN(tempload,(xforcold,xforc,xforcact,iamforc,nforc,xloadold,xload,
	      xloadact,iamload,nload,ibody,xbody,nbody,xbodyold,xbodyact,
	      t1old,t1,t1act,iamt1,nk,amta,
	      namta,nam,ampli,&time,&reltime,ttime,&dtime,ithermal,nmethod,
              xbounold,xboun,xbounact,iamboun,nboun,
	      nodeboun,ndirboun,nodeforc,ndirforc,istep,&iinc,
	      co,vold,itg,&ntg,amname,ikboun,ilboun,nelemload,sideload,mi,
              ntrans,trab,inotr,veold,integerglob,doubleglob,tieset,istartset,
              iendset,ialset,ntie,nmpc,ipompc,ikmpc,ilmpc,nodempc,coefmpc));

  /* determining the internal forces and the stiffness coefficients */

  NNEW(f,double,*neq);

  /* allocating a field for the stiffness matrix */

  NNEW(xstiff,double,(long long)27*mi[0]**ne);

  iout=-1;
  NNEW(v,double,mt**nk);
  NNEW(fn,double,mt**nk);
  NNEW(df,double,ndesi**neq);
  NNEW(stx,double,6*mi[0]**ne);
  NNEW(inum,ITG,*nk);
  results_se(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	  elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	  ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
	  prestr,iprestr,filab,eme,emn,een,iperturb,
	  f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	  ndirboun,xbounact,nboun,ipompc,
	  nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
	  &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	  xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	  &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	  emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	  iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	  fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	  &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
	  sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	  mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	  islavsurf,ielprop,prop,energyini,energy,df,&distmin,
	  &ndesi,nodedesi,ndirdesi);
  SFREE(v);SFREE(fn);SFREE(stx);SFREE(inum);
  iout=1;
  
  /* determining the system matrix and the external forces */

  NNEW(ad,double,*neq);
  NNEW(fext,double,*neq);
  NNEW(dfextminds,double,ndesi**neq);

  NNEW(au,double,*nzs);
  nmethodl=*nmethod;
  

  mafillsmmain_se(co,nk,kon,ipkon,lakon,ne,nodeboun,ndirboun,xbounact,nboun,
          ipompc,nodempc,coefmpc,nmpc,nodeforc,ndirforc,xforcact,
          nforc,nelemload,sideload,xloadact,nload,xbodyact,ipobody,
          nbody,cgr,ad,au,nactdof,icol,jq,irow,neq,nzl,&nmethodl,
          ikmpc,ilmpc,ikboun,ilboun,
          elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
          ielorien,norien,orab,ntmat_,
          t0,t1act,ithermal,prestr,iprestr,vold,iperturb,sti,
          nzs,stx,adb,aub,iexpl,plicon,nplicon,plkcon,nplkcon,
          xstiff,npmat_,&dtime,matname,mi,
	  ncmat_,mass,&stiffness,&buckling,&rhsi,&intscheme,physcon,
	  shcon,nshcon,cocon,ncocon,ttime,&time,istep,&iinc,&coriolis,
          ibody,xloadold,&reltime,veold,springarea,nstate_,
	  xstateini,xstate,thicke,integerglob,doubleglob,
          tieset,istartset,iendset,ialset,ntie,&nasym,pslavsurf,
          pmastsurf,mortar,clearini,ielprop,prop,&ne0,fnext,
	  &distmin,&ndesi,nodedesi,ndirdesi,dfextminds);


  /* determining the right hand side */

  NNEW(b,double,*neq);
  for(k=0;k<*neq;++k){
      b[k]=fext[k]-f[k];
  }
  
    /*  for(k=0;k<neq[1];++k){printf("f2=%" ITGFORMAT ",%e\n",k,f[k]);}
      for(k=0;k<neq[1];++k){printf("fext2=%" ITGFORMAT ",%e\n",k,fext[k]);}
      for(k=0;k<neq[1];++k){printf("ad2=%" ITGFORMAT ",%e\n",k,ad[k]);}
      for(k=0;k<nzs[1];++k){printf("au2=%" ITGFORMAT ",%e\n",k,au[k]);}
      for(k=0;k<neq[1];++k){printf("b2=%" ITGFORMAT ",%e\n",k,b[k]);}  */

  SFREE(fext);SFREE(f);

  /* generation of a substructure stiffness matrix */

  if(*nmethod!=0){

    if(*isolver==0){
#ifdef SPOOLES
      spooles(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,&symmetryflag,
              &inputformat,&nzs[2]);
#else
            printf("*ERROR in linstatic: the SPOOLES library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }
    else if((*isolver==2)||(*isolver==3)){
      if(nasym>0){
	  printf(" *ERROR in nonlingeo: the iterative solver cannot be used for asymmetric matrices\n\n");
	  FORTRAN(stop,());
      }
      preiter(ad,&au,b,&icol,&irow,neq,nzs,isolver,iperturb);
    }
    else if(*isolver==4){
#ifdef SGI
      if(nasym>0){
	  printf(" *ERROR in nonlingeo: the SGI solver cannot be used for asymmetric matrices\n\n");
	  FORTRAN(stop,());
      }
      token=1;
      sgi_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,token);
#else
            printf("*ERROR in linstatic: the SGI library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }
    else if(*isolver==5){
#ifdef TAUCS
      if(nasym>0){
	  printf(" *ERROR in nonlingeo: the TAUCS solver cannot be used for asymmetric matrices\n\n");
	  FORTRAN(stop,());
      }
      tau(ad,&au,adb,aub,&sigma,b,icol,&irow,neq,nzs);
#else
            printf("*ERROR in linstatic: the TAUCS library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }
    else if(*isolver==7){
#ifdef PARDISO
      pardiso_main(ad,au,adb,aub,&sigma,b,icol,irow,neq,nzs,
		   &symmetryflag,&inputformat,jq,&nzs[2]);
#else
            printf("*ERROR in linstatic: the PARDISO library is not linked\n\n");
            FORTRAN(stop,());
#endif
    }

    SFREE(ad);SFREE(au);

    /* calculating the displacements and the stresses and storing */
    /* the results in frd format for each valid eigenmode */

    NNEW(v,double,mt**nk);
    NNEW(fn,double,mt**nk);
    NNEW(df,double,ndesi**neq);    
    NNEW(stn,double,6**nk);
    NNEW(inum,ITG,*nk);
    NNEW(stx,double,6*mi[0]**ne);
  
    if(strcmp1(&filab[261],"E   ")==0) NNEW(een,double,6**nk);
    if(strcmp1(&filab[2697],"ME  ")==0) NNEW(emn,double,6**nk);
    if(strcmp1(&filab[522],"ENER")==0) NNEW(enern,double,*nk);
    if(strcmp1(&filab[2175],"CONT")==0) NNEW(cdn,double,6**nk);

    NNEW(eei,double,6*mi[0]**ne);
    if(*nener==1){
	NNEW(stiini,double,6*mi[0]**ne);
	NNEW(enerini,double,mi[0]**ne);}

    results_se(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
	    elcon,nelcon,rhcon,nrhcon,alcon,nalcon,alzero,ielmat,
	    ielorien,norien,orab,ntmat_,t0,t1act,ithermal,
            prestr,iprestr,filab,eme,emn,een,iperturb,
	    f,fn,nactdof,&iout,qa,vold,b,nodeboun,
	    ndirboun,xbounact,nboun,ipompc,
	    nodempc,coefmpc,labmpc,nmpc,nmethod,cam,neq,veold,accold,
	    &bet,&gam,&dtime,&time,ttime,plicon,nplicon,plkcon,nplkcon,
	    xstateini,xstiff,xstate,npmat_,epn,matname,mi,&ielas,
	    &icmd,ncmat_,nstate_,stiini,vini,ikboun,ilboun,ener,enern,
	    emeini,xstaten,eei,enerini,cocon,ncocon,set,nset,istartset,
	    iendset,ialset,nprint,prlab,prset,qfx,qfn,trab,inotr,ntrans,
	    fmpc,nelemload,nload,ikmpc,ilmpc,istep,&iinc,springarea,
	    &reltime,&ne0,xforc,nforc,thicke,shcon,nshcon,
	    sideload,xloadact,xloadold,&icfd,inomat,pslavsurf,pmastsurf,
	    mortar,islavact,cdn,islavnode,nslavnode,ntie,clearini,
	    islavsurf,ielprop,prop,energyini,energy,df,&distmin,
	    &ndesi,nodedesi,ndirdesi);

    SFREE(eei);
    if(*nener==1){
	SFREE(stiini);SFREE(enerini);}

    memcpy(&vold[0],&v[0],sizeof(double)*mt**nk);
    memcpy(&sti[0],&stx[0],sizeof(double)*6*mi[0]*ne0);

    ++*kode;

    /* for cyclic symmetric sectors: duplicating the results */

    if(*mcs>0){
	ptime=*ttime+time;
      frdcyc(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,kode,filab,een,t1act,
		   fn,&ptime,epn,ielmat,matname,cs,mcs,nkon,enern,xstaten,
                   nstate_,istep,&iinc,iperturb,ener,mi,output,ithermal,
                   qfn,ialset,istartset,iendset,trab,inotr,ntrans,orab,
	           ielorien,norien,sti,veold,&noddiam,set,nset,emn,thicke,
	           jobnamec,&ne0,cdn,mortar,nmat);
    }
    else{
	if(strcmp1(&filab[1044],"ZZS")==0){
	    NNEW(neigh,ITG,40**ne);
	    NNEW(ipneigh,ITG,*nk);
	}
	ptime=*ttime+time;
	frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1act,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,stx,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat);
	if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
    }

    SFREE(v);SFREE(stn);SFREE(inum);
    SFREE(b);SFREE(stx);SFREE(fn);

    if(strcmp1(&filab[261],"E   ")==0) SFREE(een);
    if(strcmp1(&filab[2697],"ME  ")==0) SFREE(emn);
    if(strcmp1(&filab[522],"ENER")==0) SFREE(enern);
    if(strcmp1(&filab[2175],"CONT")==0) SFREE(cdn);

  }
  else {

    /* error occurred in mafill: storing the geometry in frd format */

    ++*kode;
    NNEW(inum,ITG,*nk);for(k=0;k<*nk;k++) inum[k]=1;
    if(strcmp1(&filab[1044],"ZZS")==0){
	NNEW(neigh,ITG,40**ne);
	NNEW(ipneigh,ITG,*nk);
    }
    ptime=*ttime+time;
    frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
	    kode,filab,een,t1,fn,&ptime,epn,ielmat,matname,enern,xstaten,
	    nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
	    ntrans,orab,ielorien,norien,description,ipneigh,neigh,
	    mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
	    cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
	    thicke,jobnamec,output,qfx,cdn,mortar,cdnr,cdni,nmat);
    if(strcmp1(&filab[1044],"ZZS")==0){SFREE(ipneigh);SFREE(neigh);}
    SFREE(inum);FORTRAN(stop,());

  }

  /* updating the loading at the end of the step; 
     important in case the amplitude at the end of the step
     is not equal to one */

  for(k=0;k<*nboun;++k){xbounold[k]=xbounact[k];}
  for(k=0;k<*nforc;++k){xforcold[k]=xforcact[k];}
  for(k=0;k<2**nload;++k){xloadold[k]=xloadact[k];}
  for(k=0;k<7**nbody;k=k+7){xbodyold[k]=xbodyact[k];}
  if(*ithermal==1){
    for(k=0;k<*nk;++k){t1old[k]=t1act[k];}
    for(k=0;k<*nk;++k){vold[mt*k]=t1act[k];}
  }

  SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(t1act);SFREE(ampli);
  SFREE(xbodyact);if(*nbody>0) SFREE(ipobody);SFREE(xstiff);

  if(iglob==1){SFREE(integerglob);SFREE(doubleglob);}

  *icolp=icol;*irowp=irow;

  *konp=kon;*ipkonp=ipkon;*lakonp=lakon;*ielmatp=ielmat;*enerp=ener;
  *xstatep=xstate;

  (*ttime)+=(*tper);
 
  return;
}
