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
             ITG *mortar,ITG *mpcinfo,double *tietol,ITG *ics,ITG *icontact,
	     ITG *nobject,char *objectset,ITG *istat){
  
  char description[13]="            ",*lakon=NULL,stiffmatrix[132]="";

  ITG *inum=NULL,k,*icol=NULL,*irow=NULL,ielas=0,icmd=0,iinc=1,nasym=0,
      mass[2]={0,0}, stiffness=1, buckling=0, rhsi=1, intscheme=0,*ncocon=NULL,
      *nshcon=NULL,mode=-1,noddiam=-1,*ipobody=NULL,inewton=0,coriolis=0,iout,
      ifreebody,*itg=NULL,ntg=0,ngraph=1,mt=mi[1]+1,ne0,*integerglob=NULL,      
      icfd=0,*inomat=NULL,*islavact=NULL,*islavnode=NULL,*nslavnode=NULL,
      *islavsurf=NULL,nmethodl,*kon=NULL,*ipkon=NULL,*ielmat=NULL;
      
  double *stn=NULL,*v=NULL,*een=NULL,cam[5],*xstiff=NULL,*stiini=NULL,*tper,
         *f=NULL,*fn=NULL,qa[3],*fext=NULL,*epn=NULL,*xstateini=NULL,
         *vini=NULL,*stx=NULL,*enern=NULL,*xbounact=NULL,*xforcact=NULL,
         *xloadact=NULL,*t1act=NULL,*ampli=NULL,*xstaten=NULL,*eei=NULL,
         *enerini=NULL,*cocon=NULL,*shcon=NULL,*physcon=NULL,*qfx=NULL,
         *qfn=NULL,*cgr=NULL,*xbodyact=NULL,*springarea=NULL,*emn=NULL,         
         *clearini=NULL,ptime,*emeini=NULL,*doubleglob=NULL,*au=NULL,
         *ad=NULL,*b=NULL,*aub=NULL,*adb=NULL,*pslavsurf=NULL,
         *pmastsurf=NULL,*cdn=NULL,*xstate=NULL,*fnext=NULL,*energyini=NULL,
         *energy=NULL,*ener=NULL;
            
   /* variables introduced for sensitivity calculation */
   
  ITG ndesi,numobject,*ndirdesi=NULL,*nodedesi=NULL;
   
  double distmin,*dfextminds=NULL,*df=NULL,*g0=NULL,*dgdx=NULL,*dgdv=NULL,
         *dgdxtot=NULL,*dgdxtotglob=NULL;
         

  FILE *f1;
  
#ifdef SGI
  ITG token;
#endif
  
  /* dummy arguments for the results call */

  double *veold=NULL,*accold=NULL,bet,gam,dtime,time,reltime=1.;

  icol=*icolp;irow=*irowp;

  kon=*konp;ipkon=*ipkonp;lakon=*lakonp;ielmat=*ielmatp;xstate=*xstatep;

  tper=&timepar[1];

  time=*tper;
  dtime=*tper;

  ne0=*ne;
  
  /* allocating a field for the sensitivity analysis operations */
   
  NNEW(ndirdesi,ITG,*nk);
  NNEW(nodedesi,ITG,*nk);
  NNEW(ad,double,neq[1]);
  NNEW(au,double,nzs[2]);
  NNEW(enerini,double,mi[0]**ne);
  NNEW(emeini,double,6*mi[0]**ne); 
   
  ener=*enerp;emeini=eme;
  
  for(k=0;k<mi[0]*ne0;++k){enerini[k]=ener[k];}
  
  /* reading the stiffness matrix from previous step for sensitivity analysis */
  /* matrix stored in <jobname>.stima file */
  
  strcpy(stiffmatrix,jobnamec);
  strcat(stiffmatrix,".stm");

  if((f1=fopen(stiffmatrix,"rb"))==NULL){
    printf("*ERROR in sensitivity: cannot open stiffness-matrix file for reading");
    exit(0);
  }
  
  if(fread(ad,sizeof(double),neq[1],f1)!=neq[1]){
	  printf("*ERROR in sensitivity reading the diagonal of the stiffness matrix in the .stm-file");
	  exit(0);
  }
      
  if(fread(au,sizeof(double),nzs[2],f1)!=nzs[2]){
	  printf("*ERROR in sensitivity reading the off-diagonals of the stiffness matrix in the .stm-file");
	  exit(0);
  }  
  
  fclose(f1);   
    
  /* determining the information of the designvariable set */
   
  FORTRAN(getdesiinfo,(set,istartset,iendset,ialset,nset,
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
  NNEW(stiini,double,6*mi[0]**ne);
  
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
          &ndesi,nodedesi,ndirdesi,sti);
	  
  iout=1;
  
  /* determining the system matrix and the external forces */

  NNEW(fext,double,*neq);
  NNEW(dfextminds,double,ndesi**neq);
    
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

  /* determining the values and the derivatives of the objective functions */

  NNEW(g0,double,*nobject);
  NNEW(dgdx,double,ndesi**nobject);
  NNEW(dgdxtot,double,ndesi**nobject);
  NNEW(dgdv,double,*neq**nobject);
  
  iout=-1; 
  objectivesmmain_se(co,nk,kon,ipkon,lakon,ne,v,stn,inum,stx,
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
          islavsurf,ielprop,prop,energyini,energy,&distmin,
          &ndesi,nodedesi,ndirdesi,nobject,objectset,g0,dgdx,dgdv,sti,
	  dgdxtot,dfextminds,df); 
  iout=1;
  
  /* calculating the normal direction for every designvariable */
  
  
  /* preparing the sensitivities for the output in the frd-file */
  
  NNEW(dgdxtotglob,double,4**nk**nobject);
  
  FORTRAN(sensitivity_glob,(dgdxtot,dgdxtotglob,nobject,&ndesi,nodedesi,
           ndirdesi,nk));
	       
  /* writing the sensitivities in the frd-file */
  
  for(numobject=0;numobject<*nobject;numobject++){
    frd_se(co,nk,stn,inum,nmethod,kode,filab,fn,&ptime,nstate_,istep,
        &iinc,&mode,&noddiam,description,mi,&ngraph,ne,cs,set,nset,
	istartset,iendset,ialset,thicke,jobnamec,output,
	&dgdxtotglob[numobject*4**nk],numobject);
  }  	
  
  
  SFREE(xbounact);SFREE(xforcact);SFREE(xloadact);SFREE(t1act);SFREE(ampli);
  SFREE(xbodyact);SFREE(v);SFREE(fn);SFREE(stx);SFREE(inum);SFREE(fext);
  SFREE(f);if(*nbody>0) SFREE(ipobody);SFREE(xstiff);SFREE(enerini);

  SFREE(dfextminds);SFREE(df);SFREE(g0);SFREE(dgdx);SFREE(ndirdesi);
  SFREE(nodedesi);SFREE(au);SFREE(ad);SFREE(stiini);SFREE(dgdv);
  SFREE(dgdxtotglob);
  
  *icolp=icol;*irowp=irow;

  *konp=kon;*ipkonp=ipkon;*lakonp=lakon;*ielmatp=ielmat;*enerp=ener;
  *xstatep=xstate;

  (*ttime)+=(*tper);
 
  return;
}
