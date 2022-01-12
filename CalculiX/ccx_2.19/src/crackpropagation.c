/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2021 Guido Dhondt                          */

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
#include <string.h>
#include <stdlib.h>
#include "CalculiX.h"

#define max(a,b) ((a) >= (b) ? (a) : (b))

void crackpropagation(ITG **ipkonp,ITG **konp,char **lakonp,ITG *ne,ITG *nk,
		      char *jobnamec,ITG *nboun,ITG *iamboun,double *xboun,
		      ITG *nload,char *sideload,ITG *iamload,ITG *nforc,
		      ITG *iamforc,double *xforc,ITG *ithermal,double *t1,
		      ITG *iamt1,double **cop,ITG *nkon,ITG *mi,ITG **ielmatp,
		      char *matname,char *output,ITG *nmat,char *set,ITG *nset,
		      ITG *istartset,ITG *iendset,ITG *ialset,ITG *jmax,
		      double *timepar,ITG *nelcon,double *elcon,ITG *ncmat_,
		      ITG *ntmat_,ITG *istep,char *filab,ITG *nmethod,
		      ITG *mei){

  char *lakon=NULL,description[13]="            ",
    *param=NULL;
  
  ITG *kontri=NULL,ntri,*iedg=NULL,*ipoed=NULL,*ieled=NULL,nedg,
    *ibounedg=NULL,nbounedg,*ibounnod=NULL,nbounnod,*iedno=NULL,
    *integerglob=NULL,iglob,*ifront=NULL,nfront,*ifrontrel=NULL,
    *ifront2=NULL,*ifrontrel2=NULL,*istartfront=NULL,*iendfront=NULL,
    nnfront,*istartcrackfro=NULL,*iendcrackfro=NULL,ncrack,i,j,
    *ibounnod2=NULL,*istartcrackbou=NULL,*iendcrackbou=NULL,imat,
    *isubsurffront=NULL,*iedno2=NULL,nkold,*ipkon=NULL,*kon=NULL,
    *iresort=NULL,*inum=NULL,kode=1,nstate_=0,iinc,ier=0,nstepf,
    mode=-1,noddiam=-1,*inotr=NULL,ntrans,*ielorien=NULL,norien,*ipneigh=NULL,
    *neigh=NULL,ngraph=1,mortar=0,*ielprop=NULL,*ielmat=NULL,icritic=0,
    *ifrontprop=NULL,*ifronteq=NULL,*istartfronteq=NULL,*iendfronteq=NULL,
    nfronteq,ncyc,*idist=NULL,ncrconst,nstep,nproc,ncrtem,law,nstepf2,
    *iincglob=NULL,nparam,ncyctot=0,ieqspace,*integerglobf=NULL,lcf;

  double *doubleglob=NULL,*stress=NULL,*xt=NULL,*xn=NULL,*xa=NULL,
    *acrack=NULL,*xk1=NULL,*xk2=NULL,*xk3=NULL,*doubleglobf=NULL,
    *xkeq=NULL,*phi=NULL,*psi=NULL,*co=NULL,*da=NULL,*stress2=NULL,*v=NULL,
    *stn=NULL,*een=NULL,*fn=NULL,time=0.,*epn=NULL,*enern=NULL,*sti=NULL,
    *xstaten=NULL,*qfn=NULL,*trab=NULL,*orab=NULL,*vr=NULL,*vi=NULL,
    *stnr=NULL,*stni=NULL,*vmax=NULL,*stnmax=NULL,*veold=NULL,*ener=NULL,
    *cs=NULL,*eenmax=NULL,*fnr=NULL,*fni=NULL,*emn=NULL,*thicke=NULL,*qfx=NULL,
    *cdn=NULL,*cdnr=NULL,*cdni=NULL,*prop=NULL,*costruc=NULL,*costruc2=NULL,
    *dadn=NULL,*wk1glob=NULL,*wk2glob=NULL,*wk3glob=NULL,*dkeqglob=NULL,
    *phiglob=NULL,*dadnglob=NULL,*dnglob=NULL,*charlen,damax,*crcon=NULL,
    *angle=NULL,*posfront=NULL,*dist=NULL,*tempf=NULL,*stressf=NULL,
    *a=NULL,*amin=NULL,*shape=NULL,*dk=NULL,*p=NULL,*temp=NULL,*temp2=NULL,
    *crconloc=NULL,datarget,phimax,*acrackglob=NULL,scaling,*hcfstress=NULL,
    *wk1=NULL,*wk2=NULL,*wk3=NULL,*xkeqmin=NULL,*xkeqmax=NULL,*domstep=NULL,
    *dkeq=NULL,*domphi=NULL,*xkeqminglob=NULL,*xkeqmaxglob=NULL,
    *domstepglob=NULL,*hcftemp=NULL,*xk1f=NULL,*xk2f=NULL,*xk3f=NULL,
    *xkeqf=NULL,*phif=NULL,*psif=NULL,*rglob=NULL,*r=NULL;

  co=*cop;lakon=*lakonp;ipkon=*ipkonp;kon=*konp;ielmat=*ielmatp;

  damax=timepar[0];
  phimax=timepar[1]*atan(1.)/45.;
  imat=(ITG)(floor(timepar[2]));
  nproc=(ITG)(floor(timepar[3]));
  scaling=timepar[4];

  /* storing the crack propagation data in crcon */
  
  crackpropdata(jobnamec,nelcon,elcon,&crcon,&ncrconst,&ncrtem,&imat,
		matname,ntmat_,ncmat_,&param,&nparam,&law);
      
  /* reading the LCF master file
     iglob>0: static results */

  nstep=0;
  iglob=1;
  getuncrackedresults(&jobnamec[396],&integerglob,&doubleglob,
		      &iglob,&nstep);
      
  /* reading the HCF master file
     iglob<0: frequency results */

  if(mei[0]>0){
    nstepf=0;
    iglob=-1;
    getuncrackedresults(&jobnamec[528],&integerglobf,&doubleglobf,
			&iglob,&nstepf);
  }
    
  NNEW(wk1glob,double,3**nk);
  NNEW(wk2glob,double,3**nk);
  NNEW(wk3glob,double,3**nk);
  NNEW(xkeqminglob,double,3**nk);
  NNEW(xkeqmaxglob,double,3**nk);
  NNEW(dkeqglob,double,3**nk);
  NNEW(rglob,double,3**nk);
  NNEW(phiglob,double,3**nk);
  NNEW(dadnglob,double,3**nk);
  NNEW(dnglob,double,3**nk);
  NNEW(acrackglob,double,3**nk);
  NNEW(iincglob,ITG,3**nk);
  NNEW(domstepglob,double,3**nk);

  /* loop over crack propagation increments */
  
  for(iinc=0;iinc<jmax[0];iinc++){

    printf("Increment %d\n\n",iinc+1);
  
    /* catalogue all triangles belonging to the crack(s)
       they are supposed to be of type S3 */

    NNEW(kontri,ITG,3**ne);

    FORTRAN(cattri,(ne,lakon,ipkon,kon,kontri,&ntri));

    RENEW(kontri,ITG,3*ntri);

    /* catalogue all edges of the crack(s) */

    NNEW(ipoed,ITG,*nk);
    NNEW(iedg,ITG,9*ntri+3);
    NNEW(ieled,ITG,6*ntri);

    FORTRAN(catedges_crackprop,(ipoed,iedg,&ntri,ieled,kontri,&nedg,&ier));
    if(ier==1) break;

    RENEW(iedg,ITG,3*nedg);
    RENEW(ieled,ITG,2*nedg);

    /* check for the boundary edges and boundary nodes of the
       shell mesh (= crack shape(s)) */

    NNEW(ibounedg,ITG,nedg);
    NNEW(ibounnod,ITG,nedg);
    NNEW(iedno,ITG,2*nedg);

    FORTRAN(extern_crackprop,(ieled,&nedg,ibounedg,&nbounedg,ibounnod,
			      &nbounnod,iedg,iedno));

    RENEW(ibounedg,ITG,nbounedg);
    RENEW(ibounnod,ITG,nbounnod);
    RENEW(iedno,ITG,2*nbounnod);

    /* interpolating the stress in all boundary nodes
       determining the boundary crack nodes inside the structure */

    NNEW(temp,double,nstep*nbounnod);
    NNEW(stress,double,6*nstep*nbounnod);
    NNEW(ifront,ITG,nbounnod);
    NNEW(ifrontrel,ITG,nbounnod);
    NNEW(costruc,double,3*nbounnod);

    FORTRAN(interpolextnodes,(ibounnod,&nbounnod,co,doubleglob,integerglob,
			      stress,ifront,&nfront,ifrontrel,costruc,
			      temp,&nstep));

    RENEW(ifront,ITG,nfront);
    RENEW(ifrontrel,ITG,nfront);

    /* interpolating the hcf stress in all boundary nodes */

    if(mei[0]>0){
      
      NNEW(hcftemp,double,nbounnod);
      NNEW(hcfstress,double,6*nbounnod);
      
      FORTRAN(interpolextnodesf,(ibounnod,&nbounnod,co,doubleglobf,integerglobf,
				 hcfstress,hcftemp,&nstepf,&mei[0]));

      /* scaling the hcf stress */
      
      for(i=0;i<6*nbounnod;i++){hcfstress[i]*=scaling;}
      
      SFREE(hcftemp);
    }

    /* determining the internal boundary nodes
       sorting them according to adjacency
       add the nodes immediately outside the structure */

    NNEW(ifront2,ITG,3*nfront);
    NNEW(ifrontrel2,ITG,3*nfront);
    NNEW(istartfront,ITG,nfront);
    NNEW(iendfront,ITG,nfront);
    RENEW(ifront,ITG,3*nfront);
    RENEW(ifrontrel,ITG,3*nfront);
    NNEW(istartcrackfro,ITG,nfront);
    NNEW(iendcrackfro,ITG,nfront);
    NNEW(ibounnod2,ITG,nbounnod);
    NNEW(iresort,ITG,nbounnod);
    NNEW(istartcrackbou,ITG,nfront);
    NNEW(iendcrackbou,ITG,nfront);
    NNEW(isubsurffront,ITG,nfront);
    NNEW(iedno2,ITG,2*nbounnod);
    NNEW(temp2,double,nstep*6*nbounnod);
    NNEW(stress2,double,nstep*6*nbounnod);
    NNEW(costruc2,double,3*nbounnod);

    /* determine the order of the nodes along the crack front(s) */
  
    FORTRAN(adjacentbounodes,(ifront,ifrontrel,&nfront,iedno,iedg,&nnfront,
			      ifront2,ifrontrel2,ibounnod,&nbounnod,istartfront,
			      iendfront,ibounedg,istartcrackfro,iendcrackfro,
			      &ncrack,ibounnod2,istartcrackbou,iendcrackbou,
			      isubsurffront,iedno2,stress,stress2,iresort,
			      ieled,kontri,costruc,costruc2,temp,temp2,&nstep,
			      &ier));
    if(ier==1) break;

    SFREE(iresort);
    RENEW(ifront,ITG,nfront);
    RENEW(ifrontrel,ITG,nfront);
    RENEW(istartfront,ITG,nnfront);
    RENEW(iendfront,ITG,nnfront);
    if(iinc==0) NNEW(charlen,double,ncrack);
    RENEW(istartcrackfro,ITG,ncrack);
    RENEW(iendcrackfro,ITG,ncrack);
    RENEW(istartcrackbou,ITG,ncrack);
    RENEW(iendcrackbou,ITG,ncrack);
    RENEW(isubsurffront,ITG,nfront);
    SFREE(ifront2);SFREE(ifrontrel2);SFREE(ibounnod2);SFREE(iedno2);
    SFREE(stress2);SFREE(temp2);SFREE(costruc2);

    NNEW(xt,double,3*nfront);
    NNEW(xn,double,3*nfront);
    NNEW(xa,double,3*nfront);

    NNEW(angle,double,nnfront);

    /* determine a local coordinate system in each front node based on
       the local stress tensor */
  
    FORTRAN(createlocalsys,(&nnfront,istartfront,iendfront,ifront,co,xt,xn,xa,
			    &nfront,ifrontrel,stress,iedno,ibounedg,ieled,
			    kontri,isubsurffront,istartcrackfro,iendcrackfro,
			    &ncrack,angle,&nstep,&ier));
    if(ier==1) break;

    /* determine a crack length for each front node*/

    NNEW(acrack,double,nfront);
    NNEW(posfront,double,nfront);
    FORTRAN(cracklength,(&ncrack,istartcrackfro,iendcrackfro,co,istartcrackbou,
			 iendcrackbou,costruc,ibounnod,xt,acrack,istartfront,
			 iendfront,&nnfront,isubsurffront,ifrontrel,
			 ifront,posfront,doubleglob,integerglob,
			 &nproc,&iinc,acrackglob,&ier));
    if(ier==1) break;

    NNEW(dist,double,nfront);
    NNEW(idist,ITG,nfront);
    NNEW(a,double,nfront);
    //    NNEW(amin,double,ncrack);

    FORTRAN(cracklength_smoothing,(&nnfront,isubsurffront,istartfront,
				   iendfront,ifrontrel,costruc,dist,a,
				   istartcrackfro,iendcrackfro,acrack,amin,
				   &nfront,&ncrack,idist,&nstep));

    SFREE(dist);SFREE(idist);SFREE(a);

    /* calculating the shape factors F=K/(sigma*sqrt(pi*a)) */

    NNEW(shape,double,3*nfront);

    FORTRAN(crackshape,(&nnfront,ifront,istartfront,iendfront,isubsurffront,
			angle,posfront,shape));
    
    SFREE(angle);SFREE(posfront);
    
    /* calculating the stress intensity factors (LCF) */

    NNEW(xk1,double,nstep*nfront);
    NNEW(xk2,double,nstep*nfront);
    NNEW(xk3,double,nstep*nfront);
    NNEW(xkeq,double,nstep*nfront);
    NNEW(phi,double,nstep*nfront);
    NNEW(psi,double,nstep*nfront);

    FORTRAN(stressintensity,(&nfront,ifrontrel,stress,xt,xn,xa,xk1,xk2,
			     xk3,xkeq,phi,psi,acrack,shape,&nstep));

    SFREE(psi);
    
    /* calculating the stress intensity factors (HCF) */

    if(mei[0]>0){
      
      if(mei[1]>nstep){
	printf("*ERROR in crackpropagation.c: HCF mission step %d\n",mei[1]);
	printf("       exceeds the number of LCF steps %d\n\n",nstep);
	FORTRAN(stop,());
      }

      /* creating LCF+HCF and LCF-HCF */
      
      NNEW(tempf,double,2*nbounnod);
      NNEW(stressf,double,6*2*nbounnod);

      FORTRAN(combilcfhcf,(tempf,stressf,stress,hcfstress,temp,&nbounnod,
			   mei,&nstep));
      
      SFREE(hcfstress);
      
      NNEW(xk1f,double,2*nfront);
      NNEW(xk2f,double,2*nfront);
      NNEW(xk3f,double,2*nfront);
      NNEW(xkeqf,double,2*nfront);
      NNEW(phif,double,2*nfront);
      NNEW(psif,double,2*nfront);

      nstepf2=2;
      FORTRAN(stressintensity,(&nfront,ifrontrel,stressf,xt,xn,xa,xk1f,xk2f,
			       xk3f,xkeqf,phif,psif,acrack,shape,&nstepf2));

      SFREE(psif);
    }

    SFREE(shape);
    
    /* smoothing the equivalent k-factors and deflection angles */
    
    NNEW(dist,double,nfront);
    NNEW(idist,ITG,nfront);
    NNEW(dk,double,nstep*nfront);
    NNEW(p,double,nstep*nfront);

    FORTRAN(stressintensity_smoothing,(&nnfront,isubsurffront,istartfront,
				       iendfront,ifrontrel,costruc,dist,
				       istartcrackfro,iendcrackfro,xkeq,phi,
				       &nfront,&ncrack,dk,p,idist,&phimax,
				       &nstep));

    SFREE(dist);SFREE(idist);SFREE(dk);SFREE(p);
    
    /* smoothing the HCF equivalent k-factors and deflection angles */

    if(mei[0]>0){
    
      NNEW(dist,double,nfront);
      NNEW(idist,ITG,nfront);
      NNEW(dk,double,2*nfront);
      NNEW(p,double,2*nfront);

      FORTRAN(stressintensity_smoothing,(&nnfront,isubsurffront,istartfront,
					 iendfront,ifrontrel,costruc,dist,
					 istartcrackfro,iendcrackfro,xkeqf,phif,
					 &nfront,&ncrack,dk,p,idist,&phimax,
					 &nstepf2));

      SFREE(dist);SFREE(idist);SFREE(dk);SFREE(p);
      
    }

    /* determine the target crack propagation increment */
    
    FORTRAN(calcdatarget,(ifront,co,&nnfront,istartfront,iendfront,
			  isubsurffront,&damax,&datarget,acrack,&nstep));
    
    /* propagate the nodes at the crack front(s) 
       (generate new nodes) */

    lcf=1;
    
    if(mei[0]>0){

      /* check whether propagation due to hcf */
    
      NNEW(dadn,double,nfront);
      NNEW(wk1,double,nfront);
      NNEW(wk2,double,nfront);
      NNEW(wk3,double,nfront);
      NNEW(xkeqmin,double,nfront);
      NNEW(xkeqmax,double,nfront);
      NNEW(dkeq,double,nfront);
      NNEW(r,double,nfront);
      NNEW(domstep,double,nfront);
      NNEW(domphi,double,nfront);
      NNEW(crconloc,double,ncrconst);
    
      FORTRAN(crackrate,(&nfront,ifrontrel,xkeqf,phif,ifront,dadn,&ncyc,
			 &icritic,&datarget,crcon,tempf,
			 &ncrtem,crconloc,&ncrconst,xk1f,xk2f,xk3f,&nstepf2,
			 acrack,wk1,wk2,wk3,xkeqmin,xkeqmax,
			 dkeq,domstep,domphi,param,&nparam,&law,&ier,r));
      
      if((icritic>0)||(ier==1)){
	if(icritic>0){
	  printf("WARNING: Exit because a critical value for Keq is reached along the crack front due to the modal loading \n");
	}else{
	  printf("WARNING: Exit because an error occurred in crackrate.f \n");
	}
	printf("Number of iteration= %d\n\n" ,(iinc+1));
	
	SFREE(xkeq);SFREE(phi);SFREE(xk1);SFREE(xk2);SFREE(xk3);
	SFREE(xkeqf);SFREE(phif);SFREE(xk1f);SFREE(xk2f);SFREE(xk3f);
	SFREE(crconloc);
	SFREE(kontri);SFREE(ipoed);SFREE(iedg);SFREE(ieled);SFREE(ibounedg);
	SFREE(ibounnod);SFREE(iedno);SFREE(stress);SFREE(ifront);
	SFREE(ifrontrel);SFREE(istartfront);SFREE(iendfront);SFREE(xt);
	SFREE(stressf);SFREE(tempf);
	SFREE(istartcrackfro);SFREE(iendcrackfro);SFREE(temp);
	SFREE(istartcrackbou);SFREE(iendcrackbou);SFREE(r);
	SFREE(isubsurffront);SFREE(acrack);SFREE(dkeq);SFREE(costruc);
	SFREE(dadn);SFREE(wk1);SFREE(wk2);SFREE(wk3);SFREE(xkeqmin);
	SFREE(xkeqmax);SFREE(domstep);SFREE(domphi);
	
	break;
      }

      if(icritic<0){

	/* no crack propagation due to hcf: replace step mei[1]
           of xkeq by the first step of xkeqf (=lcf+hcf) */

	for(i=0;i<nfront;i++){
	  xkeq[nstep*i+mei[1]-1]=xkeqf[2*i];
	}
	
	SFREE(xkeqf);SFREE(phif);SFREE(xk1f);SFREE(xk2f);SFREE(xk3f);
	SFREE(dadn);SFREE(wk1);SFREE(wk2);SFREE(wk3);SFREE(xkeqmin);
	SFREE(xkeqmax);SFREE(dkeq);SFREE(domstep);SFREE(domphi);SFREE(r);
	SFREE(crconloc);
      
      }else{

	/* crack propagation due to hcf: 
           if not allowed: stop
           else: accumulate hcf cycles */

	/* a negative dominant step points to HCF crack propagation */
	
	for(i=0;i<nfront;i++){
	  domstep[i]*=-1;
	}
	
	SFREE(xkeq);SFREE(phi);SFREE(xk1);SFREE(xk2);SFREE(xk3);
	SFREE(xkeqf);SFREE(phif);SFREE(xk1f);SFREE(xk2f);SFREE(xk3f);
	SFREE(crconloc);
	
	if(mei[2]==0){
	  printf("HCF crack propagation takes place; the program stops!\n\n");
	  icritic=1;
	}else{
	  printf("HCF crack propagation takes place!\n\n");
	  ncyctot+=ncyc;
	}
	lcf=0;
	  
      }
    }
    
    if(lcf==1){

      /* reset icritic to the default */
      
      icritic=0;
      
      NNEW(dadn,double,nfront);
      NNEW(wk1,double,nfront);
      NNEW(wk2,double,nfront);
      NNEW(wk3,double,nfront);
      NNEW(xkeqmin,double,nfront);
      NNEW(xkeqmax,double,nfront);
      NNEW(dkeq,double,nfront);
      NNEW(r,double,nfront);
      NNEW(domstep,double,nfront);
      NNEW(domphi,double,nfront);
      NNEW(crconloc,double,ncrconst);
    
      FORTRAN(crackrate,(&nfront,ifrontrel,xkeq,phi,ifront,dadn,&ncyc,
			 &icritic,&datarget,crcon,temp,
			 &ncrtem,crconloc,&ncrconst,xk1,xk2,xk3,&nstep,
			 acrack,wk1,wk2,wk3,xkeqmin,xkeqmax,
			 dkeq,domstep,domphi,param,&nparam,&law,&ier,r));
      
      if((icritic>0)||(ier==1)){
	if(icritic>0){
	  printf("WARNING: Exit because a critical value for Keq is reached along the crack front due to the static loading \n");
	}else{
	  printf("WARNING: Exit because an error occurred in crackrate.f \n");
	}
	printf("Number of iteration= %d\n\n" ,(iinc+1));
	
	SFREE(xkeq);SFREE(phi);SFREE(xk1);SFREE(xk2);SFREE(xk3);
	SFREE(crconloc);
	SFREE(kontri);SFREE(ipoed);SFREE(iedg);SFREE(ieled);SFREE(ibounedg);
	SFREE(ibounnod);SFREE(iedno);SFREE(stress);SFREE(ifront);
	SFREE(ifrontrel);SFREE(istartfront);SFREE(iendfront);SFREE(xt);
	SFREE(stressf);SFREE(tempf);
	SFREE(istartcrackfro);SFREE(iendcrackfro);SFREE(temp);
	SFREE(istartcrackbou);SFREE(iendcrackbou);SFREE(r);
	SFREE(isubsurffront);SFREE(acrack);SFREE(dkeq);SFREE(costruc);
	SFREE(dadn);SFREE(wk1);SFREE(wk2);SFREE(wk3);SFREE(xkeqmin);
	SFREE(xkeqmax);SFREE(domstep);SFREE(domphi);
      
	break;
      }

      ncyctot+=ncyc;

      SFREE(xkeq);SFREE(phi);SFREE(xk1);SFREE(xk2);SFREE(xk3);
      SFREE(crconloc);
    }
    
    nkold=*nk;
    NNEW(da,double,nfront);
    RENEW(co,double,3*(*nk+nfront));
    NNEW(ifrontprop,ITG,nfront);
    
    FORTRAN(crackprop,(ifrontrel,ibounnod,domphi,da,co,costruc,nk,xa,
		       xn,&nnfront,istartfront,iendfront,doubleglob,
		       integerglob,isubsurffront,dadn,&ncyc,
		       ifrontprop,&nstep,acrack,acrackglob,&datarget,
		       &ieqspace,iincglob,&iinc));

    SFREE(xn);SFREE(xa);

    /* calculate the mesh characteristic length for each front
       and the new equally spread nodes on the propagated front(s) */

    NNEW(ifronteq,ITG,2*nfront);
    NNEW(istartfronteq,ITG,nnfront);
    NNEW(iendfronteq,ITG,nnfront);

    /* the assumption is made, that no more than 2*nfront equally
       spaced nodes are necessary for the propagated front */
    
    RENEW(co,double,3*(*nk+2*nfront));
    RENEW(acrackglob,double,*nk+2*nfront);
    RENEW(iincglob,ITG,*nk+2*nfront);
    
    if(iinc==0){
      FORTRAN(characteristiclength,(co,istartcrackfro,iendcrackfro,&ncrack,
				    ifront,charlen,&datarget));
    }

    /* next lines of no equal spacing is desired */

    if(ieqspace==0){
    
      RENEW(ifronteq,ITG,nfront);
      memcpy(&ifronteq[0],&ifrontprop[0],sizeof(ITG)*nfront);
      memcpy(&istartfronteq[0],&istartfront[0],sizeof(ITG)*nnfront);
      memcpy(&iendfronteq[0],&iendfront[0],sizeof(ITG)*nnfront);
      nfronteq=nfront;

    }else{
      
    /* next lines of equal spacing is desired */
    
      FORTRAN(eqspacednodes,(co,istartfront,iendfront,&nnfront,
			     ifrontprop,nk,&nfront,ifronteq,charlen,
			     istartfronteq,iendfronteq,&nfronteq,
			     acrackglob,&ier,iendcrackfro,iincglob,
			     &iinc,dnglob,&ncyctot));
      if(ier==1) break;
    }

    RENEW(acrackglob,double,3*(*nk+nfront));
    RENEW(iincglob,ITG,3*(*nk+nfront));
    RENEW(ifronteq,ITG,nfronteq);
    
    /* update the crack mesh topology (generate new crack elements) */
    
    RENEW(lakon,char,8*(*ne+(nfronteq+nfront)));
    RENEW(ipkon,ITG,*ne+(nfronteq+nfront));
    RENEW(ielmat,ITG,mi[2]*(*ne+(nfronteq+nfront)));
    for(i=mi[2]**ne;i<mi[2]*(*ne+(nfronteq+nfront));i++){ielmat[i]=0.;}
    RENEW(kon,ITG,*nkon+9*(*ne+(nfronteq+nfront)));

    FORTRAN(extendmesh,(&nnfront,istartfront,iendfront,ifront,&nkold,ne,
			nkon,lakon,ipkon,kon,isubsurffront,co,ifronteq,
			istartfronteq,iendfronteq,&nfront,&nfronteq));

    RENEW(lakon,char,8**ne);
    RENEW(ipkon,ITG,*ne);
    RENEW(ielmat,ITG,mi[2]**ne);
    RENEW(kon,ITG,*nkon);

    /* storing the local crack results into global nodal fields */
    
    RENEW(wk1glob,double,3*(*nk+nfront));
    RENEW(wk2glob,double,3*(*nk+nfront));
    RENEW(wk3glob,double,3*(*nk+nfront));
    RENEW(xkeqminglob,double,3*(*nk+nfront));
    RENEW(xkeqmaxglob,double,3*(*nk+nfront));
    RENEW(dkeqglob,double,3*(*nk+nfront));
    RENEW(rglob,double,3*(*nk+nfront));
    RENEW(phiglob,double,3*(*nk+nfront));
    RENEW(dadnglob,double,3*(*nk+nfront));
    RENEW(dnglob,double,3*(*nk+nfront));
    RENEW(iincglob,ITG,3*(*nk+nfront));
    RENEW(domstepglob,double,3*(*nk+nfront));
    
    FORTRAN(globalcrackresults,(&nfront,ifront,wk1,wk2,wk3,dkeq,domphi,dadn,
				&ncyctot,wk1glob,wk2glob,wk3glob,dkeqglob,
				phiglob,dadnglob,dnglob,acrack,acrackglob,
				&nstep,xkeqmin,xkeqmax,xkeqminglob,xkeqmaxglob,
				&iinc,iincglob,domstep,domstepglob,r,rglob));
  
    SFREE(kontri);SFREE(ipoed);SFREE(iedg);SFREE(ieled);SFREE(ibounedg);
    SFREE(ibounnod);SFREE(iedno);SFREE(stress);SFREE(ifront);SFREE(ifrontrel);
    SFREE(istartfront);SFREE(iendfront);SFREE(xt);SFREE(stressf);SFREE(tempf);
    SFREE(istartcrackfro);SFREE(iendcrackfro);SFREE(temp);
    SFREE(istartcrackbou);SFREE(iendcrackbou);
    SFREE(isubsurffront);SFREE(acrack);
    SFREE(dkeq);SFREE(da);SFREE(costruc);SFREE(r);
    SFREE(dadn);SFREE(ifrontprop);SFREE(ifronteq);SFREE(istartfronteq);
    SFREE(iendfronteq);SFREE(wk1);SFREE(wk2);SFREE(wk3);SFREE(xkeqmin);
    SFREE(xkeqmax);SFREE(domstep);SFREE(domphi);

    if(icritic>0){
	printf("WARNING: Exit because crack propagation occurred due to the modal loading \n");
	printf("Number of iteration= %d\n\n" ,(iinc+1));
	break;
    }else if(icritic<0){
	printf("WARNING: Exit because no crack propagation anywhere along the crack front \n");
	printf("Number of iteration= %d\n\n" ,(iinc+1));
	break;
    }

  }

  SFREE(crcon);SFREE(integerglob);SFREE(doubleglob);
  if(mei[0]>0){SFREE(integerglobf);SFREE(doubleglobf);}

  for(j=0;j<*ne;j++){
    if(strcmp1(&lakon[8*j+6],"L")==0) lakon[8*j+6]='A';
  }

  NNEW(inum,ITG,*nk);
  for(j=0;j<*nk;j++){
    inum[j]=1;
  }
  
  /* storing the global mesh and the triangulation of the crack(s):
     shells must be stored as 2D-triangles */

  filab[4]='I';
  frd(co,nk,kon,ipkon,lakon,ne,v,stn,inum,nmethod,
      &kode,filab,een,t1,fn,&time,epn,ielmat,matname,enern,xstaten,
      &nstate_,istep,&iinc,ithermal,qfn,&mode,&noddiam,trab,inotr,
      &ntrans,orab,ielorien,&norien,description,ipneigh,neigh,
      mi,sti,vr,vi,stnr,stni,vmax,stnmax,&ngraph,veold,ener,ne,
      cs,set,nset,istartset,iendset,ialset,eenmax,fnr,fni,emn,
      thicke,jobnamec,output,qfx,cdn,&mortar,cdnr,cdni,nmat,
      ielprop,prop,sti);

  /* storing the crack propagation fields in frd-format */

  if(strcmp1(&filab[4611],"KEQ ")==0){
    crackfrd(nk,&ngraph,&noddiam,cs,&kode,inum,nmethod,&time,istep,&iinc,&mode,
	     description,set,nset,istartset,iendset,ialset,jobnamec,output,
	     dkeqglob,wk1glob,wk2glob,wk3glob,phiglob,dadnglob,dnglob,
	     acrackglob,xkeqminglob,xkeqmaxglob,iincglob,domstepglob,
	     rglob);
  }
  
  SFREE(wk1glob);SFREE(wk2glob);SFREE(wk3glob);SFREE(dkeqglob);
  SFREE(phiglob);SFREE(dadnglob);SFREE(dnglob);SFREE(xkeqminglob);
  SFREE(inum);SFREE(charlen);SFREE(acrackglob);SFREE(xkeqmaxglob);
  SFREE(iincglob);SFREE(domstepglob);SFREE(rglob);
  
  *cop=co;*lakonp=lakon;*ipkonp=ipkon;*konp=kon;*ielmatp=ielmat;
  
  return;
}
