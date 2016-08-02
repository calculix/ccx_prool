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

void checkconvnet(ITG *icutb, ITG *iin,
		  double *uamt, double *uamf, double *uamp, 
		  double *cam1t, double *cam1f, double *cam1p,
		  double *cam2t, double *cam2f, double *cam2p,
		  double *camt, double *camf, double *camp,
		  ITG *icntrl, double *dtheta, double *ctrl,
                  double *uama,double *cam1a,double *cam2a,double *cama,
                  double *vamt, double *vamf, double *vamp, double *vama,
                  double *qa){
  
  ITG i0,ir,ip,ic,il,ig,ia,idivergence;
  
  double c1t,c1f,c1p,c1a;
  double df,dc,db,dd,ran,can,rap,ea,cae,ral;

  i0=ctrl[0];ir=ctrl[1];ip=ctrl[2];ic=ctrl[3];il=ctrl[4];ig=ctrl[5];ia=ctrl[7];
  df=ctrl[10];dc=ctrl[11];db=ctrl[12];dd=ctrl[16];
  ran=ctrl[18];can=ctrl[19];rap=ctrl[22];
  ea=ctrl[23];cae=ctrl[24];ral=ctrl[25];
  
  /* temperature */
  
  if(*iin<=ip){c1t=0.0001*ran;}
  else{c1t=0.0001*rap;}
  
  /* mass flow */
  
  if(*iin<=ip){c1f=0.0001*ran;}
  else{c1f=0.0001*rap;}
  
  /* pressure */
  
  if(*iin<=ip){c1p=0.0001*ran;}
  else{c1p=0.0001*rap;}
  
  /* geometry */
  
  if(*iin<=ip){c1a=0.0001*ran;}
  else{c1a=0.0001*rap;}
  
  if(*cam1t<*cam2t) {*cam2t=*cam1t;}
  if(*cam1f<*cam2f) {*cam2f=*cam1f;}
  if(*cam1p<*cam2p) {*cam2p=*cam1p;}
  if(*cam1a<*cam2a) {*cam2a=*cam1a;}
  
  /* check for convergence or divergence; 
     each call of radflowload.c corresponds to a thermomechanical
     iteration, let's call the iterations within radflowload.c
     subiterations;
     the convergence check consists of a comparison of the change in
     the latest subiteration with
        - the largest change in the present thermomechanical iteration
        - the largest value in the present thermomechanical iteration */

/*  if(((*camt<=c1t**uamt)||(*camt<1.e-8**vamt))&&
     ((*camf<=c1f**uamf)||(*camf<1.e-8**vamf))&&
     ((*camp<=c1p**uamp)||(*camp<1.e-8**vamp))&&
     ((*cama<=c1a**uama)||(*cama<1.e-8**vama))&&*/
  if(((*camt<=c1t**uamt)||(*camt<c1t**vamt))&&
     ((*camf<=c1f**uamf)||(*camf<c1f**vamf))&&
     ((*camp<=c1p**uamp)||(*camp<c1p**vamp))&&
     ((*cama<=c1a**uama)||(*cama<c1a**vama))&&
     (*iin>3)){
      
      /* increment convergence reached */
      
      printf("      flow network: convergence in gas iteration %" ITGFORMAT " \n\n",*iin);
      *icntrl=1;
      *icutb=0;
  }
  
  else {

      idivergence=0;

      /* divergence based on temperatures */
      
      if((*iin>=20*i0)||(fabs(*camt)>1.e20)){
	  if((*cam1t>=*cam2t)&&(*camt>=*cam2t)&&(*camt>c1t**uamt)){
	      idivergence=1;
	  }
      }

      /* divergence based on the mass flux */
      
      if((*iin>=20*i0)||(fabs(*camf)>1.e20)){
	  if((*cam1f>=*cam2f)&&(*camf>=*cam2f)&&(*camf>c1f**uamf)){
	      idivergence=1;
	  }
      }

      /* divergence based on pressures */
      
      if((*iin>=20*i0)||(fabs(*camp)>1.e20)){
	  if((*cam1p>=*cam2p)&&(*camp>=*cam2p)&&(*camp>c1p**uamp)){
	      idivergence=1;
	  }
      }

      /* divergence based on geometry */
      
      if((*iin>=20*i0)||(fabs(*cama)>1.e20)){
	  if((*cam1a>=*cam2a)&&(*cama>=*cam2a)&&(*cama>c1a**uama)){
	      idivergence=1;
	  }
      }

      /* divergence based on the number of iterations */

      if(*iin>20*ic) idivergence=1;

      /* divergence based on singular matrix or negative pressures */

      if(*iin==0) idivergence=1;
      
      if(idivergence==1){
	  *dtheta=*dtheta*df;
	  printf("\n divergence; the increment size is decreased to %e\n",*dtheta);
	  printf(" the increment is reattempted\n\n");
	  *iin=0;
	  (*icutb)++;
	  if(*icutb>ia){
	      qa[2]=0.25;
	      *icntrl=1;
//	    printf("\n *ERROR: too many cutbacks\n");
//	    FORTRAN(stop,());
	  }
      }else{
	 	  printf("      no convergence\n\n"); 
      }
  }
  return;
}
