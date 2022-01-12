/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2021 Guido Dhondt                     */

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
#ifdef PASTIX
#include "pastix.h"
#endif

/**
 * @brief       *MASSLESS DYNAMIC CONTACT*: Main computation of step for dynamic massless contact
 *
 * calling also functions:
 * + createcontactdofs():
 * + expand_auw:
 * + extract_matrices:
 * + resforccont
 * + detectactivecont
 *
 * @param      au         LOWER triangle of STIFFNESS matrix of size: neq[0]
 * @param      ad         DIAGONAL       of STIFFNESS matrix of size: neq[0]
 * @param      aub        UPPER triangle of MASS      matrix of size: neq[0]
 * @param      adb        DIAGONAL       of MASS      matrix of size: neq[0]
 * @param      jq         Location in field **irow** of the first subdiagonal nonzero in column i (only for symmetric matrices)
 * @param      irow       Row of element i in field au (i.e. au(i))
 * @param      neq        NTOT of equations:
 *                        + neq[0]: mechanical TOTDOF,
 *                        + neq[1]: TOTDOF_mech + TOTDOF_thermal
 *                        + neq[2]: neq[1] + # of single point constraints (only for modal calculations)
 * @param      nzs        projected nonzeros:
 *                        + nzs[0]: sum of projected nonzero mechanical off-diagonal terms
 *                        + nzs[1]: nzs[0]+sum of projected nonzero thermal off-diagonal terms
 *                        + nzs[2]: nzs[1] + sum of nonzero coefficients of SPC DOFs (only for modal calculations)+
 * @param      auw        The auw   for W_b matrix??
 * @param      jqw        The jqw   for W_b matrix??
 * @param      iroww      The iroww for W_b matrix??
 * @param      nzsw       The nzsw   for W_b matrix??
 * @param      islavnode  Slave nodes are stored tie by tie
 * @param      nslavnode  Position in islavnode before the first slave node of tie i.
 * @param      nslavs     The total number of slave nodes
 * @param      imastnode  Master nodes catalogued as in **islavnode** of size nmastnode(ntie+1)
 * @param      nmastnode  The nmastnode(i) points towards the master node within tie i with the largest number
 * @param      ntie       Total number of tie constraints
 * @param      nactdof    Active degrees of freedom are stored in a two-dimensional field.
 *                        + It has as many rows as there are nodes in the model and four columns since each node has
 *                        one temperature degree of freedom and three translational degrees.
 *                        + In C this field is mapped into a one-dimensional field starting
 *                        with the degrees of freedom of node 1, then those of node 2, and so on.
 * @param      mi         Maximums for:
 *                        + mi(1): max # of integration points per element (max over all elements)
 *                        + mi(2): max degree of freedom per node (max over all nodes)
 * @param      vold      vold(0..4,1..nk)   solution field in all nodes
 *                        + 0: temperature
 *                        + 1: displacement in global x-direction
 *                        + 2: displacement in global y-direction
 *                        + 3: displacement in global z-direction
 *                        + 4: static pressure
 * @param      nk         highest node number
 * @param      fext       external mechanical forces in DOF i (due to point loads and
 distributed loads, including centrifugal and gravity loads,
 but excluding temperature loading and displacement loading)
 * @param      isolver    Sparse linear solver selector
 *                        + 'SPOOLES'           0
 *                        + 'ITERATIVESCALING'  2
 *                        + 'ITERATIVECHOLESKY' 3
 *                        + 'SGI'               4
 *                        + 'TAUCS'             5
 *                        + 'PARDISO'           7
 *                        + 'PASTIX'            8
 * @param      iperturb   Perturbation?
 *                        + iperturb(1)
 *                          + = -1     : linear iteration in a nonlinear calculation
 *                          + = 0      : linear
 *                          + = 1      : second order theory for frequency/buckling/Green calculations following a static step (PERTURBATION selected)
 *                          + $ \ge$ 2 : Newton-Raphson iterative procedure is active
 *                          + = 3      : nonlinear material (linear or nonlinear geometric and/or heat transfer)
 *                        + iperturb(2)
 *                          +0 : linear geometric (NLGEOM not selected)
 *                          +1 : nonlinear geometric (NLGEOM
 *
 */

void massless(ITG *kslav,ITG *lslav,ITG *ktot,ITG *ltot,double *au,double *ad,
	      double *auc,double *adc,
	      ITG *jq,ITG *irow,ITG *neq,ITG *nzs,double *auw,
	      ITG *jqw,ITG *iroww,ITG *nzsw,ITG *islavnode,ITG *nslavnode,
	      ITG *nslavs,ITG *imastnode,ITG *nmastnode,ITG *ntie,ITG *nactdof,
	      ITG *mi,double *vold,double *volddof, double *veold,ITG *nk,
	      double *fext,ITG *isolver,ITG *iperturb,double *co,
	      double *springarea,ITG *neqslav,ITG *neqtot,double *qb,
	      double *b ){

  printf("\033[0;32m"); // green
  //   printf("===================> DEBUG: massless.c\n");

  /* determining the RHS of the global system for massless contact */

  ITG *jqwnew=NULL,*irowwnew=NULL,symmetryflag=0,mt=mi[1]+1,ncdim,
    inputformat=0,*iacti=NULL,nacti=0,itranspose,index,i,j,k,kitermax,
    *jqbb=NULL,*irowbb=NULL,*icolbb=NULL,nzsbb,*jqbi=NULL,*irowbi=NULL,
    nzsbi,*jqib=NULL,*irowib=NULL,nzsib;

  double *auwnew=NULL,sigma=0.0,*gapdisp=NULL,*gapnorm=NULL,*cvec=NULL,sum,
    *adbbb=NULL,*aubbb=NULL,*gvec=NULL,*gmatrix=NULL,*qi_kbi=NULL,
    *veolddof=NULL,*pkglob=NULL,mufric,atol,rtol,*aubb=NULL,*adbb=NULL,
    *al=NULL,*alnew=NULL,*eps_al=NULL,*rhs=NULL,*aubi=NULL,
    *auib=NULL;

  /* expanding the matrix Wb according to the number of degrees
     of freedom */

  NNEW(jqwnew,ITG,3**nslavs+1);
  NNEW(auwnew,double,*nzsw);
  NNEW(irowwnew,ITG,*nzsw);

  /* Rearrange the row entries in the Wb matrix, column by column c
     from the order in islavnode and imastnode to the order as
     dictated by nactdof */
  
  FORTRAN(expand_auw,(auw,jqw,iroww,nslavs,auwnew,jqwnew,irowwnew,neqslav,
		      nactdof,mi,ktot,neqtot,islavnode,imastnode));

  memcpy(jqw,jqwnew,sizeof(ITG)*(3**nslavs+1));
  memcpy(auw,auwnew,sizeof(double)**nzsw);
  memcpy(iroww,irowwnew,sizeof(ITG)**nzsw);

  SFREE(jqwnew);SFREE(auwnew);SFREE(irowwnew);

  /* extracting Kbb,Kbi,Kib,Kii from the stiffness matrix */

  NNEW(jqbb,ITG,*neqtot+1);
  NNEW(aubb,double,nzs[0]);
  NNEW(adbb,double,*neqtot);
  NNEW(irowbb,ITG,nzs[0]);
  NNEW(icolbb,ITG,*neqtot);

  NNEW(jqbi,ITG,neq[0]+1);
  NNEW(aubi,double,nzs[0]);
  NNEW(irowbi,ITG,nzs[0]);

  NNEW(jqib,ITG,*neqtot+1);
  NNEW(auib,double,nzs[0]);
  NNEW(irowib,ITG,nzs[0]);

  NNEW(icolbb,ITG,*neqtot);

  /* extracting the submatrices from the global stiffness matrix */
  
  FORTRAN(extract_matrices,(au,ad,jq,irow,neq,aubb,adbb,jqbb,irowbb,neqtot,
			    &nzsbb,aubi,jqbi,irowbi,&nzsbi,auib,jqib,irowib,
			    &nzsib,ktot,icolbb));

  RENEW(aubb,double,nzsbb);
  RENEW(irowbb,ITG,nzsbb);
  RENEW(aubi,double,nzsbi);
  RENEW(irowbi,ITG,nzsbi);
  RENEW(auib,double,nzsib);
  RENEW(irowib,ITG,nzsib);
  SFREE(jqbb);

  /* calculate the residual force in the contact area and store in gapdisp */
  
  NNEW(gapdisp,double,*neqtot);
  NNEW(qi_kbi,double,*neqtot);
  
  FORTRAN(resforccont,(vold,nk,mi,aubi,irowbi,jqbi,neqtot,ktot,fext,gapdisp,
		       auib,irowib,jqib,nactdof,volddof,neq,qi_kbi));

  /* factorize Kbb and premultiply gapdisp with Kbb^{-1} */

  if(*isolver==0){
#ifdef SPOOLES
    spooles_factor_rad(adbb,aubb,adbbb,aubbb,&sigma,icolbb,irowbb,
		       neqtot,&nzsbb,&symmetryflag,&inputformat);

    spooles_solve_rad(gapdisp,neqtot);
#else
    printf("*ERROR in linstatic: the SPOOLES library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }else if(*isolver==7){
#ifdef PARDISO
#else
    printf("*ERROR in linstatic: the PARDISO library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }else if(*isolver==8){
#ifdef PASTIX
#else
    printf("*ERROR in linstatic: the PASTIX library is not linked\n\n");
    FORTRAN(stop,());
#endif
  }
  SFREE(aubb);SFREE(adbb);SFREE(irowbb);SFREE(icolbb);

  NNEW(gapnorm,double,*nslavs);
  NNEW(iacti,ITG,*nslavs);

  /* premultiply g by Wb^T and add g0 => determine active degrees => 
     reduce g to c */
  
  FORTRAN(detectactivecont,(gapnorm,gapdisp,auw,iroww,jqw,nslavs,springarea,
			    iacti,&nacti));

  SFREE(gapdisp);

  printf("---> %i ACTIVE CONTACTS!<---\n",nacti);

  /* reduced the gap dimension to the active dofs */
  
  if (nacti>0){

    /* only normal contact so far */
    
    ncdim=1;

    /* constructing the c-vector of the inclusion equation */
    
    NNEW(cvec,double,nacti*ncdim);
    for(i=0;i<*nslavs;i++){
      if(iacti[i]!=0){
	cvec[iacti[i]-1]=gapnorm[i];
      }
    }

    /* constructing the g-matrix of the inclusion equation */
    
    NNEW(gmatrix,double,nacti*nacti);

    /* calculate G = Wb^T.Kbb^(-1).Wb 
       only for active slave degrees of freedom */

    /* loop over the columns of Wb */
    
    for(i=0;i<*nslavs;i++){
      
      if(iacti[i]!=0) {
	index=3*i;// normal contact index

	NNEW(gvec,double,*neqtot);

	/* Filling the vector of Wb column */
	
	for(j=jqw[index]-1;j<jqw[index+1]-1;j++){
	  gvec[iroww[j]-1]=auw[j];
	}
 
	/* Solving the linear system per column */

	if(*isolver==0){
#ifdef SPOOLES
	  spooles_solve_rad(gvec,neqtot);
#endif
	}else if(*isolver==7){
#ifdef PARDISO
#endif
	}else if(*isolver==8){
#ifdef PASTIX
#endif
	}
	
	/* premultiplying per Wb^T */

	for(j=0;j<*nslavs;++j){
	  if(iacti[j]!=0){
	    index=3*j;
	    sum=0.0;
	    for(k=jqw[index]-1;k<jqw[index+1]-1;k++){
	      sum+=auw[k]*gvec[iroww[k]-1];
	    }
	    gmatrix[(iacti[i]-1)*nacti+(iacti[j]-1)]=sum;
	  }
	}
	SFREE(gvec);
      }
    }

    /* solve the inclusion problem (augmented Lagrange) */
    
    mufric=0.0;
    atol=1.0e-8;
    rtol=1.0e-8;
    kitermax=1000;

    NNEW(al,double,nacti*ncdim);
    NNEW(alnew,double,nacti*ncdim);
    NNEW(eps_al,double,nacti*ncdim);
    NNEW(pkglob,double,*neqtot);
    FORTRAN(auglag_inclusion,
	    (gmatrix,cvec,iacti,&nacti,&ncdim,&mufric,&atol,&rtol,
	     pkglob,&kitermax,auw,jqw,iroww,nslavs,al,
	     alnew,eps_al));
    SFREE(al);SFREE(alnew);SFREE(eps_al);SFREE(gmatrix);SFREE(cvec);

  }else{
    NNEW(pkglob,double,*neqtot);
  }

  SFREE(gapnorm);SFREE(iacti);
  
  /* compute q_b^k = Kbb^{-1}*(Wb*al-qi_kbi+fexb) */

  for(i=0;i<*neqtot;i++){
    qb[i]=pkglob[i]-qi_kbi[i]+fext[ktot[i]-1];
  }

  SFREE(qi_kbi);SFREE(pkglob);

  if(*isolver==0){
#ifdef SPOOLES
    spooles_solve_rad(qb,neqtot);
#endif
  }else if(*isolver==7){
#ifdef PARDISO
#endif
  }else if(*isolver==8){
#ifdef PASTIX
#endif
  }

  /* compute the right hand side of the global equation system */
  
  /* calculate Kii*qi */
  
  NNEW(rhs,double,neq[0]); 
  FORTRAN(op,(&neq[0],volddof,rhs,ad,au,jq,irow));

  /* calculate Kii*qi+Kib*qb */

  itranspose=0;
  FORTRAN(mulmatvec_asym,(auib,jqib,irowib,neqtot,qb,rhs,&itranspose));
  itranspose=1;
  FORTRAN(mulmatvec_asym,(aubi,jqbi,irowbi,&neq[0],qb,rhs,&itranspose));
  
  SFREE(jqbi);SFREE(aubi);SFREE(irowbi);SFREE(jqib);SFREE(auib);SFREE(irowib);

  /* calculate (Mii/Delta_t-Dii/2)*u_i^{k-1/2} */

  /* switch for the velocity from a nodal to a dof representation */
  
  NNEW(veolddof,double,neq[0]);
  for(k=0;k<*nk;++k){
    for(j=0;j<mt;j++){
      if(nactdof[mt*k+j]>0){
        veolddof[nactdof[mt*k+j]-1]=veold[mt*k+j];
      }
    }
  }

  FORTRAN(op,(&neq[0],veolddof,b,adc,auc,jq,irow));

  for(i=0;i<neq[0];++i){b[i]=fext[i]-rhs[i]+b[i];}
  SFREE(rhs);SFREE(veolddof);

  /* clearing memory for the equation solver */
  
  if(*isolver==0){
#ifdef SPOOLES
    spooles_cleanup_rad();
#elif defined(PARDISO)
      pardiso_cleanup_as(neqtot,&symmetryflag);
#elif defined(PASTIX)
#endif
  }

  printf("\033[0m"); // reset color
  return;
}
