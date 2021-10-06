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
#include "spooles.h"

/**
 * @brief       *MASSLESS DYNAMIC CONTACT*: Main computation of step for dynamic massless contact
 *
 * calling also functions:
 * + createcontactdofs():
 * + expand_auw:
 * + extract_matrices:
 * + detectactivecont1
 * + detectactivecont2
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
 * @param      auwp       The auwp   for W_b matrix??
 * @param      jqwp       The jqwp   for W_b matrix??
 * @param      irowwp     The irowwp for W_b matrix??
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
void massless(double *au, double *ad, double *aub, double *adb, ITG *jq,
              ITG *irow, ITG *neq, ITG *nzs, double **auwp, ITG **jqwp,
              ITG **irowwp, ITG *nzsw, ITG *islavnode, ITG *nslavnode,
              ITG *nslavs, ITG *imastnode, ITG *nmastnode, ITG *ntie,
              ITG *nactdof, ITG *mi, double *vold, ITG *nk,
              double *fext, ITG *isolver, ITG *iperturb, double *co, double *springarea) {

    ITG nzswnew, *jqwnew = NULL, *irowwnew = NULL, *jqw = NULL, *iroww = NULL, nmasts,
            neqbb, *jqbb = NULL, *irowbb = NULL, *kslav = NULL, *lslav = NULL, *ktot = NULL,
            *ltot = NULL, neqslav, neqtot, *jqbi = NULL, *jqib = NULL, *irowbi = NULL,
            *irowib = NULL, nzsbi, nzsib, symmetryflag = 0, inputformat = 0,nrhs=1,
            nzsbb, *icolbb = NULL,*iacti = NULL, nacti=0;

    double *auwnew = NULL, *auw = NULL, *aubb = NULL, *adbb = NULL, *aubi = NULL, *auib = NULL,
            sigma = 0.0, *gapdof = NULL, *gapnorm = NULL, *qtmp = NULL, *cvec = NULL, sum,
            *adbbb = NULL, *aubbb = NULL , *gcontvec=NULL, *gcontfull=NULL;

    auw = *auwp;
    jqw = *jqwp;
    iroww = *irowwp;
    nmasts = nmastnode[*ntie];

    /* catalogue slave dofs in kslav,lslav and
       slave and master dofs in ktot,ltot */

    NNEW(kslav, ITG, 3 * *nslavs);
    NNEW(lslav, ITG, 3 * *nslavs);
    NNEW(ktot, ITG, 3 * *nslavs + 3 * nmasts);
    NNEW(ltot, ITG, 3 * *nslavs + 3 * nmasts);

//  Create set of slave and slave+master contact DOFS. Sorted.
    FORTRAN(create_contactdofs, (kslav, lslav, ktot, ltot, nslavs, islavnode, &nmasts,
            imastnode, nactdof, mi, &neqslav, &neqtot)); // TODO: nactdof, mi, to declare

    /* expanding the matrix Wb according to the number of degrees
       of freedom */

    NNEW(jqwnew, ITG, 3 * *nslavs + 1);
    NNEW(auwnew, double, 3 * *nzsw);
    NNEW(irowwnew, ITG, 3 * *nzsw); //TODO: declaring irowwnew?? typo?

// Expand contact direction matrix W_b from nodal size to DOF size (3x)
    FORTRAN(expand_auw, (auw, jqw, iroww, nslavs, nzsw, auwnew, jqwnew, irowwnew,
            &nzswnew, &neqslav, lslav, ntie, nactdof, mi, ktot, &neqtot,
            islavnode, nslavnode, imastnode));

    *nzsw = nzswnew;

    RENEW(jqw, ITG, 3 * *nslavs + 1);
    RENEW(auw, double, *nzsw);
    RENEW(iroww, ITG, *nzsw);

    memcpy(jqw, jqwnew, sizeof(ITG) * (3 * *nslavs + 1));
    memcpy(auw, auwnew, sizeof(double) * *nzsw); // TODO: aunew auwnew?? declare
    memcpy(iroww, irowwnew, sizeof(ITG) * *nzsw);

    SFREE(jqwnew);
    SFREE(auwnew);
    SFREE(irowwnew);

    /* extracting Kbb,Kbi,Kib,Kii from the stiffness matrix */

    NNEW(jqbb, ITG, neqtot + 1);
    NNEW(aubb, double, nzs[0]); //  QUESTION: here we are allocation for the size of the  whole K matrix?
    NNEW(adbb, double, neqtot);
    NNEW(irowbb, ITG, nzs[0]);

    NNEW(jqbi, ITG, neq[0] + 1);
    NNEW(aubi, double, nzs[0]);
    NNEW(irowbi, ITG, nzs[0]);

    NNEW(jqib, ITG, neqtot + 1);// QUESTION: do we need also the IB partition? can we just use the transposed BI?
    NNEW(auib, double, nzs[0]);
    NNEW(irowib, ITG, nzs[0]);

    NNEW(icolbb, ITG, neqtot);
    // NNEW(icontactv, ITG, neqtot + 1); // is this size good?

// extracting the submatrices from the global stiffness matrix
    FORTRAN(extract_matrices, (au, ad, jq, irow, neq, nzs, aubb, adbb, jqbb, irowbb, &neqtot,
            &nzsbb, islavnode, nslavs, imastnode, &nmasts, aubi, jqbi,
            irowbi, &nzsbi, auib, jqib, irowib, &nzsib, kslav, ktot, icolbb));

    RENEW(aubb, double, nzsbb);
    RENEW(irowbb, ITG, nzsbb);
    // RENEW(icolbb, ITG, neqtot);  /* length of each column  icolbb(i)=jqbb(i+1)-jqbb(i) */
    RENEW(aubi, double, nzsbi);
    RENEW(irowbi, ITG, nzsbi);
    RENEW(auib, double, nzsib);
    RENEW(irowib, ITG, nzsib);

    /* factorize kbb */
    NNEW(gapdof, double, neqtot);
    NNEW(gapnorm, double, *nslavs);

    MNEW(qtmp, double, neq[0]);

        FORTRAN(detectactivecont1, (vold, nk, mi, aubi, irowbi, jqbi, &neqtot, fext, aubb, adbb, ltot,
            irowbb, jqbb, auw, iroww, jqw, &neqslav, gapdof,
            auib, irowib, jqib, icolbb, nactdof, qtmp, neq, co));

    /* call spooles_solve..... */

    if (*isolver == 0) {
#ifdef SPOOLES
        // spooles(adbb,aubb,adbbb,aubbb,&sigma,gapdof,icolbb,irowbb,&neqtot,&nzsbb,&symmetryflag,
        //         &inputformat,&nzsbb);
        spooles_factor(adbb,aubb,adbbb,aubbb,&sigma,icolbb,irowbb,
                       &neqtot,&nzsbb,&symmetryflag, &inputformat,&nzsbb);

        spooles_solve(gapdof,&neqtot);
#else
        printf("*ERROR in linstatic: the SPOOLES library is not linked\n\n");
        FORTRAN(stop, ());
#endif
    } else if ((*isolver == 2) || (*isolver == 3)) {
        preiter(ad, &au, gapnorm, &icolbb, &irow, neq, nzs, isolver, iperturb);
    } else if (*isolver == 4) {
#ifdef SGI
        token=1;
            sgapdof_main(ad,au,adb,aub,&sigma,g,icolbb,irow,neq,nzs,token);
#else
        printf("*ERROR in linstatic: the SGI library is not linked\n\n");
        FORTRAN(stop, ());
#endif
    } else if (*isolver == 5) {
#ifdef TAUCS
        tau(ad,&au,adb,aub,&sigma,g,icolbb,&irow,neq,nzs);
#else
        printf("*ERROR in linstatic: the TAUCS library is not linked\n\n");
        FORTRAN(stop, ());
#endif
    } else if (*isolver == 7) {
#ifdef PARDISO
        pardiso_main(adbb,aubb,adbbb,aubbb,&sigma,gapdof,icolbb,irowbb,&neqtot,nzsbb,
                 &symmetryflag,&inputformat,jqbb,&nzs[2],&nrhs); // TODO check PARDISO!
#else
        printf("*ERROR in linstatic: the PARDISO library is not linked\n\n");
        FORTRAN(stop, ());
#endif
    } else if (*isolver == 8) {
#ifdef PASTIX
        pastix_main(ad,au,adb,aub,&sigma,gapdof,icolbb,irow,neq,nzs,
                &symmetryflag,&inputformat,jq,&nzs[2],&nrhs);
#else
        printf("*ERROR in linstatic: the PASTIX library is not linked\n\n");
        FORTRAN(stop, ());
#endif
    }

    NNEW(iacti, ITG, *nslavs);
    /* premultiply g by Wb^T and add g0 => determine active degrees => reduce g to c */
    FORTRAN(detectactivecont2, (gapnorm, gapdof, auw, iroww, jqw, &neqtot, nslavs, springarea, iacti,&nacti));


// ***************** G ("DELASSUS") MATRIX + Cvector  NORMAL CONTACT ***********************++
int ncdim = 1;


// initializing cvec with values of gnorm
  NNEW(cvec,double,nacti);
  for(int i=0;i<*nslavs;i++){
    if(iacti[i]!=0){
      cvec[iacti[i]-1]=gapnorm[i];
    }
  }

   NNEW(gcontfull, double, nacti *  nacti);

//   calculate Wb^T.Kbb^(-1).Wb
    int indexnorm   ;
    for (int i = 0; i < *nslavs; i++) {
        if(iacti[i]!=0){
            indexnorm= 3*i ; // normal contact index
            //initialize in zeros vector
            NNEW(gcontvec, double, neqtot);
            for(int j=jqw[indexnorm]-1;j<jqw[indexnorm+1]-1;j++){
                gcontvec[iroww[j]-1]=auw[j];
            }
        }

        spooles_solve(gcontvec,&neqtot); // TODO add other solvers!!!!

        for (int j = 0; j < *nslavs; ++j) {
            if(iacti[j]!=0){
                indexnorm= 3*j ;
                sum=0.0;
                for (int k = jqw[indexnorm]-1; k < jqw[indexnorm+1]-1 ; k++) {
                    sum+=auw[k]*gcontvec[iroww[k]-1];
                }
                gcontfull[i*nacti+j]=sum;
            }
        }
    }

    FORTRAN(relaxval_al, (gcontfull,&nacti, &ncdim )); // last parameter: nCdim: 1:NORMALonly, 3:NTT


// conttype: contact type
int conttype = 1 ;//'NORMAL'
double mufric = 0.0;
double atol = 1.0e-8;
double rtol = 1.0e-8;
double *pkvec = NULL;
int kitermax = 1000;
int timek = 0.0;

FORTRAN( auglag_inclusion ,(conttype, gcontfull,&nacti, &ncdim, mufric, atol,rtol, pkvec, kitermax, timek  ));



    if (*isolver == 0) {
#ifdef SPOOLES
        spooles_cleanup();
#endif
    }

    SFREE(kslav);
    SFREE(lslav);
    SFREE(ktot);
    SFREE(ltot);

    *auwp = auw;
    *jqwp = jqw;
    *irowwp = iroww; // QUESTION: why ?

    SFREE(gcontvec);
    SFREE(gcontfull);

    return;
}
