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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include "CalculiX.h"

static ITG *irow1,*jq1,*num_cpus1,*neq1;

static double *aub1,*adl1,*b1,*sol1,*aux1;

void solveeqmain(double *aub,double *adl,double *b,double *sol,ITG *irow,
		 ITG *jq,ITG *neq,ITG *maxit,ITG *num_cpus){
  
  ITG i,j;
  
  double *aux=NULL;
  
  /* variables for multithreading procedure */
  
  ITG *ithread=NULL;

  pthread_t tid[*num_cpus];

  for(i=0;i<*neq;i++){
    sol[i]*b[i]*adl[i];
  }

  if(*maxit==1) return;

  for(j=2;j<=*maxit;j++){
  
    /* inter/extrapolation of v at the center of the elements
       to the center of the faces */

    aub1=aub;adl1=adl;b1=b;sol1=sol;aux1=aux;irow1=irow;jq1=jq;neq1=neq;
  
    /* create threads and wait */
  
    NNEW(ithread,ITG,*num_cpus);
    for(i=0; i<*num_cpus; i++)  {
      ithread[i]=i;
      pthread_create(&tid[i], NULL, (void *)solveeqparmt, (void *)&ithread[i]);
    }
    for(i=0; i<*num_cpus; i++)  pthread_join(tid[i], NULL);
  
    SFREE(ithread);
  }

}

/* subroutine for multithreading of solveeqpar */

void *solveeqparmt(ITG *i){

  ITG neqa,neqb,neqdelta;
    
  neqdelta=(ITG)ceil(*neq1/(double)*num_cpus1);
  neqa=*i*neqdelta+1;
  neqb=(*i+1)*neqdelta;
  if(neqb>*neq1) neqb=*neq1;
  
  FORTRAN(solveeqpar,(aub1,adl1,b1,sol1,aux1,irow1,jq1,
			  &neqa,&neqb));

  return NULL;
}
