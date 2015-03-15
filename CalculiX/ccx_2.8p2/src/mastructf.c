/*     CalculiX - A 3-dimensional finite element program                 */
/*              Copyright (C) 1998-2015 Guido Dhondt                          */

/*     This program is free software; you can redistribute it and/or     */
/*     modify it under the terms of the GNU General Public License as    */
/*     published by the Free Software Foundation(version 2);    */
/*                                                                       */

/*     This program is distributed in the hope that it will be useful,   */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of    */ 
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the      */
/*     GNU General Public License for more details.                      */

/*     You should have received a copy of the GNU General Public License */
/*     along with this program; if not, write to the Free Software       */
/*     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.         */

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "CalculiX.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

void mastructf(ITG *nk,ITG *kon,ITG *ipkon,char *lakon,ITG *ne,
	       ITG *icol,ITG *jq, ITG **mast1p, ITG **irowp,
	       ITG *isolver, ITG *neq,ITG *ipointer, ITG *nzs,
               ITG *ipnei,ITG *neiel,ITG *mi){

  ITG i,j,k,l,index,idof1,idof2,node1,isubtract,nmast,ifree=0,istart,istartold,
      nzs_,kflag,isize,*mast1=NULL,*irow=NULL,neighbor,mt=mi[1]+1,numfaces;

  /* the indices in the comments follow FORTRAN convention, i.e. the
     fields start with 1 */

  mast1=*mast1p;irow=*irowp;

  kflag=2;
  nzs_=*nzs;

  *neq=*ne;

  /* determining the nonzero locations */

  for(i=0;i<*ne;i++){
      idof1=i+1;
      if(strcmp1(&lakon[8*i+3],"8")==0){
	  numfaces=6;
      }else if(strcmp1(&lakon[8*i+3],"6")==0){
	  numfaces=5;
      }else{
	  numfaces=4;
      }

      index=ipnei[i];
      insert(ipointer,&mast1,&irow,&idof1,&idof1,&ifree,&nzs_);
      for(j=0;j<numfaces;j++){
	  neighbor=neiel[index+j];
	  if(neighbor==0) continue;
	  idof2=neighbor;
	  insert(ipointer,&mast1,&irow,&idof1,&idof2,&ifree,&nzs_);
      }

  }
  
  /*   storing the nonzero nodes in the SUPERdiagonal columns:
       mast1 contains the row numbers,
       irow the column numbers  */
  
  for(i=0;i<*neq;++i){
      if(ipointer[i]==0){
	  printf("*ERROR in mastructf: zero column\n");
	  printf("       element=%" ITGFORMAT "\n",i+1);
	  FORTRAN(stop,());
      }
      istart=ipointer[i];
      while(1){
	  istartold=istart;
	  istart=irow[istart-1];
	  irow[istartold-1]=i+1;
	  if(istart==0) break;
      }
  }
  
  if(*neq==0){
      printf("\n*WARNING: no degrees of freedom in the model\n\n");
  }
  
  nmast=ifree;
  
  /* summary */
  
  printf(" number of equations\n");
  printf(" %" ITGFORMAT "\n",*neq);
  printf(" number of nonzero lower triangular matrix elements\n");
  printf(" %" ITGFORMAT "\n",nmast-(*neq));
  printf("\n");
  
/* switching from a SUPERdiagonal inventory to a SUBdiagonal one:
   since the nonzeros are located in symmetric positions mast1
   can be considered to contain the column numbers and irow the
   row numbers; after sorting mast1 the following results:
   
   - irow contains the row numbers of the SUBdiagonal
   nonzero's, column per column
   - mast1 contains the column numbers
   
   Furthermore, the following fields are determined:       
   
   - icol(i)=# SUBdiagonal nonzero's in column i
   - jq(i)= location in field irow of the first SUBdiagonal
   nonzero in column i  */
  
/* ordering the column numbers in mast1 */
  
  FORTRAN(isortii,(mast1,irow,&nmast,&kflag));
  
/* filtering out the diagonal elements and generating icol and jq */
  
  isubtract=0;
  for(i=0;i<*neq;++i){icol[i]=0;}
  k=0;
  for(i=0;i<nmast;++i){
      if(mast1[i]==irow[i]){++isubtract;}
      else{
	  mast1[i-isubtract]=mast1[i];
	  irow[i-isubtract]=irow[i];
	  if(k!=mast1[i]){
	      for(l=k;l<mast1[i];++l){jq[l]=i+1-isubtract;}
	      k=mast1[i];
	  }
	  ++icol[k-1];
      }
  }
  nmast=nmast-isubtract;
  for(l=k;l<*neq+1;++l){jq[l]=nmast+1;}
  
/* sorting the row numbers within each column */
  
  for(i=0;i<*neq;++i){
      if(jq[i+1]-jq[i]>0){
	  isize=jq[i+1]-jq[i];
	  FORTRAN(isortii,(&irow[jq[i]-1],&mast1[jq[i]-1],&isize,&kflag));
      }
  }
  
  *nzs=jq[*neq]-1;
  
  *mast1p=mast1;*irowp=irow;
  
  return;
  
}

