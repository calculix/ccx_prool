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

void frd_se(double *co,ITG *nk,double *stn,ITG *inum,ITG *nmethod,
         ITG *kode,char *filab,double *fn,double *time,ITG *nstate_,
	 ITG *istep,ITG *iinc,ITG *mode,ITG *noddiam,char *description,
	 ITG *mi,ITG *ngraph,ITG *ne,double *cs,char *set,ITG *nset,
	 ITG *istartset,ITG *iendset,ITG *ialset,double *thicke,
	 char *jobnamec,char *output,double *dgdxglob,ITG *iobject,
	 char *objectset){
	 
     /* stores the results in frd format

     iselect selects which nodes are to be stored:
          iselect=-1 means only those nodes for which inum negative
                     ist, i.e. network nodes
          iselect=+1 means only those nodes for which inum positive
                     ist, i.e. structural nodes
          iselect=0  means both of the above */
  
  FILE *f1;
  
  char c[2]="C",m1[4]=" -1",m2[4]=" -2",m3[4]=" -3",
    p0[6]="    0",p1[6]="    1",p2[6]="    2",p3[6]="    3",p4[6]="    4",
    p5[6]="    5",p6[6]="    6",p7[6]="    7",p8[6]="    8",p9[6]="    9",
    p10[6]="   10",p11[6]="   11",
    p12[6]="   12", fneig[132]="",date[8],clock[10],newdate[21],newclock[9],
    material[59]="                                                          ",
    text[2]=" ";

  static ITG icounter=0,nkcoords,iaxial;

  ITG null,one,i,j,k,indexe,nemax,nlayer,noutloc,iset,iselect,ncomp,nope,
      nodes,ifield[7],nfield[2],icomp[7],ifieldstate[*nstate_],two,three,
      icompstate[*nstate_],ip0=0,ip1=1,ip2=2,ip3=3,ip4=4,ip5=5,ip6=6,ip7=7,
      ip8=8,ip9=9,ip10=10,ip11=11,ip12=12,imat,nelout,
      nterms,nout,noutplus,noutmin,mt=mi[1]+1;

  ITG ncomptensor=2,ifieldtensor[4]={1,1},icomptensor[2]={0,1},
      nfieldtensor[2]={2,0};
      
  int iw;

  float ifl;

  double pi,oner;

  strcpy(fneig,jobnamec);
  strcat(fneig,".frd");

  if((f1=fopen(fneig,"ab"))==NULL){
    printf("*ERROR in frd: cannot open frd file for writing...");
    exit(0);
  }

  pi=4.*atan(1.);
  null=0;
  one=1;two=2;three=3;
  oner=1.;

  /* determining nout, noutplus and noutmin 
              nout: number of structural and network nodes
              noutplus: number of structural nodes
              noutmin: number of network nodes */

  if(*nmethod!=0){
      nout=0;
      noutplus=0;
      noutmin=0;
      for(i=0;i<*nk;i++){
	  if(inum[i]==0) continue;
	  nout++;
	  if(inum[i]>0) noutplus++;
	  if(inum[i]<0) noutmin++;   
      }
  }else{
      nout=*nk;
  }  
  
  nkcoords=*nk;
  
  /* storing the sensitivities in the nodes */
  
  iselect=1;
  
  frdset(&filab[4002],set,&iset,istartset,iendset,ialset,
	 inum,&noutloc,&nout,nset,&noutmin,&noutplus,&iselect,
	 ngraph);
  
  frdheader(&icounter,&oner,time,&pi,noddiam,cs,&null,mode,
	    &noutloc,description,kode,nmethod,f1,output,istep,iinc); 
  
  if(strcmp1(&objectset[*iobject*243],"SHAPEENERGY")==0){
      fprintf(f1," -4  SENENER     2    1\n");
  }else if(strcmp1(&objectset[*iobject*243],"MASS")==0){
      fprintf(f1," -4  SENMASS     2    1\n");
  }else if(strcmp1(&objectset[*iobject*243],"DISPLACEMENT")==0){
      fprintf(f1," -4  SENDISP     2    1\n");
  }else if(strcmp1(&objectset[*iobject*243],"EIGENFREQUENCY")==0){
      fprintf(f1," -4  SENFREQ     2    1\n");
  }

  fprintf(f1," -5  DFDN        1    1    1    0\n");
  fprintf(f1," -5  DFDNFIL     1    1    2    0\n");
  
  frdselect(&dgdxglob[2**nk**iobject],dgdxglob,&iset,&nkcoords,inum,m1,istartset,
	    iendset,ialset,ngraph,&ncomptensor,ifieldtensor,icomptensor,
	    nfieldtensor,&iselect,m2,f1,output,m3);
  
  fclose(f1);
  return;
  
}
