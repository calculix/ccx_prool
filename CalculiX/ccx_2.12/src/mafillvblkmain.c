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

#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include "CalculiX.h"

static char *lakonf1;

static ITG num_cpus,*nef1,*ipnei1,*neifa1,*neiel1,*ielfa1,*ielblk1,
    *ifabou1,*nbody1,*nactdohinv1,*icyclic1,*ifatie1;

static double *auv61,*adv1,*bv1,*vfa1,*xxn1,*area1,*vel1,
       *cosa1,*umfa1,*xlet1,*xle1,*gradvfa1,*xxi1,*body1,*volume1,*dtimef1,
       *velo1,*sel1,*xrlfa1,*gamma1,*xxj1,*a11,*a21,*flux1,
       *c1,*adv61,*vela1,*velaa1,*auv31,*bv31;

void mafillvblkmain(ITG *nef,ITG *ipnei,ITG *neifa,ITG *neiel,
             double *vfa,double *xxn,double *area,double *auv,double *adv,
             double *bv,double *vel,double *cosa,
             double *umfa,double *xlet,double *xle,double *gradvfa,
	     double *xxi,double *body,double *volume,
	     ITG *ielfa,char *lakonf,ITG *ifabou,ITG *nbody,
	     double *dtimef,double *velo,
	     double *sel,double *xrlfa,double *gamma,double *xxj,
	     ITG *nactdohinv,double *a1,double *a2,
	     double *flux,ITG *icyclic,double *c,ITG *ifatie,
	     double *adv6,double *auv6,double *auv3,double *vela,
	     double *velaa,double *bv3,ITG *ielblk,ITG *nblk,
	     ITG *istartblk,ITG *iendblk,ITG *nblket,ITG *nblkze){

    ITG i,nrhs=3,im,nstart,m,n,j;
      
    /* variables for multithreading procedure */
    
    ITG sys_cpus,*ithread=NULL;
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

    envloc = getenv("CCX_NPROC_CFD");
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

    if(*nef<num_cpus) num_cpus=*nef;
    
    pthread_t tid[num_cpus];

    /* setting fields to zero */

    DMEMSET(adv,0,*nef,0.);
    DMEMSET(adv6,0,6**nef,0.);
    DMEMSET(auv6,0,6**nef,0.);
    DMEMSET(bv,0,3**nef,0.);


    /* calculating the stiffness and/or mass matrix 
       (symmetric part) */


    nef1=nef;ipnei1=ipnei;neifa1=neifa;neiel1=neiel;vfa1=vfa;xxn1=xxn;
    area1=area;auv61=auv6;adv1=adv;bv1=bv;vel1=vel;cosa1=cosa;umfa1=umfa;
    xlet1=xlet;xle1=xle;gradvfa1=gradvfa;xxi1=xxi;body1=body;volume1=volume;
    ielfa1=ielfa;ifabou1=ifabou;nbody1=nbody;
    velo1=velo;sel1=sel;xxj1=xxj;nactdohinv1=nactdohinv;a11=a1;a21=a2;
    flux1=flux;icyclic1=icyclic;c1=c;ifatie1=ifatie;adv61=adv6;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillv0mt, (void *)&ithread[i]);
    }
    for(i=0;i<num_cpus;i++)pthread_join(tid[i],NULL);
    
    SFREE(ithread);


    /* solving the equations in the xi-direction */


    DMEMSET(auv3,0,3**nef,0.);

    adv1=adv;auv61=auv6;bv1=bv;auv31=auv3;bv31=bv3;ielblk1=ielblk;
    ipnei1=ipnei;neiel1=neiel;vel1=vel;nef1=nef;adv61=adv6;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillv1mt, (void *)&ithread[i]);
    }
    for(i=0;i<num_cpus;i++)pthread_join(tid[i],NULL);
    
    SFREE(ithread);

    m=1;
    for(i=0;i<*nblk;i++){
	n=iendblk[i]-istartblk[i]+1;
	nstart=3*(istartblk[i]-1);
	FORTRAN(tridiagonal_nrhs,(&auv3[nstart],&bv3[nstart],&n,&m,&nrhs));
    }

    for(i=0;i<*nef;i++){
	for(j=0;j<3;j++){
	    vela[*nef+j**nef+i]=(bv3[i*3+j]+vel[*nef+j**nef+i])/2.;
	}
    }


    /* solving the equations in the eta-direction */


    DMEMSET(auv3,0,3**nef,0.);

    adv1=adv;auv61=auv6;bv1=bv;auv31=auv3;bv31=bv3;ielblk1=ielblk;
    ipnei1=ipnei;neiel1=neiel;vel1=vel;nef1=nef;adv61=adv6;vela1=vela;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillv2mt, (void *)&ithread[i]);
    }
    for(i=0;i<num_cpus;i++)pthread_join(tid[i],NULL);
    
    SFREE(ithread);

    for(i=0;i<*nblk;i++){
	m=nblket[i];
	n=iendblk[i]-istartblk[i]+1;
	nstart=3*(istartblk[i]-1);
	FORTRAN(tridiagonal_nrhs,(&auv3[nstart],&bv3[nstart],&n,&m,&nrhs));
    }

    for(i=0;i<*nef;i++){
	for(j=0;j<3;j++){
	    velaa[*nef+j**nef+i]=(bv3[i*3+j]+vel[*nef+j**nef+i])/2.;
	}
    }


    /* solving the equations in the zeta-direction */


    DMEMSET(auv3,0,3**nef,0.); 

    adv1=adv;auv61=auv6;bv1=bv;auv31=auv3;bv31=bv3;ielblk1=ielblk;
    ipnei1=ipnei;neiel1=neiel;vel1=vel;nef1=nef;adv61=adv6;vela1=vela;
    velaa1=velaa;
    
    /* create threads and wait */
    
    NNEW(ithread,ITG,num_cpus);
    for(i=0; i<num_cpus; i++)  {
	ithread[i]=i;
	pthread_create(&tid[i], NULL, (void *)mafillv3mt, (void *)&ithread[i]);
    }
    for(i=0;i<num_cpus;i++)pthread_join(tid[i],NULL);
    
    SFREE(ithread);

    for(i=0;i<*nblk;i++){
	m=nblkze[i];
	n=iendblk[i]-istartblk[i]+1;
	nstart=3*(istartblk[i]-1);
	FORTRAN(tridiagonal_nrhs,(&auv3[nstart],&bv3[nstart],&n,&m,&nrhs));
    }

    for(i=0;i<*nef;i++){
	for(j=0;j<3;j++){
	    vel[*nef+j**nef+i]=0.8*bv3[i*3+j]+0.2*vel[*nef+j**nef+i];
//	    vel[*nef+j**nef+i]=0.0000*bv3[i*3+j]+vel[*nef+j**nef+i];
	}
    }
  
    return;

}

/* subroutine for multithreading of mafillv0 */

void *mafillv0mt(ITG *i){

    ITG nefa,nefb,nefdelta;
    
// ceil -> floor

    nefdelta=(ITG)floor(*nef1/(double)num_cpus);
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
// next line! -> all parallel sections
    if((*i==num_cpus-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(mafillv0,(nef1,ipnei1,neifa1,neiel1,vfa1,xxn1,area1,
	    auv61,adv1,bv1,vel1,cosa1,umfa1,xlet1,xle1,gradvfa1,xxi1,
	    body1,volume1,ielfa1,ifabou1,nbody1,
	    velo1,sel1,xxj1,nactdohinv1,a11,
	    a21,flux1,&nefa,&nefb,icyclic1,c1,ifatie1,adv61));

    return NULL;
}

/* subroutine for multithreading of mafillv1 */

void *mafillv1mt(ITG *i){

    ITG nefa,nefb,nefdelta;
    
// ceil -> floor

    nefdelta=(ITG)floor(*nef1/(double)num_cpus);
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
// next line! -> all parallel sections
    if((*i==num_cpus-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(mafillv1,(adv1,auv61,bv1,auv31,bv31,ielblk1,
	    ipnei1,neiel1,vel1,nef1,&nefa,&nefb,adv61));

    return NULL;
}

/* subroutine for multithreading of mafillv2 */

void *mafillv2mt(ITG *i){

    ITG nefa,nefb,nefdelta;
    
// ceil -> floor

    nefdelta=(ITG)floor(*nef1/(double)num_cpus);
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
// next line! -> all parallel sections
    if((*i==num_cpus-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(mafillv2,(adv1,auv61,bv1,auv31,bv31,ielblk1,
    ipnei1,neiel1,vel1,nef1,&nefa,&nefb,adv61,vela1));

    return NULL;
}

/* subroutine for multithreading of mafillv3 */

void *mafillv3mt(ITG *i){

    ITG nefa,nefb,nefdelta;
    
// ceil -> floor

    nefdelta=(ITG)floor(*nef1/(double)num_cpus);
    nefa=*i*nefdelta+1;
    nefb=(*i+1)*nefdelta;
// next line! -> all parallel sections
    if((*i==num_cpus-1)&&(nefb<*nef1)) nefb=*nef1;

    FORTRAN(mafillv3,(adv1,auv61,bv1,auv31,bv31,ielblk1,
    ipnei1,neiel1,vel1,nef1,&nefa,&nefb,adv61,vela1,velaa1));

    return NULL;
}
