!
!     CalculiX - A 3-dimensional finite element program
!     Copyright (C) 1998-2015 Guido Dhondt
!     
!     This program is free software; you can redistribute it and/or
!     modify it under the terms of the GNU General Public License as
!     published by the Free Software Foundation(version 2);
!     
!     
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of 
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
!     GNU General Public License for more details.
!     
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!    
!
!     center of gravity of the projection of the vertices for
!     visibility purposes
!     exact integration for one triangle: routine cubtri
!     if the surfaces are far enough away, one-point integration
!     is used
! 
      subroutine calcview(sideload,vold,co,pmid,e1,e2,e3,
     &     kontri,nloadtr,adview,auview,dist,idist,area,
     &     ntrit,mi,jqrad,irowrad,nzsrad,sidemean,ntria,ntrib)
!     
      implicit none
!     
      character*1 covered(160,160)
!
      character*20 sideload(*)
!     
      integer ntr,i,j,k,l,mi(*),ntria,ntrib,
     &     ncovered,kontri(4,*),nloadtr(*),
     &     i1,j1,istart,iend,jstart,jend,imin,imid,imax,
     &     k1,kflag,idist(*),ndist,i2,i3,ng,idi,idj,ntri,
     &     ix,iy,ntrit,limev,ier,nw,idata(1),ncalls,
     &     irowrad(*),jqrad(*),nzsrad,i0,nzsradv(3)
!     
      real*8 w(239),vold(0:mi(2),*),co(3,*),xn(3),xxn,
     &     pmid(3,*),e3(4,*),e1(3,*),e2(3,*),p1(3),p2(3),p3(3),
     &     x,y,porigin(3),yymin,yymax,xxmin,
     &     xxmid,xxmax,dummy,a(3,3),b(3,3),c(3,3),ddd(3),p31(3),
     &     xx(3),yy(3),ftij,adview(*),auview(*),dint,dir(3),
     &     dirloc(3),dist(*),area(*),dd,p21(3),sidemean,p(3,3),
     &     fform,ver(2,3),epsabs,epsrel,abserr,vj(3,3),unitvec(3,3),
     &     rdata(21),vertex(3,3),vertexl(2,3),factor,argument
!
      external fform
!     
      nzsradv(3)=nzsrad
!     
      ng=160
      dint=2.d0/ng
!
!     factor determines when the numerical integration using cubtri
!     is replaced by a simplified formula using only the center
!     of gravity of one of the triangles. The integration over the
!     other triangle is exact (analytical formula, see
!     "Radiosity: a Programmer's Perspective", by Ian Ashdown, Wiley, 1994)
!     If the distance between the center of gravity of the triangles
!     is larger then factor*the projected sqrt(area) of the triangle on the 
!     hemisphere, the simplified formula is taken
!
      factor=0.d0
!     
      do i=ntria,ntrib
         if(area(i).lt.1.d-20) cycle
!     
!     vertices of triangle i in local coordinates
!     
         i1=kontri(1,i)
         if(i1.eq.0) cycle
         i2=kontri(2,i)
         i3=kontri(3,i)
         do j=1,3
            porigin(j)=co(j,i1)+vold(j,i1)
            p2(j)=co(j,i2)+vold(j,i2)
            p3(j)=co(j,i3)+vold(j,i3)
            p21(j)=p2(j)-porigin(j)
            p31(j)=p3(j)-porigin(j)
         enddo
         ver(1,1)=0.d0
         ver(2,1)=0.d0
         ver(1,2)=dsqrt(p21(1)**2+p21(2)**2+p21(3)**2)
         ver(2,2)=0.d0
         ver(1,3)=p31(1)*e1(1,i)+p31(2)*e1(2,i)+p31(3)*e1(3,i)
         ver(2,3)=p31(1)*e2(1,i)+p31(2)*e2(2,i)+p31(3)*e2(3,i)
!     
         do k=1,3
            unitvec(k,1)=e1(k,i)
            unitvec(k,2)=e2(k,i)
            unitvec(k,3)=e3(k,i)
         enddo
!     
!     checking which triangles face triangle i
!     
         ndist=0
         idi=kontri(4,i)
         do j=1,ntrit
            if((kontri(1,j).eq.0).or.(area(j).lt.1.d-20)) cycle
            if(pmid(1,j)*e3(1,i)+pmid(2,j)*e3(2,i)+
     &           pmid(3,j)*e3(3,i)+e3(4,i).le.sidemean/800.d0) cycle
            if(pmid(1,i)*e3(1,j)+pmid(2,i)*e3(2,j)+
     &           pmid(3,i)*e3(3,j)+e3(4,j).le.sidemean/800.d0) cycle
!     
            idj=kontri(4,j)
            if(sideload(nloadtr(idi))(18:20).ne.
     &           sideload(nloadtr(idj))(18:20)) cycle
!     
            ndist=ndist+1
            dist(ndist)=dsqrt((pmid(1,j)-pmid(1,i))**2+
     &           (pmid(2,j)-pmid(2,i))**2+
     &           (pmid(3,j)-pmid(3,i))**2)
            idist(ndist)=j
         enddo
         if(ndist.eq.0) cycle
!     
!     ordering the triangles
!     
         kflag=2
         call dsort(dist,idist,ndist,kflag)
!     
!     initializing the coverage matrix
!     
         ncovered=0
         do i1=1,ng
            x=((i1-0.5d0)*dint-1.d0)**2
            do j1=1,ng
               y=((j1-0.5d0)*dint-1.d0)**2
               if(x+y.gt.1.d0) then
                  covered(i1,j1)='T'
                  ncovered=ncovered+1
               else
                  covered(i1,j1)='F'
               endif
            enddo
         enddo
!     
         do k1=1,ndist
            j=idist(k1)
!     
!     determining the 2-D projection of the vertices
!     of triangle j
!     
            do l=1,3
               do k=1,3
                  vertex(k,l)=co(k,kontri(l,j))-pmid(k,i)
               enddo
               dd=dsqrt(vertex(1,l)**2+vertex(2,l)**2+
     &              vertex(3,l)**2)
               do k=1,3
                  vertex(k,l)=vertex(k,l)/dd
               enddo
               vertexl(1,l)=vertex(1,l)*e1(1,i)+
     &              vertex(2,l)*e1(2,i)+
     &              vertex(3,l)*e1(3,i)
               vertexl(2,l)=vertex(1,l)*e2(1,i)+
     &              vertex(2,l)*e2(2,i)+
     &              vertex(3,l)*e2(3,i)
            enddo
!     
!     determining the center of gravity of the projected
!     triangle
!     
            do k=1,2
               dirloc(k)=(vertexl(k,1)+vertexl(k,2)+
     &              vertexl(k,3))/3.d0
            enddo
!     
!     determine the direction vector in global coordinates
!     
            do k=1,3
               dir(k)=(pmid(k,j)-pmid(k,i))/dist(k1)
            enddo
!     
!     direction vector in local coordinates of triangle i
!     
            dirloc(3)=dir(1)*e3(1,i)+dir(2)*e3(2,i)+dir(3)*e3(3,i)
!     
!     check whether this direction was already covered
!     
            ix=int((dirloc(1)+1.d0)/dint)+1
            iy=int((dirloc(2)+1.d0)/dint)+1
            if(covered(ix,iy).eq.'T') then
               cycle
            endif
!     
!     if surfaces are close, numerical integration with
!     cubtri is performed
!     
            if(dist(k1).le.factor*dsqrt(area(i))*dirloc(3)) then
!     
!     vertices of triangle j
!     
               do k=1,3
                  do l=1,3
                     vj(l,k)=co(l,kontri(k,j))+vold(l,kontri(k,j))
                  enddo
               enddo
!     
!     formfactor contribution
!     
               epsrel=0.01d0
               epsabs=0.d0
               limev=100
               nw=239
               ncalls=0
!
!              storing common data into field rdata
!
               rdata(1)=vj(1,1)
               rdata(2)=vj(1,2)
               rdata(3)=vj(1,3)
               rdata(4)=vj(2,1)
               rdata(5)=vj(2,2)
               rdata(6)=vj(2,3)
               rdata(7)=vj(3,1)
               rdata(8)=vj(3,2)
               rdata(9)=vj(3,3)
               rdata(10)=unitvec(1,1)
               rdata(11)=unitvec(1,2)
               rdata(12)=unitvec(1,3)
               rdata(13)=unitvec(2,1)
               rdata(14)=unitvec(2,2)
               rdata(15)=unitvec(2,3)
               rdata(16)=unitvec(3,1)
               rdata(17)=unitvec(3,2)
               rdata(18)=unitvec(3,3)
               rdata(19)=porigin(1)
               rdata(20)=porigin(2)
               rdata(21)=porigin(3)
!     
!     max 1000 evaluations for nw=239
!     
               call cubtri(fform,ver,epsrel,limev,ftij,abserr,ncalls,
     &              w,nw,idata,rdata,ier)
               ftij=ftij/2.d0
            endif
!     
!     updating the coverage matrix
!     
            do k=1,3
               p(k,1)=co(k,kontri(1,j))+vold(k,kontri(1,j))-pmid(k,i)
            enddo
            ddd(1)=dsqrt(p(1,1)*p(1,1)+p(2,1)*p(2,1)+p(3,1)*p(3,1))
            do k=1,3
               p1(k)=p(k,1)/ddd(1)
            enddo
            xx(1)=p1(1)*e1(1,i)+p1(2)*e1(2,i)+p1(3)*e1(3,i)
            yy(1)=p1(1)*e2(1,i)+p1(2)*e2(2,i)+p1(3)*e2(3,i)
!     
            do k=1,3
               p(k,2)=co(k,kontri(2,j))+vold(k,kontri(2,j))-pmid(k,i)
            enddo
            ddd(2)=dsqrt(p(1,2)*p(1,2)+p(2,2)*p(2,2)+p(3,2)*p(3,2))
            do k=1,3
               p2(k)=p(k,2)/ddd(2)
            enddo
            xx(2)=p2(1)*e1(1,i)+p2(2)*e1(2,i)+p2(3)*e1(3,i)
            yy(2)=p2(1)*e2(1,i)+p2(2)*e2(2,i)+p2(3)*e2(3,i)
!     
            do k=1,3
               p(k,3)=co(k,kontri(3,j))+vold(k,kontri(3,j))-pmid(k,i)
            enddo
            ddd(3)=dsqrt(p(1,3)*p(1,3)+p(2,3)*p(2,3)+p(3,3)*p(3,3))
            do k=1,3
               p3(k)=p(k,3)/ddd(3)
            enddo
            xx(3)=p3(1)*e1(1,i)+p3(2)*e1(2,i)+p3(3)*e1(3,i)
            yy(3)=p3(1)*e2(1,i)+p3(2)*e2(2,i)+p3(3)*e2(3,i)
!     
            if(dabs(xx(2)-xx(1)).lt.1.d-5) xx(2)=xx(1)+1.d-5
            if(dabs(xx(2)-xx(1)).lt.1.d-5) xx(2)=xx(1)+1.d-5
!     
!     if the surfaces are far enough away, one-point
!     integration is used
!     
            if(dist(k1).gt.factor*dsqrt(area(i))*dirloc(3)) then
               ftij=0.d0
               do k=1,3
                  l=k-1
                  if(l.lt.1) l=3
                  xn(1)=p(2,k)*p(3,l)-p(2,l)*p(3,k)
                  xn(2)=p(3,k)*p(1,l)-p(3,l)*p(1,k)
                  xn(3)=p(1,k)*p(2,l)-p(1,l)*p(2,k)
                  xxn=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
!     
!     argument of dacos must have an absolute value
!     smaller than or equal to 1.d0; due to
!     round-off the value can slightly exceed one
!     and has to be cut-off.
!     
                  argument=
     &                 (p(1,k)*p(1,l)+p(2,k)*p(2,l)+p(3,k)*p(3,l))/
     &                 (ddd(k)*ddd(l))
                  if(dabs(argument).gt.1.d0) then
                     if(argument.gt.0.d0) then
                        argument=1.d0
                     else
                        argument=-1.d0
                     endif
                  endif
                  ftij=ftij+
     &                 (e3(1,i)*xn(1)
     &                 +e3(2,i)*xn(2)
     &                 +e3(3,i)*xn(3))/xxn
     &                 *dacos(argument)
               enddo
               ftij=ftij*area(i)/2.d0
            endif
!     
            idj=kontri(4,j)
            i0=0
            call add_sm_st_as(auview,adview,jqrad,irowrad,
     &           idi,idj,ftij,i0,i0,nzsradv)          
!     
!     determining maxima and minima
!     
            xxmin=2.d0
            xxmax=-2.d0
            do k=1,3
               if(xx(k).lt.xxmin) then
                  xxmin=xx(k)
                  imin=k
               endif
               if(xx(k).gt.xxmax) then
                  xxmax=xx(k)
                  imax=k
               endif
            enddo
!     
            if(((imin.eq.1).and.(imax.eq.2)).or.
     &           ((imin.eq.2).and.(imax.eq.1))) then
               imid=3
               xxmid=xx(3)
            elseif(((imin.eq.2).and.(imax.eq.3)).or.
     &              ((imin.eq.3).and.(imax.eq.2))) then
               imid=1
               xxmid=xx(1)
            else
               imid=2
               xxmid=xx(2)
            endif
!     
!     check for equal x-values
!     
            if(xxmid-xxmin.lt.1.d-5) then
               xxmin=xxmin-1.d-5
               xx(imin)=xxmin
            endif
            if(xxmax-xxmid.lt.1.d-5) then
               xxmax=xxmax+1.d-5
               xx(imax)=xxmax
            endif
!     
!     equation of the straight lines connecting the
!     triangle vertices in the local x-y plane
!     
            a(1,2)=yy(2)-yy(1)
            b(1,2)=xx(1)-xx(2)
            c(1,2)=yy(1)*xx(2)-xx(1)*yy(2)
!     
            a(2,3)=yy(3)-yy(2)
            b(2,3)=xx(2)-xx(3)
            c(2,3)=yy(2)*xx(3)-xx(2)*yy(3)
!     
            a(3,1)=yy(1)-yy(3)
            b(3,1)=xx(3)-xx(1)
            c(3,1)=yy(3)*xx(1)-xx(3)*yy(1)
!     
            a(2,1)=a(1,2)
            b(2,1)=b(1,2)
            c(2,1)=c(1,2)
            a(3,2)=a(2,3)
            b(3,2)=b(2,3)
            c(3,2)=c(2,3)
            a(1,3)=a(3,1)
            b(1,3)=b(3,1)
            c(1,3)=c(3,1)
!     
            istart=int((xxmin+1.d0+dint/2.d0)/dint)+1
            iend=int((xxmid+1.d0+dint/2.d0)/dint)
            do i1=istart,iend
               x=dint*(i1-0.5d0)-1.d0
               yymin=-(a(imin,imid)*x+c(imin,imid))/b(imin,imid)
               yymax=-(a(imin,imax)*x+c(imin,imax))/b(imin,imax)
               if(yymin.gt.yymax) then
                  dummy=yymin
                  yymin=yymax
                  yymax=dummy
               endif
               jstart=int((yymin+1.d0+dint/2.d0)/dint)+1
               jend=int((yymax+1.d0+dint/2.d0)/dint)
               do j1=jstart,jend
                  covered(i1,j1)='T'
               enddo
               ncovered=ncovered+jend-jstart+1
            enddo
!     
            istart=int((xxmid+1.d0+dint/2.d0)/dint)+1
            iend=int((xxmax+1.d0+dint/2.d0)/dint)
            do i1=istart,iend
               x=dint*(i1-0.5d0)-1.d0
               yymin=-(a(imid,imax)*x+c(imid,imax))/b(imid,imax)
               yymax=-(a(imin,imax)*x+c(imin,imax))/b(imin,imax)
               if(yymin.gt.yymax) then
                  dummy=yymin
                  yymin=yymax
                  yymax=dummy
               endif
               jstart=int((yymin+1.d0+dint/2.d0)/dint)+1
               jend=int((yymax+1.d0+dint/2.d0)/dint)
               do j1=jstart,jend
                  covered(i1,j1)='T'
               enddo
               ncovered=ncovered+jend-jstart+1
            enddo
            if(ncovered.eq.ng*ng)exit
!     
         enddo
      enddo
!
      return
      end
!     
!     function to be integrated
!     
      real*8 function fform(x,y,idata,rdata)
!
      implicit none
!     
      integer k,l,idata(1)
!     
      real*8 pint(3),ddd(3),xn(3),vj(3,3),
     &   unitvec(3,3),p(3,3),xxn,x,y,porigin(3),rdata(21)
!
!     retrieving common data from field rdata
!
      vj(1,1)=rdata(1)
      vj(1,2)=rdata(2)
      vj(1,3)=rdata(3)
      vj(2,1)=rdata(4)
      vj(2,2)=rdata(5)
      vj(2,3)=rdata(6)
      vj(3,1)=rdata(7)
      vj(3,2)=rdata(8)
      vj(3,3)=rdata(9)
      unitvec(1,1)=rdata(10)
      unitvec(1,2)=rdata(11)
      unitvec(1,3)=rdata(12)
      unitvec(2,1)=rdata(13)
      unitvec(2,2)=rdata(14)
      unitvec(2,3)=rdata(15)
      unitvec(3,1)=rdata(16)
      unitvec(3,2)=rdata(17)
      unitvec(3,3)=rdata(18)
      porigin(1)=rdata(19)
      porigin(2)=rdata(20)
      porigin(3)=rdata(21)
!
      do k=1,3
         pint(k)=porigin(k)+x*unitvec(k,1)+y*unitvec(k,2)
      enddo
!
      do k=1,3
         p(k,1)=vj(k,1)-pint(k)
      enddo
      ddd(1)=dsqrt(p(1,1)*p(1,1)+p(2,1)*p(2,1)+p(3,1)*p(3,1))
!
      do k=1,3
         p(k,2)=vj(k,2)-pint(k)
      enddo
      ddd(2)=dsqrt(p(1,2)*p(1,2)+p(2,2)*p(2,2)+p(3,2)*p(3,2))
!
      do k=1,3
         p(k,3)=vj(k,3)-pint(k)
      enddo
      ddd(3)=dsqrt(p(1,3)*p(1,3)+p(2,3)*p(2,3)+p(3,3)*p(3,3))
!     
!     calculating the contribution
!     
      fform=0.d0
      do k=1,3
         l=k-1
         if(l.lt.1) l=3
         xn(1)=p(2,k)*p(3,l)-p(2,l)*p(3,k)
         xn(2)=p(3,k)*p(1,l)-p(3,l)*p(1,k)
         xn(3)=p(1,k)*p(2,l)-p(1,l)*p(2,k)
         xxn=dsqrt(xn(1)**2+xn(2)**2+xn(3)**2)
         fform=fform+
     &        (unitvec(1,3)*xn(1)
     &        +unitvec(2,3)*xn(2)
     &        +unitvec(3,3)*xn(3))/xxn
     &        *dacos(
     &        (p(1,k)*p(1,l)+p(2,k)*p(2,l)+p(3,k)*p(3,l))/
     &        (ddd(k)*ddd(l)))
      enddo
!
      return
      end
      
