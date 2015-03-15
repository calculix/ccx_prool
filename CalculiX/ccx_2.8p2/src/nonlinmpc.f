!
!     CalculiX - A 3-dimensional finite element program
!              Copyright (C) 1998-2015 Guido Dhondt
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
      subroutine nonlinmpc(co,vold,ipompc,nodempc,coefmpc,labmpc,
     &  nmpc,ikboun,ilboun,nboun,xbounact,aux,iaux,maxlenmpc,ikmpc,
     &  ilmpc,icascade,kon,ipkon,lakon,ne,reltime,newstep,xboun,fmpc,
     &  iit,idiscon,ncont,trab,ntrans,ithermal,mi)
!
!     updates the coefficients in nonlinear MPC's
!
      implicit none
!
c      logical isochoric
!
      character*8 lakon(*)
      character*20 labmpc(*),label
!
      integer ipompc(*),nodempc(3,*),irefnode,irotnode,idir,
     &  nmpc,index,ii,inode,node,id,ikboun(*),ilboun(*),nboun,
     &  i,j,k,idof,na,nb,nc,np,i1,i2,i3,iaux(*),maxlenmpc,n,
     &  l,m,lmax,mmax,ikmpc(*),ilmpc(*),icascade,neigh(7,8),
     &  mpc,kon(*),ipkon(*),indexe,ne,idofrem,idofins,nmpc0,nmpc01,
     &  newstep,iit,idiscon,ncont,iexpnode,indexexp,nmpcdif,ntrans,
     &  nodei,noded,lathyp(3,6),inum,ndir,number,ithermal,mi(*),
     &  newknot,indexexp1,indexexp2,indexexp3,idim
!
      real*8 co(3,*),coefmpc(*),vold(0:mi(2),*),c(3,3),dc(3,3,3),ww,
     &  e(3,3,3),d(3,3),w(3),f(3,3),c1,c2,c3,c4,c5,c6,xbounact(*),
     &  xboun(*),fmpc(*),expan,dd,a11,a12,a13,a21,a22,a23,a31,a32,a33,
     &  b11,b12,b13,b21,b22,b23,b31,b32,b33,aux(*),const,e1(3),e2(3),
     &  ddmax,a(3,3),b(3,3),xj,xi,et,ze,xlag(3,20),xeul(3,20),t1(3),
     &  coloc(3,8),reltime,csab(7),trab(7,*),pd(3),pi(3),e11(3,3),
     &  ad(3,3),ai(3,3),e22(3,3),e12(3,3),ru(3,3),ru1(3,3),
     &  ru2(3,3),ru3(3,3),u1(3,3),u2(3,3),u3(3,3),dcu(3,3,3),u(3,3),
     &  xi1,xi2,xi3,dco,dsi,dco2,dsi2
!
      data d /1.,0.,0.,0.,1.,0.,0.,0.,1./
      data e /0.,0.,0.,0.,0.,-1.,0.,1.,0.,
     &        0.,0.,1.,0.,0.,0.,-1.,0.,0.,
     &        0.,-1.,0.,1.,0.,0.,0.,0.,0./
      data neigh /1,9,2,12,4,17,5,2,9,1,10,3,18,6,
     &            3,11,4,10,2,19,7,4,11,3,12,1,20,8,
     &            5,13,6,16,8,17,1,6,13,5,14,7,18,2,
     &            7,15,8,14,6,19,3,8,15,7,16,5,20,4/
      data coloc /-1.,-1.,-1.,1.,-1.,-1.,1.,1.,-1.,-1.,1.,-1.,
     &            -1.,-1.,1.,1.,-1.,1.,1.,1.,1.,-1.,1.,1./
!     
!     latin hypercube positions in a 3 x 3 matrix
!     
      data lathyp /1,2,3,1,3,2,2,1,3,2,3,1,3,1,2,3,2,1/
!
c      irotnode=0
      irotnode=0
      irefnode=0
      if((icascade.eq.1).and.(newstep.ne.1).and.(ncont.eq.0)) icascade=0
c      isochoric=.false.
!
      ii=0
      loop: do
         ii=ii+1
         if(ii.gt.nmpc) exit
         if(labmpc(ii)(1:5).eq.'RIGID') then
!
            index=ipompc(ii)
            inode=nodempc(1,index)
            idir=nodempc(2,index)
            coefmpc(index)=1.d0
!
            index=nodempc(3,index)
            irefnode=nodempc(1,index)
            coefmpc(index)=-1.d0
!
            index=nodempc(3,index)
            node=nodempc(1,index)
!
!           check whether the rotational node is the same as in
!           the last rigid body MPC
!
            if(node.ne.irotnode) then
               irotnode=node
               w(1)=vold(1,node)
               w(2)=vold(2,node)
               w(3)=vold(3,node)
c               write(*,*) 'w   ',w(1),w(2),w(3)
               ww=dsqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
!
               c1=dcos(ww)
               if(ww.gt.1.d-10) then
                  c2=dsin(ww)/ww
               else
                  c2=1.d0
               endif
               if(ww.gt.1.d-5) then
                  c3=(1.d0-c1)/ww**2
               else
                  c3=0.5d0
               endif
!
!              rotation matrix c
!
               do i=1,3
                  do j=1,3
                     c(i,j)=c1*d(i,j)+
     &                  c2*(e(i,1,j)*w(1)+e(i,2,j)*w(2)+e(i,3,j)*w(3))+
     &                  c3*w(i)*w(j)
                  enddo
               enddo
!
               c4=-c2
               if(ww.gt.0.00464159) then
                  c5=(ww*dcos(ww)-dsin(ww))/ww**3
               else
                  c5=-1.d0/3.d0
               endif
               if(ww.gt.0.0031623) then
                  c6=(ww*dsin(ww)-2.d0+2.d0*dcos(ww))/ww**4
               else
                  c6=-1.d0/12.d0
               endif
!
!              derivative of the rotation matrix c with respect to
!              the rotation vector w
!
               do i=1,3
                  do j=1,3
                     do k=1,3
                        dc(i,j,k)=c4*w(k)*d(i,j)+
     &                            c5*w(k)*(e(i,1,j)*w(1)+
     &                                     e(i,2,j)*w(2)+e(i,3,j)*w(3))+
     &                            c2*e(i,k,j)+
     &                            c6*w(k)*w(i)*w(j)+
     &                            c3*(d(i,k)*w(j)+d(j,k)*w(i))
                     enddo
                  enddo
               enddo
!
!              dummy variable
!
               do i=1,3
                  do j=1,3
c                     f(i,j)=c(i,j)-d(i,j)-dc(i,j,1)*w(1)-dc(i,j,2)*w(2)-
c     &                                    dc(i,j,3)*w(3)
                     f(i,j)=c(i,j)-d(i,j)
                  enddo
               enddo
            endif
!
!           determining the coefficients of the rotational degrees
!           of freedom
!
            coefmpc(index)=dc(idir,1,1)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,1)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,1)*(co(3,irefnode)-co(3,inode))
!
            index=nodempc(3,index)
            coefmpc(index)=dc(idir,1,2)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,2)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,2)*(co(3,irefnode)-co(3,inode))
!
            index=nodempc(3,index)
            coefmpc(index)=dc(idir,1,3)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,3)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,3)*(co(3,irefnode)-co(3,inode))
!
!           determining the nonhomogeneous part
!
            index=nodempc(3,index)
            coefmpc(index)=1.d0
!
!           old value of the nonhomogeneous term must be zero
!
            vold(nodempc(2,index),nodempc(1,index))=0.d0
            idof=8*(nodempc(1,index)-1)+nodempc(2,index)
            call nident(ikboun,idof,nboun,id)
            xbounact(ilboun(id))=f(idir,1)*(co(1,irefnode)-co(1,inode))+
     &           f(idir,2)*(co(2,irefnode)-co(2,inode))+
     &           f(idir,3)*(co(3,irefnode)-co(3,inode))-
     &           vold(idir,irefnode)+vold(idir,inode)
!
         elseif(labmpc(ii)(1:4).eq.'KNOT') then
!
!           dependent node
!
            index=ipompc(ii)
            inode=nodempc(1,index)
            idir=nodempc(2,index)
            coefmpc(index)=1.d0
!
!           translation node
!
            index=nodempc(3,index)
            node=nodempc(1,index)
            coefmpc(index)=-1.d0
!
!           check whether knot is the same as in the previous MPC
!
            if(node.ne.irefnode) then
               newknot=1
               irefnode=node
            else
               newknot=0
            endif
!
            read(labmpc(ii)(5:5),'(i1)') idim
!
!           expansion node
!
            index=nodempc(3,index)
            iexpnode=nodempc(1,index)
!
            if((idim.eq.1).or.(idim.eq.3)) then
!
!              nodes of knot lie on a straight line (1 term in MPC)
!
               indexexp=index
            elseif(idim.eq.2) then
!
!              node of knot lie in a plane (3 terms in MPC)
!
               indexexp1=index
               index=nodempc(3,index)
!
               indexexp2=index
               index=nodempc(3,index)
!
               indexexp3=index
!
            endif
!
            if(newknot.eq.1) then
               if((idim.eq.1).or.(idim.eq.3)) then
                  expan=1.d0+vold(1,iexpnode)
               elseif(idim.eq.2) then
                  xi1=vold(1,iexpnode)
                  xi2=1.d0+vold(2,iexpnode)
                  xi3=1.d0+vold(3,iexpnode)
                  dco=dcos(xi1)
                  dsi=dsin(xi1)
                  dco2=dcos(2.d0*xi1)
                  dsi2=dsin(2.d0*xi1)
                  dd=xi2**2-xi3**2
               endif
            endif
!
!           rotation node
!
            index=nodempc(3,index)
            irotnode=nodempc(1,index)
!
            if(newknot.eq.1) then
               w(1)=vold(1,irotnode)
               w(2)=vold(2,irotnode)
               w(3)=vold(3,irotnode)
               ww=dsqrt(w(1)*w(1)+w(2)*w(2)+w(3)*w(3))
!
               c1=dcos(ww)
               if(ww.gt.1.d-10) then
                  c2=dsin(ww)/ww
               else
                  c2=1.d0
               endif
               if(ww.gt.1.d-5) then
                  c3=(1.d0-c1)/ww**2
               else
                  c3=0.5d0
               endif
!
!              rotation matrix c
!
               do i=1,3
                  do j=1,3
                     c(i,j)=c1*d(i,j)+
     &                  c2*(e(i,1,j)*w(1)+e(i,2,j)*w(2)+e(i,3,j)*w(3))+
     &                  c3*w(i)*w(j)
                  enddo
               enddo
!
               c4=-c2
               if(ww.gt.0.00464159) then
                  c5=(ww*dcos(ww)-dsin(ww))/ww**3
               else
                  c5=-1.d0/3.d0
               endif
               if(ww.gt.0.0031623) then
                  c6=(ww*dsin(ww)-2.d0+2.d0*dcos(ww))/ww**4
               else
                  c6=-1.d0/12.d0
               endif
!
!              derivative of the rotation matrix c with respect to
!              the rotation vector w
!
               do i=1,3
                  do j=1,3
                     do k=1,3
                        dc(i,j,k)=c4*w(k)*d(i,j)+
     &                            c5*w(k)*(e(i,1,j)*w(1)+
     &                                     e(i,2,j)*w(2)+e(i,3,j)*w(3))+
     &                            c2*e(i,k,j)+
     &                            c6*w(k)*w(i)*w(j)+
     &                            c3*(d(i,k)*w(j)+d(j,k)*w(i))
                     enddo
                  enddo
               enddo
!
               if((idim.eq.1).or.(idim.eq.3)) then
!
!                 dummy variable for constant term
!
                  do i=1,3
                     do j=1,3
                        f(i,j)=expan*c(i,j)-d(i,j)
                     enddo
                  enddo
!
!                 derivative of the rotation matrix w.r.t the rotation
!                 vector multiplied by the expansion coefficient
!
                  do i=1,3
                     do j=1,3
                        do k=1,3
                           dc(i,j,k)=dc(i,j,k)*expan
                        enddo
                     enddo
                  enddo
               elseif(idim.eq.2) then
!
!                 local unit vectors
!
                  do i=1,3
                     t1(i)=co(i,irotnode)
                     e1(i)=co(i,iexpnode)
                  enddo
                  e2(1)=t1(2)*e1(3)-t1(3)*e1(2)
                  e2(2)=t1(3)*e1(1)-t1(1)*e1(3)
                  e2(3)=t1(1)*e1(2)-t1(2)*e1(1)
!
                  do i=1,3
                     do j=1,3
                        e11(i,j)=e1(i)*e1(j)
                        e22(i,j)=e2(i)*e2(j)
                        e12(i,j)=e1(i)*e2(j)+e2(i)*e1(j)
!
                        u(i,j)=t1(i)*t1(j)
     &                        +((xi2*dco)**2+(xi3*dsi)**2)*e11(i,j)
     &                        +((xi2*dsi)**2+(xi3*dco)**2)*e22(i,j)
     &                        +dd*dco*dsi*e12(i,j)
                        u1(i,j)=dd*(dco2*e12(i,j)
     &                              -dsi2*(e11(i,j)-e22(i,j)))
                        u2(i,j)=2.d0*xi2*(dco*dco*e11(i,j)
     &                         +dsi*dsi*e22(i,j))+xi2*dsi2*e12(i,j)
                        u3(i,j)=2.d0*xi3*(dsi*dsi*e11(i,j)
     &                         +dco*dco*e22(i,j))-xi3*dsi2*e12(i,j)
                     enddo
                  enddo
!
!                 calculating r.u, r.u1, r.u2 and r.u3
!
                  do i=1,3
                     do j=1,3
                        ru(i,j)=0.d0
                        ru1(i,j)=0.d0
                        ru2(i,j)=0.d0
                        ru3(i,j)=0.d0
!
                        do k=1,3
                           ru(i,j)=ru(i,j)+c(i,k)*u(k,j)
                           ru1(i,j)=ru1(i,j)+c(i,k)*u1(k,j)
                           ru2(i,j)=ru2(i,j)+c(i,k)*u2(k,j)
                           ru3(i,j)=ru3(i,j)+c(i,k)*u3(k,j)
                        enddo
                        f(i,j)=ru(i,j)-d(i,j)
                     enddo
                  enddo
!
!                 calculating dc.u
!
                  do i=1,3
                     do j=1,3
                        do k=1,3
                           dcu(i,j,k)=0.d0
                           do l=1,3
                              dcu(i,j,k)=dcu(i,j,k)+dc(i,l,k)*u(l,j)
                           enddo
                           dc(i,j,k)=dcu(i,j,k)
                        enddo
                     enddo
                  enddo
               endif
!
            endif
!
!           determining the coefficients of the expansion degrees 
!           of freedom
!
            if((idim.eq.1).or.(idim.eq.3)) then
               coefmpc(indexexp)=c(idir,1)*(co(1,irefnode)-co(1,inode))+
     &           c(idir,2)*(co(2,irefnode)-co(2,inode))+
     &           c(idir,3)*(co(3,irefnode)-co(3,inode))
            elseif(idim.eq.2) then
!
!              if xi2=xi3 xi1 cannot be determined, since its coefficient
!              is always zero (cf. definition of u1)
!
               if(dabs(xi2-xi3).lt.1.d-10) then
                  coefmpc(indexexp1)=0.d0
c                  if(nodempc(2,indexexp1).ne.2) then
                     if(icascade.lt.1) icascade=1
                     nodempc(2,indexexp1)=2
c                  endif
               else
                  coefmpc(indexexp1)=ru1(idir,1)*
     &                      (co(1,irefnode)-co(1,inode))+
     &                 ru1(idir,2)*(co(2,irefnode)-co(2,inode))+
     &                 ru1(idir,3)*(co(3,irefnode)-co(3,inode))
c                  if(nodempc(2,indexexp1).ne.1) then
                     if(icascade.lt.1) icascade=1
                     nodempc(2,indexexp1)=1
c                  endif
               endif
               coefmpc(indexexp2)=ru2(idir,1)*
     &                      (co(1,irefnode)-co(1,inode))+
     &           ru2(idir,2)*(co(2,irefnode)-co(2,inode))+
     &           ru2(idir,3)*(co(3,irefnode)-co(3,inode))
               coefmpc(indexexp3)=ru3(idir,1)*
     &                      (co(1,irefnode)-co(1,inode))+
     &           ru3(idir,2)*(co(2,irefnode)-co(2,inode))+
     &           ru3(idir,3)*(co(3,irefnode)-co(3,inode))
            endif
!
!           determining the coefficients of the rotational degrees
!           of freedom
!
            coefmpc(index)=(dc(idir,1,1)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,1)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,1)*(co(3,irefnode)-co(3,inode)))
!
            index=nodempc(3,index)
            coefmpc(index)=(dc(idir,1,2)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,2)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,2)*(co(3,irefnode)-co(3,inode)))
!
            index=nodempc(3,index)
            coefmpc(index)=(dc(idir,1,3)*(co(1,irefnode)-co(1,inode))+
     &           dc(idir,2,3)*(co(2,irefnode)-co(2,inode))+
     &           dc(idir,3,3)*(co(3,irefnode)-co(3,inode)))
!
!           determining the nonhomogeneous part
!
            index=nodempc(3,index)
            coefmpc(index)=1.d0
!
!           old value of the nonhomogeneous term must be zero
!
            vold(nodempc(2,index),nodempc(1,index))=0.d0
            idof=8*(nodempc(1,index)-1)+nodempc(2,index)
            call nident(ikboun,idof,nboun,id)
            xbounact(ilboun(id))=f(idir,1)*(co(1,irefnode)-co(1,inode))+
     &           f(idir,2)*(co(2,irefnode)-co(2,inode))+
     &           f(idir,3)*(co(3,irefnode)-co(3,inode))-
     &           vold(idir,irefnode)+vold(idir,inode)
!
         elseif(labmpc(ii)(1:8).eq.'STRAIGHT') then
!
!           determining nodes and directions involved in MPC
!            
            index=ipompc(ii)
            np=nodempc(1,index)
            j=nodempc(2,index)
            index=nodempc(3,index)
            i=nodempc(2,index)
            index=nodempc(3,index)
            na=nodempc(1,index)
            index=nodempc(3,nodempc(3,index))
            nb=nodempc(1,index)
!
!           determining the coefficients
!
            index=ipompc(ii)
            c2=co(i,na)+vold(i,na)-co(i,nb)-vold(i,nb)
            if(dabs(c2).lt.1.d-5) then
               write(*,*) '*WARNING in nonlinmpc: coefficient of'
               write(*,*)
     &     '          dependent node in STRAIGHT MPC is zero'
               idofrem=8*(np-1)+j
!
!              determining a new dependent term       
!               
               ddmax=abs(c2)
               l=i
               m=j
               do k=1,2
                  l=l+1
                  m=m+1
                  if(l.gt.3) l=l-3
                  if(m.gt.3) m=m-3
                  dd=dabs(co(l,na)+vold(l,na)-co(l,nb)-vold(l,nb))
                  if(dd.gt.ddmax) then
                     ddmax=dd
                     lmax=l
                     mmax=m
                  endif
               enddo
               i=lmax
               j=mmax
               idofins=8*(np-1)+j
!
               call changedepterm(ikmpc,ilmpc,nmpc,ii,idofrem,idofins)
!
               index=ipompc(ii)
               nodempc(2,index)=j
               index=nodempc(3,index)
               nodempc(2,index)=i
               index=nodempc(3,index)
               nodempc(2,index)=j
               index=nodempc(3,index)
               nodempc(2,index)=i
               index=nodempc(3,index)
               nodempc(2,index)=j
               index=nodempc(3,index)
               nodempc(2,index)=i
               index=nodempc(3,index)
               nodempc(2,index)=j
               if(icascade.eq.0) icascade=1
               c2=co(i,na)+vold(i,na)-co(i,nb)-vold(i,nb)
            endif
            coefmpc(index)=c2
            index=nodempc(3,index)
            c3=co(j,nb)+vold(j,nb)-co(j,na)-vold(j,na)
            coefmpc(index)=c3
            index=nodempc(3,index)
            c5=co(i,nb)+vold(i,nb)-co(i,np)-vold(i,np)
            coefmpc(index)=c5
            index=nodempc(3,index)
            c6=co(j,np)+vold(j,np)-co(j,nb)-vold(j,nb)
            coefmpc(index)=c6
            index=nodempc(3,index)
            c4=co(i,np)+vold(i,np)-co(i,na)-vold(i,na)
            coefmpc(index)=c4
            index=nodempc(3,index)
            c1=co(j,na)+vold(j,na)-co(j,np)-vold(j,np)
            coefmpc(index)=c1
            index=nodempc(3,index)
!
!           nonhomogeneous term
!
            coefmpc(index)=1.d0
!
!           old value of the nonhomogeneous term must be zero
!
            idof=8*(nodempc(1,index)-1)+nodempc(2,index)
            call nident(ikboun,idof,nboun,id)
            xbounact(ilboun(id))=-c1*c2+c3*c4
            if(newstep.eq.1) xboun(ilboun(id))=xbounact(ilboun(id))
            vold(nodempc(2,index),nodempc(1,index))=
     &          (1.d0-reltime)*xboun(ilboun(id))
         elseif(labmpc(ii)(1:5).eq.'PLANE') then
!
!           determining nodes and directions involved in MPC
!            
            index=ipompc(ii)
            np=nodempc(1,index)
            i1=nodempc(2,index)
            index=nodempc(3,index)
            i2=nodempc(2,index)
            index=nodempc(3,index)
            i3=nodempc(2,index)
            index=nodempc(3,index)
            na=nodempc(1,index)
            index=nodempc(3,nodempc(3,nodempc(3,index)))
            nb=nodempc(1,index)
            index=nodempc(3,nodempc(3,nodempc(3,index)))
            nc=nodempc(1,index)
!
!           determining the coefficients
!
            a11=co(i1,np)+vold(i1,np)-co(i1,nc)-vold(i1,nc)
            a12=co(i2,np)+vold(i2,np)-co(i2,nc)-vold(i2,nc)
            a13=co(i3,np)+vold(i3,np)-co(i3,nc)-vold(i3,nc)
            a21=co(i1,na)+vold(i1,na)-co(i1,nc)-vold(i1,nc)
            a22=co(i2,na)+vold(i2,na)-co(i2,nc)-vold(i2,nc)
            a23=co(i3,na)+vold(i3,na)-co(i3,nc)-vold(i3,nc)
            a31=co(i1,nb)+vold(i1,nb)-co(i1,nc)-vold(i1,nc)
            a32=co(i2,nb)+vold(i2,nb)-co(i2,nc)-vold(i2,nc)
            a33=co(i3,nb)+vold(i3,nb)-co(i3,nc)-vold(i3,nc)
!
            b11=a22*a33-a23*a32
            b12=a31*a23-a21*a33
            b13=a21*a32-a31*a22
            b21=a32*a13-a12*a33
            b22=a11*a33-a31*a13
            b23=a31*a12-a11*a32
            b31=a12*a23-a22*a13
            b32=a21*a13-a11*a23
            b33=a11*a22-a12*a21
!
            index=ipompc(ii)
            if(dabs(b11).lt.1.d-5) then
               write(*,*) '*WARNING in nonlinmpc: coefficient of'
               write(*,*) '         dependent node in PLANE MPC is zero'
!
               idofrem=8*(nodempc(1,index)-1)+i1
!
               if(dabs(b12).gt.dabs(b13)) then
                  idofins=8*(nodempc(1,index)-1)+i2
                  call changedepterm
     &                   (ikmpc,ilmpc,nmpc,ii,idofrem,idofins)
                  coefmpc(index)=b12
                  nodempc(2,index)=i2
                  index=nodempc(3,index)
                  coefmpc(index)=b11
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b13
                  index=nodempc(3,index)
                  coefmpc(index)=b22
                  nodempc(2,index)=i2
                  index=nodempc(3,index)
                  coefmpc(index)=b21
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b23
                  index=nodempc(3,index)
                  coefmpc(index)=b32
                  nodempc(2,index)=i2
                  index=nodempc(3,index)
                  coefmpc(index)=b31
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b33
                  index=nodempc(3,index)
                  coefmpc(index)=-b12-b22-b32
                  nodempc(2,index)=i2
                  index=nodempc(3,index)
                  coefmpc(index)=-b11-b21-b31
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=-b13-b23-b33
                  if(icascade.eq.0) icascade=1
               else
                  idofins=8*(nodempc(1,index)-1)+i3
                  call changedepterm
     &                  (ikmpc,ilmpc,nmpc,ii,idofrem,idofins)
                  coefmpc(index)=b13
                  nodempc(2,index)=i3
                  index=nodempc(3,index)
                  coefmpc(index)=b12
                  index=nodempc(3,index)
                  coefmpc(index)=b11
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b23
                  nodempc(2,index)=i3
                  index=nodempc(3,index)
                  coefmpc(index)=b22
                  index=nodempc(3,index)
                  coefmpc(index)=b21
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=b33
                  nodempc(2,index)=i3
                  index=nodempc(3,index)
                  coefmpc(index)=b32
                  index=nodempc(3,index)
                  coefmpc(index)=b31
                  nodempc(2,index)=i1
                  index=nodempc(3,index)
                  coefmpc(index)=-b13-b23-b33
                  nodempc(2,index)=i3
                  index=nodempc(3,index)
                  coefmpc(index)=-b12-b22-b32
                  index=nodempc(3,index)
                  coefmpc(index)=-b11-b21-b31
                  nodempc(2,index)=i1
                  if(icascade.eq.0) icascade=1
               endif
            else
               coefmpc(index)=b11
               index=nodempc(3,index)
               coefmpc(index)=b12
               index=nodempc(3,index)
               coefmpc(index)=b13
               index=nodempc(3,index)
               coefmpc(index)=b21
               index=nodempc(3,index)
               coefmpc(index)=b22
               index=nodempc(3,index)
               coefmpc(index)=b23
               index=nodempc(3,index)
               coefmpc(index)=b31
               index=nodempc(3,index)
               coefmpc(index)=b32
               index=nodempc(3,index)
               coefmpc(index)=b33
               index=nodempc(3,index)
               coefmpc(index)=-b11-b21-b31
               index=nodempc(3,index)
               coefmpc(index)=-b12-b22-b32
               index=nodempc(3,index)
               coefmpc(index)=-b13-b23-b33
            endif
            index=nodempc(3,index)
            coefmpc(index)=1.d0
            idof=8*(nodempc(1,index)-1)+nodempc(2,index)
!
!           old value of the nonhomogeneous term must be zero
!
            call nident(ikboun,idof,nboun,id)
            xbounact(ilboun(id))=a11*b11+a12*b12+a13*b13
            if(newstep.eq.1) xboun(ilboun(id))=xbounact(ilboun(id))
            vold(nodempc(2,index),nodempc(1,index))=0.d0
c         elseif(labmpc(ii)(1:9).eq.'ISOCHORIC') then
c            isochoric=.true.
c!
c!           next segment is deactivated (CYCLID instead of CYCLIC):
c!           cylic MPC's are considered to be linear
c!
c         elseif((labmpc(ii)(1:6).eq.'CYCLID').and.(ithermal.ne.2)) then
c            index=ipompc(ii)
c            noded=nodempc(1,index)
c!
c!           check for thermal MPC
c!
c            if(nodempc(2,index).eq.0) cycle loop
c!
c!           check whether the next two MPC's are cyclic MPC's
c!           applied to the same dependent node
c!            
c            if((nodempc(1,ipompc(ii+1)).ne.noded).or.
c     &         (labmpc(ii+1)(1:6).ne.'CYCLIC').or.
c     &         (nodempc(1,ipompc(ii+2)).ne.noded).or.
c     &         (labmpc(ii+2)(1:6).ne.'CYCLIC')) then
c               write(*,*) '*WARNING in nonlinmpc: no three'
c               write(*,*) '         cyclic MPCs pertaining'
c               write(*,*) '         to the same dependent node;'
c               write(*,*) '         no update'
c               cycle loop
c            endif
c!
c!           finding the cyclic symmetry axis
c!
c            do i=1,ntrans
c               if(trab(7,i).eq.2) exit
c            enddo
c            if(i.gt.ntrans) then
c               write(*,*) '*ERROR in nonlinmpc: cyclic symmetry'
c               write(*,*) '       axis not found'
c               call exit(201)
c            endif
c            do j=1,6
c               csab(j)=trab(j,i)
c            enddo
c            csab(7)=-1
c!
c!           determining the independent node
c!                     
c            nodei=0
c            do
c               if(nodempc(1,index).ne.noded) then
c                  if(nodei.eq.0) then
c                     nodei=nodempc(1,index)
c                  elseif(nodei.ne.nodempc(1,index)) then
c                     write(*,*) '*WARNING in nonlinmpc:'
c                     write(*,*) '          cyclic symmetry conditions'
c                     write(*,*) '          between unequal meshes'
c                     write(*,*) '          no update'
c                     cycle loop
c                  endif
c               endif
c               index=nodempc(3,index)
c               if(index.eq.0) then
c                  if(nodei.eq.0) then
c                     write(*,*) '*ERROR in nonlinmpc:'
c                     write(*,*) '       no independent node found'
c                     call exit(201)
c                  else
c                     exit
c                  endif
c               endif
c            enddo
c!
c!           actual location of dependent and independent node
c!
c            do i=1,3
c               pd(i)=co(i,noded)+vold(i,noded)
c               pi(i)=co(i,nodei)+vold(i,nodei)
c            enddo
c!
c!           update transformation matrix
c!
c            call transformatrix(csab,pd,ad)
c            call transformatrix(csab,pi,ai)
c!     
c!     checking for latin hypercube positions in matrix al none of
c!     which are zero
c!     
c            do inum=1,6
c               if((dabs(ad(lathyp(1,inum),1)).gt.1.d-3).and.
c     &            (dabs(ad(lathyp(2,inum),2)).gt.1.d-3).and.
c     &            (dabs(ad(lathyp(3,inum),3)).gt.1.d-3)) exit
c            enddo
c!
c!           remove old DOFs
c!
c            do j=1,3
c               idof=8*(noded-1)+j
c               call nident(ikmpc,idof,nmpc,id)
c               if(id.lt.0) then
c                  write(*,*) '*ERROR in nonlinmpc: error in'
c                  write(*,*) '       MPC database'
c                  call exit(201)
c               elseif(ikmpc(id).ne.idof) then
c                  write(*,*) '*ERROR in nonlinmpc: error in'
c                  write(*,*) '       MPC database'
c                  call exit(201)
c               endif
c!
c               do k=id,nmpc-1
c                  ikmpc(k)=ikmpc(k+1)
c                  ilmpc(k)=ilmpc(k+1)
c               enddo
c            enddo
c!
c!           add new MPCs
c!
c            ii=ii-1
c            do ndir=1,3
c               ii=ii+1
c               number=lathyp(ndir,inum)
c               idof=8*(noded-1)+number
c               call nident(ikmpc,idof,nmpc-1,id)
c               if(id.gt.0) then
c                  if(ikmpc(id).eq.idof) then
c                     write(*,*) '*WARNING in nonlinmpc: cyclic MPC
c     & in node'
c                     write(*,*) '         ',noded,' and direction ',ndir
c                     write(*,*) '         cannot be created: the'
c                     write(*,*) '         DOF in this node is already us
c     &ed'
c                     cycle
c                  endif
c               endif
c               number=number-1
c!     
c!     updating ikmpc and ilmpc
c!     
c               do j=nmpc,id+2,-1
c                  ikmpc(j)=ikmpc(j-1)
c                  ilmpc(j)=ilmpc(j-1)
c               enddo
c               ikmpc(id+1)=idof
c               ilmpc(id+1)=nmpc
c!
c!              update the MPC coefficients
c!
c               index=ipompc(ii)
c               do j=1,3
c                  number=number+1
c                  if(number.gt.3) number=1
c                  if(dabs(ad(number,ndir)).lt.1.d-5) cycle
c                  if(index.eq.0) then
c                     write(*,*)'*ERROR in nonlinmpc: index=0'
c                     call exit(201)
c                  endif
c                  nodempc(1,index)=noded
c                  nodempc(2,index)=number
c                  coefmpc(index)=ad(number,ndir)
c                  index=nodempc(3,index)
c               enddo
c               do j=1,3
c                  number=number+1
c                  if(number.gt.3) number=1
c                  if(dabs(ai(number,ndir)).lt.1.d-5) cycle
c                  if(index.eq.0) then
c                     write(*,*)'*ERROR in nonlinmpc: index=0'
c                     call exit(201)
c                  endif
c                  nodempc(1,index)=nodei
c                  nodempc(2,index)=number
c                  coefmpc(index)=-ai(number,ndir)
c                  index=nodempc(3,index)
c               enddo
c            enddo
         elseif((labmpc(ii)(1:20).ne.'                    ').and.
     &          (labmpc(ii)(1:10).ne.'PRETENSION').and.
     &          (labmpc(ii)(1:7).ne.'CONTACT').and.
     &          (labmpc(ii)(1:7).ne.'NETWORK').and.
     &          (labmpc(ii)(1:5).ne.'FLUID').and.
     &          (labmpc(ii)(1:6).ne.'CYCLIC').and.
     &          (labmpc(ii)(1:9).ne.'SUBCYCLIC')) then
            index=ipompc(ii)
            i=0
            do
               if(index.eq.0) exit
               node=nodempc(1,index)
               i=i+1
               iaux(i)=nodempc(2,index)
               iaux(maxlenmpc+i)=node
               aux(6*maxlenmpc+i)=coefmpc(index)
               do j=1,3
                  aux(3*(i-1)+j)=co(j,node)
                  aux(3*(maxlenmpc+i-1)+j)=vold(j,node)
               enddo
               index=nodempc(3,index)
            enddo
            n=i-1
            if((labmpc(ii)(1:7).eq.'MEANROT').or.
     &         (labmpc(ii)(1:1).eq.'1')) then
               call umpc_mean_rot(aux,aux(3*maxlenmpc+1),const,
     &            aux(6*maxlenmpc+1),iaux,n,fmpc(ii),iit,idiscon,
     &            iaux(maxlenmpc+1),ikmpc,nmpc,ikboun,nboun,
     &            labmpc(ii))
            elseif(labmpc(ii)(1:4).eq.'DIST') then
               call umpc_dist(aux,aux(3*maxlenmpc+1),const,
     &            aux(6*maxlenmpc+1),iaux,n,fmpc(ii),iit,idiscon)
            elseif(labmpc(ii)(1:3).eq.'GAP') then
               call umpc_gap(aux,aux(3*maxlenmpc+1),const,
     &            aux(6*maxlenmpc+1),iaux,n,fmpc(ii),iit,idiscon)
            elseif(labmpc(ii)(1:4).eq.'USER') then
               call umpc_user(aux,aux(3*maxlenmpc+1),const,
     &            aux(6*maxlenmpc+1),iaux,n,fmpc(ii),iit,idiscon)
            else
               write(*,*) '*ERROR in nonlinmpc: mpc of type ',labmpc(ii)
               write(*,*) '       is unknown'
               call exit(201)
            endif
            index=ipompc(ii)
!
            if((iaux(1).ne.nodempc(2,index)).or.
     &         (iaux(maxlenmpc+1).ne.nodempc(1,index))) then
!
!              dependent MPC has changed
!  
               idofrem=8*(nodempc(1,index)-1)+nodempc(2,index)
               idofins=8*(nodempc(1,index)-1)+iaux(1)
               call changedepterm(ikmpc,ilmpc,nmpc,ii,idofrem,idofins)
               if(icascade.eq.0) icascade=1
            endif
!
            i=0
            do
               if(index.eq.0) exit
               i=i+1
               if(i.le.n) then
!
!                 check whether any directions have changed:
!                 necessitates calling of remastruct
!
                  if(iaux(i).ne.nodempc(2,index)) then
                     if(icascade.eq.0) icascade=1
                  endif
                  nodempc(2,index)=iaux(i)
                  coefmpc(index)=aux(6*maxlenmpc+i)
               else
                  coefmpc(index)=1.d0
!
!                 old value of the nonhomogeneous term must be zero
!
                  vold(nodempc(2,index),nodempc(1,index))=0.d0
                  idof=8*(nodempc(1,index)-1)+nodempc(2,index)
                  call nident(ikboun,idof,nboun,id)
                  xbounact(ilboun(id))=const
               endif
               index=nodempc(3,index)
            enddo
         endif
      enddo loop
c!
c!     incompressible material
c!
c      if(.not.isochoric) return
c!
c!     initialization of the mpc's
c!
c      nmpc01=0
c      nmpcdif=0
c      do i=1,nmpc
c         if(labmpc(i)(1:9).eq.'ISOCHORIC') then
c            if(nmpc01.eq.0) nmpc01=i
c            nmpcdif=i
c            index=ipompc(i)
c            do
c               if(nodempc(3,index).eq.0) then
c                  idof=8*(nodempc(1,index)-1)+nodempc(2,index)
c                  call nident(ikboun,idof,nboun,id)
c                  xbounact(ilboun(id))=0.d0
c                  exit
c               endif
c               coefmpc(index)=0.d0
c               index=nodempc(3,index)
c            enddo
c         endif
c      enddo
c      nmpc0=nmpc01-1
c      nmpcdif=nmpcdif-nmpc0
c!
c      do i=1,ne
c         if(ipkon(i).lt.0) cycle
c         if(lakon(i)(1:7).eq.'C3D20RI') then
c            indexe=ipkon(i)
c!
c            do j=1,20
c               node=kon(indexe+j)
c               do k=1,3
c                  xlag(k,j)=co(k,node)
c                  xeul(k,j)=xlag(k,j)+vold(k,node)
c               enddo
c            enddo
c!
c            do j=1,8
c               mpc=0
c               node=kon(indexe+j)
c               label(1:9)='ISOCHORIC'
c               write(label(10:20),'(i11)') node
cc               write(*,*) 'nonlinmpclab ',label
c               call cident20(labmpc(nmpc01),label,nmpcdif,id)
c               id=id+nmpc0
cc               write(*,*) 'nonlinmpclab ',id,label,labmpc(id)
c               if(id.gt.0) then
c                  if(labmpc(id).eq.label) then
c                     mpc=id
c                  endif
c               endif
c               if(mpc.eq.0) cycle
c!
c               xi=coloc(1,j)
c               et=coloc(2,j)
c               ze=coloc(3,j)
c!
c               call deuldlag(xi,et,ze,xlag,xeul,xj,a)
c!
c               b(1,1)=a(2,2)*a(3,3)-a(2,3)*a(3,2)
c               b(1,2)=a(3,1)*a(2,3)-a(2,1)*a(3,3)
c               b(1,3)=a(2,1)*a(3,2)-a(3,1)*a(2,2)
c               b(2,1)=a(3,2)*a(1,3)-a(1,2)*a(3,3)
c               b(2,2)=a(1,1)*a(3,3)-a(3,1)*a(1,3)
c               b(2,3)=a(3,1)*a(1,2)-a(1,1)*a(3,2)
c               b(3,1)=a(1,2)*a(2,3)-a(2,2)*a(1,3)
c               b(3,2)=a(2,1)*a(1,3)-a(1,1)*a(2,3)
c               b(3,3)=a(1,1)*a(2,2)-a(1,2)*a(2,1)
c!
c               index=ipompc(mpc)
c               do
c                  if(nodempc(3,index).eq.0) then
c                     coefmpc(index)=1.d0
c                     idof=8*(nodempc(1,index)-1)+nodempc(2,index)
c                     call nident(ikboun,idof,nboun,id)
c                     xbounact(ilboun(id))=xbounact(ilboun(id))+
c     &                    a(1,1)*b(1,1)+a(1,2)*b(1,2)+a(1,3)*b(1,3)
c     &                    -1.d0/xj
cc                     write(*,*) 'nonlinmpcboun ',nodempc(1,index),
cc     &                        nodempc(2,index),ilboun(id),
cc     &                       xbounact(ilboun(id))
c                     exit
c                  else
c                     node=nodempc(1,index)
c                     idir=nodempc(2,index)
c                     do k=1,7
c                        if(kon(indexe+neigh(k,j)).eq.node) then
c                           if(k.eq.1) then
c                              if(idir.eq.1) then
c                                 coefmpc(index)=coefmpc(index)+
c     &                             1.5d0*(xi*b(1,1)+et*b(1,2)+ze*b(1,3))
c                              elseif(idir.eq.2) then
c                                 coefmpc(index)=coefmpc(index)+
c     &                             1.5d0*(xi*b(2,1)+et*b(2,2)+ze*b(2,3))
c                              elseif(idir.eq.3) then
c                                 coefmpc(index)=coefmpc(index)+
c     &                             1.5d0*(xi*b(3,1)+et*b(3,2)+ze*b(3,3))
c                              endif
c                           elseif(k.eq.2) then
c                              if(idir.eq.1) then
c                               coefmpc(index)=coefmpc(index)-
c     &                                        2.d0*xi*b(1,1)
c                              elseif(idir.eq.2) then
c                               coefmpc(index)=coefmpc(index)-
c     &                                        2.d0*xi*b(2,1)
c                              elseif(idir.eq.3) then
c                               coefmpc(index)=coefmpc(index)-
c     &                                        2.d0*xi*b(3,1)
c                              endif
c                           elseif(k.eq.3) then
c                              if(idir.eq.1) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*xi*b(1,1)
c                              elseif(idir.eq.2) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*xi*b(2,1)
c                              elseif(idir.eq.3) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*xi*b(3,1)
c                              endif
c                           elseif(k.eq.4) then
c                              if(idir.eq.1) then
c                               coefmpc(index)=coefmpc(index)-
c     &                                        2.d0*et*b(1,2)
c                              elseif(idir.eq.2) then
c                               coefmpc(index)=coefmpc(index)-
c     &                                        2.d0*et*b(2,2)
c                              elseif(idir.eq.3) then
c                               coefmpc(index)=coefmpc(index)-
c     &                                        2.d0*et*b(3,2)
c                              endif
c                           elseif(k.eq.5) then
c                              if(idir.eq.1) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*et*b(1,2)
c                              elseif(idir.eq.2) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*et*b(2,2)
c                              elseif(idir.eq.3) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*et*b(3,2)
c                              endif
c                           elseif(k.eq.6) then
c                              if(idir.eq.1) then
c                              coefmpc(index)=coefmpc(index)-
c     &                                       2.d0*ze*b(1,3)
c                              elseif(idir.eq.2) then
c                              coefmpc(index)=coefmpc(index)-
c     &                                       2.d0*ze*b(2,3)
c                              elseif(idir.eq.3) then
c                              coefmpc(index)=coefmpc(index)-
c     &                                       2.d0*ze*b(3,3)
c                              endif
c                           elseif(k.eq.7) then
c                              if(idir.eq.1) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*ze*b(1,3)
c                              elseif(idir.eq.2) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*ze*b(2,3)
c                              elseif(idir.eq.3) then
c                              coefmpc(index)=coefmpc(index)+
c     &                                       0.5d0*ze*b(3,3)
c                              endif
c                           endif
c                           exit
c                        endif
c                     enddo
c                  endif
c                  index=nodempc(3,index)
c               enddo
c!
c            enddo
c         endif
c      enddo
!
      return
      end






