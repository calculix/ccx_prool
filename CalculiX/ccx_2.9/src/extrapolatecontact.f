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
      subroutine extrapolatecontact(yi,yn,ipkon,inum,kon,lakon,nfield,
     &  nk,ne,mi,ndim,co,cflag,vold,force,pslavsurf,islavact,islavnode,
     &  nslavnode,ntie,islavsurf,ielprop,prop,ielmat,ne0)
!
!     extrapolates contact values at the integration points to the 
!     nodes (for face-to-face penalty contact)
!
!     the C3D20RB element has 50 integration points, however, the
!     first 8 integration points coincide with those of a C3D20R
!     element. In this routine and in errorestimator.f the C3D20RB
!     element is treated as an ordinary C3D20R element
!
      implicit none
!
      logical force
!
      character*1 cflag
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,indexe,nope,
     &  nonei20(3,12),nfield,nonei10(3,6),nk,i,j,k,l,ndim,
     &  nonei15(3,9),m,iflag,jj,indexc,islavsurf(2,*),ll,
     &  indexcj,ifacew1(4,5),islavnode(*),nslavnode(*),ntie,
     &  nlayer,nopeexp,ifacew2(8,5),ifaceq(9,6),ifacet(7,4),
     &  nopespring,ifaces,nopespringj,ifacej,jfaces,n,nelems,
     &  idgn,idgnr,idglda,idgip(4),idgldb,idginfo,igauss,islavact(*),
     &  konl(26),node,nopes,ielprop(*),ielmat(mi(3),*),ne0,
     &  nodes(8)
!
      real*8 yi(ndim,mi(1),*),yn(nfield,*),field(nfield,20*mi(3)),
     &  co(3,*),xi,et,vold(0:mi(2),*),xs2(3,7),xsj2(3),shp2(7,8),
     &  aa(4,4),bb(4,6),pl(3,8),pslavsurf(3,*),xslavnor(3,nk),
     &  cc(3,4),dd(3,6),c_limit(2,nfield),nodepos(4,2),xn(3),
     &  t1(3),t2(3),trac(3),xquad(2,9),xtri(2,7),xl2s(3,9),
     &  stn(6,nk),dt1,dl,a2(6,2),a4(4,4),a27(20,27),a9(6,9),
     &  a8(8,8),prop(*)
!
      data nonei10 /5,1,2,6,2,3,7,3,1,8,1,4,9,2,4,10,3,4/
!
      data nonei15 /7,1,2,8,2,3,9,3,1,10,4,5,11,5,6,12,6,4,
     &  13,1,4,14,2,5,15,3,6/
!
      data nonei20 /9,1,2,10,2,3,11,3,4,12,4,1,
     &  13,5,6,14,6,7,15,7,8,16,8,5,
     &  17,1,5,18,2,6,19,3,7,20,4,8/
!
!     nodes per face for hex elements
!
      data ifaceq /4,3,2,1,11,10,9,12,21,
     &            5,6,7,8,13,14,15,16,22,
     &            1,2,6,5,9,18,13,17,23,
     &            2,3,7,6,10,19,14,18,24,
     &            3,4,8,7,11,20,15,19,25,
     &            4,1,5,8,12,17,16,20,26/
!
!     nodes per face for tet elements
!
      data ifacet /1,3,2,7,6,5,11,
     &             1,2,4,5,9,8,12,
     &             2,3,4,6,10,9,13,
     &             1,4,3,8,10,7,14/
!
!     nodes per face for linear wedge elements
!
      data ifacew1 /1,3,2,0,
     &             4,5,6,0,
     &             1,2,5,4,
     &             2,3,6,5,
     &             3,1,4,6/
!
!     nodes per face for quadratic wedge elements
!
      data ifacew2 /1,3,2,9,8,7,0,0,
     &             4,5,6,10,11,12,0,0,
     &             1,2,5,4,7,14,10,13,
     &             2,3,6,5,8,15,11,14,
     &             3,1,4,6,9,13,12,15/
!
      data iflag /2/
c!     
c!     local coordinates of the element nodes
c!
c      data xquad /-1.d0,-1.d0,
c     &             1.d0,-1.d0,
c     &             1.d0,1.d0,
c     &            -1.d0,1.d0,
c     &             0.d0,-1.d0,
c     &             1.d0,0.d0,
c     &             0.d0,1.d0,
c     &            -1.d0,0.d0,
c     &             0.d0,0.d0/
c!
c      data xtri /0.d0,0.d0,
c     &           1.d0,0.d0,
c     &           0.d0,1.d0,
c     &           .5d0,0.d0,
c     &           .5d0,.5d0,
c     &           0.d0,.5d0,
c     &           0.333333333333333d0,0.333333333333333d0/
c!
c      data a2 /  1.1455,-0.1455,1.1455,-0.1455,1.1455,-0.1455,
c     &           -0.1455,1.1455,-0.1455,1.1455,-0.1455,1.1455/
c      data a4 /  1.92705, -0.30902, -0.30902, -0.30902,
c     &          -0.30902,  1.92705, -0.30902, -0.30902,
c     &          -0.30902, -0.30902,  1.92705, -0.30902,
c     &          -0.30902, -0.30902, -0.30902,  1.92705/
c      data a9 / 1.63138,-0.32628,-0.32628,-0.52027, 0.10405, 0.10405,
c     &         -0.32628, 1.63138,-0.32628, 0.10405,-0.52027, 0.10405,
c     &         -0.32628,-0.32628, 1.63138, 0.10405, 0.10405,-0.52027,
c     &          0.55556,-0.11111,-0.11111, 0.55556,-0.11111,-0.11111,
c     &         -0.11111, 0.55556,-0.11111,-0.11111, 0.55556,-0.11111,
c     &         -0.11111,-0.11111, 0.55556,-0.11111,-0.11111, 0.55556,
c     &         -0.52027, 0.10405, 0.10405, 1.63138,-0.32628,-0.32628,
c     &          0.10405,-0.52027, 0.10405,-0.32628, 1.63138,-0.32628,
c     &          0.10405, 0.10405,-0.52027,-0.32628,-0.32628, 1.63138/
c!
c!     extrapolation from a 2x2x2=8 integration point scheme in a hex to
c!     the vertex nodes
c!    
c      data a8 /2.549,-.683,.183,-.683,-.683,.183,
c     &        -.04904,.183,-.683,2.549,-.683,.183,
c     &        .183,-.683,.183,-.04904,-.683,.183,
c     &        -.683,2.549,.183,-.04904,.183,-.683,
c     &        .183,-.683,2.549,-.683,-.04904,.183,
c     &        -.683,.183,-.683,.183,-.04904,.183,
c     &        2.549,-.683,.183,-.683,.183,-.683,
c     &        .183,-.04904,-.683,2.549,-.683,.183,
c     &        .183,-.04904,.183,-.683,-.683,.183,
c     &        -.683,2.549,-.04904,.183,-.683,.183,
c     &        .183,-.683,2.549,-.683/  
c!
c!     extrapolation from a 3x3x3=27 integration point scheme in a hex to
c!     the all nodes in a 20-node element
c!    
c      data a27 /
c     &  2.37499,-0.12559,-0.16145,-0.12559,-0.12559,-0.16145, 0.11575,
c     & -0.16145, 0.32628, 0.11111, 0.11111, 0.32628, 0.11111,-0.10405,
c     & -0.10405, 0.11111, 0.32628, 0.11111,-0.10405, 0.11111,-0.31246,
c     & -0.31246, 0.31481, 0.31481, 0.31481, 0.31481,-0.16902,-0.16902,
c     &  1.28439,-0.27072,-0.19444,-0.27072,-0.19444, 0.15961,-0.00661,
c     &  0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.12559, 2.37499,
c     & -0.12559,-0.16145,-0.16145,-0.12559,-0.16145, 0.11575, 0.32628,
c     &  0.32628, 0.11111, 0.11111, 0.11111, 0.11111,-0.10405,-0.10405,
c     &  0.11111, 0.32628, 0.11111,-0.10405,-0.31246, 0.31481, 0.31481,
c     & -0.31246, 0.31481,-0.16902,-0.16902, 0.31481,-0.27072,-0.19444,
c     & -0.27072, 1.28439, 0.15961,-0.00661, 0.15961,-0.19444,-0.27072,
c     &  0.15961, 0.15961,-0.27072,-0.48824,-0.48824,-0.48824,-0.48824,
c     &  0.22898, 0.22898, 0.22898, 0.22898, 0.05556, 0.05556, 0.05556,
c     &  0.05556, 0.05556, 0.05556, 0.05556, 0.05556,-0.22222,-0.22222,
c     & -0.22222,-0.22222, 0.31481,-0.31246,-0.31246, 0.31481,-0.16902,
c     &  0.31481, 0.31481,-0.16902,-0.27072, 1.28439,-0.27072,-0.19444,
c     &  0.15961,-0.19444, 0.15961,-0.00661, 0.15961,-0.27072,-0.27072,
c     &  0.15961,-0.12559,-0.16145,-0.12559, 2.37499,-0.16145, 0.11575,
c     & -0.16145,-0.12559, 0.11111, 0.11111, 0.32628, 0.32628,-0.10405,
c     & -0.10405, 0.11111, 0.11111, 0.11111,-0.10405, 0.11111, 0.32628,
c     &  0.31481, 0.31481,-0.31246,-0.31246,-0.16902,-0.16902, 0.31481,
c     &  0.31481,-0.19444,-0.27072, 1.28439,-0.27072,-0.00661, 0.15961,
c     & -0.19444, 0.15961, 0.15961, 0.15961,-0.27072,-0.27072,-0.16145,
c     & -0.12559, 2.37499,-0.12559, 0.11575,-0.16145,-0.12559,-0.16145,
c     &  0.11111, 0.32628, 0.32628, 0.11111,-0.10405, 0.11111, 0.11111,
c     & -0.10405,-0.10405, 0.11111, 0.32628, 0.11111,-0.31246, 0.31481,
c     & -0.16902, 0.31481,-0.31246, 0.31481,-0.16902, 0.31481,-0.27072,
c     &  0.15961, 0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.27072,
c     &  1.28439,-0.19444,-0.00661,-0.19444,-0.48824,-0.48824, 0.22898,
c     &  0.22898,-0.48824,-0.48824, 0.22898, 0.22898, 0.05556,-0.22222,
c     &  0.05556,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,
c     &  0.05556, 0.05556, 0.05556, 0.31481,-0.31246, 0.31481,-0.16902,
c     &  0.31481,-0.31246, 0.31481,-0.16902,-0.27072,-0.27072, 0.15961,
c     &  0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.19444, 1.28439,
c     & -0.19444,-0.00661,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824,
c     &  0.22898, 0.22898,-0.48824,-0.22222, 0.05556,-0.22222, 0.05556,
c     & -0.22222, 0.05556,-0.22222, 0.05556, 0.05556, 0.05556, 0.05556,
c     &  0.05556,-0.29630,-0.29630,-0.29630,-0.29630,-0.29630,-0.29630,
c     & -0.29630,-0.29630,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,
c     & -0.11111,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,
c     &  0.22898,-0.48824,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824,
c     &  0.22898,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,
c     & -0.22222, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.31481,
c     & -0.16902, 0.31481,-0.31246, 0.31481,-0.16902, 0.31481,-0.31246,
c     &  0.15961, 0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.27072,
c     & -0.27072,-0.19444,-0.00661,-0.19444, 1.28439, 0.22898, 0.22898,
c     & -0.48824,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824, 0.05556,
c     & -0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222,
c     &  0.05556, 0.05556, 0.05556, 0.05556,-0.16902, 0.31481,-0.31246,
c     &  0.31481,-0.16902, 0.31481,-0.31246, 0.31481, 0.15961,-0.27072,
c     & -0.27072, 0.15961, 0.15961,-0.27072,-0.27072, 0.15961,-0.00661,
c     & -0.19444, 1.28439,-0.19444,-0.12559,-0.16145, 0.11575,-0.16145,
c     &  2.37499,-0.12559,-0.16145,-0.12559, 0.11111,-0.10405,-0.10405,
c     &  0.11111, 0.32628, 0.11111, 0.11111, 0.32628, 0.32628, 0.11111,
c     & -0.10405, 0.11111, 0.31481, 0.31481,-0.16902,-0.16902,-0.31246,
c     & -0.31246, 0.31481, 0.31481,-0.19444, 0.15961,-0.00661, 0.15961,
c     &  1.28439,-0.27072,-0.19444,-0.27072,-0.27072,-0.27072, 0.15961,
c     &  0.15961,-0.16145,-0.12559,-0.16145, 0.11575,-0.12559, 2.37499,
c     & -0.12559,-0.16145, 0.11111, 0.11111,-0.10405,-0.10405, 0.32628,
c     &  0.32628, 0.11111, 0.11111, 0.11111, 0.32628, 0.11111,-0.10405,
c     &  0.31481,-0.16902,-0.16902, 0.31481,-0.31246, 0.31481, 0.31481,
c     & -0.31246, 0.15961,-0.00661, 0.15961,-0.19444,-0.27072,-0.19444,
c     & -0.27072, 1.28439,-0.27072, 0.15961, 0.15961,-0.27072, 0.22898,
c     &  0.22898, 0.22898, 0.22898,-0.48824,-0.48824,-0.48824,-0.48824,
c     &  0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556,
c     &  0.05556,-0.22222,-0.22222,-0.22222,-0.22222,-0.16902, 0.31481,
c     &  0.31481,-0.16902, 0.31481,-0.31246,-0.31246, 0.31481, 0.15961,
c     & -0.19444, 0.15961,-0.00661,-0.27072, 1.28439,-0.27072,-0.19444,
c     &  0.15961,-0.27072,-0.27072, 0.15961,-0.16145, 0.11575,-0.16145,
c     & -0.12559,-0.12559,-0.16145,-0.12559, 2.37499,-0.10405,-0.10405,
c     &  0.11111, 0.11111, 0.11111, 0.11111, 0.32628, 0.32628, 0.11111,
c     & -0.10405, 0.11111, 0.32628,-0.16902,-0.16902, 0.31481, 0.31481,
c     &  0.31481, 0.31481,-0.31246,-0.31246,-0.00661, 0.15961,-0.19444,
c     &  0.15961,-0.19444,-0.27072, 1.28439,-0.27072, 0.15961, 0.15961,
c     & -0.27072,-0.27072, 0.11575,-0.16145,-0.12559,-0.16145,-0.16145,
c     & -0.12559, 2.37499,-0.12559,-0.10405, 0.11111, 0.11111,-0.10405,
c     &  0.11111, 0.32628, 0.32628, 0.11111,-0.10405, 0.11111, 0.32628,
c     &  0.11111/
!
      if(nfield.eq.0) return
!
      do i=1,nk
         inum(i)=0
      enddo
!
      do i=1,nk
         do j=1,nfield
            yn(j,i)=0.d0
         enddo
      enddo
!
      i=0
      do
         i=i+1
         if(i.gt.ne) exit
!     
         if(ipkon(i).lt.0) cycle
         indexc=ipkon(i)
!     
         if((lakon(i)(1:1).ne.'E').or.(lakon(i)(7:7).ne.'C')) cycle
         nopespring=kon(indexc)
         ifaces=islavsurf(1,kon(indexc+nopespring+2))
!     
         nelems=int(ifaces/10.d0)
         lakonl=lakon(nelems)
!     
!        determining all contact elements belonging to slave face
!        "ifaces"
!
         n=i
!     
         do
            n=n+1
            if(n.gt.ne) exit
            indexcj=ipkon(n)
            nopespringj=kon(indexcj)
            ifacej=islavsurf(1,kon(indexcj+nopespringj+2))
            if(ifaces.ne.ifacej) exit
         enddo
         n=n-1
         jfaces=ifaces-10*int(ifaces/10.d0)
         indexe=ipkon(nelems)
!     
         if(lakonl(4:4).eq.'2') then
            nope=20
         elseif(lakonl(4:4).eq.'8') then
            nope=8
         elseif((lakonl(4:5).eq.'10').or.(lakonl(4:5).eq.'14')) then
            nope=10
         elseif(lakonl(4:4).eq.'4') then
            nope=4
         elseif(lakonl(4:5).eq.'15') then
            nope=15
         elseif(lakonl(4:4).eq.'6') then
            nope=6
         elseif((lakon(i)(1:1).eq.'E').and.(lakon(i)(7:7).eq.'A')) then
            inum(kon(indexe+1))=inum(kon(indexe+1))+1
            inum(kon(indexe+2))=inum(kon(indexe+2))+1
            cycle
         else
            cycle
         endif
!     
         if((lakonl(4:5).eq.'20').or.(lakonl(4:4).eq.'8').or.
     &        (((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')).and.
     &        (jfaces.gt.2))) then
            if(lakonl(7:8).ne.'LC') then
               field(1:nfield,1:20)=0.d0
               do k=1,4
                  do l=1,4
                     aa(k,l)=0.d0
                  enddo
                  do l=1,nfield
                     bb(k,l)=0.d0
                  enddo
               enddo
               nodepos(1,1)=-1.d0
               nodepos(1,2)=-1.d0
               nodepos(2,1)=1.d0
               nodepos(2,2)=-1.d0
               nodepos(3,1)=1.d0
               nodepos(3,2)=1.d0
               nodepos(4,1)=-1.d0
               nodepos(4,2)=1.d0
               do j=i,n
                  do k=1,nfield
                     if(j.eq.i) then
                        c_limit(1,k)=yi(k,1,j)
                        c_limit(2,k)=yi(k,1,j)
                     endif
                     if(c_limit(1,k).lt.yi(k,1,j)) then
                        c_limit(1,k)=yi(k,1,j)
                     endif
                     if(c_limit(2,k).gt.yi(k,1,j)) then
                        c_limit(2,k)=yi(k,1,j)
                     endif
                  enddo
                  indexcj=ipkon(j)
                  nopespringj=kon(indexcj)
                  igauss=kon(indexcj+nopespringj+1)
                  xi=pslavsurf(1,igauss)
                  et=pslavsurf(2,igauss)
                  if((n-i).lt.2) exit
                  if((n-i).eq.2) then
                     aa(j+1-i,1)=xi
                     aa(j+1-i,2)=et
                     do k=1,nfield
                        dd(j+1-i,k)=yi(k,1,j)
                     enddo
                  elseif((n-i).eq.3) then
                     call shape4q(xi,et,pl,xsj2,xs2,shp2,iflag)
                     do k=1,4
                        aa(j+1-i,k)=shp2(4,k)
                     enddo
                     do k=1,nfield
                        bb(j+1-i,k)=yi(k,1,j)
                     enddo
                  else
                     call shape4q(xi,et,pl,xsj2,xs2,shp2,iflag)
                     do k=1,4
                        do l=k,4
                           aa(k,l)=aa(k,l)+shp2(4,k)*shp2(4,l)
                        enddo
                     enddo
                     do k=1,4
                        do l=1,nfield
                           bb(k,l)=bb(k,l)+shp2(4,k)*yi(l,1,j)
                        enddo
                     enddo
                  endif
               enddo
               if(n.eq.i) then
                  do k=1,4
                     do l=1,nfield
                        bb(k,l)=bb(k,l)+yi(l,1,n)
                     enddo
                  enddo   
               elseif((n-i).eq.1) then
                  do j=i,n
                     do k=1,4
                        do l=1,nfield
                           bb(k,l)=bb(k,l)+yi(l,1,j)*0.5d0
                        enddo
                     enddo
                  enddo
               elseif((n-i).eq.2) then
                  do k=1,4
                     do l=1,nfield
                        call plane_eq(aa(1,1),aa(1,2),dd(1,l),
     &                       aa(2,1),aa(2,2),dd(2,l),aa(3,1),
     &                       aa(3,2),dd(3,l),nodepos(k,1),
     &                       nodepos(k,2),bb(k,l))
                     enddo
                  enddo
               else
                  if((n-i).ne.3) then
                     do k=1,4
                        do l=1,k-1
                           aa(k,l)=aa(l,k)
                        enddo
                     enddo
                  endif
                  idgn=4
                  idgnr=nfield
                  idglda=4
                  idgldb=4
                  call dgesv(idgn,idgnr,aa,idglda,idgip,bb,idgldb,
     &                         idginfo)
               endif
               if((lakonl(4:4).eq.'6').or.(lakonl(4:5).eq.'15')) then
                  do j=1,4
                     do k=1,nfield
                        if((c_limit(1,k).gt.bb(j,k)).and.
     &                     (c_limit(2,k).lt.bb(j,k))) then
                           field(k,ifacew1(j,jfaces))=bb(j,k)
                        endif
                        if(c_limit(1,k).lt.bb(j,k)) then
                           field(k,ifacew1(j,jfaces))=c_limit(1,k)
                        endif
                        if(c_limit(2,k).gt.bb(j,k)) then
                           field(k,ifacew1(j,jfaces))=c_limit(2,k)
                        endif   
                     enddo
                  enddo
               else
                  do j=1,4
                     do k=1,nfield
                        if((c_limit(1,k).gt.bb(j,k)).and.
     &                     (c_limit(2,k).lt.bb(j,k))) then
                           field(k,ifaceq(j,jfaces))=bb(j,k)
                        endif
                        if(c_limit(1,k).lt.bb(j,k)) then
                           field(k,ifaceq(j,jfaces))=c_limit(1,k)
                        endif
                        if(c_limit(2,k).gt.bb(j,k)) then
                           field(k,ifaceq(j,jfaces))=c_limit(2,k)
                        endif   
                     enddo
                  enddo
               endif
            endif
         elseif((lakonl(4:5).eq.'10').or.(lakonl(4:4).eq.'4').or.
     &           (((lakonl(4:5).eq.'15').or.(lakonl(4:4).eq.'6')).and.
     &           (jfaces.le.2))) then
            field(1:nfield,1:15)=0.d0
            if(lakonl(7:8).ne.'LC') then
               do k=1,3
                  do l=1,3
                     cc(k,l)=0.d0
                  enddo
                  do l=1,nfield
                     dd(k,l)=0.d0
                  enddo
               enddo
               do j=i,n
                  do k=1,nfield
                     if(j.eq.i) then
                        c_limit(1,k)=yi(k,1,j)
                        c_limit(2,k)=yi(k,1,j)
                     endif
                     if(c_limit(1,k).lt.yi(k,1,j)) then
                        c_limit(1,k)=yi(k,1,j)
                     endif
                     if(c_limit(2,k).gt.yi(k,1,j)) then
                        c_limit(2,k)=yi(k,1,j)
                     endif
                  enddo
                  indexcj=ipkon(j)
                  nopespringj=kon(indexcj)
                  igauss=kon(indexcj+nopespringj+1)
                  xi=pslavsurf(1,igauss)
                  et=pslavsurf(2,igauss)
                  if((n-i).lt.2) exit
                  if((n-i).eq.2) then
                     call shape3tri(xi,et,pl,xsj2,xs2,shp2,iflag)
                     do k=1,3
                        cc(j+1-i,k)=shp2(4,k)
                     enddo
                     do k=1,nfield
                        dd(j+1-i,k)=yi(k,1,j)
                     enddo
                  else
                     call shape3tri(xi,et,pl,xsj2,xs2,shp2,iflag)
                     do k=1,3
                        do l=1,3
                           cc(k,l)=cc(k,l)+shp2(4,k)*shp2(4,l)
                        enddo
                     enddo
                     do k=1,3
                        do l=1,nfield
                           dd(k,l)=dd(k,l)+shp2(4,k)*yi(l,1,j)
                        enddo
                     enddo
                  endif
               enddo
               if(n.eq.i) then
                  do k=1,3
                     do l=1,nfield
                        dd(k,l)=dd(k,l)+yi(l,1,j)
                     enddo
                  enddo   
               elseif((n-i).eq.1) then
                  do j=i,n
                     do k=1,3
                        do l=1,nfield
                           dd(k,l)=dd(k,l)+yi(l,1,j)*0.5d0
                        enddo
                     enddo
                  enddo
               else
                  idgn=3
                  idgnr=nfield
                  idglda=3
                  idgldb=3
                  call dgesv(idgn,idgnr,cc,idglda,idgip,dd,idgldb,
     &                         idginfo)
               endif
               if((lakonl(4:4).eq.'6').or.(lakonl(4:5).eq.'15')) then
                  do j=1,3
                     do k=1,nfield
                        if((c_limit(1,k).gt.dd(j,k)).and.
     &                     (c_limit(2,k).lt.dd(j,k))) then
                           field(k,ifacew1(j,jfaces))=dd(j,k)
                        endif
                        if(c_limit(1,k).lt.dd(j,k)) then
                           field(k,ifacew1(j,jfaces))=c_limit(1,k)
                        endif
                        if(c_limit(2,k).gt.dd(j,k)) then
                           field(k,ifacew1(j,jfaces))=c_limit(2,k)
                        endif   
                     enddo
                  enddo
               else
                  do j=1,3
                     do k=1,nfield
                        if((c_limit(1,k).gt.dd(j,k)).and.
     &                       (c_limit(2,k).lt.dd(j,k))) then
                           field(k,ifacet(j,jfaces))=dd(j,k)
                        endif
                        if(c_limit(1,k).lt.dd(j,k)) then
                           field(k,ifacet(j,jfaces))=c_limit(1,k)
                        endif
                        if(c_limit(2,k).gt.dd(j,k)) then
                           field(k,ifacet(j,jfaces))=c_limit(2,k)
                        endif   
                     enddo
                  enddo
               endif
            endif  
         endif
!     
!     determining the field values in the midside nodes of the
!     slave face
!     
         if((lakonl(4:5).eq.'20').or.(lakonl(4:6).eq.'26R')) then
            if(lakonl(7:8).ne.'LC') then
               do j=5,8
                  do k=1,nfield
                     field(k,ifaceq(j,jfaces))=
     &                   (field(k,nonei20(2,ifaceq(j,jfaces)-8))+
     &                    field(k,nonei20(3,ifaceq(j,jfaces)-8)))/2.d0
                  enddo
               enddo
            else
               do m=1,nlayer
                  jj=20*(m-1)
                  do j=9,20
                     do k=1,nfield
                        field(k,jj+j)=(field(k,jj+nonei20(2,j-8))
     &                       +field(k,jj+nonei20(3,j-8)))/2.d0
                     enddo
                  enddo
               enddo
            endif
         elseif((lakonl(4:5).eq.'10').or.(lakonl(4:5).eq.'14')) then
            do j=4,6
               do k=1,nfield
                  field(k,ifacet(j,jfaces))=
     &                (field(k,nonei10(2,ifacet(j,jfaces)-4))+
     &                 field(k,nonei10(3,ifacet(j,jfaces)-4)))/2.d0
               enddo
            enddo
         elseif(lakonl(4:5).eq.'15') then
            if(jfaces.le.2) then
               do j=4,6
                  do k=1,nfield
                     field(k,ifacew2(j,jfaces))=
     &                   (field(k,nonei15(2,ifacew2(j,jfaces)-6))+
     &                    field(k,nonei15(3,ifacew2(j,jfaces)-6)))/2.d0
                  enddo
               enddo
            else
               do j=5,8
                  do k=1,nfield
                     field(k,ifacew2(j,jfaces))=
     &                   (field(k,nonei15(2,ifacew2(j,jfaces)-6))+
     &                    field(k,nonei15(3,ifacew2(j,jfaces)-6)))/2.d0
                  enddo
               enddo
            endif
         endif
!     
!     transferring the field values into yn
!     
         if(lakonl(7:8).ne.'LC') then
            do j=1,nope
               do k=1,nfield-2
                  yn(k,kon(indexe+j))=yn(k,kon(indexe+j))+
     &                 field(k,j)
               enddo
!     
!     interchanging positions of two last rows of yn  
!     
               yn(nfield-1,kon(indexe+j))=yn(nfield-1,kon(indexe+j))+
     &              field(nfield,j)
               yn(nfield,kon(indexe+j))=yn(nfield,kon(indexe+j))+
     &              field(nfield-1,j)
               inum(kon(indexe+j))=inum(kon(indexe+j))+1
            enddo
         else
            do j=1,nope*nlayer
               do k=1,nfield
                  yn(k,kon(indexe+nopeexp+j))=
     &                 yn(k,kon(indexe+nopeexp+j))+field(k,j)
               enddo
               inum(kon(indexe+nopeexp+j))=inum(kon(indexe+nopeexp+j))+1
            enddo
         endif
!     
c     Bernhardi start
c     incompatible modes elements
         if(lakonl(1:5).eq.'C3D8I') then
            do j=1,3
               do k=1,nfield
                  yn(k,kon(indexe+nope+j))=0.0d0
               enddo
c               inum(kon(indexe+nope+j))=inum(kon(indexe+nope+j))+1
            enddo
         endif
c     Bernhardi end
!     
         i=n
!     
      enddo
!     
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!     
      do i=1,nk
         if(inum(i).gt.0) then
            do j=1,nfield
               yn(j,i)=yn(j,i)/inum(i)
            enddo
         endif
      enddo
!     
!     zeroing nonactive nodes
!  
      do i=1,nslavnode(ntie+1)
         if(islavact(i).ne.1) then
            do j=1,nfield
               yn(j,islavnode(i))=0.d0
            enddo
         endif
      enddo
!
c      do i=1,nk
c         if(islavact(i).ne.1) then
c            do j=1,nfield
c               yn(j,i)=0.d0
c            enddo
c         endif
c      enddo
!     
!     for 1d and 2d elements only:
!     finding the solution in the original nodes
!     
      if((cflag.ne.' ').and.(cflag.ne.'E')) then
         call map3dto1d2d(yn,ipkon,inum,kon,lakon,nfield,nk,ne,cflag,co,
     &        vold,force,mi)
      endif
!     
      return
      end
c
c     in the next section the contact stresses are derived from the
c     volumetric stresses
c
c
c!
c      if(nfield.eq.0) return
c!
c      do i=1,nk
c         do j=1,3
c            xslavnor(j,i)=0.d0
c         enddo
c      enddo
c!
c      do i=1,nk
c         do j=1,6
c            stn(j,i)=0.d0
c         enddo
c      enddo
c!
c      i=ne0
c      do
c         i=i+1
c         if(i.gt.ne) exit
c!     
c         if(ipkon(i).lt.0) cycle
c         indexc=ipkon(i)
c!     
c         if((lakon(i)(1:1).ne.'E').or.(lakon(i)(7:7).ne.'C')) cycle
c         nopespring=kon(indexc)
c         ifaces=islavsurf(1,kon(indexc+nopespring+2))
c!     
c         nelems=int(ifaces/10.d0)
c         lakonl=lakon(nelems)
c!     
c!        determining all contact elements belonging to slave face
c!        "ifaces"
c!
c         n=i
c!     
c         do
c            n=n+1
c            if(n.gt.ne) exit
c            indexcj=ipkon(n)
c            nopespringj=kon(indexcj)
c            ifacej=islavsurf(1,kon(indexcj+nopespringj+2))
c            if(ifaces.ne.ifacej) exit
c         enddo
c         n=n-1
c         jfaces=ifaces-10*int(ifaces/10.d0)
c         indexe=ipkon(nelems)
c!
c         if(lakonl(7:8).eq.'LC') then
c            nlayer=0
c            do j=1,mi(3)
c               if(ielmat(j,i).gt.0) then
c                  nlayer=nlayer+1
c               else
c                  exit
c               endif
c            enddo
c!
c            if(lakonl(4:4).eq.'2') then
c               nopeexp=28
c            elseif(lakonl(4:5).eq.'15') then
c               nopeexp=21
c            endif
c         endif
c!     
c!     nopes: # of nodes in the slave face
c!     nope: # of nodes in the element
c!     
c         if(lakonl(4:5).eq.'8R') then
c            nopes=4
c            nope=8
c         elseif(lakonl(4:4).eq.'8') then
c            nopes=4
c            nope=8
c         elseif(lakonl(4:5).eq.'20') then
c            nopes=8
c            nope=20
c         elseif(lakonl(4:5).eq.'10') then
c            nopes=6
c            nope=10
c         elseif(lakonl(4:4).eq.'4') then
c            nopes=3
c            nope=4
c!     
c!     treatment of wedge faces
c!     
c         elseif(lakonl(4:4).eq.'6') then
c            nope=6
c            if(jfaces.le.2) then
c               nopes=3
c            else
c               nopes=4
c            endif
c         elseif(lakonl(4:5).eq.'15') then
c            nope=15
c            if(jfaces.le.2) then
c               nopes=6
c            else
c               nopes=8
c            endif
c         endif
c!  
c!        identifying the local node numbers belonging to the
c!        face
c!   
c         if((nope.eq.20).or.(nope.eq.8)) then
c            do m=1,nopes
c               nodes(m)=ifaceq(m,jfaces)
c            enddo
c         elseif((nope.eq.10).or.(nope.eq.4)) 
c     &           then
c            do m=1,nopes
c               nodes(m)=ifacet(m,jfaces)
c            enddo
c         elseif(nope.eq.15) then
c            do m=1,nopes
c               nodes(m)=ifacew2(m,jfaces)
c            enddo
c         else
c            do m=1,nopes
c               nodes(m)=ifacew1(m,jfaces)
c            enddo
c         endif
c!
c!        extrapolation of the stress tensor from the integration
c!        points to the nodes of the slave face
c!
c         if((lakonl(4:7).eq.'20RB').and.
c     &        (lakonl(8:8).ne.'R').and.
c     &        (lakonl(8:8).ne.'C')) then
c            call beamextscheme(yi(1,1,nelems),ndim,nfield,lakonl,
c     &              ielprop(nelems),prop,field,mi)
c!               
cc     Bernhardi start
c         elseif((lakonl(4:6).eq.'20R').or.
c     &           (lakonl(4:5).eq.'8 ').or.(lakonl(4:5).eq.'8I')) then
cc     Bernhardi end
c            if(lakonl(7:8).ne.'LC') then
c               do j=1,8
c                  do k=1,nfield
c                     field(k,j)=0.d0
c                     do l=1,8
c                        field(k,j)=field(k,j)+a8(j,l)*yi(k,l,nelems)
c                     enddo
c                  enddo
c               enddo
c            else
c               do m=1,nlayer
c                  jj=20*(m-1)
c                  ll=8*(m-1)
c                  do j=1,8
c                     do k=1,nfield
c                        field(k,jj+j)=0.d0
c                        do l=1,8
c                           field(k,jj+j)=
c     &                          field(k,jj+j)+a8(j,l)*yi(k,ll+l,nelems)
c                        enddo
c                     enddo
c                  enddo
c               enddo
c            endif
c         elseif(lakonl(4:4).eq.'8') then
c            do j=1,8
c               do k=1,nfield
c                  field(k,j)=yi(k,1,nelems)
c               enddo
c            enddo
c         elseif(lakonl(4:5).eq.'10') then
c            do j=1,4
c               do k=1,nfield
c                  field(k,j)=0.d0
c                  do l=1,4
c                     field(k,j)=field(k,j)+a4(j,l)*yi(k,l,nelems)
c                  enddo
c               enddo
c            enddo
c         elseif(lakonl(4:4).eq.'2') then
c            do j=1,20
c               do k=1,nfield
c                  field(k,j)=0.d0
c                  do l=1,27
c                     field(k,j)=field(k,j)+a27(j,l)*yi(k,l,nelems)
c                  enddo
c               enddo
c            enddo
c         elseif(lakonl(4:4).eq.'4') then
c            do j=1,4
c               do k=1,nfield
c                  field(k,j)=yi(k,1,nelems)
c               enddo
c            enddo
c         elseif(lakonl(4:4).eq.'1') then
c            do j=1,6
c               do k=1,nfield
c                  field(k,j)=0.d0
c                  do l=1,9
c                     field(k,j)=field(k,j)+a9(j,l)*yi(k,l,nelems)
c                  enddo
c               enddo
c            enddo
c         else
c            do j=1,6
c               do k=1,nfield
c                  field(k,j)=0.d0
c                  do l=1,2
c                     field(k,j)=field(k,j)+a2(j,l)*yi(k,l,nelems)
c                  enddo
c               enddo
c            enddo
c         endif
c!     
c!     determining the field values in the midside nodes
c!     
c         if(lakonl(4:6).eq.'20R') then
c            if(lakonl(7:8).ne.'LC') then
c               do j=9,20
c                  do k=1,nfield
c                     field(k,j)=(field(k,nonei20(2,j-8))+
c     &                    field(k,nonei20(3,j-8)))/2.d0
c                  enddo
c               enddo
c            else
c               do m=1,nlayer
c                  jj=20*(m-1)
c                  do j=9,20
c                     do k=1,nfield
c                        field(k,jj+j)=(field(k,jj+nonei20(2,j-8))
c     &                       +field(k,jj+nonei20(3,j-8)))/2.d0
c                     enddo
c                  enddo
c               enddo
c            endif
c         elseif(lakonl(4:5).eq.'10') then
c            do j=5,10
c               do k=1,nfield
c                  field(k,j)=(field(k,nonei10(2,j-4))+
c     &                 field(k,nonei10(3,j-4)))/2.d0
c               enddo
c            enddo
c         elseif(lakonl(4:5).eq.'15') then
c            do j=7,15
c               do k=1,nfield
c                  field(k,j)=(field(k,nonei15(2,j-6))+
c     &                 field(k,nonei15(3,j-6)))/2.d0
c               enddo
c            enddo
c         endif
c!     
c!     transferring the field values into stn
c!     
c         if(lakonl(7:8).ne.'LC') then
c            do j=1,nopes
c               m=nodes(j)
c               do k=1,nfield
c                  stn(k,kon(indexe+m))=stn(k,kon(indexe+m))+
c     &                 field(k,m)
c               enddo
c               inum(kon(indexe+m))=inum(kon(indexe+m))+1
c            enddo
c         else
c            do j=1,nope*nlayer
c               do k=1,nfield
c                  stn(k,kon(indexe+nopeexp+j))=
c     &                 stn(k,kon(indexe+nopeexp+j))+field(k,j)
c               enddo
c               inum(kon(indexe+nopeexp+j))=inum(kon(indexe+nopeexp+j))+1
c            enddo
c         endif
c!     
cc     Bernhardi start
cc     incompatible modes elements
c         if(lakonl(1:5).eq.'C3D8I') then
c            do j=1,3
c               do k=1,nfield
c                  stn(k,kon(indexe+nope+j))=0.0d0
c               enddo
c            enddo
c         endif
cc     Bernhardi end
c!     
c!     calculation of the normals at the nodes of the slave face
c!     
c!     actual position of the nodes belonging to the
c!     slave surface
c!     
c         do k=1,nope
c            konl(k)=kon(indexe+k)
c         enddo
c!     
c         do m=1,nopes
c            do k=1,3
c               xl2s(k,m)=co(k,konl(nodes(m)))+
c     &              vold(k,konl(nodes(m)))
c            enddo
c         enddo
c!         
c!     calculate the normal vector in the nodes belonging to the slave surface
c!     
c         if(nopes.eq.8) then
c            do m=1,nopes
c               xi=xquad(1,m)
c               et=xquad(2,m)
c               call shape8q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
c               dl=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2)
c     &              + xsj2(3)*xsj2(3))
c               xsj2(1)=xsj2(1)/dl
c               xsj2(2)=xsj2(2)/dl
c               xsj2(3)=xsj2(3)/dl
c!     
c               node=konl(nodes(m))
c!     
c               do k=1,3
c                  xslavnor(k,node)=xslavnor(k,node)+xsj2(k)
c               enddo
c!     
c            enddo
c         elseif(nopes.eq.4) then
c            do m=1,nopes
c               xi=xquad(1,m)
c               et=xquad(2,m)
c               call shape4q(xi,et,xl2s,xsj2,xs2,shp2,iflag)
c               dl=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
c     &              + xsj2(3)*xsj2(3))
c               xsj2(1)=xsj2(1)/dl
c               xsj2(2)=xsj2(2)/dl
c               xsj2(3)=xsj2(3)/dl
c!     
c               node=konl(nodes(m))
c!     
c               do k=1,3
c                  xslavnor(k,node)=xslavnor(k,node)+xsj2(k)
c               enddo
c!     
c            enddo
c         elseif(nopes.eq.6) then
c            do m=1,nopes
c               xi=xtri(1,m)
c               et=xtri(2,m)
c               call shape6tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
c               dl=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
c     &              + xsj2(3)*xsj2(3))
c               xsj2(1)=xsj2(1)/dl
c               xsj2(2)=xsj2(2)/dl
c               xsj2(3)=xsj2(3)/dl
c!     
c               node=konl(nodes(m))
c!     
c               do k=1,3
c                  xslavnor(k,node)=xslavnor(k,node)+xsj2(k)
c               enddo
c!     
c            enddo
c         else
c            do m=1,nopes
c               xi=xtri(1,m)
c               et=xtri(2,m)
c               call shape3tri(xi,et,xl2s,xsj2,xs2,shp2,iflag)
c               dl=dsqrt(xsj2(1)*xsj2(1) + xsj2(2)*xsj2(2) 
c     &              + xsj2(3)*xsj2(3))
c               xsj2(1)=xsj2(1)/dl
c               xsj2(2)=xsj2(2)/dl
c               xsj2(3)=xsj2(3)/dl
c!     
c               node=konl(nodes(m))
c!     
c               do k=1,3
c                  xslavnor(k,node)=xslavnor(k,node)+xsj2(k)
c               enddo
c!     
c            enddo
c         endif
c!     
c         i=n
c!     
c      enddo
c!
c!     taking the mean of nodal contributions coming from different
c!     elements having the node in common
c!
c      do i=1,nk
c         if(inum(i).gt.0) then
c            do j=1,6
c               stn(j,i)=stn(j,i)/inum(i)
c            enddo
c         endif
c      enddo
c!
c!     calculating the traction in the slave nodes
c!
c      do i=1,nk
c         if(inum(i).ne.0) then
c!
c!           determining the mean normal
c!
c
c            do j=1,3
c               xn(j)=xslavnor(j,i)
c            enddo
c            dl=dsqrt(xn(1)*xn(1)+xn(2)*xn(2)+xn(3)*xn(3))
c            do j=1,3
c               xn(j)=xn(j)/dl
c            enddo
c!
c!           determining the tangent unit vectors
c!
c            if(1.d0-dabs(xn(1)).lt.1.5231d-6) then       
c               t1(1)=-xn(3)*xn(1)
c               t1(2)=-xn(3)*xn(2)
c               t1(3)=1.d0-xn(3)*xn(3)
c            else
c               t1(1)=1.d0-xn(1)*xn(1)
c               t1(2)=-xn(1)*xn(2)
c               t1(3)=-xn(1)*xn(3)
c            endif
c            dt1=dsqrt(t1(1)*t1(1)+t1(2)*t1(2)+t1(3)*t1(3))
c            do j=1,3
c               t1(j)=t1(j)/dt1
c            enddo
c            t2(1)=xn(2)*t1(3)-xn(3)*t1(2)
c            t2(2)=xn(3)*t1(1)-xn(1)*t1(3)
c            t2(3)=xn(1)*t1(2)-xn(2)*t1(1)           
c!
c!           calculating the traction
c!
c            trac(1)=stn(1,i)*xn(1)+stn(4,i)*xn(2)+stn(5,i)*xn(3)
c            trac(2)=stn(4,i)*xn(1)+stn(2,i)*xn(2)+stn(6,i)*xn(3)
c            trac(3)=stn(5,i)*xn(1)+stn(6,i)*xn(2)+stn(3,i)*xn(3)
c!
c!           determining the contact pressure
c!
c            yn(4,i)=-(trac(1)*xn(1)+trac(2)*xn(2)+trac(3)*xn(3))
c!
c!           determining the contact shear stress components
c!
c            yn(6,i)=-(trac(1)*t1(1)+trac(2)*t1(2)+trac(3)*t1(3))
c            yn(5,i)=trac(1)*t2(1)+trac(2)*t2(2)+trac(3)*t2(3)
c!
c         endif
c      enddo
c!     
c!     zeroing nonactive nodes
c!  
c      do i=1,nslavnode(ntie+1)
c         if(islavact(i).ne.1) then
c            do j=1,nfield
c               yn(j,islavnode(i))=0.d0
c            enddo
c         endif
c      enddo
c!     
c!     for 1d and 2d elements only:
c!     finding the solution in the original nodes
c!     
c      if((cflag.ne.' ').and.(cflag.ne.'E')) then
c         call map3dto1d2d(yn,ipkon,inum,kon,lakon,nfield,nk,ne,cflag,co,
c     &        vold,force,mi)
c      endif
c!     
c      return
c      end
