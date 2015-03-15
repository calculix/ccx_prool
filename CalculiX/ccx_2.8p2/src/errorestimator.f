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
      subroutine errorestimator(yi,yn,ipkon,inum,kon,lakon,nk,
     &  ne,mi,ielmat,thicke,nterms)
!
!     calculates an estimate of the error for the worst principal
!     stress and the von mises stress for mechanical calculations
!     (nterms=6).
!     For thermal calculations the error estimator is based on the
!     heat flux vector (nterms=3).
!
!     each node belongs to n different elements. The stress tensor
!     is extrapolated from each of these elements towards the node
!     and the corresponding worst principal stress is calculated. In
!     that way one gets n different worst principal stress values. The
!     error is the maximum difference between these values.
!
      implicit none
!
      character*8 lakon(*),lakonl
!
      integer ipkon(*),inum(*),kon(*),mi(*),ne,indexe,nope,
     &  nonei20(3,12),nfield,nonei10(3,6),nk,i,j,k,l,
     &  nonei15(3,9),konl,nopev,nterms,ishift,
     &  mint3d,m,iflag,node,jj,ll,ielmat(mi(3),*),
     &  nlayer,nopeexp,ilayer,kk,mint2d,nopes,kl,ki,kflag,three
!
      real*8 yi(nterms,mi(1),*),yn(nterms,*),field(999,20*mi(3)),
     &  a8(8,8),size,
     &  a4(4,4),a27(20,27),a9(6,9),a2(6,2),al(3),worstprin,
     &  xi,et,ze,xl(3,20),xsj,shp(4,20),
     &  yiloc(6,27),a(3,3),b(3,3),c(3,3),tlayer(4),
     &  dlayer(4),xlayer(mi(3),4),thickness,xs2(3,7),xl2(3,8),
     &  xsj2(3),shp2(7,8),thicke(mi(3),*),v1,vmises
!
      include "gauss.f"
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
      data a2 /  1.1455,-0.1455,1.1455,-0.1455,1.1455,-0.1455,
     &           -0.1455,1.1455,-0.1455,1.1455,-0.1455,1.1455/
      data a4 /  1.92705, -0.30902, -0.30902, -0.30902,
     &          -0.30902,  1.92705, -0.30902, -0.30902,
     &          -0.30902, -0.30902,  1.92705, -0.30902,
     &          -0.30902, -0.30902, -0.30902,  1.92705/
      data a9 / 1.63138,-0.32628,-0.32628,-0.52027, 0.10405, 0.10405,
     &         -0.32628, 1.63138,-0.32628, 0.10405,-0.52027, 0.10405,
     &         -0.32628,-0.32628, 1.63138, 0.10405, 0.10405,-0.52027,
     &          0.55556,-0.11111,-0.11111, 0.55556,-0.11111,-0.11111,
     &         -0.11111, 0.55556,-0.11111,-0.11111, 0.55556,-0.11111,
     &         -0.11111,-0.11111, 0.55556,-0.11111,-0.11111, 0.55556,
     &         -0.52027, 0.10405, 0.10405, 1.63138,-0.32628,-0.32628,
     &          0.10405,-0.52027, 0.10405,-0.32628, 1.63138,-0.32628,
     &          0.10405, 0.10405,-0.52027,-0.32628,-0.32628, 1.63138/
      data a8 /2.549,-.683,.183,-.683,-.683,.183,
     &        -.04904,.183,-.683,2.549,-.683,.183,
     &        .183,-.683,.183,-.04904,-.683,.183,
     &        -.683,2.549,.183,-.04904,.183,-.683,
     &        .183,-.683,2.549,-.683,-.04904,.183,
     &        -.683,.183,-.683,.183,-.04904,.183,
     &        2.549,-.683,.183,-.683,.183,-.683,
     &        .183,-.04904,-.683,2.549,-.683,.183,
     &        .183,-.04904,.183,-.683,-.683,.183,
     &        -.683,2.549,-.04904,.183,-.683,.183,
     &        .183,-.683,2.549,-.683/      
      data a27 /
     &  2.37499,-0.12559,-0.16145,-0.12559,-0.12559,-0.16145, 0.11575,
     & -0.16145, 0.32628, 0.11111, 0.11111, 0.32628, 0.11111,-0.10405,
     & -0.10405, 0.11111, 0.32628, 0.11111,-0.10405, 0.11111,-0.31246,
     & -0.31246, 0.31481, 0.31481, 0.31481, 0.31481,-0.16902,-0.16902,
     &  1.28439,-0.27072,-0.19444,-0.27072,-0.19444, 0.15961,-0.00661,
     &  0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.12559, 2.37499,
     & -0.12559,-0.16145,-0.16145,-0.12559,-0.16145, 0.11575, 0.32628,
     &  0.32628, 0.11111, 0.11111, 0.11111, 0.11111,-0.10405,-0.10405,
     &  0.11111, 0.32628, 0.11111,-0.10405,-0.31246, 0.31481, 0.31481,
     & -0.31246, 0.31481,-0.16902,-0.16902, 0.31481,-0.27072,-0.19444,
     & -0.27072, 1.28439, 0.15961,-0.00661, 0.15961,-0.19444,-0.27072,
     &  0.15961, 0.15961,-0.27072,-0.48824,-0.48824,-0.48824,-0.48824,
     &  0.22898, 0.22898, 0.22898, 0.22898, 0.05556, 0.05556, 0.05556,
     &  0.05556, 0.05556, 0.05556, 0.05556, 0.05556,-0.22222,-0.22222,
     & -0.22222,-0.22222, 0.31481,-0.31246,-0.31246, 0.31481,-0.16902,
     &  0.31481, 0.31481,-0.16902,-0.27072, 1.28439,-0.27072,-0.19444,
     &  0.15961,-0.19444, 0.15961,-0.00661, 0.15961,-0.27072,-0.27072,
     &  0.15961,-0.12559,-0.16145,-0.12559, 2.37499,-0.16145, 0.11575,
     & -0.16145,-0.12559, 0.11111, 0.11111, 0.32628, 0.32628,-0.10405,
     & -0.10405, 0.11111, 0.11111, 0.11111,-0.10405, 0.11111, 0.32628,
     &  0.31481, 0.31481,-0.31246,-0.31246,-0.16902,-0.16902, 0.31481,
     &  0.31481,-0.19444,-0.27072, 1.28439,-0.27072,-0.00661, 0.15961,
     & -0.19444, 0.15961, 0.15961, 0.15961,-0.27072,-0.27072,-0.16145,
     & -0.12559, 2.37499,-0.12559, 0.11575,-0.16145,-0.12559,-0.16145,
     &  0.11111, 0.32628, 0.32628, 0.11111,-0.10405, 0.11111, 0.11111,
     & -0.10405,-0.10405, 0.11111, 0.32628, 0.11111,-0.31246, 0.31481,
     & -0.16902, 0.31481,-0.31246, 0.31481,-0.16902, 0.31481,-0.27072,
     &  0.15961, 0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.27072,
     &  1.28439,-0.19444,-0.00661,-0.19444,-0.48824,-0.48824, 0.22898,
     &  0.22898,-0.48824,-0.48824, 0.22898, 0.22898, 0.05556,-0.22222,
     &  0.05556,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,
     &  0.05556, 0.05556, 0.05556, 0.31481,-0.31246, 0.31481,-0.16902,
     &  0.31481,-0.31246, 0.31481,-0.16902,-0.27072,-0.27072, 0.15961,
     &  0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.19444, 1.28439,
     & -0.19444,-0.00661,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824,
     &  0.22898, 0.22898,-0.48824,-0.22222, 0.05556,-0.22222, 0.05556,
     & -0.22222, 0.05556,-0.22222, 0.05556, 0.05556, 0.05556, 0.05556,
     &  0.05556,-0.29630,-0.29630,-0.29630,-0.29630,-0.29630,-0.29630,
     & -0.29630,-0.29630,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,
     & -0.11111,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,-0.11111,
     &  0.22898,-0.48824,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824,
     &  0.22898,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,
     & -0.22222, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.31481,
     & -0.16902, 0.31481,-0.31246, 0.31481,-0.16902, 0.31481,-0.31246,
     &  0.15961, 0.15961,-0.27072,-0.27072, 0.15961, 0.15961,-0.27072,
     & -0.27072,-0.19444,-0.00661,-0.19444, 1.28439, 0.22898, 0.22898,
     & -0.48824,-0.48824, 0.22898, 0.22898,-0.48824,-0.48824, 0.05556,
     & -0.22222, 0.05556,-0.22222, 0.05556,-0.22222, 0.05556,-0.22222,
     &  0.05556, 0.05556, 0.05556, 0.05556,-0.16902, 0.31481,-0.31246,
     &  0.31481,-0.16902, 0.31481,-0.31246, 0.31481, 0.15961,-0.27072,
     & -0.27072, 0.15961, 0.15961,-0.27072,-0.27072, 0.15961,-0.00661,
     & -0.19444, 1.28439,-0.19444,-0.12559,-0.16145, 0.11575,-0.16145,
     &  2.37499,-0.12559,-0.16145,-0.12559, 0.11111,-0.10405,-0.10405,
     &  0.11111, 0.32628, 0.11111, 0.11111, 0.32628, 0.32628, 0.11111,
     & -0.10405, 0.11111, 0.31481, 0.31481,-0.16902,-0.16902,-0.31246,
     & -0.31246, 0.31481, 0.31481,-0.19444, 0.15961,-0.00661, 0.15961,
     &  1.28439,-0.27072,-0.19444,-0.27072,-0.27072,-0.27072, 0.15961,
     &  0.15961,-0.16145,-0.12559,-0.16145, 0.11575,-0.12559, 2.37499,
     & -0.12559,-0.16145, 0.11111, 0.11111,-0.10405,-0.10405, 0.32628,
     &  0.32628, 0.11111, 0.11111, 0.11111, 0.32628, 0.11111,-0.10405,
     &  0.31481,-0.16902,-0.16902, 0.31481,-0.31246, 0.31481, 0.31481,
     & -0.31246, 0.15961,-0.00661, 0.15961,-0.19444,-0.27072,-0.19444,
     & -0.27072, 1.28439,-0.27072, 0.15961, 0.15961,-0.27072, 0.22898,
     &  0.22898, 0.22898, 0.22898,-0.48824,-0.48824,-0.48824,-0.48824,
     &  0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556, 0.05556,
     &  0.05556,-0.22222,-0.22222,-0.22222,-0.22222,-0.16902, 0.31481,
     &  0.31481,-0.16902, 0.31481,-0.31246,-0.31246, 0.31481, 0.15961,
     & -0.19444, 0.15961,-0.00661,-0.27072, 1.28439,-0.27072,-0.19444,
     &  0.15961,-0.27072,-0.27072, 0.15961,-0.16145, 0.11575,-0.16145,
     & -0.12559,-0.12559,-0.16145,-0.12559, 2.37499,-0.10405,-0.10405,
     &  0.11111, 0.11111, 0.11111, 0.11111, 0.32628, 0.32628, 0.11111,
     & -0.10405, 0.11111, 0.32628,-0.16902,-0.16902, 0.31481, 0.31481,
     &  0.31481, 0.31481,-0.31246,-0.31246,-0.00661, 0.15961,-0.19444,
     &  0.15961,-0.19444,-0.27072, 1.28439,-0.27072, 0.15961, 0.15961,
     & -0.27072,-0.27072, 0.11575,-0.16145,-0.12559,-0.16145,-0.16145,
     & -0.12559, 2.37499,-0.12559,-0.10405, 0.11111, 0.11111,-0.10405,
     &  0.11111, 0.32628, 0.32628, 0.11111,-0.10405, 0.11111, 0.32628,
     &  0.11111/
!
      data iflag /1/
      data kflag /1/
      data three /3/
!
      do i=1,nk
         do j=1,nterms
            yn(j,i)=0.d0
         enddo
      enddo
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(7:8).eq.'LC') then
            nlayer=0
            do j=1,mi(3)
               if(ielmat(j,i).gt.0) then
                  nlayer=nlayer+1
               else
                  exit
               endif
            enddo
!
            if(lakonl(4:4).eq.'2') then
               nopeexp=28
            elseif(lakonl(4:5).eq.'15') then
               nopeexp=21
            endif
         endif
!
         if(lakonl(1:1).eq.'F') then
            cycle
         elseif(lakonl(4:4).eq.'2') then
            nope=20
            nopev=8
         elseif(lakonl(4:4).eq.'8') then
            nope=8
            nopev=8
         elseif(lakonl(4:5).eq.'10') then
            nope=10
            nopev=4
         elseif(lakonl(4:4).eq.'4') then
            nope=4
            nopev=4
         elseif(lakonl(4:5).eq.'15') then
            nope=15
            nopev=6
         elseif(lakonl(4:4).eq.'6') then
            nope=6
            nopev=6
         else
            cycle
         endif
!
!
!        determining the field values in the vertex nodes
!        for C3D20R and C3D8: trilinear extrapolation (= use of the
!                             C3D8 linear interpolation functions)
!        for C3D8R: constant field value in each element
!        for C3D10: use of the C3D4 linear interpolation functions
!        for C3D4: constant field value in each element
!        for C3D15: use of the C3D6 linear interpolation functions
!        for C3D6: use of a linear interpolation function
!
c     Bernhardi start
         if((lakonl(4:6).eq.'20R').or.(lakonl(4:6).eq.'26R').or.
     &          (lakonl(4:5).eq.'8 ').or.(lakonl(4:5).eq.'8I')) then
c     Bernhardi end
            if(lakonl(7:8).ne.'LC') then
               do j=1,8
                  do k=1,nterms
                     field(k,j)=0.d0
                     do l=1,8
                        field(k,j)=field(k,j)+a8(j,l)*yi(k,l,i)
                     enddo
                  enddo
               enddo
            else
               do m=1,nlayer
                  jj=20*(m-1)
                  ll=8*(m-1)
                  do j=1,8
                     do k=1,nterms
                        field(k,jj+j)=0.d0
                        do l=1,8
                           field(k,jj+j)=
     &                          field(k,jj+j)+a8(j,l)*yi(k,ll+l,i)
                        enddo
                     enddo
                  enddo
               enddo
            endif
         elseif(lakonl(4:4).eq.'8') then
            do j=1,8
               do k=1,nterms
                  field(k,j)=yi(k,1,i)
               enddo
            enddo
         elseif(lakonl(4:5).eq.'10') then
            do j=1,4
               do k=1,nterms
                  field(k,j)=0.d0
                  do l=1,4
                     field(k,j)=field(k,j)+a4(j,l)*yi(k,l,i)
                  enddo
               enddo
            enddo
         elseif(lakonl(4:4).eq.'2') then
            do j=1,20
               do k=1,nterms
                  field(k,j)=0.d0
                  do l=1,27
                     field(k,j)=field(k,j)+a27(j,l)*yi(k,l,i)
                  enddo
               enddo
            enddo
         elseif(lakonl(4:4).eq.'4') then
            do j=1,4
               do k=1,nterms
                  field(k,j)=yi(k,1,i)
               enddo
            enddo
         elseif(lakonl(4:4).eq.'1') then
            do j=1,6
               do k=1,nterms
                  field(k,j)=0.d0
                  do l=1,9
                     field(k,j)=field(k,j)+a9(j,l)*yi(k,l,i)
                  enddo
               enddo
            enddo
         else
            do j=1,6
               do k=1,nterms
                  field(k,j)=0.d0
                  do l=1,2
                     field(k,j)=field(k,j)+a2(j,l)*yi(k,l,i)
                  enddo
               enddo
            enddo
         endif
!     
!     transferring the field values into yn
!     
         if(lakonl(7:8).ne.'LC') then
            if(nterms.eq.6) then
               do j=1,nopev
                  c(1,1)=field(1,j)
                  c(2,2)=field(2,j)
                  c(3,3)=field(3,j)
                  c(1,2)=field(4,j)
                  c(1,3)=field(5,j)
                  c(2,3)=field(6,j)
!     
!     calculate the eigenvalues
!     
                  call calceigenvalues(c,al)
!     
                  if(dabs(al(3)).gt.dabs(al(1))) then
                     worstprin=al(3)
                  else
                     worstprin=al(1)
                  endif
                  yn(1,kon(indexe+j))=yn(1,kon(indexe+j))+worstprin
                  yn(2,kon(indexe+j))=yn(2,kon(indexe+j))+worstprin**2
                  yn(3,kon(indexe+j))=yn(3,kon(indexe+j))+1
!
!     calculate the von Mises stress
!     
                  v1=(c(1,1)+c(2,2)+c(3,3))/3.d0
                  c(1,1)=c(1,1)-v1/3.d0
                  c(2,2)=c(2,2)-v1/3.d0
                  c(3,3)=c(3,3)-v1/3.d0
                  vmises=c(1,1)*c(1,1)+c(2,2)*c(2,2)+c(3,3)*c(3,3)+
     &                 2.d0*(c(1,2)*c(1,2)+c(2,3)*c(2,3)+c(1,3)*c(1,3))
                  vmises=dsqrt(3.d0*vmises/2.d0)
                  yn(4,kon(indexe+j))=yn(4,kon(indexe+j))+vmises
                  yn(6,kon(indexe+j))=yn(6,kon(indexe+j))+vmises**2
               enddo
            else
               do j=1,nopev
                  size=dsqrt(field(1,j)**2+field(2,j)**2+field(3,j)**2)
                  yn(1,kon(indexe+j))=yn(1,kon(indexe+j))+size
                  yn(2,kon(indexe+j))=yn(2,kon(indexe+j))+size**2
                  yn(3,kon(indexe+j))=yn(3,kon(indexe+j))+1
               enddo
            endif
         else
            if(nterms.eq.6) then
               do j=1,nope*nlayer
                  c(1,1)=field(1,j)
                  c(2,2)=field(2,j)
                  c(3,3)=field(3,j)
                  c(1,2)=field(4,j)
                  c(1,3)=field(5,j)
                  c(2,3)=field(6,j)
!     
                  call calceigenvalues(c,al)
!     
                  if(dabs(al(3)).gt.dabs(al(1))) then
                     worstprin=al(3)
                  else
                     worstprin=al(1)
                  endif
                  yn(1,kon(indexe+nopeexp+j))=
     &                 yn(1,kon(indexe+nopeexp+j))+worstprin
                  yn(2,kon(indexe+nopeexp+j))=
     &                 yn(2,kon(indexe+nopeexp+j))+worstprin**2
                  yn(3,kon(indexe+nopeexp+j))=
     &                 yn(3,kon(indexe+nopeexp+j))+1
!     
!                 calculate the von Mises stress
!
                  v1=(c(1,1)+c(2,2)+c(3,3))/3.d0
                  c(1,1)=c(1,1)-v1/3.d0
                  c(2,2)=c(2,2)-v1/3.d0
                  c(3,3)=c(3,3)-v1/3.d0
                  vmises=c(1,1)*c(1,1)+c(2,2)*c(2,2)+c(3,3)*c(3,3)+
     &                 2.d0*(c(1,2)*c(1,2)+c(2,3)*c(2,3)+c(1,3)*c(1,3))
                  vmises=dsqrt(3.d0*vmises/2.d0)
                  yn(4,kon(indexe+nopeexp+j))=
     &                 yn(4,kon(indexe+nopeexp+j))+vmises
                  yn(6,kon(indexe+nopeexp+j))=
     &                 yn(6,kon(indexe+nopeexp+j))+vmises**2
               enddo
            else
               do j=1,nope*nlayer
                  size=dsqrt(field(1,j)**2+field(2,j)**2+field(3,j)**2)
                  yn(1,kon(indexe+j))=yn(1,kon(indexe+j))+size
                  yn(2,kon(indexe+j))=yn(2,kon(indexe+j))+size**2
                  yn(3,kon(indexe+j))=yn(3,kon(indexe+j))+1
               enddo
            endif
         endif
!
      enddo
!
!     taking the mean of nodal contributions coming from different
!     elements having the node in common
!
      if(nterms.eq.6) then
         do i=1,nk
            if(yn(3,i).gt.1) then
               yn(5,i)=(yn(6,i)-(yn(4,i)**2)/yn(3,i))/
     &              (yn(3,i)-1.d0)
               if(yn(5,i).gt.0.d0) then
                  yn(5,i)=dsqrt(yn(5,i))
               else
                  yn(5,i)=0.d0
               endif
               yn(3,i)=(yn(2,i)-(yn(1,i)**2)/yn(3,i))/
     &              (yn(3,i)-1.d0)
               if(yn(3,i).gt.0.d0) then
                  yn(3,i)=dsqrt(yn(3,i))
               else
                  yn(3,i)=0.d0
               endif
            else
               yn(3,i)=0.d0
               yn(5,i)=0.d0
            endif
            yn(1,i)=0.d0
            yn(2,i)=0.d0
            yn(4,i)=0.d0
            yn(6,i)=0.d0
         enddo
      else
         do i=1,nk
            if(yn(3,i).gt.1) then
               yn(3,i)=(yn(2,i)-(yn(1,i)**2)/yn(3,i))/
     &              (yn(3,i)-1.d0)
               if(yn(3,i).gt.0.d0) then
                  yn(3,i)=dsqrt(yn(3,i))
               else
                  yn(3,i)=0.d0
               endif
            else
               yn(3,i)=0.d0
            endif
            yn(1,i)=0.d0
            yn(2,i)=0.d0
         enddo
      endif
!
!        determining the field values in the midside nodes
!
      if(nterms.eq.6) then
         ishift=2
      else
         ishift=4
      endif
!
      do i=1,ne
!
         if(ipkon(i).lt.0) cycle
         indexe=ipkon(i)
         lakonl=lakon(i)
!
         if(lakonl(7:8).eq.'LC') then
            nlayer=0
            do j=1,mi(3)
               if(ielmat(j,i).gt.0) then
                  nlayer=nlayer+1
               else
                  exit
               endif
            enddo
!
            if(lakonl(4:4).eq.'2') then
               nopeexp=28
            elseif(lakonl(4:5).eq.'15') then
               nopeexp=21
            endif
         endif
!
         if(lakonl(4:5).eq.'20') then
            if(lakonl(7:8).ne.'LC') then
               do j=9,20
                  do k=3,5,ishift
                     yn(k,kon(indexe+j))=(
     &                    yn(k,kon(indexe+nonei20(2,j-8)))+
     &                    yn(k,kon(indexe+nonei20(3,j-8))))/2.d0
                  enddo
               enddo
            else
               do m=1,nlayer
                  jj=20*(m-1)
                  do j=9,20
                     do k=3,5,ishift
                        yn(k,kon(indexe+nopeexp+jj+j))=(
     &                      yn(k,kon(indexe+nopeexp+jj+nonei20(2,j-8)))+
     &                      yn(k,kon(indexe+nopeexp+jj+nonei20(3,j-8))))
     &                      /2.d0
                     enddo
                  enddo
               enddo
            endif
         elseif(lakonl(4:5).eq.'10') then
            do j=5,10
               do k=3,5,ishift
                  yn(k,kon(indexe+j))=(yn(k,kon(indexe+nonei10(2,j-4)))+
     &                 yn(k,kon(indexe+nonei10(3,j-4))))/2.d0
               enddo
            enddo
         elseif(lakonl(4:5).eq.'15') then
            do j=7,15
               do k=3,5,ishift
                  yn(k,kon(indexe+j))=(yn(k,kon(indexe+nonei15(2,j-6)))+
     &                 yn(k,kon(indexe+nonei15(3,j-6))))/2.d0
               enddo
            enddo
         endif
      enddo
!     
      return
      end
