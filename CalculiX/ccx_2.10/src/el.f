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
      subroutine el(elas,eloc,alph,tt,ithermal,icmd,beta,stre,mattyp)
!
!     calculates stiffness and stresses for linear elastic materials
!
!     icmd=1: calculation of the stress corresponding to the
!             thermal strain
!     icmd=2: calculation of the stress at the mechanical strain
!
      implicit none
!
      integer ithermal,mattyp,icmd,m1,m2,m3,m4
!
      real*8 elas(21),alph(6),eloc(*),beta(*),stre(*),e,um,un,
     &  al,am1,anisox(3,3,3,3),tt,ekl(3,3),tkl(3,3),exx,eyy,ezz,
     &  exy,exz,eyz
!
      if((icmd.eq.1).and.(ithermal.eq.0)) return
!
!     reformulating the elastic properties
!
      if(mattyp.eq.1) then
         e=elas(1)
         un=elas(2)
         al=un*e/(1.d0+un)/(1.d0-2.d0*un)
         um=e/2.d0/(1.d0+un)
         am1=al+2.d0*um
      elseif(mattyp.eq.3) then
         call anisotropic(elas,anisox)
      endif
!
!     calculating the stress corresponding to the thermal strain
!     and adding it to the initial stresses
!
      if((ithermal.eq.1).or.(ithermal.eq.3)) then
         if(mattyp.eq.1) then
            beta(1)=((3.d0*al+2.d0*um)*alph(1))*tt+
     &           beta(1)
            beta(2)=((3.d0*al+2.d0*um)*alph(1))*tt+
     &           beta(2)
            beta(3)=((3.d0*al+2.d0*um)*alph(1))*tt+
     &           beta(3)
         elseif(mattyp.eq.2) then
            beta(1)=(elas(1)*alph(1)+
     &           elas(2)*alph(2)+
     &           elas(4)*alph(3))*tt+beta(1)
            beta(2)=(elas(2)*alph(1)+
     &           elas(3)*alph(2)+
     &           elas(5)*alph(3))*tt+beta(2)
            beta(3)=(elas(4)*alph(1)+
     &           elas(5)*alph(2)+
     &           elas(6)*alph(3))*tt+beta(3)
         elseif(mattyp.eq.3) then
            beta(1)=(elas(1)*alph(1)+
     &           elas(2)*alph(2)+
     &           elas(4)*alph(3)+
     &           elas(7)*alph(4)+
     &           elas(11)*alph(5)+
     &           elas(16)*alph(6))*tt+beta(1)
            beta(2)=(elas(2)*alph(1)+
     &           elas(3)*alph(2)+
     &           elas(5)*alph(3)+
     &           elas(6)*alph(4)+
     &           elas(12)*alph(5)+
     &           elas(17)*alph(6))*tt+beta(2)
            beta(3)=(elas(4)*alph(1)+
     &           elas(5)*alph(2)+
     &           elas(6)*alph(3)+
     &           elas(9)*alph(4)+
     &           elas(13)*alph(5)+
     &           elas(18)*alph(6))*tt+beta(3)
            beta(4)=(elas(7)*alph(1)+
     &           elas(8)*alph(2)+
     &           elas(9)*alph(3)+
     &           elas(10)*alph(4)+
     &           elas(14)*alph(5)+
     &           elas(19)*alph(6))*tt+beta(4)
            beta(5)=(elas(11)*alph(1)+
     &           elas(12)*alph(2)+
     &           elas(13)*alph(3)+
     &           elas(14)*alph(4)+
     &           elas(15)*alph(5)+
     &           elas(20)*alph(6))*tt+beta(5)
            beta(6)=(elas(16)*alph(1)+
     &           elas(17)*alph(2)+
     &           elas(18)*alph(3)+
     &           elas(19)*alph(4)+
     &           elas(20)*alph(5)+
     &           elas(21)*alph(6))*tt+beta(6)
         endif
      endif
!
!     calculating the stress at mechanical strain
!
      if(icmd.eq.2) then
         if(mattyp.eq.1) then
            exx=eloc(1)
            eyy=eloc(2)
            ezz=eloc(3)
            exy=2.d0*eloc(4)
            exz=2.d0*eloc(5)
            eyz=2.d0*eloc(6)
            stre(1)=am1*exx+al*(eyy+ezz)-beta(1)
            stre(2)=am1*eyy+al*(exx+ezz)-beta(1)
            stre(3)=am1*ezz+al*(exx+eyy)-beta(1)
            stre(4)=um*exy
            stre(5)=um*exz
            stre(6)=um*eyz
         elseif(mattyp.eq.2) then
            exx=eloc(1)
            eyy=eloc(2)
            ezz=eloc(3)
            exy=2.d0*eloc(4)
            exz=2.d0*eloc(5)
            eyz=2.d0*eloc(6)
            stre(1)=elas(1)*exx+elas(2)*eyy+
     &           elas(4)*ezz-beta(1)
            stre(2)=elas(2)*exx+elas(3)*eyy+
     &           elas(5)*ezz-beta(2)
            stre(3)=elas(4)*exx+elas(5)*eyy+
     &           elas(6)*ezz-beta(3)
            stre(4)=elas(7)*exy
            stre(5)=elas(8)*exz
            stre(6)=elas(9)*eyz
         elseif(mattyp.eq.3) then
            ekl(1,1)=eloc(1)
            ekl(2,2)=eloc(2)
            ekl(3,3)=eloc(3)
            ekl(1,2)=eloc(4)
            ekl(1,3)=eloc(5)
            ekl(2,3)=eloc(6)
            ekl(2,1)=eloc(4)
            ekl(3,1)=eloc(5)
            ekl(3,2)=eloc(6)
            do m1=1,3
               do m2=1,m1
                  tkl(m1,m2)=0.d0
                  do m3=1,3
                     do m4=1,3
                        tkl(m1,m2)=tkl(m1,m2)+
     &                       anisox(m1,m2,m3,m4)*ekl(m3,m4)
                     enddo
                  enddo
               enddo
            enddo
            stre(1)=tkl(1,1)-beta(1)
            stre(2)=tkl(2,2)-beta(2)
            stre(3)=tkl(3,3)-beta(3)
            stre(4)=tkl(2,1)-beta(4)
            stre(5)=tkl(3,1)-beta(5)
            stre(6)=tkl(3,2)-beta(6)
         endif
      endif
!
      return
      end
