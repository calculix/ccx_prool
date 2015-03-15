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
      subroutine controlss(inpc,textpart,ctrl,istep,istat,n,iline,ipol,
     &  inl,ipoinp,inp,ipoinpc)
!
!     reading the input deck: *STEP
!
      implicit none
!
      character*1 inpc(*)
      character*132 textpart(16)
!
      integer i,j,k,istep,istat,n,key,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ipoinpc(0:*)
!
      real*8 ctrl(*)
!
      do i=2,n
         if(textpart(i)(1:5).eq.'RESET') then
            ctrl(1)=4.5d0
            ctrl(2)=8.5d0
            ctrl(3)=9.5d0
            ctrl(4)=16.5d0
            ctrl(5)=10.5d0
            ctrl(6)=4.5d0
            ctrl(7)=0.
            ctrl(8)=5.5d0
            ctrl(9)=0.
            ctrl(10)=0.
            ctrl(11)=0.25d0
            ctrl(12)=0.5d0
            ctrl(13)=0.75d0
            ctrl(14)=0.85d0
            ctrl(15)=0.
            ctrl(16)=0.
            ctrl(17)=1.5d0
            ctrl(18)=0.
            ctrl(19)=0.005d0
            ctrl(20)=0.01d0
            ctrl(21)=0.d0
            ctrl(22)=0.d0
            ctrl(23)=0.02d0
            ctrl(24)=1.d-5
            ctrl(25)=1.d-3
            ctrl(26)=1.d-8
            write(*,*)
            write(*,*) 
     &         '*INFO: control parameters reset to default'
            exit
!            
         elseif(textpart(i)(1:29).eq.'PARAMETERS=TIMEINCREMENTATION') 
     &      then
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            do j=1,n
               read(textpart(j)(1:10),'(i10)',iostat=istat) k
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*CONTROLS%")
               ctrl(j)=dble(k)+0.5d0
            enddo
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            do j=1,n
               read(textpart(j)(1:20),'(f20.0)',iostat=istat) ctrl(j+10)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*CONTROLS%")
            enddo
            write(*,*) '*INFO: time control parameters set to:'
            write(*,*) '       i0 = ',int(ctrl(1))
            write(*,*) '       ir = ',int(ctrl(2))
            write(*,*) '       ip = ',int(ctrl(3))
            write(*,*) '       ic = ',int(ctrl(4))
            write(*,*) '       il = ',int(ctrl(5))
            write(*,*) '       ig = ',int(ctrl(6))
            write(*,*) '       is = ',int(ctrl(7))
            write(*,*) '       ia = ',int(ctrl(8))
            write(*,*) '       ij = ',int(ctrl(9))
            write(*,*) '       it = ',int(ctrl(10))
            write(*,*) '       df = ',ctrl(11)
            write(*,*) '       dc = ',ctrl(12)
            write(*,*) '       db = ',ctrl(13)
            write(*,*) '       da = ',ctrl(14)
            write(*,*) '       ds = ',ctrl(15)
            write(*,*) '       dh = ',ctrl(16)
            write(*,*) '       dd = ',ctrl(17)
            write(*,*) '       wg = ',ctrl(18)
            exit
!
         elseif(textpart(i)(1:16).eq.'PARAMETERS=FIELD') then
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            do j=1,n
               read(textpart(j)(1:20),'(f20.0)',iostat=istat) ctrl(j+18)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*CONTROLS%")
            enddo
            write(*,*) '*INFO: field control parameters set to:'
            write(*,*) '       ran = ',ctrl(19)
            write(*,*) '       can = ',ctrl(20)
            write(*,*) '       qa0 = ',ctrl(21)
            write(*,*) '       qau = ',ctrl(22)
            write(*,*) '       rap = ',ctrl(23)
            write(*,*) '        ea = ',ctrl(24)
            write(*,*) '       cae = ',ctrl(25)
            write(*,*) '       ral = ',ctrl(26)
            exit
!
         elseif(textpart(i)(1:21).eq.'PARAMETERS=LINESEARCH') then
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            do j=1,n
               read(textpart(j)(1:20),'(f20.0)',iostat=istat) ctrl(j+27)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*CONTROLS%")
            enddo
            write(*,*) '*INFO: field control parameters set to:'
            write(*,*) '       nls = ',ctrl(28)
            write(*,*) '       smaxls = ',ctrl(29)
            write(*,*) '       sminls = ',ctrl(30)
            write(*,*) '       fls = ',ctrl(31)
            write(*,*) '       etls = ',ctrl(32)
            exit
         else
            write(*,*) 
     &        '*WARNING in controlss: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*CONTROLS%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end








