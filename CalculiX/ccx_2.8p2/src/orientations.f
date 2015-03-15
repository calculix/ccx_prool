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
      subroutine orientations(inpc,textpart,orname,orab,norien,
     &  norien_,istep,istat,n,iline,ipol,inl,ipoinp,inp,ipoinpc)
!
!     reading the input deck: *ORIENTATION
!
      implicit none
!
      character*1 inpc(*)
      character*80 orname(*)
      character*132 textpart(16)
!
      integer norien,norien_,istep,istat,n,key,i,iline,ipol,inl,
     &  ipoinp(2,*),inp(3,*),ipoinpc(0:*)
!
      real*8 orab(7,*)
!
      if(istep.gt.0) then
         write(*,*) '*ERROR in orientations: *ORIENTATION should be'
         write(*,*) '  placed before all step definitions'
         call exit(201)
      endif
!
      norien=norien+1
      if(norien.gt.norien_) then
         write(*,*) '*ERROR in orientations: increase norien_'
         call exit(201)
      endif
!
!     rectangular coordinate system: orab(7,norien)=1
!     cylindrical coordinate system: orab(7,norien)=-1
!     default is rectangular
!
      orab(7,norien)=1.d0
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            orname(norien)=textpart(i)(6:85)
            if(textpart(i)(86:86).ne.' ') then
               write(*,*) '*ERROR in orientations: name too long'
               write(*,*) '       (more than 80 characters)'
               write(*,*) '       orientation name:',textpart(i)(1:132)
               call exit(201)
            endif
         elseif(textpart(i)(1:7).eq.'SYSTEM=') then
            if(textpart(i)(8:8).eq.'C') then
               orab(7,norien)=-1.d0
            endif
         else
            write(*,*) 
     &        '*WARNING in orientations: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*ORIENTATION%")
         endif
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)'*ERROR in orientations: definition of the following'
         write(*,*) '  orientation is not complete: ',orname(norien)
         call exit(201)
      endif
!
      do i=1,6
         read(textpart(i)(1:20),'(f20.0)',iostat=istat) orab(i,norien)
         if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*ORIENTATION%")
      enddo
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end

