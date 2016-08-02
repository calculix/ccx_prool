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
      subroutine designvariabless(inpc,textpart,tieset,tietol,istep,
     &       istat,n,iline,ipol,inl,ipoinp,inp,ntie,ntie_,ipoinpc)
!
!     reading the input deck: *DESIGNVARIABLES
!
      implicit none
!
      character*1 inpc(*)
      character*81 tieset(3,*)
      character*132 textpart(16)
!
      integer istep,istat,n,i,key,ipos,iline,ipol,inl,ipoinp(2,*),
     &  inp(3,*),ntie,ntie_,ipoinpc(0:*)
!
      real*8 tietol(3,*)
!
!     Check of correct position in Inputdeck
!
      if(istep.gt.0) then
         write(*,*) '*ERROR in designvariabless: *DESIGNVARIABLES'
         write(*,*) ' should be placed before all step definitions'
         call exit(201)
      endif
!
!     Check of correct number of ties
!
      ntie=ntie+1
      if(ntie.gt.ntie_) then
         write(*,*) '*ERROR in designvaribaless: increase ntie_'
         call exit(201)
      endif
!
!     Read in *DESIGNVARIABLES
!
      do i=2,n
         if(textpart(i)(1:5).eq.'NAME=') then
            read(textpart(i)(6:85),'(a80)',iostat=istat) 
     &           tieset(1,ntie)(1:80)
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &           "*DESIGNVARIABLE%")
         endif
       enddo  
!
      if(tieset(1,ntie)(1:1).eq.' ') then
         write(*,*) 
     &'*ERROR in designvaribaless: designvariables name is lacking'
         call inputerror(inpc,ipoinpc,iline,
     &"*DESIGNVARIABLE%")
         call exit(201)
      endif
!
!     Add "D" at the end of the name of the designvariable keyword
!      
      tieset(1,ntie)(81:81)='D' 
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
      if((istat.lt.0).or.(key.eq.1)) then
         write(*,*)
     &'*ERROR in designvaribaless: definition of the designvariables'
         write(*,*) '      is not complete.'
         call exit(201)
      endif
!
!     Read the name of the designvariable node set
!
      tieset(2,ntie)(1:81)=textpart(1)(1:81)
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &     ipoinp,inp,ipoinpc)
!
      return
      end



