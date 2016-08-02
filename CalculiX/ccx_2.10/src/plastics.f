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
      subroutine plastics(inpc,textpart,nelcon,nmat,ntmat_,npmat_,
     &        plicon,nplicon,plkcon,nplkcon,iplas,iperturb,nstate_,
     &        ncmat_,elcon,matname,irstrt,istep,istat,n,iline,ipol,
     &        inl,ipoinp,inp,ipoinpc)
!
!     reading the input deck: *PLASTIC
!
      implicit none
!
      logical iso
!
      character*1 inpc(*)
      character*80 matname(*)
      character*132 textpart(16)
!
      integer nelcon(2,*),nmat,ntmat_,ntmat,npmat_,npmat,istep,
     &  n,key,i,nplicon(0:ntmat_,*),nplkcon(0:ntmat_,*),ncmat_,
     &  iplas,iperturb(*),istat,nstate_,kin,itemp,ndata,ndatamax,id,
     &  irstrt,iline,ipol,inl,ipoinp(2,*),inp(3,*),ipoinpc(0:*)
!
      real*8 plicon(0:2*npmat_,ntmat_,*),plkcon(0:2*npmat_,ntmat_,*),
     & temperature,plconloc(802),t1l,elcon(0:ncmat_,ntmat_,*)
!
      iso=.true.
!
      ntmat=0
      npmat=0
!
      if((istep.gt.0).and.(irstrt.ge.0)) then
         write(*,*) '*ERROR reading *PLASTIC: *PLASTIC should be placed'
         write(*,*) '  before all step definitions'
         call exit(201)
      endif
!
      if(nmat.eq.0) then
         write(*,*) 
     &      '*ERROR reading *PLASTIC: *PLASTIC should be preceded'
         write(*,*) '  by a *MATERIAL card'
         call exit(201)
      endif
!
      if((nelcon(1,nmat).ne.2).and.(nelcon(1,nmat).ne.9)) then
         write(*,*) 
     &        '*ERROR reading *PLASTIC: *PLASTIC should be preceded'
         write(*,*) '  by an *ELASTIC,TYPE=ISO card or'
         write(*,*) '  by an *ELASTIC,TYPE=ORTHO card'
         call exit(201)
      endif
!
      iperturb(1)=3
!
      if(nelcon(1,nmat).eq.2) then
         iplas=1
         nelcon(1,nmat)=-51
         nstate_=max(nstate_,13)
      else
         write(*,*) '*ERROR reading *PLASTIC: *PLASTIC cannot '
         write(*,*) '       be used with an anisotropic material'
         write(*,*) '       use a material user subroutine instead'
         call exit(201)
      endif
!
      do i=2,n
         if(textpart(i)(1:10).eq.'HARDENING=') then
            if(textpart(i)(11:19).eq.'KINEMATIC') then
               iso=.false.
            elseif(textpart(i)(11:18).eq.'COMBINED') then
               iso=.false.
            elseif(textpart(i)(11:14).eq.'USER') then
               call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &              ipoinp,inp,ipoinpc)
               return
            endif
            exit
         else
            write(*,*) 
     &        '*WARNING reading *PLASTIC: parameter not recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*PLASTIC%")
         endif
      enddo
!
      if(iso) then
!
!        isotropic hardening coefficients
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*PLASTIC%")
!
!           first temperature
!
            if(ntmat.eq.0) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *PLASTIC: increase ntmat_'
                  call exit(201)
               endif
               nplicon(0,nmat)=ntmat
               plicon(0,ntmat,nmat)=temperature
!
!           new temperature
!
            elseif(plicon(0,ntmat,nmat).ne.temperature) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *PLASTIC: increase ntmat_'
                  call exit(201)
               endif
               nplicon(0,nmat)=ntmat
               plicon(0,ntmat,nmat)=temperature
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &              plicon(2*npmat+i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*PLASTIC%")
            enddo
            npmat=npmat+1
            if(npmat.gt.npmat_) then
               write(*,*) '*ERROR reading *PLASTIC: increase npmat_'
               call exit(201)
            endif
            nplicon(ntmat,nmat)=npmat
         enddo
      else
!
!        kinematic hardening coefficients
!
         do
            call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &           ipoinp,inp,ipoinpc)
            if((istat.lt.0).or.(key.eq.1)) exit
            read(textpart(3)(1:20),'(f20.0)',iostat=istat) temperature
            if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*PLASTIC%")
!
!           first temperature
!
            if(ntmat.eq.0) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *PLASTIC: increase ntmat_'
                  call exit(201)
               endif
               nplkcon(0,nmat)=ntmat
               plkcon(0,ntmat,nmat)=temperature
!
!           new temperature
!
            elseif(plkcon(0,ntmat,nmat).ne.temperature) then
               npmat=0
               ntmat=ntmat+1
               if(ntmat.gt.ntmat_) then
                  write(*,*) '*ERROR reading *PLASTIC: increase ntmat_'
                  call exit(201)
               endif
               nplkcon(0,nmat)=ntmat
               plkcon(0,ntmat,nmat)=temperature
            endif
            do i=1,2
               read(textpart(i)(1:20),'(f20.0)',iostat=istat) 
     &              plkcon(2*npmat+i,ntmat,nmat)
               if(istat.gt.0) call inputerror(inpc,ipoinpc,iline,
     &"*PLASTIC%")
            enddo
            npmat=npmat+1
            if(npmat.gt.npmat_) then
               write(*,*) '*ERROR reading *PLASTIC: increase npmat_'
               call exit(201)
            endif
            nplkcon(ntmat,nmat)=npmat
         enddo
      endif
!
      if(ntmat.eq.0) then
         write(*,*) 
     &       '*ERROR reading *PLASTIC: *PLASTIC card without data'
         call exit(201)
      endif
!
      return
      end

