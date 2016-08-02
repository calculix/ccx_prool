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
      subroutine sensitivitys(inpc,textpart,nmethod,iperturb,isolver,
     &  istep,istat,n,tinc,tper,tmin,tmax,idrct,iline,ipol,inl,ipoinp,
     &  inp,ithermal,cs,ics,tieset,istartset,
     &  iendset,ialset,ipompc,nodempc,coefmpc,nmpc,nmpc_,ikmpc,
     &  ilmpc,mpcfree,mcs,set,nset,labmpc,ipoinpc,iexpl,cfd,ttime,
     &  iaxial,nelcon,nmat,tincf)
!
!     reading the input deck: *SENSITIVITY
!
!     isolver=0: SPOOLES
!             2: iterative solver with diagonal scaling
!             3: iterative solver with Cholesky preconditioning
!             4: sgi solver
!             5: TAUCS
!             7: pardiso
!
!      iexpl==0:  structure:implicit, fluid:incompressible
!
      implicit none
!
      logical timereset
!
      character*1 inpc(*)
      character*20 labmpc(*),solver
      character*81 set(*),tieset(3,*)
      character*132 textpart(16)
!
      integer nmethod,iperturb,isolver,istep,istat,n,key,i,idrct,
     &  iline,ipol,inl,ipoinp(2,*),inp(3,*),ithermal,ics(*),iexpl,
     &  istartset(*),iendset(*),ialset(*),ipompc(*),nodempc(3,*),
     &  nmpc,nmpc_,ikmpc(*),ilmpc(*),mpcfree,nset,mcs,ipoinpc(0:*),
     &  cfd,iaxial,nelcon(2,*),nmat
!
      real*8 tinc,tper,tmin,tmax,cs(17,*),coefmpc(*),ttime,tincf
!
      idrct=0
      tper=1.d0
      tmin=0.d0
      tmax=0.d0
      timereset=.false.
!
      if(istep.lt.1) then
         write(*,*) '*ERROR reading *SENSITIVITY: *SENSITIVITY can
     &only be used within a STEP'     
         call exit(201)
      endif
!
!     default solver
!
      solver='                    '
      if(isolver.eq.0) then
         solver(1:7)='SPOOLES'
      elseif(isolver.eq.2) then
         solver(1:16)='ITERATIVESCALING'
      elseif(isolver.eq.3) then
         solver(1:17)='ITERATIVECHOLESKY'
      elseif(isolver.eq.4) then
         solver(1:3)='SGI'
      elseif(isolver.eq.5) then
         solver(1:5)='TAUCS'
      elseif(isolver.eq.7) then
         solver(1:7)='PARDISO'
      endif
!
      do i=2,n
         if(textpart(i)(1:7).eq.'SOLVER=') then
            read(textpart(i)(8:27),'(a20)') solver
         elseif((textpart(i)(1:6).eq.'DIRECT').and.
     &          (textpart(i)(1:9).ne.'DIRECT=NO')) then
            idrct=1
         elseif(textpart(i)(1:9).eq.'TIMERESET') then
            timereset=.true.
         elseif(textpart(i)(1:17).eq.'TOTALTIMEATSTART=') then
            read(textpart(i)(18:37),'(f20.0)',iostat=istat) ttime
         else
            write(*,*) 
     &        '*WARNING reading *SENSITIVITY: parameter not 
     &recognized:'
            write(*,*) '         ',
     &                 textpart(i)(1:index(textpart(i),' ')-1)
            call inputwarning(inpc,ipoinpc,iline,
     &"*STATIC%")
         endif
      enddo
!
      if(solver(1:7).eq.'SPOOLES') then
         isolver=0
      elseif(solver(1:16).eq.'ITERATIVESCALING') then
         isolver=2
      elseif(solver(1:17).eq.'ITERATIVECHOLESKY') then
         isolver=3
      elseif(solver(1:3).eq.'SGI') then
         isolver=4
      elseif(solver(1:5).eq.'TAUCS') then
         isolver=5
      elseif(solver(1:7).eq.'PARDISO') then
         isolver=7
      else
         write(*,*) '*WARNING reading *SENSITIVITY: unknown solver;'
         write(*,*) '         the default solver is used'
      endif
!
      nmethod=12
!
      call getnewline(inpc,textpart,istat,n,key,iline,ipol,inl,
     &        ipoinp,inp,ipoinpc)
!
!
      return
      end

