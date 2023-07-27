!f2py

! =========== Potential-specific initializations =============
subroutine init_linvib(mass,omega,Vconst,Vlin,nf,ns)
   use types
   use linvib, only : init
   integer  :: nf,ns
   real(dp) :: mass(nf), omega(nf)
   real(dp) :: Vconst(ns,ns), Vlin(nf,ns,ns)
!
!  Initialize linear vibronic potential
!
   call init(nf,ns,mass,omega,Vconst,Vlin)
end subroutine

subroutine init_bilinvib(mass,omega,Vconst,Vlin,Vbilin,nf,ns)
   use types
   use bilinvib, only : init
   integer  :: nf,ns
   real(dp) :: mass(nf), omega(nf)
   real(dp) :: Vconst(ns,ns), Vlin(nf,ns,ns), Vbilin(nf,nf,ns,ns)
!
!  Initialize bilinear vibronic potential
!
   call init(nf,ns,mass,omega,Vconst,Vlin,Vbilin)
end subroutine

subroutine init_frexc(mass,omega,Vconst,kappa,nf,nf_bath,ns)
   use types
   use frexc, only : init
   integer  :: nf,ns
   real(dp) :: mass(nf), omega(nf)
   real(dp) :: Vconst(ns,ns), kappa(nf_bath)
!
!  Initialize Frenkel-exciton potential
!
   call init(nf,ns,mass,omega,Vconst,kappa)
end subroutine

subroutine init_tully(model,mass)
   use types
   use tully, only : init
   integer :: model
   real(dp) :: mass(1)
!
!  Initialize Tully model
!
   call init(model,mass)
end subroutine

! ================ General initializations =====================

subroutine init_mash()
   use types
   use mash, only : init
!
!  Initialize MASH module.
!
   call init()
end subroutine

! ================ Useful wrappers for potentials etc. =================

subroutine get_vad(q,nf,ns,Vad,U)
   use types
   use pes, only : potad
   integer, intent(in) :: nf,ns
   real(dp), intent(in) :: q(nf)
   real(dp), intent(out) :: Vad(ns), U(ns,ns)
!
!  Wrapper for adiabatic potential
!
   call potad(q,Vad,U)
end subroutine

subroutine pot_mat(q,nf,ns,V,dVdq)
   use types
   use pes, only : pot, grad
   integer, intent(in)  :: nf,ns
   real(dp), intent(in) :: q(nf)
   real(dp), intent(out):: v(ns,ns), dVdq(nf,ns,ns)
!
!  Wrapper for diabatic potential matrix
!
   call pot(q,v)
   call grad(q,dVdq)
end subroutine

subroutine mashpot(q,qe,pe,nf,ns,v)
   use types
   use mash, only : mash_pot
   integer, intent(in)  :: nf
   real(dp), intent(in) :: q(nf),qe(ns),pe(ns)
   real(dp), intent(out):: v
!
!  Wrapper for mash potential
!
   call mash_pot(q,qe,pe,v)
end subroutine


subroutine mashgrad(q,qe,pe,nf,ns,dvdq)
   use types
   use pes, only : potad, grad_a
   use mash, only : cstate
   integer, intent(in) :: nf, ns
   real(dp), intent(in) :: q(nf), qe(ns), pe(ns)
   real(dp), intent(out) :: dvdq(nf)
!
!  Wrapper for gradient
!
   real(dp), allocatable :: Vad(:), U(:,:)
   integer :: a
   allocate(Vad(ns),U(ns,ns))
   call potad(q, Vad, U)
   call cstate(qe,pe,U,a)
   call grad_a(q, U, a, dvdq)
   deallocate(Vad, U)
end subroutine

subroutine nac(q,d,nf,ns)
   use types
   use pes, only : potad, mynac => nac
   integer, intent(in) :: nf, ns
   real(dp), intent(in) :: q(nf)
   real(dp), intent(out) :: d(nf,ns,ns)
!
!   Wrapper for nonadiabatic coupling vector
!
   real(dp), allocatable :: Vad(:), U(:,:)
   allocate(Vad(ns),U(ns,ns))
   call potad(q,Vad,U)
   call mynac(q,Vad,U,d)
   deallocate(Vad,U)
end subroutine

subroutine mashnacdir(q,cad,Vad,U,n,m,d,nf,ns)
   use types
   use pes, only : nacdir
   integer, intent(in) :: nf, ns
   real(dp), intent(in) :: q(nf), Vad(ns), U(ns,ns)
   complex(dpc), intent(in) :: cad(ns)
   integer, intent(in) :: n,m
   real(dp), intent(out) :: d(nf)
!
!   Wrapper for direction of momentum rescaling/reversal
!
   call nacdir(q,cad,Vad,U,n,m,d)
end subroutine

subroutine dia2ad(q,qe,pe,qa,pa,nf,ns)
   use types
   use pes, only : potad
   real(dp), intent(in) :: q(nf), qe(ns), pe(ns)
   real(dp), intent(out) :: qa(ns), pa(ns)
   integer, intent(in) :: nf, ns
!
!  Convert diabatic amplitudes to adiabatic
!
   real(dp), allocatable :: Vad(:), U(:,:)
   allocate(Vad(ns),U(ns,ns))
   call potad(q,Vad,U)
   qa(:) = matmul(qe, U)
   pa(:) = matmul(pe, U)
   deallocate(Vad,U)
end subroutine

subroutine mash_pops(q, qe, pe, pop, rep, typ, nf, ns)
   use types
   use mash, only : pops
   real(dp), intent(in) :: q(nf), qe(ns), pe(ns)
   integer, intent(in) :: nf, ns, typ
   character, intent(in) :: rep
   real(dp), intent(out) :: pop(ns)
!
!  Wrapper for mash.pops
!
   call pops(q, qe, pe, pop, rep, typ)
end subroutine



! =============== Main function for running a trajectory ===============
subroutine runtrj(q, p, qe, pe, qt, pt, qet, pet, Et, &
   dt, ierr, nt, nf, ns)
   use types
   use pes, only : potad, grad_a
   use mash, only : store, evolve, ham, cstate
   integer :: nt, nf, ns
   real(dp), intent(in) :: dt
   real(dp), intent(inout) :: q(nf), p(nf), qe(ns), pe(ns)
   real(dp), intent(out) :: qt(nt+1,nf), pt(nt+1,nf), &
      qet(nt+1,ns), pet(nt+1,ns), Et(nt+1)
   integer, intent(out) :: ierr
!
!  Run a trajectory
!
   real(dp), allocatable :: Vad(:), U(:,:), dvdq(:)
   ! logical :: fromfile = .false. ! Debugging option: 
   integer :: a !, atmp

   allocate(Vad(ns),U(ns,ns),dvdq(nf))

   call potad(q,Vad,U)
   call cstate(qe,pe,U,a)

   ! if (fromfile) then
   !    print*, 'Loading from state.out'
   !    open (unit=4,file='state.out')
   !    read (4,*) q
   !    read (4,*) p
   !    read (4,*) qe
   !    read (4,*) pe
   !    read (4,*) a
   !    close (unit=4)
   !    call potad(q,Vad,U)
   ! end if

   call grad_a(q,U,a,dvdq)
   do it = 1, nt
      call store(q, p, qe, pe, it, qt, pt, qet, pet)
      Et(it) = ham(q,p,qe,pe)
      ! atmp = a

      call evolve(q, p, qe, pe, Vad, U, dvdq, a, dt, ierr)

      ! if (ierr.gt.0 .and. .not.fromfile) then
      !    open (unit=4,file='state.out')
      !    write (4,*) qt(it,:)
      !    write (4,*) pt(it,:)
      !    write (4,*) qet(it,:)
      !    write (4,*) pet(it,:)
      !    write (4,*) atmp
      !    close (unit=4)
      ! end if

      if (ierr.gt.0) exit ! Don't waste time if trajectory will be discarded

      ! if (q(1).ne.q(1)) then
      !    ! NaN
      !    stop 'Trajectory gives Not A Number!'
      ! end if
      ! oldh = newh
   end do
   if (ierr.eq.0) then
      call store(q, p, qe, pe, nt+1, qt, pt, qet, pet)
      Et(nt+1) = ham(q,p,qe,pe)
   end if

   deallocate(Vad,U,dvdq)
end subroutine


! ============= Parallelized functions ==============

subroutine runpar_poponly(q, p, qe, pe, bdt, Et, ierr, &
   dt, nt, nf, ns, np)
   use types
   use mash, only : pops
   use omp_lib
   integer :: nt, nf, ns, np
   real(dp), intent(in) :: dt
   real(dp), intent(inout) :: q(nf,np), p(nf,np), qe(ns,np), pe(ns,np)
   real(dp), intent(out) :: bdt(nt+1,ns)
   real(dp),allocatable :: dbdt(:,:,:)
   real(dp),allocatable :: qt(:,:,:), pt(:,:,:), qet(:,:,:), pet(:,:,:)
   real(dp), intent(out) :: Et(nt+1)
   integer, intent(out) :: ierr(np)
   real(dp), allocatable :: dEt(:,:)
!
!  Run a batch of trajectories in parallel
!  np = number of parallel threads
!
   allocate(dbdt(nt+1,ns,np),qt(nt+1,nf,np),pt(nt+1,nf,np), &
            qet(nt+1,ns,np),pet(nt+1,ns,np),dEt(nt+1,np))
!$omp parallel do default(shared) private(j)
   do j=1,np
      call runtrj(q(:,j), p(:,j), qe(:,j), pe(:,j), qt(:,:,j), &
         pt(:,:,j), qet(:,:,j), pet(:,:,j), dEt(:,j),  &
         dt, ierr(j), nt, nf, ns)

      do it=1,nt+1
         call pops(qt(it,:,j), qet(it,:,j), pet(it,:,j), &
            dbdt(it,:,j), 'd', 2)
      end do

      if (ierr(j).ne.0) then
         ! Error -> discard trajectory
         dbdt(:,:,j) = 0.d0
      end if
   end do
!$omp end parallel do
   bdt = sum(dbdt,dim=3)
   Et  = sum(dEt,dim=2)
   deallocate(dbdt,dEt,qt,pt,qet,pet)
end subroutine

subroutine runpar_all(q, p, qe, pe, bdt, Et, ierr, &
   dt, nt, nf, ns, np)
   use types
   use mash, only : obsbls
   use omp_lib
   integer :: nt, nf, ns, np
   real(dp), intent(in) :: dt
   real(dp), intent(inout) :: q(nf,np), p(nf,np), qe(ns,np), pe(ns,np)
   complex(dpc), intent(out) :: bdt(nt+1,ns,ns)
   real(dp), intent(out) :: Et(nt+1)
   integer, intent(out) :: ierr(np)
   complex(dpc),allocatable :: dbt(:,:,:,:)
   real(dp), allocatable :: qt(:,:,:), pt(:,:,:), qet(:,:,:), pet(:,:,:)
   real(dp), allocatable :: dEt(:,:)
!
!  Run a batch of trajectories in parallel
!  np = number of parallel threads
!
   allocate(dbt(nt+1,ns,ns,np),qt(nt+1,nf,np),pt(nt+1,nf,np), &
   qet(nt+1,ns,np),pet(nt+1,ns,np),dEt(nt+1,np))
   allocate(dEt(nt+1,np))
!$omp parallel do default(shared) private(j)
   do j=1,np
      call runtrj(q(:,j), p(:,j), qe(:,j), pe(:,j), qt(:,:,j), &
         pt(:,:,j), qet(:,:,j), pet(:,:,j), dEt(:,j),  &
         dt, ierr(j), nt, nf, ns)
      do it=1,nt+1
         call obsbls(qt(it,:,j), qet(it,:,j), pet(it,:,j), & 
            dbt(it,:,:,j), 'd', 2, .false.)
      end do
      if (ierr(j).ne.0) then
         ! Error -> discard trajectory
         dbt(:,:,:,j) = 0.d0
      end if
   end do
!$omp end parallel do
   bdt = sum(dbt,dim=4)
   Et  = sum(dEt,dim=2)
   deallocate(dbt,qt,pt,qet,pet,dEt)
end subroutine
