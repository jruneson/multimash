module mash

   ! Module for running mapping approach to surface hopping

   use types
   implicit double precision (a-h,o-z)

   real(dp) :: alpha  ! ns-dependent constant in observables

   real(dp), allocatable :: Umeas(:,:)

contains

! =============== Initialization =============
   subroutine init()
      use pes, only : ns
!
!    Initialize module
!
      Hn = 0.d0
      do n=1,ns
         Hn = Hn + 1.d0/n
      end do
      alpha = (ns-1.d0)/(Hn-1.d0)
   end subroutine

   subroutine initbasis(Umeas_)
      use pes, only : ns
      real(dp) :: Umeas_(ns,ns)
!
!     Initialize basis to do measurements in
!
      allocate(Umeas(ns,ns))
      Umeas = Umeas_
   end subroutine

! ============= Potential, Hamiltonian =============

   subroutine mash_pot(q, qe, pe, vtot)
      use pes, only : nf, ns, potad
      real(dp), intent(in) :: q(nf), qe(ns), pe(ns)
      real(dp), intent(out) :: vtot
!
!     Compute MASH potential at (q,qe,pe)
!
      integer :: a
      real(dp), allocatable :: U(:,:), Vad(:)
      allocate(U(ns,ns),Vad(ns))
      call potad(q,Vad,U)
      call cstate(qe,pe,U,a)
      vtot = Vad(a)
      deallocate(U,Vad)
   end subroutine


   real(dp) function ham(q, p, qe, pe)
      use pes, only : nf, ns, mass
      real(dp), intent(in) :: q(nf), p(nf), qe(ns), pe(ns)
!
!     Function to compute Hamiltonian at a phase-space point
!
      call mash_pot(q, qe, pe, vtot)
      ham = vtot + 0.5d0*sum(p**2/mass)
   end function


! =============== Observable-related subroutines =============

   subroutine cstate(qe, pe, U, a)
      use pes, only : ns
      real(dp), intent(in) :: qe(:), pe(:), U(:,:)
      integer, intent(out) :: a
!
!     Get current adiabatic state with the highest population |c_a|^2
!
      real(dp), allocatable :: qa(:), pa(:)
      complex(dpc), allocatable :: ca(:)
      allocate(qa(ns),pa(ns),ca(ns))
      ! Convert to adiabatic rep.
      qa = matmul(qe,U)
      pa = matmul(pe,U)
      ca = dcmplx(qa,pa)
      call cstate_ad(ca,a)
      deallocate(qa,pa,ca)
   end subroutine

   subroutine cstate_ad(ca, a)
      use pes, only : ns
      complex(dpc), intent(in) :: ca(:)
      integer, intent(out) :: a
!
!     Get current adiabatic state with the highest population |c_a|^2
!
      integer :: b
      real(dp), allocatable :: pop(:)
      allocate(pop(ns))
      ! Convert to adiabatic rep.
      pop = abs(ca)**2
      a = 1
      do b = 2,ns
         if (pop(b)>pop(a)) a=b
      end do
      deallocate(pop)
   end subroutine

   subroutine cstate2_ad(ca, a, b)
      use pes, only : ns
      complex(dpc), intent(in) :: ca(:)
      integer :: a, b
! 
!     Get state b with second-highest population excluding a
!
      real(dp), allocatable :: pop(:)
      allocate(pop(ns))
      pop = abs(ca)**2
      pop(a) = 0.d0
      b = 1
      do n = 2,ns
         if (pop(n)>pop(b)) b=n
      end do
      deallocate(pop)
   end subroutine

   subroutine obsbls_mash(c, obs, typ, poponly)
      use pes, only : ns
      complex(dpc), intent(inout) :: c(:)
      complex(dpc), intent(inout) :: obs(:,:)
      integer :: typ
      logical :: poponly
!
!  Calculate binned observables in given representation. 
!     obs(k,k) -- population
!     obs(k,l) -- coherence 
!
!     typ==1 -- binned (Theta function)
!     typ==2 -- weighted in the Phin way
!
      overn = 1.d0/ns
      call cstate_ad(c,l) !highest populated state
      obs = 0.d0
      if (typ.eq.1) then
         obs(l,l) = 1.d0
      else if (typ.eq.2) then
         do k=1,ns
            obs(k,k) = overn + alpha*(abs(c(k))**2 - overn)
         end do
      end if
      if (poponly) then
         return
      end if

      ! Coherences
      do k=1,ns
         do l=k+1,ns
            obs(k,l) = alpha*conjg(c(k))*c(l)
            obs(l,k) = conjg(obs(k,l))
         end do
      end do
   end subroutine

   subroutine obsbls(q, qe, pe, obs, rep, typ, poponly)
      use pes, only : ns, potad
      real(dp), intent(in) :: q(:), qe(:), pe(:)
      complex(dpc), intent(out) :: obs(:,:)
      character, intent(in) :: rep ! representation ('d' or 'a')
      integer, intent(in) :: typ ! projection type theta_n (1) or Phi_n (2)
      logical, intent(in) :: poponly ! Only compute populations
!
!  Observables in diabatic or adiabatic representation
!
      real(dp), allocatable :: Vad(:), U(:,:)
      complex(dpc), allocatable :: c(:)
      allocate(c(ns))
      if (rep.eq.'a') then
         ! Adiabatic observables
         allocate(Vad(ns),U(ns,ns))
         call potad(q,Vad,U)
         c = dcmplx(matmul(qe,U),matmul(pe,U))
         call obsbls_mash(c, obs, typ, poponly)
         deallocate(Vad,U)
      else if (rep.eq.'d') then
         ! Diabatic observables
         c = dcmplx(matmul(qe,Umeas),matmul(pe,Umeas))
         call obsbls_mash(c, obs, typ, poponly)
      else
         stop 'obsbls: Undefined representation'
      end if
      deallocate(c)
   end subroutine

   subroutine pops(q, qe, pe, pop, rep, typ)
      use pes, only : ns
      real(dp), intent(in) :: q(:), qe(:), pe(:)
      real(dp), intent(inout) :: pop(:)
      character, intent(in) :: rep
      integer, intent(in) :: typ
!
!  All populations |n><n|
!
      complex(dpc), allocatable :: obs(:,:)
      allocate(obs(ns,ns))
      call obsbls(q, qe, pe, obs, rep, typ, .true.)
      do n=1,ns
         pop(n) = real(obs(n,n))
      end do
      deallocate(obs)
   end subroutine

! =============== Dynamics-related subroutines ===============
   subroutine evolve(q, p, qe, pe, Vad, U, dvdq, a, dtbase, ierr)
      use pes, only : nf, ns, grad_a
      real(dp), intent(inout) :: q(:), p(:), qe(:), pe(:), Vad(:), &
                                 U(:,:), dvdq(:)
      integer, intent(inout) :: a, ierr
      real(dp), intent(in) :: dtbase
      real(dp), allocatable :: q0(:),p0(:),qe0(:),pe0(:), &
                               Vad0(:), U0(:,:), dvdq0(:)
      complex(dpc), allocatable :: ca0(:), ca1(:)
      integer :: b
      logical :: accepted
!
!     Perform a time step
!
      allocate(q0(nf),p0(nf),qe0(ns),pe0(ns),Vad0(ns),U0(ns,ns),dvdq0(nf))
      allocate(ca0(ns),ca1(ns))

      ierr = 0
      maxhop = 30
      dt = dtbase
      do ihop=1,maxhop ! Limit number of hops to look for in a timestep dt
         ! Store initial values
         call savetmp(q,p,qe,pe,Vad,U,dvdq,q0,p0,qe0,pe0,Vad0,U0,dvdq0)
      
         ! Store initial adiabatic wavefunction
         ca0 = dcmplx(matmul(qe,U),matmul(pe,U))
         
         ! Attempt full step
         call verlet(q,p,qe,pe,Vad,U,dvdq,a,dt)

         ! Calculate new active state
         ca1 = dcmplx(matmul(qe,U),matmul(pe,U))
         call cstate_ad(ca1,b)

         if (a.eq.b) then
            ! Stayed on state - we're done
            exit
         else
            ! States have changed - find crossing time with bisection root search
            tl = 0.0
            tr = dt
            fl = deltaP(ca0,a)
            fr = deltaP(ca1,a)
            do iter=1,10
               tm = (tl+tr)/2 ! Mid point
               call savetmp(q0,p0,qe0,pe0,Vad0,U0,dvdq0,q,p,qe,pe,Vad,U,dvdq)
               call verlet(q,p,qe,pe,Vad,U,dvdq,a,tm)
               ca1 = dcmplx(matmul(qe,U),matmul(pe,U))
               fm = deltaP(ca1,a)
               if (fm.gt.0) then
                  tl = tm
               else
                  tr = tm
               end if
               ! print*, iter, tm, dt, fm
            end do
            call cstate2_ad(ca1,a,b)
            call cross(q, p, ca1, a, b, Vad, U, accepted)
            if (accepted) then
               tx = tr
            else
               tx = tl
            end if
            call savetmp(q0,p0,qe0,pe0,Vad0,U0,dvdq0,q,p,qe,pe,Vad,U,dvdq)
            call verlet(q,p,qe,pe,Vad,U,dvdq,a,tx)
            ca1 = dcmplx(matmul(qe,U),matmul(pe,U))
            call cross(q, p, ca1, a, b, Vad, U, accepted)
            if (accepted) then
               call grad_a(q,U,b,dvdq)
               a = b
            end if
            dt = dt - tx
         end if
         if (ihop.eq.maxhop) then
            ! Too many hops -- highlight trajectory
            ierr=1
         end if
      end do
      deallocate(ca0,ca1)
      deallocate(q0,p0,qe0,pe0,Vad0,U0)
   end subroutine

   real(dp) function deltaP(ca,a)
      complex(dpc) :: ca(:)
      integer :: a,b
!
!     Difference in population between state a and highest populated state apart from a
!
      call cstate2_ad(ca,a,b)
      deltaP = abs(ca(a))**2 - abs(ca(b))**2
   end function

   subroutine savetmp(q,p,qe,pe,Vad,U,dvdq,q0,p0,qe0,pe0,Vad0,U0,dvdq0)
      real(dp), intent(inout) :: q(:),p(:),qe(:),pe(:),Vad(:),U(:,:),dvdq(:),&
            q0(:),p0(:),qe0(:),pe0(:), Vad0(:),U0(:,:),dvdq0(:)   
!
!     Store temporary variables and potential information
!
      q0 = q
      p0 = p
      qe0 = qe
      pe0 = pe
      Vad0 = Vad
      U0 = U
      dvdq0 = dvdq
   end subroutine


   subroutine verlet(q, p, qe, pe, Vad, U, dvdq, a, dt)
      use pes, only : potad, grad_a
      real(dp), intent(inout) :: q(:), p(:), qe(:), pe(:), Vad(:), &
                                 U(:,:), dvdq(:)
      integer, intent(in) :: a
      real(dp), intent(in) :: dt
!
!     Perform a time step on state n
!
      dt2 = dt/2
      call step_e(qe,pe,Vad,U,dt2)
      call step_p(p,dvdq,dt2)
      call step_q(q,p,dt)
      call potad(q, Vad, U)
      call grad_a(q,U,a,dvdq)
      call step_p(p,dvdq,dt2)
      call step_e(qe,pe,Vad,U,dt2)
   end subroutine


   subroutine step_p(p,dvdq,dt)
      real(dp), intent(inout) :: p(:)
      real(dp), intent(in) :: dvdq(:), dt
!
!    Evolve nuclear momenta
!
      p = p - dvdq * dt
   end subroutine

   subroutine step_q(q, p, dt)
      use pes, only : mass
      real(dp), intent(inout) :: q(:), p(:)
      real(dp), intent(in) :: dt
!
!     Evolve q and check for switch of adiabatic potential
!
      q = q + dt*p/mass
   end subroutine

   subroutine step_e(qe, pe, Vad, U, dt)
      use pes, only : ns
      real(dp), intent(inout) :: qe(:), pe(:)
      real(dp), intent(in) :: Vad(:), U(:,:)
      real(dp), intent(in) :: dt
!
!     Evolve electronic coefficients c = qe + i*pe
!
      complex(dpc), allocatable :: c(:)
      allocate(c(ns))
      ! U transforms from adiabatic to diabatic basis
      qe = matmul(qe, U) ! Same as qe = transpose(U).qe
      pe = matmul(pe, U)
      c = dcmplx(qe, pe)
      c = exp(- iu * dt * Vad) * c
      qe = matmul(U, real(c)) ! qe = U.qe
      pe = matmul(U, aimag(c))
      deallocate(c)
   end subroutine


   subroutine cross(q,p,cad,n,m,Vad,U,accepted)
      use pes, only : nf, ns, mass, nacdir
      real(dp), intent(inout) :: p(:)
      complex(dpc), intent(inout) :: cad(:)
      real(dp), intent(in) :: q(:), Vad(:), U(:,:)
      logical, intent(out) :: accepted
!
!     If energetically allowed, hop from n to m and rescale momentum
!     If not, reverse propagation direction along given component of the NAC vector
!
      real(dp), allocatable :: d(:,:,:), dj(:), pnac(:), porth(:)
      allocate(d(nf,ns,ns),dj(nf),pnac(nf),porth(nf))
      call nacdir(q,cad,Vad,U,n,m,dj)
      dj = dj/sqrt(mass)
      
      ! Use mass-scaled momenta
      p = p/sqrt(mass)
      if (nf.eq.1) then
         pnac = p
      else
         pnac = dot_product(p,dj)/dot_product(dj,dj)*dj
      end if
      porth = p - pnac
      Ekin = 0.5d0*sum(pnac**2) !/mass)
      Vdiff = Vad(m)-Vad(n)
      if ((Ekin-Vdiff).gt.0.d0) then
         ! Rescale momentum 
         pnac = sqrt(2.d0*(Ekin-Vdiff)) & ! (no masses since p is mass-scaled)
             * pnac/sqrt(dot_product(pnac,pnac))
         p = porth + pnac
         accepted = .true.
      else
         ! Reverse momentum along NAC vector
         pnac = -pnac
         p = porth + pnac
         accepted = .false.
      end if
      p = p*sqrt(mass)
      deallocate(d,dj,pnac,porth)
   end subroutine

   subroutine store(q, p, qe, pe, it, qt, pt, qet, pet)
      use types
      real(dp), intent(in) :: q(:), p(:), qe(:), pe(:)
      integer, intent(in) :: it
      real(dp), intent(inout) :: qt(:,:), pt(:,:), qet(:,:), pet(:,:)
!
!  Store state at a given timestep
!
      qt(it,:) = q
      pt(it,:) = p
      qet(it,:) = qe
      pet(it,:) = pe
   end subroutine 


end module
