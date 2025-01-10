module mash

   ! Module for running mapping approach to surface hopping

   use types
   implicit double precision (a-h,o-z)

   real(dp) :: alpha   ! ns-dependent constant used in observables
   real(dp) :: beta    ! Reciprocal temperature

contains

! =============== Initialization =============
   subroutine init(beta_)
      use pes, only : ns
      real(dp) :: beta_
!
!    Initialize module
!
      Hn = 0.d0
      do n=1,ns
         Hn = Hn + 1.d0/n
      end do
      alpha = (ns-1.d0)/(Hn-1.d0)

      beta = beta_
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
   
   real(dp) function ham_a(q, p, a)
      use pes, only : ns, mass, potad
      real(dp), intent(in) :: q(:), p(:)
      integer :: a
!
!     Function to compute Hamiltonian at a phase-space point
!
      real(dp), allocatable :: U(:,:), Vad(:)
      allocate(U(ns,ns),Vad(ns))
      call potad(q,Vad,U)
      ham_a = Vad(a) + 0.5d0*sum(p**2/mass)
      deallocate(U,Vad)
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

   subroutine pops_phi(c, pop)
      use pes, only : ns
      complex(dpc), intent(inout) :: c(:)
      real(dp), intent(inout) :: pop(:)
!
!  Calculate Phi_n observable (|n><m| becomes alpha_N c_n^*c_m + beta_N delta_{nm})
!
      overn = 1.d0/ns
      pop = overn + alpha*(abs(c)**2 - overn)
   end subroutine

   subroutine pops(q, qe, pe, pop, rep)
      use pes, only : ns, potad
      real(dp), intent(in) :: q(:), qe(:), pe(:)
      real(dp), intent(out) :: pop(:)
      character, intent(in) :: rep ! representation ('d' for diabatic/site, 'e' for exciton, 'a' for adiabatic)
!
!  Observables in diabatic or adiabatic representation
!
      real(dp), allocatable :: Vad(:), U(:,:)
      complex(dpc), allocatable :: c(:)
      if (.not. (rep.eq.'d' .or. rep.eq.'e' .or. rep.eq.'a')) then
         stop 'obsbls: Undefined representation'
      end if
      allocate(c(ns))
      if (rep.eq.'d') then
         ! Diabatic (original basis) observables
         c = dcmplx(qe,pe)
         call pops_phi(c, pop)
      else if (rep.eq.'e') then
         ! Diabatic (exciton) observables
         allocate(Vad(ns),U(ns,ns))
         call potad(q*0.d0,Vad,U)
         c = dcmplx(matmul(qe,U),matmul(pe,U))
         call pops_phi(c, pop)
         deallocate(Vad,U)
      else if (rep.eq.'a') then
         ! Adiabatic observables
         allocate(Vad(ns),U(ns,ns))
         call potad(q,Vad,U)
         c = dcmplx(matmul(qe,U),matmul(pe,U))
         call pops_phi(c, pop)
         deallocate(Vad,U)
      end if
      deallocate(c)
   end subroutine


! =============== Dynamics-related subroutines ===============
   subroutine evolve(q, p, qe, pe, Vad, U, dvdq, a, dtbase)
      use pes, only : nf, ns, grad_a
      real(dp), intent(inout) :: q(:), p(:), qe(:), pe(:), Vad(:), &
                                 U(:,:), dvdq(:)
      integer, intent(inout) :: a
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

      maxhop = 10
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

         if (b.eq.a) then
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
            ! Too many hops -- just finish it off
            call verlet(q,p,qe,pe,Vad,U,dvdq,a,dt)
            ca1 = dcmplx(matmul(qe,U),matmul(pe,U))
            call cstate_ad(ca1,b)
            if (b.ne.a) then
               call cross(q, p, ca1, a, b, Vad, U, accepted)
               if (accepted) then
                  call grad_a(q,U,b,dvdq)
                  a = b
               end if
            end if                 
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
      use pes, only : nf,mass,omega,cayley
      real(dp), intent(inout) :: q(nf), p(nf)
      real(dp), intent(in) :: dt
!
!     Cayley evolution through a time interval dt.
!
      if (cayley) then
         do i = 1,nf
            em = mass(i)
            dtm = dt/em
            w2 = omega(i)**2
            x = w2*dt**2/4
            pp = (1-x)/(1+x)
            pq = -em*w2*dt/(1+x)
            qp = dtm/(1+x)
            qq = (1-x)/(1+x)
            
            pnew = pq*q(i) + pp*p(i)  
            q(i) = qq*q(i) + qp*p(i)
            p(i) = pnew
         end do
      else
         q = q + dt*p/mass
      end if
   end subroutine

   subroutine step_e(qe, pe, Vad, U, dt)
      use pes, only : ns
      use maths, only : iu, symevp
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
      Ekin = 0.5d0*sum(pnac**2)
      Vdiff = Vad(m)-Vad(n)
      if ((Ekin-Vdiff).gt.0.d0) then
         ! Rescale momentum 
         pnac = sqrt(2.d0*(Ekin-Vdiff)) & ! (no masses since p is mass-scaled)
             * pnac/sqrt(dot_product(pnac,pnac))
         accepted = .true.
      else
         ! Reverse momentum along NAC vector
         pnac = -pnac
         accepted = .false.
      end if      
      p = porth + pnac
      p = p*sqrt(mass)
      deallocate(d,dj,pnac,porth)
   end subroutine

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

   subroutine store(q, p, qe, pe, a, it, qt, pt, qet, pet, at)
      use types
      real(dp), intent(in) :: q(:), p(:), qe(:), pe(:)
      integer, intent(in) :: a, it
      real(dp), intent(inout) :: qt(:,:), pt(:,:), qet(:,:), pet(:,:)
      integer, intent(inout) :: at(:)
!
!  Store state at a given timestep
!
      qt(it,:) = q
      pt(it,:) = p
      qet(it,:) = qe
      pet(it,:) = pe
      at(it) = a
   end subroutine 

end module
