! Frenkel-exciton potential with sites coupled to identical independent baths
! This contains faster routines for computing adiabatic gradients than the linvib module

module frexc

   use types

   real(dp), allocatable :: mw2(:), kappa(:,:), Vconst(:,:)
   integer :: nf_bath

contains

   subroutine init(nf_,ns_,mass_,omega_,Vconst_,kappa_)
      use pes, only : pesinit=>init, potptr=>pot, gradptr=>grad, nf, ns, omega, cayley
      use pes, only : grad_a,grad_ab,pes_get_vconst=>get_vconst,pes_get_vlin=>get_vlin
      use pes, only : grad_av
      integer :: nf_, ns_
      real(dp) :: mass_(nf_), omega_(nf_)
      real(dp) :: Vconst_(ns_,ns_), kappa_(:,:)
!
!     Initialize module
!     omega_ should contain nf elements, consisting of ns identical blocks of nf_bath elements each
!     kappa_ should have shape (nf_bath,ns)
!
      call pesinit(nf_,ns_,mass_)
      potptr => pot
      gradptr => grad
      pes_get_vconst => get_vconst
      pes_get_vlin => get_vlin
      nf_bath = nf/ns
      if (ns*nf_bath.ne.nf) stop 'frexc.f90: nf has do be divisible by ns'
      if (.not.allocated(mw2)) allocate(mw2(nf_bath),Vconst(ns,ns),kappa(nf_bath,ns),omega(nf))
      mw2 = mass_(:nf_bath)*omega_(:nf_bath)**2
      omega = omega_
      kappa = kappa_
      Vconst = Vconst_
      ! Override generic adiabatic gradients with faster ones
      grad_a => gradad_diag
      grad_ab => gradad
      grad_av => grad_vector
      cayley = .true.
   end subroutine


   subroutine pot(q, V)
      use pes, only : ns,nf
      real(dp), intent(in) :: q(nf)
      real(dp), intent(out) :: V(ns,ns)
!
!     Diabatic potential matrix
!
      V = Vconst
      call addbath(q(1), V)
    end subroutine

   subroutine addbath(q, V)
      use pes, only : ns
      real(dp), intent(in) :: q(nf_bath,ns)
      real(dp), intent(inout) :: V(ns,ns)
!
!     Potential for Frenkel-exciton model. 
!     Reshaping of q is done automatically when called with shape q(nf)
!
      real(dp) :: V0
      real(dp), allocatable :: q2(:,:)
      allocate(q2(nf_bath,ns))
      q2 = q**2
      V0 = sum(matmul(mw2,q2))/2
      do n=1,ns
         V(n,n) = V(n,n) + V0 + sum(kappa(:,n)*q(:,n))
      end do
      deallocate(q2)
   end subroutine

   subroutine grad(q, G)
      use pes, only : nf,ns
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: G(nf,ns,ns)
!
!     Analytic gradient of diabatic potential
!
      call grad_all(q, G)
   end subroutine

   subroutine grad_all(q, G)
      use pes, only : ns, cayley
      real(dp), intent(in) :: q(nf_bath,ns)
      real(dp), intent(out) :: G(nf_bath,ns,ns,ns)
!
!     Auxiliary function for grad to help reshaping
!
      G = 0.d0
      do n=1,ns
         do m=1,ns
            if (.not.cayley) G(:,n,m,m) = mw2*q(:,n)
         end do
         G(:,n,n,n) = G(:,n,n,n) + kappa(:,n)
      end do
   end subroutine

   subroutine gradad_diag(q,U,a,dvdq)
      use pes, only : nf
      integer :: a
      real(dp), intent(in) :: q(:), U(:,:)
      real(dp), intent(out) :: dvdq(nf)
      call gradad_diag_reshaped(q,U,a,dvdq)
   end subroutine

   subroutine gradad_diag_reshaped(q,U,a,dvdq)
      use pes, only : ns, cayley
      integer :: a
      real(dp), intent(in) :: q(nf_bath,ns), U(:,:)
      real(dp), intent(out) :: dvdq(nf_bath,ns)
!
!     Compute gradient of adiabatic potential a for special case of a Frenkel-exciton model
!     U -- transformation matrix from adia to dia basis
!     dvdq -- gradient of adiabatic state a, automatically reshaped to dvdq(nf)
!
      if (cayley) then
         do n=1,ns
            dvdq(:,n) = U(n,a)**2 * kappa(:,n)
         end do
      else
         do n=1,ns
            dvdq(:,n) = mw2*q(:,n) + U(n,a)**2 * kappa(:,n)
         end do
      end if
   end subroutine

   subroutine gradad(q,U,Gad)
      use pes, only : ns,nf
      real(dp), intent(in) :: q(:), U(:,:)
      real(dp), intent(out) :: Gad(nf,ns,ns)
      call gradad_reshaped(q,U,Gad)
   end subroutine

   subroutine gradad_reshaped(q,U,Gad)
      use pes, only : ns
      real(dp), intent(in) :: q(nf_bath,ns), U(:,:)
      real(dp), intent(out) :: Gad(nf_bath,ns,ns,ns)
!
!     Compute full gradient of adiabatic potential for special case of a Frenkel-exciton model
!     U -- transformation matrix from adia to dia basis
!     G -- gradient of adiabatic state a, automatically reshaped to G(nf,ns,ns)
!     NOTE: includes gradient of harmonic term no matter if cayley is true or false
!           because this function is used for hops in both cases
!
      integer :: a,b
      Gad = 0.d0
      do n=1,ns
         do a=1,ns
            Gad(:,n,a,a) = mw2*q(:,n)
            do b=1,ns
               Gad(:,n,a,b) = Gad(:,n,a,b) + U(n,a)*U(n,b)*kappa(:,n)
            end do
         end do
      end do
   end subroutine

   subroutine grad_vector(q,U,Gad,b)
      use pes, only : nf,ns
      integer :: b
      real(dp), intent(in) :: q(:), U(:,:)
      real(dp), intent(out) :: Gad(nf,ns)
      call grad_vector_reshaped(q,U,Gad,b)
   end subroutine

   subroutine grad_vector_reshaped(q,U,Gad,b)
      use pes, only : ns
      integer :: a,b
      real(dp), intent(in) :: q(nf_bath,ns), U(:,:)
      real(dp), intent(out) :: Gad(nf_bath,ns,ns)
!
!     Compute gradient of adiabatic potential for special case of a Frenkel-exciton model
!     U -- transformation matrix from adia to dia basis
!     Gad -- gradient of adiabatic state a, automatically reshaped to Gad(nf,ns)
!
      Gad = 0.d0
      do n=1,ns
         Gad(:,n,b) = mw2*q(:,n)
         do a=1,ns
            Gad(:,n,a) = Gad(:,n,a) + U(n,a)*U(n,b)*kappa(:,n)
         end do
      end do
   end subroutine

   subroutine get_vconst(Vconst_)
      use pes, only : ns
      real(dp), intent(out) :: Vconst_(ns,ns)
      Vconst_ = Vconst
   end subroutine

   subroutine get_vlin(vlin)
      use pes, only : ns, nf
      real(dp), intent(out) :: vlin(nf,ns,ns)
      real(dp), allocatable :: vlin4d(:,:,:,:)
      allocate(vlin4d(nf_bath,ns,ns,ns))
      vlin4d = 0.d0
      do n=1,ns
         vlin4d(:,n,n,n) = kappa(:,n)
      end do
      vlin = reshape(vlin4d,(/nf,ns,ns/))
   end subroutine

end module
