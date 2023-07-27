! Frenkel-exciton potential with sites coupled to identical independent baths
! This contains faster routines for computing adiabatic gradients than the linvib module

module frexc

   use types

   real(dp), allocatable :: mw2(:), kappa(:), Vconst(:,:)
   integer :: nf_bath

contains

   subroutine init(nf_,ns_,mass_,omega_,Vconst_,kappa_)
      use pes, only : pesinit=>init, potptr=>pot, gradptr=>grad, nf, ns
      use pes, only : grad_a,grad_ab
      integer :: nf_, ns_
      real(dp) :: mass_(nf_), omega_(nf_)
      real(dp) :: Vconst_(ns_,ns_), kappa_(:)
!
!     Initialize module
!     kappa_ should have length nf/ns
!
      call pesinit(nf_,ns_,mass_)
      potptr => pot
      gradptr => grad
      nf_bath = nf/ns
      if (ns*nf_bath.ne.nf) stop 'frexc.f90: nf has do be divisible by ns'
      allocate(mw2(nf_bath),Vconst(ns,ns),kappa(nf_bath))
      mw2 = mass_(:nf_bath)*omega_(:nf_bath)**2
      kappa = kappa_
      Vconst = Vconst_
      ! Override generic adiabatic gradients with faster ones
      grad_a => gradad_diag
      grad_ab => gradad
   end subroutine


   subroutine pot(q, V)
      use pes, only : ns
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: V(ns,ns)
!
!     Diabatic potential matrix
!
      V = Vconst
      call addbath(q, V)
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
      real(dp), allocatable :: V1(:)
      allocate(V1(ns))
      V0 = sum(matmul(mw2,q**2))/2      
      V1 = matmul(kappa,q)
      do n=1,ns
         V(n,n) = V(n,n) + V0 + V1(n)
      end do
      deallocate(V1)
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
      use pes, only : ns
      real(dp), intent(in) :: q(nf_bath,ns)
      real(dp), intent(out) :: G(nf_bath,ns,ns,ns)
!
!     Auxiliary function for grad to help reshaping
!
      G = 0.d0
      do n=1,ns
         do m=1,ns
            G(:,n,m,m) = mw2*q(:,n)
         end do
         G(:,n,n,n) = G(:,n,n,n) + kappa
      end do
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
!
      integer :: a,b
      Gad = 0.d0
      do n=1,ns
         do a=1,ns
            Gad(:,n,a,a) = mw2*q(:,n)
            do b=1,ns
               Gad(:,n,a,b) = Gad(:,n,a,b) + U(n,a)*U(n,b)*kappa
            end do
         end do
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
      use pes, only : ns
      integer :: a
      real(dp), intent(in) :: q(nf_bath,ns), U(:,:)
      real(dp), intent(out) :: dvdq(nf_bath,ns)
!
!     Compute gradient of adiabatic potential a for special case of a Frenkel-exciton model
!     U -- transformation matrix from adia to dia basis
!     dvdq -- gradient of adiabatic state a, automatically reshaped to dvdq(nf)
!
      do n=1,ns
         dvdq(:,n) = mw2*q(:,n) + U(n,a)**2 * kappa
      end do
   end subroutine

end module
