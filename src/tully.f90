! Tully canon potentials

module tully

   use types

   integer :: model
   ! Choices:
   ! 1 : Miller's modified version of Tully 1
   ! 2 : Tully 2
   ! 3 : Not currently implemented
   real(dp) :: A, B, C, D
   real(dp) :: E

contains
 
   subroutine init(model_,mass_)
      use pes, only : pesinit=>init, potptr=>pot, gradptr=>grad
      integer :: model_
      real(dp) :: mass_(1)
!
!     Initialize module
!
      call pesinit(nf_,ns_,mass_)
      potptr => pot
      gradptr => grad
      call pesinit(1,2,mass_)
      model = model_
      if (model.eq.1) then
         A = 0.01d0
         B = 1.6d0
         C = 0.005d0
         D = 1.d0
      else if (model.eq.2) then
         A = 0.1d0
         B = 0.28d0
         C = 0.015d0
         D = 0.06d0
         E = 0.05d0
      else
         stop 'Model not implemented'
      end if
   end subroutine


   subroutine pot(q, V)
      use pes, only : ns
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: V(ns,ns)
!
!     Diabatic potential matrix
!
      if (model.eq.1) then
         V(1,1) = A*tanh(B*q(1))
         V(2,2) = -V(1,1)
         V(1,2) = C*exp(-D*q(1)**2)
         V(2,1) = V(1,2)
      else if (model.eq.2) then
         V(1,1) = 0.d0
         V(2,2) = E - A*exp(-B*q(1)**2)
         V(1,2) = C*exp(-D*q(1)**2)
         V(2,1) = V(1,2)
      end if
   end subroutine
 
 
   subroutine grad(q, G)
      use pes, only : nf,ns
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: G(nf,ns,ns)
!
!     Analytic gradient of diabatic potential 
!
      if (model.eq.1) then
         G(1,1,1) = A*B/cosh(B*q(1))**2
         G(1,2,2) = -G(1,1,1)
         G(1,1,2) = -2.d0*C*D*q(1)*exp(-D*q(1)**2)
         G(1,2,1) = G(1,1,2)
      else if (model.eq.2) then
         G(1,1,1) = 0.d0
         G(:,2,2) = 2.d0 * A * B * q * exp(-B*q**2)
         G(:,1,2) = -2.d0 * C * D * q * exp(-D*q**2)
         G(:,2,1) = G(:,1,2) 
      end if
   end subroutine
 
 
 
end module
 