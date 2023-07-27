! Linear vibronic coupling potential

module bilinvib

   use types

   real(dp), allocatable :: omega(:), Vconst(:,:), Vlin(:,:,:)
   real(dp), allocatable :: Vbilin(:,:,:,:)

contains

   subroutine init(nf_,ns_,mass_,omega_,Vconst_,Vlin_,Vbilin_)
      use pes, only : pesinit=>init, potptr=>pot, gradptr=>grad, nf, ns
      integer :: nf_, ns_
      real(dp) :: mass_(nf_), omega_(nf_)
      real(dp) :: Vconst_(ns_,ns_), Vlin_(nf_,ns_,ns_), &
                  Vbilin_(nf_,nf_,ns_,ns_)
!
!     Initialize module
!
      call pesinit(nf_,ns_,mass_)
      potptr => pot
      gradptr => grad
      allocate(omega(nf),Vconst(ns,ns),Vlin(nf,ns,ns),&
               Vbilin(nf,nf,ns,ns))
      omega = omega_
      Vconst = Vconst_
      Vlin = Vlin_
      Vbilin = Vbilin_
   end subroutine


   subroutine pot(q, V)
      use pes, only : mass,nf,ns
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: V(ns,ns)
!
!     Diabatic potential matrix
!
      real(dp) :: V0
      V0 = 0.5d0 * sum(mass*(omega*q)**2)
      do n=1,ns
         do m=1,ns
            V(n,m) = Vconst(n,m) + sum(Vlin(:,n,m)*q)
            do j=1,nf
               V(n,m) = V(n,m) + sum(Vbilin(:,j,n,m)*q(j)*q(:))
            end do
         end do
         V(n,n) = V(n,n) + V0
      end do
   end subroutine


   subroutine grad(q, G)
      use pes, only : mass,nf,ns
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: G(nf,ns,ns)
!
!     Analytic gradient of diabatic potential 
!
      real(dp), allocatable :: G0(:)
      allocate(G0(nf))
      G0 = mass*omega**2 * q
      G = Vlin
      do j=1,nf
         G = G + 2.d0*Vbilin(:,j,:,:)*q(j)
      end do
      do n=1,ns
         G(:,n,n) = G(:,n,n) + G0
      end do
      deallocate(G0)
   end subroutine



end module
