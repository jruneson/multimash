! Generic pes module

module pes
   use types
   integer :: nf ! number of nuclear degrees of freedom
   integer :: ns ! number of electronic states
   real(dp), allocatable :: mass(:)
   procedure (pot_interface), pointer :: pot
   procedure (grad_interface), pointer :: grad
   procedure (grad_a_interface), pointer :: grad_a
   procedure (grad_ab_interface), pointer :: grad_ab

   interface
      subroutine pot_interface(q,V)
         ! Diabatic potential matrix
         import ns,dp
         real(dp), intent(in) :: q(:)
         real(dp), intent(out) :: V(ns,ns)
      end subroutine pot_interface

      subroutine grad_interface(q,dVdq)
         ! Diabatic gradient matrix
         import nf,ns,dp
         real(dp), intent(in) :: q(:)
         real(dp), intent(out) :: dVdq(nf,ns,ns)
      end subroutine grad_interface

      subroutine grad_ab_interface(q,U,dVdq)
         ! Adiabatic gradient matrix
         import nf,ns,dp
         real(dp), intent(in) :: q(:), U(:,:)
         real(dp), intent(out) :: dVdq(nf,ns,ns)
      end subroutine grad_ab_interface

      subroutine grad_a_interface(q,U,a,dVdq)
         ! Adiabatic gradient of single state
         import nf,dp
         real(dp), intent(in) :: q(:), U(:,:)
         integer :: a
         real(dp), intent(out) :: dVdq(nf)
      end subroutine grad_a_interface
   end interface

   contains
   subroutine init(nf_,ns_,mass_)
      real(dp) :: mass_(:)
      nf = nf_
      ns = ns_
      if (.not. allocated(mass)) allocate(mass(nf))
      mass = mass_
      ! Default adiabatic routines: compute from diabatic 
      grad_a => gradad_diag
      grad_ab => gradad
   end subroutine

   subroutine potad(q,Vad,U)
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: Vad(ns), U(ns,ns)
!
!     Compute full set of adiabatic potentials and eigenstates
!
      call pot(q,U)
      call symevp(U,ns,ns,Vad,ierr)
   end subroutine

   subroutine gradad(q,U,Gad)
      real(dp), intent(in) :: q(:), U(:,:)
      real(dp), intent(out) :: Gad(nf,ns,ns)
!
!     Gradient in adiabatic representation
!     U -- rotation from adiabatic to diabatic representation (output of potad)
!
      real(dp), allocatable :: Gdia(:,:,:)

      allocate(Gdia(nf,ns,ns))
      
      ! get H' in diabatic basis
      call grad(q,Gdia)

      ! convert H' to adia
      do i=1,nf
         Gad(i,:,:) = matmul(transpose(U),matmul(Gdia(i,:,:),U))
      end do
      deallocate(Gdia)
   end subroutine

   subroutine gradad_diag(q, U, a, dvdq)
      real(dp), intent(in) :: q(:), U(:,:)
      integer :: a
      real(dp), intent(out) :: dvdq(nf)
!
!     Compute gradient of potential for a given adiabatic state.
!     (This is faster than "gradad", which computes all elements of the adiabatic gradient matrix)
!
      real(dp), allocatable :: Gdia(:,:,:)

      allocate(Gdia(nf,ns,ns))
      
      ! get H' in diabatic basis
      call grad(q,Gdia)

      ! convert H' to adia. Only compute needed bits of Gad
      do i=1,nf
         dvdq(i) = dot_product(U(:,a),matmul(Gdia(i,:,:),U(:,a)))
      end do

      deallocate(Gdia)
   end subroutine

   subroutine nac(q,Vad,U,d)
      real(dp), intent(in) :: q(:), Vad(:), U(:,:)
      real(dp), intent(out) :: d(nf,ns,ns)
!
!     Nonadiabatic coupling vector using Hellman-Feynman theorem
!
      real(dp), allocatable :: Gad(:,:,:)
      allocate(Gad(nf,ns,ns))
      call gradad(q,U,Gad)
      do k=1,ns
         d(:,k,k) = 0.d0
         do l=k+1,ns
            d(:,k,l) = Gad(:,k,l)/(Vad(l)-Vad(k))
            d(:,l,k) = -d(:,k,l)
         end do
      end do
      deallocate(Gad)
   end subroutine

   subroutine nacdir(q,cad,Vad,U,n,m,dj)
      real(dp), intent(in) :: q(:), Vad(:), U(:,:)
      complex(dpc), intent(in) :: cad(:)
      real(dp), intent(out) :: dj(nf)
!
!     Direction of momentum rescaling/reversal.
!     Note: does not include masses
!
      real(dp), allocatable :: d(:,:,:)
      allocate(d(nf,ns,ns))
      call nac(q,Vad,U,d)
      dj = 0.d0
      do k=1,ns
         dj = dj + d(:,k,n)*real(conjg(cad(k))*cad(n))
         dj = dj - d(:,k,m)*real(conjg(cad(k))*cad(m))
      end do
      deallocate(d)
   end subroutine

   subroutine symevp (a,lda,n,d,ierr)
      implicit double precision (a-h,o-z)
!
!     ------------------------------------------------------------------
!     This subroutine uses LAPACK DSYEVD to
!     diagonalise a real symmetric matrix.
!     ------------------------------------------------------------------
!
      dimension a(lda,n),d(n)
      allocatable :: work(:),iwork(:)
!
      lwork = 1+6*n+2*n*n
      liwork = 3+5*n
      allocate (work(lwork),iwork(liwork))
      call dsyevd ('V','U',n,a,lda,d,work,lwork,iwork,liwork,ierr)
      deallocate (work,iwork)
      return
   end

end module