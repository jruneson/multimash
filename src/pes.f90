! Generic pes module

module pes
   use types
   integer :: nf ! number of nuclear degrees of freedom
   integer :: ns ! number of electronic states
   real(dp), allocatable :: mass(:), omega(:)
   procedure (pot_interface), pointer :: pot
   procedure (grad_interface), pointer :: grad
   procedure (grad_a_interface), pointer :: grad_a
   procedure (grad_ab_interface), pointer :: grad_ab
   procedure (grad_av_interface), pointer :: grad_av
   procedure (get_vlin_interface), pointer :: get_vlin
   procedure (get_vconst_interface), pointer :: get_vconst
   logical :: cayley = .false.

   interface
      subroutine pot_interface(q,V)
         ! Diabatic potential matrix
         import ns,nf,dp
         real(dp), intent(in) :: q(nf)
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

      subroutine grad_av_interface(q,U,Gad,a)
         ! Adiabatic gradient vector
         import nf,ns,dp
         real(dp), intent(in) :: q(:), U(:,:)
         integer :: a
         real(dp), intent(out) :: Gad(nf,ns)
      end subroutine grad_av_interface

      subroutine grad_a_interface(q,U,a,dVdq)
         ! Adiabatic gradient of single state
         import nf,dp
         real(dp), intent(in) :: q(:), U(:,:)
         integer :: a
         real(dp), intent(out) :: dVdq(nf)
      end subroutine grad_a_interface

      subroutine get_vconst_interface(Vconst)
         ! Constant part of potential
         import ns,dp
         real(dp), intent(out) :: Vconst(ns,ns)
      end subroutine get_vconst_interface

      subroutine get_vlin_interface(Vlin)
         ! Linear part of potential
         import nf,ns,dp
         real(dp), intent(out) :: Vlin(nf,ns,ns)
      end subroutine get_vlin_interface
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
      grad_av => grad_vector
   end subroutine

   subroutine potad(q,Vad,U)
      use maths, only : symevp
      real(dp), intent(in) :: q(:)
      real(dp), intent(out) :: Vad(ns), U(ns,ns)
!
!     Compute full set of adiabatic potentials and eigenstates
!
      call pot(q,U)
      call symevp(U,ns,ns,Vad,ierr)
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

   subroutine grad_vector(q,U,Gad,b)
      real(dp), intent(in) :: q(:), U(:,:)
      integer :: b
      real(dp), intent(out) :: Gad(nf,ns)
!
!     Returns Gad(i,a,b) for a given adiabatic state b
!
      real(dp), allocatable :: Gdia(:,:,:)

      allocate(Gdia(nf,ns,ns))
      
      ! get H' in diabatic basis
      call grad(q,Gdia)

      ! convert H' to adia.
      do i=1,nf
         Gad(i,:) = matmul(transpose(U),matmul(Gdia(i,:,:),U(:,b)))
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
      call grad_ab(q,U,Gad)
      do k=1,ns
         d(:,k,k) = 0.d0
         do l=k+1,ns
            d(:,k,l) = Gad(:,k,l)/(Vad(l)-Vad(k))
            d(:,l,k) = -d(:,k,l)
         end do
      end do
      deallocate(Gad)
   end subroutine

   subroutine nac_a(q,b,Vad,U,d)
      real(dp), intent(in) :: q(:), Vad(:), U(:,:)
      integer :: a,b
      real(dp), intent(out) :: d(nf,ns)
!
!     Returns components d(i,a,b) for a given adiabatic state b
!
      real(dp), allocatable :: Gad(:,:)
      allocate(Gad(nf,ns))
      call grad_av(q,U,Gad,b)
      do a=1,ns
         if (a.eq.b) then
            d(:,a) = 0.d0
         else
            d(:,a) = Gad(:,a)/(Vad(b)-Vad(a))
         end if
      end do
      deallocate(Gad)
   end subroutine


   subroutine nacdir(q,cad,Vad,U,a,b,dj)
      real(dp), intent(in) :: q(:), Vad(:), U(:,:)
      complex(dpc), intent(in) :: cad(:)
      real(dp), intent(out) :: dj(nf)
      integer :: a,b
!
!     Direction of momentum rescaling/reversal.
!     Note: does not include masses
!
      real(dp), allocatable :: d(:,:)
      allocate(d(nf,ns))
      call nac_a(q,a,Vad,U,d)
      dj = 0.d0
      do k=1,ns
         dj = dj + d(:,k)*real(conjg(cad(k))*cad(a))
      end do
      call nac_a(q,b,Vad,U,d)
      do k=1,ns
         dj = dj - d(:,k)*real(conjg(cad(k))*cad(b))
      end do
      deallocate(d)
   end subroutine

end module