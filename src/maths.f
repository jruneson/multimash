      module maths
         use types
         implicit double precision (a-h,o-z)

         complex(DPC), parameter :: iu = (0.0_dp,1.0_dp)
         complex(DPC), parameter :: c_zero = (0.0_dp,0.0_dp)
         complex(DPC), parameter :: c_one = (1.0_dp,0.0_dp)
         contains

         subroutine symevp (a,lda,n,d,ierr)
c
c     ------------------------------------------------------------------
c     This subroutine uses LAPACK DSYEVD to
c     diagonalise a real symmetric matrix.
c     ------------------------------------------------------------------
c
            dimension a(lda,n),d(n)
            allocatable :: work(:),iwork(:)
c
            lwork = 1+6*n+2*n*n
            liwork = 3+5*n
            allocate (work(lwork),iwork(liwork))
            call dsyevd ('V','U',n,a,lda,d,work,lwork,iwork,liwork,ierr)
            deallocate (work,iwork)
            return
         end


         subroutine gasdev (g,n)
c
c     ------------------------------------------------------------------
c     Generates an array of n normal deviates.
c     ------------------------------------------------------------------
c
            dimension g(n)
c
            twopi = 2*dacos(-1.d0)
            do j = 1,n,2
               ! x = 1.d0-random_number()! rand() is in [0,1) and can be exactly 0
               call random_number(x)
               x = 1.d0 - x
               call random_number(y)
               r = dsqrt(-2*dlog(x))
               phi = twopi*y
               g(j) = r*dcos(phi)
               if (j .eq. n) exit
               g(j+1) = r*dsin(phi)
            enddo
         end
      
      function trace (A)
         dimension A(:,:)
         trace = 0.d0
         n = size(A,1)
         do j=1,n
            trace = trace + A(j,j)
         end do
      end function

      end module
