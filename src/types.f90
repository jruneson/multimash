MODULE types

! Symbolic names for kind types of 4-, 2-, and 1-byte integers
INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)

! Symbolic names for kind types of single- and double-precision reals
INTEGER, PARAMETER :: SP = KIND(1.0)
INTEGER, PARAMETER :: DP = KIND(1.0D0)

! Symbolic names for kind types of single- and double-precision complex
INTEGER, PARAMETER :: SPC = KIND((1.0,1.0))
INTEGER, PARAMETER :: DPC = KIND((1.0D0,1.0D0))

! Symbolic name for kind type of default logical
INTEGER, PARAMETER :: LGT = KIND(.true.)

real(dp) :: pi = dacos(-1.d0)
complex(DPC), parameter :: iu = (0.0_dp,1.0_dp)

END MODULE types
