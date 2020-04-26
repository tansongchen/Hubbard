include 'lapack.f90'

module Hubbard
use lapack95
use f95_precision
implicit none

real(8), parameter :: alpha = 1.5d0, beta = 1d0, unit = 1d0
integer, parameter :: n = 8, n2 = n * n
integer, parameter :: nOcc = 4, nVir = n - nOcc, nEH = nOcc * nVir

contains

pure function HoppingMatrix() result(T)
    real(8), dimension(N,N) :: T
    integer :: l

    do l = 1, N - 1
        if (mod(l, 2) == 1) then
            ! l = 1, 3, 5, 7
            T(l,l+1) = alpha
            T(l+1,l) = alpha
        else
            ! l = 2, 4, 6
            T(l,l+1) = beta
            T(l+1,l) = beta
        end if
    end do
end function HoppingMatrix
end module
    