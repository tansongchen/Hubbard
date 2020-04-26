include '0_Hubbard.f90'

module HF
    use Hubbard
    implicit none
contains

pure function FockMatrix(T, P) result(F)
    real(8), dimension(:,:), intent(in) :: T(:,:), P(:,:)
    real(8), dimension(N,N) :: F
    integer :: l

    F = T
    forall (l=1:N) 
        F(l,l) = F(l,l) + unit * P(l,l) / 2
    end forall
end function FockMatrix

subroutine EigenSolver(F, C, e)
    real(8), dimension(:,:), intent(in) :: F
    real(8), dimension(:,:), intent(out) :: C
    real(8), dimension(:), intent(out) :: e
    integer :: i, j

    C = F
    call syev(C, e, 'v')
    ! In HF calculation, matrix C must be sorted according to eigenvalue,
    ! so that P = CθC^T is correct.
    do i = 1, N
        do j = 1, N - 1
            if (e(j) > e(j+1)) then
                e(j:j+1) = e(j+1:j:-1)
                C(:,j:j+1) = C(:,j+1:j:-1)
            end if
        end do
    end do
end subroutine EigenSolver
end module

program main
    use HF
    implicit none
    real(8), dimension(N,N) :: T = 0, P = 0, P_ = 0, F = 0, C = 0, C_ = 0
    real(8), dimension(N) :: energy = 0
    integer, parameter :: maxIter = 1e4
    real(8), parameter :: tolerance = 1d-10
    integer :: iter, l
    real t1, t2

    T = HoppingMatrix()
    print *, 'Self-consistent field equation solver start...'
    call cpu_time(t1)
    P = 0
    P_ = 1
    iter = 0
    do while (sum((P - P_)**2) > tolerance .and. iter <= maxIter)
        iter = iter + 1
        F = FockMatrix(T, P)
        call EigenSolver(F, C, energy)
        C_(:,1:Nocc) = C(:,1:Nocc)
        P_ = P
        P = 2 * matmul(C_, transpose(C_))
    end do
    if (iter > maxIter) then
        print *, 'Failed to reach self-consistency!'
        stop
    else
        open(10, file='output/bse/HF.dat')
        call cpu_time(t2)
        print *, 'Successful to reach self-consistency, using', iter, 'cycles and', t2-t1, 'seconds.'
        write (10, '(16e25.17)') energy
        do l = 1, N
            write (10, '(16e25.17)') C(:,l)
        end do
        F = FockMatrix(T, P)
        do l = 1, N
            write (10, '(16e25.17)') F(:,l)
        end do
        close(10)
    end if
end program