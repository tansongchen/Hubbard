include 'parse.f90'

module Eigenvalue
    use parse
    implicit none
contains

pure function calculateKx(C, e) result(Kx)
    real(8), dimension(nEH, nEH) :: Kx 
    real(8), intent(in), dimension(:,:) :: C
    real(8), intent(in), dimension(:) :: e
    integer :: i, a, j, b, ia, jb

    do ia = 1, nEH
        do jb = 1, nEH
            call one_to_two(ia, i, a)
            call one_to_two(jb, j, b)
            Kx(ia,jb) = sum(C(:,i) * C(:,j) * C(:,a) * C(:,b))
        end do
    end do
end function calculateKx

function calculateKd(C, e, mode, W, omega) result(Kd)
    real(8), dimension(nEH, nEH) :: Kd
    real(8), intent(in), dimension(:,:) :: C
    real(8), intent(in), dimension(:) :: e
    character(*), intent(in) :: mode
    real(8), intent(in), dimension(:,:,:,:), optional :: W
    real(8), intent(in), optional :: omega

    real(8) :: omega_0, EA, EB, factor
    integer :: a, b, i, j, l, m, p, q, ia, jb

    Kd = 0
    select case (mode)
        case ('Hartree')
            Kd = 0
        case ('Hartree-Fock')
            do ia = 1, nEH
            do jb = 1, nEH
                call one_to_two(ia, i, a)
                call one_to_two(jb, j, b)
                Kd(ia,jb) = sum(C(:,i) * C(:,j) * C(:,a) * C(:,b))
            end do
            end do
        case default
            do ia = 1, nEH
            do jb = 1, nEH
                call one_to_two(ia, i, a)
                call one_to_two(jb, j, b)
                do l = 1, n
                do m = 1, n
                do p = 1, n
                do q = 1, n
                    Kd(ia,jb) = Kd(ia,jb) + C(l,a) * C(m,i) * C(p,b) * C(q,j) * W(l,m,p,q)
                end do
                end do
                    ! Kd(ia,jb) = Kd(ia,jb) + C(l,a) * C(m,i) * C(l,b) * C(m,j) * WS(l,m)
                end do
                end do
            end do
            end do
    end select
    if (mode == 'BSE_dynamic') then
        omega_0 = 1.747d0
        do ia = 1, nEH
        do jb = 1, nEH
            call one_to_two(ia, i, a)
            call one_to_two(jb, j, b)
            EA = omega - (e(b) - e(i))
            EB = omega - (e(a) - e(j))
            factor = (omega_0 / (omega_0 - EA) + omega_0 / (omega_0 - EB)) / 2
            ! print *, factor
            Kd(ia,jb) = Kd(ia,jb) * min(factor, 1d0)
        end do
        end do
    end if
end function calculateKd

pure function prepareCasida(D, Kx, Kd) result(M)
    real(8), dimension(2*nEH, 2*nEH) :: M
    real(8), intent(in), dimension(:,:) :: D, Kx, Kd
    real(8), dimension(nEH, nEH) :: A, B

    A = D + 2 * Kx - Kd
    B = 2 * Kx - Kd
    M(1:nEH, 1:nEH) = A
    M(1:nEH, nEH+1:2*nEH) = B
    M(nEH+1:2*nEH, 1:nEH) = -B
    M(nEH+1:2*nEH, nEH+1:2*nEH) = -A
end function prepareCasida

subroutine Casida(M, omega, V)
    ! LAPACK eigensolver requires complex vector and matrix to store the real and imaginary part,
    ! but the imaginary part should be 0, so we discard it
    real(8), intent(in), dimension(:,:) :: M
    real(8), intent(out), dimension(:) :: omega
    real(8), intent(out), dimension(:,:), optional :: V
    complex(8), dimension(2*nEH, 2*nEH) :: M_c
    complex(8), dimension(2*nEH, 2*nEH) :: Vl_c, Vr_c
    complex(8), dimension(2*nEH) :: omega_c
    integer :: index, index_

    M_c = M
    call geev(M_c, omega_c, Vl_c, Vr_c)
    omega = real(omega_c)
    if (present(V)) V = abs(Vr_c)
    ! We check that the imaginary part should be very small
    ! print *, maxval(aimag(omega_c))
    ! Sort the eigenvalues
    do index = 1, 2 * nEH
        do index_ = 1, 2 * nEH - 1
            if (omega(index_) > omega(index_ + 1)) then
                omega(index_:index_+1) = omega(index_+1:index_:-1)
                if (present(V)) V(:,index_:index_+1) = V(:,index_+1:index_:-1)
            end if
        end do
    end do
end subroutine Casida

subroutine CasidaTDA(M, omega, V)
    real(8), intent(in), dimension(:,:) :: M
    real(8), intent(out), dimension(:) :: omega
    real(8), intent(out), dimension(:,:), optional :: V

    real(8), dimension(nEH, nEH) :: M_
    integer :: index, index_
    
    M_ = M(1:nEH,1:nEH)
    call syev(M_, omega(1:nEH))
    if (present(V)) then
        V = 0
        V(1:nEH,1:nEH) = M_
    end if
    do index = 1, nEH
        do index_ = 1, nEH - 1
            if (omega(index_) > omega(index_ + 1)) then
                omega(index_:index_+1) = omega(index_+1:index_:-1)
                if (present(V)) V(:,index_:index_+1) = V(:,index_+1:index_:-1)
            end if
        end do
    end do
end subroutine
end module Eigenvalue

program main
    use Eigenvalue
    implicit none
    character(15), dimension(4), parameter :: modes = ['Hartree', 'Hartree-Fock', 'BSE_static', 'BSE_dynamic']
    character(15) :: mode
    integer :: i, iMode
    real(8), dimension(2*nEH, 2*nEH) :: M, V
    real(8), dimension(2*nEH) :: omegas
    real(8), dimension(nEH,nEH) :: D, Kx, Kd
    real(8), dimension(n,n,n,n) :: Q, W
    real(8), dimension(n,n) :: C, WS
    real(8), dimension(n) :: e
    real(8), parameter :: tolerance = 1d-4
    real(8) :: omegaL = 1.7d0, omegaR = 1.8d0, omega = 0d0
    real(8), dimension(4) :: results

    call readOrbitals(C, e)
    D = calculateD(e)
    Q = calculateQ(C, e)
    W = calculateW(Q)
    Kx = calculateKx(C, e)
    do iMode = 1, size(modes)
        mode = modes(iMode)
        if (mode == 'BSE_dynamic') then
            do while (omegaR - omegaL > tolerance)
                omega = (omegaL + omegaR) / 2
                Kd = calculateKd(C, e, mode, W, omega)
                M = prepareCasida(D, Kx, Kd)
                call Casida(M, omegas)
                if (omegas(nEH + 1) > omega) then
                    omegaL = omega
                else
                    omegaR = omega
                end if
            end do
            results(4) = omegas(nEH + 1)
        else
            Kd = calculateKd(C, e, mode, W)
            M = prepareCasida(D, Kx, Kd)
            call Casida(M, omegas)
            results(iMode) = omegas(nEH + 1)
        end if
    end do
    open(11, file='data/casida.dat')
    write(11, '(4f10.5)') results
    close(11)
end program