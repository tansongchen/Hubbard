include '0_Hubbard.f90'

module parse
    use Hubbard
    implicit none
contains

subroutine readOrbitals(C, e)
    real(8), intent(out), dimension(:,:) :: C
    real(8), intent(out), dimension(:) :: e
    integer :: i

    open(10, file='result/HF.dat')
    read(10, *) e
    do i = 1, n
        read(10, *) C(:,i)
    end do
    close(10)
end subroutine readOrbitals

subroutine readFockOperator(F)
    real(8), intent(out), dimension(:,:) :: F
    integer :: i

    open(10, file='result/HF.dat')
    do i = 1, n + 1
        read(10, *) F(:,1)
    end do
    do i = 1, n
        read(10, *) F(:,i)
    end do
    close(10)
end subroutine readFockOperator

pure function calculateZ(C) result(Z)
    real(8), dimension(n,n) :: Z
    real(8), intent(in), dimension(:,:) :: C
    integer :: i
    real(8), dimension(n,n) :: D

    D = 0
    forall (i=1:n)
        D(i,i) = i
    end forall
    Z = matmul(transpose(C), matmul(D, C))
end function

pure subroutine one_to_two(ia, i, a)
    integer, intent(in) :: ia
    integer, intent(out) :: i, a
    i = (ia - 1) / nVir + 1
    a = ia - (i - 1) * nVir + nOcc
end subroutine

pure subroutine two_to_one(ia, i, a)
    integer, intent(in) :: i, a
    integer, intent(out) :: ia
    ia = (i - 1) * nVir + (a - nOcc)
end subroutine

pure subroutine one_to_two_w(lp, l, p)
    integer, intent(in) :: lp
    integer, intent(out) :: l, p
    l = (lp - 1) / n + 1
    p = lp - (l - 1) * n
end subroutine

pure subroutine two_to_one_w(lp, l, p)
    integer, intent(in) :: l, p
    integer, intent(out) :: lp
    lp = (l - 1) * n + p
end subroutine

pure function calculateQ(C, e) result(Q)
    real(8), dimension(n,n,n,n) :: Q
    real(8), intent(in), dimension(:,:) :: C
    real(8), intent(in), dimension(:) :: e
    integer :: l, s, p, t, i, j, a, b

    Q = 0
    do l = 1, n
    do s = 1, n
    do p = 1, n
    do t = 1, n
        if (l == p) then
            do i = 1, nOcc
            do a = 1, nVir
                Q(l,s,p,t) = Q(l,s,p,t) - 4 / (e(nOcc + a) - e(i)) * C(l,i) * C(p,a) * C(s,a) * C(t,i)
            end do
            end do
        end if
    end do
    end do
    end do
    end do
end function calculateQ

pure function calculateWS(C, e) result(WS)
    real(8), dimension(n,n) :: WS
    real(8), intent(in), dimension(:,:) :: C
    real(8), intent(in), dimension(:) :: e
    integer :: l, s, i, a
    real(8), dimension(n,n) :: Q
    integer, dimension(n) :: ipiv

    Q = 0
    do l = 1, n
    do s = 1, n
        do i = 1, nOcc
        do a = nOcc + 1, nOcc + nVir
            Q(l,s) = Q(l,s) - 4 / (e(a) - e(i)) * C(l,i) * C(l,a) * C(s,a) * C(s,i)
        end do
        end do
    end do
    end do
    WS = 0
    forall (l=1:n)
        WS(l,l) = 1
    end forall
    WS = WS - Q
    call getrf(WS, ipiv)
    call getri(WS, ipiv)
end function calculateWS

function calculateW(Z) result(W)
    real(8), dimension(n,n,n,n) :: W
    real(8), intent(in), dimension(:,:,:,:) :: Z

    real(8), dimension(n2,n2) :: W_matrix
    real(8), dimension(n2,n2) :: Z_matrix
    integer, dimension(n2) :: ipiv
    integer :: iW, l, m, p, q, lp, mq

    Z_matrix = 0
    do l = 1, n
    do m = 1, n
    do p = 1, n
    do q = 1, n
        call two_to_one_w(lp, l, p)
        call two_to_one_w(mq, m, q)
        Z_matrix(lp,mq) = Z(l,m,p,q)
    end do
    end do
    end do
    end do
    W_matrix = 0
    do iW = 1, n2
        W_matrix(iW,iW) = 1
    end do
    W_matrix = W_matrix - Z_matrix
    call getrf(W_matrix, ipiv)
    call getri(W_matrix, ipiv)
    ! print '(64f13.10)', W_matrix
    ! left-multiply it by U
    do iW = 1, n2
        call one_to_two_w(iW, l, p)
        if (l /= p) W_matrix(iW,:) = 0
    end do
    W = 0
    do l = 1, n
    do m = 1, n
    do p = 1, n
    do q = 1, n
        call two_to_one_w(lp, l, p)
        call two_to_one_w(mq, m, q)
        W(l,m,p,q) = W_matrix(lp,mq)
    end do
    end do
    end do
    end do
end function calculateW

pure function calculateD(e) result(D)
    real(8), dimension(nEH, nEH) :: D
    real(8), intent(in), dimension(:) :: e
    integer :: i, a, ia

    do ia = 1, nEH
        call one_to_two(ia, i, a)
        D(ia,ia) = e(a) - e(i)
    end do
end function calculateD
end module