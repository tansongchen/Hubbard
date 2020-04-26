include 'parse.f90'

module CI
    use parse
    implicit none
    integer, parameter :: nFrontier = 4
    real(8), dimension(nFrontier,nFrontier) :: A = 0
    real(8), dimension(n,n) :: H1, S1
    real(8), dimension(n,n,n,n) :: H2
    real(8), dimension(:,:), allocatable, target :: H, S
    real(8), dimension(:,:), pointer :: F => null()
    real(8), dimension(:), allocatable :: EN
    integer, dimension(:,:), allocatable :: O
    integer, parameter :: nElec = nOcc * 2, nHole = nVir * 2, nTot = nElec + nHole
    integer :: nDet = 0

contains

subroutine getHamiltonian()
    real(8), dimension(n,n) :: C
    real(8), dimension(n) :: e
    real(8), dimension(n,n2) :: CC
    real(8), dimension(n,n) :: T
    integer :: i, j, k, l, ij, kl
    
    T = HoppingMatrix()
    call readOrbitals(C, e)
    H1 = matmul(matmul(transpose(C), T), C)
    S1 = calculateZ(C)
    do i = 1, n
    do j = 1, n
    do k = 1, n
    do l = 1, n
        H2(i,j,k,l) = unit * sum(C(:,i) * C(:,j) * C(:,k) * C(:,l))
    end do
    end do
    end do
    end do
end subroutine

subroutine initialize()
    integer(8) :: i, j
    integer(1), dimension(nTot) :: bool
    integer(1), dimension(nTot, 10000) :: buffer

    nDet = 0
    do i = 0, 2**nTot - 1
        do j = 0, nTot - 1
            bool(j+1) = ibits(i,j,1)
        end do
        if (sum(bool) == nElec .and. sum(bool(1::2)) == (sum(bool(2::2)) + mod(nElec, 2))) then
            nDet = nDet + 1
            buffer(:,nDet) = bool
        end if
    end do
    allocate(H(nDet, nDet), S(nDet, nDet), EN(nDet), O(nTot,nDet))
    H = 0
    S = 0
    F => H(:,1:nFrontier)
    EN = 0
    O = buffer(:, 1:nDet)
end subroutine

subroutine finalize()
    F => null()
    deallocate(H, S, EN, O)
end subroutine

pure real(8) function oneBody(i1, j1) result(o1)
    ! calculate <i1|h|j1>
    integer, intent(in) :: i1, j1
    integer :: spatialI1, spatialJ1

    spatialI1 = (i1 + 1) / 2
    spatialJ1 = (j1 + 1) / 2
    if (mod(i1 - j1, 2) == 0) then
        o1 = H1(spatialI1, spatialJ1)
    else
        o1 = 0
    end if
end function

pure real(8) function opticalExcitation(i1, j1) result(o1)
    ! calculate <i1|z|j1>
    integer, intent(in) :: i1, j1
    integer :: spatialI1, spatialJ1

    spatialI1 = (i1 + 1) / 2
    spatialJ1 = (j1 + 1) / 2
    if (mod(i1 - j1, 2) == 0) then
        o1 = S1(spatialI1, spatialJ1)
    else
        o1 = 0
    end if
end function

pure real(8) function twoBody(i1, i2, j1, j2) result(o2)
    ! calculate <i1 i2||j1 j2>, physicist's notation
    integer, intent(in) :: i1, i2, j1, j2
    integer :: spatialI1, spatialI2, spatialJ1, spatialJ2
    spatialI1 = (i1 + 1) / 2
    spatialI2 = (i2 + 1) / 2
    spatialJ1 = (j1 + 1) / 2
    spatialJ2 = (j2 + 1) / 2
    o2 = 0
    if (mod(i1 - j1, 2) == 0 .and. mod(i2 - j2, 2) == 0) then
        o2 = o2 + H2(spatialI1,spatialI2,spatialJ1,spatialJ2)
    end if
    if (mod(i1 - j2, 2) == 0 .and. mod(i2 - j1, 2) == 0) then
        o2 = o2 - H2(spatialI1,spatialI2,spatialJ2,spatialJ1)
    end if
end function

subroutine parse1(boolI, boolJ, i1, j1, phase)
    integer, dimension(:), intent(in) :: boolI, boolJ
    integer, intent(out) :: i1, j1, phase
    integer :: switch = 0, i
    phase = -1
    do i = 1, nTot
        if (boolI(i) /= boolJ(i)) then
            switch = 1 - switch
            if (boolI(i) == 1) then
                i1 = i
            else
                j1 = i
            end if
        else
            if (switch == 1 .and. boolI(i) == 1) phase = -phase
        end if
    end do
end subroutine

subroutine parse2(boolI, boolJ, i1, i2, j1, j2, phase)
    integer, dimension(:), intent(in) :: boolI, boolJ
    integer, intent(out) :: i1, j1, i2, j2, phase
    integer :: switch = 0, i
    i1 = 0
    i2 = 0
    j1 = 0
    j2 = 0
    phase = 1
    do i = 1, nTot
        if (boolI(i) /= boolJ(i)) then
            switch = 1 - switch
            if (boolI(i) == 1) then
                if (i1 == 0) then
                    i1 = i
                else
                    i2 = i
                end if
            else
                if (j1 == 0) then
                    j1 = i
                else
                    j2 = i
                end if
            end if
        else
            if (switch == 1 .and. boolI(i) == 1) phase = -phase
        end if
    end do
end subroutine

subroutine calculateH()
    integer :: i, j, k, l, i1, i2, j1, j2, diff, phase
    integer, dimension(nTot) :: boolI, boolJ

    do i = 1, nDet
        do j = i, nDet
            boolI = O(:,i)
            boolJ = O(:,j)
            diff = sum(abs(boolJ - boolI))
            if (diff == 0) then
                ! H(i,j) = sum([(oneBody(k, k), k = 1, nTot)], mask=boolI == 1) + sum([((twoBody(k, l, k, l), k = 1, nTot), l = 1, nTot)], mask=[((boolI(k) * boolJ(l), k = 1, nTot), l = 1, nTot)] == 1) / 2
                H(i,j) = 0
                do k = 1, nTot
                    if (boolI(k) == 1) H(i,j) = H(i,j) + oneBody(k,k)
                end do
                do k = 1, nTot
                    do l = 1, nTot
                        if (boolI(k) == 1 .and. boolI(l) == 1) H(i,j) = H(i,j) + twoBody(k,l,k,l) / 2
                    end do
                end do
            else if (diff == 2) then
                call parse1(boolI, boolJ, i1, j1, phase)
                H(i,j) = phase * (oneBody(i1, j1) + sum([(twoBody(i1, k, j1, k), k = 1, nTot)], mask=boolI * boolJ == 1))
            else if (diff == 4) then
                call parse2(boolI, boolJ, i1, i2, j1, j2, phase)
                H(i,j) = phase * twoBody(i1, i2, j1, j2)
            else if (diff > 4) then
                H(i,j) = 0
            else
                print '(16I2)', boolI, boolJ
                error stop 'Illegal Difference'
            end if
        end do
    end do
    do i = 1, nDet
        do j = 1, i-1
            H(i,j) = H(j,i)
        end do
    end do
end subroutine

subroutine calculateS()
    integer :: i, j, k, l, i1, i2, j1, j2, diff, phase
    integer, dimension(nTot) :: boolI, boolJ

    do i = 1, nDet
        do j = i, nDet
            boolI = O(:,i)
            boolJ = O(:,j)
            diff = sum(abs(boolJ - boolI))
            if (diff == 0) then
                S(i,j) = sum([(opticalExcitation(k, k), k = 1, nTot)], mask=boolI == 1)
            else if (diff == 2) then
                call parse1(boolI, boolJ, i1, j1, phase)
                S(i,j) = phase * opticalExcitation(i1, j1)
            else if (diff > 2) then
                S(i,j) = 0
            else
                print '(16I2)', boolI, boolJ
                error stop 'Illegal Difference'
            end if
        end do
    end do
    do i = 1, nDet
        do j = 1, i-1
            S(i,j) = S(j,i)
        end do
    end do
end subroutine

end module

program main
    use CI
    implicit none

    call getHamiltonian()
    call initialize()
    print *, 'Initialize Finished'
    call calculateH()
    call calculateS()
    print *, 'H/S Calculation Finished'
    call syev(H, EN, 'v')
    print *, 'SYEV Finished'
    A = matmul(transpose(F), matmul(S, F))
    open(10, file='output/bse/ci.dat')
    write(10, '(4f15.7)') EN(2:nFrontier+1) - EN(1)
    write(10, '(4f15.7)') A(1,2:nFrontier+1)
    close(10)
    call finalize()
end program