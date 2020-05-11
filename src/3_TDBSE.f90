include 'parse.f90'

module tdbse
    use parse
    implicit none
    complex(8), parameter :: Im = (0, 1)
contains

pure function propagator(O, dt) result(U)
    complex(8), dimension(n,n) :: U
    real(8), intent(in), dimension(:,:) :: O
    real(8), intent(in) :: dt

    real(8), dimension(n,n) :: Q
    real(8), dimension(n) :: w
    integer :: i

    Q = O
    w = 0
    U = 0
    call syev(Q, w, 'v')
    forall (i=1:n)
        U(i,i) = exp(-Im * w(i) * dt)
    end forall
    U = matmul(matmul(Q, U), transpose(Q))
end function propagator

pure function propagatorHermitian(O, dt) result(U)
    complex(8), dimension(n,n) :: U
    complex(8), intent(in), dimension(:,:) :: O
    real(8), intent(in) :: dt

    complex(8), dimension(n,n) :: Q
    real(8), dimension(n) :: w
    integer :: i

    Q = O
    w = 0
    U = 0
    call heev(Q, w, 'v')
    forall (i=1:n)
        U(i,i) = exp(-Im * w(i) * dt)
    end forall
    U = matmul(matmul(Q, U), conjg(transpose(Q)))
    ! U = 0
    ! forall (i=1:n)
    !     U(i,i) = 1
    ! end forall
    ! U = U + Im * O * dt
end function propagatorHermitian

pure function propagatorGeneral(O, dt) result(U)
    complex(8), dimension(n,n) :: U
    complex(8), intent(in), dimension(:,:) :: O
    real(8), intent(in) :: dt
    integer :: i

    U = 0d0
    forall (i=1:n)
        U(i,i) = 1d0
    end forall
    U = U - Im * O * dt
end function propagatorGeneral

pure function hartree(C, dt) result(J)
    complex(8), dimension(n,n) :: J
    complex(8), intent(in), dimension(:,:) :: C
    real(8), intent(in) :: dt

    complex(8), dimension(n,n) :: P
    real(8), dimension(n,n) :: G
    integer :: i

    P = 2 * matmul(C, conjg(transpose(C)))
    G = 0
    forall (i=1:n)
        G(i,i) = P(i,i)
    end forall
    J = propagator(G, dt)
end function

function fockS(C, dt, mode, WS, CPrev, HQP, plasmon) result(K)
    complex(8), dimension(n,n) :: K
    complex(8), intent(in), dimension(:,:) :: C, CPrev
    real(8), intent(in) :: dt, plasmon
    character(*), intent(in) :: mode
    ! real(8), intent(in), dimension(:,:,:,:) :: W
    real(8), intent(in), dimension(:,:) :: WS
    real(8), intent(in), dimension(n,n) :: HQP

    complex(8), dimension(n,n) :: P, PPrev, PDiff, PComm, Gc
    real(8), dimension(n,n) :: G
    integer :: i, l, q, m, pp

    K = 0
    if (mode == 'Hartree') then
        forall (i=1:n)
            K(i,i) = 1
        end forall
    else if (mode == 'Hartree-Fock') then
        P = 2 * matmul(C, conjg(transpose(C)))
        G = 0
        forall (i=1:n)
            G(i,i) = -P(i,i) / 2
        end forall
        K = propagator(G, dt)
    else if (mode == 'BSE_static') then
        P = 2 * matmul(C, conjg(transpose(C)))
        Gc = -WS * P / 2
        K = propagatorHermitian(Gc, dt)
    else if (mode == 'BSE_dynamic') then
        P = 2 * matmul(C, conjg(transpose(C)))
        Gc = -WS * P / 2
        K = propagatorHermitian(Gc, dt)
    ! else if (mode == 'BSE_dynamic') then
    !     P = 2 * matmul(C, conjg(transpose(C)))
    !     PPrev = 2 * matmul(CPrev, conjg(transpose(CPrev)))
    !     Gc = -WS * P / 2
    !     PDiff = (P - PPrev) / dt * Im - (matmul(P, HQP) - matmul(HQP, P))
    !     Gc = Gc + (matmul(Gc, HQP) - matmul(HQP, Gc)) * Im / plasmon
    !     Gc = Gc - WS * (PDiff + PComm) / plasmon / 2
    !     K = propagatorHermitian(Gc, dt / 2)
    else
        K = 0
    end if
end function

pure function dynamic(C, e, W) result(WS)
    real(8), dimension(n,n) :: WS
    real(8), intent(in), dimension(:,:) :: C
    real(8), intent(in), dimension(:) :: e
    real(8), intent(in), dimension(:,:,:,:) :: W

    real(8), dimension(n,n,n,n) :: WF, WC
    real(8) :: omega_0, omega, EA, EB, factor
    integer :: a, b, i, j, l, m, p, q, ia, jb, l1, l2

    WF = 0
    ! do ia = 1, nEH
    ! do jb = 1, nEH
    !     call one_to_two(ia, i, a)
    !     call one_to_two(jb, j, b)
    do i = 1, n
    do a = 1, n
    do j = 1, n
    do b = 1, n
        do l = 1, n
        do m = 1, n
        do p = 1, n
        do q = 1, n
            WF(i,a,j,b) = WF(i,a,j,b) + C(l,a) * C(m,i) * C(p,b) * C(q,j) * W(l,m,p,q)
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do

    omega_0 = 1.747d0
    omega = 1.747d0
    do ia = 1, nEH
    do jb = 1, nEH
        call one_to_two(ia, i, a)
        call one_to_two(jb, j, b)
        EA = omega - (e(b) - e(i))
        EB = omega - (e(a) - e(j))
        factor = (omega_0 / (omega_0 - EA) + omega_0 / (omega_0 - EB)) / 2
        ! print *, factor
        WF(i,a,j,b) = WF(i,a,j,b) * min(factor, 1d0)
    end do
    end do

    WC = 0
    do l = 1, n
    do m = 1, n
    do p = 1, n
    do q = 1, n
        ! do ia = 1, nEH
        ! do jb = 1, nEH
        ! call one_to_two(ia, i, a)
        ! call one_to_two(jb, j, b)
        do i = 1, n
        do a = 1, n
        do j = 1, n
        do b = 1, n
            WC(l,m,p,q) = WC(l,m,p,q) + C(l,a) * C(m,i) * C(p,b) * C(q,j) * WF(i,a,j,b)
        end do
        end do
        end do
        end do
    end do
    end do
    end do
    end do

    WS = 0
    do l1 = 1, n
    do l2 = 1, n
        ! WS(l1,l2) = (sum(W(l1,:,l1,l2)) + sum(W(l1,l2,l1,:))) / n
        WS(l1,l2) = (sum(WC(l1,l2,:,:)) + sum(WC(l1,:,:,l2)) + sum(WC(:,l2,l1,:)) + sum(WC(:,:,l1,l2))) / 2 / n
    end do
    end do
end function

subroutine perturb(C1, gamma)
    complex(8), intent(inout), dimension(n,nOcc) :: C1
    real(8), intent(in) :: gamma
    integer :: l

    do l = 1, n
        C1(l,:) = C1(l,:) * exp(Im * gamma * l)
    end do
end subroutine perturb
end module tdbse

program main
    use tdbse
    implicit none
    integer i, l, l1, l2, l3, l4, iMode
    character(15), dimension(4), parameter :: modes = ['Hartree', 'Hartree-Fock', 'BSE_static', 'BSE_dynamic']
    character(15) :: mode
    ! matrices
    real(8), parameter :: gamma = 0.0001
    real(8), dimension(n,n) :: T = 0, C_ = 0, HQP = 0, WS = 0
    real(8), dimension(n,nOcc) :: C = 0
    complex(8), dimension(n,nOcc) :: C0 = 0, C1 = 0, C0Prev = 0, C1Prev = 0
    complex(8), dimension(n,n) :: H = 0, J0 = 0, J1 = 0, K0 = 0, K1 = 0, P0 = 0, P1 = 0
    real(8), dimension(n) :: e
    real(8), dimension(n,n,n,n) :: Q, W = 0
    ! propagation
    real(8), parameter :: omegaMin = 1.65d0, omegaMax = 1.8d0, omegaInterval = 1d-3
    integer, parameter :: nOmega = nint((omegaMax - omegaMin) / omegaInterval)
    real(8), parameter, dimension(0:nOmega) :: omegas = [(omegaMin + omegaInterval * i, i = 0, nOmega)]
    real(8), parameter :: tMax = 2 / omegaInterval, dt = 1d-1
    integer, parameter :: nsteps = nint(tMax / dt)
    real(8), dimension(nsteps) :: ts = [(dt * i, i = 1, nsteps)]
    ! fourier transform
    real(8), dimension(nsteps) :: dacf = 0d0
    complex(8), dimension(nsteps) :: fourier = 0d0
    real(8), dimension(0:nOmega,4) :: spectrum = 0d0
    integer :: iTStep, iOmega
    real(8) :: omega, plasmon

    T = HoppingMatrix()
    call readOrbitals(C_, e)
    call readFockOperator(HQP)
    ! plasmon = e(nOcc + 1) - e(nOcc)
    plasmon = 10000
    Q = calculateQ(C_, e)
    W = calculateW(Q)
    ! WS = calculateWS(C_, e)
    C = C_(:,1:nOcc)
    H = propagator(T, dt / 2)
    do iMode = 1, 4
        mode = modes(iMode)
        if (mode == 'BSE_static') then
            do l1 = 1, n
            do l2 = 1, n
                WS(l1,l2) = (sum(W(l1,:,l1,l2)) + sum(W(l1,l2,l1,:))) / n
            end do
            end do
        else if (mode == 'BSE_dynamic') then
            WS = dynamic(C_, e, W)
        end if
        C0 = C
        C1 = C
        call perturb(C1, gamma)
        do iTStep = 1, nsteps
            ! J
            J0 = hartree(C0, dt / 2)
            J1 = hartree(C1, dt / 2)
            C0 = matmul(J0, C0)
            C1 = matmul(J1, C1)
            ! H
            C0 = matmul(H, C0)
            C1 = matmul(H, C1)
            ! K
            K0 = fockS(C0, dt, mode, WS, C0Prev, HQP, plasmon)
            K1 = fockS(C1, dt, mode, WS, C1Prev, HQP, plasmon)
            C0Prev = C0
            C1Prev = C1
            C0 = matmul(K0, C0)
            C1 = matmul(K1, C1)
            ! H
            C0 = matmul(H, C0)
            C1 = matmul(H, C1)
            ! J
            J0 = hartree(C0, dt / 2)
            J1 = hartree(C1, dt / 2)
            C0 = matmul(J0, C0)
            C1 = matmul(J1, C1)
            ! Autocorrelation Function
            P0 = 2 * matmul(C0, conjg(transpose(C0)))
            P1 = 2 * matmul(C1, conjg(transpose(C1)))
            dacf(iTStep) = sum([((P1(l,l) - P0(l,l)) * l, l = 1, n)]) / gamma
        end do
        do iOmega = 0, nOmega
            omega = omegas(iOmega)
            fourier = exp(Im * omega * ts) * exp(-ts**2*omegaInterval**2*2) * dt
            spectrum(iOmega, iMode) = abs(aimag(dot_product(fourier, dacf)))
        end do
    end do
    open(10, file='data/spectrum.dat')
    do iOmega = 0, nOmega
        write(10, '(5f10.3)') omegas(iOmega), spectrum(iOmega, :)
    end do
    close(10)
end program main
