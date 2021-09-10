! Compute Monsoon index using Berger (1978) orbital solution
! for a selected transient time series. User defines initial
! and final time and a time step size. Output on stdout.
!
! Author: Maduvi
! Date: 28-07-2021
PROGRAM monsindex

    USE calormod

    IMPLICIT NONE

    INTEGER             :: i, niter, ll
    INTEGER             :: operr
    INTEGER             :: ktyp, kuni
    INTEGER             :: Ts
    REAL                :: t, tini, tfin, tinc
    REAL                :: ecc, xob, perh
    REAL                :: lat
    REAL                :: perh0, ecc0, xob0, prec0
    REAL                :: latT, latE
    REAL                :: Rs, Rw, RsT0, RsE0, Rsols
    REAL                :: grad0, grad, mi
    REAL                :: QsT0, QsE0, QsT, QsE, Qs
    REAL                :: qlat
    REAL                :: accum
    CHARACTER(LEN=100)  :: arg1

    NAMELIST /INPUT_INFO/ ktyp, kuni, tini, tfin, tinc, perh, ecc, xob, qlat

    ! Read namelist file name with input variables from command line
    CALL GET_COMMAND_ARGUMENT(1, arg1)
    
    ! read namelist file
    OPEN (unit=1, file=TRIM(arg1), status='old', iostat=operr)
    IF (operr>0) THEN
        WRITE(*,'(A)') "mindex_b78: error: could not find namelist file."
        CALL EXIT(0)
    END IF
    READ (1, nml=INPUT_INFO)
    
    ! latitudes
    latT = 23.45
    latE = 0.0

    ! get present day values of orbparams
    CALL BERGER(0.0, perh0, ecc0, xob0)
    prec0  = ecc0 * SIN(perh0 * pir)

    ! get present day caloric means
    CALL PRESENT_CALORIC_SUMMER(latT, ecc0, xob0, perh0, QsT0)
    CALL PRESENT_CALORIC_SUMMER(latE, ecc0, xob0, perh0, QsE0)
    grad0 = 2 * QsT0 - QsE0
    
    IF (ktyp .eq. 1) THEN
        ! number of iterations
        niter = INT((tfin - tini) / tinc)
        WRITE(*, 100) 'YEAR', 'QSTROP', 'QSEQUA', 'MINDEX',&
            & 'ECC', 'XOB', 'PERH'
        
        ! iterate for every kyear
        DO i = 0, niter        

            ! get time
            t = tini + tinc * i
            
            ! orbital parameters
            CALL BERGER(t, perh, ecc, xob)
            
            ! caloric summer
            CALL CALORIC_SUMMER(t, latT, perh, ecc, xob, prec0, xob0, QsT0, QsT)
            CALL CALORIC_SUMMER(t, latE, perh, ecc, xob, prec0, xob0, QsE0, QsE)
            
            ! compute monsoon index
            grad = 2 * QsT - QsE
            mi = grad - grad0

            IF (KUNI .eq. 0) THEN
                WRITE(*, 101) t, QsT, QsE, mi, ecc, xob, perh
            ELSE
                WRITE(*, 101) t, QsT * conv, QsE * conv, mi * conv, ecc, xob, perh
            END IF
        END DO

    ELSE IF (ktyp .eq. 0) THEN
        ! caloric summer
        CALL CALORIC_SUMMER(t, latT, perh, ecc, xob, prec0, xob0, QsT0, QsT)
        CALL CALORIC_SUMMER(t, latE, perh, ecc, xob, prec0, xob0, QsE0, QsE)
        
        ! compute monsoon index
        grad = 2 * QsT - QsE
        mi = grad - grad0
        
        WRITE(*, '(A)') 'Input:'
        WRITE(*, '(A)') ''
        WRITE(*, '(A, F12.8)') 'PERH: ', perh
        WRITE(*, '(A, F12.8)') 'ECC: ', ecc
        WRITE(*, '(A, F12.8)') 'XOB: ', xob
        WRITE(*, '(A)') '===================='
        WRITE(*, '(A)') 'Results:'
        WRITE(*, '(A)') ''

        IF (KUNI .eq. 0) THEN
            WRITE(*, '(A, F12.3)') 'QsT:         ', QsT
            WRITE(*, '(A, F12.3)') 'QsE:         ', QsE
            WRITE(*, '(A, F12.3)') 'Mons. index: ', mi
        ELSE
            WRITE(*, '(A, F12.3)') 'QsT:         ', QsT * conv
            WRITE(*, '(A, F12.3)') 'QsE:         ', QsE * conv
            WRITE(*, '(A, F12.3)') 'Mons. index: ', mi * conv
        END IF

    ELSE

        ! number of iterations
        niter = INT((tfin - tini) / tinc)
        WRITE(*, 102) 'YEAR', 'RASTRO', 'QCALOR', 'RSOLST'        
        ! iterate for every kyear
        DO i = 0, niter        

            ! get time
            t = tini + tinc * i
            
            ! orbital parameters
            CALL BERGER(t, perh, ecc, xob)

            ! different insolation values
            CALL ASTRO_SUMMER(qlat, ecc, xob, perh, Rs)
            CALL CALORIC_SUMMER(t, qlat, perh, ecc, xob, prec0, xob0, QsE0, Qs)
            CALL SUMMER_SOLSTICE(qlat, ecc, xob, perh, Rsols)

            WRITE(*, 103)   t, Rs / hyear, Qs, Rsols
            
        END DO
        
    END IF
    
100 FORMAT (A12, A12, A12, A10, A10, A10, A10, A10, A10)
101 FORMAT (F12.2, F12.2, F12.2, F10.5, F10.6, F10.5, F10.5)
102 FORMAT (A12, A12, A12, A12)
103 FORMAT (F12.2, F12.2, F12.2, F12.2)
    
END PROGRAM monsindex
