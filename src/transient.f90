! Compute Monsoon index using Berger (1978) orbital solution
! for a selected transient time series. User defines initial
! and final time and a time step size. Output on stdout.
!
! Author: Maduvi
! Date: 28-07-2021
PROGRAM transient_mindex

    USE calormod

    IMPLICIT NONE

    INTEGER    :: i, niter
    INTEGER    :: operr
    REAL       :: t, tini, tfin, tinc
    REAL       :: ecc, xob, perh, prec
    REAL       :: perh0, ecc0, xob0, prec0
    REAL       :: latT, latE
    REAL       :: RsT0, RsE0, RsT, RsE, Raux
    REAL       :: dRsT, dRsE
    REAL       :: mT, mE, dQsT, dQsE, grad0, grad, mi
    REAL       :: QsT0, QsE0, QsT, QsE, deps
    CHARACTER(LEN=50)  :: arg1

    NAMELIST /INPUT_INFO/ tini, tfin, tinc
    
    ! latitudes
    latT = 23.45
    latE = 0

    ! get present day values of orbparams
    CALL BERGER(0.0, perh0, ecc0, xob0)
    prec0  = ecc0 * SIN(perh0 * pir)

    ! get present day caloric means
    CALL CALORIC_SUMMER(latT, ecc0, xob0, perh0, QsT0)
    CALL CALORIC_SUMMER(latE, ecc0, xob0, perh0, QsE0)
    grad0 = 2 * QsT0 - QsE0
    
    ! compute delta per degree in obliquity
    CALL ASTRO_SUMMER(latT, ecc0, xob0, perh0, RsT0)
    CALL ASTRO_SUMMER(latT, ecc0, xob0 + 1, perh0, Raux)
    dRsT = Raux - RsT0
    
    CALL ASTRO_SUMMER(latE, ecc0, xob0, perh0, RsE0)
    CALL ASTRO_SUMMER(latE, ecc0, xob0 + 1, perh0, Raux)
    dRsE = Raux - RsE0

    ! Read namelist file name with input variables from command line
    CALL GET_COMMAND_ARGUMENT(1, arg1)
    
    ! read namelist file
    OPEN (unit=1, file=TRIM(arg1), status='old', iostat=operr)
    IF(operr>0) THEN
        WRITE(*,'(A)') "mindex_b78: error: could not find namelist file."
        CALL EXIT(0)
    END iF
    READ (1, nml=INPUT_INFO)

    ! number of iterations
    niter = INT((tfin - tini) / tinc)

    WRITE(*, '(A10, A12, A12, A10, A10, A10, A10, A10, A10)') 'YEAR', 'QSTROP', &
        &'QSEQUA', 'MINDEX', 'PRECC', 'ECC', 'XOB', 'PERH'
    ! iterate for every kyear
    DO i = 0, niter        

        ! get time
        t = tini + tinc * i
        
        ! orbital parameters
        CALL BERGER(t, perh, ecc, xob)

        ! delta epsilon
        deps = xob - xob0

        ! delta precession index
        prec  = (ecc * SIN(perh * pir)) - prec0

        ! coefficient m
        mT = 2 * Tmin * Slang * COS(latT * pir) / (pi**2 * SQRT(1 - ecc**2))
        mE = 2 * Tmin * Slang * COS(latE * pir) / (pi**2 * SQRT(1 - ecc**2))

        ! change in Q following milankovitch
        dQsT = (dRsT * deps - prec * mT) / hyear
        dQsE = (dRsE * deps - prec * mE) / hyear

        ! absolute caloric summer
        QsT = QsT0 + dQsT
        QsE = QsE0 + dQsE

        ! compute monsoon index
        grad = 2 * QsT - QsE
        mi = grad - grad0

        WRITE(*, '(F10.2, F12.2, F12.2, F10.5, F10.5, F10.6, F10.5, F10.5)') t,&
            &QsT, QsE, mi, prec, ecc, xob, perh
    END DO
END PROGRAM transient_mindex
