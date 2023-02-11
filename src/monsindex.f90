!
! Compute monsoon forcing index (Rossignol-Strick, 1983)
! using Berger (1978) orbital solution for selected transient
! time series. User defines initial and final time and a time
! step size. Output in console stdout.
!
! Author: Maduvi
! Date: 28-07-2021
!
PROGRAM monsindex

    USE calormod

    IMPLICIT NONE

    INTEGER             :: i, niter
    INTEGER             :: operr
    INTEGER             :: ktyp, kuni
    REAL(KIND=8)        :: Ts, Tw
    REAL(KIND=8)        :: t, tini, tfin, tinc
    REAL(KIND=8)        :: ecc, xob, perh
    REAL(KIND=8)        :: lat
    REAL(KIND=8)        :: perh0, ecc0, xob0, prec0
    REAL(KIND=8)        :: Rs, Rw, RsT0, RsE0, Ss, Sw
    REAL(KIND=8)        :: RsT, RwT, SsT, SwT
    REAL(KIND=8)        :: RsE, RwE, SsE, SwE
    REAL(KIND=8)        :: grad0, grad, mi
    REAL(KIND=8)        :: QsT0, QsE0, QsT, QsE, Qs
    REAL(KIND=8)        :: qlat, s0
    REAL(KIND=8)        :: accum
    CHARACTER(LEN=256)  :: arg, cwd, nlf
    REAL(KIND=8), PARAMETER      :: latT=23.45, latE=0.0
    REAL(KIND=8), DIMENSION(360) :: solarm

    NAMELIST /INPUT_INFO/ ktyp, kuni, tini, tfin, tinc, perh, ecc, xob, qlat, s0

    ! Read namelist file name with input variables from command line
    CALL GETCWD(cwd)
    CALL GET_COMMAND_ARGUMENT(1, arg)
    nlf = TRIM(cwd) // '/' // TRIM(arg)
    
    ! read namelist file
    OPEN (unit=1, file=nlf, status='old', iostat=operr)
    IF (operr>0) THEN
        WRITE(*,'(A)') "mindex: error: missing namelist file"
        CALL EXIT(0)
    END IF
    READ (1, nml=INPUT_INFO)
    
    ! get present day values of orbparams
    CALL BERGER(0.0D0, perh0, ecc0, xob0)
    prec0 = ecc0 * SIN(perh0 * pir)

    ! get present day caloric means
    CALL TODAY_CALSUMMER(s0, latT, ecc0, xob0, perh0, QsT0)
    CALL TODAY_CALSUMMER(s0, latE, ecc0, xob0, perh0, QsE0)
    grad0 = (2.0 * QsT0) - QsE0

    ! the following depends on namelist input
    IF (ktyp .eq. 0) THEN
        ! number of iterations
        niter = INT((tfin - tini) / tinc)
        WRITE(*, 100) 'year', 'qst', 'qse', 'mindex', 'ecc', 'xob', 'perh'
        
        ! iterate for every kyear
        DO i = 0, niter        
            ! get time
            t = tini + (tinc * i)
            
            ! orbital parameters
            CALL BERGER(t, perh, ecc, xob)
            
            ! caloric summer
            CALL CALSUMMER(s0, latT, perh, ecc, xob, prec0, xob0, QsT0, QsT)
            CALL CALSUMMER(s0, latE, perh, ecc, xob, prec0, xob0, QsE0, QsE)
            
            ! compute monsoon index
            grad = (2.0 * QsT) - QsE
            mi = grad - grad0

            IF (KUNI .eq. 1) THEN
                QsT = QsT * conv
                QsE = QsE * conv
                mi = mi * conv
            END IF
            
            WRITE(*, 101) t, QsT, QsE, mi, ecc, xob, perh
            
        END DO
        
    ELSE IF (ktyp .eq. 1) THEN
        ! number of iterations
        niter = INT((tfin - tini) / tinc)

        WRITE(*, 102) 'year', 'sastro', 'wastro', 'ssols', 'wsols', 'Ts', 'Tw', 'ecc', 'perh', 'xob'

        DO i = 0, niter
            ! get time
            t = tini + (tinc * i)
            
            CALL BERGER(t, perh, ecc, xob) ! orbital parameters

            CALL SEASONS_LENGTH(ecc, perh, Ts, Tw)  ! to average astroseasons

            ! different insolation values
            CALL SINSOL(s0, qlat, ecc, xob, perh, solarm, Rs, Rw, Ss, Sw)

            IF (KUNI .eq. 1) THEN
                Rs = Rs * conv
                Rw = Rw * conv
                Ss = Ss * conv
                Sw = Sw * conv
            END IF
            
            WRITE(*, 103) t, Rs / Ts, Rw / Tw, Ss, Sw, Ts, Tw, ecc, perh, xob
        END DO

    ELSE IF (ktyp .eq. 2) THEN

        CALL BERGER(tini, perh, ecc, xob) ! orbital parameters        

        ! caloric summer
        CALL CALSUMMER(s0, latT, perh, ecc, xob, prec0, xob0, QsT0, QsT)
        CALL CALSUMMER(s0, latE, perh, ecc, xob, prec0, xob0, QsE0, QsE)
        
        ! compute monsoon index
        grad = (2.0 * QsT) - QsE
        mi = grad - grad0

        ! also extra insolation data, seasons lengths

        CALL SEASONS_LENGTH(ecc, perh, Ts, Tw)  ! to average astroseasons
        
        CALL SINSOL(s0, latT, ecc, xob, perh, solarm, RsT, RwT, SsT, SwT)
        CALL SINSOL(s0, latE, ecc, xob, perh, solarm, RsE, RwE, SsE, SwE)
        CALL SINSOL(s0, qlat, ecc, xob, perh, solarm, Rs, Rw, Ss, Sw)

        IF (KUNI .eq. 1) THEN
            QsT = QsT * conv
            QsE = QsE * conv
            mi = mi * conv
            RsT = RsT * conv
            RsE = RsE * conv
            RwT = RwT * conv
            RwE = RwE * conv
            SsT = SsT * conv
            SsE = SsE * conv
            SwT = SwT * conv
            SwE = SwE * conv
            Rs = Rs * conv
            Rw = Rw * conv
            Ss = Ss * conv
            Sw = Sw * conv
        END IF
        
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Input (time and lat):        '
        WRITE(*, '(A, F12.2)') 'tini:            ',        tini
        WRITE(*, '(A, F12.2)') 'qlat:            ', qlat
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Orbital parameters:          '
        WRITE(*, '(A, F12.8)') 'perh:            ', perh
        WRITE(*, '(A, F12.8)') 'ecc :            ', ecc
        WRITE(*, '(A, F12.8)') 'xob :            ', xob
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Astronomical seasons qlat:   '
        WRITE(*, '(A, F12.3)') 'Summer (avg):    ', Rs / Ts
        WRITE(*, '(A, F12.3)') 'Winter (avg):    ', Rw / Tw
        WRITE(*, '(A, F12.3)') 'Summer solst.:   ', Ss
        WRITE(*, '(A, F12.3)') 'Winter solst.:   ', Sw
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Results monsoon index:       '
        WRITE(*, '(A, F12.3)') 'QsTro (yr/2 avg):', QsT
        WRITE(*, '(A, F12.3)') 'QsEqu (yr/2 avg):', QsE
        WRITE(*, '(A, F12.3)') 'Mindex:          ', mi
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Astronomical seasons tropic: '
        WRITE(*, '(A, F12.3)') 'Summer (total):  ', RsT
        WRITE(*, '(A, F12.3)') 'Winter (total):  ', RwT
        WRITE(*, '(A, F12.3)') 'Summer solst.:   ', SsT
        WRITE(*, '(A, F12.3)') 'Winter solst.:   ', SwT
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Astronomical seasons equator:'
        WRITE(*, '(A, F12.3)') 'Summer (total):  ', RsE
        WRITE(*, '(A, F12.3)') 'Winter (total):  ', RwE
        WRITE(*, '(A, F12.3)') 'Summer solst.:   ', SsE
        WRITE(*, '(A, F12.3)') 'Winter solst.:   ', SwE        
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Astronomical seasons lengths:'
        WRITE(*, '(A, F12.3)') 'Tsummer:         ', Ts
        WRITE(*, '(A, F12.3)') 'Twinter:         ', Tw

    ELSE IF (ktyp .eq. 3) THEN

        IF (ecc == 0 .and. perh == 0 .and. xob == 0) THEN
            WRITE(*, *) "mindex: warning: orbital parameters all 0"
        END IF
        
        ! caloric summer
        CALL CALSUMMER(s0, latT, perh, ecc, xob, prec0, xob0, QsT0, QsT)
        CALL CALSUMMER(s0, latE, perh, ecc, xob, prec0, xob0, QsE0, QsE)
        
        ! compute monsoon index
        grad = (2.0 * QsT) - QsE
        mi = grad - grad0

        ! also extra insolation data, seasons lengths

        CALL SEASONS_LENGTH(ecc, perh, Ts, Tw)  ! to average astroseasons
        
        CALL SINSOL(s0, latT, ecc, xob, perh, solarm, RsT, RwT, SsT, SwT)
        
        CALL SINSOL(s0, latE, ecc, xob, perh, solarm, RsE, RwE, SsE, SwE)
        
        CALL SINSOL(s0, latE, ecc, xob, perh, solarm, Rs, Rw, Ss, Sw)

        IF (KUNI .eq. 1) THEN
            QsT = QsT * conv
            QsE = QsE * conv
            mi = mi * conv
            RsT = RsT * conv
            RsE = RsE * conv
            RwT = RwT * conv
            RwE = RwE * conv
            SsT = SsT * conv
            SsE = SsE * conv
            SwT = SwT * conv
            SwE = SwE * conv
            Rs = Rs * conv
            Rw = Rw * conv
            Ss = Ss * conv
            Sw = Sw * conv
        END IF
        
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Input (orbparams and lat):   '
        WRITE(*, '(A, F12.8)') 'perh:            ',        perh
        WRITE(*, '(A, F12.8)') 'ecc :            ', ecc
        WRITE(*, '(A, F12.8)') 'xob :            ', xob
        WRITE(*, '(A, F12.2)') 'qlat:            ', qlat
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Astronomical seasons qlat:   '
        WRITE(*, '(A, F12.3)') 'Summer (avg):    ', Rs / Ts
        WRITE(*, '(A, F12.3)') 'Winter (avg):    ', Rw / Tw
        WRITE(*, '(A, F12.3)') 'Summer solst.:   ', Ss
        WRITE(*, '(A, F12.3)') 'Winter solst.:   ', Sw
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Results monsoon index:       '
        WRITE(*, '(A, F12.3)') 'QsTro (yr/2 avg):', QsT
        WRITE(*, '(A, F12.3)') 'QsEqu (yr/2 avg):', QsE
        WRITE(*, '(A, F12.3)') 'Mindex:          ', mi
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Astronomical seasons tropic: '
        WRITE(*, '(A, F12.3)') 'Summer (total):  ', RsT
        WRITE(*, '(A, F12.3)') 'Winter (total):  ', RwT
        WRITE(*, '(A, F12.3)') 'Summer solst.:   ', SsT
        WRITE(*, '(A, F12.3)') 'Winter solst.:   ', SwT
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Astronomical seasons equator:'
        WRITE(*, '(A, F12.3)') 'Summer (total):  ', RsE
        WRITE(*, '(A, F12.3)') 'Winter (total):  ', RwE
        WRITE(*, '(A, F12.3)') 'Summer solst.:   ', SsE
        WRITE(*, '(A, F12.3)') 'Winter solst.:   ', SwE        
        WRITE(*, '(A)')        '-----------------------------'
        WRITE(*, '(A)')        'Astronomical seasons lengths:'
        WRITE(*, '(A, F12.3)') 'Tsummer:         ', Ts
        WRITE(*, '(A, F12.3)') 'Twinter:         ', Tw
        
    ELSE
        WRITE(*, *) "mindex: ktyp must be 0, 1, 2 or 3"
        STOP
    END IF
    
100 FORMAT (A12, A12, A12, A10, A10, A10, A10, A10, A10)
101 FORMAT (F12.2, F12.5, F12.5, F10.5, F10.6, F10.5, F10.5)
102 FORMAT (A12,   A12,   A12,   A12,   A12,   A12,   A12,   A12,   A12,   A12)
103 FORMAT (F12.2, F12.4, F12.4, F12.4, F12.4, F12.2, F12.2, F12.4, F12.4, F12.4)
    
END PROGRAM monsindex
