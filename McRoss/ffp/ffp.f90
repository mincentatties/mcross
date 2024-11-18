!------------------------------------------------------------------------
!************************************************************************

PROGRAM ffp

! Initialization (parameter setting) of FFP for unix and windows (gfortran)
!
!                                     JRS 16/1/20
!------------------------------------------------------------------------

IMPLICIT NONE

REAL(8), DIMENSION(500)  :: q(500), ff(500)
REAL(8)                  :: g

LOGICAL               :: pltflg = .TRUE.

INTEGER               :: n, ii
CHARACTER (LEN = 1)   :: ans, pathstring, slash
CHARACTER (LEN = 2)   :: rare(11)
CHARACTER (LEN = 3)   :: ion, atom, os, delete
CHARACTER (LEN = 255) :: readdata, wdr

DATA rare /'Ce','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb'/

!------------------------------------------------------------------------

    os = 'lin'; delete='rm '; slash='/'
    CALL get_environment_variable("PATH",pathstring)
	IF (pathstring /= '/') THEN
        os ='win'
        delete = 'del'
        slash = '\\'
	END IF

    CALL getcwd(wdr); wdr=TRIM(wdr)//slash

    q = 0.0                                !initiallize q-array
    n = 100                                !number of points
    
    PRINT '(A,$)', 'Input ion <Mn2+>'
    READ(*,'(A3)',IOSTAT=ii), ion

    IF (ion == '   ') ion = 'Mn2'
    IF (ion(3:3) == ' ') ion = TRIM(ion)//'+'
    atom = ion
    IF (ion(3:3) == '+') atom = ion(1:2) 

    IF (ANY(ion(1:2) == rare)) THEN
        PRINT*, 'Rare-earth: g-factor will be calculated'
        g = 0
    ELSE 
        PRINT '(A,$)', 'Input g-factor <2.0>'
        READ(*,'(A255)'), readdata
        IF (readdata(1:4) == '    ') THEN 
            g = 2.0 
        ELSE 
            READ(readdata,*) g
        END IF
    END IF

    PRINT '(A,$)', 'Plot data? <y>'
    READ(*,'(A1)') ans
    IF (ans == 'n'.OR.ans == 'N') pltflg = .FALSE.

    CALL formfactor(q,n,g,ion,pltflg,wdr,ff)

    PRINT '(A,$)', 'Save data? <n>'
    READ(*,'(A1)') ans
    IF (ans == 'y'.OR.ans == 'Y') THEN 
        PRINT*, 'Data saved to '//TRIM(wdr)//TRIM(atom)//'.ff'
    ELSE
        PRINT*, 'Data not saved'
        CALL SYSTEM(delete//' "'//TRIM(wdr)//TRIM(ion)//'.ff"')
    END IF

    IF (pltflg) CALL system(delete//' "'//TRIM(wdr)//'main.gnu"')
    

END PROGRAM ffp
!-------------------------------------------------------------------------
!*************************************************************************