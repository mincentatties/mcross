!---------------------------------------------------------------------
    PROGRAM cor_calc
!
! Calls mag_correl to calculation spin correlations from magsro mos files
! 
!                                               JRS 21/3/19
!---------------------------------------------------------------------
        USE McRoss
        IMPLICIT NONE

        INTEGER, DIMENSION(13) :: buff
        INTEGER                :: ii=0, isim=1, status, shnum
        CHARACTER(LEN=255)     :: simstr,f_try

        CALL get_os()
        CALL getcwd(wdr); wdr=TRIM(wdr)//slash

! Check for input file as argument
        ii = IARGC()
        IF (ii > 0) THEN
            CALL GETARG(1,ufile)
        ELSE
            PRINT '(/,A,/)', ' COR_CALC: Please supply filestem of simulations'
            STOP
        ENDIF

! Find and count .mos files 
        DO
            WRITE(simstr,'(A1,I0.2)') '_',isim
            f_try = TRIM(wdr)//TRIM(ufile)//TRIM(simstr)//'.mos'
            CALL STAT(f_try, buff, status)
            IF (status .NE. 0) EXIT
            isim = isim + 1
        END DO
        numsims = isim - 1

        IF (numsims == 0) THEN
            PRINT '(/,A,/)', ' COR_CALC: Error - no simulations found'
            STOP
        END IF

        PRINT '(/,A,I2,A,/)', ' COR_CALC: calculating spin correlations from ', numsims, ' configurations'

! Calculate correlations
        shnum = 0
        CALL mag_correl(shnum)
        CALL SYSTEM(delete//' "'//TRIM(wdr)//'main.gnu"')
        

    END PROGRAM cor_calc
!---------------------------------------------------------------------
