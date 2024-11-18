!------------------------------------------------------------------------
!************************************************************************

SUBROUTINE formfactor(q_in,nq,g,ion,pltflg,path,ff)

!------------------------------------------------------------------
! Function to calculate the magnetic form factor (in the dipole
! approximation) from Brown values
!------------------------------------------------------------------

! called by ffp and MagSRO
!                                              JRS 29/1/19

!------------------------------------------------------------------------
!************************************************************************

IMPLICIT NONE

REAL(8), DIMENSION(*) :: q_in, ff
REAL(8), ALLOCATABLE  :: s_in(:), j0(:), j2(:), term(:,:)
REAL(8), PARAMETER    :: pi = acos(-1.0) 
REAL(8)               :: g, l, s, j, ffo(19), gg

INTEGER              :: atomtype = 0, iq, nq,ii, iplot

CHARACTER            :: ion*(*),atom*2,ans*1
CHARACTER (LEN = 2)  :: three(12),fourd(12),rare(12),actin(12)
CHARACTER (LEN = 3)  :: delete
CHARACTER (LEN = 255):: path, fname
LOGICAL              :: pltflg 
         
DATA three/'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','  ','  ','  '/
DATA fourd/'Y ','Zr','Nb','Mo','Ru','Rh','Pd','  ','  ','  ','  ','  '/
DATA rare /'Ce','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Pr'/
DATA actin/'U ','Np','Pu','Am','  ','  ','  ','  ','  ','  ','  ','  '/

ALLOCATE(s_in(nq), j0(nq), j2(nq), term(6,nq))

!------------------------------------------------------------------------

! Make Q value array if none is input 
    IF (q_in(2) == 0.0) q_in(:100) = (/ (FLOAT(iq)/5.0, iq = 1, nq) /)

! Deduce atom type

    IF (ion(3:3) == '+') THEN
        atom = ion(1:1)//' '
        ion  = ion(1:2)//' '
    ELSE
        atom = ion(1:2)
    END IF

    IF (ANY(atom == three) .OR. ANY(atom == fourd)) atomtype = 1
    IF (ANY(atom == rare )) atomtype = 2
    IF (ANY(atom == actin)) atomtype = 3

    IF (atomtype == 0) THEN
        PRINT*,'FORMFACTOR: Unknown atom type'
        RETURN
    END IF

    CALL ff_table(ion,ffo)

    IF (ffo(1).EQ.0.0) THEN
        PRINT*,'ERROR: no formfactor available for '//ion
        RETURN
    END IF

!------------------------------------------------------------------------
! Calculate formfactor

    s_in(:nq) = q_in(:nq)/(4.0 * pi)
    DO ii = 1,3
        term(ii,:nq)   = ffo(ii*2-1)*EXP(-ffo(ii*2)*(s_in(:nq)**2))
        term(ii+3,:nq) = ffo(ii*2+7)*(s_in(:nq)**2)*EXP(-ffo(ii*2+8)*(s_in(:nq)**2))
    END DO
    j0 = SUM(term(1:3,:nq), DIM = 1) + ffo(7)
    j2 = SUM(term(4:6,:nq), DIM = 1) + ffo(15)*(s_in(:nq)**2)
    IF (atomtype == 1 .OR. atomtype == 3) THEN
         IF (g == 0.0) g = 2
         ff(:nq) = j0(:nq) + ((1-(2.0/g))*j2(:nq))
    ELSE
        l = ffo(17)
        s = ffo(18)
        j = ffo(19)
        IF (l == 0 .AND. j == 0.0 .AND. s == 0.0) THEN
            PRINT*, 'FORMFACTOR: L, S and J all zero...'
            RETURN
        END IF
        IF (g == 0) g  = 1 + ((j*(j+1)) - (l*(l+1)) + (s*(s+1))) / (2*j*(j+1))
        ff(:nq)=(j0(:nq)*(2. - 2./g)) + ((j0(:nq) + j2(:nq))*((2./g) - 1))
    END IF
!------------------------------------------------------------------------
!  Write out formfactor data

    fname = TRIM(path)//TRIM(ion)//'.ff'
    OPEN(10, FILE = fname)
    WRITE(10, '(''# '',''Ion: '',A4,1X,''g = '',F6.4,1X)') TRIM(ion)//'+', g
    WRITE(10,'(A38)') '#   Q       S       j0      j2      ff'
    WRITE(10,'(A41)') '#----------------------------------------'
    WRITE(10,'(5(F8.4))') (q_in(ii),s_in(ii),j0(ii),j2(ii),ff(ii),ii = 1,nq)


!------------------------------------------------------------------------
!  Plot formfactor data

    IF (pltflg) THEN
        DO iplot = 0,1
            OPEN(10,FILE = TRIM(path)//'main.gnu')
            WRITE(10,*) 'set encoding iso_8859_1'
            WRITE(10,*) 'cd '//''''//TRIM(path)//''''
            IF (iplot == 1) WRITE(10,*) 'set terminal pdf enhanced color size 6, 6 fsize 14'
            IF (iplot == 1) WRITE(10,*) 'set output '//''''//TRIM(fname)//'.pdf'//''''
            WRITE(10,*) 'set tics in'
            WRITE(10,*) 'set size square 1,1'
            WRITE(10,*) 'set xlabel "s=4{/Symbol p}/{/Symbol l} (1/\305)" offset 0,0'
            WRITE(10,*) 'set ylabel "F(s), j1(s), j2(s)" offset 0,0'
            WRITE(10,*) 'set mxtics 5'
            WRITE(10,*) 'set mytics 5'
            WRITE(10,*) 'set title "Formfactor for '//TRIM(ion)//'+ ion"'
            WRITE(10,*) 'set style data lines'
            WRITE(10,'(A)') &
               'plot "'  //TRIM(ion)//'.ff'//'" using 2:5 t ''F(s) '' w lines ls 1 lw 3 lc "red"' // &
                    ', "'//TRIM(ion)//'.ff'//'" using 2:3 t ''j0(s)'' w lines ls 1 lw 3 lc "blue"' // &
                    ', "'//TRIM(ion)//'.ff'//'" using 2:4 t ''j2(s)''w lines ls 1 lw 3 lc "green"'            
            IF (iplot == 0) WRITE(10,*) 'pause -1 "Hit <CR> to end"'
            CLOSE(10)
            IF (iplot == 0) CALL system('gnuplot "'//TRIM(path)//'main.gnu"')
            IF (iplot == 1) CALL system('gnuplot "'//TRIM(path)//'main.gnu" 2>/dev/null')
            IF (iplot == 1) EXIT
            PRINT '(A,$)', 'Keep Harcopy? <n>'
            READ(*,'(A1)') ans
            IF (ans == 'y' .OR. ans == 'Y') CYCLE
            EXIT
        END DO
    ENDIF

    END SUBROUTINE formfactor
!------------------------------------------------------------------------
!************************************************************************