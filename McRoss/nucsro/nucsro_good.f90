!------------------------------------------------------------------------
!************************************************************************

PROGRAM NucSRO 

! RMC procedure to fit normalised, corrected, and radially averaged nuclear 
! diffuse structure factors in binary systems.
! At present, only deals with orthogonal structures.

! input files required:
! 1) simulation parameter file as argument
! 2) Raw data file        -> q,xs,err
!                            (xs is cross-section in barns/st/atom)
! output files:
! 1) .fit file             -> the computed structure factor
! 2) .chi file             -> the variation of ChiSq as a function of 
!                             the number of attempted moves
! 3) .pos file             -> output of the atoms in the box with host 
!                             atoms labelled "0" and impurity atoms .
!                             labelled "1". Used as input for magsro  
!                             program
! 4) .cow file             -> Warren-Cowley parameters as a function of
!                             shell distance
!
!                        Ross Stewart: Created       12/5/98
!                                      Last Modified 22/4/20
!
! EXTERNALS:  McRoss
!------------------------------------------------------------------------
!************************************************************************

USE McRoss
IMPLICIT NONE

REAL(8), ALLOCATABLE  :: prob(:),sumsinc(:),cal(:),wc(:),ft_compt(:,:)
REAL,    DIMENSION(2) :: tarray
INTEGER, ALLOCATABLE  :: hits(:), goodimp(:),newnb(:), oldnb(:), nb(:)

REAL(8)            :: sumdsq, anchor, conc, c, laue, olddsq, r, q, ii_level
INTEGER            :: shnum, nimp, stop_time, mcount, maxatom, isim, ishell
INTEGER            :: ii, iatom, iseed, icount, nx, ny, nz, iq
INTEGER            :: new_imp, new_host, n, numd, check, acceptance=0
LOGICAL            :: fflag, aflg, simflg, accept
CHARACTER(LEN=1)   :: ans
CHARACTER(LEN=50)  :: atime
CHARACTER(LEN=255) :: simstr, xtitle, ytitle, title, ffile

!------------------------------------------------------------------------
!  Title

    CALL get_os()
    CALL SYSTEM(clear)

    PRINT*
    PRINT*, '   *******************************************************'
    PRINT*, '   *  NN     NN UU    UU CCCCCCC SSSSSSS RRRRRR  OOOOOOO *'
    PRINT*, '   *  NNNN   NN UU    UU CC      SS      RR   RR OO   OO *'
    PRINT*, '   *  NN NN  NN UU    UU CC      SSSSSSS RRRRR   OO   OO *'
    PRINT*, '   *  NN  NN NN UU    UU CC           SS RR RR   OO   OO *'
    PRINT*, '   *  NN    NNN UUUUUUUU CCCCCCC SSSSSSS RR   RR OOOOOOO *'
    PRINT*, '   *******************************************************'
    PRINT*
    PRINT*, '        RMC Modelling Of Binary Diffuse Scattering'
    PRINT*, '                     By Ross Stewart'
    PRINT*
    PRINT*, '         ISIS, Rutherford Appleton Lab, March 2020'
    PRINT*

    CALL getcwd(wdr); wdr=TRIM(wdr)//slash
    CALL nucsro_run_params(conc,aflg,anchor,simflg,n,ii_level)

!------------------------------------------------------------------------
! Read In Data and include anchor point
    IF (.NOT. simflg) CALL read_data(aflg,anchor,numd)

! calculate initial guess of the Laue scattering
! average the scattering level and then subtract the ISO_INC
!    laue = SUM(cod) / numd - ii_level

! initialise arrays
    ALLOCATE(cal(numd), sumsinc(numd))
    cal=0.0; sumsinc=0.0
    ALLOCATE(imp(n), goodimp(n), hits(n))
    goodimp=0

! generate model by copying unit cell in x y and z directions
    icount = 0
    DO nx = 0, 2*nunitc(1)-1
        DO ny = 0, 2*nunitc(2)-1
            DO nz = 0, 2*nunitc(3)-1
                DO iatom = 1, numat
                    icount = icount + 1
                    pos(:,icount) = pos(:,iatom) + 1d0*[nx,ny,nz]
                END DO
            END DO 
        END DO
    END DO

    PRINT '(1X,A,I2,A,I2,A,I2,A)', ' NUCSRO: ', nunitc(1)*2,'*',nunitc(2)*2,'*',nunitc(3)*2,' unit cells'
    PRINT '(1X,A,3(1X,F6.3))',     ' NUCSRO: Lattice constants:',latt(1),latt(2),latt(3)                                   
    PRINT '(1X,A,I2,A,I5,A)',      ' NUCSRO: ',numat,' atoms per unit cell, and ',n, ' atoms in total'
    PRINT '(1X,A,F5.2,A1)',        ' NUCSRO: Impurity concentration of ',conc*100.,'%'
    PRINT '(1X,A,I5,A)',           ' NUCSRO: There are ',NINT(conc * FLOAT(n)),' impurity atoms'
    PRINT '(1X,A,F11.2)',          ' NUCSRO: RMC weight = ',weight
    PRINT*
    
!------------------------------------------------------------------------
! sort into shells, assign neighbours, and calculate shell sums
    shnum  = 0                    ! tell shellsort that shells are not determined
    CALL shellsort(n,shnum)
    ALLOCATE (prob(shnum),wc(shnum),newnb(shnum),oldnb(shnum),nb(shnum),ft_compt(shnum,numd))

! calculate sinc terms in advance
    IF (.NOT.simflg) THEN
        DO iq = 1, numd
            q = qd(iq)
            ft_compt(:,iq)  = SIN(q*sdist) / q / sdist
        END DO
    END IF
!------------------------------------------------------------------------
! start of multiple RMC simulations loop

    DO isim = 1, numsims
        mcount = 0                    ! initialise monte-carlo count
        nb=0; prob=0.0; wc=0.0; newnb=0; oldnb=0; hits=0
        olddsq = 0.0

        WRITE(simstr,'(A1,I0.2)') '_',isim
        simstr = TRIM(ufile)//simstr

        CALL init_random_seed()            ! initialise random numbers
        
    ! assign and count impurity atoms
        imp = 0
        nimp = 0
        DO
            CALL RANDOM_NUMBER(r)
            iatom = INT(r * FLOAT(n+1))
            IF (iatom == 0) CYCLE
            IF (imp(iatom) /= 1) THEN
                IF (nimp + 1 > NINT(conc * FLOAT(n))) EXIT
                imp(iatom) = 1
                nimp = nimp + 1
            END IF
        END DO
        hits = PACK( (/ (ii, ii=1,n) /), MASK = imp == 1)

        c = DBLE(nimp) / DBLE(n)

        IF (simflg) THEN
            goodimp=imp
            GO TO 90
        END IF

        PRINT '(/,A)', '*****************************************************************'
        PRINT '(1X,A,I2,A4,I2)',       ' NUCSRO: Simulation #',isim, ' of ', numsims
        PRINT '(A)', '*****************************************************************'

        CALL amplitude_sum(n,hits,nimp,imp,shnum,nb)   ! Nimp*N_B - total number of B atoms around B 

        IF (idebug) PAUSE

    !------------------------------------------------------------------------
    
        CALL ctrlc                    ! trap ctrlc command to interrupt

        check = maxmoves * nimp / 4
        ALLOCATE (chisq(maxmoves*nimp))

        DO mcount = 0, maxmoves*nimp        ! Main Monte-Carlo Loop

    ! Calculate scattering
            wc = (2d0 * nb / DBLE(nimp) - coord * c) / (1. - c)     !Z_n * alpha_n
            DO iq = 1,numd
                    cal(iq) = 1d0 + SUM(wc * ft_compt(:,iq))
            END DO
            laue   = SUM(cal*(cod-ii_level)/errd**2)/SUM(cal**2/errd**2)       ! from Joe's paper
            sumdsq = SUM(((laue*cal - (cod-ii_level))**2) / (errd**2))

    ! Take appropriate action depending on moves and chisq
            accept = .TRUE.
            CALL RANDOM_NUMBER(r)
            IF (EXP(-weight*(sumdsq-olddsq)/2.0) < r) accept = .FALSE.
            IF(accept .OR. mcount == 0) THEN
                goodimp  = imp
                oldnb = 0.0; newnb = 0.0
                olddsq   = sumdsq
                IF (mcount /= 0) chisq(mcount) = sumdsq
                acceptance = acceptance + 1
            END IF

    ! Report progress and update plot
            IF (MOD(mcount,check) == 0) THEN
                PRINT '(/,A,I3,A)', '---- ',NINT(100d0 * mcount / (maxmoves * nimp)) ,'% complete ----'
                PRINT '(A,I4,A,F7.1)', 'After ',mcount / nimp,' moves per impurity, Chisq = ',olddsq
                PRINT '(A,F7.4)', 'Laue level = ', laue
                IF (mcount /= 0) PRINT '(A,F5.2,A)', 'Acceptance rate =', 100.0 * DBLE(acceptance) / DBLE(check), '%'

    ! write out fit file
                ffile = TRIM(simstr)//'.fit'
                OPEN(10,FILE = TRIM(wdr)//ffile)
                WRITE(10,'(2(F8.5,1X))') (qd(ishell), laue*cal(ishell)+ii_level, ishell=1,numd)
                CLOSE(10)
    ! plot fit
                xtitle = 'Momentum Transfer (\305^{-1})'
                ytitle = 'Cross-Section (barns st^{-1} atom^{-1})'
                title  = 'Data and fit from '//ffile
                CALL gnuplot_fitmonitor(ffile,xtitle,ytitle,title)

                acceptance = 0
            END IF

    ! set impurity positions to the most recent accepted values and reset nb
            imp = goodimp
            nb = nb + oldnb - newnb

            IF (stflg) EXIT                             ! Stop program on ctrlc

    ! choose trial impurity and host atoms randomly
30          CALL RANDOM_NUMBER(r)
            new_host = INT(r * DBLE(n)) + 1
            IF (imp(new_host) /= 1) GO TO 30
40          CALL RANDOM_NUMBER(r)
            new_imp = INT(r * DBLE(n)) + 1
            IF (imp(new_imp) /= 0) GO TO 40

    ! calculate differences in amplitude sum due to swapped atoms
            CALL amplitude_sum(n,(/new_host/),1,imp,shnum,oldnb)    
            imp(new_host) = 0; imp(new_imp)  = 1
            CALL amplitude_sum(n,(/new_imp/),1,imp,shnum,newnb)    
            nb = nb - oldnb + newnb

        END DO                                        ! End of main Monte-Carlo Loop

!------------------------------------------------------------------------
        IF (stflg) EXIT                             ! Stop program on ctrlc

!       write out chi and position files

        chisq = PACK(chisq, MASK=chisq>0.0)
        OPEN(10,FILE = TRIM(wdr)//TRIM(simstr)//'.chi')
        WRITE(10,'(I7,1X,F7.1)') (ii, chisq(ii), ii=1,SIZE(chisq))
        CLOSE(10)

90      OPEN(10,FILE = TRIM(wdr)//TRIM(simstr)//'.pos')
        IF (simflg) WRITE(10,*) "! Nuclear position generator for MagSRO"
        IF (.NOT.simflg) WRITE(10,'(A14,A,A17,I3)') " ! RMC fit of ",TRIM(dfile)," -  Simulation # ",isim
        WRITE(10,*) "!------------------------------------------------"
        WRITE(10,'(4(I3),A)') numat, nunitc*2, '            !  # atoms in unit cell; # unit cells in each direction' 
        WRITE(10,'(3(F7.4,1X),A)') latt(1), latt(2), latt(3), '   ! lattice constants'
        WRITE(10,'(A2,I4,A,F6.3,A)') " !", nimp, ' impurities in model - concentration of ',c*100,'%'
        WRITE(10,'(A,F7.4)') ' ! "Laue" scale factor - ', laue
        WRITE(10,*) "!"
        WRITE(10,*) "!   x          y          z     imp"
        WRITE(10,*) "!--------------------------------------"
        WRITE(10,'(3(F10.5,1X),I2)') (pos(1,iatom),pos(2,iatom),pos(3,iatom),goodimp(iatom),iatom = 1, n)
        CLOSE(10)

        IF (simflg) THEN 
            PRINT*,'NUCSRO: positions generated' 
            CYCLE
        END IF

        DEALLOCATE(chisq)

    END DO                                            ! End of simulations loop

!------------------------------------------------------------------------

    IF (stflg) THEN
        PRINT '(/,A,/)','"Stop" command given' 
        numsims = isim - 1
    END IF

! Calculate WC parameters from outputed .pos files
    CALL nuc_correl(shnum)
    CALL SYSTEM(delete//' "'//TRIM(wdr)//'main.gnu"')

! Determine and display run-time
    IF (.NOT.simflg) PRINT '(/,A,I3,A)', ' NUCSRO: ',numsims,' simulations completed'
    stop_time = ETIME(tarray)
    PRINT '(A,I2,A,I2,A)', ' Running time: ', stop_time / 60, ' mins ', MOD(stop_time, 60), ' secs'

END PROGRAM NucSRO

!------------------------------------------------------------------------
!************************************************************************

SUBROUTINE NucSRO_run_params(conc,aflg,anchor,simflg,numtot,ii_level)

! Get input run parameters for NucSRO program

!------------------------------------------------------------------------

USE McRoss
IMPLICIT NONE

INTEGER             :: ii, box(3), numtot, k = 0
REAL(8)             :: laue, conc, anchor, ii_level
LOGICAL             :: aflg, simflg
CHARACTER(LEN=300)  :: line, ans

!------------------------------------------------------------------------

    aflg   = .FALSE.
    simflg = .FALSE.
    laue   = 0.0
    anchor = 0.0
    numat  = 0
    ii_level = 0

! Check for input file as argument
    ii = IARGC()
    IF (ii > 0) THEN
        CALL GETARG(1,parfile)
        OPEN(8,FILE = TRIM(wdr)//parfile, STATUS='OLD')
    ELSE
        PRINT '(A,$)', ' Input parfile name + ext: '
        READ(*,*) parfile
        OPEN(8,FILE = TRIM(wdr)//parfile, STATUS='OLD')
    ENDIF

    DO
        READ(8,'(A)',END=23) line
        ii = INDEX(line,'SITE')
        IF (ii /= 0) numat = numat + 1
    END DO
23  CLOSE(8)

    OPEN(8,FILE = TRIM(wdr)//parfile, STATUS='OLD')
    DO
        READ(8,'(A)',END=22) line
        line = TRIM(line) 
        IF (line(1:1) == '!') CYCLE

        ii = INDEX(line,'DATA')
        IF (ii /= 0) THEN
            READ(line(ii+4:),*) dfile
            dfile = TRIM(dfile)//'.dat'
            CYCLE
        END IF
        ii = INDEX(line,'SIMULATE')
        IF (ii /= 0) THEN
            simflg = .TRUE.
            CYCLE
        END IF
        ii = INDEX(line,'B_FRACTION')
        IF (ii /= 0) THEN
            READ(line(ii+10:),*) conc
            CYCLE
        END IF
        ii = INDEX(line, 'CELL')
        IF (ii /= 0) THEN
            READ(line(ii+4:),*) latt(1), latt(2), latt(3)
            CYCLE
        ENDIF
        ii = INDEX(line,'BOX')
        IF (ii/= 0) THEN
            READ(line(ii+3:),*) box(1), box(2), box(3)
            nunitc =INT(box/2)
            CYCLE
        ENDIF
        ii = INDEX(line,'ANCHOR')
        IF (ii/= 0) THEN
            aflg = .TRUE.
            READ(line(ii+6:),*) anchor
            CYCLE
        ENDIF
        ii = INDEX(line,'WEIGHT')
        IF (ii/= 0) THEN
            READ(line(ii+6:),*) weight
            CYCLE
        ENDIF        
        ii = INDEX(line,'RUNS')
        IF (ii/= 0) THEN
            READ(line(ii+4:),*) numsims
            CYCLE
        ENDIF
        ii = INDEX(line,'OUTPUT')
        IF (ii/= 0) THEN
            READ(line(ii+6:),*) ufile
            CYCLE
        ENDIF
        ii = INDEX(line,'SWAPS')
        IF (ii/= 0) THEN
            READ(line(ii+5:),*) maxmoves
            CYCLE
        ENDIF
        ii = INDEX(line,'ISO_INC')
        IF (ii/= 0) THEN
            READ(line(ii+7:),*) ii_level
            CYCLE
        ENDIF
        
    END DO
22  CLOSE(8)

    numtot = numat * PRODUCT(nunitc) * 8
    ALLOCATE(pos(3,numtot))
    pos=0.0
    
    OPEN(8,FILE = TRIM(wdr)//parfile, STATUS='OLD')
    DO
        READ(8,'(A)',END=24) line
        ii = INDEX(line,'SITE')
        IF (ii /= 0) THEN
            k = k+1
            READ(line(ii+4:),*) pos(1,k),pos(2,k),pos(3,k)
        END IF
    END DO
24  CLOSE(8)

END SUBROUTINE NucSRO_run_params

!------------------------------------------------------------------------
!************************************************************************
