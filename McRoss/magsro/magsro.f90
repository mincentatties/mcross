!------------------------------------------------------------------------
!************************************************************************

PROGRAM MagSRO 

! RMC procedure to fit normalised, corrected, and radially
! averaged magnetic diffuse structure factors in binary systems
! At present, only deals with orthogonal structures

! input files required:
! 1) Raw data file        -> xye format
! 2) Atomic position file -> .pos file from NucSRO
! 3) Input parameters file

! output files:
! 1) .fit file             -> the computed structure factor
! 2) .mos file             -> magnetic atoms in the box with moment components
! 3) spinvert "spins" file -> for use with SCATTY and SPINCORREL
! 4) .cor file             -> spin correlations vs. distance 

!                           Ross Stewart: Created       20/5/98
!                                         Last Modified 26/8/21

! gfortran version

!      EXTERNALS: McRoss, formfactor, ff_table

!------------------------------------------------------------------------

USE McRoss
IMPLICIT NONE

REAL(8), PARAMETER    :: e_charge=1.602176487e-19 
REAL(8), PARAMETER    :: e_mass = 9.10938215e-31, gamma_n = -1.9130427
REAL(8), ALLOCATABLE  :: ft_compt(:,:), nft_compt(:,:), sumsinc(:), ff(:)
REAL(8), ALLOCATABLE  :: s_acc(:,:), pos_0(:,:)
REAL(8), ALLOCATABLE  :: a_ij(:), b_ij(:), d_aij(:), d_bij(:), cal(:), goodcal(:)

REAL(8)               :: sumdsq, olddsq, r_0, r, anchor(2), g, conc, sfac, j0(19), q
REAL(8)               :: max_spin_move = 0.2d0, norm, ds(3)
REAL                  :: tarray(2)

INTEGER, ALLOCATABLE      :: centre(:)
INTEGER                   :: stop_time,acceptance=0,mcount,n,nimp,nmag,numd,isim 
INTEGER                   :: ishell,icentre,imag,ispin,shnum,dim,ii,iq,check
INTEGER                   :: total_accept, cellpos, iatom
 
LOGICAL                   :: aflg  = .FALSE.
LOGICAL, ALLOCATABLE      :: mask(:)

CHARACTER (LEN = 1)       :: maga, ans
CHARACTER (LEN = 3)       :: ion
CHARACTER (LEN = 255)     :: simstr, spinstr, ffile, xtitle, ytitle, title

!------------------------------------------------------------------------
!  Title

    CALL get_os()
    CALL SYSTEM(clear)

    PRINT*
    PRINT*, '   *****************************************************'
    PRINT*, '   *  MM    MM   AAA   GGGGGGG SSSSSSS RRRRRR  OOOOOOO *'
    PRINT*, '   *  MM MM MM AA   AA GG      SS      RR   RR OO   OO *'
    PRINT*, '   *  MM    MM AAAAAAA GG  GGG SSSSSSS RRRRRR  OO   OO *'
    PRINT*, '   *  MM    MM AA   AA GG   GG      SS RR RR   OO   OO *'
    PRINT*, '   *  MM    MM AA   AA GGGGGGG SSSSSSS RR   RR OOOOOOO *'
    PRINT*, '   *****************************************************'
    PRINT*
    PRINT*, '       RMC Modelling of Magnetic Diffuse Scattering'
    PRINT*, '                     By Ross Stewart'
    PRINT*
    PRINT*, '           MacOSX GFortran version using GNUplot'
    PRINT*, '          ISIS, Rutherford Appleton Lab, Apr 2020' 
    PRINT*

!------------------------------------------------------------------------

    CALL getcwd(wdr); wdr=TRIM(wdr)//slash
    CALL magsro_run_params(ion,g,maga,aflg,anchor,dim)

    IF (maga == 'A') tmp = 'all'
    IF (maga == 'H') tmp = 'host'
    IF (maga == 'I') tmp = 'impurity' 

! Classical electron radius in 10^-12 cm units 
    r_0   = 1.e7 * (e_charge**2) / e_mass
    
! Read In Data and include anchor point
    CALL read_data(aflg,numd,anchor)
    ALLOCATE(ff(numd), sumsinc(numd), cal(numd))

! read in host and imp atomic position file, and allocate arrays
    CALL read_pos(n)
    ALLOCATE(s(3,n), s_acc(3,n), centre(n), mask(n), pos_0(3,n))
    pos_0 = pos                         ! keep all atoms for SPINVERT output

! assign spins to host, impurity or all atoms
    mask = imp > 0
    nimp = SIZE(PACK(imp,mask))
    conc = FLOAT(nimp)/FLOAT(n)
    mask = imp < 2                      ! set mask = .true. for all atoms
    IF (maga == 'H') mask = imp == 0    ! set mask = .true. for host atoms
    IF (maga == 'I') mask = imp == 1    ! set mask = .true. for impurity atoms
    nmag = UBOUND(PACK(imp,mask), DIM = 1)
    pos(1,1:nmag) = PACK(pos_0(1,:), mask)  ! remove all non-magnetic atoms
    pos(2,1:nmag) = PACK(pos_0(2,:), mask)
    pos(3,1:nmag) = PACK(pos_0(3,:), mask)

! Generate Form Factor
    CALL formfactor(qd,numd,g,ion,.FALSE.,wdr,ff)

! number of moves 
    maxmoves = maxmoves * nmag     ! convert to total moves
    check = maxmoves / 10

!------------------------------------------------------------------------
! Print out debugging messages
    IF (.NOT. aflg)    PRINT '(A)',     ' MagSRO: No Anchor at Q=0'
    IF (aflg) PRINT '(A,F6.3,A)',       ' MagSRO: Data anchored to ',anchor(1),' at Q = 0'
    PRINT '(A9,A3,A,F3.1)',             ' MagSRO: ',ion,'+ ions - g-factor = ',g
    IF (dim == 3) PRINT '(A)',          ' MagSRO: Using 3d-Heisenberg spins'
    IF (dim == 2) PRINT '(A)',          ' MagSRO: Using 2d-XY spins'                       
    PRINT '(A,F6.2)',                   ' MagSRO: Impurity concentration of ',conc
    PRINT '(A19,I5,A)',                 ' MagSRO: There are ',nimp,' impurity nuclei'
    PRINT '(A8,I2,A1,I2,A1,I2,A)',      ' MagSRO: ', nunitc(1)*2,'*',nunitc(2)*2,'*', & 
                                                     nunitc(3)*2,' unit cells'
    PRINT '(A)',                        ' MagSRO: Moments assigned to '//TRIM(tmp)//' atoms' 
    PRINT '(A,I5,A24)',                 ' MagSRO: ', nmag,' magnetic atoms in model'
    PRINT '(A,I3,A)',                   ' MagSRO: ', maxmoves/nmag,' moves per spin'
    PRINT '(A,F11.2,/)',                ' MagSRO: RMC weight = ',weight

!-----------------------------------------------------------------------
! sort magnetic sites into shells and calculate distances and r vectors
    CALL shellsort(nmag,shnum)
    ALLOCATE(ft_compt(shnum,numd), nft_compt(shnum,numd))
    ALLOCATE(a_ij(shnum), b_ij(shnum), d_aij(shnum), d_bij(shnum))

! calculate sinc terms in advance
    DO iq = 1, numd
        q = qd(iq)
        ft_compt(:,iq)  = SIN(q*sdist) / q / sdist
        nft_compt(:,iq) = SIN(q*sdist) / (q*sdist)**3 - COS(q*sdist) / (q*sdist)**2 
    END DO

    CALL ctrlc                     ! trap ctrlc command   

    !-----------------------LOOP OVER SIMULATIONS---------------------------
   
    DO isim = 1, numsims

        total_accept = 0
        ALLOCATE (chisq(maxmoves+1))

        PRINT '(/,A)', '*****************************************************************'
        PRINT '(1X,A,I2,A4,I2)',       ' MAGSRO: Simulation #',isim, ' of ', numsims
        PRINT '(A)',   '*****************************************************************'
    
        mcount = 0
        WRITE(simstr,'(A1,I0.2)') '_',isim
        spinstr = TRIM(ufile)//'_spins'//simstr 
        simstr  = TRIM(ufile)//simstr

        CALL init_random_seed()            ! initialise random numbers

    ! randomly assign spin directions    
        DO ispin = 1, nmag
            CALL ran_orient(s(:,ispin))
        END DO

        s_acc  = s
        olddsq = 1d10

        IF (idebug) PAUSE

        centre = [(icentre, icentre=1,nmag)]
        CALL blech_averbach(nmag,centre,nmag,shnum,a_ij,b_ij)

    !-----------------------LOOP OVER MC MOVES-------------------------------
        DO

    ! Calculate scattering
            DO iq = 1,numd
                    sumsinc(iq) = SUM(a_ij * ft_compt(:,iq) + b_ij* nft_compt(:,iq))
            END DO
            sumsinc = 2d0 * sumsinc / DBLE(nmag)
            cal = (gamma_n * r_0 / 2.)**2 * (ff**2) * (g**2) * (2./3. + sumsinc)      

    ! calculate final chisq, and model scattering function
            sfac   = SUM(cal(2:)*cod(2:)/errd(2:)**2) / SUM(cal(2:)**2/errd(2:)**2)
            ! avoids possible anchor point with small error
            sumdsq = SUM(((sfac*cal - cod)**2) / (errd**2))

    ! Take appropriate action depending on moves and chisq
            CALL RANDOM_NUMBER(r)
            IF (EXP(-weight*(sumdsq-olddsq)/2d0) > r) THEN
                s_acc(:,icentre) = s_acc(:,icentre) + s(:,icentre)
                d_aij = 0.0; d_bij = 0.0
                olddsq = sumdsq
                goodcal = cal
                chisq(mcount+1) = sumdsq
                acceptance   = acceptance + 1
                total_accept = total_accept + 1
            END IF

    ! Report progress and update plot     
            IF (MOD(mcount,check) == 0 .AND. mcount /= 0) THEN
                PRINT '(/,A,I3,A)', '---- ',NINT(100d0 * mcount / maxmoves),'% complete ----'
                PRINT '(A,I3,A,F7.1)', 'After ',mcount / nmag,' moves per spin, Chisq = ',olddsq
                PRINT '(A,F7.4)', 'S(S+1) = ', sfac
                IF (mcount /= 0) PRINT '(A,F6.2,A)', 'Acceptance rate = ', 100.0 * REAL(acceptance) / REAL(check), '%'
    
    ! write out fit and res file            
                ffile   = TRIM(simstr)//'.fit'
                OPEN(10,FILE = TRIM(wdr)//ffile)
                WRITE(10,'(3(F8.5,1X))') (qd(ii), sfac*goodcal(ii), cod(ii) - sfac*goodcal(ii), ii=1,numd)
                WRITE(10,'(A13,F8.4)') 'S(S+1) = ',sfac
                CLOSE(10)

    ! plot fit
                xtitle = 'Momentum Transfer (\305^{-1})'
                ytitle = 'Cross-Section (barns st^{-1} atom^{-1})'
                title  = 'Data and fit from '//ffile
                CALL gnuplot_fitmonitor(ffile,xtitle,ytitle,title,aflg,anchor)

                acceptance = 0
            END IF


    ! set spin orientations and coefficients to the most recently accepted values
            s = s_acc
            a_ij = a_ij - d_aij
            b_ij = b_ij - d_bij

            IF (stflg .OR. mcount == maxmoves) EXIT                                ! Stop program

    ! Rotate a random spin, to a random angle
            CALL RANDOM_NUMBER(r)
            icentre = INT(r * nmag) + 1                                            ! Choose spin
            CALL ran_orient(ds)                                                    ! find ds (spin nudge vector)
            s(:,icentre) = s(:,icentre) + max_spin_move * ds
            norm = NORM2(s(:,icentre))
            s(:,icentre) = s(:,icentre) / norm                                     ! new spin
            s(:,icentre) = s(:,icentre) - s_acc(:,icentre)                         ! calculate deltaS for randomly chosen spin
            CALL blech_averbach(nmag,[icentre],1,shnum,d_aij,d_bij)
            a_ij = a_ij + d_aij
            b_ij = b_ij + d_bij
            mcount = mcount + 1

        END DO                                       ! End of main Monte-Carlo Loop
       
        IF (stflg) EXIT                              ! Stop program before writing incomplete files

!------------------------------------------------------------------------
! calculate R-factor and find form factor cooefficients
        rfactor = 100.*SQRT(SUM(((sfac*goodcal - cod) / errd)**2) / SUM((cod / errd)**2))
        CALL ff_table(ion,j0)

! output chisq file
        chisq = PACK(chisq, MASK=chisq>0.0)
        OPEN(10,FILE = TRIM(wdr)//TRIM(simstr)//'.chi')
        WRITE(10,'(I7,1X,F7.1)') (ii, chisq(ii), ii=1,SIZE(chisq))
        CLOSE(10)

! output magnetic orientation file
        OPEN(10,FILE = TRIM(wdr)//TRIM(simstr)//'.mos')
        WRITE(10,'(A14,A,A17,I3)') " ! RMC fit of ",TRIM(dfile)," -  Simulation # ",isim
        WRITE(10,*) "!------------------------------------------------"
        WRITE(10,'(4(I3),A)') numat, nunitc*2, '            !  # atoms in unit cell; # unit cells in each direction' 
        WRITE(10,'(3(F7.4,1X),A)') latt(1), latt(2), latt(3), '   ! lattice constants'
        WRITE(10,'(A,F7.4)') ' ! S(S+1) scale factor - ', sfac
        WRITE(10,'(A,A3)') ' ! Magnetic ion: ', ion
        WRITE(10,*) "!"
        WRITE(10,*) "!    x            y           z          mux         muy         muz"
        WRITE(10,*) "!---------------------------------------------------------------------"
        WRITE(10,'(6(F11.5,1X))') (pos(1,imag),pos(2,imag),pos(3,imag),s(1,imag),s(2,imag),s(3,imag),imag=1,nmag)
        CLOSE(10)

!write out spinvert "spins" file        
        OPEN(10, FILE = TRIM(wdr)//TRIM(spinstr)//'.txt')
        WRITE(10,'(A)') 'TITLE   '//TRIM(ufile)
        WRITE(10,'(A,6(F11.5,1X)))')  'CELL   ', latt, 90.,90.,90.
        DO iatom = 1,numat
            WRITE(10,'(A,3(F11.5,1X)))') 'SITE   ', pos(1,iatom), pos(2,iatom),pos(3,iatom)
        END DO
        WRITE(10, '(A,3(I5))')       'BOX    ', nunitc*2
        WRITE(10, '(A,7(F11.5,1X))') 'FORM_FACTOR_J0   ', j0(:7)
        WRITE(10, '(A,I10)')         'PROPOSED_MOVES  ', maxmoves
        WRITE(10, '(A,I10)')         'ACCEPTED_MOVES  ', total_accept
        WRITE(10, '(A,F11.2)')       'WEIGHT  ', weight
        WRITE(10, '(A,F11.5)')       'CHI_SQUARED  ', olddsq
        WRITE(10, '(A,F11.5)')       'R_FACTOR  ', rfactor
        WRITE(10, '(A,F11.5)')       'SCALE   ', sfac * g**2
        WRITE(10, '(A,F11.5)')       'FLAT_BACKGROUND   ', 0.
        WRITE(10, '(A,F11.5)')       'LINEAR_BACKGROUND   ', 0.
        ! spinvert spin files need non-magnetic atoms ommited        
        cellpos = 0
        imag = 0
        DO iatom = 1,n                                                                              ! loop over all atoms
            cellpos = cellpos + 1                                                                   ! unit cell coordinate
            IF ((imp(iatom) == 0 .AND. maga /= 'I') .OR. (imp(iatom) == 1 .AND. maga /= 'H')) THEN  ! decide if magnetic
                imag = imag + 1                                                                     ! count magnetic atoms
                WRITE(10,'(A,4(I5),3(F11.5,1X))') 'SPIN   ', &
                    cellpos,INT(pos_0(1,iatom)),INT(pos_0(2,iatom)),INT(pos_0(3,iatom)),s(1,imag),s(2,imag),s(3,imag)
            END IF
            IF (cellpos == numat) cellpos = 0                                                       !reset cell coordinate index                           
        END DO
        CLOSE(10)

        DEALLOCATE(chisq)

    END DO                                          ! End of simulations loop

!------------------------------------------------------------------------

    IF (stflg) THEN
        PRINT '(/,A,/)','"Stop" command given' 
        numsims = isim - 1
    END IF

! calculate the spin-correlations
    IF (numsims > 0) CALL mag_correl(shnum)

! delete gnu file
    CALL SYSTEM(delete//' "'//TRIM(wdr)//'main.gnu"')

! flag the end of the calculation
    PRINT '(/,A,I3,A)', ' MAGSRO: ',numsims,' simulations completed'
    stop_time = ETIME(tarray)
    PRINT '(A,I2,A,I2,A)', ' Running time: ', stop_time / 60, ' mins ', MOD(stop_time, 60), ' secs'

END PROGRAM MagSRO


!------------------------------------------------------------------------
!************************************************************************

    SUBROUTINE ran_orient(m)

! caluclates random spin components of unit modulus for 2d or 3d
! From the method of Marsaglia, Ann.Math.Stat. 43 (1972) p645
! Called by MagSRO

    USE McRoss
    IMPLICIT NONE
    REAL(8)    :: m(3), rr, r1, r2, theta, pi

        pi = ACOS(-1.)
        CALL RANDOM_NUMBER(r1)
        CALL RANDOM_NUMBER(r2)
        m(3) = 2.*r1 -1.
        rr = SQRT(1.-m(3)*m(3))
        theta = 2.*pi*r2
        m(1) = rr*COS(theta)
        m(2) = rr*SIN(theta)
              
     END SUBROUTINE ran_orient

!------------------------------------------------------------------------
!************************************************************************

SUBROUTINE MagSRO_run_params(ion,g,maga,aflg,anchor,dim)

! Get input run parameters for MagSRO program

USE McRoss
IMPLICIT NONE

REAL(8)    :: g, anchor(*)
INTEGER             :: dim, ii, jj
CHARACTER (LEN = 1) :: maga, ans
CHARACTER (LEN = 3) :: ion
CHARACTER (LEN=300) :: line
LOGICAL             :: aflg

!------------------------------------------------------------------------
! defaults
    maga = 'A'
    dim  = 3
    g = 0.0

! Check for input file as argument
    ii = IARGC()
    IF (ii > 0) THEN
        CALL GETARG(1,parfile)
        OPEN(8,FILE = TRIM(wdr)//parfile, STATUS='OLD')
    ELSE
        PRINT '(A,$)', ' Input parfile name: '
        READ(*,'(A255)') parfile
        OPEN(8,FILE = TRIM(wdr)//parfile, STATUS='OLD')
    END IF

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
        ii = INDEX(line,'SUPERCELL')
        IF (ii /= 0) THEN
            READ(line(ii+9:),*) pfile
            pfile = TRIM(pfile)//'.pos'
            CYCLE
        END IF
        ii = INDEX(line,'MAG_ION')
        IF (ii /= 0) THEN
            READ(line(ii+7:),*) ion
            CYCLE
        END IF
        ii = INDEX(line,'G_FACTOR')
        IF (ii /= 0) THEN
            READ(line(ii+8:),*) g
            CYCLE
        END IF
        ii = INDEX(line,'MAGNETIC')
        IF (ii /= 0) THEN
            READ(line(ii+8:),*) maga
            IF (maga == 'h') maga = 'H'
            IF (maga == 'i') maga = 'I'
            IF (maga == 'a') maga = 'A'
            CYCLE
        END IF
        ii = INDEX(line,'DIM')
        IF (ii /= 0) THEN
            READ(line(ii+3:),*) ans
            IF (ans == 'H' .OR. ans == 'h') dim = 3
            IF (ans == 'X' .OR. ans == 'X') dim = 2
            IF (ans == 'I' .OR. ans == 'i') dim = 1
            CYCLE
        END IF
        ii = INDEX(line,'ANCHOR')
        IF (ii/= 0) THEN
            aflg = .TRUE.
            READ(line(ii+6:),*) anchor(:2)
            CYCLE
        ENDIF
        ii = INDEX(line,'WEIGHT')
        IF (ii/= 0) THEN
            READ(line(ii+6:),*) weight
            CYCLE
        ENDIF
        ii = INDEX(line,'MOVES')
        IF (ii/= 0) THEN
            READ(line(ii+5:),*) maxmoves
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
    END DO

22  CLOSE(8)

END SUBROUTINE MagSRO_run_params
!------------------------------------------------------------------------
!************************************************************************
