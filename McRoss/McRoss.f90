!------------------------------------------------------------------------
!************************************************************************
MODULE McRoss

! Module containing common variables and procedures for McRoss programs
!                                                    JRS 26/8/21

!------------------------------------------------------------------------

IMPLICIT NONE

LOGICAL              :: iprint = .TRUE., idebug = .FALSE.
! iprint  -> turns on extra print statements
! idebug  -> outputs shells and distances debugging files, and inserts pauses
REAL(8), ALLOCATABLE :: pos(:,:),s(:,:),chisq(:),dist(:,:),r_unit(:,:,:)
REAL(8), ALLOCATABLE :: qd(:),cod(:),errd(:),sdist(:),coord(:)
INTEGER, ALLOCATABLE :: imp(:),shell_from_r(:)
REAL(8)              :: shellbin = 0.0001, latt(3), weight, rfactor
INTEGER              :: numat, maxmoves, numsims, nunitc(3)
CHARACTER (LEN = 255):: dfile, pfile, ufile, parfile, wdr, tmp, fmt
CHARACTER (LEN = 3)  :: os, delete, ptsize
CHARACTER (LEN = 5)  :: clear
CHARACTER (LEN = 1)  :: slash
LOGICAL              :: stflg

CONTAINS

!------------------------------------------------------------------------
    SUBROUTINE read_data(aflg,numd,anchor)

! Reads data file (X Y E) format.  Called by NucSRO and MagSRO.

    IMPLICIT NONE
    REAL(8), OPTIONAL :: anchor(*)
    REAL(8) :: dummy
    INTEGER :: idata, ii, k = 0, numd
    LOGICAL :: aflg

! count lines
        numd = 0
        OPEN(10, FILE = TRIM(wdr)//dfile, STATUS = 'OLD')
        DO
            READ(10,*, IOSTAT = ii) dummy, dummy, dummy
            IF (ii /= 0) EXIT
            numd = numd + 1
        END DO
        CLOSE(10)

        IF (aflg) k = 1
        ALLOCATE(qd(numd + k), cod(numd + k), errd(numd + k))

!read data
        OPEN(10, FILE = TRIM(wdr)//dfile, STATUS = 'OLD')
        DO idata = 1 + k, numd + k
            READ(10,*) qd(idata), cod(idata), errd(idata)
        END DO
        CLOSE(10)

        IF (errd(numd) == 0.0) errd(numd) = errd(numd-1)    !avoids odd zeroes if no carriage returns

! add anchor point at Q=0 if required
        IF (aflg) THEN 
            qd(1)   = 1.D-6
            cod(1)  = anchor(1)
            errd(1) = anchor(1) * anchor(2)
            numd = numd + 1
        END IF
    
    IF (iprint) PRINT*, 'READ_DATA: finished'

    END SUBROUTINE read_data

!------------------------------------------------------------------------
    SUBROUTINE read_pos(n)

! Reads atomic positions from .pos file generated from NucSRO.  
! Called by MagSRO and wc_calc

    IMPLICIT NONE
    INTEGER :: n, iatom, ii, shnum, ishell

        OPEN(10,FILE = TRIM(wdr)//pfile, STATUS = 'OLD')
        READ(10,*);READ(10,*)
        READ(10,*) numat, nunitc
        nunitc = nunitc / 2
        n = numat * PRODUCT(nunitc) * 8
        READ(10,*) latt
        READ(10,*);READ(10,*);READ(10,*);READ(10,*);READ(10,*)

        IF (.NOT.ALLOCATED(pos)) ALLOCATE(pos(3,n), imp(n))
        DO iatom = 1, n
            READ(10,*,IOSTAT=ii) pos(1,iatom),pos(2,iatom),pos(3,iatom),imp(iatom)
        END DO

        CLOSE(10)
        IF (idebug) PRINT*, 'READ_POS: finished'

    END SUBROUTINE read_pos

!------------------------------------------------------------------------
    SUBROUTINE read_mos(n)

! Reads atomic positions from .mos file generated from MagSRO, called by mag_correl()  

    IMPLICIT NONE
    INTEGER :: n, iatom, ii, shnum, ishell 

        OPEN(10,FILE = TRIM(wdr)//pfile, STATUS = 'OLD')
        READ(10,*);READ(10,*)
        READ(10,*) numat, nunitc
        nunitc = nunitc / 2
        n = numat * PRODUCT(nunitc) * 8
        READ(10,*) latt
        READ(10,*);READ(10,*);READ(10,*);READ(10,*);READ(10,*)

        IF (.NOT.ALLOCATED(pos)) ALLOCATE(pos(3,n), s(3,n))
        DO iatom = 1, n
            READ(10,*,IOSTAT=ii) pos(1,iatom),pos(2,iatom),pos(3,iatom),s(1,iatom),s(2,iatom),s(3,iatom)
            IF (ii /= 0) EXIT
        END DO
        n = iatom - 1

        CLOSE(10)
        IF (idebug) PRINT*, 'READ_MOS: finished'

    END SUBROUTINE read_mos

!-------------------------------------------------------------------------
    SUBROUTINE shellsort(n,shnum)

! Sorts atoms into shells and creates distance and lookup arrays
! Called by MagSRO, NucSRO, nuc_correl and mag_correl

    IMPLICIT NONE
    INTEGER, ALLOCATABLE :: n_atoms_r(:)
    REAL(8), ALLOCATABLE :: d1(:)
    REAL(8) :: r_max, a_i(3), a_j(3), r_ij(3), d_ij, min_d=0.0
    REAL    :: eps = EPSILON(1.0)
    INTEGER :: n, n_rbins, i, j, shnum, rbin, ishell=0, iatom, counter=0

        r_max = MINVAL(nunitc * latt) - eps
        n_rbins = NINT(r_max / shellbin) + 2

        ALLOCATE(dist(n,n), r_unit(3,n,n), n_atoms_r(n_rbins), d1(numat))

        n_atoms_r = 0; dist = 0d0; r_unit=0d0

! double loop through atoms - find distances, radial bins and count atoms
       DO i = 1, n
            a_i = pos(:,i)
            DO j = 1, i
                IF (i == j) CYCLE
                a_j = pos(:,j)
                r_ij = a_j - a_i

                ! boundary conditions
                WHERE (r_ij >   1d0*nunitc) r_ij = r_ij - 2d0*nunitc
                WHERE (r_ij <= -1d0*nunitc) r_ij = r_ij + 2d0*nunitc

                ! real distances
                r_ij = r_ij * latt
                d_ij = SQRT(DOT_PRODUCT(r_ij,r_ij))
                dist(j,i) = d_ij
                dist(i,j) = d_ij
                r_unit(:,i,j) = r_ij / d_ij
                r_unit(:,j,i) = -1d0 * r_unit(:,i,j)
                IF (d_ij > r_max) CYCLE

                rbin = NINT(d_ij / shellbin) + 1                ! assign bin for distance
                n_atoms_r(rbin) = n_atoms_r(rbin) + 2           ! two atoms ij = ji
            END DO
        END DO

        shnum = COUNT(n_atoms_r /= 0)

! check for non-equivalent sites and inform user
        DO iatom = 1, numat
            d1(iatom) = MINVAL(dist(iatom,:), MASK = dist(iatom,:) > eps)
        END DO
        DO iatom = 1,4
            min_d   = MINVAL(d1,MASK = d1 > min_d + eps)
            PRINT '(1X,A,I1,A,F7.4,A,I2,A)', &
            'SHELLSORT: Site ',iatom,', NN distance ',min_d,'Å - ',COUNT(ABS(d1-min_d)<eps), '-fold'
            counter = counter + COUNT(ABS(d1-min_d) < eps)
            IF (counter == numat) EXIT
            IF (counter > numat) THEN
                PRINT*, 'SHELLSORT: ERROR - problem with non-equivalent positions in cell'
                STOP
            END IF
        END DO
        IF (iatom > 1) PRINT*, 'SHELLSORT: more than one non-equivalent found in cell - expect non-integer coordination'

! write shell distance and coordiation arrays, and shell from r array
        ALLOCATE(coord(shnum),sdist(shnum),shell_from_r(n_rbins))

        IF (idebug) OPEN(10,FILE='shells.dat')
        DO rbin = 1, n_rbins
            IF (n_atoms_r(rbin) /= 0) THEN
                ishell = ishell + 1
                sdist(ishell) = (DBLE(rbin) - 1d0) * shellbin
                shell_from_r(rbin) = ishell
                coord(ishell) = DBLE(n_atoms_r(rbin)) / DBLE(n)
                IF (idebug) WRITE(10,*) sdist(ishell), coord(ishell)
            END IF
        END DO
        IF (idebug) CLOSE(10)
        
        DEALLOCATE(n_atoms_r, d1)
        PRINT '(1X,A,I4,A)',   'SHELLSORT: Using ',shnum,' shells'

    END SUBROUTINE shellsort

!------------------------------------------------------------------------
    SUBROUTINE amplitude_sum(atoms,centre,ncentres,imp,shnum,nb)

!  Calculates the number (nb) of B atoms around B centre.
!  Calculates half the total number since p_ij = p_ji.

    IMPLICIT NONE

    INTEGER, DIMENSION(*) :: nb, imp, centre
    INTEGER :: shnum, j, i, ish, iimp, atoms, rbin, ncentres

        nb(:shnum) = 0
        DO iimp = 1, ncentres
            i = centre(iimp)
            DO j = 1, atoms
                IF (i == j) CYCLE                    ! don't count i=j
                IF (j < i .AND. ncentres > 1) CYCLE  ! ij == ji
                IF (dist(i,j) > sdist(shnum) + shellbin) CYCLE

                rbin = NINT(dist(i,j) / shellbin) + 1   
                ish = shell_from_r(rbin)
                IF (imp(j) == 1) THEN
                    nb(ish) = nb(ish) + 1
                END IF
            END DO
        END DO

    END SUBROUTINE amplitude_sum

!------------------------------------------------------------------------
    SUBROUTINE blech_averbach(atoms,centre,ncentres,shnum,a_ij,b_ij,dot,nat)

! Calculates spin correlation functions for shnum near neighbour shells.
! Calculates half the total number since a_ij = a_ji, etc.

    IMPLICIT NONE
    REAL(8), DIMENSION(*), INTENT(out)  :: a_ij, b_ij
    REAL(8), DIMENSION(*), INTENT(out), OPTIONAL :: dot
    REAL(8)  :: rij(3), si(3), sj(3), si_sj, si_sj_z
    INTEGER, DIMENSION(*), INTENT(in) :: centre
    INTEGER, DIMENSION(*), OPTIONAL, INTENT(out) :: nat
    INTEGER, INTENT(in)  :: shnum, ncentres, atoms
    INTEGER  i, j, imag, ish, rbin

        a_ij(:shnum) = 0.0; b_ij(:shnum) = 0.0
        IF (PRESENT(dot)) THEN
            dot(:shnum) = 0.0; nat(:shnum)=0.0
        END IF
        
        DO imag = 1, ncentres
            i = centre(imag)
            DO j = 1, atoms
                IF (j == i) CYCLE                    ! don't count i=j
                IF (j < i .AND. ncentres > 1) CYCLE  ! ij == ji
                IF (dist(i,j) > sdist(shnum) + shellbin) CYCLE       ! maximum shell to consider

                rbin = NINT(dist(i,j) / shellbin) + 1
                ish  = shell_from_r(rbin) 

                si = s(:,i); sj = s(:,j); rij = r_unit(:,i,j)
                si_sj   = DOT_PRODUCT(si,sj)
                si_sj_z = DOT_PRODUCT(si,rij) * DOT_PRODUCT(sj,rij)

                a_ij(ish) = a_ij(ish) + si_sj - si_sj_z
                b_ij(ish) = b_ij(ish) + 3d0*si_sj_z - si_sj

                IF (PRESENT(dot)) THEN
                    dot(ish)  = dot(ish) + si_sj
                    nat(ish)  = nat(ish) + 1
                END IF
            END DO
        END DO

    END SUBROUTINE blech_averbach

!------------------------------------------------------------------------
    SUBROUTINE nuc_correl(shnum)

        IMPLICIT NONE
        REAL(8), ALLOCATABLE :: prob_av(:),prob_sdev(:),wc(:),wcerr(:),conc(:),prob_diff(:),prob_dist(:,:)
        INTEGER, ALLOCATABLE :: hits(:), nb(:)
        REAL(8)  :: c
        INTEGER  :: n, isim, shnum, ii, ishell, nimp
        CHARACTER(LEN=255) :: simstr, cfile, xtitle, ytitle, title

! read in .pos files and calculate distribution of shell dependent porbabilities
        DO isim = 1, numsims
            WRITE(simstr,'(A1,I0.2)') '_',isim
            simstr = TRIM(ufile)//simstr
            pfile = TRIM(simstr)//'.pos'
            CALL read_pos(n)
            IF (shnum == 0 .AND. isim == 1) CALL shellsort(n, shnum)
            IF (.NOT.ALLOCATED(hits)) THEN 
                ALLOCATE(hits(n),nb(shnum),prob_dist(numsims,shnum),prob_av(shnum),conc(numsims))
                ALLOCATE(prob_sdev(shnum),wc(shnum),wcerr(shnum),prob_diff(numsims))
            END IF
            hits = PACK( [(ii, ii=1,n)], MASK = imp == 1)
            nimp = SIZE(hits)
            conc(isim) = DBLE(nimp) / DBLE(n)
            CALL amplitude_sum(n,hits,nimp,imp,shnum,nb)
            prob_dist(isim,:) = 2. * DBLE(nb) / DBLE(nimp) / coord 
            IF (nimp == 0) prob_dist(isim,:) = 0
        END DO
        c = SUM(conc) / DBLE(numsims)

! calculate average and standard deviation 
        prob_av = SUM(prob_dist, DIM=1) / DBLE(numsims)
        DO ishell = 1, shnum
            prob_diff = prob_dist(:,ishell) - prob_av(ishell)
            prob_sdev(ishell) = 0.0
            IF (numsims > 1) prob_sdev(ishell) = SQRT( SUM(prob_diff**2) / DBLE(numsims-1) ) !see wikipedia!!
        END DO

! calcalate WC parameters from average probabilities
        wc   = (prob_av - c) / (1.0 - c)
        wcerr = prob_sdev / (1.0 - c)

! output Warren-Cowley parameters
        cfile = TRIM(ufile)//'.cow'
        title = 'Warren-Cowley parameters from '//TRIM(ufile)//'.cow'
        OPEN(10,FILE = TRIM(wdr)//cfile)
        WRITE(10,'(A)') "  Shell r/Å  WC param   error      Zn"
        WRITE(10,'(A)') "------------------------------------------"
        WRITE(10,'(3(F10.5,1X),F8.2)') (sdist(ishell),wc(ishell),wcerr(ishell),coord(ishell),ishell=1,shnum)
        CLOSE(10)

! plot Warren-Cowley parameters
        xtitle = 'Radial Distance (\305)'
        ytitle = '{/Symbol a}(0n)'
        CALL gnuplot_fitmonitor(cfile,xtitle,ytitle,title)

        PRINT*, 'NUC_CORREL: WC parameters calculated'

    END SUBROUTINE nuc_correl

!------------------------------------------------------------------------
    SUBROUTINE mag_correl(shnum)

    IMPLICIT NONE
    REAL(8), ALLOCATABLE :: dot_dist(:,:),dot(:),a_ij(:),b_ij(:),dot_av(:),dot_sdev(:),dot_diff(:)
    INTEGER, ALLOCATABLE :: centre(:), nat(:)
    INTEGER :: nmag, shnum, isim, icentre, ishell
    CHARACTER(LEN=255) :: simstr, cfile, title, xtitle, ytitle  

    DO isim = 1,numsims
        WRITE(simstr,'(A1,I0.2)') '_',isim
        simstr = TRIM(ufile)//simstr
        pfile = TRIM(simstr)//'.mos' 
        CALL read_mos(nmag)
        IF (shnum == 0 .AND. isim == 1) CALL shellsort(nmag, shnum)
        IF (.NOT.ALLOCATED(centre)) THEN
            ALLOCATE(centre(nmag),a_ij(shnum),b_ij(shnum),dot(shnum),dot_dist(numsims,shnum))
            ALLOCATE(dot_av(shnum),dot_sdev(shnum),dot_diff(numsims),nat(shnum))
        END IF
        centre = (/(icentre, icentre=1,nmag)/)
        CALL blech_averbach(nmag,centre,nmag,shnum,a_ij,b_ij,dot,nat)
        dot_dist(isim,:) = dot / REAL(nat)     ! normalised to nmag/2 since only i < j counted and Si.Sj = Sj.Si
    END DO

! calculate average and standard deviation 
    dot_av = SUM(dot_dist, DIM=1) / REAL(numsims)
    dot_sdev = 0.0
    DO ishell = 1, shnum
        dot_diff = dot_dist(:,ishell) - dot_av(ishell)
        IF (numsims > 1) dot_sdev(ishell) = SQRT( SUM(dot_diff**2) / REAL(numsims-1) )
    END DO

! output correlations
        cfile = TRIM(ufile)//'.cor'
        title = 'Spin-correlations from '//TRIM(ufile)//'.cor'
        OPEN(10,FILE = TRIM(wdr)//cfile)
        WRITE(10,'(A)') "  Shell r/Å  (S0.S1)   error"
        WRITE(10,'(A)') "--------------------------------"
        WRITE(10,'(3(F10.5,1X))') (sdist(ishell),dot_av(ishell),dot_sdev(ishell),ishell=1,shnum)
        CLOSE(10)

! plot Spin-correlations
        xtitle = 'Radial Distance (\305)'
        ytitle = '< S_{0}.S_{i} > / S(S+1)'
        CALL gnuplot_fitmonitor(cfile,xtitle,ytitle,title)

        PRINT*, 'MAG_CORREL: Magnetic correlations calculated'

    END SUBROUTINE mag_correl
!----------------------------------------------------------------------
!      
    SUBROUTINE gnuplot_fitmonitor(ffile,xtitle,ytitle,title,aflg,anchor)

    !  Generic gnuplot batch file operation routine developed for use with
    !  nucsro and magsro

    IMPLICIT NONE
    LOGICAL, OPTIONAL :: aflg
    REAL(8), OPTIONAL :: anchor(*)
    CHARACTER (LEN=255) :: ffile,xtitle,ytitle,title,input_pfile,ofile,anch_str=""

        ofile = TRIM(ffile)//'.pdf'

        IF(PRESENT(aflg)) THEN
            IF (aflg) WRITE(anch_str,*)', "<echo 0 ',anchor(1),'" ls 2 pt 7 ps '//TRIM(ptsize)//' lc "red"'
        END IF

        input_pfile = 'plot "'//TRIM(dfile)//'" using 1:2:3 ls 2 pt 7 ps '//TRIM(ptsize)//' lc "black", "' & 
                   //TRIM(ffile)//'" using 1:2 w lines ls 1 lw 3 lc "red", "'&
                   //TRIM(ffile)//'" using 1:3 w lines ls 1 lw 3 lc "blue"'//anch_str
        IF (INDEX(ffile,'.co') > 0) THEN 
            input_pfile = 'plot "'//TRIM(ffile)//'" using 1:2:3 pt 7 ps '//TRIM(ptsize)//' lc "black"'
        END IF

        OPEN(10,file = TRIM(wdr)//'main.gnu')
        WRITE(10,*) 'set terminal pdf enhanced color size 6, 6 fsize 14'
        WRITE(10,*) 'set encoding iso_8859_1'
        WRITE(10,*) 'set output "'//TRIM(ofile)//'"'
        WRITE(10,*) 'cd '//''''//wdr(:INDEX(wdr,'  ') - 2)//''''
        WRITE(10,*) 'set tics in'
        WRITE(10,*) 'set size square 1,1'
        WRITE(10,*) 'set xlabel "'//TRIM(xtitle)//'" offset 0,0'
        WRITE(10,*) 'set ylabel "'//TRIM(ytitle)//'" offset 0,0'
        WRITE(10,*) 'set mxtics 5'
        WRITE(10,*) 'set mytics 5'
        WRITE(10,*) 'set title  "'//TRIM(title) //'" offset 0,0 noenhanced'
        WRITE(10,*) 'unset key'
        WRITE(10,*) 'set border lw 2'
        WRITE(10,*) 'set errorbars small'
        WRITE(10,*) 'set style data error'
        WRITE(10,*) 'set xzeroaxis lw 2'
        WRITE(10,'(A)') TRIM(input_pfile)
        CLOSE(10)

        CALL system('gnuplot "'//TRIM(wdr)//'main.gnu" 2>/dev/null')

    END SUBROUTINE gnuplot_fitmonitor

!------------------------------------------------------------------------

    SUBROUTINE ctrlc

! Traps the ctrl-C event

    EXTERNAL :: ctrlc_ast
    INTEGER  :: signal, sigint, iret
        
        sigint = 2
        iret = SIGNAL(sigint, ctrlc_ast)

    END SUBROUTINE ctrlc
!------------------------------------------------------------------------

    SUBROUTINE init_random_seed()

    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
              
        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))
        CALL SYSTEM_CLOCK(COUNT=clock)
        seed = clock + 37 * [(i - 1, i = 1, n)]
        CALL RANDOM_SEED(PUT = seed)
        DEALLOCATE(seed)

    END SUBROUTINE init_random_seed

    !------------------------------------------------------------------------
    SUBROUTINE get_os()

    IMPLICIT NONE
    CHARACTER (LEN=1) :: pathstring
    
    os = 'lin'; delete='rm '; slash='/'; clear='clear'; ptsize='1.5'
    CALL get_environment_variable("PATH",pathstring)
    IF (pathstring /= '/') THEN
        os ='win'
        delete = 'del'
        slash  = '\\'
        clear  = 'cls'
        ptsize = '0.3'
    END IF

    END SUBROUTINE get_os
!------------------------------------------------------------------------
    
END MODULE McRoss

!------------------------------------------------------------------------
!************************************************************************

SUBROUTINE ctrlc_ast

USE McRoss
CHARACTER(LEN = 1) :: ans

    stflg = .TRUE.

END SUBROUTINE ctrlc_ast

!-------------------------------------------------------------------------


