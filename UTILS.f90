!------------------------------------------------------------------------------------------------------------------------------------------
MODULE  UTILS
!------------------------------------------------------------------------------------------------------------------------------------------
    USE MPI
    IMPLICIT NONE
    !--------------------------------------------------------------------------------------------------------------------------------------
    CHARACTER*500                               ::   info_str  
    !--------------------------------------------------------------------------------------------------------------------------------------
    ! constant
    !--------------------------------------------------------------------------------------------------------------------------------------
    REAL*8                                      ::   bohr
    REAL*8                                      ::   pi
    !--------------------------------------------------------------------------------------------------------------------------------------
    ! input setting
    !--------------------------------------------------------------------------------------------------------------------------------------
    INTEGER                                     ::   flag_xc
    INTEGER                                     ::   n_periodic
    REAL*8                                      ::   mbd_supercell_cutoff  
    REAL*8                                      ::   mbd_cfdm_dip_cutoff 
    REAL*8                                      ::   mbd_scs_dip_cutoff
    LOGICAL,                  DIMENSION(3)      ::   mbd_scs_vacuum_axis
    !--------------------------------------------------------------------------------------------------------------------------------------
    ! geometry input
    !--------------------------------------------------------------------------------------------------------------------------------------
    CHARACTER*2, ALLOCATABLE, DIMENSION(:)      ::   atom_name
    REAL*8,      ALLOCATABLE, DIMENSION(:)      ::   hirshfeld_volume             ! hirshfeld_volume from hirshfeld partioning
    REAL*8,      ALLOCATABLE, DIMENSION(:,:)    ::   coords
    REAL*8,                   DIMENSION(3,3)    ::   lattice_vector
    !--------------------------------------------------------------------------------------------------------------------------------------
    INTEGER,     ALLOCATABLE, DIMENSION(:)      ::   task_list
    REAL*8,      ALLOCATABLE, DIMENSION(:,:)    ::   coupled_atom_pol
    REAL*8,      ALLOCATABLE, DIMENSION(:,:)    ::   relay_matrix                 ! 3N x 3N (imaginary) frequency dependent matrix
    REAL*8,      ALLOCATABLE, DIMENSION(:,:)    ::   relay_matrix_periodic
    REAL*8,      ALLOCATABLE, DIMENSION(:)      ::   alpha_omega
    REAL*8,      ALLOCATABLE, DIMENSION(:)      ::   R_p                          ! effective dipole spread of QHO, [PRB,75,045407(2007)]
    REAL*8,      ALLOCATABLE, DIMENSION(:)      ::   Rvdw_iso                             
    REAL*8,      ALLOCATABLE, DIMENSION(:)      ::   alpha_eff
    REAL*8,      ALLOCATABLE, DIMENSION(:)      ::   C6_eff
    REAL*8,                   DIMENSION(3,3)    ::   mol_pol_tensor
    REAL*8,                   DIMENSION(0:20)   ::   casimir_omega
    REAL*8,                   DIMENSION(0:20)   ::   casimir_omega_weight
    !--------------------------------------------------------------------------------------------------------------------------------------
    ! mpi
    !--------------------------------------------------------------------------------------------------------------------------------------
    INTEGER                                     ::   nproc                        ! number of processes
    INTEGER                                     ::   iproc                        ! process id
    INTEGER                                     ::   ierr                         ! error handle
    INTEGER                                     ::   n_atoms                      ! number of atoms
    !--------------------------------------------------------------------------------------------------------------------------------------
    ! scalapack
    !--------------------------------------------------------------------------------------------------------------------------------------
    INTEGER                                     ::   nprow                        ! processes grid along row
    INTEGER                                     ::   npcol                        ! processes grid along col
    INTEGER                                     ::   ipcol                        ! process id along row
    INTEGER                                     ::   iprow                        ! process id along col
    INTEGER                                     ::   icontxt                      ! blacs related
    INTEGER                                     ::   desc_A(9)                    ! description of A matrix
    INTEGER                                     ::   desc_B(9)                    ! description of B matrix
    !--------------------------------------------------------------------------------------------------------------------------------------

CONTAINS

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE init_constants()
        !----------------------------------------------------------------------------------------------------------------------------------
        flag_xc                = 1                                                ! XC= PBE
        mbd_supercell_cutoff   =  25.0D0                                          ! Angstrom
        mbd_cfdm_dip_cutoff    = 100.0D0                                          ! Angstrom
        mbd_scs_dip_cutoff     = 120.0D0                                          ! Angstrom
        n_periodic             = 0  
        bohr                   = 0.52917721D0
        pi                     = 3.14159265358979323846d0
        mbd_scs_vacuum_axis(1) = .FALSE. 
        mbd_scs_vacuum_axis(2) = .FALSE. 
        mbd_scs_vacuum_axis(3) = .FALSE. 
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE init_constants
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE start_mpi()
        !----------------------------------------------------------------------------------------------------------------------------------
        CALL MPI_INIT (ierr)
        CALL MPI_COMM_RANK (MPI_COMM_WORLD, iproc, ierr)
        CALL MPI_COMM_SIZE (MPI_COMM_WORLD, nproc, ierr)
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE stop_mpi()
        !----------------------------------------------------------------------------------------------------------------------------------
        CALL MPI_FINALIZE(ierr)
        IF(ierr.ne.0) WRITE(*,*) "ERROR in terminating MPI"
        IF(ierr.ne.0) WRITE(*,*) "Normal Termination"
        STOP
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE allocate_task()
        !----------------------------------------------------------------------------------------------------------------------------------
        INTEGER :: i_row
        !----------------------------------------------------------------------------------------------------------------------------------
        IF(.NOT. ALLOCATED(task_list)) ALLOCATE(task_list(n_atoms))
        DO i_row=1,n_atoms
            task_list(i_row) = MOD(i_row,nproc)
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE allocate_task
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE init_blacs()
        CALL BLACS_GET(0, 0, icontxt)
        CALL BLACS_GRIDINIT(icontxt, "Row", nprow, npcol)
        CALL BLACS_GRIDINFO(icontxt, nprow, npcol, iprow, ipcol)
        WRITE(*,"(A,I1,A,I1,A,I1,A,I1,A,I1)") "nproc=",nproc,"   iproc=",iproc,"   iprow=",iprow,"   ipcol=",ipcol,"   icontxt=",icontxt
    END SUBROUTINE init_blacs
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE sync_tensors(temp_matrix,num_element)
        !----------------------------------------------------------------------------------------------------------------------------------
        INTEGER                                      :: num_element
        REAL*8, DIMENSION(num_element,num_element)   :: temp_matrix
        REAL*8, DIMENSION(num_element,num_element)   :: temp_matrix_mpi
        !----------------------------------------------------------------------------------------------------------------------------------
        call MPI_ALLREDUCE(temp_matrix, temp_matrix_mpi, num_element*num_element, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, ierr)
        temp_matrix = temp_matrix_mpi
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE MBD_at_rsSCS(ene_mbd_rsSCS)
        !----------------------------------------------------------------------------------------------------------------------------------
        ! local variable 
        !----------------------------------------------------------------------------------------------------------------------------------
        INTEGER                :: i_freq
        INTEGER                :: i_myatom
        INTEGER                :: errorflag
        INTEGER                :: LWORK
        REAL*8                 :: C6_free
        REAL*8                 :: alpha_free
        REAL*8                 :: R_free
        REAL*8                 :: mol_c6_coeff
        REAL*8                 :: C6_atom
        REAL*8                 :: R_vdw_free
        REAL*8                 :: ene_mbd_rsSCS
        REAL*8, DIMENSION(9)   :: WORK(9)
        REAL*8, DIMENSION(3)   :: mol_pol_eigen(3)
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Input
        !----------------------------------------------------------------------------------------------------------------------------------
        ! o hirshfeld_volume from hirshfeld partioning 
        ! o relay_matrix is 3N x 3N (imaginary) frequency dependent matrix
        ! o R_p is effective dipole spread of quantum harmonic oscilator
        ! o alpha_omega is single pole polarizability of atom in molecule
        ! o coupled_atom_pol contains atom resolved screened polarizabilities at every imaginary frequency 
        ! o alpha_eff is atom resolved screened polarizabilty at ZERO frequency 
        ! o C6_eff is atom resolved screened C6 coefficients calculated by integrating coupled_atom_pol over entire frequency range 
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Allocation
        !----------------------------------------------------------------------------------------------------------------------------------
        ! JJ: as distributed computing, you can not allocate the same matrix everywhere 3N*3N
        !----------------------------------------------------------------------------------------------------------------------------------
        IF(.NOT. ALLOCATED(relay_matrix))              ALLOCATE(relay_matrix(3*n_atoms,3*n_atoms))
        IF(.NOT. ALLOCATED(R_p))                       ALLOCATE(R_p(n_atoms))
        IF(.NOT. ALLOCATED(Rvdw_iso))                  ALLOCATE(Rvdw_iso(n_atoms))
        IF(.NOT. ALLOCATED(alpha_omega))               ALLOCATE(alpha_omega(n_atoms))
        IF(.NOT. ALLOCATED(coupled_atom_pol))          ALLOCATE(coupled_atom_pol(20,n_atoms))
        IF(.NOT. ALLOCATED(alpha_eff))                 ALLOCATE(alpha_eff(n_atoms))
        IF(.NOT. ALLOCATED(C6_eff))                    ALLOCATE(C6_eff(n_atoms))
        IF(n_periodic .GT. 0) THEN
            IF(.NOT. ALLOCATED(relay_matrix_periodic)) ALLOCATE(relay_matrix_periodic(3*n_atoms,3*n_atoms))
        END IF
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Initialization
        !----------------------------------------------------------------------------------------------------------------------------------
        C6_eff               = 0.0d0 
        alpha_eff            = 0.0d0
        mol_c6_coeff         = 0.0d0
        coupled_atom_pol     = 0.0d0
        casimir_omega        = 0.0d0
        casimir_omega_weight = 0.0d0 
        CALL gauss_legendre_grid()
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Write Title
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)')"| Many-Body Dispersion (MBD@rsSCS) energy "
        CALL output_string(info_str)
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)')"| System dynamic molecular polarizability alpha(iw) (bohr^3)"
        CALL output_string(info_str)
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)')"| ---------------------------------------------------------------------------"
        CALL output_string(info_str)
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)')"|  omega(Ha)   alpha_iso      alpha_xx       alpha_yy       alpha_zz"
        CALL output_string(info_str)
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Loop over Casimir-Polder frequencies
        !----------------------------------------------------------------------------------------------------------------------------------
        DO i_freq=0,20
            !------------------------------------------------------------------------------------------------------------------------------
            ! zeroout array before getting actual frequency dependent parameters
            !------------------------------------------------------------------------------------------------------------------------------
            R_p            = 0.0d0
            Rvdw_iso       = 0.0d0 
            alpha_omega    = 0.0d0
            relay_matrix   = 0.0d0
            mol_pol_tensor = 0.0d0
            !------------------------------------------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------------------------------------------
            ! compute frequency dependent proprieties (R_p)
            !------------------------------------------------------------------------------------------------------------------------------
            DO i_myatom=1,n_atoms
                CALL get_vdw_param(atom_name(i_myatom),C6_free,alpha_free,R_vdw_free)
                CALL get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),C6_free,alpha_free,casimir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
                Rvdw_iso(i_myatom)= (hirshfeld_volume(i_myatom)**0.333333333333333333333333333d0)*R_vdw_free
            END DO
            !------------------------------------------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------------------------------------------
            ! compute screened polarizability using isotropically damped tensor, (OUT: relay_matrix^{-1})
            !------------------------------------------------------------------------------------------------------------------------------
            CALL calculate_scs_matrix() 
            !------------------------------------------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------------------------------------------
            ! Contract relay_matrix^-1 for polarizability (3*3 block)
            !------------------------------------------------------------------------------------------------------------------------------
            CALL contract_matrix(mol_pol_tensor)                                                                            ! equation (18)
            !------------------------------------------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------------------------------------------
            ! get coupled atomic "local" polarizability by summing over row or columns of relay tensor
            !------------------------------------------------------------------------------------------------------------------------------
            IF(i_freq .EQ. 0) THEN
               CALL calculate_screened_polarizability(alpha_eff)                                                    ! static polarizability
            ELSE
               CALL calculate_screened_polarizability(coupled_atom_pol(i_freq,:))                      ! frequency dependent polarizability
            END IF 
            !------------------------------------------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------------------------------------------
            ! Casimir-Polder integral for molecular C6-coefficent 
            !------------------------------------------------------------------------------------------------------------------------------
            LWORK=9 
            CALL DSYEV('V','U',3,mol_pol_tensor,3,mol_pol_eigen,WORK,LWORK,errorflag)                ! JJ: this one is not symmetric either
            IF(i_freq .GE. 1) THEN
                mol_c6_coeff = mol_c6_coeff + (casimir_omega_weight(i_freq)* (sum(mol_pol_eigen)/3.0)**2)
            END IF
            
            WRITE(info_str,'(2x,"| ",F10.6,4(e15.6))') casimir_omega(i_freq), SUM(mol_pol_eigen)/3.0, mol_pol_eigen(1), mol_pol_eigen(2), mol_pol_eigen(3)
            CALL output_string(info_str)  
            !------------------------------------------------------------------------------------------------------------------------------
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)') "| ---------------------------------------------------------------------------"
        CALL output_string(info_str)  
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A,f20.8,2x,A)') "| System C6 coefficient    :  ",0.954929658551372d0*mol_c6_coeff,"  hartree bohr^6"
        CALL output_string(info_str)  
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)') "| ---------------------------------------------------------------------------"
        CALL output_string(info_str)  
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)') "| Partitioned C6 (hartree bohr^6) | alpha (bohr^3)  of atom in system"
        CALL output_string(info_str)  
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)') "| ---------------------------------------------------------------------------"
        CALL output_string(info_str)  
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        DO i_myatom=1,n_atoms,1
            !------------------------------------------------------------------------------------------------------------------------------
            C6_atom = 0.0d0 
            !------------------------------------------------------------------------------------------------------------------------------
            DO i_freq=1,20,1
                C6_atom = C6_atom + (casimir_omega_weight(i_freq)*coupled_atom_pol(i_freq,i_myatom)**2)
            END DO
            !------------------------------------------------------------------------------------------------------------------------------
            WRITE(info_str,'(2x,A,I4,2x,A,f12.6,6x,f12.6)') "| ATOM ",i_myatom,atom_name(i_myatom), C6_atom*0.954929658551372d0,alpha_eff(i_myatom)
            CALL output_string(info_str)  
            !------------------------------------------------------------------------------------------------------------------------------
            C6_eff(i_myatom)=C6_atom*0.954929658551372d0                                                        ! effctive C6 for each atom
            !------------------------------------------------------------------------------------------------------------------------------
        END DO 
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        WRITE(info_str,'(2x,A)') "| ---------------------------------------------------------------------------"
        CALL output_string(info_str)
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Compute many-body-dispersion energy new alpha_eff and C6_eff
        !----------------------------------------------------------------------------------------------------------------------------------
        call get_mbd_rSCS(ene_mbd_rsSCS)
        !----------------------------------------------------------------------------------------------------------------------------------

        !----------------------------------------------------------------------------------------------------------------------------------
        write(info_str,'(2X,A,F20.9,A,F20.9,A)') "| MBD@rsSCS energy   :", ene_mbd_rsSCS," Ha     ",ene_mbd_rsSCS*27.211," eV"
        call output_string(info_str)
        !----------------------------------------------------------------------------------------------------------------------------------


        !----------------------------------------------------------------------------------------------------------------------------------
        ! clean up
        !----------------------------------------------------------------------------------------------------------------------------------
        IF(ALLOCATED(R_p))                      DEALLOCATE(R_p)
        IF(ALLOCATED(Rvdw_iso))                 DEALLOCATE(Rvdw_iso)
        IF(ALLOCATED(alpha_omega))              DEALLOCATE(alpha_omega)
        IF(ALLOCATED(relay_matrix))             DEALLOCATE(relay_matrix)
        IF(ALLOCATED(relay_matrix_periodic))    DEALLOCATE(relay_matrix_periodic)
        IF(ALLOCATED(coupled_atom_pol))         DEALLOCATE(coupled_atom_pol)
        IF(ALLOCATED(alpha_eff))                DEALLOCATE(alpha_eff)
        IF(ALLOCATED(C6_eff))                   DEALLOCATE(C6_eff)
        !----------------------------------------------------------------------------------------------------------------------------------
        RETURN
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE MBD_at_rsSCS
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE calculate_scs_matrix()
    !--------------------------------------------------------------------------------------------------------------------------------------
        INTEGER                       :: i_row
        INTEGER                       :: i_col
        INTEGER                       :: i_index
        INTEGER                       :: j_index
        INTEGER                       :: i_lattice
        INTEGER                       :: j_lattice
        INTEGER                       :: k_lattice
        INTEGER                       :: periodic_cell_i
        INTEGER                       :: periodic_cell_j
        INTEGER                       :: periodic_cell_k
        INTEGER                       :: errorflag
        INTEGER, DIMENSION(3*n_atoms) :: IPIV                  ! LAPACK related 
        REAL*8,  DIMENSION(3*n_atoms) :: WORK                  ! LAPACK related
        REAL*8,  DIMENSION(3,3)       :: TPP
        REAL*8,  DIMENSION(3)         :: dxyz
        REAL*8,  DIMENSION(3)         :: coord_curr
        REAL*8                        :: r_ij
        REAL*8                        :: r_pp
        REAL*8                        :: Rvdw12
        REAL*8                        :: beta
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! initio values
        !----------------------------------------------------------------------------------------------------------------------------------
        relay_matrix=0.0d0
        !----------------------------------------------------------------------------------------------------------------------------------
        SELECT CASE (flag_xc)
            CASE (1) ! PBE
                beta=0.83
            CASE (2) !PBE0
                beta=0.85
            CASE (3) !HSE 
                beta=0.85
            CASE DEFAULT
                beta=1.0
                WRITE(info_str,'(A)')"***  WARNING range seperated parameter beta is not defined for"
                CALL output_string(info_str)
                WRITE(info_str,'(A)')"this fxc defaulting it to 0.0 , which will give MBD energy"
                CALL output_string(info_str)
                WRITE(info_str,'(A)')"WITHOUT short range damping"
                CALL output_string(info_str)
        END SELECT
        !----------------------------------------------------------------------------------------------------------------------------------
        
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! compute relay matrix of cluster or unit cell
        !----------------------------------------------------------------------------------------------------------------------------------
        DO i_row=1,n_atoms                                 ! loop through block index i
            IF(iproc .EQ. task_list(i_row)) THEN      
                DO i_col=i_row,n_atoms                     ! loop through block index j
                    TPP=0.0d0
                    IF(i_row .EQ. i_col) THEN              ! diagonal block
                        DO i_index=1,3
                            DO j_index=1,3
                                IF(i_index .EQ. j_index) THEN
                                    relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=1.d0/alpha_omega(i_row)
                                ELSE 
                                    relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=0.d0
                                END IF   
                            END DO
                        END DO
                    ELSE
                        dxyz(:) = coords(:,i_col)-coords(:,i_row)
                        r_ij    = DSQRT((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0)
                        r_pp    = DSQRT(R_p(i_row)**2 + R_p(i_col)**2)                         
                        Rvdw12  = Rvdw_iso(i_row) +  Rvdw_iso(i_col)                           ! ??
                        CALL SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)
                        DO i_index=1,3,1
                           DO j_index=1,3,1
                            relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=TPP(i_index,j_index)
                            relay_matrix(3*i_col-3+j_index,3*i_row-3+i_index)=TPP(i_index,j_index)
                           END DO
                        END DO
                    END IF
                END DO
            END IF
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
        
        CALL sync_tensors(relay_matrix,3*n_atoms)
        
        !----------------------------------------------------------------------------------------------------------------------------------
        IF (n_periodic .GT. 0) THEN
            relay_matrix_periodic=0.0d0 
            !------------------------------------------------------------------------------------------------------------------------------
            IF(.NOT. mbd_scs_vacuum_axis(1)) THEN
                periodic_cell_i = CEILING((mbd_scs_dip_cutoff/bohr) / SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 + lattice_vector(3,1)**2))
            ELSE
                periodic_cell_i = 0
            END IF 
            !------------------------------------------------------------------------------------------------------------------------------
            IF(.NOT. mbd_scs_vacuum_axis(2)) THEN
                periodic_cell_j = CEILING((mbd_scs_dip_cutoff/bohr) / SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 + lattice_vector(3,2)**2))
            ELSE
                periodic_cell_j = 0
            END IF 
            !------------------------------------------------------------------------------------------------------------------------------
            IF(.NOT. mbd_scs_vacuum_axis(3)) THEN
                periodic_cell_k = CEILING((mbd_scs_dip_cutoff/bohr) / SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 + lattice_vector(3,3)**2))
            ELSE
                periodic_cell_k = 0
            END IF 
            !------------------------------------------------------------------------------------------------------------------------------
            ! sum over images for long range interactions
            !------------------------------------------------------------------------------------------------------------------------------
            DO i_lattice = -periodic_cell_i, periodic_cell_i,1
                DO j_lattice = -periodic_cell_j, periodic_cell_j,1
                    DO k_lattice = -periodic_cell_k, periodic_cell_k,1
                        IF((ABS(i_lattice) .NE. 0) .OR. (ABS(j_lattice) .NE. 0) .OR. (ABS(k_lattice) .NE. 0)) THEN !#1
                            DO i_row = 1, n_atoms
                                IF(iproc .EQ. task_list(i_row)) THEN
                                    DO i_col = i_row, n_atoms
                                        !-------------------------------------------------------------------------------------------------
                                        coord_curr = 0.d0
                                        dxyz = 0.d0
                                        r_ij = 0.d0
                                        TPP  = 0.d0
                                        Rvdw12= 0.d0  
                                        !-------------------------------------------------------------------------------------------------
                                        ! find the coordinate of images
                                        !-------------------------------------------------------------------------------------------------
                                        coord_curr(:) = coords(:,i_col) + i_lattice*lattice_vector(:,1) + &
                                                                          j_lattice*lattice_vector(:,2) + &
                                                                          k_lattice*lattice_vector(:,3)
                                        !-------------------------------------------------------------------------------------------------
                                        dxyz(:) = coords(:,i_row) - coord_curr(:)
                                        r_ij    = SQRT((dxyz(1))**2.0d0 + (dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
                                        r_pp    = SQRT(R_p(i_row)**2 + R_p(i_col)**2)
                                        Rvdw12  = Rvdw_iso(i_row) + Rvdw_iso(i_col)
                                        !-------------------------------------------------------------------------------------------------
                                        IF(r_ij .LE. mbd_scs_dip_cutoff/bohr) THEN
                                            CALL SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)
                                            DO i_index=1,3
                                                DO j_index=1,3
                                                    relay_matrix_periodic(3*i_row-3+i_index,3*i_col-3+j_index) = &
                                                    relay_matrix_periodic(3*i_row-3+i_index,3*i_col-3+j_index) + TPP(i_index,j_index)
                                                    IF(i_col .NE. i_row) THEN
                                                        relay_matrix_periodic(3*i_col-3+j_index,3*i_row-3+i_index) = &
                                                        relay_matrix_periodic(3*i_col-3+j_index,3*i_row-3+i_index) + TPP(i_index,j_index)
                                                    END IF
                                                END DO
                                            END DO
                                        END IF
                                        !-------------------------------------------------------------------------------------------------
                                    END DO
                                END IF
                            END DO
                        END IF  !#1
                    END DO
                END DO
            END DO
            !------------------------------------------------------------------------------------------------------------------------------
            CALL sync_tensors(relay_matrix_periodic,3*n_atoms)
            relay_matrix = relay_matrix + relay_matrix_periodic  
            !------------------------------------------------------------------------------------------------------------------------------
        END IF
        
        !----------------------------------------------------------------------------------------------------------------------------------
        CALL DGETRF(3*n_atoms, 3*n_atoms, relay_matrix, 3*n_atoms, IPIV, errorflag)            ! compute LU factorization of A
        !----------------------------------------------------------------------------------------------------------------------------------
        IF(errorflag .NE. 0) THEN
            WRITE(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
            CALL output_string(info_str)
        END IF
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        CALL DGETRI(3*n_atoms, relay_matrix, 3*n_atoms, IPIV, WORK, 3*n_atoms, errorflag )     ! compute inverse of A (solve LU A^{-1} = I)
        !----------------------------------------------------------------------------------------------------------------------------------
        IF(errorflag .NE. 0) THEN
            WRITE(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
            CALL output_string(info_str)
        END IF
        RETURN
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE calculate_scs_matrix
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE SCS_TENSOR_MBD_rsSCS(dxyz, r_ij, r_pp, Rvdw12, beta, TPP)                                                    ! equation (14)
    !--------------------------------------------------------------------------------------------------------------------------------------
        REAL*8, INTENT(IN)                   :: r_ij
        REAL*8, INTENT(IN)                   :: r_pp
        REAL*8, INTENT(IN)                   :: Rvdw12
        REAL*8, INTENT(IN)                   :: beta
        REAL*8, INTENT(IN),  DIMENSION(3)    :: dxyz
        REAL*8, INTENT(OUT), DIMENSION(3,3)  :: TPP           ! dipole-dipole interaction tensor between two QHO Tp-p damped with Fermi damping 
        REAL*8                               :: derf, erf
        !----------------------------------------------------------------------------------------------------------------------------------
        ! local vars
        !----------------------------------------------------------------------------------------------------------------------------------
        REAL*8,              DIMENSION(3,3)  :: TPQ
        REAL*8,              DIMENSION(3,3)  :: r_tensor
        REAL*8                               :: zeta_l
        REAL*8                               :: zeta_r
        REAL*8                               :: ratio
        REAL*8                               :: d_param       ! Fermi damping parameter 
        REAL*8                               :: fermi_param
        integer                              :: i
        integer                              :: j
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        TPP(:,:)    = 0.0d0
        d_param     = 6.0d0
        fermi_param = r_ij/(beta*Rvdw12)
        ratio       = r_ij/r_pp
        zeta_l      = (derf(ratio) - (dsqrt(4.0d0/pi) * ratio * dexp(-(ratio**2.0d0)))) / (r_ij**5.0d0)
        zeta_r      = (dsqrt(16.0d0/pi) * dexp(-(ratio**2.0d0))) / (r_pp * r_pp * r_pp * r_ij * r_ij)
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! R_{i}R_{j} tensor
        !----------------------------------------------------------------------------------------------------------------------------------
        DO i=1,3
            DO j=1,3
                r_tensor(i,j)=dxyz(i)*dxyz(j)
            END DO
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        TPQ = r_tensor*3.0d0
        DO i=1,3
           TPQ(i,i)=TPQ(i,i)-(r_ij*r_ij)
        end do
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        TPP = -TPQ * zeta_l + r_tensor * zeta_r
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Fermi type damping  
        !----------------------------------------------------------------------------------------------------------------------------------
        TPP = TPP*(1.0-(1.0/(1.0+exp(-d_param*(fermi_param-1.0)))))
        !----------------------------------------------------------------------------------------------------------------------------------
        RETURN
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE SCS_TENSOR_MBD_rsSCS
    !--------------------------------------------------------------------------------------------------------------------------------------

subroutine MBD_TENSOR_MBD_rsSCS(dxyz,r_ij,Rvdw12,beta,Tij)
      ! o Tij- Fermi type*Grad_i X Grad_j (1/r)
      ! tensor derived from bare coulomb potential damped with
      ! Fermi type damping

      real*8,dimension(3),intent(in) :: dxyz
      real*8,intent(in) :: r_ij
      real*8,intent(in) :: Rvdw12
      real*8,intent(in) :: beta
      real*8,dimension(3,3),intent(out) :: Tij
      ! This needs declarding for the PGI Architecture
      real*8 :: erf
      real*8 :: d_param,fermi_param
      ! local vars
      integer :: i_index, j_index

      d_param = 6.0
      fermi_param =r_ij/(beta*Rvdw12)

             do i_index=1,3,1
               do j_index=1,3,1
                Tij(i_index,j_index)=3.d0*dxyz(i_index)*dxyz(j_index)
               enddo
             enddo

             do i_index=1,3,1
               Tij(i_index,i_index) = Tij(i_index,i_index) - r_ij**2
             enddo
             Tij=Tij/r_ij**5
             Tij=-1.d0*Tij

             Tij=Tij*(1.0/(1.0+exp(-d_param*(fermi_param-1.0))))
      return
endsubroutine MBD_TENSOR_MBD_rsSCS


    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE get_mbd_rSCS(ene_mbd_rsSCS)
    !--------------------------------------------------------------------------------------------------------------------------------------
        REAL*8,  DIMENSION(:,:), ALLOCATABLE   :: cfdm_hamiltonian 
        REAL*8,  DIMENSION(:,:), ALLOCATABLE   :: cfdm_hamiltonian_periodic
        REAL*8,  DIMENSION(:,:), ALLOCATABLE   :: coords_SL
        REAL*8,  DIMENSION(:),   ALLOCATABLE   :: cfdm_eigenvalues
        REAL*8,  DIMENSION(:),   ALLOCATABLE   :: omega_cfdm_SL
        REAL*8,  DIMENSION(:),   ALLOCATABLE   :: R_vdw_SL
        REAL*8,  DIMENSION(:),   ALLOCATABLE   :: alpha_eff_SL
        REAL*8,  DIMENSION(:),   ALLOCATABLE   :: WORK
        REAL*8,  DIMENSION(:),   ALLOCATABLE   :: task_list_SL  
        REAL*8,  DIMENSION(3,3)                :: TPP
        REAL*8,  DIMENSION(3,3)                :: lattice_vector_SL
        REAL*8,  DIMENSION(3)                  :: dxyz,coord_curr
        REAL*8                                 :: r_ij
        REAL*8                                 :: C6_free
        REAL*8                                 :: alpha_free
        REAL*8                                 :: R_vdw_free
        REAL*8                                 :: Rvdw_12
        REAL*8                                 :: beta
        REAL*8                                 :: CFDM_prefactor
        REAL*8                                 :: E_int
        REAL*8                                 :: E_nonint
        REAL*8                                 :: sigma
        REAL*8                                 :: ene_mbd_rsSCS
        INTEGER                                :: errorflag
        INTEGER                                :: i_atom
        INTEGER                                :: j_atom
        INTEGER                                :: SL_i                      ! MBD super cell index's
        INTEGER                                :: SL_j                      ! MBD super cell index's
        INTEGER                                :: SL_k                      ! MBD super cell index's
        INTEGER                                :: i_index                   ! 3x3 block index
        INTEGER                                :: j_index                   ! 3x3 block index
        INTEGER                                :: i_lattice                 ! lattice increament index
        INTEGER                                :: j_lattice                 ! lattice increament index
        INTEGER                                :: k_lattice                 ! lattice increament index
        INTEGER                                :: periodic_cell_i
        INTEGER                                :: periodic_cell_j
        INTEGER                                :: periodic_cell_k
        INTEGER                                :: NIMAG
        INTEGER                                :: n_atoms_SL
        INTEGER                                :: LWORK
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! set beta
        !----------------------------------------------------------------------------------------------------------------------------------
        SELECT CASE (flag_xc)
            CASE (1) ! PBE
                beta=0.83
            CASE (2) !PBE0
                beta=0.85
            CASE (3) !HSE 
                beta=0.85
            CASE DEFAULT
                beta=1.0
                WRITE(info_str,'(A)') "***  WARNING range seperated parameter beta is not defined for"
                CALL output_string(info_str)
                WRITE(info_str,'(A)') "this fxc defaulting it to 0.0 , which will give MBD energy"
                CALL output_string(info_str)
                WRITE(info_str,'(A)') "WITHOUT short range damping"
                CALL output_string(info_str)
        END SELECT
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        IF ( n_periodic .EQ. 0) THEN
            !------------------------------------------------------------------------------------------------------------------------------
            ! For Cluster calculation NO super cell needed intiate this index so as to avoide any invalid 
            ! floting point when normalizing energy in case of cluster/molecule
            !------------------------------------------------------------------------------------------------------------------------------
            SL_i = 1 
            SL_j = 1
            SL_k = 1
            !------------------------------------------------------------------------------------------------------------------------------
            ! number of atoms in cluster /molecule
            !------------------------------------------------------------------------------------------------------------------------------
            n_atoms_SL = n_atoms
            !------------------------------------------------------------------------------------------------------------------------------
            IF(.NOT. ALLOCATED(cfdm_hamiltonian))          ALLOCATE(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
            IF(.NOT. ALLOCATED(cfdm_eigenvalues))          ALLOCATE(cfdm_eigenvalues(3*n_atoms_SL)) 
            IF(.NOT. ALLOCATED(coords_SL))                 ALLOCATE(coords_SL(3,n_atoms_SL)) 
            IF(.NOT. ALLOCATED(omega_cfdm_SL))             ALLOCATE(omega_cfdm_SL(n_atoms_SL))
            IF(.NOT. ALLOCATED(R_vdw_SL))                  ALLOCATE(R_vdw_SL(n_atoms_SL))
            IF(.NOT. ALLOCATED(alpha_eff_SL))              ALLOCATE(alpha_eff_SL(n_atoms_SL))
            IF(.NOT. ALLOCATED(WORK))                      ALLOCATE(WORK((3*n_atoms_SL)*(3+(3*n_atoms_SL)/2)))
            IF(.NOT. ALLOCATED(task_list_SL))              ALLOCATE(task_list_SL(n_atoms_SL))
            !------------------------------------------------------------------------------------------------------------------------------
            ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
            !------------------------------------------------------------------------------------------------------------------------------
            cfdm_hamiltonian=0.0d0 
            cfdm_eigenvalues=0.0d0
            !------------------------------------------------------------------------------------------------------------------------------
            ! RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN CLUSTER/MOLECULE
            !------------------------------------------------------------------------------------------------------------------------------
            DO i_atom=1,n_atoms_SL
                task_list_SL(i_atom) = MOD(i_atom,nproc)  
            END DO
            !------------------------------------------------------------------------------------------------------------------------------
            DO i_atom =1,n_atoms
                R_vdw_free=0.0d0
                alpha_free=0.0d0
                C6_free=0.d0
                CALL get_vdw_param(atom_name(i_atom),C6_free,alpha_free,R_vdw_free)
                !--------------------------------------------------------------------------------------------------------------------------
                R_vdw_SL(i_atom)= R_vdw_free*((alpha_eff(i_atom)/alpha_free)**0.333333333333333333333333333d0)
                omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                coords_SL(:,i_atom)   = coords(:,i_atom)
                alpha_eff_SL(i_atom) = alpha_eff(i_atom)
            END DO
            !------------------------------------------------------------------------------------------------------------------------------
        ELSE
            !------------------------------------------------------------------------------------------------------------------------------
            IF(.NOT. mbd_scs_vacuum_axis(1)) THEN
                SL_i = CEILING((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 +lattice_vector(3,1)**2))
            ELSE
                SL_i = 1 
            END IF
            !------------------------------------------------------------------------------------------------------------------------------
            IF(.NOT. mbd_scs_vacuum_axis(2)) THEN
                SL_j = CEILING((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 +lattice_vector(3,2)**2)) 
            ELSE
                SL_j = 1 
            END IF
            !------------------------------------------------------------------------------------------------------------------------------
            IF(.NOT. mbd_scs_vacuum_axis(3)) THEN
                SL_k = CEILING((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 +lattice_vector(3,3)**2)) 
            ELSE
                SL_k = 1 
            END IF
            !------------------------------------------------------------------------------------------------------------------------------
          
            !------------------------------------------------------------------------------------------------------------------------------
            ! number of atoms in super cell 
            !------------------------------------------------------------------------------------------------------------------------------
            n_atoms_SL = n_atoms*SL_i*SL_j*SL_k 
            !------------------------------------------------------------------------------------------------------------------------------
            WRITE(info_str,'(2X,A,i3,A,i3,A,i3,A)') "| Creating super cell of dimension",  SL_i," X ", SL_j," X ", SL_k, " in MBD calculation" 
            CALL output_string(info_str)
            WRITE(info_str,'(2x,"| containing", I6,A)') n_atoms_SL, "  atom"
            CALL output_string(info_str)
            !------------------------------------------------------------------------------------------------------------------------------
            IF(.NOT. ALLOCATED(cfdm_hamiltonian))             ALLOCATE(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
            IF(.NOT. ALLOCATED(cfdm_hamiltonian_periodic))    ALLOCATE(cfdm_hamiltonian_periodic(3*n_atoms_SL,3*n_atoms_SL)) 
            IF(.NOT. ALLOCATED(cfdm_eigenvalues))             ALLOCATE(cfdm_eigenvalues(3*n_atoms_SL)) 
            IF(.NOT. ALLOCATED(coords_SL))                    ALLOCATE(coords_SL(3,n_atoms_SL)) 
            IF(.NOT. ALLOCATED(omega_cfdm_SL))                ALLOCATE(omega_cfdm_SL(n_atoms_SL))
            IF(.NOT. ALLOCATED(R_vdw_SL))                     ALLOCATE(R_vdw_SL(n_atoms_SL))
            IF(.NOT. ALLOCATED(alpha_eff_SL))                 ALLOCATE(alpha_eff_SL(n_atoms_SL))
            IF(.NOT. ALLOCATED(WORK))                         ALLOCATE(WORK(3*n_atoms_SL*(3+(3*n_atoms_SL)/2)))
            IF(.NOT. ALLOCATED(task_list_SL))                 ALLOCATE(task_list_SL(n_atoms_SL))
            !------------------------------------------------------------------------------------------------------------------------------
            ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
            !------------------------------------------------------------------------------------------------------------------------------
            cfdm_hamiltonian = 0.0d0 
            cfdm_eigenvalues = 0.0d0
            cfdm_hamiltonian_periodic = 0.d0 
            !------------------------------------------------------------------------------------------------------------------------------
            ! RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN SUPERCELL/MOLECULE
            !------------------------------------------------------------------------------------------------------------------------------
            DO i_atom=1,n_atoms_SL 
                task_list_SL(i_atom) = MOD(i_atom,nproc) 
            END DO
            !------------------------------------------------------------------------------------------------------------------------------
            ! Get all the parameter required for MBD with super cell calculation
            !------------------------------------------------------------------------------------------------------------------------------
            j_atom = 0
            DO i_lattice = 0, SL_i-1
                DO j_lattice = 0, SL_j-1
                    DO k_lattice = 0, SL_k-1
                        DO i_atom = 1, n_atoms             ! <--- atom index
                            j_atom = j_atom +1             ! <--- dummy index supercell
                            coords_SL(:,j_atom) = coords(:,i_atom) + i_lattice*lattice_vector(:,1) + j_lattice*lattice_vector(:,2) + k_lattice*lattice_vector(:,3)   
                            R_vdw_free=0.0d0
                            alpha_free=0.0d0  
                            CALL get_vdw_param(atom_name(i_atom),C6_free,alpha_free,R_vdw_free) !FIXME
                            
                            R_vdw_SL(j_atom)=R_vdw_free*((alpha_eff(i_atom)/alpha_free)**0.333333333333333333333333333d0)
                            alpha_eff_SL(j_atom) = alpha_eff(i_atom)
                            omega_cfdm_SL(j_atom) =(4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                        END DO 
                    END DO            
                END DO            
            END DO            
            lattice_vector_SL(:,1)  = lattice_vector(:,1)*SL_i 
            lattice_vector_SL(:,2)  = lattice_vector(:,2)*SL_j 
            lattice_vector_SL(:,3)  = lattice_vector(:,3)*SL_k 
            !------------------------------------------------------------------------------------------------------------------------------
        END IF 
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Construct CFDM hamiltonian matrix 
        !----------------------------------------------------------------------------------------------------------------------------------
        DO i_atom=1,n_atoms_SL
            IF(iproc .EQ. task_list_SL(i_atom)) THEN
                DO j_atom=i_atom,n_atoms_SL
                    IF(i_atom .EQ. j_atom) THEN
                        DO i_index=1,3
                            DO j_index=1,3
                                IF(i_index .EQ. j_index) THEN
                                    cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = omega_cfdm_SL(i_atom)**2.0
                                END IF  
                            END DO
                        END DO
                    ELSE
                        r_ij = 0.0d0
                        TPP  = 0.0d0
                        dxyz = 0.d0
                        dxyz(:) = coords_SL(:,i_atom)-coords_SL(:,j_atom)
                        r_ij=DSQRT((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
                        Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
                        sigma=(r_ij/Rvdw_12)**beta
                        CALL MBD_TENSOR_MBD_rsSCS(dxyz,r_ij,Rvdw_12,beta,TPP)
                        CFDM_prefactor=omega_cfdm_SL(i_atom)*omega_cfdm_SL(j_atom)*DSQRT(alpha_eff_SL(i_atom)*alpha_eff_SL(j_atom))
                        
                        ! Transfer each dipole matrix to CFDM hamiltonian i.j accordingly         
                        DO i_index=1,3
                            DO j_index=1,3
                                cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = TPP(i_index,j_index)*CFDM_prefactor
                                cfdm_hamiltonian(3*j_atom-3+j_index,3*i_atom-3+i_index) = TPP(i_index,j_index)*CFDM_prefactor
                            END DO
                        END DO
                    END IF
                END DO
            END IF
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
        
        CALL sync_tensors(cfdm_hamiltonian,3*n_atoms_SL)
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Adds dipole field due to image cells based spherical cutoff mbd_cfdm_dip_cutoff
        !----------------------------------------------------------------------------------------------------------------------------------
        IF (n_periodic .GT. 0) THEN 
            IF(.NOT. mbd_scs_vacuum_axis(1)) THEN
                periodic_cell_i = CEILING((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,1)**2 + lattice_vector_SL(2,1)**2 + lattice_vector_SL(3,1)**2))
            ELSE
                periodic_cell_i = 0 
            END IF   
            
            IF (.NOT. mbd_scs_vacuum_axis(2)) THEN
                periodic_cell_j = CEILING((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,2)**2 + lattice_vector_SL(2,2)**2 + lattice_vector_SL(3,2)**2)) 
            ELSE
                periodic_cell_j = 0 
            END IF   
            
            IF (.NOT.mbd_scs_vacuum_axis(3)) THEN
                periodic_cell_k = CEILING((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,3)**2 + lattice_vector_SL(2,3)**2 + lattice_vector_SL(3,3)**2)) 
            ELSE
                periodic_cell_k = 0 
            END IF
         
            DO i_lattice = -periodic_cell_i, periodic_cell_i
                DO j_lattice = -periodic_cell_j, periodic_cell_j
                    DO k_lattice = -periodic_cell_k, periodic_cell_k
                        IF((ABS(i_lattice) .NE. 0) .OR. (ABS(j_lattice) .NE. 0) .OR. (ABS(k_lattice) .NE. 0)) THEN
                            DO i_atom=1,n_atoms_SL
                                IF(iproc .EQ. task_list_SL(i_atom)) THEN
                                    ! LOOP GOES OVER UPPPER TRIANGLE OF HAMILTONIAN 
                                    DO j_atom=i_atom,n_atoms_SL
                                        r_ij=0.0d0
                                        TPP=0.0d0
                                        dxyz=0.d0
                                        coord_curr(:) = coords_SL(:,i_atom) + i_lattice*lattice_vector_SL(:,1) + &
                                                                              j_lattice*lattice_vector_SL(:,2) + &
                                                                              k_lattice*lattice_vector_SL(:,3)
                                        dxyz(:) = coords_SL(:,j_atom) - coord_curr(:)
                                        r_ij = DSQRT((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
                                        
                                        IF(r_ij .LE. mbd_cfdm_dip_cutoff/bohr) THEN
                                            Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
                                            sigma=(r_ij/Rvdw_12)**beta
                                            CALL MBD_TENSOR_MBD_rsSCS(dxyz,r_ij,Rvdw_12,beta,TPP) 
                                            CFDM_prefactor=omega_cfdm_SL(i_atom)*omega_cfdm_SL(j_atom) * SQRT(alpha_eff_SL(i_atom)*alpha_eff_SL(j_atom))
                                            DO i_index=1,3
                                                DO j_index=1,3
                                                    ! FILL UPPER BLOCK
                                                    cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index) = &
                                                    cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index) + (TPP(i_index,j_index)*CFDM_prefactor)
                                                    ! FILL LOWER BLOCK ALSO MAKE SURE THAT YOU DONT ADD FIELD
                                                    ! CONTRIBUTION TWO TIMES IN DIAGONAL BLOCK below is the check
                                                    IF (i_atom .NE. j_atom) THEN
                                                        cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index) = &
                                                        cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index) + (TPP(i_index,j_index)*CFDM_prefactor)
                                                    END IF
                                                END DO
                                            END DO
                                        END IF 
                                    END DO
                                END IF
                            END DO
                        END IF
                    END DO
                END DO
            END DO
            CALL sync_tensors(cfdm_hamiltonian_periodic,3*n_atoms_SL)    
            cfdm_hamiltonian = cfdm_hamiltonian + cfdm_hamiltonian_periodic
        END IF
        
        errorflag=0
        LWORK=3*n_atoms_SL*(3+(3*n_atoms_SL)/2)
        CALL DSYEV('V','U',3*n_atoms_SL,cfdm_hamiltonian,3*n_atoms_SL,cfdm_eigenvalues,WORK,LWORK,errorflag)
        
        E_int=0.0d0
        E_nonint =0.0d0
        ene_mbd_rsSCS = 0.0
        NIMAG=0 
        IF(errorflag .EQ. 0) THEN
            DO i_atom = 1,n_atoms_SL
                E_nonint = E_nonint + omega_cfdm_SL(i_atom)    
            END DO
            
            DO i_atom =1,3*n_atoms_SL
                IF(cfdm_eigenvalues(i_atom) .GE. 0.0d0) THEN
                    E_int = E_int + DSQRT(cfdm_eigenvalues(i_atom))
                ELSE
                    NIMAG= NIMAG +1  
                END IF 
            END DO
            
            ene_mbd_rsSCS = ((0.5*E_int)-(1.5* E_nonint))/(SL_i*SL_j*SL_k)
            IF(NIMAG .GT. 0) THEN
                WRITE(info_str,'(A,I4,A)')"***WARNING: found ",NIMAG," negative eigenvalues in MBD energy calculation."
                CALL output_string(info_str)  
            END IF
        ELSE
            WRITE(info_str,*)"***Error:- Digonalization of CFDM hamiltonian failed"
            CALL output_string(info_str)
        END IF
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! clean up
        !----------------------------------------------------------------------------------------------------------------------------------
        IF(ALLOCATED(cfdm_hamiltonian))             DEALLOCATE(cfdm_hamiltonian)
        IF(ALLOCATED(cfdm_hamiltonian_periodic))    DEALLOCATE(cfdm_hamiltonian_periodic) 
        IF(ALLOCATED(cfdm_eigenvalues))             DEALLOCATE(cfdm_eigenvalues)
        IF(ALLOCATED(coords_SL))                    DEALLOCATE(coords_SL)
        IF(ALLOCATED(omega_cfdm_SL))                DEALLOCATE(omega_cfdm_SL)
        IF(ALLOCATED(R_vdw_SL))                     DEALLOCATE(R_vdw_SL)
        IF(ALLOCATED(alpha_eff_SL))                 DEALLOCATE(alpha_eff_SL)  
        !----------------------------------------------------------------------------------------------------------------------------------
        RETURN
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE get_mbd_rSCS
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE get_alpha_omega_and_Rp(HF,C6_free,alpha_free,w,a_w,Rpw)
        implicit none
        real*8,intent(in)::HF,C6_free,alpha_free,w
        real*8,intent(out)::a_w,Rpw
        real*8::eff_freq,w_eff_freq_ratio
        ! alpha(iomega)
        eff_freq=0.0d0
        w_eff_freq_ratio=0.0d0
        eff_freq= ((4.d0/3.d0)*C6_free/(alpha_free**2.0))**2
        w_eff_freq_ratio=(w*w)/eff_freq
        a_w=(HF*alpha_free)/(1.0+w_eff_freq_ratio)
     
        !Rp(iomega) ! Rp_eff definition from   A. Mayer, Phys. Rev. B, 75,045407(2007)
        Rpw= ((a_w/3.d0)*dsqrt(2.0d0/pi))**0.333333333333333333333333333d0
        
        RETURN
    END SUBROUTINE get_alpha_omega_and_Rp
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE contract_matrix(tensor)
    !--------------------------------------------------------------------------------------------------------------------------------------
        IMPLICIT NONE
        !----------------------------------------------------------------------------------------------------------------------------------
        INTEGER                  :: ir
        INTEGER                  :: ic
        INTEGER                  :: i
        INTEGER                  :: j
        INTEGER                  :: i_row
        INTEGER                  :: i_col
        REAL*8, DIMENSION(3,3)   :: tensor
        !----------------------------------------------------------------------------------------------------------------------------------
        tensor(:,:)=0.0d0
        !----------------------------------------------------------------------------------------------------------------------------------
        DO ir=1,n_atoms
            DO ic=1,n_atoms
                i_row=0
                DO i=3*ir-2,3*ir
                    i_row=i_row+1
                    i_col=0
                    DO j=3*ic-2,3*ic
                        i_col=i_col+1
                        tensor(i_row,i_col)=tensor(i_row,i_col) + relay_matrix(i,j)
                    END DO
                END DO
            END DO
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
        RETURN
    !--------------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE contract_matrix
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE calculate_screened_polarizability(iso_polar_coupled)
    !--------------------------------------------------------------------------------------------------------------------------------------
        IMPLICIT NONE
        !----------------------------------------------------------------------------------------------------------------------------------
        INTEGER                                  :: i_row
        INTEGER                                  :: i_col
        INTEGER                                  :: i_index
        INTEGER                                  :: j_index
        INTEGER                                  :: errorflag
        INTEGER                                  :: LWORK
        REAL*8, DIMENSION(9)                     :: WORK
        REAL*8, DIMENSION(3)                     :: eigen
        REAL*8, DIMENSION(3,3)                   :: matrix
        REAL*8, DIMENSION(n_atoms), INTENT(OUT)  :: iso_polar_coupled
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        ! polar_coupled is formed by summing either rows 'blocks' or column 'blocks' of 
        ! relay_matrix_screened and contains full anisotropic atom resolved polarizability 
        !----------------------------------------------------------------------------------------------------------------------------------
        iso_polar_coupled=0.d0
        !----------------------------------------------------------------------------------------------------------------------------------
        
        !----------------------------------------------------------------------------------------------------------------------------------
        DO i_row=1,n_atoms
            matrix=0.0
            DO i_col=1,n_atoms,1
                DO i_index=1,3,1
                    DO j_index=1,3,1
                        matrix(i_index,j_index) = matrix(i_index,j_index) + relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)
                    END DO
                END DO
            END DO
            !------------------------------------------------------------------------------------------------------------------------------
            ! iso_polar_coupled contains average of trances of atom resolved polarizabilities in polar_coupled
            !------------------------------------------------------------------------------------------------------------------------------
            LWORK=9  
            CALL DSYEV('V','U',3,matrix,3,eigen,WORK,LWORK,errorflag)   ! JJ: why is this symmetric ??
            iso_polar_coupled(i_row) = sum(eigen)/3.d0
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
        RETURN
    END SUBROUTINE calculate_screened_polarizability
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE gauss_legendre_grid()
        casimir_omega( 0) = 0.0000000; casimir_omega_weight( 0) = 0.0000000
        casimir_omega( 1) = 0.0392901; casimir_omega_weight( 1) = 0.0786611
        casimir_omega( 2) = 0.1183580; casimir_omega_weight( 2) = 0.0796400
        casimir_omega( 3) = 0.1989120; casimir_omega_weight( 3) = 0.0816475
        casimir_omega( 4) = 0.2820290; casimir_omega_weight( 4) = 0.0847872
        casimir_omega( 5) = 0.3689190; casimir_omega_weight( 5) = 0.0892294
        casimir_omega( 6) = 0.4610060; casimir_omega_weight( 6) = 0.0952317
        casimir_omega( 7) = 0.5600270; casimir_omega_weight( 7) = 0.1031720
        casimir_omega( 8) = 0.6681790; casimir_omega_weight( 8) = 0.1136050
        casimir_omega( 9) = 0.7883360; casimir_omega_weight( 9) = 0.1273500
        casimir_omega(10) = 0.9243900; casimir_omega_weight(10) = 0.1456520
        casimir_omega(11) = 1.0817900; casimir_omega_weight(11) = 0.1704530
        casimir_omega(12) = 1.2684900; casimir_omega_weight(12) = 0.2049170
        casimir_omega(13) = 1.4966100; casimir_omega_weight(13) = 0.2544560
        casimir_omega(14) = 1.7856300; casimir_omega_weight(14) = 0.3289620
        casimir_omega(15) = 2.1691700; casimir_omega_weight(15) = 0.4480920
        casimir_omega(16) = 2.7106200; casimir_omega_weight(16) = 0.6556060
        casimir_omega(17) = 3.5457300; casimir_omega_weight(17) = 1.0659600
        casimir_omega(18) = 5.0273400; casimir_omega_weight(18) = 2.0635700
        casimir_omega(19) = 8.4489600; casimir_omega_weight(19) = 5.6851000
        casimir_omega(20) = 25.451700; casimir_omega_weight(20) = 50.955800
        RETURN
    END SUBROUTINE gauss_legendre_grid
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE get_vdw_param(atom,C6, alpha, R0)
        !----------------------------------------------------------------------------------------------------------------------------------
        ! Free-atom C6, polarizability and vdW radius
        !----------------------------------------------------------------------------------------------------------------------------------
        IMPLICIT NONE
        !----------------------------------------------------------------------------------------------------------------------------------
        CHARACTER*2 :: atom
        REAL*8      :: C6
        REAL*8      :: alpha
        REAL*8      :: R0
        !----------------------------------------------------------------------------------------------------------------------------------
        
        SELECT CASE (atom)
        
        case('H')
        alpha=4.500000
        C6=6.500000
        R0=3.100000
        
        case('He')
        alpha=1.380000
        C6=1.460000
        R0=2.650000
        
        case('Li')
        alpha=164.200000
        C6=1387.000000
        R0=4.160000
        
        case('Be')
        alpha=38.000000
        C6=214.000000
        R0=4.170000
        
        case('B')
        alpha=21.000000
        C6=99.500000
        R0=3.890000
        
        case('C')
        alpha=12.000000
        C6=46.600000
        R0=3.590000
        
        case('N')
        alpha=7.400000
        C6=24.200000
        R0=3.340000
        
        case('O')
        alpha=5.400000
        C6=15.600000
        R0=3.190000
        
        case('F')
        alpha=3.800000
        C6=9.520000
        R0=3.040000
        
        case('Ne')
        alpha=2.670000
        C6=6.380000
        R0=2.910000
        
        case('Na')
        alpha=162.700000
        C6=1556.000000
        R0=3.730000
        
        case('Mg')
        alpha=71.000000
        C6=627.000000
        R0=4.270000
        
        case('Al')
        alpha=60.000000
        C6=528.000000
        R0=4.330000
        
        case('Si')
        alpha=37.000000
        C6=305.000000
        R0=4.200000
        
        case('P')
        alpha=25.000000
        C6=185.000000
        R0=4.010000
        
        case('S')
        alpha=19.600000
        C6=134.000000
        R0=3.860000
        
        case('Cl')
        alpha=15.000000
        C6=94.600000
        R0=3.710000
        
        case('Ar')
        alpha=11.100000
        C6=64.300000
        R0=3.550000
        
        case('K')
        alpha=292.900000
        C6=3897.000000
        R0=3.710000
        
        case('Ca')
        alpha=160.000000
        C6=2221.000000
        R0=4.650000
        
        case('Sc')
        alpha=120.000000
        C6=1383.000000
        R0=4.590000
        
        case('Ti')
        alpha=98.000000
        C6=1044.000000
        R0=4.510000
        
        case('V')
        alpha=84.000000
        C6=832.000000
        R0=4.440000
        
        case('Cr')
        alpha=78.000000
        C6=602.000000
        R0=3.990000
        
        case('Mn')
        alpha=63.000000
        C6=552.000000
        R0=3.970000
        
        case('Fe')
        alpha=56.000000
        C6=482.000000
        R0=4.230000
        
        case('Co')
        alpha=50.000000
        C6=408.000000
        R0=4.180000
        
        case('Ni')
        alpha=48.000000
        C6=373.000000
        R0=3.820000
        
        case('Cu')
        alpha=42.000000
        C6=253.000000
        R0=3.760000
        
        case('Zn')
        alpha=40.000000
        C6=284.000000
        R0=4.020000
        
        case('Ga')
        alpha=60.000000
        C6=498.000000
        R0=4.190000
        
        case('Ge')
        alpha=41.000000
        C6=354.000000
        R0=4.200000
        
        case('As')
        alpha=29.000000
        C6=246.000000
        R0=4.110000
        
        case('Se')
        alpha=25.000000
        C6=210.000000
        R0=4.040000
        
        case('Br')
        alpha=20.000000
        C6=162.000000
        R0=3.930000
        
        case('Kr')
        alpha=16.800000
        C6=129.600000
        R0=3.820000
        
        case('Rb')
        alpha=319.200000
        C6=4691.000000
        R0=3.720000
        
        case('Sr')
        alpha=199.000000
        C6=3170.000000
        R0=4.540000
        
        case('Y')
        alpha=126.7370
        C6=1968.580
        R0=4.81510
        
        case('Zr')
        alpha=119.97
        C6=1677.91
        R0=4.53
        
        case('Nb')
        alpha=101.603
        C6=1263.61
        R0=4.2365
        
        case('Mo')
        alpha=88.4225785
        C6=1028.73
        R0=4.099
        
        case('Tc')
        alpha=80.083
        C6=1390.87
        R0=4.076
        
        case('Ru')
        alpha=65.8950
        C6=609.754
        R0=3.99530
        
        case('Rh')
        alpha=56.1
        C6=469.0
        R0=3.95
        
        case('Pd')
        alpha=23.680000
        C6=157.500000
        R0=3.66000
        
        case('Ag')
        alpha=50.600000
        C6=339.000000
        R0=3.820000
        
        case('Cd')
        alpha=39.7
        C6=452.0
        R0=3.99
        
        case('In')
        alpha=70.22000
        C6=707.046000
        R0=4.23198000
        
        case('Sn')
        alpha=55.9500
        C6=587.41700
        R0=4.303000
        
        case('Sb')
        alpha=43.671970 
        C6=459.322
        R0=4.2760
        
        case('Te')
        alpha=37.65
        C6=396.0
        R0=4.22
        
        case('I')
        alpha=35.000000
        C6=385.000000
        R0=4.170000
        
        case('Xe')
        alpha=27.300000
        C6=285.900000
        R0=4.080000
        
        case('Cs')
        alpha=427.12
        C6=6582.08
        R0=3.78
        
        case('Ba')
        alpha=275.0
        C6=5727.0
        R0=4.77
        
        case('Hf')
        alpha=99.52
        C6=1274.8
        R0=4.21
        
        case('Ta')
        alpha=82.53
        C6=1019.92
        R0=4.15
        
        case('W')
        alpha=71.041
        C6=847.93
        R0=4.08
        
        case('Re')
        alpha=63.04
        C6=710.2
        R0=4.02
        
        case('Os')
        alpha=55.055
        C6=596.67
        R0=3.84
        
        case('Ir')
        alpha=42.51
        C6=359.1
        R0=4.00
        
        case('Pt')
        alpha=39.68
        C6=347.1
        R0=3.92
        
        case('Au')
        alpha=36.5
        C6=298.0
        R0=3.86
        
        case('Hg')
        alpha=33.9
        C6=392.0
        R0=3.98
        
        case('Tl')
        alpha=69.92
        C6=717.44
        R0=3.91
        
        case('Pb')
        alpha=61.8
        C6=697.0
        R0=4.31
        
        case('Bi')
        alpha=49.02
        C6=571.0
        R0=4.32
        
        case('Po')
        alpha=45.013
        C6=530.92
        R0=4.097
        
        case('At')
        alpha=38.93
        C6=457.53
        R0=4.07
        
        case('Rn')
        alpha=33.54
        C6=390.63
        R0=4.23
        
        case default
        C6=0.0 
        alpha=1.0
        R0=1.0
        write(info_str,'(1X,4A)') '*** WARNING: VdW parameters not defined for atom: ', atom
        call output_string(info_str)
        stop
        
        END SELECT
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE get_vdw_param
    !--------------------------------------------------------------------------------------------------------------------------------------


    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE output_string(info_str)
        !----------------------------------------------------------------------------------------------------------------------------------
        CHARACTER*500 :: info_str  
        !----------------------------------------------------------------------------------------------------------------------------------
        IF(iproc .EQ. 0) WRITE(*,*)TRIM(info_str)
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE output_string
    !--------------------------------------------------------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------------------------------------------------------------
END MODULE UTILS
!------------------------------------------------------------------------------------------------------------------------------------------
