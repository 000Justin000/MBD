PROGRAM MBD_rsSCS
!-------------------------------------------------------------------------
    USE UTILS
    IMPLICIT NONE
    !---------------------------------------------------------------------
    ! --- local variables ---
    !---------------------------------------------------------------------
    CHARACTER*200                ::   geometry_file
    CHARACTER*200                ::   settings_file
    CHARACTER*200                ::   string
    CHARACTER*200                ::   dster
    CHARACTER*200                ::   temp
    INTEGER                      ::   io_file
    INTEGER                      ::   io_line
    INTEGER                      ::   i_index
    INTEGER                      ::   j_index
    REAL*8                       ::   ene_mbd_rsSCS
    REAL*8                       ::   start_time
    REAL*8                       ::   stop_time
    INTEGER, DIMENSION(7)        ::   input_settings
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! --- external functions ---
    !---------------------------------------------------------------------
    INTEGER                      ::   IARGC
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! initialize
    !---------------------------------------------------------------------
    geometry_file=""
    settings_file=""
    input_settings(:)=0
    CALL init_constants()
    CALL start_mpi() 
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! print parallelization info
    !---------------------------------------------------------------------
    IF (iproc==0) THEN
        WRITE(*,'(2X,A,I8,1X,A)') "Using  ", nproc,   "parallel tasks."
        WRITE(*,'(2X,A,I8,1X,A)') "Reading", IARGC(), "files."
    END IF
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! read geometry and settings filenames from command line arguments
    !---------------------------------------------------------------------
    IF (IARGC() .EQ. 0) THEN
        WRITE(*,*)"ERROR-Missing input file "
        WRITE(*,*)"Please specify input file"
        WRITE(*,*)"... Aborting current run"
        CALL stop_mpi()
        STOP
    END IF
    !---------------------------------------------------------------------
    IF (IARGC() .GE. 1) CALL GETARG(1,geometry_file)
    IF (IARGC() .GE. 2) CALL GETARG(2,settings_file)
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! read geometry from geometry_file
    !---------------------------------------------------------------------
    io_file=0
    OPEN(999,FILE=geometry_file(1:LEN_TRIM(geometry_file)),STATUS="old",IOSTAT=io_file,ERR=100)
    READ(999,*)n_atoms
    READ(999,*)
    !---------------------------------------------------------------------
    IF(.NOT. ALLOCATED(coords))                 ALLOCATE(coords(3,n_atoms))
    IF(.NOT. ALLOCATED(atom_name))              ALLOCATE(atom_name(n_atoms))
    IF(.NOT. ALLOCATED(hirshfeld_volume))       ALLOCATE(hirshfeld_volume(n_atoms))
    io_line=0  
    DO i_index=1,n_atoms,1
        READ(999,*,IOSTAT=io_line) atom_name(i_index), coords(:,i_index), hirshfeld_volume(i_index)
    END DO
    !---------------------------------------------------------------------
    io_line=0 
    READ(999,*,IOSTAT=io_line) string
    IF(INDEX(string,"lattice_vector") .GE. 1) THEN
        BACKSPACE(999)
        n_periodic     = 3       !n_periodic>0 indicate extended system  
        lattice_vector = 0.0d0
        DO i_index=1,3,1
            READ(999,*,IOSTAT=io_line) dster, lattice_vector(:,i_index)
        END DO
    END IF  
    CLOSE(999)
    !---------------------------------------------------------------------

    WRITE(*,*) "lattice vector:"
    WRITE(*,*) lattice_vector

    !---------------------------------------------------------------------
    ! read settings from settings_file
    !---------------------------------------------------------------------
    io_file=0
    IF(IARGC() .GE. 2) THEN
        OPEN(111,FILE=settings_file(1:LEN_TRIM(settings_file)),STATUS="old",IOSTAT=io_file,ERR=100)
        io_line=0  
        DO WHILE (io_line .EQ. 0)
            READ(111,'(A)',IOSTAT=io_line) string
            IF(INDEX(string,"xc") .GE. 1) THEN
                BACKSPACE(111)
                input_settings(1)=1
                READ(111,*,IOSTAT=io_line) dster, flag_xc
            ELSE IF(INDEX(string,"mbd_supercell_cutoff") .GE. 1) THEN
                BACKSPACE(111)
                input_settings(2)=1
                READ(111,*,IOSTAT=io_line) dster, mbd_supercell_cutoff
            ELSE IF(INDEX(string,"mbd_cfdm_dip_cutoff") .GE. 1) THEN
                BACKSPACE(111)
                input_settings(3)=1
                READ(111,*,IOSTAT=io_line) dster, mbd_cfdm_dip_cutoff
            ELSE IF(INDEX(string,"mbd_scs_dip_cutoff") .GE. 1) THEN
                BACKSPACE(111)
                input_settings(4)=1
                READ(111,*,IOSTAT=io_line) dster, mbd_scs_dip_cutoff
            ELSE IF(INDEX(string,"mbd_scs_vacuum_axis") .GE. 1) THEN
                BACKSPACE(111)
                input_settings(5)=1
                READ(111,*,IOSTAT=io_line) dster, mbd_scs_vacuum_axis(:)
            ELSE IF(INDEX(string,"nprow") .GE. 1) THEN
                BACKSPACE(111)
                input_settings(6)=1
                READ(111,*,IOSTAT=io_line) dster, nprow
            ELSE IF(INDEX(string,"npcol") .GE. 1) THEN
                BACKSPACE(111)
                input_settings(7)=1
                READ(111,*,IOSTAT=io_line) dster, npcol
            END IF
        END DO
        CLOSE(111)
    END IF
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! print geometry and settings info
    !---------------------------------------------------------------------
    IF (iproc==0) THEN  
        WRITE(*,'(3x,A)')"|----------------------------------Input Geometry----------------------------------" 
        WRITE(*,'(3x,A,I5)') "Number of atoms",n_atoms
        WRITE(*,*) 
        DO i_index=1,n_atoms,1
            WRITE(*,'(4x,A,4(2x,f15.6))') atom_name(i_index), coords(:,i_index), hirshfeld_volume(i_index)
        END DO
        WRITE(*,'(3x,A)')"|----------------------------------------------------------------------------------"
    END IF
    !---------------------------------------------------------------------
    IF (iproc==0)THEN
        WRITE(*,'(3x,A)')"|--------------------------------MBD@rSCS Settings---------------------------"
        IF(input_settings(1) .EQ. 1)THEN
            IF(flag_xc .EQ. 1) THEN
                WRITE(temp,'(A)') "  (PBE)"     
            ELSE IF(flag_xc .EQ. 2) THEN
                WRITE(temp,'(A)') "  (PBE0)"     
            ELSE IF(flag_xc .EQ. 3) THEN
                WRITE(temp,'(A)') "  (HSE)"     
            END IF 
            WRITE(*,'(4x,A,I4,A)') "Exchange-Correlation functional:", flag_xc, TRIM(temp)
        ELSE IF(flag_xc .EQ. 1) THEN
            WRITE(temp,'(A)') "(PBE)"
        ELSE IF(flag_xc .EQ. 2) THEN
            WRITE(temp,'(A)') "(PBE0)"
        ELSE IF(flag_xc.eq.3) THEN
            WRITE(temp,'(A)') "(HSE)"
        END IF
        WRITE(*,'(4x,A)') "xc type not specified in setting file, using default for "
        WRITE(*,'(4x,A,I4,A)') "Exchange-Correlation functional:", flag_xc, TRIM(temp)
    END IF  
    !---------------------------------------------------------------------
    IF(n_periodic .GT. 0) THEN 
        IF(input_settings(2) .EQ. 1) THEN
            WRITE(*,'(4x,A,f10.5,A)') "Radius used to construct the supercell in periodic MBD calculations", mbd_supercell_cutoff," Angstrom"
        ELSE
            WRITE(*,'(4x,A)'        ) "mbd_supercell_cutoff  not specified in setting file, using default for "
            WRITE(*,'(4x,A,f10.5,A)') "Radius used to construct the supercell in periodic MBD calculations", mbd_supercell_cutoff," Angstrom"
        END IF

        IF(input_settings(3) .EQ. 1) THEN
            WRITE(*,'(4x,A,f10.5,A)') "Radius used to integrate dipole field in periodic MBD calculations", mbd_cfdm_dip_cutoff," Angstrom"
        ELSE
            WRITE(*,'(4x,A)'        ) "mbd_cfdm_dip_cutoff  not specified in setting file, using default for "
            WRITE(*,'(4x,A,f10.5,A)') "Radius used to integrate dipole field in periodic MBD calculations", mbd_cfdm_dip_cutoff," Angstrom"
        END IF

        IF(input_settings(4) .EQ. 1) THEN
            WRITE(*,'(4x,A,f10.5,A)') "Radius used to integrate dipole field in periodic SCS calculations", mbd_scs_dip_cutoff," Angstrom"
        ELSE
            WRITE(*,'(4x,A)'        ) "mbd_scs_dip_cutoff  not specified in setting file, using default for "
            WRITE(*,'(4x,A,f10.5,A)') "Radius used to integrate dipole field in periodic SCS calculations", mbd_scs_dip_cutoff," Angstrom"
        END IF

        IF(input_settings(5).eq.1) THEN
            WRITE(*,'(4x,A)') " Found mbd_scs_vacuum_axis keyword in setting file" 
            WRITE(*,'(4x,A)') " This keyword specifies direction of vacuum region in case of low-dimensional systems"
            WRITE(*,'(4x,A)') "Default: No vacuum mbd_scs_vacuum_axis .false. .false. .false."
            WRITE(*,'(4x,A)') "For example in the case of periodic slab along xy directions the keyword"
            WRITE(*,'(4x,A)') "mbd_scs_vacuum_axis .false. .false. .true. is needed to specify z as"
            WRITE(*,'(4x,A)') "vacuum direction in MBD/SCS calculation."
        END IF
        
        WRITE(*,'(3x,A)')"|---------------------------------------------------------------------------"
    END IF  
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! scale the coordinates/lattice_vector to atomic unit
    !---------------------------------------------------------------------
    coords = coords/bohr
    lattice_vector = lattice_vector/bohr
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! major calculation
    !---------------------------------------------------------------------
    CALL init_blacs()
    !---------------------------------------------------------------------
    start_time = MPI_Wtime()
    CALL MBD_at_rsSCS(ene_mbd_rsSCS)
    stop_time  = MPI_Wtime()
    !---------------------------------------------------------------------
    IF (iproc==0) WRITE(*,'(A,F18.10)') "Total Calculation time: ", stop_time-start_time
    !---------------------------------------------------------------------

    !---------------------------------------------------------------------
    ! clean up
    !---------------------------------------------------------------------
    IF(ALLOCATED(coords))             DEALLOCATE(coords)
    IF(ALLOCATED(atom_name))          DEALLOCATE(atom_name)
    IF(ALLOCATED(hirshfeld_volume))   DEALLOCATE(hirshfeld_volume)
    CALL stop_mpi()
    !---------------------------------------------------------------------
    STOP
    !---------------------------------------------------------------------
100 IF (io_file .NE. 0) WRITE(*,*) "Error in reading input file"
!-------------------------------------------------------------------------
END PROGRAM MBD_rsSCS
