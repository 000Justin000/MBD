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
    REAL*8                                      ::   three_by_pi 
    !--------------------------------------------------------------------------------------------------------------------------------------
    ! input setting
    !--------------------------------------------------------------------------------------------------------------------------------------
    INTEGER                                     ::   flag_xc
    INTEGER                                     ::   n_periodic
    REAL*8                                      ::   mbd_cfdm_dip_cutoff 
    REAL*8                                      ::   mbd_supercell_cutoff  
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
    INTEGER                                     ::   myid     
    INTEGER                                     ::   mpiierror                    ! mpi_error_handler 
    INTEGER                                     ::   n_tasks                      ! number of theads    
    INTEGER                                     ::   n_atoms                      ! number of atoms
    !--------------------------------------------------------------------------------------------------------------------------------------

CONTAINS

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE init_constants()
        !----------------------------------------------------------------------------------------------------------------------------------
        flag_xc                = 1                                                ! XC= PBE
        mbd_cfdm_dip_cutoff    = 100.d0                                           ! Angstrom
        mbd_supercell_cutoff   = 25.d0                                            ! Angstrom
        mbd_scs_dip_cutoff     = 120.0                                            ! Angstrom
        n_periodic             = 0  
        bohr                   = 0.52917721d0
        pi                     = 3.14159265358979323846d0
        three_by_pi            = 0.954929658551372d0
        mbd_scs_vacuum_axis(1) = .false. 
        mbd_scs_vacuum_axis(2) = .false. 
        mbd_scs_vacuum_axis(3) = .false. 
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE init_constants
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE start_mpi()
        !----------------------------------------------------------------------------------------------------------------------------------
        CALL MPI_INIT (mpiierror)
        CALL MPI_COMM_RANK (MPI_COMM_WORLD, myid,    mpiierror)
        CALL MPI_COMM_SIZE (MPI_COMM_WORLD, n_tasks, mpiierror)
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE stop_mpi()
        !----------------------------------------------------------------------------------------------------------------------------------
        CALL MPI_FINALIZE(mpiierror)
        IF(mpiierror.ne.0) WRITE(*,*) "ERROR in terminating MPI"
        IF(mpiierror.ne.0) WRITE(*,*) "Normal Termination"
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
            task_list(i_row) = MOD(i_row,n_tasks)
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
    END SUBROUTINE allocate_task
    !--------------------------------------------------------------------------------------------------------------------------------------

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE sync_tensors(temp_matrix,num_element)
        !----------------------------------------------------------------------------------------------------------------------------------
        INTEGER                                      :: num_element
        REAL*8, DIMENSION(num_element,num_element)   :: temp_matrix
        REAL*8, DIMENSION(num_element,num_element)   :: temp_matrix_mpi
        !----------------------------------------------------------------------------------------------------------------------------------
        call MPI_ALLREDUCE(temp_matrix, temp_matrix_mpi, num_element*num_element, MPI_DOUBLE_PRECISION, MPI_SUM,MPI_COMM_WORLD, mpiierror)
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
        IF(n_periodic .eq. 0) THEN
            WRITE(info_str,'(2x,A)')"| Dynamic molecular polarizability alpha(iw) (bohr^3)"
            CALL output_string(info_str)
        ELSE
            WRITE(info_str,'(2x,A)')"| Dynamic polarizability of unit cell alpha(iw) (bohr^3)"
            CALL output_string(info_str)
        END IF
        WRITE(info_str,'(2x,A)')"| ---------------------------------------------------------------------------"
        CALL output_string(info_str)
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
            ! 
            !------------------------------------------------------------------------------------------------------------------------------
            DO i_myatom=1,n_atoms
                CALL get_vdw_param(atom_name(i_myatom),C6_free,alpha_free,R_vdw_free)
                CALL get_alpha_omega_and_Rp(hirshfeld_volume(i_myatom),C6_free,alpha_free,casimir_omega(i_freq),alpha_omega(i_myatom),R_p(i_myatom))
                Rvdw_iso(i_myatom)= (hirshfeld_volume(i_myatom)**0.333333333333333333333333333d0)*R_vdw_free
            END DO
            !------------------------------------------------------------------------------------------------------------------------------
            
            !------------------------------------------------------------------------------------------------------------------------------
            ! compute screened polarizability using isotropically damped tensor matrix(INPUT:-alpha_omega,R_p,Rvdw_iso), Ambrosetti et al 
            !------------------------------------------------------------------------------------------------------------------------------
            CALL calculate_scs_matrix() 
            
            ! Contract relay_matrix^-1 for polarizability
            call contract_matrix(mol_pol_tensor)  
            
            ! get coupled atomic "local" polarizability by summing over row or columns of
            ! relay tensor
            if (i_freq.eq.0) then
               ! Store static polarizability for post SCS routines  
               call calculate_screened_polarizability(alpha_eff)
            endif
            
            if (i_freq.ge.1) then
               call calculate_screened_polarizability(coupled_atom_pol(i_freq,:))
            endif 
            
            ! Casimir-Polder integral for molecular C6-coefficent 
            LWORK=9 
            call DSYEV('V','U',3,mol_pol_tensor,3,mol_pol_eigen,WORK,LWORK,errorflag) 
            if(i_freq.ge.1)  then
                mol_c6_coeff = mol_c6_coeff + (casimir_omega_weight(i_freq)* (sum(mol_pol_eigen)/3.0)**2)
            endif
            
            write(info_str,'(2x,"| ",F10.6,4(e15.6))')casimir_omega(i_freq),&
                   sum(mol_pol_eigen)/3.0,mol_pol_eigen(1),mol_pol_eigen(2),mol_pol_eigen(3)
            call output_string(info_str)  
        !----------------------------------------------------------------------------------------------------------------------------------
        END DO
        !----------------------------------------------------------------------------------------------------------------------------------
        
        
              write(info_str,'(2x,A)')&
              "| ---------------------------------------------------------------------------"
              call output_string(info_str)  

              if(n_periodic .eq. 0) then   
              Write(info_str,'(2x,A,f20.8,2x,A)')&
              "| Molecular C6 coefficient    :  ",0.954929658551372d0*mol_c6_coeff,"  hartree bohr^6"
              call output_string(info_str)  
              else
              Write(info_str,'(2x,A,f20.8,2x,A)')&
              "| Unit cell C6 coefficient    :  ",0.954929658551372d0*mol_c6_coeff,"  hartree bohr^6"
              call output_string(info_str)  
              endif

              write(info_str,'(2x,A)')&
              "| ---------------------------------------------------------------------------"
              call output_string(info_str)  
              if(n_periodic .eq. 0) then
              write(info_str,'(2x,A)')"| Partitioned C6 (hartree bohr^6) | alpha (bohr^3)  of atom in molecule"
              call output_string(info_str)  
              else
              write(info_str,'(2x,A)')"| Partitioned C6 (hartree bohr^6) | alpha (bohr^3)  of atom in unit cell"
              call output_string(info_str)  
           
              endif  
              write(info_str,'(2x,A)')&
              "| ---------------------------------------------------------------------------"
              call output_string(info_str)  

              do i_myatom=1,n_atoms,1
                 C6_atom = 0.0d0 
                 do i_freq=1,20,1
                 C6_atom = C6_atom + (casimir_omega_weight(i_freq)*coupled_atom_pol(i_freq,i_myatom)**2)
                 enddo
                 write(info_str,'(2x,A,I4,2x,A,f12.6,6x,f12.6)')"| ATOM ",i_myatom,&
                             atom_name(i_myatom), C6_atom*0.954929658551372d0,alpha_eff(i_myatom)
                 call output_string(info_str)  

                 !store effctive C6 for post SCS routines 
                  C6_eff(i_myatom)=C6_atom*0.954929658551372d0
              enddo 


              write(info_str,'(2x,A)')&
              "| ---------------------------------------------------------------------------"
              call output_string(info_str)

              !  
              ! Compute many-body-dispersion energy new alpha_eff and C6_eff

              call get_mbd_rSCS(ene_mbd_rsSCS)
              write(info_str,'(2X,A,F20.9,A,F20.9,A)') &
              "| MBD@rsSCS energy   :", ene_mbd_rsSCS," Ha     ",ene_mbd_rsSCS*27.211," eV"
              call output_string(info_str)


 
              if(allocated(R_p))                   deallocate(R_p)
              if(allocated(Rvdw_iso))              deallocate(Rvdw_iso)
              if(allocated(alpha_omega))           deallocate(alpha_omega)
              if(allocated(relay_matrix))          deallocate(relay_matrix)
              if(allocated(relay_matrix_periodic)) deallocate(relay_matrix_periodic)
              if(allocated(coupled_atom_pol))      deallocate(coupled_atom_pol)
              if(allocated(alpha_eff))             deallocate(alpha_eff)
              if(allocated(C6_eff))                deallocate(C6_eff)



      return
    END SUBROUTINE MBD_at_rsSCS

    !--------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE calculate_scs_matrix()
    !--------------------------------------------------------------------------------------------------------------------------------------
        integer ::i_index,j_index
        real*8,dimension(3,3)::TPP
        real*8,dimension(3) :: dxyz
        real*8,dimension(3) :: coord_curr
        real*8 :: r_ij
        real*8 :: r_pp
        real*8 :: Rvdw12
        real*8 :: beta
        integer :: i_row, i_col
        integer :: i_lattice, j_lattice, k_lattice
        integer :: errorflag,periodic_cell_i,periodic_cell_j,periodic_cell_k

    !For LAPACK
    integer,dimension(3*n_atoms):: IPIV
    real*8,dimension(3*n_atoms):: WORK

    ! initio values
    relay_matrix=0.0d0

        select case (flag_xc)
          case (1) ! PBE
               beta=0.83
          case (2) !PBE0
               beta=0.85
          case (3) !HSE 
               beta=0.85

          case default
               beta=1.0
               write(info_str,'(A)')"***  WARNING range seperated parameter beta is not defined for"
               call output_string(info_str)
               write(info_str,'(A)')"this fxc defaulting it to 0.0 , which will give MBD energy"
               call output_string(info_str)
               write(info_str,'(A)')"WITHOUT short range damping"
               call output_string(info_str)
        end select



     ! compute relay matrix of  cluster or unit cell
       do i_row=1,n_atoms,1 !#1
         if(myid.eq.task_list(i_row)) then      
         do i_col=i_row,n_atoms,1 !#2
         TPP=0.d0
            if(i_row.eq.i_col) then  !$1
               do i_index=1,3,1
                  do j_index=1,3,1
                     if(i_index.eq.j_index) then
                        relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=1.d0/alpha_omega(i_row)
                     else 
                        relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=0.d0
                     endif   
                  enddo
               enddo

            else
               dxyz(:) = coords(:,i_col)-coords(:,i_row)
               r_ij = dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
               r_pp = dSqrt(R_p(i_row)**2 + R_p(i_col)**2)
               Rvdw12 = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
               call SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)
               do i_index=1,3,1
                  do j_index=1,3,1
                   relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)=TPP(i_index,j_index)
                   relay_matrix(3*i_col-3+j_index,3*i_row-3+i_index)=TPP(i_index,j_index)
                  enddo
               enddo

            endif !$1
         enddo   !#2
        endif ! task 
       enddo  !#1
  call sync_tensors(relay_matrix,3*n_atoms)
  if (n_periodic .gt. 0) then
   
      relay_matrix_periodic=0.d0 

      if(.NOT.mbd_scs_vacuum_axis(1)) then     
      periodic_cell_i = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,1)**2 +& 
                                                               lattice_vector(2,1)**2 +&
                                                               lattice_vector(3,1)**2)) 
      else
      periodic_cell_i = 0
      endif 
      if(.NOT.mbd_scs_vacuum_axis(2)) then     
      periodic_cell_j = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,2)**2 +& 
                                                               lattice_vector(2,2)**2 +&
                                                               lattice_vector(3,2)**2)) 
      else
      periodic_cell_j = 0
      endif 
      if(.NOT.mbd_scs_vacuum_axis(3)) then     
      periodic_cell_k = ceiling((mbd_scs_dip_cutoff/bohr)/SQRT(lattice_vector(1,3)**2 +& 
                                                               lattice_vector(2,3)**2 +&
                                                               lattice_vector(3,3)**2))
      else
      periodic_cell_k = 0
      endif 
            
      do i_lattice = -periodic_cell_i, periodic_cell_i,1
         do j_lattice = -periodic_cell_j, periodic_cell_j,1
            do k_lattice = -periodic_cell_k, periodic_cell_k,1
               if((abs(i_lattice).ne.0).or.(abs(j_lattice).ne.0).or.(abs(k_lattice).ne.0))then !#1
                  do i_row = 1, n_atoms, 1 ! atom1 loop
                    if(myid.eq.task_list(i_row)) then
                     do i_col = i_row, n_atoms, 1 ! atom2 loop
                          coord_curr = 0.d0
                          dxyz = 0.d0
                                r_ij = 0.d0
                                TPP  = 0.d0
                               Rvdw12= 0.d0  
                          ! find the coordinate of images
                          coord_curr(:) = coords(:,i_col) + i_lattice*lattice_vector(:,1) + &
                                                            j_lattice*lattice_vector(:,2) + &
                                                            k_lattice*lattice_vector(:,3)

                          dxyz(:) = coords(:,i_row)- coord_curr(:)
                          r_ij    = sqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
                          r_pp    = sqrt(R_p(i_row)**2 + R_p(i_col)**2)
                          Rvdw12  = Rvdw_iso(i_row) +  Rvdw_iso(i_col)
                          if(r_ij.le.mbd_scs_dip_cutoff/bohr) then
                            call SCS_TENSOR_MBD_rsSCS(dxyz,r_ij,r_pp,Rvdw12,beta,TPP)
                            do i_index=1,3,1
                            do j_index=1,3,1
                            relay_matrix_periodic(3*i_row-3+i_index,3*i_col-3+j_index)=&
                            relay_matrix_periodic(3*i_row-3+i_index,3*i_col-3+j_index) + TPP(i_index,j_index)
                             if(i_col.NE.i_row) then
                               relay_matrix_periodic(3*i_col-3+j_index,3*i_row-3+i_index)=&
                               relay_matrix_periodic(3*i_col-3+j_index,3*i_row-3+i_index) + TPP(i_index,j_index)
                              endif
                            enddo
                            enddo
                          endif
                     enddo !atom2 loop
                   endif! parallel
                  enddo !atom1 loop
               endif  !#1
            enddo
         enddo
      enddo
   call sync_tensors(relay_matrix_periodic,3*n_atoms)
   relay_matrix = relay_matrix + relay_matrix_periodic  
   endif

   call DGETRF(3*n_atoms, 3*n_atoms, relay_matrix, 3*n_atoms, IPIV, errorflag)
   if(errorflag.ne.0) then
   write(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
   call output_string(info_str)
   endif

   call DGETRI(3*n_atoms, relay_matrix, 3*n_atoms, IPIV, WORK,3*n_atoms,errorflag )
   if(errorflag.ne.0) then
   write(info_str,'(A)')"Error** Matrix inversion failed in SCS module"
   call output_string(info_str)
   endif
   return
   endsubroutine  calculate_scs_matrix

   subroutine SCS_TENSOR_MBD_rsSCS(dxyz, r_ij, r_pp,Rvdw12,beta, TPP)
      real*8,dimension(3),intent(in) :: dxyz
      real*8,intent(in) :: r_ij
      real*8,intent(in) :: r_pp
      real*8,intent(in) :: Rvdw12
      real*8,intent(in) :: beta
      real*8,dimension(3,3),intent(out) :: TPP
      ! This needs declarding for the PGI Architecture
      real*8 :: derf, erf

      ! local vars
      real*8,dimension(3,3) :: TPQ,r_tensor
      real*8:: zeta_l,zeta_r,ratio,d_param,fermi_param
      integer :: i,j

      !o dipole-dipole interaction tensor between two quantum harmonic
      !  oscilators Tp-p damped with Fermi damping 
      TPP(:,:)=0.d0
      ! Fermi damping parameter 
      d_param = 6.0
      fermi_param = r_ij/(beta*Rvdw12)
      !!!!!!!!!!!!!!!!!!!!!!!!!  
      ratio=r_ij/r_pp
      zeta_l=(derf(ratio)-(dsqrt(4.0d0/pi) * ratio*dexp(-(ratio**2.0d0))))/r_ij**5.0d0
      zeta_r=(dsqrt(16.0d0/pi) * dexp(-(ratio**2.0d0)))/(r_pp * r_pp * r_pp*r_ij * r_ij)

                         !!Tensor product
                          do i=1,3,1
                            do j=1,3,1
                              r_tensor(i,j)=dxyz(i)*dxyz(j)
                            enddo
                          enddo

                          TPQ=r_tensor*3.0d0
                          do i=1,3,1
                             TPQ(i,i)=TPQ(i,i)-(r_ij*r_ij)
                          enddo

                          TPQ=TPQ*zeta_l
                          r_tensor=r_tensor * zeta_r
                          TPP=TPQ-r_tensor
                          TPP=-1.d0*TPP
                          ! Fermi type damping  
                          TPP=TPP*(1.0-(1.0/(1.0+exp(-d_param*(fermi_param-1.0)))))
      return
endsubroutine SCS_TENSOR_MBD_rsSCS

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


subroutine get_mbd_rSCS(ene_mbd_rsSCS)

     ! local vars  
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian 
     real*8,dimension(:,:),allocatable:: cfdm_hamiltonian_periodic
     real*8,dimension(:,:),allocatable:: coords_SL
     real*8,dimension(:),allocatable:: cfdm_eigenvalues
     real*8,dimension(:),allocatable:: omega_cfdm_SL
     real*8,dimension(:),allocatable:: R_vdw_SL
     real*8,dimension(:),allocatable::alpha_eff_SL
     real*8,dimension(:),allocatable:: WORK
     real*8,dimension(:),allocatable:: task_list_SL  
     real*8,dimension(3,3):: TPP,lattice_vector_SL
     real*8,dimension(3)::dxyz,coord_curr
     real*8:: r_ij
     real*8:: C6_free
     real*8:: alpha_free
     real*8:: R_vdw_free
     real*8:: Rvdw_12
     real*8:: beta
     real*8:: CFDM_prefactor
     real*8:: E_int
     real*8:: E_nonint
     real*8:: sigma
     real*8:: ene_mbd_rsSCS
     integer :: errorflag,i_atom,j_atom
     integer :: SL_i,SL_j,SL_k ! MBD super cell index's
     integer :: i_index, j_index ! 3x3 block index
     integer :: i_lattice,j_lattice,k_lattice ! lattice increament index
     integer :: periodic_cell_i,periodic_cell_j,periodic_cell_k
     integer :: NIMAG
     integer ::n_atoms_SL
     integer :: LWORK
!     Begin work       
        select case (flag_xc)
          case (1) ! PBE
               beta=0.83
          case (2) !PBE0
               beta=0.85
          case (3) !HSE 
               beta=0.85
          case default
               beta=1.0
               write(info_str,'(A)')"***  WARNING range seperated parameter beta is not defined for"
               call output_string(info_str)
               write(info_str,'(A)')"this fxc defaulting it to 0.0 , which will give MBD energy"
               call output_string(info_str)
               write(info_str,'(A)')"WITHOUT short range damping"
               call output_string(info_str)
        end select

      if ( n_periodic .eq. 0) then
      ! For Cluster calculation NO super cell needed 
      ! intiate this index so as to avoide any invalid floting point when
      ! normalizing energy in case of cluster/molecule
      SL_i = 1 
      SL_j = 1
      SL_k = 1

!     number of atoms in cluster /molecule
      n_atoms_SL = n_atoms
 
      if(.NOT.allocated(cfdm_hamiltonian))          allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_eigenvalues))          allocate(cfdm_eigenvalues(3*n_atoms_SL)) 
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))
      if(.NOT.allocated(WORK))                      allocate(WORK((3*n_atoms_SL)*(3+(3*n_atoms_SL)/2)))
      if(.NOT.allocated(task_list_SL))              allocate(task_list_SL(n_atoms_SL))

      ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
      cfdm_hamiltonian=0.0d0 
      cfdm_eigenvalues=0.0d0

!     RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN CLUSTER/MOLECULE
      do i_atom=1,n_atoms_SL
        task_list_SL(i_atom)   = MOD(i_atom,n_tasks)  
      enddo

          do i_atom =1,n_atoms
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0
                  C6_free=0.d0
                  call get_vdw_param(atom_name(i_atom),C6_free,alpha_free,R_vdw_free)

                  R_vdw_SL(i_atom)= R_vdw_free*((alpha_eff(i_atom)/alpha_free)**0.333333333333333333333333333d0)
                  omega_cfdm_SL(i_atom) = (4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
                  coords_SL(:,i_atom)   = coords(:,i_atom)
                  alpha_eff_SL(i_atom) = alpha_eff(i_atom)

         enddo
      else

      if(.NOT.mbd_scs_vacuum_axis(1)) then
      SL_i = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,1)**2 + lattice_vector(2,1)**2 +lattice_vector(3,1)**2)) 
      else
      SL_i = 1 
      endif

      if(.NOT.mbd_scs_vacuum_axis(2)) then
      SL_j = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,2)**2 + lattice_vector(2,2)**2 +lattice_vector(3,2)**2)) 
      else
      SL_j = 1 
      endif

      if(.NOT.mbd_scs_vacuum_axis(3)) then
      SL_k = ceiling((mbd_supercell_cutoff/bohr)/SQRT(lattice_vector(1,3)**2 + lattice_vector(2,3)**2 +lattice_vector(3,3)**2)) 
      else
      SL_k = 1 
      endif
       
!     number of atoms in super cell 
      n_atoms_SL = n_atoms*SL_i*SL_j*SL_k 
      write(info_str,'(2X,A,i3,A,i3,A,i3,A)') &
         "| Creating super cell of dimension",  SL_i," X ", SL_j," X ", &
         SL_k, " in MBD calculation" 
      call output_string(info_str)
      write(info_str,'(2x,"| containing", I6,A)')n_atoms_SL, "  atom"
      call output_string(info_str)
        


      if(.NOT.allocated(cfdm_hamiltonian))          allocate(cfdm_hamiltonian(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_hamiltonian_periodic)) allocate(cfdm_hamiltonian_periodic(3*n_atoms_SL,3*n_atoms_SL)) 
      if(.NOT.allocated(cfdm_eigenvalues))          allocate(cfdm_eigenvalues(3*n_atoms_SL)) 
      if(.NOT.allocated(coords_SL))                 allocate(coords_SL(3,n_atoms_SL)) 
      if(.NOT.allocated(omega_cfdm_SL))             allocate(omega_cfdm_SL(n_atoms_SL))
      if(.NOT.allocated(R_vdw_SL))                  allocate(R_vdw_SL(n_atoms_SL))
      if(.NOT.allocated(alpha_eff_SL))              allocate(alpha_eff_SL(n_atoms_SL))
      if(.NOT.allocated(WORK))                      allocate(WORK(3*n_atoms_SL*(3+(3*n_atoms_SL)/2)))
      if(.NOT.allocated(task_list_SL))              allocate(task_list_SL(n_atoms_SL))
      ! MUST BE INTIATED TO ZERO TO AVOID ANY POSSIBLE noise in LAPACK call
      cfdm_hamiltonian = 0.0d0 
      cfdm_eigenvalues = 0.0d0
      cfdm_hamiltonian_periodic = 0.d0 
!     RECOMPUTE TASK LIST DEPENDING ON NUMBER OF ATOMS IN SUPERCELL/MOLECULE
      do i_atom=1,n_atoms_SL 
        task_list_SL(i_atom)   = MOD(i_atom,n_tasks) 
      enddo

      ! Get all the parameter  required for MBD with super cell calculation
       j_atom = 0
       do i_lattice = 0, SL_i-1,1
         do j_lattice = 0, SL_j-1,1
            do k_lattice = 0, SL_k-1,1
              do i_atom =1,n_atoms !<--- atom index
                  j_atom = j_atom +1  !<--- dummy index supercell
                  coords_SL(:,j_atom) = coords(:,i_atom) + i_lattice*lattice_vector(:,1) + &
                                                           j_lattice*lattice_vector(:,2) + &
                                                           k_lattice*lattice_vector(:,3)   
                  R_vdw_free=0.0d0
                  alpha_free=0.0d0  
                  call get_vdw_param(atom_name(i_atom),C6_free,alpha_free,R_vdw_free) !FIXME

                  R_vdw_SL(j_atom)=R_vdw_free*((alpha_eff(i_atom)/alpha_free)**0.333333333333333333333333333d0)
                  alpha_eff_SL(j_atom) = alpha_eff(i_atom)
                  omega_cfdm_SL(j_atom) =(4.d0/3.d0)*C6_eff(i_atom)/(alpha_eff(i_atom)**2.0)
              enddo 
            enddo            
         enddo            
       enddo            
      lattice_vector_SL(:,1)  = lattice_vector(:,1)*SL_i 
      lattice_vector_SL(:,2)  = lattice_vector(:,2)*SL_j 
      lattice_vector_SL(:,3)  = lattice_vector(:,3)*SL_k 
      endif 

      ! Construct CFDM hamiltonian matrix 
      do i_atom=1,n_atoms_SL,1 !$1
         if(myid.eq.task_list_SL(i_atom)) then
         do j_atom=i_atom,n_atoms_SL,1 !$2

         if(i_atom.eq.j_atom) then  !#1

               do i_index=1,3,1
                  do j_index=1,3,1
                  if(i_index.eq.j_index) then
                   cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = omega_cfdm_SL(i_atom)**2.0
                  endif  
                  enddo
               enddo

         else
            r_ij=0.0d0
            TPP=0.0d0
            dxyz=0.d0
            dxyz(:)= coords_SL(:,i_atom)-coords_SL(:,j_atom)
            r_ij=dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )
            Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
            sigma=(r_ij/Rvdw_12)**beta
            call MBD_TENSOR_MBD_rsSCS(dxyz,r_ij,Rvdw_12,beta,TPP)
            CFDM_prefactor=omega_cfdm_SL(i_atom)*omega_cfdm_SL(j_atom)*dsqrt(alpha_eff_SL(i_atom)*alpha_eff_SL(j_atom))

                   ! Transfer each dipole matrix to CFDM hamiltonian i.j accordingly         
               do i_index=1,3,1
                  do j_index=1,3,1
                   cfdm_hamiltonian(3*i_atom-3+i_index,3*j_atom-3+j_index) = TPP(i_index,j_index)*CFDM_prefactor
                   cfdm_hamiltonian(3*j_atom-3+j_index,3*i_atom-3+i_index) = TPP(i_index,j_index)*CFDM_prefactor
                  enddo
               enddo

         endif !#1
         enddo !$2
       endif ! tasks 
      enddo !$1
      call sync_tensors(cfdm_hamiltonian,3*n_atoms_SL)

! Adds dipole field due to image cells based spherical cutoff mbd_cfdm_dip_cutoff
if (n_periodic .gt. 0) then 
      if(.NOT.mbd_scs_vacuum_axis(1)) then
      periodic_cell_i = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,1)**2 +&
                                                                lattice_vector_SL(2,1)**2 +&
                                                                lattice_vector_SL(3,1)**2))
      else
      periodic_cell_i = 0 
      endif   

      if(.NOT.mbd_scs_vacuum_axis(2)) then
      periodic_cell_j = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,2)**2 +&
                                                                lattice_vector_SL(2,2)**2 +&
                                                                lattice_vector_SL(3,2)**2)) 
      else
      periodic_cell_j = 0 
      endif   

      if(.NOT.mbd_scs_vacuum_axis(3)) then
      periodic_cell_k = ceiling((mbd_cfdm_dip_cutoff/bohr)/SQRT(lattice_vector_SL(1,3)**2 +&
                                                                lattice_vector_SL(2,3)**2 +&
                                                                lattice_vector_SL(3,3)**2)) 
      else
      periodic_cell_k = 0 
      endif  
 
      do i_lattice = -periodic_cell_i, periodic_cell_i,1          !$7
         do j_lattice = -periodic_cell_j, periodic_cell_j,1         !$6  
            do k_lattice = -periodic_cell_k, periodic_cell_k,1        !$5 
               if((abs(i_lattice).ne.0).or.(abs(j_lattice).ne.0).or.(abs(k_lattice).ne.0))then!$4
                do i_atom=1,n_atoms_SL,1 !$3
                 if(myid.eq.task_list_SL(i_atom)) then !$ tasks
                  ! LOOP GOES OVER UPPPER TRIANGLE OF HAMILTONIAN 
                  do j_atom=i_atom,n_atoms_SL,1 !$2

                      r_ij=0.0d0
                      TPP=0.0d0
                      dxyz=0.d0
                      coord_curr(:) = coords_SL(:,i_atom) + i_lattice*lattice_vector_SL(:,1) + &
                                                            j_lattice*lattice_vector_SL(:,2) + &
                                                            k_lattice*lattice_vector_SL(:,3)

                      dxyz(:)= coords_SL(:,j_atom)- coord_curr(:)
                      r_ij=dsqrt((dxyz(1))**2.0d0 +(dxyz(2))**2.0d0 + (dxyz(3))**2.0d0 )

                      if(r_ij.le.mbd_cfdm_dip_cutoff/bohr)  then
                        Rvdw_12=R_vdw_SL(i_atom)+R_vdw_SL(j_atom)
                        sigma=(r_ij/Rvdw_12)**beta
                        call MBD_TENSOR_MBD_rsSCS(dxyz,r_ij,Rvdw_12,beta,TPP) 
                        CFDM_prefactor=omega_cfdm_SL(i_atom)*omega_cfdm_SL(j_atom)*&
                                        sqrt(alpha_eff_SL(i_atom)*alpha_eff_SL(j_atom))
                        do i_index=1,3,1
                          do j_index=1,3,1
                          !FILL UPPER BLOCK
                          cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index)=&
                          cfdm_hamiltonian_periodic(3*i_atom-3+i_index,3*j_atom-3+j_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          !FILL LOWER BLOCK ALSO MAKE SURE THAT YOU DONT ADD FIELD
                          !CONTRIBUTION TWO TIMES IN DIAGONAL BLOCK below is the
                          !check
                          if(i_atom.NE.j_atom) then
                          cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index)=&
                          cfdm_hamiltonian_periodic(3*j_atom-3+j_index,3*i_atom-3+i_index)+ (TPP(i_index,j_index)*CFDM_prefactor)
                          endif
                         enddo
                        enddo
                      endif 

                  enddo !$2
                 endif !$tasks
                enddo !$3
               endif !$4
            enddo !$5
         enddo !$6
      enddo !$7
   call sync_tensors(cfdm_hamiltonian_periodic,3*n_atoms_SL)    
   cfdm_hamiltonian = cfdm_hamiltonian + cfdm_hamiltonian_periodic
endif

    errorflag=0
    LWORK=3*n_atoms_SL*(3+(3*n_atoms_SL)/2)
    call DSYEV('V','U',3*n_atoms_SL,cfdm_hamiltonian,3*n_atoms_SL,cfdm_eigenvalues,WORK,LWORK,errorflag)

    E_int=0.0d0
    E_nonint =0.0d0
    ene_mbd_rsSCS =   0.0
    NIMAG=0 
    if (errorflag.eq.0) then
             do i_atom =1,n_atoms_SL
             E_nonint = E_nonint + omega_cfdm_SL(i_atom)    
             enddo

             do i_atom =1,3*n_atoms_SL
              if(cfdm_eigenvalues(i_atom).ge.0.d0) then
                E_int = E_int + dsqrt(cfdm_eigenvalues(i_atom))
              else
                NIMAG= NIMAG +1  
              endif 
             enddo

         ene_mbd_rsSCS = ((0.5*E_int)-(1.5* E_nonint))/(SL_i*SL_j*SL_k)
         if(NIMAG.gt.0) then
         write(info_str,'(A,I4,A)')"***WARNING: found ",NIMAG," negative eigenvalues in MBD energy calculation."
         call output_string(info_str)  
         endif
    else
           write(info_str,*)"***Error:- Digonalization of CFDM hamiltonian failed"
           call output_string(info_str)
    endif
   !
      if(allocated(cfdm_hamiltonian))          deallocate(cfdm_hamiltonian)
      if(allocated(cfdm_hamiltonian_periodic)) deallocate(cfdm_hamiltonian_periodic) 
      if(allocated(cfdm_eigenvalues))          deallocate(cfdm_eigenvalues)
      if(allocated(coords_SL))                 deallocate(coords_SL)
      if(allocated(omega_cfdm_SL))             deallocate(omega_cfdm_SL)
      if(allocated(R_vdw_SL))                  deallocate(R_vdw_SL)
      if(allocated(alpha_eff_SL))              deallocate(alpha_eff_SL)  
  return
  endsubroutine get_mbd_rSCS



subroutine get_alpha_omega_and_Rp(HF,C6_free,alpha_free,w,a_w,Rpw)
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

       return
endsubroutine get_alpha_omega_and_Rp

subroutine contract_matrix(tensor)
        implicit none

        integer::ir,ic,i,j,i_row,i_col
        real*8,dimension(3,3)::tensor
        tensor(:,:)=0.0d0

        do ir=1,n_atoms,1
             do ic=1,n_atoms,1
                 i_row=0
                    do i=3*ir-2,3*ir,1
                       i_row=i_row+1
                       i_col=0
                          do j=3*ic-2,3*ic,1
                          i_col=i_col+1
                 tensor(i_row,i_col)=tensor(i_row,i_col) + relay_matrix(i,j)
                          enddo
                    enddo
             enddo
        enddo

        return
endsubroutine contract_matrix

      subroutine calculate_screened_polarizability(iso_polar_coupled)
        implicit none

        integer::i_row,i_col,i_index,j_index
        real*8,dimension(3,3)::matrix
        real*8,dimension(n_atoms),intent(OUT)::iso_polar_coupled
        real*8 :: WORK(9),eigen(3)
        integer ::   LWORK,errorflag
        iso_polar_coupled=0.d0

        ! o polar_coupled is formed by summing either rows 'blocks' or column
        ! 'blocks' of relay_matrix_screened and contains full anisotropic atom
        ! resolved polarizability 

        do i_row=1,n_atoms,1
        matrix=0.0
          do i_col=1,n_atoms,1
             do i_index=1,3,1
                do j_index=1,3,1
                   matrix(i_index,j_index) = matrix(i_index,j_index) + relay_matrix(3*i_row-3+i_index,3*i_col-3+j_index)
                enddo
             enddo
          enddo
       !o iso_polar_coupled contains average of trances of atom resolved
       !polarizabilities in polar_coupled
        LWORK=9  
        call DSYEV('V','U',3,matrix,3,eigen,WORK,LWORK,errorflag)  
        iso_polar_coupled(i_row) =(sum(eigen))/3.d0
        enddo
        return
      endsubroutine calculate_screened_polarizability



subroutine gauss_legendre_grid()
         casimir_omega(0) = 0.0000000 ;  casimir_omega_weight(0) = 0.0000000
         casimir_omega(1) = 0.0392901 ;  casimir_omega_weight(1) = 0.0786611
         casimir_omega(2) = 0.1183580 ;  casimir_omega_weight(2) = 0.0796400
         casimir_omega(3) = 0.1989120 ;  casimir_omega_weight(3) = 0.0816475
         casimir_omega(4) = 0.2820290 ;  casimir_omega_weight(4) = 0.0847872
         casimir_omega(5) = 0.3689190 ;  casimir_omega_weight(5) = 0.0892294
         casimir_omega(6) = 0.4610060 ;  casimir_omega_weight(6) = 0.0952317
         casimir_omega(7) = 0.5600270 ;  casimir_omega_weight(7) = 0.1031720
         casimir_omega(8) = 0.6681790 ;  casimir_omega_weight(8) = 0.1136050
         casimir_omega(9) = 0.7883360 ;  casimir_omega_weight(9) = 0.1273500
        casimir_omega(10) = 0.9243900 ; casimir_omega_weight(10) = 0.1456520
        casimir_omega(11) = 1.0817900 ; casimir_omega_weight(11) = 0.1704530
        casimir_omega(12) = 1.2684900 ; casimir_omega_weight(12) = 0.2049170
        casimir_omega(13) = 1.4966100 ; casimir_omega_weight(13) = 0.2544560
        casimir_omega(14) = 1.7856300 ; casimir_omega_weight(14) = 0.3289620
        casimir_omega(15) = 2.1691700 ; casimir_omega_weight(15) = 0.4480920
        casimir_omega(16) = 2.7106200 ; casimir_omega_weight(16) = 0.6556060
        casimir_omega(17) = 3.5457300 ; casimir_omega_weight(17) = 1.0659600
        casimir_omega(18) = 5.0273400 ; casimir_omega_weight(18) = 2.0635700
        casimir_omega(19) = 8.4489600 ; casimir_omega_weight(19) = 5.6851000
        casimir_omega(20) = 25.451700 ; casimir_omega_weight(20) = 50.955800
return
endsubroutine gauss_legendre_grid

! Free-atom C6, polarizability and vdW radius
subroutine get_vdw_param(atom,C6, alpha, R0)
implicit none

! local variables
character*2 :: atom
real*8 :: C6
real*8 :: alpha
real*8 :: R0

select case (atom)
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

end select

end subroutine get_vdw_param


subroutine output_string(info_str)
character*400 :: info_str  
if(myid.eq.0) write(*,*)trim(info_str)
endsubroutine output_string


!-------------------------------------------------------------------------
END MODULE UTILS

