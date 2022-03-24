program HartreeFock

   ! Demonstration program that can be used as a starting point
   ! Lucas Visscher, March 2022

   use molecular_structure
   use ao_basis
   use compute_integrals
   use diagonalization

     implicit none

     ! Variable containing the molecular structure
     type(molecular_structure_t) :: molecule
     ! Variable containing the atomic orbital basis
     type(basis_set_info_t) :: ao_basis

     ! Variable naming as in the description of the exercise
     integer  :: n_AO, n_occ, it
     integer  :: kappa, lambda
     real(8)  :: E_HF_1, E_HF_2, treshold
     real(8), allocatable :: F(:,:),V(:,:),T(:,:),S(:,:), C(:,:), eps(:), D(:,:), H_core(:,:)

     ! The following large array can be eliminated when Fock matrix contruction is implemented
     real(8), allocatable :: ao_integrals (:,:,:,:)

     treshold = 1e-9
     E_HF_1   = 1
     E_HF_2   = 0
     it       = 0
   
     ! Definition of the molecule
     call define_molecule(molecule)

     ! Definition of the GTOs
     call define_basis(ao_basis)
     n_AO = ao_basis%nao
   
     ! Definition of the number of occupied orbitals
     n_occ = 3 ! hardwired for this demonstration program, should be set via input

     ! Compute the overlap matrix
     allocate (S(n_AO,n_AO))
     call   compute_1e_integrals ("OVL",ao_basis,ao_basis,S)

     ! Compute the kinetic matrix
     allocate (T(n_AO,n_AO))
     call   compute_1e_integrals ("KIN",ao_basis,ao_basis,T)

     ! Compute the potential matrix
     allocate (V(n_AO,n_AO))
     call   compute_1e_integrals ("POT",ao_basis,ao_basis,V,molecule)
     
     ! Calculating the atom orbital integrals
     allocate (ao_integrals(n_AO,n_AO,n_AO,n_AO))
     call generate_2int (ao_basis,ao_integrals)

     ! Compute the core Hamiltonian matrix (the potential is positive, we scale with -e = -1 to get to the potential energy matrix)
     allocate (F(n_AO,n_AO))
     allocate (H_core(n_AO,n_AO))
     H_core = T - V
     print *, 'Core Hamiltonian matrix:'
     write (*,'(3f15.5)') H_core
     allocate (D(n_AO,n_AO))

     ! Diagonalize the core hamiltonian, creating a dummy set of coefficients
     allocate (C(n_AO,n_AO))
     allocate (eps(n_AO))
     call solve_genev (F,S,C,eps)
     print*, "Orbital energies for the core Hamiltonian:",eps

     
      do
         if (abs(E_HF_2 - E_HF_1) < treshold) then
           exit
         else if (it > 100) then
         exit
         else
            it = it +1
           ! Form the density matrix
     
            do lambda = 1, n_ao
             do kappa = 1, n_ao
               D(kappa,lambda) = sum(C(kappa,1:n_occ)*C(lambda,1:n_occ))
             end do
            end do

           ! Saving the Hartree Fock energy of the previous iteration
            E_HF_1 = E_HF_2

           !Create the Fock matrix with hcore and 2-electron integrals
            do lambda = 1, n_AO
              do kappa = 1, n_AO
                F = H_core + 2.D0 * ao_integrals(:,:,kappa,lambda) * D
                F = F - 1.D0 * ao_integrals(:,lambda,kappa,:)* D
              end do
            end do
           ! Compute the Hartree-Fock energy
                E_HF_2 = sum((H_core+F)*D)
           ! Diagonalize the Fock matrix
           call solve_genev (F,S,C,eps)
          end if
       end do 
     
     print*, "The Hartree-Fock energy:    ", E_HF_2
     print*, 'The energy was found after', it, 'iterations'  
   end
   subroutine define_molecule(molecule)
     ! This routine should be improved such that an arbitrary molecule can be given as input
     ! the coordinates below are for a be-he dimer oriented along the x-axis with a bond length of 2 au
     use molecular_structure
     type(molecular_structure_t), intent(inout) :: molecule
     real(8) :: charge(2),coord(3,2)
     charge(1)   = 4.D0
     charge(2)   = 2.D0
     coord       = 0.D0
     coord(1,2)  = 2.D0
     call add_atoms_to_molecule(molecule,charge,coord)
   end subroutine

   subroutine define_basis(ao_basis)
    ! This routine can be extended to use better basis sets 
    ! The coordinates of the shell centers are the nuclear coordinates
    ! Think of a refactoring of define_molecule and define_basis to ensure consistency 
     use ao_basis
     type(basis_set_info_t), intent(inout) :: ao_basis
     type(basis_func_info_t) :: gto
     ! Be:  2 uncontracted s-funs:    l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),4.D0)
     call add_shell_to_basis(ao_basis,0,(/0.D0,0.D0,0.D0/),1.D0)
     ! He:  1 uncontracted s-fun:     l      coord          exp      
     call add_shell_to_basis(ao_basis,0,(/2.D0,0.D0,0.D0/),1.D0)
   end subroutine

   
