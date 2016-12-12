  subroutine td_write_floquet(out_floquet, hm, gr, st, ks, iter)
    type(c_ptr),       intent(inout)   :: out_floquet
    type(hamiltonian_t), intent(inout) :: hm
    type(grid_t),      intent(inout)   :: gr
    type(states_t),    intent(inout)   :: st !< at iter=0 this is the
    groundstate                                                                                                                                       
    type(v_ks_t),      intent(in)      :: ks
    integer,           intent(in)      :: iter 

    CMPLX, allocatable :: hmss(:,:), psi(:,:,:), hpsi(:,:,:),
    temp_state1(:,:), temp_state2(:,:)
    CMPLX, allocatable :: HFloquet(:,:,:), HFloq_eff(:,:), temp(:,:)
    FLOAT, allocatable :: eigenval(:), bands(:,:)
    character(len=80) :: filename
    integer :: it, nT, ik, ist, jst, in, im, inm, file, idim, nik, ik_count
    integer :: Forder, Fdim, m0, n0, n1, nst, ii, jj, lim_nst
    logical :: downfolding = .false.
    type(mesh_t) :: mesh
    type(states_t) :: hm_st

    FLOAT :: dt, Tcycle, omega

    PUSH_SUB(td_write_floquet)

    ! this does not depend on propagation, so we do it only once                                                                                                                                                        
    if(.not. iter == 0) then
       POP_SUB(td_write_floquet)
       return
    end if

    mesh = gr%der%mesh
    nst = st%nst

    !for now no domain distributionallowed                                                                                                                                                                              
    ASSERT(mesh%np == mesh%np_global)

   ! this is used to initialize the hpsi (more effiecient ways?)                                                                                                                                                        
    call states_copy(hm_st, st)

    !%Variable TDFloquetFrequency                                                                                                                                                                                       
    !%Type float                                                                                                                                                                                                        
    !%Default 0                                                                                                                                                                                                         
    !%Section Time-Dependent::TD Output                                                                                                                                                                                 
    !%Description                                                                                                                                                                                                       
    !% Frequency for the Floquet analysis, this should be the carrier
    !%frequency or integer multiples of it.                                                                                                             
    !% Other options will work, but likely be nonsense.                                                                                                                                                                 
    !%                                                                                                                                                                                                                  
    !%End                                                                                                                                                                                                               
    call parse_variable('TDFloquetFrequency', M_ZERO, omega, units_inp%energy)
    call messages_print_var_value(stdout,'Frequency used for Floquet
    !%analysis', omega)
    if(omega==M_ZERO) then
       message(1) = "Please give a non-zero value for TDFloquetFrequency"
       call messages_fatal(1)
    endif

    ! get time of one cycle                                                                                                                                                                                             
    Tcycle=M_TWO*M_PI/omega

    !%Variable TDFloquetSample                                                                                                                                                                                          
    !%Type integer                                                                                                                                                                                                      
    !%Default 20                                                                                                                                                                                                        
    !%Section Time-Dependent::TD Output                                                                                                                                                                                 
    !%Description                                                                                                                                                                                                       
    !% Number of points on which one Floquet cycle is sampled in the
    !%time-integral of the Floquet analysis.                                                                                                             
    !%                                                                                                                                                                                                                  
    !%End                                                                                                                                                                                                               
    call parse_variable('TDFloquetSample',20 ,nt)
    call messages_print_var_value(stdout,'Number of Floquet time-sampling
    !%points', nT)
    dt = Tcycle/real(nT)

    !%Variable TDFloquetDimension                                                                                                                                                                                       
    !%Type integer                                                                                                                                                                                                      
    !%Default -1                                                                                                                                                                                                        
    !%Section Time-Dependent::TD Output                                                                                                                                                                                 
    !%Description                                                                                                                                                                                                       
    !% Order of Floquet Hamiltonian. If negative number is given, downfolding
    !%is performed.                                                                                                                             
    !%End                                                                                                                                                                                                               
    call parse_variable('TDFloquetDimension',-1,Forder)
    if(Forder.ge.0) then
       call messages_print_var_value(stdout,'Order of multiphoton
    !%Floquet-Hamiltonian', Forder)
       !Dimension of multiphoton Floquet-Hamiltonian                                                                                                                                                                    
       Fdim = 2*Forder+1
    else
       message(1) = 'Floquet-Hamiltonian is downfolded'
       call messages_info(1)
       downfolding = .true.
       Forder = 1
       Fdim = 3
    endif

    dt = Tcycle/real(nT)


   POP_SUB(td_write_floquet)

  end subroutine td_write_floquet

