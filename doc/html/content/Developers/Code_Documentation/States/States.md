---
Title: States
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

Wavefunctions in Octopus
========================

The wave functions in Octopus are referred to as the states.

They are handled by the module `states_oct_m`, which is defined in [states/states.F90](https://gitlab.com/octopus-code/octopus/blob/develop/src/states/states.F90).

The exact way how states are stored in memory is flexible and depends on optimization  (and accelerator) settings in the input file. In particular, the order of inices depends on the `PACKED` setting.

States are stored in a hierarchy of 'containers'. Concepts of this hierarchy 
include _groups_, _batches_ and _blocks_.



The top level data structure, describing states is:

```Fortran
#include_type_def states_abst_t
```

This structure contains mainly metadata about states, describing _how_ states are represented, as well as pointers to other quantities, which are common to all states, such as the desnity, the current, etc.

The dimensions object contains a number of variables. The most relevant for this discussion is the `dim` variable, which denotes the dimension of one state, being `1` for spin-less states and `2` for spinors.

```Fortran
#include_type_def states_dim_t
```


The `states_t` object is initialized in `states_init()`, which is called from `system_init()`.

The wave funtions themselves are stored in
```Fortran
    type(states_group_t)     :: group
```
which, in turn, is defined in the module `states_group_oct_m` in [`src/states/states_group.F90`](https://gitlab.com/octopus-code/octopus/blob/develop/src/states/states_group.F90):

```Fortran
#include_type_def states_group_t
```

The `group` contains all wave functions, grouped together in blocks or batches.
They are organised in an array of `batch_t` structures.
```Fortran
    type(batch_t), pointer   :: psib(:, :)            !< A set of wave-functions
```
The indexing is as follows: `psib(ib,iqb)` where `ib` is the block index, and `iqn` the **k**-point. See below for the routine `states_init_block(st, mesh, verbose)` which creates the `group` object. On a given node, only wave functions of local blocks are available.
The `group` object does contain all information on how the batches are distributed over nodes.






Creating the wave functions
---------------------------
A number of steps in initializing the `states_t` object are called from the `system_init()` routine:

`states_init()`:  
parses states-related input variables, and allocates memory for some book keeping variables. It does not allocate any memory for the states themselves.

`states_distribute_nodes()`:  
...


`states_density_init()`:  
allocates memory for the density (`rho`) and the core density (`rho_core`).

`states_exec_init()`:  
1. Fills in the block size (`st\%d\%block_size`);
2. Finds out whether or not to pack the states (`st\%d\%pack_states`);
3. Finds out the orthogonalization method (`st\%d\%orth_method`).
  


Memory for the actual wave functions is allocated in `states_allocate_wfns()` which is called from the corresponding `*_run()` routines, such as `scf_run()` or `td_run()`, etc.



```Fortran
  ! ---------------------------------------------------------
  !> Allocates the KS wavefunctions defined within a states_t structure.
  subroutine states_allocate_wfns(st, mesh, wfs_type, alloc_Left)
    type(states_t),         intent(inout)   :: st
    type(mesh_t),           intent(in)      :: mesh
    type(type_t), optional, intent(in)      :: wfs_type
    logical,      optional, intent(in)      :: alloc_Left !< allocate an additional set of wfs to store left eigenstates

    PUSH_SUB(states_allocate_wfns)

    if (present(wfs_type)) then
      ASSERT(wfs_type == TYPE_FLOAT .or. wfs_type == TYPE_CMPLX)
      st%priv%wfs_type = wfs_type
    end if

    call states_init_block(st, mesh)  !< Allocate memory for blocks and states
                                      !! Determine the block structure and set up index arrays
    call states_set_zero(st)          !< Initialize the states to zero.

    POP_SUB(states_allocate_wfns)
  end subroutine states_allocate_wfns
```


The routine `states_init_block` initializes the data components in `st` that describe how the states are distributed in blocks:

`st%nblocks`: this is the number of blocks in which the states are divided. 
Note that this number is the total number of blocks,
regardless of how many are actually stored in each node.  
`block_start`: in each node, the index of the first block.  
`block_end`: in each node, the index of the last block.
If the states are not parallelized, then `block_start` is 1 and `block_end` is `st%nblocks`.  
`st%iblock(1:st%nst, 1:st%d%nik)`: it points, for each state, to the block that contains it.  
`st%block_is_local()`: `st%block_is_local(ib)` is `.true.` if block `ib` is stored in the running node.  
`st%block_range(1:st%nblocks, 1:2)`: Block ib contains states fromn st\%block_range(ib, 1) to st\%block_range(ib, 2)  
`st%block_size(1:st%nblocks)`: Block ib contains a number st\%block_size(ib) of states.  
`st%block_initialized`: it should be .false. on entry, and .true. after exiting this routine.  

The set of batches `st%psib(1:st%nblocks)` contains the `block`s themselves.
  
```Fortran  
  subroutine states_init_block(st, mesh, verbose)
    type(states_t),           intent(inout) :: st
    type(mesh_t),             intent(in)    :: mesh
    logical, optional,        intent(in)    :: verbose

    integer :: ib, iqn, ist
    logical :: same_node, verbose_
    integer, allocatable :: bstart(:), bend(:)

    PUSH_SUB(states_init_block)

    ! Allocate memory for block indices:

    SAFE_ALLOCATE(bstart(1:st%nst))
    SAFE_ALLOCATE(bend(1:st%nst))
    SAFE_ALLOCATE(st%group%iblock(1:st%nst, 1:st%d%nik))

    st%group%iblock = 0

    verbose_ = optional_default(verbose, .true.)

    ! count and assign blocks
    ib = 0
    st%group%nblocks = 0
    bstart(1) = 1 

    do ist = 1, st%nst     ! Loop over all states:
      INCR(ib, 1)

      st%group%iblock(ist, st%d%kpt%start:st%d%kpt%end) = st%group%nblocks + 1

      same_node = .true.
      if(st%parallel_in_states .and. ist /= st%nst) then
        ! We have to avoid that states that are in different nodes end
        ! up in the same block
        same_node = (st%node(ist + 1) == st%node(ist))
      end if

      if(ib == st%d%block_size .or. ist == st%nst .or. .not. same_node) then
        ib = 0
        INCR(st%group%nblocks, 1)
        bend(st%group%nblocks) = ist
        if(ist /= st%nst) bstart(st%group%nblocks + 1) = ist + 1
      end if

    end do ! ist

    ! Allocate memory for the blocks:

    SAFE_ALLOCATE(st%group%psib(1:st%group%nblocks, st%d%kpt%start:st%d%kpt%end))

    ! Mark the local blocks:

    SAFE_ALLOCATE(st%group%block_is_local(1:st%group%nblocks, st%d%kpt%start:st%d%kpt%end))
    st%group%block_is_local = .false.
    st%group%block_start  = -1
    st%group%block_end    = -2  ! this will make that loops block_start:block_end do not run if not initialized

    ! Allocate memory for the states in the local blocks:

    do ib = 1, st%group%nblocks         ! loop over blocks
      if(bstart(ib) >= st%st_start .and. bend(ib) <= st%st_end) then
        if(st%group%block_start == -1) st%group%block_start = ib
        st%group%block_end = ib
        do iqn = st%d%kpt%start, st%d%kpt%end
          st%group%block_is_local(ib, iqn) = .true.

          if (states_are_real(st)) then
            call batch_init(st%group%psib(ib, iqn), st%d%dim, bend(ib) - bstart(ib) + 1)
            call dbatch_allocate(st%group%psib(ib, iqn), bstart(ib), bend(ib), mesh%np_part)
          else
            call batch_init(st%group%psib(ib, iqn), st%d%dim, bend(ib) - bstart(ib) + 1)
            call zbatch_allocate(st%group%psib(ib, iqn), bstart(ib), bend(ib), mesh%np_part)
          end if
          ! The above batch_init calls resolve to batch_init_empty(this, dim, nst)

        end do
      end if
    end do

    ! copy block indices to `group` data structure:

    SAFE_ALLOCATE(st%group%block_range(1:st%group%nblocks, 1:2))
    SAFE_ALLOCATE(st%group%block_size(1:st%group%nblocks))
    
    st%group%block_range(1:st%group%nblocks, 1) = bstart(1:st%group%nblocks)
    st%group%block_range(1:st%group%nblocks, 2) = bend(1:st%group%nblocks)
    st%group%block_size(1:st%group%nblocks) = bend(1:st%group%nblocks) - bstart(1:st%group%nblocks) + 1

    st%group%block_initialized = .true.

    SAFE_ALLOCATE(st%group%block_node(1:st%group%nblocks))

    ASSERT(associated(st%node))
    ASSERT(all(st%node >= 0) .and. all(st%node < st%mpi_grp%size))
    
    do ib = 1, st%group%nblocks
      st%group%block_node(ib) = st%node(st%group%block_range(ib, 1))
      ASSERT(st%group%block_node(ib) == st%node(st%group%block_range(ib, 2)))
    end do
        
    SAFE_DEALLOCATE_A(bstart)
    SAFE_DEALLOCATE_A(bend)
    POP_SUB(states_init_block)
  end subroutine states_init_block
```

The allocation of memory for the actual wave functions is performed in `batch_init_empty()` and `X(batch_allocate)()`. This routine, and the related `X(batch_add_state)()` show most clearly how the different memory blocks are related.

[`batch_init_empty()`](https://gitlab.com/octopus-code/octopus/blob/develop/src/grid/batch.F90) allocates the memory for `batch_state_t` `states` and `batch_states_l_t` `states_linear` and _nullifies_ the pointers within this types. Note that no memory for the actual wave functions has been allocated yet.

```Fortran
  subroutine batch_init_empty (this, dim, nst)
    type(batch_t), intent(out)   :: this
    integer,       intent(in)    :: dim    !< dimension of the state (see st%d%dim), either 1 or 2.
    integer,       intent(in)    :: nst    !< number of states in the batch
    
    integer :: ist

    PUSH_SUB(batch_init_empty)

    this%is_allocated = .false.
    this%nst = nst
    this%dim = dim
    this%current = 1
    nullify(this%dpsicont, this%zpsicont, this%spsicont, this%cpsicont)
    
    SAFE_ALLOCATE(this%states(1:nst))
    do ist = 1, nst
      nullify(this%states(ist)%dpsi)
      nullify(this%states(ist)%zpsi)
      nullify(this%states(ist)%spsi)
      nullify(this%states(ist)%cpsi)
    end do
    
    this%nst_linear = nst*dim
    SAFE_ALLOCATE(this%states_linear(1:this%nst_linear))
    do ist = 1, this%nst_linear
      nullify(this%states_linear(ist)%dpsi)
      nullify(this%states_linear(ist)%zpsi)
      nullify(this%states_linear(ist)%spsi)
      nullify(this%states_linear(ist)%cpsi)
    end do
    
    this%in_buffer_count = 0
    this%status = BATCH_NOT_PACKED

    this%ndims = 2
    SAFE_ALLOCATE(this%ist_idim_index(1:this%nst_linear, 1:this%ndims))

    POP_SUB(batch_init_empty)
  end subroutine batch_init_empty
```


Now, memory for a block of states is allocated in the continuous storage block `X(psicont)`. 
This is done using the `X(batch_allocate)`m defined in [`src/grid/batch_inc.F90`](https://gitlab.com/octopus-code/octopus/blob/develop/src/grid/batch_inc.F90) routine.


```Fortran
subroutine X(batch_allocate)(this, st_start, st_end, np, fill_zeros)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: st_start
  integer,        intent(in)    :: st_end
  integer,        intent(in)    :: np
  logical, optional, intent(in) :: fill_zeros

  integer :: ist

  PUSH_SUB(X(batch_allocate))


  !> allocate memory in the rank-3 object `X(psicont)(:,:,:)`:
  SAFE_ALLOCATE(this%X(psicont)(1:np, 1:this%dim, 1:st_end - st_start + 1))

  !> if requested, initialize to zero:
  if (optional_default(fill_zeros, .true.)) this%X(psicont) = R_TOTYPE(M_ZERO)

  this%is_allocated = .true.

  !> add the states to "states" and "states_linear" (one by one)
  do ist = st_start, st_end
    call X(batch_add_state)(this, ist, this%X(psicont)(:, :, ist - st_start + 1))
  end do

  POP_SUB(X(batch_allocate))
end subroutine X(batch_allocate)
```
Sub-blocks of this continuous memory (representing one state) are then added to `states` and `states_linear`, using:

```Fortran
subroutine X(batch_add_state)(this, ist, psi)
  type(batch_t),  intent(inout) :: this
  integer,        intent(in)    :: ist
  R_TYPE, target, intent(in)    :: psi(:, :)

  integer :: idim, ii

  PUSH_SUB(X(batch_add_state))

  ASSERT(this%current <= this%nst)

  ! populate the 'states' array:

  this%states(this%current)%ist    =  ist
  this%states(this%current)%X(psi) => psi

  ! now we also populate the linear array
  do idim = 1, this%dim
    ii = this%dim*(this%current - 1) + idim
    this%states_linear(ii)%X(psi) => psi(:, idim)
    this%ist_idim_index(ii, 1) = ist
    this%ist_idim_index(ii, 2) = idim
  end do

  this%current = this%current + 1

  POP_SUB(X(batch_add_state))
end subroutine X(batch_add_state)
```


```Fortran
subroutine X(batch_add_state_linear)(this, psi)
  type(batch_t),  intent(inout) :: this
  R_TYPE, target, intent(in)    :: psi(:)

  PUSH_SUB(X(batch_add_state_linear))

  ASSERT(this%current <= this%nst_linear)
  this%states_linear(this%current)%X(psi) => psi
  this%ist_idim_index(this%current, 1) = this%current

  this%current = this%current + 1

  POP_SUB(X(batch_add_state_linear))
end subroutine X(batch_add_state_linear)
```

Questions:
----------

How are the different objects pointing to states related?

* The usual storage for states is in `states_linear` which can be shadowed by `pack` in case 
  packed states are used.

What is the difference between `batch_add_state` and `batch_add_state_linear`?
