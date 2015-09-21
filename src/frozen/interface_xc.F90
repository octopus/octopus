#include "global.h"

module interface_xc_m

  use derivatives_m
  use global_m
  use grid_m
  use kinds_m
  use mesh_m
  use messages_m
  use profiling_m
  use XC_F90(lib_m)

  private

  public ::     &
    XC_UNKNOWN, &
    XC_NONE

  public ::         &
    interface_xc_t

  public ::                      &
    interface_xc_init,           &
    interface_xc_start,          &
    interface_xc_get,            &
    interface_xc_lda_exc,        &
    interface_xc_lda_vxc,        &
    interface_xc_lda_exc_vxc,    &
    interface_xc_gga_exc,        &
    interface_xc_gga_vxc,        &
    interface_xc_gga_exc_vxc,    &
    interface_xc_copy,           &
    interface_xc_end

  integer, parameter :: XC_UNKNOWN =-1
  integer, parameter :: XC_NONE    = 0

  integer, parameter :: BLOCK_SIZE = 1000

  integer, dimension(3), parameter :: xc_flags_nd=(/XC_FLAGS_1D, XC_FLAGS_2D, XC_FLAGS_3D/)

  type :: interface_xc_t
    private
    type(mesh_t),        pointer :: mesh   =>null()
    type(derivatives_t), pointer :: der    =>null()
    integer                      :: nblock = 0
    integer                      :: ndim   = 0
    integer                      :: nspin  = 0
    integer                      :: id     = XC_NONE ! identifier
    integer                      :: family = XC_NONE ! LDA, GGA, etc.
    integer                      :: spin   = XC_NONE ! XC_UNPOLARIZED | XC_POLARIZED
    integer                      :: kind   = XC_NONE ! exchange, correlation, or exchange-correlation
    integer                      :: flags  = XC_NONE ! XC_FLAGS_HAVE_EXC + XC_FLAGS_HAVE_VXC + ...
    type(XC_F90(pointer_t))      :: conf             ! the pointer used to call the library
    type(XC_F90(pointer_t))      :: info             ! information about the functional
  end type interface_xc_t

  interface interface_xc_get
    module procedure interface_xc_get_info
    module procedure interface_xc_get_mesh
    module procedure interface_xc_get_derivatives
  end interface interface_xc_get

contains

  ! ---------------------------------------------------------
  subroutine interface_xc_init(this, id, nspin, blocksize)
    type(interface_xc_t), intent(out) :: this
    integer,              intent(in)  :: id
    integer,              intent(in)  :: nspin
    integer,    optional, intent(in)  :: blocksize

    PUSH_SUB(interface_xc_init)

    this%id = id
    ASSERT(this%id>XC_UNKNOWN)
    if(this%id>XC_NONE)then
      this%family = XC_F90(family_from_id)(this%id)
      ASSERT(this%family>XC_NONE)
      this%nblock = BLOCK_SIZE
      if(present(blocksize)) this%nblock = blocksize
      ASSERT(this%nblock>0)
      this%nspin = nspin
      ASSERT(this%nspin>0)
      ASSERT(this%nspin<3)
      this%spin = XC_UNPOLARIZED
      if(this%nspin>1) this%spin = XC_POLARIZED
      call XC_F90(func_init)(this%conf, this%info, this%id, this%spin)
      this%kind = XC_F90(info_kind)(this%info)
      this%flags = XC_F90(info_flags)(this%info)
    end if

    POP_SUB(interface_xc_init)
  end subroutine interface_xc_init

  ! ---------------------------------------------------------
  subroutine interface_xc_start(this, mesh, der)
    type(interface_xc_t),        intent(inout) :: this
    type(mesh_t),        target, intent(in)    :: mesh
    type(derivatives_t), target, intent(in)    :: der

    PUSH_SUB(interface_xc_start)

    ASSERT(.not.associated(this%mesh))
    ASSERT(.not.associated(this%der))
    this%mesh => mesh
    this%der => der
    this%ndim = this%mesh%sb%dim
    if(this%id>XC_NONE)then
      if(iand(this%flags,xc_flags_nd(this%ndim))==0)then
        write(unit=message(1), fmt="(a,i1.1,a2)") "Cannot use the specified functional in ", this%ndim, "D."
        call messages_fatal(1)
      end if
    end if

    POP_SUB(interface_xc_start)
  end subroutine interface_xc_start

  ! ---------------------------------------------------------
  subroutine interface_xc_get_info(this, id, family, kind, nspin, polarized, ndim)
    type(interface_xc_t), intent(in)  :: this
    integer,    optional, intent(out) :: id
    integer,    optional, intent(out) :: family
    integer,    optional, intent(out) :: kind
    integer,    optional, intent(out) :: nspin
    logical,    optional, intent(out) :: polarized
    integer,    optional, intent(out) :: ndim

    PUSH_SUB(interface_xc_get_info)

    if(present(id)) id = this%id
    if(present(family)) family = this%family
    if(present(kind)) kind = this%kind
    if(present(nspin)) nspin = this%nspin
    if(present(polarized)) polarized = (this%spin==XC_POLARIZED)
    if(present(ndim)) ndim = this%ndim

    POP_SUB(interface_xc_get_info)
  end subroutine interface_xc_get_info

  ! ---------------------------------------------------------
  subroutine interface_xc_get_mesh(this, that)
    type(interface_xc_t), target, intent(in) :: this
    type(mesh_t),        pointer             :: that

    PUSH_SUB(interface_xc_get_mesh)

    nullify(that)
    if(associated(this%mesh)) that => this%mesh

    POP_SUB(interface_xc_get_mesh)
  end subroutine interface_xc_get_mesh

  ! ---------------------------------------------------------
  subroutine interface_xc_get_derivatives(this, that)
    type(interface_xc_t), target, intent(in) :: this
    type(derivatives_t), pointer             :: that

    PUSH_SUB(interface_xc_get_derivatives)

    nullify(that)
    if(associated(this%der)) that => this%der

    POP_SUB(interface_xc_get_derivatives)
  end subroutine interface_xc_get_derivatives

  ! ---------------------------------------------------------
  elemental subroutine interface_xc_density(irho, orho)
    real(kind=wp), intent(in)  :: irho
    real(kind=wp), intent(out) :: orho

    if(irho<0.0_wp)then
      orho = 0.0_wp
    else
      orho = irho
    end if

  end subroutine interface_xc_density

  ! ---------------------------------------------------------
  pure subroutine interface_xc_block_density(nb, irho, orho)
    integer,                       intent(in)  :: nb
    real(kind=wp), dimension(:,:), intent(in)  :: irho
    real(kind=wp), dimension(:,:), intent(out) :: orho

    integer :: ib

    orho = 0.0_wp
    do ib = 1, nb
      call interface_xc_density(irho(ib,:),orho(:,ib))
    end do

  end subroutine interface_xc_block_density

  ! ---------------------------------------------------------
  pure subroutine interface_xc_block_gradient(nb, ns, gradient, sigma)
    integer,                         intent(in)  :: nb
    integer,                         intent(in)  :: ns
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    real(kind=wp), dimension(:,:),   intent(out) :: sigma

    integer :: ip

    sigma = 0.0_wp
    do ip = 1, nb
      sigma(1,ip) = dot_product(gradient(ip,:,1),gradient(ip,:,1))
      if(ns==3)then
        sigma(2,ip) = dot_product(gradient(ip,:,1),gradient(ip,:,2))
        sigma(3,ip) = dot_product(gradient(ip,:,2),gradient(ip,:,2))
      end if
    end do

  end subroutine interface_xc_block_gradient

  ! ---------------------------------------------------------
  pure subroutine interface_xc_deblock_vxc(nb, ipot, opot)
    integer,                       intent(in)  :: nb
    real(kind=wp), dimension(:,:), intent(in)  :: ipot
    real(kind=wp), dimension(:,:), intent(out) :: opot

    integer :: ip

    do ip = 1, nb
      opot(ip,:) = ipot(:,ip)
    end do

  end subroutine interface_xc_deblock_vxc
  
  ! ---------------------------------------------------------
  pure subroutine interface_xc_deblock_dedgd(nb, ns, vsigma, grad, dedgd)
    integer,                         intent(in)  :: nb
    integer,                         intent(in)  :: ns
    real(kind=wp), dimension(:,:),   intent(in)  :: vsigma
    real(kind=wp), dimension(:,:,:), intent(in)  :: grad
    real(kind=wp), dimension(:,:,:), intent(out) :: dedgd

    integer :: ip

    do ip = 1, nb
      dedgd(ip,:,1) = 2.0_wp * vsigma(1,ip) * grad(ip,:,1)
      if(ns==3)then
        dedgd(ip,:,1) = dedgd(ip,:,1) + vsigma(2,ip) * grad(ip,:,2)
        dedgd(ip,:,2) = 2.0_wp * vsigma(3,ip) * grad(ip,:,2) + vsigma(2,ip) * grad(ip,:,1)
      end if
    end do

  end subroutine interface_xc_deblock_dedgd

  ! ---------------------------------------------------------
  subroutine interface_xc_lda_exc(this, density, exc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc

    real(kind=wp), dimension(this%spin,this%nblock) :: l_dens
    real(kind=wp), dimension(          this%nblock) :: l_zk
    integer                                         :: ip, nb

    PUSH_SUB(interface_xc_lda_exc)

    exc = 0.0_wp
    if((this%id>XC_NONE).and.(iand(this%flags, XC_FLAGS_HAVE_EXC)/=0))then
      do ip = 1, this%mesh%np, this%nblock
        nb = min(this%mesh%np-ip+1, this%nblock)
        call interface_xc_block_density(nb, density(ip:,:), l_dens)
        call XC_F90(lda_exc)(this%conf, nb, l_dens(1,1), l_zk(1))
        exc(ip:ip+nb-1) = l_zk(1:nb)
      end do
    end if

    POP_SUB(interface_xc_lda_exc)
  end subroutine interface_xc_lda_exc

  ! ---------------------------------------------------------
  subroutine interface_xc_lda_vxc(this, density, vxc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: vxc

    real(kind=wp), dimension(this%spin,this%nblock) :: l_dens, l_dedd
    integer                                         :: ip, nb

    PUSH_SUB(interface_xc_lda_vxc)

    vxc = 0.0_wp
    if((this%id>XC_NONE).and.(iand(this%flags, XC_FLAGS_HAVE_VXC)/=0))then
      do ip = 1, this%mesh%np, this%nblock
        nb = min(this%mesh%np-ip+1, this%nblock)
        call interface_xc_block_density(nb, density(ip:,:), l_dens)
        call XC_F90(lda_vxc)(this%conf, nb, l_dens(1,1), l_dedd(1,1))
        call interface_xc_deblock_vxc(nb, l_dedd, vxc(ip:,:))
      end do
    end if

    POP_SUB(interface_xc_lda_vxc)
  end subroutine interface_xc_lda_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_lda_exc_vxc(this, density, exc, vxc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc

    real(kind=wp), dimension(this%spin,this%nblock) :: l_dens, l_dedd
    real(kind=wp), dimension(          this%nblock) :: l_zk
    integer                                         :: ip, nb

    PUSH_SUB(interface_xc_lda_exc_vxc)

    exc = 0.0_wp
    vxc = 0.0_wp
    if((this%id>XC_NONE)                           &
      .and.(iand(this%flags,XC_FLAGS_HAVE_EXC)/=0) &
      .and.(iand(this%flags,XC_FLAGS_HAVE_VXC)/=0))then
      do ip = 1, this%mesh%np, this%nblock
        nb = min(this%mesh%np-ip+1, this%nblock)
        call interface_xc_block_density(nb, density(ip:,:), l_dens)
        call XC_F90(lda_exc_vxc)(this%conf, nb, l_dens(1,1), l_zk(1), l_dedd(1,1))
        call interface_xc_deblock_vxc(nb, l_dedd, vxc(ip:,:))
        exc(ip:ip+nb-1) = l_zk(1:nb)
      end do
    end if

    POP_SUB(interface_xc_lda_exc_vxc)
  end subroutine interface_xc_lda_exc_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_block_gga_exc(this, density, gradient, exc)
    type(interface_xc_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:),   intent(in)  :: density
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    real(kind=wp), dimension(:),     intent(out) :: exc

    real(kind=wp), dimension(this%spin,    this%nblock) :: l_dens
    real(kind=wp), dimension(2*this%spin-1,this%nblock) :: l_sigma
    real(kind=wp), dimension(              this%nblock) :: l_zk
    integer                                             :: ip, ns, nb

    PUSH_SUB(interface_xc_block_gga_exc)

    exc = 0.0_wp
    if((this%id>XC_NONE).and.(iand(this%flags, XC_FLAGS_HAVE_EXC)/=0))then
      ns = 2 * this%spin - 1
      do ip = 1, this%mesh%np, this%nblock
        nb = min(this%mesh%np-ip+1, this%nblock)
        call interface_xc_block_density(nb, density(ip:,:), l_dens)
        call interface_xc_block_gradient(nb, ns, gradient(ip:,:,:), l_sigma)
        call XC_F90(gga_exc)(this%conf, nb, l_dens(1,1), l_sigma(1,1), l_zk(1))
        exc(ip:ip+nb-1) = l_zk(1:nb)
      end do
    end if

    POP_SUB(interface_xc_block_gga_exc)
  end subroutine interface_xc_block_gga_exc

  ! ---------------------------------------------------------
  subroutine interface_xc_block_gga_vxc(this, density, gradient, vxc, dedg)
    type(interface_xc_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:),   intent(in)  :: density
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    real(kind=wp), dimension(:,:),   intent(out) :: vxc
    real(kind=wp), dimension(:,:,:), intent(out) :: dedg

    real(kind=wp), dimension(this%spin,    this%nblock) :: l_dens, l_dedd
    real(kind=wp), dimension(2*this%spin-1,this%nblock) :: l_sigma
    real(kind=wp), dimension(2*this%spin-1,this%nblock) :: l_vsigma
    integer                                             :: ip, ns, nb

    PUSH_SUB(interface_xc_block_gga_vxc)

    vxc = 0.0_wp
    if((this%id>XC_NONE).and.(iand(this%flags, XC_FLAGS_HAVE_VXC)/=0))then
      ns = 2 * this%spin - 1
      do ip = 1, this%mesh%np, this%nblock
        nb = min(this%mesh%np-ip+1, this%nblock)
        call interface_xc_block_density(nb, density(ip:,:), l_dens)
        call interface_xc_block_gradient(nb, ns, gradient(ip:,:,:), l_sigma)
        call XC_F90(gga_vxc)(this%conf, nb, l_dens(1,1), l_sigma(1,1), l_dedd(1,1), l_vsigma(1,1))
        call interface_xc_deblock_dedgd(nb, ns, l_vsigma, gradient(ip:,:,:), dedg(ip:,:,:))
        call interface_xc_deblock_vxc(nb, l_dedd, vxc(ip:,:))
      end do
    end if

    POP_SUB(interface_xc_block_gga_vxc)
  end subroutine interface_xc_block_gga_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_block_gga_exc_vxc(this, density, gradient, exc, vxc, dedg)
    type(interface_xc_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:),   intent(in)  :: density
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    real(kind=wp), dimension(:),     intent(out) :: exc
    real(kind=wp), dimension(:,:),   intent(out) :: vxc
    real(kind=wp), dimension(:,:,:), intent(out) :: dedg

    real(kind=wp), dimension(this%spin,    this%nblock) :: l_dens, l_dedd
    real(kind=wp), dimension(2*this%spin-1,this%nblock) :: l_sigma
    real(kind=wp), dimension(2*this%spin-1,this%nblock) :: l_vsigma
    real(kind=wp), dimension(              this%nblock) :: l_zk
    integer                                             :: ip, ns, nb

    PUSH_SUB(interface_xc_block_gga_exc_vxc)

    exc = 0.0_wp
    vxc = 0.0_wp
    if((this%id>XC_NONE)                           &
      .and.(iand(this%flags,XC_FLAGS_HAVE_EXC)/=0) &
      .and.(iand(this%flags,XC_FLAGS_HAVE_VXC)/=0))then
      ns = 2 * this%spin - 1
      do ip = 1, this%mesh%np, this%nblock
        nb = min(this%mesh%np-ip+1, this%nblock)
        call interface_xc_block_density(nb, density(ip:,:), l_dens)
        call interface_xc_block_gradient(nb, ns, gradient(ip:,:,:), l_sigma)
        call XC_F90(gga_exc_vxc)(this%conf, nb, l_dens(1,1), l_sigma(1,1), l_zk(1), l_dedd(1,1), l_vsigma(1,1))
        call interface_xc_deblock_dedgd(nb, ns, l_vsigma, gradient(ip:,:,:), dedg(ip:,:,:))
        call interface_xc_deblock_vxc(nb, l_dedd, vxc(ip:,:))
        exc(ip:ip+nb-1) = l_zk(1:nb)
      end do
    end if

    POP_SUB(interface_xc_block_gga_exc_vxc)
  end subroutine interface_xc_block_gga_exc_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_gga_vxc_correction(this, dedg, vxc)
    type(interface_xc_t),            intent(in)    :: this
    real(kind=wp), dimension(:,:,:), intent(inout) :: dedg !Should be intent in only
    real(kind=wp), dimension(:,:),   intent(inout) :: vxc

    real(kind=wp), dimension(this%mesh%np_part) :: div
    integer                                     :: is

    PUSH_SUB(interface_xc_gga_vxc_correction)

    ! subtract the divergence of the functional derivative of Exc with respect to
    ! the gradient of the density.
    do is = 1, this%spin
      call dderivatives_div(this%der, dedg(:,:,is), div)
      vxc(:,is) = vxc(:,is) - div
    end do

    POP_SUB(interface_xc_gga_vxc_correction)
  end subroutine interface_xc_gga_vxc_correction

  ! ---------------------------------------------------------
  subroutine interface_xc_gga_exc(this, density, exc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc

    real(kind=wp), dimension(this%mesh%np_part,          this%spin) :: dens
    real(kind=wp), dimension(this%mesh%np_part,this%ndim,this%spin) :: grad
    integer                                                         :: isp

    PUSH_SUB(interface_xc_gga_exc)

    call interface_xc_density(density, dens)
    do isp = 1, this%spin
      call dderivatives_grad(this%der, dens(:,isp), grad(:,:,isp))
    end do
    call interface_xc_block_gga_exc(this, dens, grad, exc)

    POP_SUB(interface_xc_gga_exc)
  end subroutine interface_xc_gga_exc

  ! ---------------------------------------------------------
  subroutine interface_xc_gga_vxc(this, density, vxc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: vxc

    real(kind=wp), dimension(this%mesh%np_part,          this%spin) :: dens
    real(kind=wp), dimension(this%mesh%np_part,this%ndim,this%spin) :: grad
    real(kind=wp), dimension(this%mesh%np_part,this%ndim,this%spin) :: dedg
    integer                                                         :: isp

    PUSH_SUB(interface_xc_gga_vxc)

    call interface_xc_density(density, dens)
    do isp = 1, this%spin
      call dderivatives_grad(this%der, dens(:,isp), grad(:,:,isp))
    end do
    call interface_xc_block_gga_vxc(this, dens, grad, vxc, dedg)
    dedg(this%mesh%np+1:,:,:) = 0.0_wp
    call interface_xc_gga_vxc_correction(this, dedg, vxc)

    POP_SUB(interface_xc_gga_vxc)
  end subroutine interface_xc_gga_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_gga_exc_vxc(this, density, exc, vxc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc

    real(kind=wp), dimension(this%mesh%np_part,          this%spin) :: dens
    real(kind=wp), dimension(this%mesh%np_part,this%ndim,this%spin) :: grad
    real(kind=wp), dimension(this%mesh%np_part,this%ndim,this%spin) :: dedg
    integer                                                         :: isp

    PUSH_SUB(interface_xc_gga_exc_vxc)

    call interface_xc_density(density, dens)
    do isp = 1, this%spin
      call dderivatives_grad(this%der, dens(:,isp), grad(:,:,isp))
    end do
    call interface_xc_block_gga_exc_vxc(this, dens, grad, exc, vxc, dedg)
    dedg(this%mesh%np+1:,:,:) = 0.0_wp
    call interface_xc_gga_vxc_correction(this, dedg, vxc)

    POP_SUB(interface_xc_gga_exc_vxc)
  end subroutine interface_xc_gga_exc_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_copy(this, that)
    type(interface_xc_t), intent(inout) :: this
    type(interface_xc_t), intent(in)    :: that

    PUSH_SUB(interface_xc_copy)

    call interface_xc_end(this)
    if(that%id>XC_NONE)then
      call interface_xc_init(this, that%id, that%nspin, that%nblock)
      if(that%ndim>0)then
        call interface_xc_start(this, that%mesh, that%der)
      end if
    end if

    POP_SUB(interface_xc_copy)
  end subroutine interface_xc_copy

  ! ---------------------------------------------------------
  subroutine interface_xc_end(this)
    type(interface_xc_t), intent(inout) :: this

    PUSH_SUB(interface_xc_end)

    nullify(this%mesh,this%der)
    this%nblock = 0
    this%ndim = 0
    this%nspin = 0
    if(this%id>XC_NONE)call XC_F90(func_end)(this%conf)
    this%id = XC_NONE
    this%family = XC_NONE
    this%spin = XC_NONE
    this%kind = XC_NONE
    this%flags = XC_NONE

    POP_SUB(interface_xc_end)
  end subroutine interface_xc_end

end module interface_xc_m

!! Local Variables:
!! mode: f90
!! End:
