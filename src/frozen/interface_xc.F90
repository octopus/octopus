#include "global.h"

module interface_xc_m

  use global_m
  use messages_m
  use profiling_m

  use XC_F90(lib_m)

  use json_m,  only: JSON_OK, json_object_t, json_get
  use kinds_m, only: wp

  implicit none

  private
  public ::                     &
    interface_xc_init,          &
    interface_xc_set_parameter, &
    interface_xc_get_kind,      &
    interface_xc_density,       &
    interface_xc_lda_exc,       &
    interface_xc_lda_exc_vxc,   &
    interface_xc_lda_vxc,       &
    interface_xc_gga_exc,       &
    interface_xc_gga_exc_vxc,   &
    interface_xc_gga_vxc,       &
    interface_xc_end

  integer, public, parameter :: XC_UNKNOWN =-1
  integer, public, parameter :: XC_NONE    = 0

  integer, public, parameter :: INT_EXP_SCREENED = 0
  integer, public, parameter :: INT_SOFT_COULOMB = 1

  integer, parameter :: n_block = 1000

  integer, dimension(3), parameter :: xc_flags_nd=(/XC_FLAGS_1D, XC_FLAGS_2D, XC_FLAGS_3D/)


  type, public :: interface_xc_t
    private
    integer                 :: kind  = XC_NONE ! exchange, correlation, or exchange-correlation
    integer                 :: id    = XC_NONE ! identifier
    integer                 :: spin  = XC_NONE ! XC_UNPOLARIZED | XC_POLARIZED
    integer                 :: flags = XC_NONE ! XC_FLAGS_HAVE_EXC + XC_FLAGS_HAVE_VXC + ...
    integer                 :: nspin = 0
    integer                 :: ndim  = 0
    type(XC_F90(pointer_t)) :: conf            ! the pointer used to call the library
    type(XC_F90(pointer_t)) :: info            ! information about the functional
  end type interface_xc_t

contains

  ! ---------------------------------------------------------
  subroutine interface_xc_init(this, id, ndim, nspin)
    type(interface_xc_t), intent(out) :: this
    integer,              intent(in)  :: id
    integer,              intent(in)  :: ndim
    integer,              intent(in)  :: nspin
    !
    integer :: family
    !
    ASSERT(id>XC_UNKNOWN)
    if(id>XC_NONE)then
      this%id=id
      ASSERT(nspin>0)
      this%nspin=nspin
      ASSERT(ndim>0)
      this%ndim=ndim
      this%spin=XC_UNPOLARIZED
      if(nspin>1)this%spin=XC_POLARIZED
      call XC_F90(func_init)(this%conf, this%info, this%id, this%spin)
      this%kind=XC_F90(info_kind)(this%info)
      this%flags=XC_F90(info_flags)(this%info)
      if(iand(this%flags,xc_flags_nd(this%ndim))==0)then
        write(unit=message(1), fmt="(a,i1.1,a2)") "Cannot use the specified functional in ", this%ndim, "D."
        call messages_fatal(1)
      end if
    end if
    return
  end subroutine interface_xc_init

  ! ---------------------------------------------------------
  elemental function interface_xc_is_unpolarized(this) result(that)
    type(interface_xc_t), intent(in) :: this
    !
    logical :: that
    !
    that=(this%spin==XC_UNPOLARIZED)
    return
  end function interface_xc_is_unpolarized

  ! ---------------------------------------------------------
  elemental function interface_xc_is_polarized(this) result(that)
    type(interface_xc_t), intent(in) :: this
    !
    logical :: that
    !
    that=(this%spin==XC_POLARIZED)
    return
  end function interface_xc_is_polarized

  ! ---------------------------------------------------------
  elemental function interface_xc_get_kind(this) result(kind)
    type(interface_xc_t), intent(in) :: this
    !
    integer :: kind
    !
    kind=this%kind
    return
  end function interface_xc_get_kind

  ! ---------------------------------------------------------
  pure subroutine interface_xc_density(this, rho_out, rho_in)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(out) :: rho_out
    real(kind=wp), dimension(:,:), intent(in)  :: rho_in
    !
    real(kind=wp) :: a, b, c, d, dtot, dpol
    integer       :: i, np
    !
    np=size(rho_in, dim=1)
    do i = 1, np
      select case(this%nspin)
      case(1,2)
        where(rho_in(i,:)<0.0_wp)
          rho_out(i,:)=0.0_wp
        elsewhere
          rho_out(i,:)=rho_in(i,:)
        end where
      case(4) !SPINORS
        a=rho_in(i,1)
        b=rho_in(i,2)
        dtot=a+b
        c=rho_in(i,3)
        d=rho_in(i,4)
        dpol=sqrt((a-b)**2+4.0_wp*(c**2+d**2))
        rho_out(i,1)=max(0.5_wp*(dtot+dpol),0.0_wp)
        rho_out(i,2)=max(0.5_wp*(dtot-dpol),0.0_wp)
      case default
        rho_out(i,:)=-1.0_wp
      end select
    end do
    return
  end subroutine interface_xc_density

  ! ---------------------------------------------------------
  pure subroutine interface_xc_lda_in_density(rho_out, rho_in)
    real(kind=wp), dimension(:,:), intent(out) :: rho_out
    real(kind=wp), dimension(:,:), intent(in)  :: rho_in
    !
    real(kind=wp) :: a, b, c, d, dtot, dpol
    integer       :: i, np, ns
    !
    np=min(size(rho_out, dim=2), size(rho_in, dim=1))
    ns=size(rho_in,dim=2)
    do i = 1, np
      select case(ns)
      case(1,2)
        where(rho_in(i,:)<0.0_wp)
          rho_out(:,i)=0.0_wp
        elsewhere
          rho_out(:,i)=rho_in(i,:)
        end where
      case(4) !SPINORS
        a=rho_in(i,1)
        b=rho_in(i,2)
        dtot=a+b
        c=rho_in(i,3)
        d=rho_in(i,4)
        dpol=sqrt((a-b)**2+4.0_wp*(c**2+d**2))
        rho_out(1,i)=max(0.5_wp*(dtot+dpol),0.0_wp)
        rho_out(2,i)=max(0.5_wp*(dtot-dpol),0.0_wp)
      case default
        rho_out(:,i)=-1.0_wp
      end select
    end do
    return
  end subroutine interface_xc_lda_in_density

  ! ---------------------------------------------------------
  pure subroutine interface_xc_gga_in_density(rho_out, rho_in)
    real(kind=wp), dimension(:,:), intent(out) :: rho_out
    real(kind=wp), dimension(:,:), intent(in)  :: rho_in
    !
    integer :: i, np, ns
    !
    np=min(size(rho_out, dim=2), size(rho_in, dim=1))
    ns=size(rho_in, dim=2)
    select case(ns)
    case(1)
      rho_out(1,1:np)=rho_in(1:np,1)
    case(2)
      rho_out(:,1:np)=transpose(rho_in(1:np,:))
    case default
      rho_out(:,1:np)=-1.0_wp
    end select
    return
  end subroutine interface_xc_gga_in_density

  ! ---------------------------------------------------------
  pure subroutine interface_xc_in_gradient(sigma, gradient)
    real(kind=wp), dimension(:,:),   intent(out) :: sigma
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    !
    real(kind=wp) :: a, b, c, d, dtot, dpol
    integer       :: i, np, ns
    !
    np=min(size(sigma, dim=2), size(gradient, dim=1))
    ns=size(gradient,dim=3)
    do i = 1, np
      sigma(1,i)=dot_product(gradient(i,:,1),gradient(i,:,1))
      if(ns>1)then
        sigma(2,i)=dot_product(gradient(i,:,1),gradient(i,:,2))
        sigma(3,i)=dot_product(gradient(i,:,2),gradient(i,:,2))
      end if
    end do
    return
  end subroutine interface_xc_in_gradient

  ! ---------------------------------------------------------
  pure subroutine interface_xc_out_vxc(pot_out, rho_in, pot_in)
    real(kind=wp), dimension(:,:), intent(out) :: pot_out
    real(kind=wp), dimension(:,:), intent(in)  :: rho_in
    real(kind=wp), dimension(:,:), intent(in)  :: pot_in
    !
    real(kind=wp) :: va, vb, vtot, vdlt, vpol, ra, rb, rc, rd, rtot, rdlt, rpol
    integer       :: i, np, ns
    !
    np=min(size(pot_out, dim=1), size(pot_in, dim=2))
    ns=size(pot_in, dim=1)
    do i = 1, np
      select case(ns)
      case(1,2)
        pot_out(i,:)=pot_in(:,i)
      case(4) !SPINORS
        ra=rho_in(i,1)
        rb=rho_in(i,2)
        rtot=ra+rb
        rdlt=ra-rb
        rc=rho_in(i,3)
        rd=rho_in(i,4)
        rpol=sqrt(rdlt**2+4.0_wp*(rc**2+rd**2))
        va=pot_in(1,i)
        vb=pot_in(2,i)
        vtot=va+vb
        vdlt=(va-vb)/max(rpol,tiny(rpol))
        vpol=vdlt*rdlt
        pot_out(i,1)=0.5_wp*(vtot+vpol)
        pot_out(i,2)=0.5_wp*(vtot-vpol)
        pot_out(i,3)=vdlt*rc
        pot_out(i,4)=vdlt*rd
      end select
    end do
    return
  end subroutine interface_xc_out_vxc
  
  ! ---------------------------------------------------------
  pure subroutine interface_xc_out_dedgd(dedgd, gradient, vsigma)
    real(kind=wp), dimension(:,:,:), intent(out) :: dedgd
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    real(kind=wp), dimension(:,:),   intent(in)  :: vsigma
    !
    real(kind=wp) :: a, b, c, d, dtot, dpol
    integer       :: i, np, ns
    !
    np=min(size(vsigma, dim=2), size(dedgd, dim=1))
    ns=size(gradient,dim=3)
    do i = 1, np
      if(ns==1)then
        dedgd(i,:,1)=2.0_wp*vsigma(1,i)*gradient(i,:,1)
      else
        dedgd(i,:,1)=2.0_wp*vsigma(1,i)*gradient(i,:,1)+vsigma(2,i)*gradient(i,:,2)
        dedgd(i,:,2)=2.0_wp*vsigma(3,i)*gradient(i,:,2)+vsigma(2,i)*gradient(i,:,1)
      end if
    end do
    return
  end subroutine interface_xc_out_dedgd

  ! ---------------------------------------------------------
  subroutine interface_xc_set_parameter(this, config)
    type(interface_xc_t), intent(inout) :: this
    type(json_object_t),  intent(in)    :: config
    !
    real(kind=wp) :: rtmp
    integer       :: itmp, ierr
    !
    select case(this%id)
    case(XC_LDA_C_XALPHA)
      call json_get(config, "xalpha", rtmp, ierr)
      if(ierr/=JSON_OK)rtmp=1.0_wp
      call XC_F90(lda_c_xalpha_set_par)(this%conf, rtmp)
    case(XC_LDA_X_1D)
      call json_get(config, "Interaction1D", itmp, ierr)
      if(ierr/=JSON_OK)itmp=INT_SOFT_COULOMB
      call json_get(config, "Interaction1DScreening", rtmp, ierr)
      if(ierr/=JSON_OK)rtmp=1.0_wp
      call XC_F90(lda_x_1d_set_par)(this%conf, itmp, rtmp)
    case(XC_LDA_C_1D_CSC)
      call json_get(config, "Interaction1D", itmp, ierr)
      if(ierr/=JSON_OK)itmp=INT_SOFT_COULOMB
      call json_get(config, "Interaction1DScreening", rtmp, ierr)
      if(ierr/=JSON_OK)rtmp=1.0_wp
      call XC_F90(lda_c_1d_csc_set_par)(this%conf, itmp, rtmp)
    case(XC_LDA_C_2D_PRM)
      call json_get(config, "nel", rtmp, ierr)
      ASSERT(ierr==JSON_OK)
      call XC_F90(lda_c_2d_prm_set_par)(this%conf, rtmp)
    end select
    return
  end subroutine interface_xc_set_parameter

  ! ---------------------------------------------------------
  subroutine interface_xc_lda_exc(this, density, exc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    !
    real(kind=wp), dimension(this%spin,n_block) :: l_dens
    real(kind=wp), dimension(n_block)           :: l_zk
    integer                                     :: i, n, np
    !
    exc=0.0_wp
    np=size(exc)
    if(iand(this%flags, XC_FLAGS_HAVE_EXC)/=0)then
      do i = 1, np, n_block
        n=min(np-i+1, n_block)
        call interface_xc_lda_in_density(l_dens, density(i:,:))
        call XC_F90(lda_exc)(this%conf, n, l_dens(1,1), l_zk(1))
        exc(i:i+n-1)=l_zk(1:n)
      end do
    end if
    return
  end subroutine interface_xc_lda_exc

  ! ---------------------------------------------------------
  subroutine interface_xc_lda_exc_vxc(this, density, exc, vxc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:),   intent(out) :: exc
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    real(kind=wp), dimension(this%spin,n_block) :: l_dens, l_dedd
    real(kind=wp), dimension(n_block)           :: l_zk
    integer                                     :: i, n, np
    !
    exc=0.0_wp
    vxc=0.0_wp
    np=size(exc)
    if((iand(this%flags,XC_FLAGS_HAVE_EXC)/=0).and.(iand(this%flags,XC_FLAGS_HAVE_VXC)/=0))then
      do i = 1, np, n_block
        n=min(np-i+1, n_block)
        call interface_xc_lda_in_density(l_dens, density(i:,:))
        call XC_F90(lda_exc_vxc)(this%conf, n, l_dens(1,1), l_zk(1), l_dedd(1,1))
        call interface_xc_out_vxc(vxc(i:,:), density(i:,:), l_dedd)
        !vxc(i:i+n-1,1)=l_dedd(1,1:n)
        exc(i:i+n-1)=l_zk(1:n)
      end do
    end if
    return
  end subroutine interface_xc_lda_exc_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_lda_vxc(this, density, vxc)
    type(interface_xc_t),          intent(in)  :: this
    real(kind=wp), dimension(:,:), intent(in)  :: density
    real(kind=wp), dimension(:,:), intent(out) :: vxc
    !
    real(kind=wp), dimension(this%spin,n_block) :: l_dens, l_dedd
    integer                                     :: i, n, np
    !
    vxc=0.0_wp
    np=size(vxc, dim=1)
    if(iand(this%flags, XC_FLAGS_HAVE_VXC)/=0)then
      do i = 1, np, n_block
        n=min(np-i+1, n_block)
        call interface_xc_lda_in_density(l_dens, density(i:,:))
        call XC_F90(lda_vxc)(this%conf, n, l_dens(1,1), l_dedd(1,1))
        call interface_xc_out_vxc(vxc(i:,:), density(i:,:), l_dedd)
      end do
    end if
    return
  end subroutine interface_xc_lda_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_gga_exc(this, density, gradient, exc)
    type(interface_xc_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:),   intent(in)  :: density
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    real(kind=wp), dimension(:),     intent(out) :: exc
    !
    real(kind=wp), dimension(this%spin,n_block)     :: l_dens
    real(kind=wp), dimension(2*this%spin-1,n_block) :: l_sigma
    real(kind=wp), dimension(n_block)               :: l_zk
    integer                                         :: i, n, np
    !
    exc=0.0_wp
    np=size(exc)
    if(iand(this%flags, XC_FLAGS_HAVE_EXC)/=0)then
      do i = 1, np, n_block
        n=min(np-i+1,n_block)
        call interface_xc_gga_in_density(l_dens, density(i:,:))
        call interface_xc_in_gradient(l_sigma, gradient(i:,:,:))
        call XC_F90(gga_exc)(this%conf, n, l_dens(1,1), l_sigma(1,1), l_zk(1))
        exc(i:i+n-1)=l_zk(1:n)
      end do
    end if
    return
  end subroutine interface_xc_gga_exc

  ! ---------------------------------------------------------
  subroutine interface_xc_gga_exc_vxc(this, density, gradient, exc, vxc, dedg)
    type(interface_xc_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:),   intent(in)  :: density
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    real(kind=wp), dimension(:),     intent(out) :: exc
    real(kind=wp), dimension(:,:),   intent(out) :: vxc
    real(kind=wp), dimension(:,:,:), intent(out) :: dedg
    !
    real(kind=wp), dimension(this%spin,n_block)     :: l_dens, l_dedd
    real(kind=wp), dimension(2*this%spin-1,n_block) :: l_sigma
    real(kind=wp), dimension(2*this%spin-1,n_block) :: l_vsigma
    real(kind=wp), dimension(n_block)               :: l_zk
    integer                                         :: i, n, np
    !
    exc=0.0_wp
    np=size(exc)
    if((iand(this%flags,XC_FLAGS_HAVE_EXC)/=0).and.(iand(this%flags,XC_FLAGS_HAVE_VXC)/=0))then
      do i = 1, np, n_block
        n=min(np-i+1, n_block)
        call interface_xc_gga_in_density(l_dens, density(i:,:))
        call interface_xc_in_gradient(l_sigma, gradient(i:,:,:))
        call XC_F90(gga_exc_vxc)(this%conf, n, l_dens(1,1), l_sigma(1,1), l_zk(1), l_dedd(1,1), l_vsigma(1,1))
        call interface_xc_out_dedgd(dedg(i:,:,:), gradient(i:,:,:), l_vsigma)
        call interface_xc_out_vxc(vxc(i:,:), density(i:,:), l_dedd)
        exc(i:i+n-1)=l_zk(1:n)
      end do
    end if
    return
  end subroutine interface_xc_gga_exc_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_gga_vxc(this, density, gradient, vxc, dedg)
    type(interface_xc_t),            intent(in)  :: this
    real(kind=wp), dimension(:,:),   intent(in)  :: density
    real(kind=wp), dimension(:,:,:), intent(in)  :: gradient
    real(kind=wp), dimension(:,:),   intent(out) :: vxc
    real(kind=wp), dimension(:,:,:), intent(out) :: dedg
    !
    real(kind=wp), dimension(this%spin,n_block)     :: l_dens, l_dedd
    real(kind=wp), dimension(2*this%spin-1,n_block) :: l_sigma
    real(kind=wp), dimension(2*this%spin-1,n_block) :: l_vsigma
    integer                                         :: i, n, np
    !
    vxc=0.0_wp
    np=size(vxc, dim=1)
    if(iand(this%flags, XC_FLAGS_HAVE_VXC)/=0)then
      do i = 1, np, n_block
        n=min(np-i+1, n_block)
        call interface_xc_gga_in_density(l_dens, density(i:,:))
        call interface_xc_in_gradient(l_sigma, gradient(i:,:,:))
        call XC_F90(gga_vxc)(this%conf, n, l_dens(1,1), l_sigma(1,1), l_dedd(1,1), l_vsigma(1,1))
        call interface_xc_out_dedgd(dedg(i:,:,:), gradient(i:,:,:), l_vsigma)
        call interface_xc_out_vxc(vxc(i:,:), density(i:,:), l_dedd)
      end do
    end if
    return
  end subroutine interface_xc_gga_vxc

  ! ---------------------------------------------------------
  subroutine interface_xc_end(this)
    type(interface_xc_t), intent(inout) :: this
    !
    if(this%id>XC_NONE)&
      call XC_F90(func_end)(this%conf)
    this%id=XC_NONE
    return
  end subroutine interface_xc_end

end module interface_xc_m

!! Local Variables:
!! mode: f90
!! End:
