
! elem
!   metdet
!     real_kind (np, np)
!   Dinv
!     real_kind (np, np, 2, 2)
!   rmetdet
!     real_kind (np, np)


! deriv
!   Dvv
!     real_kind (np, np)

module divergence

  use iso_c_binding
  implicit none

  integer, parameter, public :: np = 4
  integer, parameter, public :: dim = 2

  type, bind(c), public :: element_t
     real (kind=c_double) :: metdet(np,np)
     real (kind=c_double) :: Dinv(dim,dim,np,np)
     real (kind=c_double) :: rmetdet(np,np)
  end type element_t

  type, bind(c), public :: derivative_t
     real (kind=c_double) :: Dvv(np,np)
  end type derivative_t

  real (kind=c_double), parameter, private :: rrearth = 1.5683814303638645D-7

  interface
     subroutine puts(s) bind(c)
       ! For Debugging: Use as follows:
       ! character (len=*), parameter :: strout="Printing from Fortran"
       ! call puts(strout)
       use iso_c_binding
       character :: s(*)
     end subroutine puts
  end interface

contains
  subroutine divergence_sphere_fortran(v,Dvv,metdet,Dinv,rmetdet,div) bind(c)
    !
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  div(v)  spherical divergence of v
    !
    real(kind=c_double), intent(in) :: v(dim,np,np)  ! in lat-lon coordinates
    real(kind=c_double), intent(in) :: Dvv(np,np)
    real(kind=c_double), intent(in) :: metdet(np,np)
    real(kind=c_double), intent(in) :: Dinv(dim,dim,np,np)
    real(kind=c_double), intent(in) :: rmetdet(np,np)
    real(kind=c_double), intent(out) :: div(np,np)

    ! Local
    integer i
    integer j
    integer l

    real(kind=c_double) ::  dudx00
    real(kind=c_double) ::  dvdy00
    real(kind=c_double) ::  gv(np,np,2),vvtemp(np,np)
#if 1
    ! convert to contra variant form and multiply by g
    do j=1,np
       do i=1,np
          do l=1,2
             gv(i,j,l)=metdet(i,j)*(Dinv(1,l,i,j)*v(1,i,j) + Dinv(2,l,i,j)*v(2,i,j))
          enddo
       enddo
    enddo
#endif
    ! compute d/dx and d/dy
#if 1
    do l=1,np
       do j=1,np
          dudx00=0.0d0
          dvdy00=0.0d0
          !DIR$ UNROLL(NP)
          do i=1,np
             dudx00 = dudx00 + Dvv(i,l  )*gv(i,j  ,1)
             dvdy00 = dvdy00 + Dvv(i,l  )*gv(j  ,i,2)
          end do
          div(l  ,j  ) = dudx00
          vvtemp(j  ,l  ) = dvdy00
       end do
    end do
#endif
#if 1
    div(:,:)=(div(:,:)+vvtemp(:,:))*(rmetdet(:,:)*rrearth)
#endif
  end subroutine divergence_sphere_fortran

end module divergence
