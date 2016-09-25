
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

  type, bind(c), public :: element_t
     real (kind=c_double) :: metdet(np,np)
     real (kind=c_double) :: Dinv(2,2,np,np)
     real (kind=c_double) :: rmetdet(np,np)
  end type element_t

  type, bind(c), public :: derivative_t
     real (kind=c_double) :: Dvv(np,np)
  end type derivative_t

  real (kind=c_double), parameter, private :: rrearth = 1.5683814303638645D-7

contains
  subroutine divergence_sphere_fortran(v,deriv,elem,div) bind(c)
    !
    !   input:  v = velocity in lat-lon coordinates
    !   ouput:  div(v)  spherical divergence of v
    !
    real(kind=c_double), intent(in) :: v(2,np,np)  ! in lat-lon coordinates
    type (derivative_t), intent(in) :: deriv
    type (element_t), intent(in) :: elem
    real(kind=c_double), intent(out) :: div(np,np)

    ! Local
    integer i
    integer j
    integer l

    real(kind=c_double) ::  dudx00
    real(kind=c_double) ::  dvdy00
    real(kind=c_double) ::  gv(np,np,2),vvtemp(np,np)

    ! convert to contra variant form and multiply by g
    do j=1,np
       do i=1,np
          gv(i,j,1)=elem%metdet(i,j)*(elem%Dinv(1,1,i,j)*v(1,i,j) + elem%Dinv(2,1,i,j)*v(2,i,j))
          gv(i,j,2)=elem%metdet(i,j)*(elem%Dinv(1,2,i,j)*v(1,i,j) + elem%Dinv(2,2,i,j)*v(2,i,j))
       enddo
    enddo

    ! compute d/dx and d/dy         
    do j=1,np
       do l=1,np
          dudx00=0.0d0
          dvdy00=0.0d0
          !DIR$ UNROLL(NP)
          do i=1,np
             dudx00 = dudx00 + deriv%Dvv(i,l  )*gv(i,j  ,1)
             dvdy00 = dvdy00 + deriv%Dvv(i,l  )*gv(j  ,i,2)
          end do
          div(l  ,j  ) = dudx00
          vvtemp(j  ,l  ) = dvdy00
       end do
    end do

    div(:,:)=(div(:,:)+vvtemp(:,:))*(elem%rmetdet(:,:)*rrearth)
  end subroutine divergence_sphere_fortran

end module divergence
