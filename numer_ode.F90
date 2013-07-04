!!
!! Ordinary Differential Equations Solvers
!!
!! taken form NOSE, originally written by Tomas Mancal
!!
module numer_ode

	use std_types

	implicit none

  !!
  !! 4-th order Runge-Kutta method
  !!
  interface ode_rk4

     subroutine rk4_sp(y,dydx,x,h,yout,derivs)
       real(SP), dimension(:), intent(IN) :: y,dydx
       real(SP), intent(IN) :: x,h
       real(SP), dimension(:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx)
            real(SP), intent(IN) :: x
            real(SP), dimension(:), intent(IN) :: y
            real(SP), dimension(:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_sp

     subroutine rk4_dp(y,dydx,x,h,yout,derivs)
       real(DP), dimension(:), intent(IN) :: y,dydx
       real(DP), intent(IN) :: x,h
       real(DP), dimension(:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx)
            real(DP), intent(IN) :: x
            real(DP), dimension(:), intent(IN) :: y
            real(DP), dimension(:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_dp

     subroutine rk4_dpc(y,dydx,x,h,yout,derivs)
       complex(DPC), dimension(:), intent(IN) :: y,dydx
       real(DP), intent(IN) :: x,h
       complex(DPC), dimension(:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx)
            real(DP), intent(IN) :: x
            complex(DPC), dimension(:), intent(IN) :: y
            complex(DPC), dimension(:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_dpc

     subroutine rk4_dpc_mat(y,dydx,x,h,yout,derivs)
       complex(DPC), dimension(:,:), intent(IN) :: y,dydx
       real(DP), intent(IN) :: x,h
       complex(DPC), dimension(:,:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx)
            real(DP), intent(IN) :: x
            complex(DPC), dimension(:,:), intent(IN) :: y
            complex(DPC), dimension(:,:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_dpc_mat

     subroutine rk4_dpc_mat3(y,dydx,x,h,yout,derivs)
       complex(DPC), dimension(:,:,:), intent(IN)  :: y,dydx
       real(DP), intent(IN) :: x,h
       complex(DPC), dimension(:,:,:), intent(OUT) :: yout
       interface
          subroutine derivs(x,y,dydx)
            real(DP), intent(IN) :: x
            complex(DPC), dimension(:,:,:), intent(IN)  :: y
            complex(DPC), dimension(:,:,:), intent(OUT) :: dydx
          end subroutine derivs
       end interface
     end subroutine rk4_dpc_mat3

  end interface

contains

	! double precision complex matrix version
	subroutine rk4_dpc_mat(y,dydx,x,h,yout,derivs)
		complex(DPC), dimension(:,:), intent(IN)  :: y,dydx
		real(DP),                     intent(IN)  :: x,h
		complex(DPC), dimension(:,:), intent(OUT) :: yout
		interface
			subroutine derivs(x,y,dydx)
				real(DP), intent(IN) :: x
				complex(DPC), dimension(:,:), intent(IN)  :: y
				complex(DPC), dimension(:,:), intent(OUT) :: dydx
			end subroutine derivs
		end interface
		integer(I4B) :: ndum
		real(DP) :: h6,hh,xh
		complex(DPC), dimension(size(y,1),size(y,2)) :: dym,dyt,yt

		hh=h*0.5_dp
		h6=h/6.0_dp
		xh=x+hh
		yt=y+hh*dydx
		call derivs(xh,yt,dyt)
		yt=y+hh*dyt
		call derivs(xh,yt,dym)
		yt=y+h*dym
		dym=dyt+dym
		call derivs(x+h,yt,dyt)
		yout=y+h6*(dydx+dyt+2.0_dp*dym)

	end subroutine rk4_dpc_mat

	! double precision complex matrix with 3 indices version
	subroutine rk4_dpc_mat3(y,dydx,x,h,yout,derivs)
		complex(DPC), dimension(:,:,:), intent(IN)  :: y,dydx
		real(DP),                     intent(IN)  :: x,h
		complex(DPC), dimension(:,:,:), intent(OUT) :: yout
		interface
			subroutine derivs(x,y,dydx)
				real(DP), intent(IN) :: x
				complex(DPC), dimension(:,:,:), intent(IN)  :: y
				complex(DPC), dimension(:,:,:), intent(OUT) :: dydx
			end subroutine derivs
		end interface
		integer(I4B) :: ndum
		real(DP) :: h6,hh,xh
		complex(DPC), dimension(size(y,1),size(y,2),size(y,3)) :: dym,dyt,yt

		hh=h*0.5_dp
		h6=h/6.0_dp
		xh=x+hh
		yt=y+hh*dydx
		call derivs(xh,yt,dyt)
		yt=y+hh*dyt
		call derivs(xh,yt,dym)
		yt=y+h*dym
		dym=dyt+dym
		call derivs(x+h,yt,dyt)
		yout=y+h6*(dydx+dyt+2.0_dp*dym)

	end subroutine rk4_dpc_mat3

! double precision complex version
subroutine rk4_dpc(y,dydx,x,h,yout,derivs)
  complex(DPC), dimension(:), intent(IN)  :: y,dydx
  real(DP),                   intent(IN)  :: x,h
  complex(DPC), dimension(:), intent(OUT) :: yout
  interface
     subroutine derivs(x,y,dydx)
       real(DP), intent(IN) :: x
       complex(DPC), dimension(:), intent(IN)  :: y
       complex(DPC), dimension(:), intent(OUT) :: dydx
     end subroutine derivs
  end interface
  integer(I4B) :: ndum
  real(DP) :: h6,hh,xh
  complex(DPC), dimension(size(y)) :: dym,dyt,yt
  hh=h*0.5_sp
  h6=h/6.0_sp
  xh=x+hh
  yt=y+hh*dydx
  call derivs(xh,yt,dyt)
  yt=y+hh*dyt
  call derivs(xh,yt,dym)
  yt=y+h*dym
  dym=dyt+dym
  call derivs(x+h,yt,dyt)
  yout=y+h6*(dydx+dyt+2.0_dp*dym)
end subroutine rk4_dpc

! double precision version
subroutine rk4_dp(y,dydx,x,h,yout,derivs)
  real(DP), dimension(:), intent(IN)  :: y,dydx
  real(DP),                  intent(IN)  :: x,h
  real(DP), dimension(:), intent(OUT) :: yout
  interface
     subroutine derivs(x,y,dydx)
       real(DP), intent(IN) :: x
       real(DP), dimension(:), intent(IN)  :: y
       real(DP), dimension(:), intent(OUT) :: dydx
     end subroutine derivs
  end interface
  integer(I4B) :: ndum
  real(DP) :: h6,hh,xh
  real(DP), dimension(size(y)) :: dym,dyt,yt
  hh=h*0.5_sp
  h6=h/6.0_sp
  xh=x+hh
  yt=y+hh*dydx
  call derivs(xh,yt,dyt)
  yt=y+hh*dyt
  call derivs(xh,yt,dym)
  yt=y+h*dym
  dym=dyt+dym
  call derivs(x+h,yt,dyt)
  yout=y+h6*(dydx+dyt+2.0_dp*dym)
end subroutine rk4_dp

! single precision version
subroutine rk4_sp(y,dydx,x,h,yout,derivs)
  real(SP), dimension(:), intent(IN) :: y,dydx
  real(SP), intent(IN) :: x,h
  real(SP), dimension(:), intent(OUT) :: yout
  interface
     subroutine derivs(x,y,dydx)
       real(SP), intent(IN) :: x
       real(SP), dimension(:), intent(IN) :: y
       real(SP), dimension(:), intent(OUT) :: dydx
     end subroutine derivs
  end interface
  integer(I4B) :: ndum
  real(SP) :: h6,hh,xh
  real(SP), dimension(size(y)) :: dym,dyt,yt
  hh=h*0.5_sp
  h6=h/6.0_sp
  xh=x+hh
  yt=y+hh*dydx
  call derivs(xh,yt,dyt)
  yt=y+hh*dyt
  call derivs(xh,yt,dym)
  yt=y+h*dym
  dym=dyt+dym
  call derivs(x+h,yt,dyt)
  yout=y+h6*(dydx+dyt+2.0_sp*dym)
end subroutine rk4_sp



end module numer_ode
