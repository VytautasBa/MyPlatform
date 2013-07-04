!
! system parameters (hardcoded for now) + internal resources
!
module resources

	use std_types
	use numer_matrix

	implicit none

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!
	!! external / primary parameters
	!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! method parameters
	!!
	integer, parameter :: Ks = 1 ! number of Bose-Einstein poles, i.e. there are Ks+1 poles in total
	integer(sp), parameter :: Ls = 15 ! depth of the hierarchy: # tiers ABOVE 0'th, i.e. there are Ls+1 ADO's in total
	! these need to be 'parameters' because some non-allocatable matrices they give dimensions to

	!! system + system-bath parameters
	!!
	integer, parameter :: Ns = 2 ! number of sites

	real(dp), dimension(:), allocatable :: en
	complex(dpc), dimension(:,:), allocatable :: cpl	! cpl

	real(dp), dimension(:), allocatable :: dd			! dipole moments

	real(dp), dimension(:), allocatable :: re_en		! reorganization energies
	real(dp), dimension(:), allocatable :: gamma_i		! fluctuation/dissipation time scales

	!! control + other parameters
	!!
	real(dp) 	:: dt	! time step in fs
	integer(sp) :: gt1	! # time steps per grid step; e.g., dg1 * dt = const. controls the density of integration
	integer(sp) :: gt2	! # time steps per population grid step
	integer(sp) :: Nt1	! # grid steps; i.e., duration of an experiment = Nt1 * gt1 * dt
	integer(sp) :: Nt2	! # grid steps for population time (again, Nt*gt = control over density of points)

	integer(sp) :: NFFT_basic = 1024 ! number of points in FFT input data

	real(dp) :: rwa = 0.0_dp			! the RWA frequency (e.g. to do shifts by optical frequencies in FFT)
	real(dp) :: Temperature = 70.0_dp	! given in Kelvin

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!
	!! internal / derivative parameters
	!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	real(dp) :: beta_int

	integer :: NN ! number of ADO's

	complex(dpc), dimension(:,:), allocatable  :: HH ! Hamiltonian

	complex(dpc), dimension(:,:,:), allocatable :: ADO	! the Auxiliary Density Operator

	integer, dimension(:,:,:), allocatable :: index_n	! INDEX[I,n,k]
	integer, dimension(:,:,:), allocatable :: index_p1	! INDEX_p1[I,n,k]
	integer, dimension(:,:,:), allocatable :: index_m1

	complex(dpc), dimension(:,:), allocatable :: c_ik		! correlation function coefficients (i= site #,k= pole #)
	real(dp), dimension(:,:), allocatable     :: gamma_ik 	! gamma's = cor. f. poles
	real(dp), dimension(:), allocatable       :: delta_i	! delta's (sites)

	complex(dpc), dimension(:,:), allocatable :: corr_fun	! correlation functions
	complex(dpc), dimension(:,:), allocatable :: goft		! g(t) functions

	real(dp) :: R_Ks

	character(len=256) :: out_dir

contains

	subroutine set_values()

		print *, "Initializing MyPlatform..."

		!!
		!! general parameters:
		!!

		! 'operational'
		out_dir="/home/workspace/myplatform/out"

		! 'calculational'
		beta_int = 1/(Temperature * kB_intK)

		dt  = 1.0_dp
		gt1 = 1
		gt2 = 1
		Nt1 = 1024 ! better keep it a power of 2 when using Spectroscopy
		Nt2 = 151

		! 'system'
		call set_system_parameters()

		!!
		!! HEOM-speciffic resources
		!!
		NN = binomial(Ns*(Ks+1)+Ls,Ls)

		allocate(ADO(-1:NN,Ns,Ns))
		ADO = (0.0_dp,0.0_dp)

		allocate( index_n(0:NN-1,Ns,Ks+1))
		allocate(index_p1(0:NN-1,Ns,Ks+1))
		allocate(index_m1(0:NN-1,Ns,Ks+1))

	end subroutine set_values



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! SYSTEM PARAMETERS
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine set_system_parameters()
		integer :: i

		print *, "Setting system parameters..."

		if (.not. allocated(en)) then
			allocate(en(Ns))
		end if
		en = 0.0_dp

		if (.not. allocated(dd)) then
			allocate(dd(Ns))
		end if
		dd = 0.0_dp

		if (.not. allocated(cpl)) then
			allocate(cpl(Ns,Ns))
		end if
		cpl = 0.0_dp

		if (.not. allocated(HH)) then
			allocate(HH(Ns,Ns))
		end if
		HH = 0.0_dp

		if (.not. allocated(re_en)) then
			allocate(re_en(Ns))
		end if
		re_en = 0.0_dp

		if (.not. allocated(gamma_i)) then
			allocate(gamma_i(Ns))
		end if
		gamma_i = 0.0_dp

		!!! hardcoding for now
		en(1) = 0.0_dp * Energy_cm_to_internal
		en(2) = 400.0_dp * Energy_cm_to_internal

		dd(1) = 0.5_dp
		dd(2) = 0.5_dp

		re_en(1) = 30.0_dp * Energy_cm_to_internal !30
		re_en(2) = 30.0_dp * Energy_cm_to_internal !30

		gamma_i(1) = 1.0_dp/100.0_dp
		gamma_i(2) = 1.0_dp/100.0_dp

		cpl(1,2) = 80.0_dp * Energy_cm_to_internal !80
		cpl(2,1) = cpl(1,2)

		HH = cpl
		forall(i= 1:Ns)
			HH(i,i) = HH(i,i) + en(i)
		end forall

	end subroutine set_system_parameters



	subroutine set_corr_functions()
		real(dp), dimension(:,:), allocatable :: LL, LL1	! auxiliary matrices
		real(dp), dimension(:,:), allocatable :: SL, LLd, SL1, LL1d	! matrices for diagonalization of aux. m.
		real(dp), dimension(:), allocatable :: fi, eta		! Padé poles; Padé coefficients
		real(dp), dimension(:), allocatable :: ev, ev1		! eigenvalues of aux. matrices
		integer :: k, m, t

		if (.not. allocated(c_ik)) then
			allocate(c_ik(Ns,Ks+1))
		end if

		if (.not. allocated(gamma_ik)) then
			allocate(gamma_ik(Ns,Ks+1))
		end if

		if (.not. allocated(delta_i)) then
			allocate(delta_i(Ns))
		end if

		if (.not. allocated(corr_fun)) then
			allocate(corr_fun(Ns,Nt1*gt1))
		end if

		if (.not. allocated(goft)) then
			allocate(goft(Ns,-Nt1:Nt1+2*Nt2))
		end if

		R_Ks = 1.0_dp / (4.0_dp*(Ks+1)*(2*Ks+3))

		allocate(fi(Ks),eta(Ks))
		allocate(LL(2*Ks+1,2*Ks+1),SL(2*Ks+1,2*Ks+1),LLd(2*Ks+1,2*Ks+1))
		allocate(LL1(2*Ks,2*Ks),SL1(2*Ks,2*Ks),LL1d(2*Ks,2*Ks))
		allocate(ev(Ks),ev1(Ks))

		LL = 0.0_dp
		do m = 1, 2*Ks
			LL(m, m+1) = 1/(sqrt(2*real(m)+1)*sqrt(2*real(m)+3))
			LL(m+1, m) = LL(m, m+1)
		end do

		LL1 = 0.0_dp
		do m = 1, 2*Ks-1
			LL1(m, m+1) = 1/(sqrt(2*real(m)+3)*sqrt(2*real(m)+5))
			LL1(m+1, m) = LL1(m, m+1)
		end do

		call spec(LL,SL,LLd)		! 2*Ks+1 eigenvalues: 0 + Ks unique (with both signs)
		call spec(LL1,SL1,LL1d)		! 2*Ks eigenvalues: Ks unique (with both signs)

		! picking out the unique eigenvalues (they are sorted, which helps a lot)
		do m = 1, Ks
			ev(m) = LLd(int((2*Ks+1)/2)+1+m, int((2*Ks+1)/2)+1+m) ! taking elements above the middle one, which is 0
			ev1(m)= LL1d(Ks+m, Ks+m) ! taking elements from the "upper" part of LL1d, because of even number
		end do

		fi = 2.0_dp /(beta_int * ev)
		eta = 0.5_dp * R_Ks
		do k = 1, Ks
			do m = 1 ,Ks
				eta(k) = eta(k) * 4.0_dp * (1/ev1(m)**2 - 1/ev(k)**2)
				if (m.ne.k) then
					eta(k) = eta(k) / (   4.0_dp * (1/ev(m)**2 - 1/ev(k)**2)   )
				end if
			end do
		end do

		! prepare output values
		c_ik = (0.0_dp,0.0_dp)
		c_ik(1:Ns,1) = - (0.0_dp,1.0_dp) * gamma_i(1:Ns) * re_en(1:Ns) ! the rest of the imaginary parts are 0
		c_ik(1:Ns,1) = c_ik(1:Ns,1) + 2.0_dp*re_en(1:Ns)/beta_int - 				 &
					   2*R_Ks*beta_int*re_en(1:Ns)*gamma_i(1:Ns)**2
		do m = 1, Ns
			do k = 1, Ks
				c_ik(m,1) = c_ik(m,1) + 4.0_dp*eta(k)*re_en(m)*gamma_i(m)**2 /		 &
							(beta_int*(gamma_i(m)**2 - fi(k)**2))
				! for c_ik(m,k) of k > 1 we simply have
				c_ik(m,k+1) = 4.0_dp*eta(k)*re_en(m)*fi(k)*gamma_i(m) / 			 &
							  (beta_int*(fi(k)**2 - gamma_i(m)**2))
			end do
		end do

!		print *, "C_ik: ", real(c_ik(1,1)), aimag(c_ik(1,1))
!		print *, "C_ik: ", real(c_ik(1,2)), aimag(c_ik(1,2))
!		print *, "C_ik: ", real(c_ik(1,3)), aimag(c_ik(1,3))
!		print *, "C_ik: ", real(c_ik(2,1)), aimag(c_ik(2,1))
!		print *, "C_ik: ", real(c_ik(2,2)), aimag(c_ik(2,2))
!		print *, "C_ik: ", real(c_ik(2,3)), aimag(c_ik(2,3))

		gamma_ik(1:Ns,1) = gamma_i(1:Ns)
		forall(m = 1:Ns)
			gamma_ik(m,2:Ks+1) = fi(1:Ks)
		end forall

		delta_i = 4.0_dp * R_Ks * re_en * gamma_i * beta_int

		! create & write out correlation/g(t) functions
!		open(unit=11,file=(trim(out_dir)//"/"//trim("corr_fun.dat")), position='append')
!		open(unit=12,file=(trim(out_dir)//"/"//trim("goft.dat")), position='append')
		goft = 0.0_dp
		do t = 1, Nt1+2*Nt2
			do m = 1, Ns
				do k = 1, Ks+1
!					corr_fun(m,t) = corr_fun(m,t) + c_ik(m,k) * exp(-gamma_ik(m,k)*dt*t)
					goft(m,t) = goft(m,t) + &
							    (c_ik(m,k)/gamma_ik(m,k)) * (t*gt1*dt+( exp(-gamma_ik(m,k)*t*gt1*dt)-1 )/gamma_ik(m,k) )
				end do
			end do
!		write(11,'(5f16.8)') t*dt, real(corr_fun(1,t)),aimag(corr_fun(1,t)),real(corr_fun(2,t)),aimag(corr_fun(2,t))
!		write(12,'(5f16.8)') t*dt, real(goft(1,t)),aimag(goft(1,t)),real(goft(2,t)),aimag(goft(2,t))
		end do
		do	m = 1, Ns
		do  t = 1, Nt1
			goft(m,-t) = conjg(goft(m,t))
		end do
		end do
!		close(unit=11)
!		close(unit=12)

		! clean-up
		deallocate(LL, SL, LLd)
		deallocate(LL1, SL1, LL1d)
		deallocate(fi, eta, ev, ev1)

	end subroutine set_corr_functions



	subroutine cleanup_platform()

		print *, "Cleaning up MyPlatform..."

		if (allocated(en)) then
			deallocate(en)
		end if

		if (allocated(dd)) then
			deallocate(dd)
		end if

		if (allocated(re_en)) then
			deallocate(re_en)
		end if

		if (allocated(gamma_i)) then
			deallocate(gamma_i)
		end if

		if (allocated(cpl)) then
			deallocate(cpl)
		end if

		if (allocated(HH)) then
			deallocate(HH)
		end if

		if (allocated(index_n)) then
			deallocate(index_n)
		end if

		if (allocated(index_p1)) then
			deallocate(index_p1)
		end if

		if (allocated(index_m1)) then
			deallocate(index_m1)
		end if

		if (allocated(corr_fun)) then
			deallocate(corr_fun)
		end if

		if (allocated(goft)) then
			deallocate(goft)
		end if

		if (allocated(c_ik)) then
			deallocate(c_ik)
		end if

		if (allocated(gamma_ik)) then
			deallocate(gamma_ik)
		end if

		if (allocated(delta_i)) then
			deallocate(delta_i)
		end if

		if (allocated(ADO)) then
			deallocate(ADO)
		end if

	end subroutine cleanup_platform



end module resources
