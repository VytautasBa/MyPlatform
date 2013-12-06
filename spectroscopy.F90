module spectroscopy

	use std_types
	use resources
	use numer_interp
	use numer_fft
	use numer_matrix

	implicit none


	complex(dpc), dimension(:), allocatable :: R_trf	! Response function for TRF; R(tau) (later upgrade to R(tau,T))
	real(dp), dimension(:), allocatable     :: F_trf 	! TRF signal F(w) at T (later F(w,T))

	double precision, parameter	:: REGULARIZATION_TIME_PER_TMAX = 0.250_dp ! convolution subroutine parameter

contains

	subroutine do_spectroscopy()
		integer :: nr_R ! which diagram to calculate: 0 or 1 (all), or 2{1-6}

		call init_spectroscopy()

		print *, "Doing Spectroscopy..."

!		call forster_weighted_fluor()

		nr_R = 0

		call time_resolved_fluor_resp(nr_R)		! calculates the response function
		call time_resolved_fluor_spec(nr_R)		! calculates the spectrum from response function

		call cleanup_spectroscopy()    ! deallocation of resources

	end subroutine do_spectroscopy



	subroutine init_spectroscopy()

		call set_values() ! initializes parameters: 'operational' and 'system'

		call set_corr_functions()

		allocate(R_trf(Nt1))
		R_trf = 0.0_dp

	end subroutine init_spectroscopy



	subroutine time_resolved_fluor_resp(nr_R)
		integer, intent(in) :: nr_R
		complex(dpc), dimension(:,:,:,:), allocatable :: AUX1, AUX2, AUX3
		real(dp), dimension(Nt1+Nt2) :: time

		integer 	:: n_conv	! limit of extent for convolution; BEWARE!
		integer 	:: n_prol	! multiplier of time extention for convolution; ALSO BEWARE!
		integer(dp) :: d_ratio	! ratio of internal time step dt to internal time step dt1 during <densification>
		complex(dpc), dimension(:,:), allocatable :: conv_in_e, conv_out_e

		integer :: e, b, c, i
		integer :: tau, tau1, tau2

		!!parameters for calculations using convolution
		n_conv = 1024
		n_prol = 1 ! must be power of 2
		if (Nt1 < n_conv) then
			print *, "ERROR: Insufficient time span: Nt1 < n_conv!"
			stop
		end if
		if (Nt2 > n_conv) then
			print *, "ERROR: Insufficient time span: Nt2 > n_conv!"
			stop
		end if
		d_ratio = 4 ! must be a power of 2
		! (Nt1;2*Nt1] FT 'sees' it as t < 0 => leave zeros there
		allocate(conv_in_e(2,2*n_conv*n_prol*d_ratio),conv_out_e(1,2*n_conv*n_prol*d_ratio+1))

		!! prepare resources/functions
		!!ATENZIONE!! allocation 1:Nt1, because Nt1>Nt2 FOR NOW!
		allocate(AUX1(Nt1+Nt2,Ns,Ns,Ns)); AUX1 = 0
		allocate(AUX2(Nt1*n_prol+Nt2,Ns,Ns,Ns)); AUX2 = 0
		allocate(AUX3(Nt1*n_prol,Ns,Ns,Ns)); AUX3 = 0

		forall(i = 1:Nt1+Nt2)
			time(i) = (i-1)*gt1*dt
		end forall

!		call init_kroneker(Ns)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		select case (nr_R)

!!!!!!!!!!!!!!!
!!! CALCULATE monomeric response function R0_trf
!!!!!!
		case(0)
		print *, "Calculating diagram #0"

		do e = 1, Ns
			do tau = 1, Nt1
				R_trf(tau) = R_trf(tau) + (dd(e)**4) * &
											exp( -(0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt - conjg(goft(e,tau-1)) &
												 +(0.0_dp,2.0_dp)*aimag(goft(e,Nt2)-goft(e,Nt2+tau)) )
			end do
		end do


!!!!!!!!!!!!!!!
!!! CALCULATE R(1_X) response functions for dimer
!!!!!!
!		case(1)
!		print *, "Calculating diagrams #1.X"

		do e = 1, Ns
			do b = 1, Ns
			if (e/=b) then
				do tau = 1, Nt1
					!! ------R1.1------
					AUX1 = 0
					do tau1 = Nt2, Nt2+tau
						AUX1(tau1,b,b,e) = exp(	  (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-Nt2)*gt1*dt	&
												- (0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
												- goft(b,Nt2+tau-tau1) - conjg(goft(e,tau1-Nt2))	)
					end do
					R_trf(tau) = R_trf(tau) - (0.0_dp,1.0_dp) * dd(b)*(dd(e)**3)*cpl(e,b)	*						&
										exp( -(0.0_dp,1.0_dp)*en(b)*(tau-1)*gt1*dt	 								&
										 	 +(0.0_dp,2.0_dp)*aimag(goft(e,Nt2-1)) )*								&
										 (				 intsplin(time(Nt2:Nt2+tau), real(AUX1(Nt2:Nt2+tau,b,b,e)))	&
										+(0.0_dp,1.0_dp)*intsplin(time(Nt2:Nt2+tau),aimag(AUX1(Nt2:Nt2+tau,b,b,e)))	)
					!! ------R1.1 coh.------
                    AUX1 = 0
                    do tau1 = Nt2, Nt2+tau
                        AUX1(tau1,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-Nt2)*gt1*dt &
                                                - goft(e,tau1-1) - conjg(goft(b,tau1-1))          &
                                                - goft(b,Nt2+tau-tau1) + conjg(goft(b,tau1-Nt2))    )
                    end do
                    R_trf(tau) = R_trf(tau) - (0.0_dp,1.0_dp) * dd(e)*(dd(b)**3)*cpl(b,e)   *                 &
                                        exp( -(0.0_dp,1.0_dp)*en(b)*(tau-1)*gt1*dt                            &
                                             -(0.0_dp,1.0_dp)*(en(e)-en(b))*(Nt2-1) - conjg(goft(b,tau-1))    &
                                             +conjg(goft(b,Nt2+tau)) - conjg(goft(b,Nt2-1)) )*                &
                                         (               intsplin(time(Nt2:Nt2+tau), real(AUX1(Nt2:Nt2+tau,b,b,e))) &
                                        +(0.0_dp,1.0_dp)*intsplin(time(Nt2:Nt2+tau),aimag(AUX1(Nt2:Nt2+tau,b,b,e))) )
					!! ------R1.2------
					AUX1 = 0
					do tau1 = 1, Nt2
						AUX1(tau1,b,b,e) = exp(	  (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt	&
												- (0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
												- goft(b,Nt2+tau-tau1-1) - goft(e,Nt2-tau1)	)
					end do
					R_trf(tau) = R_trf(tau) - (0.0_dp,1.0_dp) * dd(b)*(dd(e)**3)*cpl(e,b)	*			&
										exp( -(0.0_dp,1.0_dp)*en(b)*(tau-1)*gt1*dt						&
											 -(0.0_dp,1.0_dp)*(en(b)-en(e))*(Nt2-1)*gt1*dt				&
											 +(0.0_dp,2.0_dp)*aimag(goft(e,Nt2-1)) )*					&
										 (				 intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))	&
										+(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e)))	)
                   !! ------R1.2 coh.------
                   AUX1 = 0
                   do tau1 = 1, Nt2
                       AUX1(tau1,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt   &
                                               - conjg(goft(b,tau1-1)) + goft(b,Nt2-tau1)          &
                                               - goft(b,Nt2+tau-tau1-1) - goft(e,tau1-1) )
                   end do
                   R_trf(tau) = R_trf(tau) - (0.0_dp,1.0_dp) * dd(e)*(dd(b)**3)*cpl(b,e)   *           &
                                       exp( -(0.0_dp,1.0_dp)*en(b)*(tau-1)*gt1*dt                      &
                                            -conjg(goft(b,tau-1)) + conjg(goft(b,Nt2+tau-2))           &
                                            -conjg(goft(b,Nt2-1)) ) *                                  &
                                        (               intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e))) &
                                       +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e))) )
					!! ------R1.3------
					AUX1 = 0
					do tau1 = 1, Nt2
						AUX1(tau1,b,b,e) = exp(	 (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau1-1)*gt1*dt	&
												+(0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
												- conjg(goft(e,Nt2+tau-tau1-1)) - conjg(goft(b,Nt2-tau1))	)
					end do
					R_trf(tau) = R_trf(tau) + (0.0_dp,1.0_dp) * dd(b)*(dd(e)**3)*cpl(e,b)	*			&
										exp( -(0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt 						&
											 -(0.0_dp,1.0_dp)*(en(e)-en(b))*(Nt2-1)*gt1*dt				&
											 -(0.0_dp,2.0_dp)*aimag(goft(e,Nt2+tau-2))	)*				&
										 (				 intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))	&
										+(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e)))	)
                   !! ------R1.3 coh.------
                   AUX1 = 0
                   do tau1 = 1, Nt2
                       AUX1(tau1,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt         &
                                               - conjg(goft(b,tau1-1)) - goft(e,tau1-1)                  &
                                               - conjg(goft(e,Nt2-tau1)) + conjg(goft(e,Nt2+tau-tau1-1)) )
                   end do
                   R_trf(tau) = R_trf(tau) + (0.0_dp,1.0_dp) * dd(b)*(dd(e)**3)*cpl(b,e)   *           &
                                       exp( -(0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt                      &
                                            -conjg(goft(e,tau-1)) - goft(e,Nt2+tau-2)           &
                                            +goft(e,Nt2-1) ) *                                  &
                                        (               intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e))) &
                                       +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e))) )
				end do
			end if
			end do
		end do


!!!!!!!!!!!!!!!
!!! R(1_2)_RTF for dimer
!!!!!!
		!case(21)
		print *, "Calculating diagram #2.1"

		do e = 1, Ns
			do b = 1, Ns
			if (e/=b) then
				do tau = 1, Nt1
					AUX1 = 0; AUX2 = 0
					do tau2 = Nt2, Nt2+tau
						do tau1 = 1, Nt2
							AUX1(tau1,b,b,e) = exp(	  (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt	&
													- (0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
													- goft(e,Nt2-tau1) + goft(e,Nt2+tau-tau1-1)			&
													- goft(e,tau2-tau1)-goft(b,tau2-tau1)	)
						end do
						AUX2(tau2,b,b,e) = exp(	 (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt		&
											    +(0.0_dp,2.0_dp) * aimag(goft(e,tau2-1))				&
												- goft(e,Nt2+tau-tau2) + conjg(goft(e,tau2-Nt2)) ) *	&
									 (				 intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))	&
									+(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e)))	)
					end do
					R_trf(tau) = R_trf(tau) - (dd(e)**4)*(cpl(e,b)**2)	*									&
										exp( - (0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt - conjg(goft(e,tau-1))	&
											 + (0.0_dp,2.0_dp) * aimag(goft(e,Nt2-1)-goft(e,Nt2+tau-2))	)*	&
											 ( intsplin(time(Nt2:Nt2+tau), real(AUX2(Nt2:Nt2+tau,b,b,e)))	&
							  +(0.0_dp,1.0_dp)*intsplin(time(Nt2:Nt2+tau),aimag(AUX2(Nt2:Nt2+tau,b,b,e)))	)
				end do
			end if
			end do
		end do


!!!!!!!!!!!!!!!
!!! R(1_2)_RTF for dimer COHERENT I.C.s
!!!!!!
        !case(21c)
        print *, "Calculating diagram #2.1 coherent i.c."

        do e = 1, Ns
            do b = 1, Ns
            if (e/=b) then
                do tau = 1, Nt1
                    AUX1 = 0; AUX2 = 0
                    do tau2 = Nt2, Nt2+tau
                        do tau1 = 1, Nt2
                            AUX1(tau1,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt   &
                                                    - goft(e,tau1-1) - conjg(goft(b,tau1-1))            &
                                                    + goft(b,Nt2-tau1)  + goft(e,Nt2+tau-tau1-1)        &
                                                    - goft(e,tau2-tau1) - goft(b,tau2-tau1)   )
                        end do
                        AUX2(tau2,b,b,e) = exp(  (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt        &
                                                + goft(e,tau2-1) + conjg(goft(b,tau2-1))                &
                                                - goft(e,Nt2+tau-tau2-1) - conjg(goft(b,tau2-Nt2)) ) *  &
                                     (               intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))     &
                                    +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e))) )
                    end do
                    R_trf(tau) = R_trf(tau) - (dd(b)**2)*(dd(e)**2)*(cpl(e,b)**2)  *            &
                                        exp( - (0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt             &
                                             - conjg(goft(b,Nt2-1)) - goft(e,Nt2+tau-2)         &
                                             - (0.0_dp,1.0_dp)*(en(e)-en(b))*(Nt2-1)*gt1*dt )*  &
                                             ( intsplin(time(Nt2:Nt2+tau), real(AUX2(Nt2:Nt2+tau,b,b,e)))   &
                              +(0.0_dp,1.0_dp)*intsplin(time(Nt2:Nt2+tau),aimag(AUX2(Nt2:Nt2+tau,b,b,e)))   )
                end do
            end if
            end do
        end do


!!!!!!!!!!!!!!!
!!! R(2_2)_RTF for dimer
!!!!!!
		!case(22)
		print *, "Calculating diagram #2.2"

		do e = 1, Ns
			do b = 1, Ns
			if (e/=b) then
				do tau = 1, Nt1
					AUX1 = 0; AUX2 = 0
					do tau2 = Nt2, Nt2+tau
						do tau1 = 1, Nt2
							AUX1(tau1,b,b,e) = exp(	- (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau1-1)*gt1*dt	&
													- (0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
													- goft(b,Nt2-tau1) + goft(b,Nt2+tau-tau1-1)			&
													- goft(e,tau2-tau1)- goft(b,tau2-tau1)	)
							AUX1(tau1,b,b,e) = conjg(AUX1(tau1,b,b,e))
						end do
						AUX2(tau2,b,b,e) = exp(	 (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau2-1)*gt1*dt	&
											    -(0.0_dp,2.0_dp) * aimag(goft(e,tau2-1))			&
												- goft(b,Nt2+tau-tau2) - goft(b,tau2-Nt2) ) *		&
									 (				 intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))	&
									+(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e)))	)
					end do
					R_trf(tau) = R_trf(tau) + (dd(e)**2)*(dd(b)**2)*(cpl(e,b)**2)	*						&
								exp( - (0.0_dp,1.0_dp)*en(b)*(tau-1)*gt1*dt - conjg(goft(b,tau-1))	)	*	&
											 ( intsplin(time(Nt2:Nt2+tau), real(AUX2(Nt2:Nt2+tau,b,b,e)))	&
							  +(0.0_dp,1.0_dp)*intsplin(time(Nt2:Nt2+tau),aimag(AUX2(Nt2:Nt2+tau,b,b,e)))	)
				end do
			end if
			end do
		end do


!!!!!!!!!!!!!!!
!!! R(2_2)_RTF for dimer COHERENT I.C.s
!!!!!!
        !case(22c)
        print *, "Calculating diagram #2.2 coherent i.c."

        do e = 1, Ns
            do b = 1, Ns
            if (e/=b) then
                do tau = 1, Nt1
                    AUX1 = 0; AUX2 = 0
                    do tau2 = Nt2, Nt2+tau
                        do tau1 = 1, Nt2
                            AUX1(tau1,b,b,e) = exp( - (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt   &
                                                    - goft(b,tau1-1) - conjg(goft(e,tau1-1))            &
                                                    - goft(e,Nt2-tau1) - goft(b,Nt2+tau-tau1-1)         &
                                                    + goft(e,tau2-tau1) + goft(b,tau2-tau1) )
                            AUX1(tau1,b,b,e) = conjg(AUX1(tau1,b,b,e))
                        end do
                        AUX2(tau2,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau2-1)*gt1*dt       &
                                                - goft(e,tau2-1) - conjg(goft(b,tau2-1))                &
                                                - conjg(goft(e,tau2-Nt2)) - goft(b,Nt2+tau-tau2-1) ) *  &
                                     (               intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))     &
                                    +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e)))     )
                    end do
                    R_trf(tau) = R_trf(tau) + (dd(e)**2)*(dd(b)**2)*(cpl(e,b)**2)   *           &
                                        exp( - (0.0_dp,1.0_dp)*en(b)*(tau-1)*gt1*dt             &
                                             + conjg(goft(b,Nt2+tau-2)) + goft(e,Nt2-1)         &
                                             - (0.0_dp,1.0_dp)*(en(b)-en(e))*(Nt2-1)*gt1*dt ) * &
                                          (               intsplin(time(1:Nt2), real(AUX2(1:Nt2,b,b,e))) &
                                         +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX2(1:Nt2,b,b,e))) )
                end do
            end if
            end do
        end do


!!!!!!!!!!!!!!!
!!! R(3_2)_RTF for dimer; there are no coherence transfer terms
!!!!!!
		!case(23)
		print *, "Calculating diagram #2.3"

		do e = 1, Ns
			do b = 1, Ns
			if (e/=b) then
				do tau = 1, Nt1
				    AUX1 = 0; AUX2 = 0
					do tau1 = 1, Nt2
						do tau2 = 1, Nt2
							AUX2(tau2,b,b,e) = exp(	- (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt	&
													- (0.0_dp,2.0_dp) * aimag(goft(e,tau2-1))			&
													- goft(b,Nt2-tau2) + goft(b,Nt2+tau-tau2)			&
													- goft(e,tau1-tau2)	- goft(b,tau1-tau2)	)
							AUX2(tau2,b,b,e) = conjg(AUX2(tau2,b,b,e))
						end do
						AUX1(tau1,b,b,e) = exp(	 (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt	&
											    -(0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
												+ goft(b,Nt2-tau1) - goft(b,Nt2+tau-tau1) ) *		&
									 (				 intsplin(time(1:Nt2), real(AUX2(1:Nt2,b,b,e)))	&
									+(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX2(1:Nt2,b,b,e)))	)
					end do
					R_trf(tau) = R_trf(tau) + (dd(e)**2)*(dd(b)**2)*(cpl(e,b)**2)	*							&
										exp( - (0.0_dp,1.0_dp)*en(b)*(tau-1)*gt1*dt - conjg(goft(b,tau-1)) ) *  &
											 (				 intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))		&
											+(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e)))	    )
				end do
			end if
			end do
		end do

 !(the full version of R(3_2))
!		do e = 1, Ns
!			do b = 1, Ns
!			do c = 1, Ns
!			if ((e/=b).and.(e/=c)) then
!				do tau = 1, Nt1
!					do tau1 = 1, Nt2
!						do tau2 = 1, Nt2
!							AUX2(tau2,c,b,e) = exp(	- (0.0_dp,1.0_dp)*(en(e)-en(c))*(tau2-1)*gt1*dt	&
!													- (0.0_dp,2.0_dp) * aimag(goft(e,tau2))			&
!												 	- goft(c,Nt2-tau2)								&
!													+ kroneker(c,b)*goft(c,Nt2+Nt1-tau2)			&
!													- goft(e,tau1-tau2)								&
!													- kroneker(c,b)*goft(b,tau1-tau2)	)
!							AUX2(tau2,c,b,e) = conjg(AUX2(tau2,c,b,e))
!						end do
!						AUX1(tau1,c,b,e) = exp(	(0.0_dp,1.0_dp)*(en(b)-en(e))*gt1*(tau1-1)*dt	&
!												- (0.0_dp,2.0_dp) * aimag(goft(e,tau1))			&
!												+ kroneker(c,b)*goft(b,Nt2-tau1)				&
!												- goft(b,Nt2+Nt1-tau1)	) *	&
!										(				 intsplin(time, real(AUX2(1:Nt2,c,b,e)))	&
!										+(0.0_dp,1.0_dp)*intsplin(time,aimag(AUX2(1:Nt2,c,b,e)))	)
!					end do
!					R_trf(tau) = R_trf(tau) + (dd(e)**2)*dd(b)*dd(c)*cpl(b,e)*cpl(e,c)				&
!											* exp( -(0.0_dp,1.0_dp)*( (en(b)-en(c))*Nt2+en(b)*tau )	&
!												   -kroneker(b,c)*conjg(goft(b,tau)) )	*			&
!											(	intsplin(time, real(AUX1(1:Nt2,c,b,e)))				&
!											+ (0.0_dp,1.0_dp)*intsplin(time,aimag(AUX1(1:Nt2,c,b,e)))	)
!				end do
!			end if
!			end do
!			end do
!		end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TESTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	! A = CONVOLUTION
!	conv_in_e = 0.0_dp; conv_out_e = 0.0_dp
!	do tau1 = 1, n_conv
!		conv_in_e(1,tau1) = exp((0.0_dp,1.0_dp) * (en(2)-en(1))*(tau1-1)*gt1*dt)
!		conv_in_e(2,tau1) = exp(- goft(2,tau1-1) - goft(1,tau1-1))
!	end do
!	call fft_row_convolution(conv_out_e(1:1,2:2*n_conv+1),&
!  							 conv_in_e(1:1,1:2*n_conv),	&
!  							 conv_in_e(2:2,1:2*n_conv), &
!							 REGULARIZATION_TIME_PER_TMAX)
!open(unit=11,file=( trim(out_dir) // "/" // trim("convolution.dat") ), position='append')
!        do i = 1, 2*n_conv
!        	write(11, '(3f16.8)') (i-1)*gt1*dt, real(conv_out_e(1,i)), aimag(conv_out_e(1,i))
!        end do
!close(unit=11)

	! B = INTEGRATION
!	do tau2 = 1, n_conv
!		do tau1 = 1, tau2
!		AUX1(tau1,1,1,1) = exp((0.0_dp,1.0_dp)*(en(2)-en(1))*(tau1-1)*gt1*dt - goft(2,tau2-tau1) - goft(1,tau2-tau1))
!		end do
!	AUX2(tau2,1,1,1) = ( intsplin(time(1:tau2), real(AUX1(1:tau2,1,1,1)))	&
!		+(0.0_dp,1.0_dp)*intsplin(time(1:tau2),aimag(AUX1(1:tau2,1,1,1)))	)
!	end do
!open(unit=11,file=( trim(out_dir) // "/" // trim("integral.dat") ), position='append')
!        do i = 1, Nt1
!        	write(11, '(3f16.8)') (i-1)*gt1*dt, real(AUX2(i,1,1,1)), aimag(AUX2(i,1,1,1))
!        end do
!close(unit=11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!
!!! R(3_2)_RTF for dimer COHERENT I.C.s
!!!!!!
        !case(23c)
        print *, "Calculating diagram #2.3 coherent i.c."

        do e = 1, Ns
            do b = 1, Ns
            if (e/=b) then
                do tau = 1, Nt1
                    AUX1 = 0; AUX2 = 0
                    do tau2 = 1, Nt2
                        do tau1 = 1, Nt2
                            AUX1(tau1,b,b,e) = exp( - (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt   &
                                                    - goft(b,tau1-1) - conjg(goft(e,tau1-1))            &
                                                    - goft(e,Nt2-tau1) - goft(b,Nt2+tau-tau1-1)           &
                                                    + goft(e,tau2-tau1) + goft(b,tau2-tau1) )
                            AUX1(tau1,b,b,e) = conjg(AUX1(tau1,b,b,e))
                        end do
                        AUX2(tau2,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau2-1)*gt1*dt   &
                                                - goft(e,tau2-1) - conjg(goft(b,tau2-1))            &
                                                - goft(e,Nt2-tau2) - goft(b,Nt2+tau-tau2-1) ) *       &
                                     (               intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e))) &
                                    +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e))) )
                    end do
                    R_trf(tau) = R_trf(tau) + (dd(e)**2)*(dd(b)**2)*(cpl(e,b)**2)   *           &
                                        exp( - (0.0_dp,1.0_dp)*en(b)*(tau-1)*gt1*dt             &
                                             + conjg(goft(b,Nt2+tau-2)) + goft(e,Nt2-1)         &
                                             - (0.0_dp,1.0_dp)*(en(b)-en(e))*(Nt2-1)*gt1*dt ) * &
                                          (               intsplin(time(1:Nt2), real(AUX2(1:Nt2,b,b,e))) &
                                         +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX2(1:Nt2,b,b,e))) )
                end do
            end if
            end do
        end do


!!!!!!!!!!!!!!!
!!! R(4_2)_RTF for dimer
!!!!!!
		!case(24)
		print *, "Calculating diagram #2.4"

		do e = 1, Ns
			do b = 1, Ns
			if (e/=b) then
				do tau = 1, Nt1

					conv_in_e = 0.0_dp; conv_out_e = 0.0_dp
					do tau1 = 1, n_conv
						AUX2(tau1,b,b,e) = exp((0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt	&
												-(0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))		&
												+ goft(e,Nt2+tau-tau1-1) - goft(e,Nt2-tau1))
						AUX3(tau1,b,b,e) = exp(- goft(e,tau1-1) - goft(b,tau1-1))
					end do
					! AUX evolution is longer than n_conv: 0's are 'padded' at tau1=(n_prol-1)*n_conv
					call densify_linear(AUX2(1:n_conv*n_prol,b,b,e),conv_in_e(1,1:n_conv*n_prol*d_ratio))
					call densify_linear(AUX3(1:n_conv*n_prol,b,b,e),conv_in_e(2,1:n_conv*n_prol*d_ratio))
					call fft_row_convolution(conv_out_e(1:1,2:2*n_conv*n_prol*d_ratio+1),&
				  							 conv_in_e(1:1,1:2*n_conv*n_prol*d_ratio),	&
				  							 conv_in_e(2:2,1:2*n_conv*n_prol*d_ratio), &
											 REGULARIZATION_TIME_PER_TMAX)

					AUX1 = 0; AUX2 = 0
					do tau2 = 1, Nt2
!						do tau1 = 1, tau2
!							AUX2(tau1,b,b,e) = exp(	  (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt	&
!													- (0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
!													- goft(e,Nt2-tau1) + goft(e,Nt2+tau-tau1-1) 		&
!													- goft(e,tau2-tau1) - goft(b,tau2-tau1)	)
!						end do
						AUX1(tau2,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt	&
												+ (0.0_dp,2.0_dp) * aimag(goft(e,tau2-1))			&
												+ (goft(e,Nt2-tau2)) - goft(e,Nt2+tau-tau2-1) )	*	&
											conv_out_e(1,(tau2-1)*d_ratio+1)*(gt1*dt/d_ratio)!comment out if 'directly'
!										 (				 intsplin(time(1:tau2), real(AUX2(1:tau2,b,b,e)))	&
!										+(0.0_dp,1.0_dp)*intsplin(time(1:tau2),aimag(AUX2(1:tau2,b,b,e)))	)
					end do
					R_trf(tau) = R_trf(tau) - (dd(e)**4)*(cpl(e,b)**2)	*	&
										exp( -(0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt - conjg(goft(e,tau-1))	&
											 +(0.0_dp,2.0_dp) * aimag(goft(e,Nt2-1)-goft(e,Nt2+tau-2)) )*	&
										 (				 intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))		&
										+(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e)))	)
				end do
			end if
			end do
		end do


!!!!!!!!!!!!!!!
!!! R(4_2)_RTF for dimer COHERENT I.C.s
!!!!!!
        !case(24c)
        print *, "Calculating diagram #2.4 coherent i.c."

        do e = 1, Ns
            do b = 1, Ns
            if (e/=b) then
                do tau = 1, Nt1
                    AUX1 = 0; AUX2 = 0
                    do tau2 = 1, Nt2
                        do tau1 = 1, tau2
                            AUX1(tau1,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt   &
                                                    - goft(e,tau1-1) - conjg(goft(b,tau1-1))            &
                                                    + goft(b,Nt2-tau1)  + goft(e,Nt2+tau-tau1-1)        &
                                                    - goft(e,tau2-tau1) - goft(b,tau2-tau1)   )
                        end do
                        AUX2(tau2,b,b,e) = exp(  (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt        &
                                                + goft(e,tau2-1) + conjg(goft(b,tau2-1))                &
                                                - goft(e,Nt2+tau-tau2-1) - goft(b,Nt2-tau2) ) *         &
                                     (               intsplin(time(1:tau2), real(AUX1(1:tau2,b,b,e)))   &
                                    +(0.0_dp,1.0_dp)*intsplin(time(1:tau2),aimag(AUX1(1:tau2,b,b,e))) )
                    end do
                    R_trf(tau) = R_trf(tau) - (dd(b)**2)*(dd(e)**2)*(cpl(e,b)**2)  *            &
                                        exp( - (0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt             &
                                             - conjg(goft(b,Nt2-1)) - goft(e,Nt2+tau-2)         &
                                             - (0.0_dp,1.0_dp)*(en(e)-en(b))*(Nt2-1)*gt1*dt ) * &
                                             ( intsplin(time(1:Nt2), real(AUX2(1:Nt2,b,b,e)))   &
                              +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX2(1:Nt2,b,b,e)))   )
                end do
            end if
            end do
        end do


!!!!!!!!!!!!!!!
!!! R(5_2)_RTF for dimer; there are no coherence transfer terms
!!!!!!
		!case(25)
		print *, "Calculating diagram #2.5"

		do e = 1, Ns
			do b = 1, Ns
			if (e/=b) then
				do tau = 1, Nt1

!					conv_in_e = 0.0_dp; conv_out_e = 0.0_dp
!					do tau1 = 1, n_conv
!						AUX2(tau1,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt	&
!												- (0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
!												- goft(e,Nt2+tau-tau1-1) + goft(e,Nt2-tau1)		)
!						AUX3(tau1,b,b,e) = exp(	- goft(e,tau1-1) - goft(b,tau1-1)	)
!						AUX2(tau1,b,b,e) = conjg(AUX2(tau1,b,b,e))
!						AUX3(tau1,b,b,e) = conjg(AUX3(tau1,b,b,e))
!					end do
!					! AUX evolution is longer than n_conv: 0's are 'padded' at tau1=(n_prol-1)*n_conv
!					call densify_linear(AUX2(1:n_conv*n_prol,b,b,e),conv_in_e(1,1:n_conv*n_prol*d_ratio))
!					call densify_linear(AUX3(1:n_conv*n_prol,b,b,e),conv_in_e(2,1:n_conv*n_prol*d_ratio))
!					call fft_row_convolution(conv_out_e(1:1,2:2*n_conv*n_prol*d_ratio+1),&
!				  							 conv_in_e(1:1,1:2*n_conv*n_prol*d_ratio),	&
!				  							 conv_in_e(2:2,1:2*n_conv*n_prol*d_ratio), &
!											 REGULARIZATION_TIME_PER_TMAX)

					AUX1 = 0; AUX2 = 0
					do tau2 = 1, Nt2
						do tau1 = 1, tau2
							AUX2(tau1,b,b,e) = exp(	  (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt	&
													- (0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
													+ goft(e,Nt2-tau1) - goft(e,Nt2+tau-tau1-1) 		&
													- goft(e,tau2-tau1) - goft(b,tau2-tau1)	)
							AUX2(tau1,b,b,e) = conjg(AUX2(tau1,b,b,e))
						end do
						AUX1(tau2,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt	&
												+ (0.0_dp,2.0_dp) * aimag(goft(e,tau2-1))			&
												- (goft(e,Nt2-tau2)) + goft(e,Nt2+tau-tau2-1) )
						AUX1(tau2,b,b,e) = conjg(AUX1(tau2,b,b,e)) *	&
!											conv_out_e(1,(tau2-1)*d_ratio+1)*(gt1*dt/d_ratio)!comment out if 'directly'
										 (				 intsplin(time(1:tau2), real(AUX2(1:tau2,b,b,e)))	&
										+(0.0_dp,1.0_dp)*intsplin(time(1:tau2),aimag(AUX2(1:tau2,b,b,e)))	)
					end do
					R_trf(tau) = R_trf(tau) - (dd(e)**4)*(cpl(e,b)**2)	*	&
										exp( -(0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt - conjg(goft(e,tau-1)) +	&
											(0.0_dp,2.0_dp) * aimag(goft(e,Nt2-1)-goft(e,Nt2+tau-2)) )	*	&
										 (				 intsplin(time(1:Nt2), real(AUX1(1:Nt2,b,b,e)))		&
										+(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX1(1:Nt2,b,b,e)))	)
				end do
			end if
			end do
		end do


!!!!!!!!!!!!!!!
!!! R(5_2)_RTF for dimer COHERENT I.C.s
!!!!!!
        !case(25c)
        print *, "Calculating diagram #2.5 coherent i.c."

        do e = 1, Ns
            do b = 1, Ns
            if (e/=b) then
                do tau = 1, Nt1
                    AUX1 = 0; AUX2 = 0
                    do tau2 = 1, Nt2
                        do tau1 = 1, tau2
                            AUX1(tau1,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt          &
                                                    - goft(e,tau1-1) - conjg(goft(b,tau1-1))                   &
                                                    + conjg(goft(b,Nt2-tau1))  + conjg(goft(e,Nt2+tau-tau1-1)) &
                                                    - conjg(goft(e,tau2-tau1)) - conjg(goft(b,tau2-tau1))   )
                        end do
                        AUX2(tau2,b,b,e) = exp(  (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt               &
                                                + goft(e,tau2-1) + conjg(goft(b,tau2-1))                       &
                                                - conjg(goft(e,Nt2+tau-tau2-1)) - conjg(goft(b,Nt2-tau2)) ) *  &
                                     (               intsplin(time(1:tau2), real(AUX1(1:tau2,b,b,e)))          &
                                    +(0.0_dp,1.0_dp)*intsplin(time(1:tau2),aimag(AUX1(1:tau2,b,b,e))) )
                    end do
                    R_trf(tau) = R_trf(tau) - (dd(b)**2)*(dd(e)**2)*(cpl(e,b)**2)  *            &
                                        exp( - (0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt             &
                                             - conjg(goft(b,Nt2-1)) - goft(e,Nt2+tau-2)         &
                                             - (0.0_dp,1.0_dp)*(en(e)-en(b))*(Nt2-1)*gt1*dt )*  &
                                             ( intsplin(time(1:Nt2), real(AUX2(1:Nt2,b,b,e)))   &
                              +(0.0_dp,1.0_dp)*intsplin(time(1:Nt2),aimag(AUX2(1:Nt2,b,b,e)))   )
                end do
            end if
            end do
        end do


!!!!!!!!!!!!!!!
!!! R(6_2)_RTF for dimer; there are no coherence transfer terms
!!!!!!
		!case(26)
		print *, "Calculating diagram #2.6"

		do e = 1, Ns
			do b = 1, Ns
			if (e/=b) then
				do tau = 1, Nt1
				    AUX1 = 0; AUX2 = 0
					do tau2 = Nt2, Nt2+tau
						do tau1 = Nt2, tau2
							AUX1(tau1,b,b,e) = exp(	  (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt	&
													- (0.0_dp,2.0_dp) * aimag(goft(e,tau1-1))			&
													- conjg(goft(e,tau1-1)) + goft(e,Nt2+tau-tau1)			&
													- goft(e,tau2-tau1)-goft(b,tau2-tau1)	)
						end do
						AUX2(tau2,b,b,e) = exp(	 (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt	&
											    +(0.0_dp,2.0_dp) * aimag(goft(e,tau2-1))		&
												- goft(e,Nt2+tau-tau2) + conjg(goft(e,tau2-Nt2)) ) *		&
									 (				 intsplin(time(Nt2:tau2), real(AUX1(Nt2:tau2,b,b,e)))	&
									+(0.0_dp,1.0_dp)*intsplin(time(Nt2:tau2),aimag(AUX1(Nt2:tau2,b,b,e)))	)
					end do
					R_trf(tau) = R_trf(tau) - (dd(e)**4)*(cpl(e,b)**2)	*						&
										exp( - (0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt - conjg(goft(e,tau-1))	&
											 + (0.0_dp,2.0_dp) * aimag(goft(e,Nt2-1)-goft(e,Nt2+tau-2)) )*	&
											 ( intsplin(time(Nt2:Nt2+tau), real(AUX2(Nt2:Nt2+tau,b,b,e)))	&
							  +(0.0_dp,1.0_dp)*intsplin(time(Nt2:Nt2+tau),aimag(AUX2(Nt2:Nt2+tau,b,b,e)))	)
				end do
			end if
			end do
		end do


!!!!!!!!!!!!!!!
!!! R(6_2)_RTF for dimer COHERENT I.C.s
!!!!!!
        !case(26c)
        print *, "Calculating diagram #2.6 coherent i.c."

        do e = 1, Ns
            do b = 1, Ns
            if (e/=b) then
                do tau = 1, Nt1
                    AUX1 = 0; AUX2 = 0
                    do tau2 = Nt2, Nt2+tau
                        do tau1 = Nt2, tau2
                            AUX1(tau1,b,b,e) = exp(   (0.0_dp,1.0_dp) * (en(b)-en(e))*(tau1-1)*gt1*dt   &
                                                    - goft(e,tau1-1) - conjg(goft(b,tau1-1))            &
                                                    + conjg(goft(b,tau1-Nt2))  + goft(e,Nt2+tau-tau1-1)        &
                                                    - goft(e,tau2-tau1) - goft(b,tau2-tau1)   )
                        end do
                        AUX2(tau2,b,b,e) = exp(  (0.0_dp,1.0_dp) * (en(e)-en(b))*(tau2-1)*gt1*dt        &
                                                + goft(e,tau2-1) + conjg(goft(b,tau2-1))                &
                                                - goft(e,Nt2+tau-tau2-1) - conjg(goft(b,tau2-Nt2)) ) *  &
                                     (               intsplin(time(Nt2:tau2), real(AUX1(Nt2:tau2,b,b,e)))     &
                                    +(0.0_dp,1.0_dp)*intsplin(time(Nt2:tau2),aimag(AUX1(Nt2:tau2,b,b,e))) )
                    end do
                    R_trf(tau) = R_trf(tau) - (dd(b)**2)*(dd(e)**2)*(cpl(e,b)**2)  *            &
                                        exp( - (0.0_dp,1.0_dp)*en(e)*(tau-1)*gt1*dt             &
                                             - conjg(goft(b,Nt2-1)) - goft(e,Nt2+tau-2)         &
                                             - (0.0_dp,1.0_dp)*(en(e)-en(b))*(Nt2-1)*gt1*dt )*  &
                                             ( intsplin(time(Nt2:Nt2+tau), real(AUX2(Nt2:Nt2+tau,b,b,e)))   &
                              +(0.0_dp,1.0_dp)*intsplin(time(Nt2:Nt2+tau),aimag(AUX2(Nt2:Nt2+tau,b,b,e)))   )
                end do
            end if
            end do
        end do


		end select
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!		R_trf(1) = R_trf(1) - (0.0_dp,1.0_dp)*aimag(R_trf(1))
!		open(unit=11,file=( trim(out_dir) // "/" // trim("R_TRF.dat") ), position='append')
!		do i = 1, Nt1
!			write(11, '(3f16.8)') (i-1)*gt1*dt, real(R_TRF(i)), aimag(R_TRF(i))
!		end do
!		close(unit=11)

		!! clean subroutine
		deallocate(AUX1, AUX2, AUX3, conv_in_e, conv_out_e)

!		call delete_kroneker()

	end subroutine time_resolved_fluor_resp



	subroutine time_resolved_fluor_spec(nr_R)
		integer, intent(in) :: nr_R
		complex(dpc), dimension(:),   allocatable    :: sig
        complex(dpc), dimension(:,:), allocatable    :: dat

        integer		:: i, padfac, NFFT
        real(dp)	:: oma, dom
		character(len = 20)	:: pop_time
		character(len=4)	:: method_nr
		character(len=50)	:: name

		padfac = 7 !4 <=> resolution 2cm^-1
		NFFT = (2**padfac)*NFFT_basic

		allocate(sig(NFFT))
		allocate(dat(1,NFFT))
		if(.not. allocated(F_trf)) then
			allocate(F_trf(NFFT))
		end if
        F_trf = 0.0_dp

        dom = 2.0_dp*PI_D/(NFFT*gt1*dt)

        sig = 0.0_dp
        sig(1:Nt1) = R_trf(1:Nt1)

        do i = 1, NFFT/2
            sig(NFFT/2+i) = conjg(sig(NFFT/2-i+2))
        end do

        !
        ! FFT a la numerical recipes
        !
        dat(1,:) = sig(:)
        call fft_row(dat,1_i4b)
        sig(1:NFFT/2-1)  = dat(1,NFFT/2+2:NFFT)
        sig(NFFT/2:NFFT) = dat(1,1:NFFT/2+1)

		F_trf = real(sig)

		! output spectra
		write(pop_time, '(i4)') Nt2
		write(method_nr, '(i4)') nr_R
		name = 'TRF_' // trim(method_nr) // '_' // trim(pop_time) // 'fs.dat'
		open(unit=12,file=( trim(out_dir) // "/" // trim(name) ), position='append')
        do i = 1, NFFT
        	oma = (-NFFT/2 + i)*dom * Energy_internal_to_cm + rwa
        	write(12, '(2f16.8)') oma, F_trf(i)
        end do
        close(unit=12)

		deallocate(sig)
		deallocate(dat)
		deallocate(F_trf)

	end subroutine time_resolved_fluor_spec



	subroutine forster_weighted_fluor()
		real(dp), dimension(Ns,Ns) :: f_rates
		real(dp), dimension(Ns) :: population
		real(dp), dimension(Ns) :: pop_0		! initial conditions
		integer :: tt, i
		! for a dimer
		real(dp) :: f, b, r ! forward (b <- a), backward (a <- b) rates; r = f + b, the decay rate

		pop_0(1) = dd(1)**2/(dd(1)**2+dd(2)**2); pop_0(2) = dd(2)**2/(dd(1)**2+dd(2)**2)

		call forster_rates(f_rates)
		f = f_rates(2,1); b = f_rates(1,2); r = f + b

		open(unit=12,file=( trim(out_dir) // "/" // "forster_ev_dimer.dat" ), position='append')! evolution
		do tt = 1, Nt1, 1 ! output step																		! evolution

!		tt = Nt2	! coment out for evolution
			population(1) = (pop_0(1)/r)   * (b + f*exp(-r*(tt-1)*gt1*dt)) +	&
							(pop_0(2)*b/r) * (1 -   exp(-r*(tt-1)*gt1*dt))
			population(2) = (pop_0(1)*f/r) * (1 -   exp(-r*(tt-1)*gt1*dt)) +	&
							(pop_0(2)/r)   * (f + b*exp(-r*(tt-1)*gt1*dt))
        	write(12, '(3f16.8)') (tt-1)*gt1*dt, population(1), population(2)

		end do																					! evolution
		close(unit=12)																			! evolution

		do tt = 1, Nt1
			do i = 1, Ns
				R_trf(tt) = R_trf(tt) + &
							population(i)*exp(-(0.0_dp,1.0_dp)*(en(i)-2*re_en(i))*(tt-1)*gt1*dt-conjg(goft(i,tt-1)))
			end do
		end do

	end subroutine forster_weighted_fluor



	subroutine forster_rates(rate)
		real(dp), dimension(:,:), intent(out) :: rate	! f_r(a,b) = rate b -> a !!
		real(dp)                              :: integ	!	overlap(i,j)
		real(dp),     dimension(Nt1)          :: tt
		real(dp),     dimension(1)            :: w
		complex(dpc), dimension(Nt1)          :: yy
		complex(dpc), dimension(Nt1,1)        :: pyr,pyi
		integer :: i,j,t

		rate = 0.0_dp
		forall(t = 1:Nt1)
			tt(t) = (t-1)*gt1*dt
		end forall

		do i = 1, Ns
			do j = i+1, Ns ! (i+1, Ns) calc. explicitly only the "upper" triangle

				w(1) = -(en(i)-en(j))
				do t = 1, Nt1
                	yy(t) = exp(	- goft(i,t-1) - goft(j,t-1) - (0.0_dp,1.0_dp)*(re_en(i)+re_en(j))*tt(t)	)
                end do

                call primitive(tt, real(yy),w,pyr)
                call primitive(tt,aimag(yy),w,pyi)
                integ = 2.0_dp * real(	pyr(Nt1,1) + (0.0_dp,1.0_dp)*pyi(Nt1,1)	)

                ! correction for numerical or other errors (negative rates)
                if (integ > 0.0_dp) then
                	rate(i,j) =  (cpl(i,j)**2) * integ
                 else
                	rate(i,j) = 0.0_dp
                end if
				rate(j,i) = rate(i,j) * exp(beta_int * (en(i)-en(j))) ! "lower" triangle by detailed balance

			end do
			rate(i,i) = -sum(rate(:,i))	! "diagonal" by summing the "outgoing" rates
		end do


	end subroutine forster_rates



	subroutine cleanup_spectroscopy()

		print *, "Cleaning up Spectroscopy..."

		call cleanup_platform()    ! deallocation of resources

		if (allocated(R_trf)) then
			deallocate(R_trf)
		end if

	end subroutine cleanup_spectroscopy



	!! COMPLEMENTARY FUNCTIONS/SUBROUTINES

	!
	! Densifying a function by d_ratio using simple linear interpolation
	!
	subroutine densify_linear(in,out)
 		complex(dpc), dimension(:), intent(in)		:: in
 		complex(dpc), dimension(:), intent(out)		:: out
 		integer :: i, j, ratio, size_in, size_out

 		size_in  = size(in,1)
 		size_out = size(out,1)
 		ratio = size_out/size_in

 		do i = 1, size_in
 		  do j = 1, ratio
 		  	if (i/=size_in) then
 		      out((i-1)*ratio+j) = in(i) + (in(i+1)-in(i))*(j-1)/ratio
 		      else
 		      out((i-1)*ratio+j) = in(i)
 		    end if
 		  end do
 		end do

	end subroutine



end module
