module HEOM

	use std_types
	use prep_rules
	use resources


	implicit none


contains

	subroutine do_HEOM()

		call initialize_HEOM() ! activating parameters, etc.
		call propagate_HEOM()  ! the actual calculation
		call cleanup_platform()    ! deallocation of resources

	end subroutine do_HEOM



	subroutine initialize_HEOM()

		call set_values() ! initializes parameters: 'operational' and 'system'

		!! if not created during previous runs for the system with the same (N; K; L)
		call create_index()

		call set_corr_functions()
        call create_coefficients()

	end subroutine initialize_HEOM



	subroutine propagate_HEOM()
		integer :: i, a, b
		real(dp), dimension(2*size(ADO,2)*size(ADO,2)) :: RDO_row
		character(len = 20) 	:: buffer

		print *, "Doing HEOM..."
		print *, "The number of ADO's is ", NN

		write(buffer,'(i3)') 2*Ns*Ns+1
		buffer = '(' // trim(buffer) // 'f16.8)' ! forming a descriptor for outputting row by row
		open(unit=12,file=(trim(out_dir)//"/"//trim("RDO_s.dat")), position='append')

		!!! set initial conditions
		ADO(0,1,1) = 1.0_dp

		do a = 1, Ns
			do b = 1, Ns
				RDO_row((a-1)*2*Ns+(b-1)*2 + 1) =  real(ADO(0,a,b))
				RDO_row((a-1)*2*Ns+(b-1)*2 + 2) = aimag(ADO(0,a,b))
			end do
		end do

		write(12,buffer) 0.0_dp, RDO_row(1:2*Ns*Ns) ! the first line is written down

		!!! run calculation
		do i = 1, Nt1
		    print *, 'Runing time step #', i, ' out of ', Nt1
			call make_step_HEOM(ADO,-1)
			do a = 1, Ns
				do b = 1, Ns
					RDO_row((a-1)*2*Ns+(b-1)*2 + 1) =  real(ADO(0,a,b))
					RDO_row((a-1)*2*Ns+(b-1)*2 + 2) = aimag(ADO(0,a,b))
				end do
			end do
			write(12,buffer) i*gt1*dt, RDO_row(1:2*Ns*Ns)
!			if (abs(ADO(0,1,1)) > 1) then
!				print *, "+++ Terminating. Negative values, starting at time step ", i
!				exit
!			end if
		end do

		close(unit=12)

	end subroutine propagate_HEOM



	subroutine make_step_HEOM(y,lower)
		integer, intent(in) :: lower ! a Fortran pecularity: passing over a special boundary index; i.e. -1 for ADO here
		complex(dpc), intent(inout), dimension(lower:,:,:)  :: y
		integer :: i
		!! rk4 parameters
		real(dp) :: hh, h6
		complex(dpc), dimension(lower:(size(y,1)-2),size(y,2),size(y,3)) :: yt, dy, dyt, dym

		hh = dt * 0.5_dp
		h6 = dt / 6.0_dp
        yt = 0; dy = 0; dyt = 0; dym = 0

 		do i = 1, gt1
			call derivs_HEOM(y,dy,-1)
			!! rk4 a la Numerical Recipes (without 'time', since derivs in HEOM are time-independent)
			yt = y + hh*dy
			call derivs_HEOM(yt,dyt,-1)
			yt = y + hh*dyt
			call derivs_HEOM(yt,dym,-1)
			yt = y + dt*dym ! buvo hh, t.y. klaida...
			dym = dym + dyt
			call derivs_HEOM(yt,dyt,-1)
			y = y + h6 * (dy + dyt + 2.0_dp*dym)
		end do

	end subroutine make_step_HEOM



	subroutine derivs_HEOM(y,dy,lower)
		integer, intent(in) :: lower
		complex(dpc), intent(in),  dimension(lower:,:,:) :: y
		complex(dpc), intent(out), dimension(lower:,:,:) :: dy
		integer :: I, n, k, a, b
		real(dp), dimension(size(y,2),size(y,2)) :: Q_i ! projector

		 dy = 0.0_dp

!		Q_i = 0.0_dp
!
!		do I = 0, NN-1 ! check the upper bound...
!			dy(I,:,:) = - (0.0_dp,1.0_dp) * (matmul(HH,y(I,:,:)) - matmul(y(I,:,:),HH)) &
!						- sum( gamma_ik(:,:) * index_n(I,:,:)) * y(I,:,:)
!
!			do n = 1, Ns
!				Q_i(n,n) = 1.0_dp
!				! 'residual' terms
!				dy(I,:,:) = dy(I,:,:) - 0.5_dp*delta_i(n) * & !
!								( matmul(Q_i,y(I,:,:)) - 2.0_dp*matmul(Q_i,matmul(y(I,:,:),Q_i)) + &
!								matmul(y(I,:,:),Q_i) ) ! because QQ == Q
!                do k = 1, Ks+1
!					! '-1' neighbours
!					dy(I,:,:) = dy(I,:,:) &
!								+ (0.0_dp,1.0_dp) * real(c_ik(n,k)) * index_n(I,n,k) * &
!								( matmul(Q_i,y(index_m1(I,n,k),:,:)) - matmul(y(index_m1(I,n,k),:,:),Q_i) ) &
!								+ aimag(c_ik(n,k)) * index_n(I,n,k) * & ! "-" 'correctly'...
!								( matmul(Q_i,y(index_m1(I,n,k),:,:)) + matmul(y(index_m1(I,n,k),:,:),Q_i) )
!					! '+1' neightbours
!					dy(I,:,:) = dy(I,:,:) + (0.0_dp,1.0_dp) * &
!								( matmul(Q_i,y(index_p1(I,n,k),:,:)) - matmul(y(index_p1(I,n,k),:,:),Q_i) )

					!! normalized version
					! '-1' neighbours
!                    dy(I,:,:) = dy(I,:,:) + (0.0_dp,1.0_dp) * &
!                                sqrt( index_n(I,n,k) / sqrt( c_ik(n,k) * conjg(c_ik(n,k)) ) ) * &
!                                ( c_ik(n,k) *matmul(Q_i,y(index_m1(I,n,k),:,:)) - &
!                            conjg(c_ik(n,k))*matmul(y(index_m1(I,n,k),:,:),Q_i) )
!                    ! '+1' neightbours
!                    dy(I,:,:) = dy(I,:,:) + (0.0_dp,1.0_dp) * &
!                                sqrt( (index_n(I,n,k) + 1) * sqrt( c_ik(n,k) * conjg(c_ik(n,k)) ) ) * &
!                                ( matmul(Q_i,y(index_p1(I,n,k),:,:)) - matmul(y(index_p1(I,n,k),:,:),Q_i) )

!					Q_i(n,n) = 0.0_dp
!				end do
!			end do
!		end do


		!!!!!!!!!!
		!!
		!! alternative version: matrix elements
		!!
		!!!!!!!!!!
		call init_kroneker(Ns)

        do I = 0, NN-1 ! check the upper bound...
            do a = 1, Ns
            do b = 1, Ns

                do n = 1, Ns
                    dy(I,a,b) = dy(I,a,b) - (0.0_dp,1.0_dp) * ( HH(a,n)*y(I,n,b) - y(I,a,n)*HH(n,b) )
                end do

                dy(I,a,b) = dy(I,a,b) - sum( gamma_ik(:,:) * index_n(I,:,:) ) * y(I,a,b) &
                                      - 0.5_dp * (delta_i(a) + delta_i(b) - 2.0_dp*delta_i(a)*kroneker(a,b)) * y(I,a,b)

                do k = 1, Ks+1
                    dy(I,a,b) = dy(I,a,b) &
                    ! '-1' neighbours
                              + (0.0_dp,1.0_dp) *       conjg(norma_m1(I,a,k))  * y(index_m1(I,a,k),a,b) &
                              - (0.0_dp,1.0_dp) * (norma_m1(I,b,k)) * y(index_m1(I,b,k),a,b) &
                    ! '+1' neightbours
                              + (0.0_dp,1.0_dp) * norma_p1(I,a,k) * y(index_p1(I,a,k),a,b) &
                              - (0.0_dp,1.0_dp) * norma_p1(I,b,k) * y(index_p1(I,b,k),a,b)
                end do

            end do
            end do
        end do

		call delete_kroneker()


	end subroutine derivs_HEOM



end module HEOM
