!
! prepares the indexing of ADO's (index) and their neighbors (index_[p/m]1)
!
module prep_rules

	use std_types
	use resources

	implicit none

	integer, dimension(:), allocatable 	 :: u	! number of elements in the l'th tier
	integer, dimension(:,:), allocatable :: nm	! the index vector/matrix of an ADO



contains

	subroutine create_index()
		integer :: n ! numbering sites
		integer :: k ! numbering poles
		integer :: l ! numbering tiers
		integer :: I ! numbering ADO's
		integer :: j, jj ! counting elements generated in a tier
		logical	:: vector_is_unique
		integer :: a, b, c, check

		allocate(u(0:Ls))
		allocate(nm(Ns,Ks+1))

	print *, "Creating index matrix for HEOM..."

		! the l'th tier has u(l-1) ADO's: 0'th has 1, etc.
		u(Ls) = NN
		do l = 0, Ls-1
			u(l) = binomial(Ns*(Ks+1)+l,l)
		end do

		! create the 0'th entry in index_n
		index_n = 0
		index_p1 = -1
		index_m1 = -1
		nm = 0

		! create the 1'st tier in index_n and '+/- 1' neighbour matrix
		I = 1
		do n = 1, Ns
		  do k = 1, Ks+1
			  index_n(I,n,k) = 1
			  index_p1(0,n,k) = I
			  index_m1(I,n,k) = 0
			  I = I + 1
		  end do
		end do
!print *, "Hola DEA!"
		! create subsequent tiers "recurrently"
		do l = 1, Ls-1

			j = u(l)	! counts successful tries in a tier

			index_n(j,1,1) = index_n(u(l-1),1,1) + 1 ! first element needed explicitly for later procedures
			index_p1(u(l-1),1,1) = j
			index_m1(j,1,1) = u(l-1)

			do I = u(l-1), u(l)-1 ! because u(l) is already the first element from the next tier

				nm = index_n(I,:,:) ! pick the I'th index matrix from a lower tier and try adding +1 in diff. positions

				do n = 1, Ns
					do k = 1, Ks+1
				  		nm(n,k) = nm(n,k) + 1	! generate new index matrix and check if it is unique:
				  		!***********************************************************
						check = 0 ! = matrices are differrent

						do0: do a = u(l), j

							do1: do b = 1, Ns
								 do c = 1, Ks+1
									if (nm(b,c) /= index_n(a,b,c)) then
										check = 0 ! matrices are different at least at the element j,k
										exit do1
									else
										check = 1
									end if
								end do
							end do do1

							if (check > 0) then
								vector_is_unique = .false. ! = there's at least one matching matrix
								index_p1(I,n,k) = a
								index_m1(a,n,k) = I
								exit do0
							else
								vector_is_unique = .true.
							end if

						end do do0
				  		!***********************************************************
				  		if (vector_is_unique) then	! check if the new vector is unique
							j = j + 1
							index_n(j,:,:) = nm		! record the newly generated index
							index_p1(I,n,k) = j		! register +1 neighbourhood
							index_m1(j,n,k) = I		! register -1 neighbourhood (reciprocity)
				  		end if
				  		nm(n,k) = nm(n,k) - 1	! prepare index matrix for a new try
					end do
				end do

			end do

		end do

		! writing out the index matrix into an external file; needs to have a trigger
		call write_index_to_file()

		deallocate(u, nm)

	print *, "Creating index matrix: done"

	end subroutine create_index



	subroutine create_coefficients()
	    integer :: I, n, k
        character(len = 20) :: buffer

!        open(unit=12,file=(trim(out_dir)//"/"//trim("norma_m1.txt")), position='append')
!        buffer = "(3i6,  2f16.8)"
	    do I = 0, NN-1
	        do n = 1, Ns
	            do k = 1, Ks+1
	                norma_m1(I,n,k) = sqrt(  index_n(I,n,k)      / sqrt( c_ik(n,k) * conjg(c_ik(n,k)) ) ) *  c_ik(n,k)
	                norma_p1(I,n,k) = sqrt( (index_n(I,n,k) + 1) * sqrt( c_ik(n,k) * conjg(c_ik(n,k)) ) )
!	                write(12,buffer) I, n, k, real(norma_m1(I,n,k)), aimag(norma_m1(I,n,k))
	            end do
	        end do
	    end do
!	    close(unit=12)

    end subroutine create_coefficients



	subroutine write_index_to_file()
		integer :: I, n, k, eol
		integer, dimension(size(nm,1)*size(nm,2)) :: index_row
		character(len = 20) 	:: buffer

		write(buffer,'(i4)') (Ns*(Ks+1))+1
		buffer = '(' // trim(buffer) // 'i5)' ! forming a descriptor for outputting row by row

		open(unit=12,file=(trim(out_dir)//"/"//trim("index_matrix.dat")), status='replace')!, position='append')
		do I = 0, NN-1
			do n = 1, Ns
				do k = 1, Ks+1
					index_row((n-1)*(Ks+1)+k) = index_n(I,n,k)
				end do
			end do
			write(12,buffer) I, index_row(1:Ns*(Ks+1))
		end do
		close(unit=12)

		open(unit=12,file=(trim(out_dir)//"/"//trim("index_p1_matrix.dat")), status='replace')!, position='append')
		do I = 0, NN-1
			do n = 1, Ns
				do k = 1, Ks+1
					index_row((n-1)*(Ks+1)+k) = index_p1(I,n,k)
				end do
			end do
			write(12,buffer) I, index_row(1:Ns*(Ks+1))
		end do
		close(unit=12)

		open(unit=12,file=(trim(out_dir)//"/"//trim("index_m1_matrix.dat")), status='replace')!, position='append')
		do I = 0, NN-1
			do n = 1, Ns
				do k = 1, Ks+1
					index_row((n-1)*(Ks+1)+k) = index_m1(I,n,k)
				end do
			end do
			write(12,buffer) I, index_row(1:Ns*(Ks+1))
		end do
		close(unit=12)

	end subroutine write_index_to_file



end module prep_rules
