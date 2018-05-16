!-------------------------------------------------------------
!Êý¾Ý»ã×Ü
!-------------------------------------------------------------

subroutine write_data_global(phi,phi_global)
	include 'openNS3d.h'
	integer :: i_global,ii,id_mid,ka,i,count,inum,iter
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP) :: phi
	real(kind=OCFD_REAL_KIND),allocatable,dimension(:) :: temp_phi
	real(kind=OCFD_REAL_KIND),dimension(1:nx_global) :: phi_global

if (my_id .eq. 0) then
    
    do 100 i=0,npx0-1		
		ka=i
			
		if (ka .eq. 0) then
			do ii=1,nx
					i_global=i_offset(i)-1+ii
					phi_global(i_global)=phi(ii)
			enddo

		else
			allocate(temp_phi(1:i_nn(i)))

			call MPI_RECV(temp_phi,i_nn(i),OCFD_DATA_TYPE,  &
			ka,1,MPI_COMM_WORLD,status,ierr) 
			
			do ii=1,i_nn(i)
				i_global=i_offset(i)-1+ii
				phi_global(i_global)=temp_phi(ii)
			enddo	
	
			deallocate(temp_phi)

		endif
100	continue
	
else
	allocate(temp_phi(1:i_nn(npx)))
	do ii=1,i_nn(npx)

		temp_phi(ii)=phi(ii)

	enddo
				
	call MPI_SEND(temp_phi,i_nn(npx),OCFD_DATA_TYPE,  &
	0,1,MPI_COMM_WORLD,ierr) 

	deallocate(temp_phi)

endif
end subroutine write_data_global