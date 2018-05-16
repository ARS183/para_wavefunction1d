!=========================================================
!==============�������ṩ�������ܣ����������Tecplot�ļ�==	
!=========================================================	
	
	subroutine Toplt1D(x,v,name,n_name,N)
!	n_name:�ļ��洢�ı��������ȣ������ܳ��ȣ����ֳ���Ĭ��Ϊ6λ
!	eg: call Toplt(x,y,p,'p',1,0.001,400,500)
	include 'openNS3d.h'
	integer :: n_name,N
	real(kind=OCFD_REAL_KIND),dimension(1:nx_global) :: x,v
	character(len=n_name) name
	character*20 fname

	fname =name//'.plt'
	
	open(1,file=fname,status='unknown')
	write(1,*) 'variables=x,'//name
	write(1,*) 'zone i=',N
	do 10 i=1,N
	  write(1,'(8E25.15)') x(i),v(i)
10	continue	
	close(1)

	end subroutine Toplt1D

!
!------------------------------------------------------------------
!

	subroutine Toplt2D(x,y,v,name,n_name,Nx,Ny)
	double precision x,y,v
	character(len=n_name) name
	character*20 fname
	dimension x(Nx,Ny),y(Nx,Ny),v(Nx,Ny)

	fname =name//'.plt'
	
	open(1,file=fname,status='unknown')
	write(1,*) 'variables=x,y,'//name
	write(1,*) 'zone i=',Nx,'j=',Ny
	do 10 j=1,Ny
	do 10 i=1,Nx
	  write(1,'(8E25.15)') x(i,j),y(i,j),v(i,j)
10	continue	
	close(1)

	end subroutine Toplt2D