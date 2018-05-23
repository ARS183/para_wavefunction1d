!subroutine wavefunc1d_solver(x,u,uexat,nx,hx,dt,Nt)
subroutine wavefunc1d_solver(x,u,Nt)
    include 'openNS3d.h'		
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::uexat1(1-LAP:nx+LAP),uexat2(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::uexat3(1-LAP:nx+LAP),uexat4(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::R(1-LAP:nx+LAP),x(1-LAP:nx+LAP)

    real(kind=OCFD_REAL_KIND),allocatable :: u0_global(:),u_global(:),x_global(:)
    real(kind=OCFD_REAL_KIND),allocatable :: uexat1_global(:),uexat2_global(:)
    real(kind=OCFD_REAL_KIND),allocatable :: uexat3_global(:),uexat4_global(:)

    integer :: Nt,NT1,NT2,NT3,NT4

    NT1=2.4d0/dt
    NT2=4.2d0/dt
    NT3=9.d0/dt
    NT4=15.d0/dt
    
    call check_x1d(x)

    call analyticSol_TKS_u0(x,u)
    call analyticSol_TKS_uexat(x,uexat1,dt,NT1)
    call analyticSol_TKS_uexat(x,uexat2,dt,NT2)
    call analyticSol_TKS_uexat(x,uexat3,dt,NT3)
    call analyticSol_TKS_uexat(x,uexat4,dt,NT4)

    allocate(uexat1_global(1:nx_global),uexat2_global(1:nx_global))
	allocate(uexat3_global(1:nx_global),uexat4_global(1:nx_global))
	allocate(u0_global(1:nx_global),x_global(1:nx_global))

    call write_data_global(x,x_global)
	call write_data_global(u,u0_global)
	call write_data_global(uexat1,uexat1_global)
	call write_data_global(uexat2,uexat2_global)
	call write_data_global(uexat3,uexat3_global)
	call write_data_global(uexat4,uexat4_global)

	if  (my_id .ne. 0) then
		deallocate(x_global,u0_global,uexat1_global,uexat2_global,uexat3_global,uexat4_global)
	endif

!-------------------------------------------------------------
	if (my_id .eq. 0) then
		call Toplt1D(x_global,u0_global,'u0',2,nx_global)
		call Toplt1D(x_global,uexat1_global,'uexat1',6,nx_global)
		call Toplt1D(x_global,uexat2_global,'uexat2',6,nx_global)
		call Toplt1D(x_global,uexat3_global,'uexat3',6,nx_global)
		call Toplt1D(x_global,uexat4_global,'uexat4',6,nx_global)
	endif



    do 800 i=1,Nt
        call RK4(u)
        if (i==NT1) then
            allocate(u_global(1:nx_global))
            call write_data_global(u,u_global)
            if (my_id .eq. 0) then
                call Toplt1D(x_global,u_global,'u1',2,nx_global)
            endif
            deallocate(u_global)
        elseif (i==NT2) then
            allocate(u_global(1:nx_global))
            call write_data_global(u,u_global)
            if (my_id .eq. 0) then
                call Toplt1D(x_global,u_global,'u2',2,nx_global) 
            endif
            deallocate(u_global)
        elseif (i==NT3) then
            allocate(u_global(1:nx_global))
            call write_data_global(u,u_global)
            if (my_id .eq. 0) then
            call Toplt1D(x_global,u_global,'u3',2,nx_global) 
            endif
            deallocate(u_global)
        elseif (i==NT4) then
            allocate(u_global(1:nx_global))
            call write_data_global(u,u_global)
            if (my_id .eq. 0) then
            call Toplt1D(x_global,u_global,'u4',2,nx_global) 
            endif
            deallocate(u_global)
        endif

800	continue


end subroutine wavefunc1d_solver


!subroutine RK4(nx,u,dt,hx)	!4阶RK方法
subroutine RK4(u)
    implicit none
    include 'openNS3d.h'
    real(kind=OCFD_REAL_KIND) :: u(1-LAP:nx+LAP), R(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND) :: R1(1-LAP:nx+LAP),R2(1-LAP:nx+LAP),R3(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND):: u1(1-LAP:nx+LAP),u2(1-LAP:nx+LAP),u3(1-LAP:nx+LAP)
!    real(kind=OCFD_REAL_KIND):: hx,dt
!    integer :: nx

    call check_x1d(u)
	call ComputeR(u,R)	

	u1(1:nx)=u(1:nx)+(dt/2.d0)*R(1:nx)

    call check_x1d(u1)
	call ComputeR(u1,R1)

	u2(1:nx)=u(1:nx)+(dt/2.d0)*R1(1:nx)

    call check_x1d(u2)
    call ComputeR(u2,R2)

	u3(1:nx)=u(1:nx)+dt*R2(1:nx)
    
    call check_x1d(u3)
    call ComputeR(u3,R3)

    u(1:nx)=u(1:nx)+(dt/6.d0)*(R(1:nx)+2.d0*R1(1:nx)+2.d0*R2(1:nx)+R3(1:nx))

!    deallocate(R1,R2,R3,u1,u2,u3)
    
    end subroutine RK4
        



!subroutine ComputeR(nx,u,R,hx)  !计算右端项
subroutine ComputeR(u,R)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP),R(1-LAP:nx+LAP)
!	integer :: nx

    call Du2Dx_PADE4(u,d2f)
    call DuDx_UCC_UpWind(u,d2f,df)

    R(1:nx)=-0.5d0*df(1:nx)
!    deallocate(d2f,df)

end subroutine ComputeR





subroutine analyticSol(x,u,uexat,deltat,Nt)
    include 'openNS3d.h'
    integer :: i,Nt
    real(kind=OCFD_REAL_KIND):: x(1-LAP:nx+LAP),u(1-LAP:nx+LAP),uexat(1-LAP:nx+LAP)
   real(kind=OCFD_REAL_KIND):: deltat
        	
    do i=1,nx
        u(i)=exp(-(x(i)-5.d0)**2)
        uexat(i)=exp(-(x(i)-5.d0-dble(Nt)*deltat)**2)
    enddo
    
end subroutine



subroutine analyticSol_TKS_u0(x,u)
    include 'openNS3d.h'
    integer :: i,Nt
    real(kind=OCFD_REAL_KIND):: x(1-LAP:nx+LAP),u(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND):: deltat,phi
    phi=50.d0
 !    k0=0.838242d0/hx
    k0=0.7857d0/hx
        	
    do i=1,nx
        u(i)=exp(-phi*(x(i)-1.5d0)**2)*dsin(k0*x(i))
    enddo
    
end subroutine


subroutine analyticSol_TKS_uexat(x,uexat,deltat,nt)
    include 'openNS3d.h'
    integer :: i,nt
    real(kind=OCFD_REAL_KIND):: x(1-LAP:nx+LAP),uexat(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND):: deltat,phi,c
    phi=50.d0
!    k0=0.838242d0/hx
    k0=0.7857d0/hx
    c=0.5d0
        	
    do i=1,nx
        uexat(i)=exp(-phi*(x(i)-1.5d0-c*dble(nt)*deltat)**2)*dsin(k0*(x(i)-c*dble(nt)*deltat))
    enddo
    
end subroutine