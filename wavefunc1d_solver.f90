!subroutine wavefunc1d_solver(x,u,uexat,nx,hx,dt,Nt)
subroutine wavefunc1d_solver(x,u,uexat,Nt)
    include 'openNS3d.h'		
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),uexat(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::R(1-LAP:nx+LAP),x(1-LAP:nx+LAP)
    integer :: Nt
    
    call check_x1d(x)

    call analyticSol(x,u,uexat,dt,Nt)

    do 800 i=1,Nt

        

        call RK4(u)

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

    if (npx==0) then
        call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,1)
        call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,1)
        call D2F_PADE4(u,d2f,nx,hx)
        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)
        call DF_UCC45_P(u,d2f,df,nx,hx,1)
    elseif (npx==npx0-1) then
        call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,2)
        call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,2)
        call D2F_PADE4(u,d2f,nx,hx)
        call OCFD_DF_SB_UCC45(u,df,nx,hx)
        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)
        call DF_UCC45_P(u,d2f,df,nx,hx,2)
    else
        call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,0)
        call D2F_PADE4(u,d2f,nx,hx)
        call OCFD_DF_SB_UCC45(u,df,nx,hx)
        call DF_UCC45_P(u,d2f,df,nx,hx,0)
    endif

    R(1:nx)=-df(1:nx)
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