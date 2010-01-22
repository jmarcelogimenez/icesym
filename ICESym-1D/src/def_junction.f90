module def_junction
  
  use def_simulator
  use utilities
  use gasdyn_utils
  use, intrinsic :: ISO_C_BINDING
  
  type, BIND(C) :: this
    integer(C_INT) :: nnod, ndof, modelo_junc, nnod_input
  end type this

contains
  
  subroutine state_initial_junction(myData, atm, globalData, state_ini, type_end) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(this):: myData
    type(dataSim) :: globalData
    integer(C_INT), dimension(myData%nnod) :: type_end
    real(C_DOUBLE) :: state_ini(0:((myData%nnod+myData%nnod_input)*myData%ndof)-1)
    real(C_DOUBLE) :: atm(0:2)
    
    do i=0,myData%nnod-1
       state_ini(i*myData%ndof)   = atm(0)
       state_ini(i*myData%ndof+1) = atm(1)
       state_ini(i*myData%ndof+2) = atm(2)
    enddo

  end subroutine state_initial_junction
  
  subroutine solve_junction(myData, globalData, state, new_state, type_end, dx_tube, Area_tube, &
       twall_tube, dAreax_tube) BIND(C)
    use, intrinsic :: ISO_C_BINDING
    type(this) :: myData
    type(dataSim) :: globalData 
    integer(C_INT), dimension(myData%nnod) :: type_end
    real(C_DOUBLE), dimension(myData%nnod) :: dx_tube,Area_tube, twall_tube, dAreax_tube
    real(C_DOUBLE) :: state(0:((myData%nnod+myData%nnod_input)*myData%ndof)-1)
    real(C_DOUBLE) :: new_state(0:(myData%nnod*myData%ndof)-1)
    
    integer :: i,j,viscous_flow,heat_flow
    real*8 :: ga,R_gas,dt
    real*8 :: Area_P,dAreax_P,Twall_P
    real*8, dimension(myData%ndof) :: RHS
    real*8, dimension(myData%ndof,myData%nnod) :: RHSv,U,Uref
    
    ga    = globalData%ga
    R_gas = globalData%R_gas
    dt    = globalData%dt

    viscous_flow = globalData%viscous_flow
    heat_flow    = globalData%heat_flow

    do i=0,myData%nnod-1
       do j=0,myData%ndof-1
          U(j+1,i+1)    = state(i*myData%ndof+j)
          Uref(j+1,i+1) = state((i+myData%nnod)*myData%ndof+j)
       end do
    end do

    RHSv = 0.0d0
    do i=1,myData%nnod
       Area_P   = Area_tube(i)
       dAreax_P = dAreax_tube(i)
       Twall_P  = twall_tube(i)
       
       call rhschar(U, Area_P, dAreax_P, Twall_P, ga, R_gas, dt, &
            viscous_flow, heat_flow, RHS)
       RHSv(:,i) = RHS
    end do

    call solve_equation_system(myData, Uref, Area_tube, type_end, RHSv, &
         ga, U)

    do i=0,myData%nnod-1
       do j=0,myData%ndof-1
          new_state(i*myData%ndof+j) = U(j+1,i+1)
       enddo
    enddo

  end subroutine solve_junction

  subroutine solve_equation_system(myData, Uref, Area_tube, type_end, RHS, &
       ga, U)
    !
    !
    !
    implicit none

    type(this), intent(in) :: myData
    integer, dimension(myData%nnod), intent(in) :: type_end
    real*8, intent(in) :: ga
    real*8, dimension(myData%nnod), intent(in) :: Area_tube
    real*8, dimension(myData%ndof,myData%nnod), intent(in) :: RHS,Uref
    real*8, dimension(myData%ndof,myData%nnod), intent(inout) :: U

    integer :: nnod,ndof,niter,niter_max,i
    real*8 :: tol,normres,tol_junction,relax
    real*8, dimension(myData%ndof*myData%nnod) :: R,dU
    real*8, dimension(myData%ndof,myData%nnod) :: Uini
    real*8, dimension(myData%ndof*myData%nnod,myData%ndof*myData%nnod) :: J
    logical :: continue_flag

    nnod = myData%nnod
    ndof = myData%ndof

    niter_max = 50
    tol       = 1.0d-6

    Uini = U

    niter = 0
    continue_flag = .true.
    do while(continue_flag)
       R  = 0.0d0
       J  = 0.0d0
       call junction_residue(nnod, U, Uref, Area_tube, type_end, RHS, &
            ga, R)
       call junction_jaco(nnod, U, Uref, Area_tube, type_end, RHS, &
            ga, J)
       call linear_solver(R, J, dU, nnod*ndof)
       
       if(any(isnan(dU))) then
          write(*,*) 'NANs in JUNCTION SOLVER'
          stop
       end if
       relax = tanh(0.05*(niter+1))
       do i=1,nnod
          U(:,i) = U(:,i)-relax*dU(ndof*(i-1)+1:ndof*i)
       end do
       
       normres = normvec(R, nnod*ndof, 2)
       niter   = niter+1
       if((normres.lt.tol).or.(niter.gt.niter_max)) continue_flag = .false.
    end do
    
  end subroutine solve_equation_system

  subroutine junction_residue(ntubes, U, Uref, Areas, type_end, RHS, ga, res)
    !
    !
    !
    implicit none
    
    integer, intent(in) :: ntubes
    integer, dimension(ntubes), intent(in) :: type_end
    real*8, intent(in) :: ga
    real*8, dimension(ntubes), intent(in) :: Areas
    real*8, dimension(3,ntubes), intent(in) :: Uref,RHS
    real*8, dimension(3,ntubes), intent(inout) :: U
    real*8, dimension(3*ntubes), intent(out) :: res
    
    integer :: iq,ip,j,np,nq,pipe_count
    integer, dimension(ntubes) :: type_end_junction
    real*8 :: mdot,edot,g1,p_1,h_1,delta
    real*8 :: rho_P,v_P,p_P,a_P,rho_S,v_S,p_S,a_S
    logical :: inflow
    
    g1    = ga-1.
    delta = 0.5*g1
    
    ip = 0
    iq = 0
    
    res  = 0.0d0
    
    type_end_junction = 0
    call junction_type(ntubes, U, type_end, type_end_junction)
    
    np = sum(iabs(type_end_junction), 1, type_end_junction.gt.0)
    nq = sum(iabs(type_end_junction), 1, type_end_junction.lt.0)
    
    pipe_count = -1
    
    do j=1,ntubes
       
       pipe_count = pipe_count+1
       
       rho_P = U(1,j)
       v_P   = U(2,j)
       p_P   = U(3,j)
       
       rho_S = Uref(1,j)
       v_S   = Uref(2,j)
       p_S   = Uref(3,j)
       
       a_P  = dsqrt(ga*p_P/rho_P)
       a_S  = dsqrt(ga*p_S/rho_S)
       
       mdot = rho_P*type_end(j)*v_P*Areas(j)
       edot = mdot/g1*(a_P**2 + delta*v_P**2)
       
       inflow = .false.
       if(type_end_junction(j).gt.0) inflow = .true.
       
       ! Mass balance
       res(1) = res(1) + mdot
       ! Energy balance
       res(2) = res(2) + edot
       ! Mach line characeristics
       if(type_end(j).gt.0) then
          res(j+2) = (p_P-p_S)/(rho_S*a_S)-v_S+RHS(2,j)+v_P
       else
          res(j+2) = (p_P-p_S)/(rho_S*a_S)+v_S+RHS(3,j)-v_P
       end if
       
       ! Equalize all the pressures
       if(pipe_count.gt.0) then
          res(2+ntubes+pipe_count) = p_P - p_1
       else
          p_1 = p_P
       end if
       
       ! Path line characteristics
       if(inflow) then
          ip = ip+1
          if(ip.gt.ntubes) stop ' ====== ERROR ====== @ junction_residue '
          res(2*ntubes+ip+1) = (p_P/p_S)**(1./ga)-rho_P/rho_S*dexp(RHS(1,j))
       endif
    end do
    
    ! Equalize all the outgoing enthalpy
    if(nq.gt.1) then
       do j=1,ntubes
          rho_P = U(1,j)
          v_P   = U(2,j)
          p_P   = U(3,j)
          a_P   = dsqrt(ga*p_P/rho_P)
          
          inflow = .false.
          if(type_end_junction(j).gt.0) inflow = .true.
          
          if(.not.inflow) then
             if(iq.gt.0) then
                res(2*ntubes+ip+iq+1) = (a_P**2+delta*v_P**2) - h_1
                iq = iq+1
             else
                iq = 1
                h_1 = a_P**2+delta*v_P**2
             endif
          endif
       end do
    endif
    
  end subroutine junction_residue
  
  subroutine junction_type(ntubes, U, no, type_end)
    !
    !
    !
    implicit none
    
    integer, intent(in) :: ntubes
    integer, dimension(ntubes), intent(in) :: no
    real*8, dimension(3,ntubes), intent(inout) :: U
    integer, dimension(ntubes), intent(out) :: type_end

    integer :: np,nq,k
    real*8 :: tol,pt_avg
    real*8, dimension(ntubes) :: pt,vn

    tol = 1.0d-06
    ! Normal velocity
    vn = U(2,:)*no

    where(vn.gt. tol) type_end =  1
    where(vn.lt.-tol) type_end = -1

    np = sum(iabs(type_end), 1, type_end.gt.0)
    nq = sum(iabs(type_end), 1, type_end.lt.0)

    if((np+nq).lt.ntubes) then
       ! There are undefined ends
       if(np+nq .eq. 0) then
          ! All are undefinied
          pt     = U(3,:)+0.5*U(1,:)*U(2,:)**2
          pt_avg = sum(pt)/ntubes

          where(pt.gt.pt_avg) type_end =  1
          where(pt.le.pt_avg) type_end = -1

          np = sum(iabs(type_end), 1, (type_end.gt.0))
          nq = sum(iabs(type_end), 1, (type_end.lt.0))

          if((np.eq.0).or.(nq.eq.0)) then
             ! Undefined by normal velocity or by total pressure
             type_end = no
          endif
       else
          ! There are some undefinied
          if(np.eq.0) then
             do k=1,ntubes
                if(type_end(k).eq.0) type_end(k) = 1
             end do
          else
             do k=1,ntubes
                if(type_end(k).eq.0) type_end(k) = -1
             end do
          endif
       endif
    else
       if((np.eq.0).or.(nq.eq.0)) then
          ! All are undefinied but have the same sense
          pt     = U(3,:)+0.5*U(1,:)*U(2,:)**2
          pt_avg = sum(pt)/ntubes
          where(pt.gt.pt_avg)
             type_end = 1
          elsewhere
             type_end = -1
          endwhere

          np = sum(iabs(type_end), 1, type_end.gt.0)
          nq = sum(iabs(type_end), 1, type_end.lt.0)
       endif
       if((np.eq.0).or.(nq.eq.0)) then
          ! Undefined by normal velocity or by total pressure
          type_end = no
       endif
    endif

    ! Check sign of velocities
    where(type_end*no.gt.0) U(2,:) = +dabs(U(2,:))
    where(type_end*no.lt.0) U(2,:) = -dabs(U(2,:))

  end subroutine junction_type

  subroutine junction_jaco(ntubes, U, Uref, Areas, type_end, RHS, ga, mat)
    !
    !
    !
    implicit none

    integer, intent(in) :: ntubes
    integer, dimension(ntubes), intent(in) :: type_end
    real*8, intent(in) :: ga
    real*8, dimension(ntubes), intent(in) :: Areas
    real*8, dimension(3,ntubes), intent(in) :: Uref,RHS
    real*8, dimension(3,ntubes), intent(inout) :: U
    real*8, dimension(3*ntubes,3*ntubes), intent(out) :: mat

    integer :: j,k
    real*8 :: tol
    real*8, dimension(3) :: dU
    real*8, dimension(3,ntubes) :: Up
    real*8, dimension(3*ntubes) :: res_p,res_m

    tol = 1.0d-10

    mat = 0.0d0

    do j=1,ntubes
       dU = 1.0e-7*U(:,j)
       do k=1,3
          if(dabs(dU(k)).lt.tol) dU(k) = tol
       end do
       do k=1,3
          Up      = U
          Up(k,j) = Up(k,j)+dU(k)
          call junction_residue(ntubes, Up, Uref, Areas, type_end, RHS, ga, res_p)
          Up      = U
          Up(k,j) = Up(k,j)-dU(k)
          call junction_residue(ntubes, Up, Uref, Areas, type_end, RHS, ga, res_m)
          mat(:,3*(j-1)+k) = (res_p-res_m)/(2.0*dU(k))
       end do
    end do

  end subroutine junction_jaco

end module def_junction
