Program SV
  ! SV=Steady Viscous
  use, intrinsic :: iso_c_binding
  use mod_BuildMatrices
  implicit none
  character(len=20) :: string
  ! Grid variables
  real (c_double),allocatable :: xp(:,:),xu(:,:),Soln(:)
  integer (c_int),allocatable :: tp(:,:),tu(:,:),tuplot(:,:)
  integer (c_int),allocatable :: N(:),W(:),S(:),E(:)
  integer (c_int) :: nx,nu,ne,nnz,k,k2
  ! for SLU
  integer(c_int) :: nprocs,nrhs=1,info,nxSLU,nnzSLU
  ! date_time return
  integer (c_int) :: values(8) 
  real (c_double)  :: start,finish  ! timers

  call Read_XT()
  call get_environment_variable("OMP_NUM_THREADS",string)
  if (string=="") then
     print *,"ERROR, OMP_NUM_THREADS has not been set"
     stop
  end if
  string=trim(string)
  read(string,*)nprocs
  allocate(Soln(2*nu+nx))  
  call date_and_time(VALUES=values)
  start=60*values(6)+values(7)+dble(values(8))/1000.
  call BuildMatrices(xu,tu,nx,nu,ne,nnz,N,E,S,W,Soln)
  call date_and_time(VALUES=values)
  finish=60*values(6)+values(7)+float(values(8))/1000.
  print *,'SV build matrices, wall clock:',finish-start
  ! redefine nx, nnz based on Int,Int matrix
  nx=size(LHS%colptr)-1;nnz=size(LHS%rowind)

  ! Solve 
  call c_bridge_pdgssv(nprocs,nx,nnz,nrhs,LHS%data,&
       &LHS%rowind,LHS%colptr,tildeb,nx,info)
  print *,'info=',info
  Soln(Interior)=tildeb
  ! Save Soln
  call Write_Vector_RDP("Soln",Soln)
contains

  Subroutine Read_XT()
    use, intrinsic :: iso_c_binding
    ! order of variables in XT.mat
    !save XT.mat xp tp xu tu tuplot N W S E
    implicit none
    integer (c_int) :: unitno = 10
    integer (c_int) :: nn,k
    open(unit=unitno,file="XT.mat",form="formatted"&
         &,status="old",action="read")
    call Skip(unitno,3)
    ! read xp and tp
    read(unitno,fmt="(8x,i14)")nx
    !print *,'xp nx= ',nx
    allocate (xp(nx,2))
    call Skip(unitno,1)
    do k=1,nx
       read(unitno,*)xp(k,1),xp(k,2)
    end do
    ! on to tp
    call Skip(unitno,4)
    read(unitno,fmt="(8x,I14)")ne
    !print *,'tp ne= ',ne
    allocate (tp(ne,3))
    call Skip(unitno,1)
    do k=1,ne
       read(unitno,*)tp(k,1),tp(k,2),tp(k,3)
    end do
    call Skip(unitno,4)
    ! read xu and tu, tuplot
    read(unitno,fmt="(8x,i14)")nu
    !print *,'xu  nu= ',nu
    allocate (xu(nu,2))
    call Skip(unitno,1)
    do k=1,nu
       read(unitno,*)xu(k,1),xu(k,2)
    end do
    ! on to tu
    call Skip(unitno,4)
    read(unitno,fmt="(8x,I14)")ne
    allocate (tu(ne,6))
    !print *,'tu ne= ',ne
    call Skip(unitno,1)
    do k=1,ne
       read(unitno,*)tu(k,1:6)
    end do
    ! on to tuplot
    call Skip(unitno,4)
    read(unitno,fmt="(8x,I14)")nn
    allocate (tuplot(nn,6))
     !print *,'tuplot nn= ',nn
     call Skip(unitno,1)
    do k=1,nn
       read(unitno,*)tuplot(k,1:3)
    end do
    ! on to N
    call Skip(unitno,4)
    read(unitno,fmt="(8x,I14)")nn
    allocate (N(nn))
    !print *,'N nn= ',nn
    call Skip(unitno,1)
    do k=1,nn
       read(unitno,*)N(k)
    end do
    ! on to W
    call Skip(unitno,4)
    read(unitno,fmt="(8x,I14)")nn
    !nn=1
    allocate (W(nn))
    !print *,'W nn= ',nn
    call Skip(unitno,1)
    do k=1,nn
       read(unitno,*)W(k)
    end do
    ! on to S
    call Skip(unitno,4)
    read(unitno,fmt="(8x,I14)")nn
    allocate (S(nn))
    !print *,'S nn= ',nn
    call Skip(unitno,1)
    do k=1,nn
       read(unitno,*)S(k)
    end do
    ! on to E
    call Skip(unitno,4)
    read(unitno,fmt="(8x,I14)")nn
    !nn=1
    allocate (E(nn))
    !print *,'E nn= ',nn
    call Skip(unitno,1)
    do k=1,nn
       read(unitno,*)E(k)
    end do
    close(unit=unitno)
  End Subroutine Read_XT
    
  Subroutine Skip(unitno,no)
    implicit none
    integer ( kind = 4 ), intent(in) :: unitno,no
    integer ( kind = 4 ) :: kx
    character (len=5) :: ch
    do kx=1,no
       read(unitno,fmt="(A5)")ch
    end do     
  End Subroutine Skip

  Subroutine Write_Vector_RDP(varname,variable)
    implicit none
    integer (c_int) :: kx,ky,ncol=1,nrow
    character(len=4),intent(in) :: varname
    character(:),allocatable :: wrtfmt
    character(range(ncol)) :: res
    real(c_double),intent(in) :: variable(:)
    open(unit=10,file=varname,form="formatted",&
         &status="replace",action="write")
    write(10,fmt="(A)")"# created by CDW "
    nrow=size(variable)
    write(10,fmt="(A,A)")"# name: ",varname
    write(10,fmt="(A)")"# type: matrix"
    write(10,fmt="(A,i0)")"# rows: ",nrow
    write(10,fmt="(A,i0)")"# columns: ",ncol
    write(res,'(i0)') ncol ! complex values
    wrtfmt="("  //trim(res)//"(e19.12))"
    do ky=1,nrow
       write(10,fmt=wrtfmt)variable(ky)
    end do
    write(10,*)" "
    write(10,*)" "
    close(10)
  End Subroutine Write_Vector_RDP

End Program SV
