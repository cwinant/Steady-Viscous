Module mod_SparseHB
  use, intrinsic :: iso_c_binding
  implicit none
  Type SparseHB
     integer(c_int),allocatable :: colptr(:),rowind(:)
     real(c_double),allocatable :: data(:)
  End Type SparseHB
  ! Global Sparse Matrices these seem to need to be here
  ! rather than in mod_BuildMatrices
  type(sparseHB) :: LHS,RHS,A
  ! date_time return
  integer (c_int) :: vls(8) 
  real (c_double)  :: st,fn  ! timers

contains
  
  Function AllocateSparseHB(nx,nnz) 
    ! nx=no coluns in full array
    ! nnz=no of non-zero entries in the full array
    integer(c_int), intent(in) :: nx,nnz
    type(sparseHB) :: AllocateSparseHB
    allocate(AllocateSparseHB%colptr(nx+1))
    allocate(AllocateSparseHB%rowind(nnz))
    allocate(AllocateSparseHB%data(nnz))
    AllocateSparseHB%data=0.0d0 ! initialize
  End Function AllocateSparseHB

  Function PHB(row,col)
    ! for col,row PHB return the corresponding position
    ! in the HB structure,
    implicit none
    integer (c_int),intent(in) :: col,row
    integer (c_int) :: k,PHB
       do k = LHS%colptr(col), LHS%colptr(col+1)-1
          if (LHS%rowind(k).ne. row) then
             PHB = 0
          else
             PHB = k
             exit
          end if
       end do
  End Function PHB

  Subroutine PrintSparseHB(M)
    implicit none
    type(sparseHB) :: M
    write(*,fmt="(10i4)")M%colptr
    write(*,fmt="(10i4)")M%rowind
    write(*,fmt="(5(1pe10.3))")M%data
  End Subroutine PrintSparseHB
    
  Subroutine SparseHBPrep(x,t,nx,nu,ne,nnz)
    implicit none
    real(c_double),intent(in) :: x(:,:)
    integer(c_int),intent(in) :: t(:,:),nx,ne,nu
    integer(c_int),intent(out) :: nnz
    integer(c_int) :: temp(2),ke,kx,ky,LSol
    ! to hold temporary values
    integer(c_int),allocatable :: colptr(:),rowind(:)
    integer(c_int),allocatable :: op(:),nop(:),id(:,:),list(:),order(:)
    ! This if for mixed quad-linear TF:
    ! The sparsematrices will be 2*nu+np, minus the no of boundary points
    !  Set up how sparse matrix indexing etc..
    !........................  SparsePrep ......................
    !! Step 1 make an extended id table that covers all points
    !! size(id)=Lsolby60
    LSol=2*nu+nx ! length of solution vector
    !! Step 1a 
    !! allocate and initialize id
    allocate(id(LSol,60));id=-1; ! allocate and initialize to -1
    ! this part determines how many connections to each node
    ! this does u points first, then adds v and p points
    allocate(nop(LSol))
    nop=0
    do ke=1,ne                 !for each element (triangle)
       do kx=1,6                !for each of the six nodes in t
       ! what are the other five? note op is allocated by next statement
          op=t(ke,:)!op=pack(t(ke,1:6),t(ke,1:6)/=t(ke,kx))
          do ky=1,6             ! for each of these other points
             ! do we already have this connection?
             !! op(ky)=node number
             if (all(id(op(ky),:)/=t(ke,kx))) then !no add it to list
                nop(op(ky))=nop(op(ky))+1
                id(op(ky),nop(op(ky)))=t(ke,kx)
                ! mixed quad lin, second matrix
                nop(op(ky))=nop(op(ky))+1
                id(op(ky),nop(op(ky)))=nu+t(ke,kx) 
                if (t(ke,kx).le.nx) then
                   !mixed quad lin third entry
                   nop(op(ky))=nop(op(ky))+1
                   id(op(ky),nop(op(ky)))=2*nu+t(ke,kx) !3nd
                end if
             end if
          end do
       end do
    end do
    id(nu+1:2*nu,:)=id(1:nu,:);nop(nu+1:2*nu)=nop(1:nu)
    !! modify last np entries in id table
    do kx=2*nu+1,LSol
       deallocate(op)
       ! eliminate >2*nu connections
       op=pack(id(kx-2*nu,:),id(kx-2*nu,:)<=2*nu)
       op=pack(op,op/=-1)
       nop(kx)=size(op)
       id(kx,1:nop(kx))=op
    end do
    nnz=sum(nop)
    !! Step 2 determine colptr and rowind
    allocate(colptr(Lsol+1));colptr(1)=1
    allocate(rowind(nnz))
    !! this part makes colptr
    do kx=1,LSol! for each column
       colptr(kx+1)=colptr(kx)+nop(kx)
    end do
    deallocate(op)! re-use op
    !! this part makes rowind (and colind?)
    do kx=1,LSol! for each column
       op=pack(id(kx,:),id(kx,:)/=-1)! for each nonneg entry in that col.
       allocate(order(size(op)),list(size(op)))
       order=(/(ky,ky=1,size(op))/)! instead of bubble sort
       do ke=1,size(op)
          if (ke==1)then
             list=(/(ky,ky=(op(order(1))-1)+1,op(order(1)))/)
          else
             list=(/list,(ky,ky=(op(order(ke))-1)+1,op(order(ke)))/)
          endif
       end do
       rowind(colptr(kx):colptr(kx+1)-1)=list
       deallocate(order,list)
    end do
    LHS=AllocateSparseHB(nx,nnz)
    LHS%colptr=colptr;LHS%rowind=rowind;
    RHS=AllocateSparseHB(nx,nnz)
    RHS%colptr=colptr;RHS%rowind=rowind;
    A=AllocateSparseHB(nx,nnz)
    A%colptr=colptr;A%rowind=rowind;
  End Subroutine SparseHBPrep

    Function HBSubMatrix(Big,keepcol,keeprow) 
    integer(c_int), intent(in) :: keepcol(:),keeprow(:)
    type(sparseHB) :: Big,HBSubMatrix
    
    integer(c_int) :: subnx,subnnz,subkcol,subkrow,thisrow,k,index
    integer(c_int),allocatable :: rows(:)
    real(c_double),allocatable :: data(:)
    subnx=size(keepcol)
    subnnz=0
    HBSubMatrix=AllocateSparseHB(subnx,size(Big%rowind))
    HBSubMatrix%colptr=0;HBSubMatrix%colptr(1)=1
    call date_and_time(VALUES=vls)
    st=60*vls(6)+vls(7)+dble(vls(8))/1000.

    do subkcol=1,subnx
       HBSubMatrix%colptr(subkcol+1)=HBSubMatrix%colptr(subkcol)
       do subkrow=1,size(keeprow)
          index=PHB(keeprow(subkrow),keepcol(subkcol))
          if (index.gt.0) then
             subnnz=subnnz+1
             HBSubMatrix%colptr(subkcol+1)=HBSubMatrix%colptr(subkcol+1)+1
             HBSubMatrix%rowind(subnnz)=subkrow
             HBSubMatrix%data(subnnz)=Big%data(index)
          endif
       end do
    end do
    call date_and_time(VALUES=vls)
    fn=60*vls(6)+vls(7)+float(vls(8))/1000.
    print *,'HBS-subkcol loop , wall clock:',fn-st
    !trim
    HBSubMatrix%rowind=HBSubMatrix%rowind(1:subnnz)
    HBSubMatrix%data=HBSubMatrix%data(1:subnnz)
  End Function HBSubMatrix


  !The following works but is slow
!!$  Function SparseHBSubMatrix(M,kpcol,kprow)
!!$    implicit none
!!$    integer(c_int), intent(in) :: kpcol(:),kprow(:)
!!$    type(sparseHB) :: M,SparseHBSubMatrix
!!$    integer(c_int) :: subncol,subnrow,subnnz,subkcol,subkrow,k1,k2
!!$    subncol=size(kpcol);subnrow=size(kprow);subnnz=0
!!$    ! allocate too big,then trim
!!$    SparseHBSubMatrix=AllocateSparseHB(subncol,size(M%rowind))
!!$    SparseHBSubMatrix%colptr=0;SparseHBSubMatrix%colptr(1)=1
!!$    do subkcol=1,subncol
!!$       SparseHBSubMatrix%colptr(subkcol+1)=SparseHBSubMatrix%colptr(subkcol)
!!$       do subkrow=1,subnrow
!!$          do k1=M%colptr(kpcol(subkcol)),&
!!$               &M%colptr(kpcol(subkcol)+1)-1
!!$             ! for each of the rows in Sparse M in this column
!!$             k2=M%rowind(k1)
!!$
!!$             if (kprow(subkrow)==k2) then
!!$                ! this a row we want to keep
!!$                subnnz=subnnz+1! where the values fits in SparseHBSub
!!$
!!$                SparseHBSubMatrix%colptr(subkcol+1)=&
!!$                     &SparseHBSubMatrix%colptr(subkcol+1)+1
!!$                SparseHBSubMatrix%rowind(subnnz)=subkrow
!!$                SparseHBSubMatrix%data(subnnz)=M%data(k1)
!!$             end if
!!$          end do
!!$       end do
!!$    end do
!!$    SparseHBSubMatrix%rowind=SparseHBSubMatrix%rowind(1:subnnz)
!!$    SparseHBSubMatrix%data=SparseHBSubMatrix%data(1:subnnz)
!!$
!!$  End Function SparseHBSubMatrix

  Subroutine SparseHBMul(M,F,R)
    ! nx=no coluns in full array
    ! nnz=no of non-zero entries in the full array
    real(c_double),intent(in) :: F(:)
    real(c_double),intent(out) :: R(:)
    integer(c_int) :: knz,kcol,krow,k,nnz,nx
    type(sparseHB) :: M
    nx=size(M%colptr)-1;nnz=size(M%rowind)
    R=0;knz=0
    do kcol=1,nx
       do k=M%colptr(kcol),M%colptr(kcol+1)-1
          krow=M%rowind(k)
          knz=knz+1
          R(krow)=R(krow)+M%data(knz)*F(kcol)
       end do
    end do
  End Subroutine SparseHBMul
       
  Subroutine SparseHBToFull(M,FullMat)
    ! nx=no coluns in full array
    ! nnz=no of non-zero entries in the full array
    real(c_double),intent(out) :: FullMat(:,:)
    integer(c_int) :: ncol,nrow,kcol,krow,PrintFull
    type(sparseHB) :: M
    FullMat=0
    ncol=size(M%colptr)-1
    nrow=maxval(M%rowind)
    do kcol=1,ncol
        do krow=M%colptr(kcol),M%colptr(kcol+1)-1
          FullMat(M%rowind(krow),kcol)=M%data(krow)
       end do
    end do
       
  End Subroutine SparseHBToFull

End Module mod_SparseHB
