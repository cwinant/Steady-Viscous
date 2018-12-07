Module mod_BuildMatrices
  use, intrinsic :: iso_c_binding
  use mod_SparseHB
  implicit none
  ! Element Matrices
  real(c_double) :: Me(6,6),Kexx(6,6),Kexy(6,6),Keyx(6,6),Keyy(6,6)
  real(c_double) :: Dex(3,6),Dey(3,6),Aex(6,6,6),Aey(6,6,6)
  integer(c_int),allocatable :: Dirichlet(:),Interior(:)
  real(c_double),allocatable :: tildeb(:)
  ! date_time return
  integer (c_int) :: vals(8) 
  real (c_double)  :: str,fin  ! timers


contains

  Subroutine BuildMatrices(xu,tu,nx,nu,ne,nnz,N,E,S,W,Soln)
    implicit none
    integer(c_int),intent(in) :: tu(:,:),N(:),E(:),S(:),W(:)
    integer(c_int),intent(in) :: nx,ne,nu
    integer (c_int),intent(inout) :: nnz
    real(c_double),intent(in) :: xu(:,:)
    real(c_double),intent(inout) :: Soln(:)
    !
    ! QuadLin determines element matrices
    ! Me,Kxx,Kxy,Kyx,Kyy,Ax,Ay are all 6 by 6
    ! Dx,Dy are 6 by 3
    ! begin with initial UpsideUp element
    call ElementMatrices(xu(tu(1,1:5:2),:))
    ! for random shape elements, ElementMatrices has to
    ! be called for each element from Global Matrices
    call GlobalMatrices(xu,tu,nx,nu,ne,nnz,N,E,S,W,Soln)
  End Subroutine BuildMatrices

  Subroutine ElementMatrices(x)
    implicit none
    ! Initially called once for UU (canonical) triangle
    ! ElementMatrices are
    ! Me etc..  are 6 by 6
    ! Dex,Dey are 3 by 6
    
    ! Step 1 get polynomial coefficients V(6x6)for quad TF
    ! and L(3x3) for linear TF
    real(c_double),intent(in) :: x(:,:)
    real (c_double) :: theta,a,b,c,R(2,2),R2(4,4)
    ! temp arrays needed to find Kexx etc...
    real (c_double),allocatable :: T1(:,:),T2(:,:),T3(:,:),T4(:,:)
    real (c_double) :: Kxixi(6,6),Kxieta(6,6),Ketaxi(6,6),Ketaeta(6,6)
    real (c_double) :: Dxi(3,6),Deta(3,6)
    real (c_double) :: nodes(6,2),V(6,6),L(3,3),work(6,3)
    real (c_double) :: canodes(6,2)=&
         &reshape((/0.,0.5,1.,0.5,0.,0.,0.,0.,0.,0.5,1.,0.5/),(/6,2/))
    real (c_double),allocatable :: xieta(:,:)
    integer (c_int) :: k,klin=3,kquad=6
    integer (c_int) :: p(6)=(/0,1,0,2,1,0/)
    integer (c_int) :: q(6)=(/0,0,1,0,1,2/)
    theta=atan2(x(2,2)-x(1,2),x(2,1)-x(1,1))
    R=reshape((/cos(theta),sin(theta),-sin(theta),cos(theta)/),(/2,2/))
    xieta=matmul(x,R)
    a=xieta(2,1)-xieta(3,1)
    b=xieta(3,1)-xieta(1,1)
    c=xieta(3,2)-xieta(1,2)
    ! get V matrix of quad coefficients
    work=reshape((/1-canodes(:,1)-canodes(:,2),canodes(:,1),canodes(:,2)/),(/6,3/))
    nodes=matmul(work,reshape((/-b,a,0.d0,0.d0,0.d0,c/),(/3,2/)))
    allocate(T1(6,6),T2(6,6),T3(6,6),T4(6,6))
    T1(:,1)=1.0d0;T1(:,2)=nodes(:,1);T1(:,3)=nodes(:,2)
    T1(:,4)=T1(:,2)**2;T1(:,5)=T1(:,2)*T1(:,3);T1(:,6)=T1(:,3)**2;
    V=inv(T1) ! quad coeffs
    T1(1:3,2)=nodes(1:5:2,1);T1(1:3,3)=nodes(1:5:2,2)
    L=inv(T1(1:3,1:3)) ! linear coefs
    ! Step 2 now evaluate Me
    !% Cowper: FEM2d - eqn 35
    ! This is documented in Steady/Viscous/words/Viscous.tex
    T1=spread(p,1,6)+spread(p,2,6)! what I used to call P
    T2=spread(q,1,6)+spread(q,2,6)! what I used to call Q
    ! gamma(n+1)=factorial(n)
    T3=gamma(dble(T1+T2+3))**(-1) !1/denominator
    T4=gamma(dble(T1+1))*gamma(dble(T2+1))*T3  
    T4=c**(T2+1)*(a**(T1+1)-(-b)**(T2+1))*T4 ! what I call I elsewhere
    Me=matmul(transpose(V),matmul(T4,V))
    !! evaluate Kexx,Kexy,Keyx,Keyy
    !! first evaluate Kxixi, etc then rotate into Kexx etc...
    !arrayP=spread(p,1,6)+spread(p,2,20)!:T1
    T3=T1-2!arrayPm2=arrayP-2
    where(T3.lt.0) T3=0
    Kxixi=matmul(transpose(reshape(p,(/1,6/))),reshape(p,(/1,6/)))
    Kxixi=Kxixi*gamma(dble(T3+T2+3))**(-1)
    Kxixi=gamma(dble(T3+1))*gamma(dble(T2+1))*Kxixi  
    Kxixi=c**(T2+1)*(a**(T3+1)-(-b)**(T3+1))*Kxixi! Ixx
    Kxixi=matmul(transpose(V),matmul(Kxixi,V))  
    T3=T2-2!arrayQm2=arrayQ-2
    where(T3.lt.0) T3=0
    Ketaeta=matmul(transpose(reshape(q,(/1,6/))),reshape(q,(/1,6/)))
    Ketaeta=Ketaeta*gamma(dble(T1+T3+3))**(-1)
    Ketaeta=gamma(dble(T1+1))*gamma(dble(T3+1))*Ketaeta   
    Ketaeta=c**(T3+1)*(a**(T1+1)-(-b)**(T1+1))*Ketaeta! Iyy
    Ketaeta=matmul(transpose(V),matmul(Ketaeta,V))
    T3=T1-1!arrayPm2=arrayP-1
    where(T3.lt.0) T3=0
    T4=T2-1!arrayQm2=arrayQ-1
    where(T4.lt.0) T4=0
    Kxieta=matmul(transpose(reshape(p,(/1,6/))),reshape(q,(/1,6/)))
    Kxieta=Kxieta*gamma(dble(T3+T4+3))**(-1)
    Kxieta=gamma(dble(T3+1))*gamma(dble(T4+1))*Kxieta   
    Kxieta=c**(T4+1)*(a**(T3+1)-(-b)**(T3+1))*Kxieta! Ixy
    Kxieta=matmul(transpose(V),matmul(Kxieta,V))
    Ketaxi=matmul(transpose(reshape(q,(/1,6/))),reshape(p,(/1,6/)))
    Ketaxi=Ketaxi*gamma(dble(T3+T4+3))**(-1)
    Ketaxi=gamma(dble(T3+1))*gamma(dble(T4+1))*Ketaxi   
    Ketaxi=c**(T4+1)*(a**(T3+1)-(-b)**(T3+1))*Ketaxi! Iyx
    Ketaxi=matmul(transpose(V),matmul(Ketaxi,V))
    Kexx=Kxixi*R(1,1)**2+Ketaeta*R(2,1)**2-R(2,1)*R(1,1)*(Kxieta+Ketaxi);
    Keyy=Ketaeta*R(1,1)**2+Kxixi*R(2,1)**2+R(2,1)*R(1,1)*(Kxieta+Ketaxi);
    Kexy=R(1,1)**2*Kxieta-R(2,1)**2*Ketaxi+R(2,1)*R(1,1)*(Kxixi-Ketaeta);
    Keyx=R(1,1)**2*Ketaxi-R(2,1)**2*Kxieta+R(2,1)*R(1,1)*(Kxixi-Ketaeta);
    ! Now De: 3by6
    deallocate(T1,T2,T3,T4)
    allocate(T1(klin,kquad),T2(klin,kquad),T3(klin,kquad),T4(klin,kquad))
    T1=spread(p,1,klin)+spread(p(1:klin),2,kquad)! what I used to call P
    T2=spread(q,1,klin)+spread(q(1:klin),2,kquad)! what I used to call Q
    T3=T1-1
    where(T3.lt.0) T3=0 !Pm1
    T4=T2-1
    where(T4.lt.0) T4=0 !Qm1
    Dxi=gamma(dble(T2+T3+3))**(-1) 
    Dxi=gamma(dble(T3+1))*gamma(dble(T2+1))*Dxi
    Dxi=spread(p,1,klin)*(c**(T2+1)*(a**(T3+1)-(-b)**(T3+1)))*Dxi
    Dxi=matmul(transpose(L),matmul(Dxi,V))! Dxi
    Deta=gamma(dble(T1+T4+3))**(-1) 
    Deta=gamma(dble(T4+1))*gamma(dble(T1+1))*Deta
    Deta=spread(q,1,klin)*(c**(T4+1)*(a**(T1+1)-(-b)**(T1+1)))*Deta! Ixi
    Deta=matmul(transpose(L),matmul(Deta,V))! Deta
    Dex=Dxi*R(1,1)-Deta*R(2,1);
    Dey=Dxi*R(2,1)+Deta*R(1,1)
    ! Now Aex(6,6,6),...
  End Subroutine ElementMatrices

  Subroutine GlobalMatrices(x,t,nx,nu,ne,nnz,N,E,S,W,Soln)
    implicit none
    ! x and t are quad function coord and connectivity
    ! Steps: (1) SparsePrep: set up sparse HB structure:
    ! nx,nnz,colptr and rowind for the 2*nu+np square array
    ! (2) Assembly, fill the full arrays LHS,RHS, and allocate A
    ! (3) BC modify arrays to account for BC
    real(c_double),intent(in) :: x(:,:)
    real(c_double),intent(inout) :: Soln(:)
    integer (c_int),intent(in) :: t(:,:),nx,ne,nu
    integer (c_int),intent(inout) :: nnz
    integer (c_int),intent(in) :: N(:),E(:),S(:),W(:)
    integer (c_int) :: ke,kc,kr,nodes(6),vertex(3),signpm
    type(sparseHB) :: T1

    ! (1) allocate sparse arrays
    call SparseHBPrep(x,t,nx,nu,ne,nnz)
    ! (2) Assemble into global matrix LHS etc..
    ! for now element matrices are known to within sign
    ! This is a new approach to save memory, by assembling
    ! LHS and RHS directly, without storing intermediate
    ! results Kxx,Kxy,etc...  It is still element by element
    !!call PrintSparseHB(LHS)
    

    do ke=1,ne ! for each element
       if(ke.le.ne/2)then
          signpm=+1
       else
          signpm=-1
       end if
       nodes=t(ke,:);vertex=t(ke,1:5:2)
       do kr=1,6
          do kc=1,6
             ! 2*Kexx+Keyy
             LHS%data(PHB(nodes(kr),nodes(kc)))=&
                  &LHS%data(PHB(nodes(kr),nodes(kc)))+&
                  &2*Kexx(kr,kc)+Keyy(kr,kc)
             ! Kexx+2*Keyy
             LHS%data(PHB(nu+nodes(kr),nu+nodes(kc)))=&
                  &LHS%data(PHB(nu+nodes(kr),nu+nodes(kc)))+&
                  &Kexx(kr,kc)+2*Keyy(kr,kc)
             ! Keyx
             LHS%data(PHB(nu+nodes(kr),nodes(kc)))=&
                  &LHS%data(PHB(nu+nodes(kr),nodes(kc)))+Keyx(kr,kc)
             ! Kexy
             LHS%data(PHB(nodes(kr),nu+nodes(kc)))=&
                  &LHS%data(PHB(nodes(kr),nu+nodes(kc)))+Kexy(kr,kc)
             ! bits involving Dex and Dey
             if(kr<4)then ! last np rows
             LHS%data(PHB(2*nu+vertex(kr),nodes(kc)))=&
                  &LHS%data(PHB(2*nu+vertex(kr),nodes(kc)))+&
                  &signpm*(Dex(kr,kc))
             LHS%data(PHB(2*nu+vertex(kr),nu+nodes(kc)))=&
                  &LHS%data(PHB(2*nu+vertex(kr),nu+nodes(kc)))+&
                  &signpm*(Dey(kr,kc))
             ! last np columns: -transpose(Dex) and -transpose(Dey)
             LHS%data(PHB(nodes(kc),2*nu+vertex(kr)))=&
                  &LHS%data(PHB(nodes(kc),2*nu+vertex(kr)))-&
                  &signpm*(Dex(kr,kc))
             LHS%data(PHB(nu+nodes(kc),2*nu+vertex(kr)))=&
                  &LHS%data(PHB(nu+nodes(kc),2*nu+vertex(kr)))-&
                  &signpm*(Dey(kr,kc))
            end if
          end do
       end do
    end do
    ! (3) modify LHS to account for BC
    ! apply BC on Soln(Dirichlet)
    Soln=0.d0
    Soln(W)=1.d0-(x(W,2)/maxval(x(:,2)))**2;
    Dirichlet=(/S,N,W,nu+S,nu+N,nu+W/)
    !write(*,fmt="(a9,14i3)")"Dirichlet",Dirichlet
    Interior=(/(kc,kc=1,2*nu+nx)/)
    do ke=1,size(Dirichlet)
       Interior=pack(Interior,Interior/=Dirichlet(ke))
    end do
    !write(*,fmt="(a9,14i3)")"Interior",Interior
    ! determine tildeb T1 is a temporary  array
    call date_and_time(VALUES=vals)
    str=60*vals(6)+vals(7)+dble(vals(8))/1000.
   
    T1=HBSubMatrix(LHS,Dirichlet,Interior)

    call date_and_time(VALUES=vals)
    fin=60*vals(6)+vals(7)+float(vals(8))/1000.
    print *,'HBSubMatrixDI, wall clock:',fin-str
    allocate(tildeB(size(Interior)))
    call SparseHBMul(T1,Soln(Dirichlet),tildeB)
    tildeb=-tildeb
 
    deallocate(T1%colptr,T1%rowind,T1%data) ! re-use T1
    call date_and_time(VALUES=vals)
    str=60*vals(6)+vals(7)+dble(vals(8))/1000.

    T1=HBSubMatrix(LHS,Interior,Interior)

    call date_and_time(VALUES=vals)
    fin=60*vals(6)+vals(7)+float(vals(8))/1000.
    print *,'HBSubMatrixII, wall clock:',fin-str
    ! replace LHS with Int,Int array
    deallocate(LHS%colptr,LHS%rowind,LHS%data)
    kc=size(T1%colptr)-1;kr=size(T1%rowind)
    !print *,'GB: kc=',kc,' kr=',kr
    LHS=AllocateSparseHB(kc,kr)
    LHS%colptr=T1%colptr
    LHS%rowind=T1%rowind
    LHS%data=T1%data
  End Subroutine GlobalMatrices

  Subroutine PrintArrayRDP(arrayK)
    implicit none
    real(c_double),intent(in) :: arrayK(:,:)
    integer(c_int) :: kk,temp(2)
    temp=shape(arrayK)
    do kk=1,temp(1)
       !print *,kk
       write(*,fmt="(14(1pe9.2))")arrayK(kk,:)
    end do
  End Subroutine PrintArrayRDP

  Function inv(A) result(Ainv)
    real(c_double), dimension(:,:), intent(in) :: A
    real(c_double), dimension(size(A,1),size(A,2)) :: Ainv
    
    real(c_double), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer(c_int), dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info
    
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)
    
    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if
    
    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)
    
    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  End Function inv
End Module mod_BuildMatrices
  
