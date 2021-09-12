! ------------------------------------------------------------------------------------------!
! 					Name : Naren				Roll no.: 204103110							!
!    						M. Tech in Computational Mechanics								!
!    						Department of Mechanical Engineering							!
!    								IIT- Guwahati											!
!																							!
!    Programming assignment#2																!
!    Advanced Computational Fluid dynamics													!
!																							!
!    Solving Conduction equation (Laplace eq) using conjugate gradient method  				!
!																							!
! ------------------------------------------------------------------------------------------!


subroutine analytical()
    implicit none

    INTEGER::xGrid=100,yGrid=100,i,j,itr_t,k
    double precision:: xDelta, yDelta,x,y,z,lambdaN,lambdaM,sum1,sum2,pi,t
    double precision::alpha=1,xLen=2,yLen=1, thetai=1
    character(len=30)::x3
    character(len=40)::filename

    xDelta=xLen/xGrid
    yDelta=yLen/yGrid
    pi=4.D0*DATAN(1.D0)

    open(unit=3,file='dat_files/tempx5t1.dat',ACTION='write',STATUS='unknown')
    open(unit=4,file='dat_files/tempy5t5.dat',ACTION='write',STATUS='unknown')
    open(unit=5,file='dat_files/tempx5y25.dat',ACTION='write',STATUS='unknown')

    do itr_t=1,11

    write(x3,' (I3.3)') itr_t-1
    filename='dat_files/analytical/'//'time'//trim(x3)//'.dat'
    open(UNIT=2,file=filename,ACTION='write',STATUS='unknown')
    
    t=(itr_t-1.d0)/10

    do i=1,xGrid
        do j=1,yGrid
            x=xDelta*(i-1)
            y=yDelta*(j-1)
            sum1=0
            sum2=0
            do k=0,200
                lambdaN=(2.d0*k+1)*pi/(2*xLen)
                lambdaM=(2.d0*k+1)*pi/(2*yLen)
                sum1=sum1+(((-1.d0)**k)/(lambdaN*xLen))*exp(-alpha*(lambdaN**2)*t)*cos(lambdaN*x)
                sum2=sum2+(((-1.d0)**k)/(lambdaM*yLen))*exp(-alpha*(lambdaM**2)*t)*cos(lambdaM*y)
                !write(*,*) k,lambdaN,xDelta,xLen,(((-1.d0)**k)/(lambdaN*xLen))*exp(-alpha*(lambdaN**2)*t)*cos(lambdaN*x)
                !write(*,*) x,alpha,t
            end do
            z=thetai*4*sum1*sum2
            write(2,*) x,y,z
            if(i==xGrid/4 .and. modulo(j,4)==0 .and. itr_t==2)then
                write(3,*) y,z
            end if
            if(j==yGrid/2 .and. modulo(i,4)==0 .and. itr_t==6)then
                write(4,*) x,z
            end if
            if(i==xGrid/4.and.j==yGrid/4)then
                write(5,*) t,z
            end if
            !exit
        end do
        !exit
        write(2,*)
    end do

    CLOSE(2)
    !exit
    !write(*,*) itr_t-1
    end do
end subroutine analytical

program rk4
    implicit none
    
    INTEGER::xGrid=200,yGrid=100,i,j,itr_t
    double precision,DIMENSION(:,:),ALLOCATABLE::Temp,F1,F2,F3,F4,K1,K2,K3
    double precision::alpha=1,x,y,xDelta,yDelta,pi,t
    double precision::xLen=2,yLen=1,tStep=0.00001
    character(len=30)::x3
    character(len=40)::filename
    
    xDelta=xLen/xGrid
    yDelta=yLen/yGrid
    pi=4.D0*DATAN(1.D0)

    call analytical()

    allocate(Temp(xGrid+2,yGrid+2))
    allocate(F1(xGrid+2,yGrid+2))
    allocate(K1(xGrid+2,yGrid+2))
    allocate(F2(xGrid+2,yGrid+2))
    allocate(K2(xGrid+2,yGrid+2))
    allocate(F3(xGrid+2,yGrid+2))
    allocate(K3(xGrid+2,yGrid+2))
    allocate(F4(xGrid+2,yGrid+2))

    Temp=1

    Temp(xGrid+2,:)=0
    Temp(:,yGrid+2)=0

    F1=Temp
    F2=Temp
    F3=Temp
    F4=Temp
    K1=Temp
    K2=Temp
    K3=Temp
    
    open(unit=3,file='dat_files/rk4tempx5t1.dat',ACTION='write',STATUS='unknown')
    open(unit=4,file='dat_files/rk4tempy5t5.dat',ACTION='write',STATUS='unknown')
    open(unit=5,file='dat_files/rk4tempx5y25.dat',ACTION='write',STATUS='unknown')

    do itr_t=1,100001
        
        if(mod(itr_t-1,10000)==0) then
            write(x3,' (I3.3)') (itr_t-1)/10000
            filename='dat_files/rk4/'//'time'//trim(x3)//'.dat'
            open(UNIT=2,file=filename,ACTION='write',STATUS='unknown')
            do i=2,xGrid+1
                do j=2,yGrid+1
                    x=xDelta*(i-2)
                    y=yDelta*(j-2)
                    write(2,*) x,y,Temp(i,j)
                    if(i==1+xGrid/4 .and. itr_t==10001)then
                        write(3,*) y,Temp(i,j)
                    end if
                    if(j==yGrid/2 .and. itr_t==50001)then
                        write(4,*) x,Temp(i,j)
                    end if
                    !exit
                end do
                !exit
                write(2,*)
            end do
        
            CLOSE(2)  
            !write(*,*) K1(xGrid+2,:)  
        end if

        t=(itr_t-1)*tStep
        Temp(1,1:yGrid+1)=Temp(2,1:yGrid+1)
        Temp(1:xGrid+1,1)=Temp(1:xGrid+1,2)

        write(5,*) t,Temp(1+xGrid/4,1+yGrid/4)
              
    
        F1(2:xGrid+1,2:yGrid+1)=alpha*((Temp(1:xGrid,2:yGrid+1)-2*Temp(2:xGrid+1,2:yGrid+1)+Temp(3:xGrid+2,2:yGrid+1))/xDelta**2+&
            (Temp(2:xGrid+1,1:yGrid)-2*Temp(2:xGrid+1,2:yGrid+1)+Temp(2:xGrid+1,3:yGrid+2))/yDelta**2)
        K1(2:xGrid+1,2:yGrid+1)=Temp(2:xGrid+1,2:yGrid+1)+tStep/2*F1(2:xGrid+1,2:yGrid+1)
        K1(1,1:yGrid+1)=K1(2,1:yGrid+1)
        K1(1:xGrid+1,1)=K1(1:xGrid+1,2)

        F2(2:xGrid+1,2:yGrid+1)=alpha*((K1(1:xGrid,2:yGrid+1)-2*K1(2:xGrid+1,2:yGrid+1)+K1(3:xGrid+2,2:yGrid+1))/xDelta**2+&
            (K1(2:xGrid+1,1:yGrid)-2*K1(2:xGrid+1,2:yGrid+1)+K1(2:xGrid+1,3:yGrid+2))/yDelta**2)
        K2(2:xGrid+1,2:yGrid+1)=Temp(2:xGrid+1,2:yGrid+1)+tStep/2*F2(2:xGrid+1,2:yGrid+1)
        K2(1,1:yGrid+1)=K2(2,1:yGrid+1)
        K2(1:xGrid+1,1)=K2(1:xGrid+1,2)

        F3(2:xGrid+1,2:yGrid+1)=alpha*((K2(1:xGrid,2:yGrid+1)-2*K2(2:xGrid+1,2:yGrid+1)+K2(3:xGrid+2,2:yGrid+1))/xDelta**2+&
            (K2(2:xGrid+1,1:yGrid)-2*K2(2:xGrid+1,2:yGrid+1)+K2(2:xGrid+1,3:yGrid+2))/yDelta**2)
        K3(2:xGrid+1,2:yGrid+1)=Temp(2:xGrid+1,2:yGrid+1)+tStep*F3(2:xGrid+1,2:yGrid+1)
        K3(1,1:yGrid+1)=K3(2,1:yGrid+1)
        K3(1:xGrid+1,1)=K3(1:xGrid+1,2)

        F4(2:xGrid+1,2:yGrid+1)=alpha*((K3(1:xGrid,2:yGrid+1)-2*K3(2:xGrid+1,2:yGrid+1)+K3(3:xGrid+2,2:yGrid+1))/xDelta**2+&
            (K3(2:xGrid+1,1:yGrid)-2*K3(2:xGrid+1,2:yGrid+1)+K3(2:xGrid+1,3:yGrid+2))/yDelta**2)
        Temp(2:xGrid+1,2:yGrid+1)=Temp(2:xGrid+1,2:yGrid+1)+tStep/6*(F1(2:xGrid+1,2:yGrid+1)+2*F2(2:xGrid+1,2:yGrid+1)+&
                        2*F3(2:xGrid+1,2:yGrid+1)+F4(2:xGrid+1,2:yGrid+1))
     !exit   
    end do

    CLOSE(3)
    CLOSE(4)
    CLOSE(5)
end program rk4

