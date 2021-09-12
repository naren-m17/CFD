! ------------------------------------------------------------------------------------------!
! 										Name : Naren										!
!    						M. Tech in Computational Mechanics								!
!    						Department of Mechanical Engineering							!
!    								IIT- Guwahati											!
!																							!
!    Flow in 2D lid-driven cavity															!
!																							!
! ------------------------------------------------------------------------------------------!

!create helper function 
MODULE helper

 IMPLICIT none
 
 CONTAINS
 
 !calculate sum of error between two iterationn (norm_infinity)
 FUNCTION errorCal(imax, jmax, val, val1) RESULT (error)
 	INTEGER							::imax, jmax,i,j
 	REAL							::error
 	REAL, DIMENSION(:,:),INTENT(IN)	::val, val1
 	
 	error = 0		!intilize to zero
 	DO i = 1,imax
 		DO j = 1,jmax
 			error = error + ABS(val1(i,j) - val(i,j))
 		END DO
 	END DO
 	
 END FUNCTION errorCal
 
 !print data to the data in matrix to different file
 !fresh file will be created (all previous values will be earsed)
 SUBROUTINE printMatrixFile(mat, imax, jmax, fileName, deltaX, deltaY)
 	INTEGER						:: i,j
 	INTEGER						:: imax, jmax
 	REAL						:: deltaX, deltaY
 	REAL, DIMENSION(:, :) 		:: mat				!matrix to be printed
 	CHARACTER(LEN = 20)			:: fileName			!filename into which data to be writen
 	
 	!open file with filename
 	OPEN(1, FILE = fileName, STATUS = 'unknown')
 	
 	DO i = 1,imax
 		DO j = 1, jmax
	 		WRITE(1,*) (i - 1)*deltaX, (j - 1)*deltaY, mat(i,j)
 		END DO
 		WRITE(1,*)
 	END DO
 	
 	CLOSE(1)
 	
 END SUBROUTINE printMatrixFile
 
END MODULE helper

!main program to executed
PROGRAM  mainlidCavity
	!call function explicitly
	USE helper																							
	!disable automatic variable initilization
    IMPLICIT  NONE														
    
    !intilize all variables
	REAL								:: L = 1, W = 1, U = 1, tol = 1e-6 ,deltaX, deltaY, beta, error	
	REAL, DIMENSION (:), ALLOCATABLE	:: uC, vC						!u at center and v at center
	!all matrix to store the data
    REAL, DIMENSION (:,:), ALLOCATABLE  :: phiOld, omegaOld, phiNew, omegaNew, phiTemp, omegaTemp		
    INTEGER 						    :: xmax = 129, ymax = 129, Re = 400, i, j, step = 1
    CHARACTER(LEN = 20)					:: fileName								!file name
    
    !calculate delta values and beta parameters
    deltaX = L/(xmax - 1)
    deltaY = W/(ymax - 1)
    beta = deltaX/deltaY
    
    !allocate all values to store data
    ALLOCATE(phiOld(xmax, ymax))
    ALLOCATE(omegaOld(xmax, ymax))
    ALLOCATE(phiNew(xmax, ymax))
    ALLOCATE(omegaNew(xmax, ymax))
    ALLOCATE(phiTemp(xmax-2, ymax-2))
    ALLOCATE(omegaTemp(xmax-2, ymax-2))
    
    !set all elements to zero
	phiNew = 0
	omegaNew = 0
	
	!set boundary condition
	omegaNew(1, 2:ymax-1) = -2*phiNew(2,2:ymax-1)/deltaX**2
	omegaNew(xmax, 2:ymax-1) = -2*phiNew(xmax-1,2:ymax-1)/deltaX**2
	omegaNew(2:xmax-1, 1) = -2*phiNew(2:xmax-1, 2)/deltaY**2
	omegaNew(2:xmax-1, ymax) = -1*(2*phiNew(2:xmax-1,ymax-1) + 2*deltaY*U)/deltaY**2
	
	!outer loop over inner loop of psi and omega
	outer:DO WHILE(errorCal(xmax, ymax, omegaOld, omegaNew)/(SUM(ABS(omegaNew))) > tol)
		!store previous iteration values
		phiOld = phiNew
		omegaOld = omegaNew
		
		!set i(psi step counter) and j(omega step counter) to zero
		i=0
		j=0
		
		!set different value to run the loop
		phiTemp = 10
		
		!inner loop of psi
		innerPhi: DO WHILE(errorCal(xmax-2,ymax-2, phiNew(2:xmax-1, 2:ymax-1), phiTemp)/(SUM(ABS(phiTemp))) > tol)
			!store previous iteration value
			if(i .NE. 0) then
				phiNew(2:xmax-1, 2:ymax-1) = phiTemp
			end if
			
			!calculation of new psi value
			phiTemp = (phiNew(3:xmax,2:ymax-1) + phiNew(1:xmax-2, 2:ymax-1) + beta**2*(phiNew(2:xmax-1, 3:ymax) + & 
			phiNew(2:xmax-1,1:ymax-2)) + deltaX**2*omegaNew(2:xmax-1,2:ymax-1))/(2*(1+beta**2))

			i = i +1							!increase step countor
			
			!exit the loop if max value is satisfied
			IF (i > 20) THEN
				EXIT innerPhi
			END IF
		END DO innerPhi
		
		!set temp value
		omegaTemp = 10
		
		!inner loop of omega
		innerOmega: DO WHILE(errorCal(xmax-2,ymax-2, omegaNew(2:xmax-1, 2:ymax-1), omegaTemp)/(SUM(ABS(omegaTemp))) > tol)
			!store previous iteration value
			if(j .ne. 0) then
				omegaNew(2:xmax-1, 2:ymax-1) = omegaTemp
			end if
			
			!calculation next iteration omega value		
			omegaTemp = (omegaNew(3:xmax,2:ymax-1) + omegaNew(1:xmax-2,2:ymax-1) + &
			beta**2*(omegaNew(2:xmax-1,3:ymax)  + omegaNew(2:xmax,1:ymax-2)) - &
			beta/4*Re*(omegaNew(3:xmax, 2:ymax-1) - omegaNew(1:xmax-2,2:ymax-1))*(phiNew(2:xmax-1,3:ymax) - phiNew(2:xmax-1,1:ymax-2)) &
			+ beta*Re/4*(omegaNew(2:xmax-1,3:ymax) - omegaNew(2:xmax-1,1:ymax-2))*(phiNew(3:xmax,2:ymax-1)-& 
			phiNew(1:xmax-2,2:ymax-1)))/(2*(1+beta**2))
			
			j = j + 1								!increase step countour
			
			!exit loop of max is satisfied
			IF (j > 3) THEN
				EXIT innerOmega
			END IF
		END DO innerOmega
		
		!update boundary values omega vslue
		omegaNew(1,2:ymax-1) = -2*phiNew(2,2:ymax-1)/deltaX**2
		omegaNew(xmax, 2:ymax-1) = -2*phiNew(xmax-1, 2:ymax-1)/deltaX**2
		omegaNew(2:xmax-1,1) = -2*phiNew(2:xmax-1,2)/deltaY**2
		omegaNew(2:xmax-1,ymax) = -1*(2*phiNew(2:xmax-1,ymax-1)+2*deltaY*U)/deltaY**2
		
		!print *, errorCal(xmax, ymax, omegaOld, omegaNew)/ABS(SUM(omegaNew)), i, j  !for testing
		
		!increase step countor
		step = step + 1
	END DO outer

	!calculate final relative error
	error = errorCal(xmax, ymax, omegaOld, omegaNew)/ABS(SUM(omegaNew))

	!store psi vale to phiValue.dat file
	fileName = 'phiValue.dat'
	CALL printMatrixFile(phiNew, xmax, ymax, fileName, deltaX, deltaY)
	!store omega value to file
	fileName = 'omegaValue.dat'
	CALL printMatrixFile(omegaNew, xmax, ymax, fileName, deltaX, deltaY)		
	
	!deallocate all variables that are used in loop
	DEALLOCATE(omegaOld)
	DEALLOCATE(phiOld)	
	DEALLOCATE(omegaTemp)
	DEALLOCATE(phiTemp)
	
	!store velocity value into txt file
	OPEN(2, FILE = 'velocity.txt', STATUS = 'unknown')								!open velocity.txt
	!allocate velocity vectors 
	ALLOCATE(uC(ymax))
	ALLOCATE(vC(xmax))

	!calculate first and last element using firward and backward difference
	uC(1) = (phiNew(65, 2) - phiNew(65, 1))/deltaY
	uC(ymax) = (phiNew(65, ymax) - phiNew(65, ymax - 1))/deltaY
	!start storing u- velocity values
	WRITE(2,*) "U velocity along vetical line passing through geometric center"
	WRITE(2,*) 1, uC(1)
	!calculate and write all u-values along vertical centerline
	DO i = 2, ymax -1
		uC(2:ymax) = (phiNew(65, i+1) - phiNew(65, i-1))/(2*deltaY)
		WRITE(2,*) i, uC(i)
	END DO
	WRITE(2,*) ymax, uC(ymax)
	
	
	WRITE(2,*)
	!start storing v- velocity values
	WRITE(2,*) "V velocity along horizontal line passing thorough geometric center" 
	!calculate first and last element using firward and backward difference
	vC(1) = -1*(phiNew(2, 65) - phiNew(1, 65))/deltaX
	vC(xmax) = -1*(phiNew(xmax, 65) - phiNew(xmax - 1, 65))/deltaX
	WRITE(2,*) 1, vC(1)
	!calculate and write all v-values along horizontal centerline
	DO j = 2, ymax -1
		vC(2:xmax) = -1*(phiNew(j+1, 65) - phiNew(j-1, 65))/(2*deltaX)
		WRITE(2,*) j, vC(j)
	END DO
	WRITE(2,*) xmax, vC(xmax)
	
	CLOSE(2)
	
	!deallocate all values
	DEALLOCATE(phiNew)	
	DEALLOCATE(omegaNew)
	DEALLOCATE(uC)
	DEALLOCATE(vC)
	
	PRINT*, "Program completed with final relative error ", error
	PRINT*, "With number of steps (outer loop)", step
 END PROGRAM  mainlidCavity
