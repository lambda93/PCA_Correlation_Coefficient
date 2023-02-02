!Program is working in the OpenMP system. This contains upto the step of calculation of eigen-values and eigen_vector
!This program have the additional features to calculate the eigenfuction from the desired independent variable
! new_test
program test
implicit none
real,dimension(:),allocatable::b,variable_1,M,SD
real,dimension(:,:),allocatable::variable_2,EF,U_N
integer,dimension(:),allocatable::cont,k,ek
integer::i,j,i1,l10,j1,I12,L66,l100,nr,num_param,INFO, LWORK,I2,LDA,NEW13
real::r1,dp,CHI,R,ERROR2,SUM1993,SUM2,HP,H0,H0_V,HT,OMEGAM,OMEGAX,HS,CHI2,R01,S,RS
INTEGER,parameter ::lwmax = 5000
real,DIMENSION(:,:),ALLOCATABLE::C,DC,MATRIX,MATRIX_OUT,NEW_MATRIX
real,DIMENSION(:),ALLOCATABLE::W,DW,WORK,COFF
integer,parameter::data_points=38,actual_data_points=38,total_rows=1339                !total_rows=1339
real::dat(3,actual_data_points),U(2,total_rows)
character(len=90)::filename,coffname


l100=10
nr=1339                         !nr=1339
num_param=10
omegam = 0.3
omegax = 0.7
h0_v = 68.5
h0 = 68.5
lda = NUM_PARAM


ALLOCATE(W(NUM_PARAM))
ALLOCATE(WORK(LWMAX))
ALLOCATE(DW(NUM_PARAM))

ALLOCATE(MATRIX(NUM_PARAM,NUM_PARAM+1))
ALLOCATE(MATRIX_OUT(LDA,NUM_PARAM))
ALLOCATE(NEW_MATRIX(NUM_PARAM,NUM_PARAM+1))

allocate(ek(num_param))
allocate(k(num_param))
allocate(cont(num_param))
allocate(b(l100))
allocate(variable_1(num_param+1))   
allocate(variable_2(NUM_PARAM+1,nr))

ALLOCATE(SD(NUM_PARAM))
ALLOCATE(C(NUM_PARAM,NUM_PARAM))
ALLOCATE(DC(LDA,NUM_PARAM))
ALLOCATE(M(NUM_PARAM))

open(0,file='hz_fiducial.dat',status='old')
open(1,file='input.dat',status='old')
read(0,*)dat
read(1,*)u
 close(0)
 close(1)
 
open(69,file='./created/DATA.txt',status='replace')
open(6969,file='./created/COV_MATRIX.txt',status='replace')
OPEN(56,FILE='./created/EIGENVALUES.txt',status='replace')
OPEN(65,FILE='./created/EIGENVECTORS.txt',status='replace')
OPEN(96,FILE='./created/EIGENFUNCTIONS.txt',STATUS='replace')

OPEN(9696,FILE='./created/NEIGENFNS.txt',STATUS='replace')
OPEN(13,FILE='./created/DATA_FINAL.txt',STATUS='replace')
OPEN(1993,FILE='./created/FINAL_FILE.txt',STATUS='replace')

variable_2 = 0
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(NUM_PARAM,OMEGAM,OMEGAX,H0,H0_V,L100,U,dat,variable_2,nr) NUM_THREADS(78)
!$OMP DO

DO L66 = 1,nr
r1 = U(1,L66)                                           !Initial Prior, Upper limit
dp = U(2,L66)                                                  !step size	
I12 = 0
k=0
cont = 0
ek = 0
variable_1 = 0
print*,l66
        do i=1,l100
        b(i) = r1 - i*dp
        end do

        do l10 = 1 ,l100**num_param
                do j=1,num_param
                if(j .ne. 1)then  

                        if(ek(j) .eq. l100**(j-1))then
                        cont(j) = cont(j) + 1 
                        ek(j) = 0
                        end if

        10               if((cont(j-1) - l100*cont(j)) .eq. (l100 - 1))then
                         ek(j) = ek(j) + 1
                         end if

                  k(j) = cont(j-1) - cont(j)*l100 + 1  
                  else
                           if(mod((l10 - 1),l100) .eq. 0 .and. l10 .ge. l100)then
	  	 	   cont(1) = cont(1) + 1
	  	 	   end if
		 	 k(j) = l10 - cont(1)*l100
		end if	
		end do
	
	!* Chi calculation	
	chi=0
		do i=1,data_points        !38              
		r=dat(1,i)
		error2=dat(3,i)	
		sum2=0			
				do j=1,NUM_PARAM			
				sum2=sum2+b(K(J))*(funct(r)**(j-1))
				end do	
			
		HP=sum2
		hs = dat(2,i)	
		chi=chi+((HS-HP)**2)/(error2**2)
		end do

	!* Finding the Minimum chi and the corresponding point in parameter space
	
	I12=I12+1
	 
	 if(I12 == 1)then
	 chi2=chi
	 end if
	
			 if (chi .le. chi2)then 
			 chi2=chi
			 variable_1(1)=chi
	
			  DO J=1,NUM_PARAM
			  variable_1(J+1) = b(K(J))
			  END DO 
				
			 end if
	
	END DO

	variable_2(1,l66) = chi2	
	DO J = 1,NUM_PARAM
	variable_2(J+1,l66) = variable_1(J+1)
	END DO
END DO
!$OMP END DO
!$OMP END PARALLEL

9999 FORMAT( 11(:,5X,F12.7) )

DO L66=1,nr
WRITE(69,*)(variable_2(J,l66),J=1,NUM_PARAM+1)
END DO

print*,'complete the first DR'
!*===================================Creating Covariance Matrix===============================================
 do j=1,NUM_PARAM
 R01=0
 S=0
 	do i=1,nr
 	R01=variable_2(j+1,i)+R01		!finding mean
 	end do
 m(j)=R01/nr			!value of mean

		
 		do i=1,nr
 		S=S+((variable_2(j+1,i)-m(j))**2)/(nr-1)   		!finding variance
 		end do
 		
 SD(j)=S				!value of variance
 C(j,j)=SD(j)
 end do
  
 do I1=1,NUM_PARAM-1
 	do j=1,NUM_PARAM-I1
 	RS=0
 		do i=1,nr
 		RS=RS+variable_2(j+1,i)*variable_2(j+1+I1,i)			!finding covariance matrix
 		end do
 	C(j,j+I1)=(RS-m(j)*m(j+I1))/(nr-1)
 	C(j+I1,j)=C(j,j+I1)
 	end do
 end do

 do i=1,NUM_PARAM
!write(6969,9999)(c(i,j),j=1,NUM_PARAM)
write(6969,*)(c(i,j),j=1,NUM_PARAM)
 end do
 print*,'complete calculating the covariance matrix'
!*=============================Creating The Eigenvectors and Eigenvalues======================================
!open(1893,file='cov_matrix_rd_1-a.txt',status='old')
!read(1893,*)C
! close(1893)

9901 continue

      LWORK = -1
      CALL SSYEV( 'Vectors', 'U', NUM_PARAM, C, LDA, W, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
      CALL SSYEV( 'Vectors', 'U', NUM_PARAM, C, LDA, W, WORK, LWORK, INFO )

!*     Check for convergence.

      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF

!*     Print eigenvalues.

      CALL PRINT_MATRIX( 'Eigenvalues', 1, NUM_PARAM, W, 1, 0)

!*     Print eigenvectors.

DO I2=1,NUM_PARAM	
	DO I1=1,NUM_PARAM
		MATRIX_OUT(I1,I2) = C(I2,I1)
	END DO
END DO

      CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', NUM_PARAM, NUM_PARAM, MATRIX_OUT,LDA,1)

!*================================ARRANGING EIGENVECTORS AS REQUIRED==========================================
DO I2 = 1,NUM_PARAM
	MATRIX(I2,1) = W(I2)
END DO

DO I1=1,NUM_PARAM
	DO I2=1,NUM_PARAM
		MATRIX(I2,I1+1) = C(I1,I2)
	END DO
END DO

CALL DESCENDING(NUM_PARAM,MATRIX,NEW_MATRIX)

DO I1=1,NUM_PARAM
DW(I1)=NEW_MATRIX(I1,1)
END DO

DO I1=1,NUM_PARAM
	DO I2=1,NUM_PARAM
		DC(I2,I1) = NEW_MATRIX(I2,I1+1)
	END DO
END DO

CALL PRINT_MATRIX( 'Eigenvalues', 1, NUM_PARAM, DW, 1,0)
CALL PRINT_MATRIX( 'Eigenvectors (stored column-wise)', NUM_PARAM, NUM_PARAM, DC,LDA,1)

DEALLOCATE(WORK)
DEALLOCATE(variable_2)
print*,'complete calculating the eigenvalues'
!*======================================EIGENFUNCTONS=========================================================

ALLOCATE(EF(NUM_PARAM,actual_data_points))
ALLOCATE(U_N(NUM_PARAM,actual_data_points))

do i=1,data_points
r=dat(1,i)
error2=dat(3,i)
		
		do i1=1,NUM_PARAM
		sum2=0			
			do j=1,NUM_PARAM			
			sum2=sum2+DC(i1,j)*((funct(r))**(j-1))
			end do	
		EF(i1,i)=sum2
		end do
write(96,*)(EF(j,i),j=1,num_param)	

	SUM2 = 0
	do i1=1,num_param
	SUM2 = SUM2 + EF(i1,i)*EF(i1,i)/(error2*error2)
	end do
	
	do i1=1,num_param
	u_n(i1,i) = EF(i1,i)/sqrt(SUM2)
	end do
write(9696,*)(u_n(j,i),j=1,num_param)
end do
print*,'complete calculating the eigenfunction'
!*=====================================FINAL RUN FOR THE CHI SQUARE==========================================
DEALLOCATE(C)
allocate(variable_2(NUM_PARAM+1,nr))

!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(nr,NUM_PARAM,OMEGAM,OMEGAX,H0,H0_V,L100,EF,U,dat,variable_2) NUM_THREADS(78)
!$OMP DO

DO L66 = 1,nr
r1 = U(1,L66)                                            !Initial Prior, Upper limit
dp = U(2,L66)                                            !step size	
I12 = 0
k=0
cont = 0
ek = 0
variable_1 = 0
variable_2 = 0
!print*,'see1',L66
        do i=1,l100
        b(i) = r1 - i*dp
        end do

        do l10 = 1 ,l100**num_param
                do j=1,num_param
                if(j .ne. 1)then 

                        if(ek(j) .eq. l100**(j-1))then
                        cont(j) = cont(j) + 1 
                        ek(j) = 0
                        end if

        100             if((cont(j-1) - l100*cont(j)) .eq. (l100 - 1))then
                        ek(j) = ek(j) + 1
                        end if

                   k(j) = cont(j-1) - cont(j)*l100 + 1
                    else
                           if(mod((l10 - 1),l100) .eq. 0 .and. l10 .ge. l100)then
                            cont(1) = cont(1) + 1
                              end if
                   k(j) = l10 - cont(1)*l100
                   end if
                   end do                !ending loop of j=1,num_param

       !* Chi calculation	
 chi=0
                do i=1,data_points              
                r=dat(1,i)
                error2=dat(3,i)
                 sum2=0
                                do j=1,NUM_PARAM
                                sum2=sum2+b(K(J))*EF(J,I)
                                end do
                HP=sum2
                hs = dat(2,i)
                chi=chi+((HS-HP)**2)/(error2**2)
                end do                                       !i= 1,data_points

                !* Finding the Minimum chi and the corresponding point in parameter space
                !print*,chi	
         I12=I12+1
 
          if(I12 == 1)then
          chi2=chi
          end if

                          if (chi .le. chi2)then 
                          chi2=chi
                           variable_1(1)=chi

                          DO J=1,NUM_PARAM
                          variable_1(J+1) = b(K(J))
                          END DO 
                          
                          end if
                          END DO                             !ending loop of l10

         variable_2(1,l66) = chi2
         DO J = 1,NUM_PARAM
         variable_2(J+1,l66) = variable_1(J+1)
         END DO
!print*,'see2=',l66
END DO
!$OMP END DO
!$OMP END PARALLEL

DO L66=1,nr
WRITE(13,*)(variable_2(J,l66),J=1,NUM_PARAM+1)	
END DO

print*,'complete calculating the second DR'
!*==============================MINIMISING THE FINAL DATA-SET=================================================

 sum2 = variable_2(1,1)

DO L66 = 1,nr
	IF (sum2 .GE. variable_2(1,L66)) THEN
	J=L66
	SUM2=variable_2(1,L66)
	END IF 
END DO
WRITE(1993,*)(variable_2(I,J),I=1,NUM_PARAM + 1)

DO I = 1,NUM_PARAM + 1
variable_1(I)=variable_2(I,J)
END DO

print*,'complete finding the minimum'
!*========================================FINDING THE COEFFICIENTS AND THE FUNCTOINAL FORM==================
NEW13 = 2
DEALLOCATE(variable_2)
ALLOCATE(variable_2(NEW13,data_points))
DO I1 = NUM_PARAM - 2, NUM_PARAM

write (filename, '( "./created/resultant", I1, ".txt" )' )  I1
write (coffname, '( "./created/resultant_coff", I1, ".txt" )' )  I1

open (file=filename,unit=1932,status='replace')
open (file=coffname,unit=1947,status='replace')

ALLOCATE(COFF(I1))
        
        DO J = 1,I1
        COFF(J) = 0
                DO I = 1, NUM_PARAM
                COFF(J) = variable_1(I+1)*DC(I,J) + COFF(J)
                END DO
        END DO

	write(1947,9999)(coff(i),i=1,i1)

        DO I = 1,data_points
        SUM2 = 0
	r = dat(1,I)
                DO J = 1,I1
                SUM2 = SUM2 + COFF(J)*((funct(r))**(j-1))
                END DO
        variable_2(2,I) = SUM2
        variable_2(1,I) = DAT(1,I)       
	END DO
print*,'one step forward'
        DO I=1,data_points
        WRITE(1932,*)(variable_2(J,I),J=1,NEW13)
        END DO

DEALLOCATE(COFF)
 close(1932)
 close(1947)
END DO
print*,'complete the program_calculation of functional form'
!*==============================SUBROUTINES AND FUNCTIONS=====================================================
1001 continue
CONTAINS

	FUNCTION funct(t)
	      real, intent(in) :: t
	      real::funct
	  	funct=t/(1+t)
	end function funct

!*============================================================================================================
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA,new )
      CHARACTER*(*) ::   DESC
      INTEGER       ::   M,N,LDA
      real :: A(LDA,*)

      INTEGER     ::     I93, J93, new

      WRITE(*,*)
      WRITE(*,*) DESC
      
      IF (NEW == 0) THEN
      	DO I93 = 1, M
     	WRITE(56,9998) ( A( I93, J93 ), J93 = 1, N )     !56
      	END DO
	9998 FORMAT(11(:,1X,E13.3))
      ELSE
      	DO I93 = 1, M
      	WRITE(65,9999) ( A( J93, I93 ), J93 = 1, N ) !65	!MAKING THE MATRIX COLUMN-WISE, BUT ACTUALLY IT IS STORED IN ROW-WISE
      	END DO
	9999 FORMAT( 11(:,5X,F12.7) )
      END IF
       
      END subroutine print_matrix

!=============================================================================================================
	SUBROUTINE DESCENDING(N1,X,X_OUT)
    		INTEGER,INTENT(IN)::N1
    		INTEGER:: flag,i93,j93,count
    		real:: temp,X(N1,N1+1)
    		real,INTENT(OUT)::X_OUT(N1,N1+1)

    		count=N1
    		    
    	do
        	flag=0
        	
        	do i93=1,count-1
        	    if (ABS(x(i93,1))<ABS(x(i93+1,1)))then
        	    
        	    	do j93 = 1,NUM_PARAM+1
        	        temp=x(i93,j93)
        	        x(i93,j93)=x(i93+1,j93)
        	        x(i93+1,j93)=temp 
        	    	end do
        	    
        	    flag=1    
        	    end if
        	end do
        	if (flag==0)  then 
		
		X_OUT = X		
		
        	   exit
        	  end if
    	end do		
    	
    END SUBROUTINE DESCENDING

end program test
