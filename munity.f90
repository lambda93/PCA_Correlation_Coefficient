module mpi_mod
implicit none
real::dat(3,38),U(2,1339)
integer,parameter::l100=10,num_param=7,nr=1339,nr1=1339,data_points=38,actual_data_points=38,total_rows = 1339
integer::i,j
contains
subroutine DR(vector_receive,avgrows,f_t,output)
integer,intent(in)::vector_receive(2),avgrows,f_t
real,intent(out)::output(num_param+1,1000)
real,dimension(:),allocatable::b,FINAL_MATRIX
integer,dimension(:),allocatable::ek,k,cont
integer::int11,int12,id_number,x1,y1,i12,i,j,l10,l66
real::chi,r,chi2,error2,sum2,hp,hs,r1,dp

allocate(ek(num_param))
allocate(k(num_param))
allocate(cont(num_param))
allocate(b(l100))
allocate(FINAL_MATRIX(num_param+1))   
		int11=vector_receive(1)
		int12=vector_receive(2)
		id_number = (int11 - 1)/avgrows

		if(id_number == f_t)then
  		  x1=1
		  y1= (int12 - int11) + 1
		else
		  x1 = 1
		  y1 = avgrows
		end if 

do l66=x1,y1
r1=u(1,int11 + l66 -1)
dp=u(2,int11 + l66 -1)

I12 = 0
k=0
cont = 0
ek = 0
FINAL_MATRIX = 0
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
		do i=1,data_points              
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
			 FINAL_MATRIX(1)=chi
	
			  DO J=1,NUM_PARAM
			  FINAL_MATRIX(J+1) = b(K(J))
			  END DO 
				
			 end if
	
	END DO
	OUTPUT(1,l66) = chi2
	
	DO J = 1,NUM_PARAM
	OUTPUT(J+1,l66) = FINAL_MATRIX(J+1)
	END DO
END DO
deallocate(ek)
deallocate(k)
deallocate(cont)
deallocate(b)
deallocate(FINAL_MATRIX)   

END SUBROUTINE DR
!===============================================================================================================
subroutine DR_final(vector_receive,e_vector,avgrows,f_t,output)
integer,intent(in)::vector_receive(2),avgrows,f_t
real,intent(in)::e_vector(num_param,num_param)
real,intent(out)::output(num_param+1,1000)
real,dimension(:),allocatable::b,FINAL_MATRIX
integer,dimension(:),allocatable::ek,k,cont
integer::int11,int12,id_number,x1,y1,i12,i,j,l10,l66,i1,j1
real::chi,r,chi2,error2,sum2,hp,hs,r1,dp,coff(num_param)

allocate(ek(num_param))
allocate(k(num_param))
allocate(cont(num_param))
allocate(b(l100))
allocate(FINAL_MATRIX(num_param+1))   
		int11=vector_receive(1)
		int12=vector_receive(2)
		id_number = (int11 - 1)/avgrows

		if(id_number == f_t)then
  		  x1=1
		  y1= (int12 - int11) + 1
		else
		  x1 = 1
		  y1 = avgrows
		end if 

do l66=x1,y1
r1=u(1,int11 + l66 -1)
dp=u(2,int11 + l66 -1)

I12 = 0
k=0
cont = 0
ek = 0
FINAL_MATRIX = 0
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
		do i=1,data_points              
		r=dat(1,i)
		error2=dat(3,i)	
		sum2=0

	do j1=1,num_param
	coff(j1) = 0
	  do i1=1,num_param
	  coff(j1) = coff(j1) + b(k(i1))*e_vector(i1,j1)
	  end do
	end do
			
				do j=1,NUM_PARAM			
				sum2=sum2+coff(j)*(funct(r)**(j-1))
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
			 FINAL_MATRIX(1)=chi
	
			  DO J=1,NUM_PARAM
			  FINAL_MATRIX(J+1) = b(K(J))
			  END DO 
				
			 end if
	
	END DO
	OUTPUT(1,l66) = chi2
	
	DO J = 1,NUM_PARAM
	OUTPUT(J+1,l66) = FINAL_MATRIX(J+1)
	END DO
END DO
deallocate(ek)
deallocate(k)
deallocate(cont)
deallocate(b)
deallocate(FINAL_MATRIX)   

END SUBROUTINE DR_final

!==================================================================================================================
      SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA0,new )
      CHARACTER*(*) ::   DESC
      INTEGER       ::   M,N,LDA0
      real :: A(LDA0,*)

      INTEGER     ::     I93, J93, new

      !WRITE(*,*)
      !WRITE(*,*) DESC
      
      IF (NEW == 0) THEN
      	DO I93 = 1, M
     	WRITE(56,9998) ( A( I93, J93 ), J93 = 1, N )
      	END DO
	9998 FORMAT(11(:,1X,E13.3))
      ELSE
      	DO I93 = 1, M
      	WRITE(65,9999) ( A( J93, I93 ), J93 = 1, N )	!MAKING THE MATRIX COLUMN-WISE, BUT ACTUALLY IT IS STORED IN ROW-WISE
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



function funct(t)
      real:: t,funct
  	funct=t
end function funct

end module mpi_mod


program test
use mpi
use mpi_mod
implicit none
character (len=90) :: filename,coffname

real::pkii(num_param+1,nr),FINAL_WORKSPACE(num_param+1,1000),FINAL_WORKSPACE2(num_param+1,1000)
integer,parameter::send_data_tag=1001,return_data_tag=1000,max_rows = 100000
integer:: my_id,root_process, ierr, status(MPI_STATUS_SIZE),num_process,sender,vector2(2),vector(2),vector3(2),final_thread
integer::an_id,k,end_1,end_2,avg_rows,new1,new2,new3,new4,size_of_FINAL_WORKSPACE2,size_of_the_FINAL_WORKSPACE2
integer,dimension(:),allocatable::v

integer::I1,I2,J1
real::R01,S,RS,M(NUM_PARAM),C(NUM_PARAM,NUM_PARAM),SD(NUM_PARAM)

integer,parameter::lda=num_param,lwmax=5000

INTEGER::INFO,LWORk,i9,j9
real::MATRIX(NUM_PARAM,NUM_PARAM+1),NEW_MATRIX(NUM_PARAM,NUM_PARAM+1),MATRIX_OUT(LDA,NUM_PARAM), &
                          DC(NUM_PARAM,NUM_PARAM),W(NUM_PARAM),WORK(LWMAX),DW(NUM_PARAM) 
                                
real::EF(NUM_PARAM,actual_data_points),U_N(NUM_PARAM,actual_data_points) 
REAL::SUM29,ERR
integer::size_of_eigenvector
integer,parameter::send_tag=2002,return_tag=2000
REAL::RECEIVE_VECTOR(NUM_PARAM,NUM_PARAM)

real::sum13,coefficient(num_param+1),r
real,dimension(:),allocatable::coffee

open(0,file='hz_fiducial.dat',status='old')
open(1,file='input.dat',status='old')
read(0,*)dat
read(1,*)u
 close(0)
 close(1)
 
open(69,file='MDATA_fn3.txt',status='replace')
open(691,file='MCOV_MATRIX_fn3.txt',status='replace')
open(56,file='MEVALUE_fn3.txt',status='replace')
open(65,file='MEVECTOR_fn3.txt',status='replace')
open(96,file='MEIGENFUCNTION_fn3.txt',status='replace')
open(1857,file='MFINAL_DATA_fn3.txt',status='replace')
open(1993,file='MMinimum_chi_fn3.txt',status='replace')

size_of_eigenvector = num_param*num_param

root_process=0
!############################################START OF PARALLELIZATION#######################################################
call MPI_INIT (ierr)
	call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
        call MPI_COMM_SIZE (MPI_COMM_WORLD, num_process, ierr)
allocate(v(num_process*2))
final_thread = num_process - 1 
avg_rows = total_rows/num_process       
if(my_id == root_process)then

	avg_rows = total_rows/num_process
	v(1)=1
	v(2)=avg_rows
	
   			 do an_id=1,num_process-1
			 end_1 = avg_rows*an_id + 1
			 end_2 = end_1 + avg_rows -1
 
 					if(an_id .eq. num_process-1)then
 					end_2 = total_rows
 					end if 

			v(2*an_id + 1)=end_1
			v(2*an_id + 2)=end_2

			vector(1)=v(2*an_id + 1)
			vector(2)=v(2*an_id + 2)
			call MPI_SEND(vector,2,MPI_INT,an_id,send_data_tag,MPI_COMM_WORLD,ierr)
			end do	
!################################## JOB IN THE ROOT_PROCESSOR #######################################################
vector(1)=v(1)
vector(2)=v(2)
call DR(vector,avg_rows,final_thread,FINAL_WORKSPACE)

do new2=1,(num_param+1)
	do new1=1,avg_rows
	pkii(new2,new1)=FINAL_WORKSPACE(new2,new1)
	end do
end do

!**************************************	JOB DONE BY OTHER PROCESSORS *******************************************************
do an_id = 1, num_process - 1
	
SIZE_OF_THE_FINAL_WORKSPACE2 = (1-V(2*an_id + 1)+V(2*an_id + 2))*(num_param+1)
call MPI_RECV( FINAL_WORKSPACE2, SIZE_OF_THE_FINAL_WORKSPACE2, MPI_REAL, an_id,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

do new2=1,(num_param+1)
	do new1=V(2*an_id + 1),V(2*an_id + 2)
	pkii(new2,new1)=FINAL_WORKSPACE2(new2,new1 - an_id*avg_rows)
	end do
end do

end do


!***************************************************************************************************************************
	9999 FORMAT( 11(:,5X,F12.7) )
DO NEW3=1,TOTAL_ROWs
WRITE(69,9999)(PKII(NEW4,NEW3),NEW4=1,num_param+1)
END DO
!*===================================Creating Covariance Matrix===============================================

 do j=1,NUM_PARAM
 R01=0
 S=0
	do i=1,nr1
 	R01=PKII(j+1,i)+R01		!finding mean
 	end do
 m(j)=R01/nr1			!value of mean

		
 		do i=1,NR1
 		S=S+((PKII(j+1,i)-m(j))**2)/(nr1-1)   		!finding variance
 		end do
 		
 SD(j)=S				!value of variance
 C(j,j)=SD(j)
 end do
  
 do I1=1,NUM_PARAM-1
 	do j=1,NUM_PARAM-I1
 	RS=0
 		do i=1,nr1
 		RS=RS+PKII(j+1,i)*PKII(j+1+I1,i)			!finding covariance matrix
 		end do
 	C(j,j+I1)=(RS-m(j)*m(j+I1))/(nr1-1)
 	C(j+I1,j)=C(j,j+I1)
 	end do
 end do

 do i=1,NUM_PARAM
 write(691,9999)(c(i,j),j=1,NUM_PARAM)
 end do
 !print*,'complete calculating the covariance matrix'
!##########################eigenvalue and eigenvector calculation####################################################

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

!      CALL PRINT_MATRIX( 'Eigenvalues', 1, NUM_PARAM, W, 1, 0)

!*     Print eigenvectors.

DO J9=1,NUM_PARAM	
	DO I9=1,NUM_PARAM
		MATRIX_OUT(I9,J9) = C(J9,I9)
	END DO
END DO

 !     CALL PRINT_MATRIX( 'Eigenvectors (stored columnwise)', NUM_PARAM, NUM_PARAM, MATRIX_OUT,LDA,1)

!*================================ARRANGING EIGENVECTORS AS REQUIRED==========================================
DO J9 = 1,NUM_PARAM
	MATRIX(J9,1) = W(J9)
END DO

DO I9=1,NUM_PARAM
	DO J9=1,NUM_PARAM
		MATRIX(J9,I9+1) = C(I9,J9)
	END DO
END DO

CALL DESCENDING(NUM_PARAM,MATRIX,NEW_MATRIX)

DO I9=1,NUM_PARAM
DW(I9)=NEW_MATRIX(I9,1)
END DO

DO I9=1,NUM_PARAM
	DO J9=1,NUM_PARAM
		DC(J9,I9) = NEW_MATRIX(J9,I9+1)
	END DO
END DO

CALL PRINT_MATRIX( 'Eigenvalues', 1, NUM_PARAM, DW, 1,0)
CALL PRINT_MATRIX( 'Eigenvectors (stored column-wise)', NUM_PARAM, NUM_PARAM, DC,LDA,1)

1001 continue
!######################################EIGENFUCNTIONS#################################################################

do i=1,data_points
ERR=dat(3,i)
		
		do i1=1,NUM_PARAM
		sum29=0			
			do j=1,NUM_PARAM			
			sum29=sum29+DC(i1,j)*((funct(dat(1,i)))**(j-1))
			end do	
		EF(i1,i)=sum29
		end do
write(96,9999)(EF(j,i),j=1,num_param)	

	SUM29 = 0
	do i1=1,num_param
	SUM29 = SUM29 + EF(i1,i)*EF(i1,i)/(err*err)
	end do
	
	do i1=1,num_param
	u_n(i1,i) = EF(i1,i)/sqrt(SUM29)
	end do
write(9696,9999)(u_n(j,i),j=1,num_param)
end do
!###################################################DR_FINAL#######################################################
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)

v = 0
vector = 0
	avg_rows = total_rows/num_process
	v(1)=1
	v(2)=avg_rows
	
   			 do an_id=1,num_process-1
			 end_1 = avg_rows*an_id + 1
			 end_2 = end_1 + avg_rows -1
 
 					if(an_id .eq. num_process-1)then
 					end_2 = total_rows
 					end if 

			v(2*an_id + 1)=end_1
			v(2*an_id + 2)=end_2

			vector(1)=v(2*an_id + 1)
			vector(2)=v(2*an_id + 2)
       call MPI_SEND(vector,2,MPI_INT,an_id,send_tag,MPI_COMM_WORLD,ierr)
       call MPI_SEND(dc,SIZE_OF_EIGENVECTOR,MPI_REAL,an_id,send_tag,MPI_COMM_WORLD,ierr)
                	end do	
!################################## JOB IN THE ROOT_PROCESSOR #######################################################
vector(1)=v(1)
vector(2)=v(2)

FINAL_WORKSPACE = 0

call DR_final(vector,dc,avg_rows,final_thread,FINAL_WORKSPACE)

do new2=1,(num_param+1)
	do new1=1,avg_rows
	pkii(new2,new1)=FINAL_WORKSPACE(new2,new1)
	end do
end do

!**************************************	JOB DONE BY OTHER PROCESSORS *******************************************************
do an_id = 1, num_process - 1

	SIZE_OF_THE_FINAL_WORKSPACE2 = (1-V(2*an_id + 1)+V(2*an_id + 2))*(num_param+1)
call MPI_RECV( FINAL_WORKSPACE2, SIZE_OF_THE_FINAL_WORKSPACE2, MPI_REAL, an_id,return_tag, MPI_COMM_WORLD, status, ierr)

do new2=1,(num_param+1)
	do new1=V(2*an_id + 1),V(2*an_id + 2)
	pkii(new2,new1)=FINAL_WORKSPACE2(new2,new1 - an_id*avg_rows)
	end do
end do

end do

!***************************************************************************************************************************
DO NEW3=1,TOTAL_ROWs
WRITE(1857,9999)(PKII(NEW4,NEW3),NEW4=1,num_param+1)
END DO
!*==============================MINIMISING THE FINAL DATA-SET=================================================

 sum13 = PKII(1,1)
DO i2 = 1,nr1
	IF (sum13 .GE. PKII(1,i2)) THEN
	J=i2
	SUM13=PKII(1,i2)
	END IF 
END DO
WRITE(*,*)j
WRITE(1993,*)(PKII(I,J),I=1,NUM_PARAM + 1)

DO I = 1,NUM_PARAM + 1
coefficient(I)=PKII(I,J)
END DO
!*========================================FINDING THE COEFFICIENTS AND THE FUNCTOINAL FORM==================
DO I1 = NUM_PARAM - 2, NUM_PARAM

write (filename, '( "/home/ph15038/algo/created/Mresultant_fn3", I1, ".txt" )' )  I1
write (coffname, '( "/home/ph15038/algo/created/Mresultant_fn3_coff", I1, ".txt" )' )  I1

open (file=filename,unit=1932,status='replace')
open (file=coffname,unit=1947,status='replace')

ALLOCATE(COFFEE(I1))
        
        DO J = 1,I1
        COFFEE(J) = 0
                DO I = 1, NUM_PARAM
                COFFEE(J) = COEFFICIENT(I+1)*DC(I,J) + COFFEE(J)
                END DO
        END DO

	write(1947,9999)(COFFEE(i),i=1,i1)

        DO I = 1,data_points
        SUM13 = 0
	r = dat(1,I)
                DO J = 1,I1
                SUM13 = SUM13 + COFFEE(J)*((funct(r))**(j-1))
                END DO
        PKII(2,I) = SUM13
        PKII(1,I) = DAT(1,I)       
	END DO

        DO I=1,data_points
        WRITE(1932,9999)(PKII(J,I),J=1,2)
        END DO
close(1947)
close(1932)
DEALLOCATE(COFFEE)
END DO
!################################### JOBS ASSIGNED BY THE ROOT TO OTHER PROCESSORS #############################################

else
 	 CALL MPI_RECV ( vector2,2, MPI_INT, root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
         size_of_FINAL_WORKSPACE2 = (num_param+1)*(-vector2(1) + vector2(2) + 1 ) 
                 call DR(vector2,avg_rows,final_thread,FINAL_WORKSPACE2)  	 
	 CALL MPI_SEND( FINAL_WORKSPACE2, size_of_FINAL_WORKSPACE2 , MPI_REAL, root_process, return_data_tag, MPI_COMM_WORLD, ierr)	 
!open(191719,file='/home/ph15038/Hubble/output/st.txt',status='replace')

!=============================================DR_final(second term)======================================================
CALL MPI_BARRIER(MPI_COMM_WORLD,IERR)
 	 CALL MPI_RECV ( vector2,2, MPI_INT, root_process, SEND_TAG, MPI_COMM_WORLD, status, ierr)
 	 CALL MPI_RECV ( RECEIVE_VECTOR,SIZE_OF_EIGENVECTOR, MPI_REAL, root_process, SEND_TAG, MPI_COMM_WORLD, status, ierr)
 	size_of_FINAL_WORKSPACE2 = (num_param+1)*(-vector2(1) + vector2(2) + 1 ) 
         FINAL_WORKSPACE2 = 0
        call DR_final(vector2,RECEIVE_VECTOR,avg_rows,final_thread,FINAL_WORKSPACE2)  	 

	 CALL MPI_SEND( FINAL_WORKSPACE2, size_of_FINAL_WORKSPACE2 , MPI_REAL, root_process, return_tag, MPI_COMM_WORLD, ierr)	 

end if

call MPI_FINALIZE(ierr)
!################################################## *END* ################################################################ 
STOP
end program test
