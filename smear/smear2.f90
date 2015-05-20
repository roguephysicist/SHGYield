MODULE SMEAR_MOD
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.d0)
  REAL(DP), PARAMETER :: pi = 3.14159265358979
!!!  REAL(DP), ALLOCATABLE :: arr(:), smear(:), x(:)
  REAL(DP), ALLOCATABLE :: arr(:,:), x(:)
CONTAINS
!========================================================  
  SUBROUTINE READFILE(FNAME,nMax,LNCOL)
!========================================================
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: fname
    INTEGER, INTENT(INOUT) :: nMax
    INTEGER, INTENT(IN) :: LNCOL
    REAL(DP), ALLOCATABLE :: array_tmpX(:), array_tmp(:,:)
    !REAL(DP) :: tot, energy, h, intmpx, intmpf(6)
    REAL(DP) :: tot, energy, h, intmpx, intmpf(50)
    INTEGER :: i, ioerror
    INTEGER :: oldSize, newSize
    nMax = 1
    
    ALLOCATE( arr(nMax,(LNCOL-1)), x(nMax) )
    OPEN(UNIT = 1,FILE = FNAME)
    PRINT *, 'FILE IS OPEN'
    i=0
    oldSize=0
    newSize=0
    DO
       READ(1,FMT=*,IOSTAT=ioerror) intmpx, intmpf(1:(LNCOL-1))
       ! WRITE(6,*) intmpx, intmpf(1:6)
       !
       ! PRINT *,'ioerror ', ioerror
       ! 
       IF (ioerror==0) THEN ! all the arrays should increase.
          ! First copy the old arrays into temporary arrays
          oldSize = i
          i=i+1
          newSize = i
          
          ALLOCATE(array_tmpX(oldSize))
          ALLOCATE(array_tmp(oldSize,(LNCOL-1)))

          array_tmpX(1:oldSize)=x(1:oldSize)
          DEALLOCATE(x)
          
          ! Now increase size of x to store new value
          ! WRITE(*,*) 'i = ', i
          ALLOCATE(x(newSize))

          ! Copy back temporary array into x
          x(1:oldSize)=array_tmpX(1:oldSize)

          ! And place new varable in ith location
          x(i)=intmpx

          ! Do the same with the actual function
          DEALLOCATE(array_tmp)
          ALLOCATE(array_tmp(oldSize,(LNCOL-1)))

          array_tmp(1:oldSize,1:(LNCOL-1)) = arr(1:oldSize,1:(LNCOL-1))
          DEALLOCATE(arr)
          ALLOCATE(arr(newSize,(LNCOL-1)))

          arr(1:oldSize,1:(LNCOL-1)) = array_tmp(1:oldSize,1:(LNCOL-1))
          arr(i,1:(LNCOL-1))=intmpf(1:(LNCOL-1))
         !  WRITE(6,*) 'NEW', x(i), arr(i,1:6)
          DEALLOCATE(array_tmp)
          DEALLOCATE(array_tmpX)
       ELSE IF (ioerror == -1) THEN
          IF (i == 0) THEN
             STOP 'Error:  Input file is possibly empty'
          END IF
          nMax = i
          EXIT
       ELSE
          STOP 'ERROR'
       END IF
    ENDDO
  END SUBROUTINE READFILE

!========================================================
  SUBROUTINE WRITEFILE(OUTFILE,x,res,nMax,WNCOL)
!========================================================
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: OUTFILE
    INTEGER, INTENT(IN) :: nMax
    INTEGER, INTENT(IN) :: WNCOL
    !REAL(DP), INTENT(IN) :: x(nmax), res(nMax,6)
    REAL(DP), INTENT(IN) :: x(nmax), res(nMax,50)
    INTEGER :: i
    OPEN(UNIT=2,FILE=OUTFILE)
    DO i = 1, nMax, 1
       !WRITE(2,'(F10.5,6E16.6)') x(i), res(i,1:(WNCOL-1))
       WRITE(2,'(F10.5,50E16.6)') x(i), res(i,1:(WNCOL-1))
    ENDDO
    CLOSE(2)    
  END SUBROUTINE WRITEFILE
!!!!========================================================
!!!  SUBROUTINE INITIALIZE_G(i,x,nMax,smear)
!!!!========================================================
!!!    IMPLICIT NONE
!!!    INTEGER, INTENT(IN) :: i, nMax
!!!    REAL(DP), INTENT(IN) :: x(nMax)
!!!    REAL(DP), INTENT(OUT) :: smear(nMax)
!!!    REAL(DP) :: fwhm = 0.030 ! eV
!!!    INTEGER :: j
!!!    
!!!    DO j=1,nMax
!!!       smear(j) = DEXP(-log(2.d0)*(x(j)-x(i))**2/fwhm**2)
!!!    END DO
!!!    smear(:) = SQRT(log(2.d0)/(pi*fwhm**2))*smear(:)
!!!    !    IF (i==501) THEN
!!!    !       DO j=1,nMax
!!!    !          WRITE (99,*) x(j), smear(j)
!!!    !       END DO
!!!    !    END IF
!!!  END SUBROUTINE INITIALIZE_G
!!!========================================
  SUBROUTINE INITIALIZE_G2(fwhm,x,nMax,smear,tMax)
!!!========================================
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nMax, tMax
    REAL(DP), INTENT(IN) :: x(nMax)
    REAL(DP), INTENT(OUT) :: smear(-tMax:tMax)
    REAL(DP), INTENT(IN) :: fwhm
    INTEGER :: j
    
    ! j =-N  =>  x(N+1)
    ! j =-1  =>  x(2)=h
    ! j = 0  =>  x(1)=0
    ! j = 1  =>  x(2)=h
    ! j = N  =>  x(N+1)
    
    DO j=-tMax,-1
       smear(j) = DEXP(-4.d0*log(2.d0)*(x(-j+1))**2/fwhm**2)
    END DO
    DO j=0,tMax
       smear(j) = DEXP(-4.d0*log(2.d0)*(x(j+1))**2/fwhm**2)
    END DO
    smear(-tMax:tMax) = 2.d0/fwhm*SQRT(log(2.d0)/pi)*smear(-tMax:tMax)
  END SUBROUTINE INITIALIZE_G2
!!!=================================================
  SUBROUTINE GAUSSIAN_SMEAR(nMax,x,arr,res,h,mhwf)
!!!=================================================
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nMax
    REAL(DP), INTENT(IN) :: arr(nMax), x(nMax), h
    REAL(DP), INTENT(INOUT) :: res(nMax)
    REAL(DP), INTENT(IN) :: mhwf 
!!!    REAL(DP), ALLOCATABLE, INTENT(INOUT) :: smear(:)
    REAL(DP), ALLOCATABLE :: smear(:)
    REAL(DP), ALLOCATABLE :: partToSmear(:)
!!! for simplicity we define a new array newArr which will extend
!!! beyond the limits of arr, so that the integration can be done
!!! easily for the endpoints of the resulting function.
    REAL(DP), ALLOCATABLE :: newArr(:)
    INTEGER :: tMax
    REAL(DP) :: integrand1(nMax), integrand2(nMax)
    REAL(DP) :: energy
    INTEGER :: i, j, k
    REAL(DP) :: fwhm = 0.15! eV
    
   !! write(*,*)"mhwf=",mhwf
     fwhm=mhwf
    write(*,*)"FWHM in GAUSSIAN_SMEAR subroutine=",fwhm  
    !!!!!STOP
    tMax = NINT(2.d0*fwhm/h)
    WRITE(6,*)'tMax = ', tMax
    ALLOCATE ( partToSmear(-tMax:tMax) )
    ALLOCATE ( smear(-tMax:tMax) )
    ALLOCATE ( newArr(1-tMax:nMax+tMax) )
!!! Since the function is even, and we assume that the input function
!!! extends to x(1)=0.0d0.
    k=1
    DO j = 0, 1-tMax, -1
       k = k + 1
       newArr( j ) = arr( k )
    END DO
    newArr( 1 : nMax ) = arr(1:nMax)
    newArr( nMax + 1 : nMax + tMax) = 0.d0
    
!!! write out function to look at for debugging
    DO j= -tMax + 1, nMax
       IF (j .LE. 0) write(10,*) -x(-j+2), newArr(j)
       IF (j .GE. 1) write(10,*) x( j  ), newArr(j)
    END DO
    
    CALL INITIALIZE_G2(fwhm,x,nMax,smear,tMax)
    
    DO i= 1, nMax    ! this upper bound can be safely changed I think
!!! Find part to smear
       partToSmear(-tMax : tMax) = newArr( i - tMax : i + tMax)
       
       ! For Simpson's Rule
       ! partToSmear(-tMax) = partToSmear(-tMax)
       ! partToSmear( tMax) = partToSmear( tMax)
       ! partToSmear(-tMax+1:tMax-1:2) = 4.d0*partToSmear(-tMax+1:tMax-1:2)
       ! partToSmear(-tMax+2:tMax-2:2) = 2.d0*partToSmear(-tMax+2:tMax-2:2)
       
       res(i) = SUM(partToSmear(:)*smear(:))
       res(i) = res(i)*h
       
    END DO
    DEALLOCATE ( partToSmear )
    DEALLOCATE ( smear )
    DEALLOCATE ( newArr )
  END SUBROUTINE GAUSSIAN_SMEAR
!!!<><><><><><><><><><><><><><><><><><><><><><>
!!!<><><><><><><><><><><><><><><><><><><><><><>
!!!<><><><><><><><><><><><><><><><><><><><><><>
!===========================
 subroutine NumberofRows(filein, NofR)
!! Cabellos 28 Septiembre 2008 a 16:20
implicit none
character(len=255), intent(IN) :: filein
integer, intent(out) :: NofR
!!----local variables
integer :: n, error,i,j
character(len=255) :: s
real :: value
!!!===begin to count the number of columns of file========
 ! write(*,*)"Openeing ESTE",trim(filein)
   open(99,FILE=filein,status="unknown")
    do
      read (99,'(a)',end=100)s
       !write(*,*)"s",s
        do i =1,254 ! The very maximum that the string can contain
           read( s, *, iostat=error ) ( value,j=1,i )
            if ( error .ne. 0 ) then
             NofR =i-1
             !write(*,*)NofR
              exit
            endif
        enddo
    end do
100 continue
    !write(*,*)NofR
   close(99)
end subroutine NumberofRows
!!!!<><><><><><><><><><><><><><><><><><><><><><
!!!!<><><><><><><><><><><><><><><><><><><><><><
!!!<><><><><><><><><><><><><><><><><><><><><><>
!!!<><><><><><><><><><><><><><><><><><><><><><>
!!!<><><><><><><><><><><><><><><><><><><><><><>
END MODULE SMEAR_MOD

!!!========================================================
PROGRAM SMEARING
!!!=========================================================
  ! A simple smearing implementation
  USE SMEAR_MOD
  IMPLICIT NONE
  INTEGER :: numberOfCommandLineArguments
  integer   :: IArgC
  INTEGER :: nMax=0
  REAL(DP), ALLOCATABLE :: res(:,:)
  REAL(DP) :: tot, h
  CHARACTER(1) :: flag
  !!CHARACTER(60) :: fname, ofname
  CHARACTER(255) :: fname, ofname
  CHARACTER(255) :: mhw
  REAL(DP) :: mhwf
  INTEGER :: itmp
  INTEGER :: NCOL
  numberOfCommandLineArguments = iargc()
  if (numberOfCommandLineArguments.eq.0) THEN
  WRITE (6,*) 'Smearing v0.1  -- NO GARAUNTEES!'
  WRITE (6,*) 'smear [1] [inputfile2smear] [outputfile] [fwhm]'
  WRITE (6,*) 'stopping right now ...!'
  STOP
  end if
  
  write(*,*)"numberOfCommandLineArguments ",numberOfCommandLineArguments
  
  CALL GETARG(1,flag)
  IF (flag == '1') THEN
     PRINT *,'Gaussian Smearing'
  ELSE IF(flag == '2')THEN
     PRINT *,'Lorentzian Smearing'
  ELSE
     STOP 'Wrong option for smearing type flag'
  END IF
  CALL GETARG(2,fname)
  PRINT *,'Opening :  ',trim(fname)
  CALL GETARG(3,ofname)
  PRINT *,'Output  :  ',trim(ofname)
  call NumberofRows(fname, NCOL)    
  write(*,*)"NUMERO DE COLUMNAS ",NCOL
  !!!! CALL GETARG(4,fwhm)
 if (numberOfCommandLineArguments.eq.4) THEN
    CALL GETARG(4,mhw)
    !read(mhw,"( f5.4 )" )mhwf
    read(mhw,*)mhwf
    write(6,*)"Taking the input FWHM=",mhwf
    write(6,*)"This value is written to FILE: fromSmear(for use scripts)"    
 else 
    mhwf=0.15
    write(6,*)"Taking the DEFAULT FWHM=",mhwf  
    write(6,*)"This value is written to FILE: fromSmear(for use scripts)"    
 end if
  !! Just put this value in a file called fromSmear
   call system('rm -f fromSmear')
   OPEN(UNIT=56,FILE = 'fromSmear')
   ! write(56,"(1E12.5)")mhwf
    write(56,"(2F6.3)")mhwf
   close(56)

  !!write(*,*)mhwf 
  !!stop 
  CALL READFILE(fname,nMax,NCOL)
  
  ! x and arr should be read by now. nMax should be the array length
  
  !!  ALLOCATE(smear(nMax),res(nMax))  ! the result goes in res
  !ALLOCATE(res(nMax,6))  ! the result goes in res
  ALLOCATE(res(nMax,NCOL))  ! the result goes in res
  
  ! Output some info from File
  WRITE(6,*) 'From the file, I learned that:'
  WRITE(6,*) 'o nMax is ', nMax
  WRITE(6,*) 'o xmin is ', x(1)
  WRITE(6,*) 'o xmax is ', x(nMax)
  tot = SUM(x)
  h = tot * 2.d0/(nMax)/(nMax-1)
  WRITE(6,FMT='(A28,E18.10)') 'o spacing is approximately ',h
  
  !  CALL INITIALIZE(flag,nMax,x,smear)
!!!  IF (flag == '1') CALL GAUSSIAN_SMEAR(nMax,x,arr,smear,res,h)
  
  
  IF (flag == '1') THEN
     DO itmp = 1, (NCOL-1)
        WRITE(*,*) itmp
        CALL GAUSSIAN_SMEAR(nMax,x,arr(1:nMax,itmp),res(1:nMax,itmp),h,mhwf)
     END DO
  END IF

  !  IF (flag == '2') CALL LORENTZIAN_SMEAR(nMax,x,arr,smear,res)
  
  CALL WRITEFILE(ofname,x,res,nMax,NCOL)
  
END PROGRAM SMEARING

!**************************************************
