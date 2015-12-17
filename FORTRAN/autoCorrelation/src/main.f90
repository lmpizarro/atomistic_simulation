PROGRAM autoCorrelation

    IMPLICIT NONE

    CHARACTER(LEN("acf.out")) :: outputFile = "acf.out"
    INTEGER :: Nat,i,nbTimeStepsInTraj,iostat,dt,t,dmax,d
    DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: r ! position of site i at timestep t
    DOUBLE PRECISION :: time1, time0
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: acf
    CHARACTER(LEN=300) :: arg,trajectoryFileName,algo

    CALL CPU_TIME(time0)

    ! read all arguments that MUST be given at execution
    CALL readArguments(Nat,trajectoryFileName,algo,dmax)
    i =NbOfLinesInTraj()     ! deduce the number of timesteps in the trajectory from the number of lines in the trajectory file
    CALL testConsistencyOfAtomNumber(i,Nat)
    nbTimeStepsInTraj = i/Nat
    PRINT*,'I found',nbTimeStepsInTraj,' timesteps in your trajectory called ',trim(adjustl(trajectoryFileName))
    PRINT*,'Please be patient... Everything looks fine...'

    ! read vector of all sites i at all timesteps t
    ALLOCATE( r(nbTimeStepsInTraj,dmax,Nat), SOURCE=0.d0 )
    CALL opentraj
    DO t=1,nbTimeStepsInTraj
        IF( mod(t,1000)==1 .AND. t/=1) PRINT*,"For your information, I've read ",t-1," timesteps.",nbTimeStepsInTraj-(t-1)&
          ,"more to come"
        DO i=1,Nat
            READ(10,*) r(t,1:dmax,i)
        END DO
    END DO
    CALL closetraj
    PRINT*,"All trajectories are read. Now, be patient..."

    ! compute autocorrelation function acf(dt)= <v_i(t).v_i(t+dt)>_{i,t}   where . is the scalar product
    ALLOCATE( acf(0:nbTimeStepsInTraj-1), SOURCE=0.d0 )

    IF (TRIM(ADJUSTL(algo))=="bruteforce") THEN
        CALL bruteforce(r,acf)
    ELSE IF (TRIM(ADJUSTL(algo))=="fourierspace") THEN
        CALL fourierspace(r,acf)
    ELSE
        STOP "I do not understand the algorithm you ask for. I only recognize bruteforce and fourierspace. STOP"
    END IF

    ! acf(t) will be written in file unit 11
    OPEN(11,FILE=outputfile)
        DO dt = 0, nbTimeStepsInTraj-1
            WRITE(11,*) dt, acf(dt)
        END DO
    CLOSE(11)

    CALL CPU_TIME(time1)
    IF ((time1-time0)>120.d0) THEN
        PRINT*,"-- Finished in ",NINT((time1-time0)/60.d0)," min. GGHF ;) --"
    ELSE
        PRINT*,"-- Finished in ",NINT(time1-time0)," s. GGHF ;) --"
    END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CONTAINS


        SUBROUTINE testConsistencyOfAtomNumber(i,Nat)
            INTEGER, INTENT(IN) :: i,Nat
            IF( MODULO(i,Nat)/=0 ) THEN
                PRINT*,"ERROR: Inconsistency between the number of lines in file ",TRIM(ADJUSTL(trajectoryFileName))
                PRINT*,"=====  and the number of atoms given as argument:",Nat
                PRINT*,"       Modulo(nlines,Nat)/=0"
                STOP
            END IF
        END SUBROUTINE testConsistencyOfAtomNumber

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE opentraj
            CALL inquireFileExistence(trajectoryFileName)
            ! read positions
            OPEN(10, FILE=trajectoryFileName,STATUS='old',IOSTAT=iostat)
            IF (iostat /= 0) THEN
                WRITE(*,*) 'ERROR: l.82 File open failed for',trajectoryFileName
                WRITE(*,*) '=====  Error code is ', iostat
                STOP
            END IF
        END SUBROUTINE opentraj

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE closetraj
            CLOSE(10)
        END SUBROUTINE closetraj

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE inquireFileExistence(fileName)
            CHARACTER(LEN=*), INTENT(IN) :: fileName
            INTEGER, PARAMETER :: stderr = 0
            LOGICAL :: exist
            INQUIRE(FILE=fileName, EXIST=exist)
            IF( .NOT. exist) THEN
                WRITE(stderr,*) "ERROR: File ",fileName," does not exist."
                write(stderr,*) "====="
                STOP
            END IF
        END SUBROUTINE inquireFileExistence

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        FUNCTION NbOfLinesInTraj()
            INTEGER :: NbOfLinesInTraj
            CALL opentraj
            ! computes the number of lines in traj.in and deduces the number of timesteps
            NbOfLinesInTraj = -1
            DO WHILE (iostat == 0)
                READ(10,*,IOSTAT=iostat)
                NbOfLinesInTraj = NbOfLinesInTraj + 1
            END DO
            CALL closetraj
        END FUNCTION NbOfLinesInTraj

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE readArguments(Nat,trajectoryFileName,algo,nColumns)
            IMPLICIT NONE
            INTEGER, INTENT(OUT) :: Nat,nColumns
            CHARACTER(LEN=*),INTENT(OUT) :: trajectoryFileName,algo
            CALL GET_COMMAND_ARGUMENT(1,arg,STATUS=i)
                IF (i<0) THEN
                    STOP "STOP. The length of the argument is too big for me :( "
                ELSE IF (i>0) THEN
                    STOP "Argument retrieval failed. You should execute the program with the number of atoms as 1st argument"
                END IF
                READ(arg,*) Nat
                PRINT*,"Number of atoms: ",Nat
            CALL GET_COMMAND_ARGUMENT(2,arg,STATUS=i)
                IF (i<0) THEN
                    STOP "STOP. The length of the argument is too big for me :( "
                ELSE IF (i>0) THEN
                    STOP "Argument retrieval failed. You should execute the program with the filename as 2nd argument"
                END IF
                READ(arg,*) i
                SELECT CASE (i)
                CASE (1)
                    algo='bruteforce'
                CASE (2)
                    algo="fourierspace"
                CASE DEFAULT
                    STOP "Argument for the algo must be 1 (bruteforce) or 2 (fourierspace)"
                END SELECT
                PRINT*,"I'll use algorithm ",TRIM(ADJUSTL(algo))
            CALL GET_COMMAND_ARGUMENT(3,arg,STATUS=i)
                IF (i<0) THEN
                    STOP "STOP. The length of the argument is too big for me :( "
                ELSE IF (i>0) THEN
                    STOP "Argument retrieval failed. You should execute the program with the number of atoms as 1st argument"
                END IF
                READ(arg,*) nColumns
                PRINT*,"Number of columns in file: ",nColumns
            CALL GET_COMMAND_ARGUMENT(4,arg,STATUS=i)
                IF (i<0) THEN
                    STOP "STOP. The length of the argument is too big for me :( "
                ELSE IF (i>0) THEN
                    STOP "Argument retrieval failed. You should execute the program with the algorithm as 3nd argument"
                END IF
                trajectoryFileName = TRIM(ADJUSTL(arg))
                PRINT*,"I'll read the trajectory from ",TRIM(ADJUSTL(trajectoryFileName))
        END SUBROUTINE readArguments

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE bruteforce(r,acf)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:,:,:), CONTIGUOUS, INTENT(IN) :: r
            DOUBLE PRECISION, DIMENSION(0:), CONTIGUOUS, INTENT(OUT) :: acf
            DOUBLE PRECISION, DIMENSION(SIZE(r,1),SIZE(r,2)) :: ri ! position of site i at timestep t for a given atom
            INTEGER :: i,dt,nt,Nat,nbTimeStepsInTraj,d,dmax
            DOUBLE PRECISION :: time0, time1, remainingTimeInSec
            PRINT*,"I'll use the brute force algorithm"
            CALL CPU_TIME(time0)
            nbTimeStepsInTraj = SIZE(r,1)
            dmax = SIZE(r,2)
            Nat = SIZE(r,3)   !r(nbTimeStepsInTraj,x:z,Nat)
            DO i= 1, Nat
                ri = r(:,:,i)
                DO d=1,dmax
                    DO dt = 0, nbTimeStepsInTraj-1
                        nt = nbTimeStepsInTraj-dt
                            acf(dt) =acf(dt)+SUM( ri(1:nt,d)*ri(1+dt:nt+dt,d) )/DBLE(nt)
                    END DO
                END DO

                CALL CPU_TIME(time1)
                IF(MODULO(i,MAX(INT(Nat*0.1),1))==0) THEN
                    remainingTimeInSec = DBLE(Nat-i)*(time1-time0)/DBLE(i)
                    IF (remainingTimeInSec>60.d0) THEN
                        PRINT*,'Estimated remaining time ≈ ',NINT(remainingTimeInSec/60.d0),' min'
                    ELSE
                        PRINT*,'Estimated remaining time ≈ ',NINT(remainingTimeInSec),' s'
                    END IF
                END IF
            END DO
            acf = acf/DBLE(Nat)
        END SUBROUTINE bruteforce

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        SUBROUTINE fourierspace(r,acf)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(:,:,:), CONTIGUOUS, INTENT(IN) :: r
            DOUBLE PRECISION, DIMENSION(0:), CONTIGUOUS, INTENT(OUT) :: acf
            INTEGER(KIND=KIND(1)), PARAMETER :: i2b = KIND(1) !> simple precision integer
            INTEGER(i2b), PARAMETER :: dp = KIND(0.0d0) !> double precision real
            INTEGER(i2b), PARAMETER :: sp = KIND(0.0) !> simple precision real
            INTEGER(i2b), PARAMETER :: i4b = 2_i2b * i2b !> double precision integer
            TYPE :: fftw3Needs
                INTEGER(i4b) :: plan_forward, plan_backward ! descriptors of our FFTs
                REAL(dp),    ALLOCATABLE, DIMENSION(:) :: in_forward, out_backward ! input of FFT and output (result) of inverse FFT
                COMPLEX(dp), ALLOCATABLE, DIMENSION(:) :: out_forward, in_backward ! output (result) of FFT and input of inverse FFT
            END TYPE
            TYPE( fftw3Needs ) :: fftw3
            INTEGER :: s,nt
            INCLUDE "/usr/include/fftw3.f" ! SHOULD BE REMOVED VERY FAST BUT NEEDS A MAKEFILE OR ./CONFIG TO INCLUDE THE CORRECT TREE IN MAC, UBUNTU, FEDORA, etc.

            PRINT*,"I'll use the fourier space algorithm to compute the autocrosscorrelation"
            CALL CPU_TIME(time0)

            dmax = SIZE(r,2)
            nt = SIZE(r,1)
            s = nt*2
            Nat = SIZE(r,3)

            ALLOCATE(fftw3%in_forward(s), SOURCE=0.d0)
            ALLOCATE(fftw3%out_forward(s/2+1)) ! when from real to complex then half the signal is enough
            ALLOCATE(fftw3%out_backward(s))
            ALLOCATE(fftw3%in_backward(s/2+1))

            CALL dfftw_plan_dft_r2c_1d ( fftw3%plan_forward,  s, fftw3%in_forward,  fftw3%out_forward,  FFTW_ESTIMATE)!FFTW_MEASURE )
            CALL dfftw_plan_dft_c2r_1d ( fftw3%plan_backward, s, fftw3%in_backward, fftw3%out_backward, FFTW_ESTIMATE)!FFTW_MEASURE )

            acf=0.d0
            DO s=1,Nat ! atoms
                DO i=1,dmax ! x,y,z ...
                    fftw3%in_forward(LBOUND(r,1):UBOUND(r,1)) =r(:,i,s)
                    CALL dfftw_execute ( fftw3%plan_forward )
                    fftw3%in_backward = fftw3%out_forward* CONJG(fftw3%out_forward) ! I AM SURE AN INTRINSIC FUNCTION EXISTS FOR x=z.z*
                    CALL dfftw_execute ( fftw3%plan_backward )
                    acf =acf +fftw3%out_backward(LBOUND(r,1):UBOUND(r,1))
                END DO
            END DO

            acf =acf/(DBLE(Nat*SIZE(fftw3%in_forward,1)))
            DO s=1,nt
                acf(s-1)=acf(s-1)/DBLE(nt-(s-1))
            END DO

        END SUBROUTINE fourierspace


END PROGRAM autoCorrelation
