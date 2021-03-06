! ############### CONST ###################
MODULE CONST ! math constants
	COMPLEX, PARAMETER :: Z0 = (0.,0.), Z1 = (1.,0.), ZI = (0.,1.)
	REAL, PARAMETER :: PI = 4*ATAN(1.)
END MODULE CONST
! ############### MODEL ###################
MODULE MODEL ! model settings
	USE CONST
	! physical parameters
	REAL    :: THETA = 1.*PI ! theta-term
	REAL    :: CROSS = 1.    ! crossing: 1. = allow, 0. = avoid
	REAL    :: BETA = 0.53    ! 0.440687
	! lattice parameter
	INTEGER :: L1 = 10
	INTEGER :: L2 = 20
	! ATRG parameters
	INTEGER :: DCUT = 128
END MODULE MODEL
! ############## PHYSICS ###################
MODULE PHYSICS
	USE TENSORIAL
	TYPE HUGE_TENSOR
		REAL :: LEV
		TYPE(TENSOR) :: TEN
	END TYPE HUGE_TENSOR
CONTAINS
! ------------ set MPO tensor ---------------
! square lattice MPO (new)
FUNCTION GET_MPO() RESULT(T)
! output: T - MPO tensor
!     │
!     3
! ─ 1 T 2 ─
!     4
!     │ 
	USE MODEL
	TYPE(TENSOR) :: T ! MPO tensor output
	! local variables
	TYPE(TENSOR) :: X, Y, UA, UB, S
	COMPLEX, ALLOCATABLE :: WA(:,:)
	COMPLEX :: Q, B
	
! ++++++++ set the vertex tensor here ++++++++
	Q = THETA * ZI
	B = BETA * Z1
	X =  TENSOR([2,2,2,2],[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[EXP(2*B),EXP(Q/4.),EXP(Q/4.),EXP(Z0),EXP(Q/4.),CROSS*EXP(-2*B - Q/2.),EXP(Z0),EXP(-Q/4.),EXP(Q/4.),EXP(Z0),CROSS*EXP(-2*B + Q/2.),EXP(-Q/4.),EXP(Z0),EXP(-Q/4.),EXP(-Q/4.),EXP(2*B)])
! ++++++++++++++++++++++++++++++++++++++++++++
	! SVD of X to UA - S - UB
	CALL SVD(X,[1,4],[3,2],UA,UB,S)
	S%VALS = SQRT(S%VALS) ! split S to half
	! attach the half S to U
	UA = TEN_PROD(S,UA,[2],[3])
	UB = TEN_PROD(S,UB,[2],[3])
	! set Y tensor
	Y = EYE_TEN([2,2,2])
	! contract U from both sides with Y
	UA = TEN_PROD(UA,Y,[3],[3])
	UB = TEN_PROD(UB,Y,[3],[3])
	T = TEN_TRANS(TEN_PROD(UB,UA,[2,4],[4,2]),[1,3,2,4])
END FUNCTION GET_MPO
! square lattice MPO
FUNCTION GET_MPO0() RESULT(T)
! output: T - MPO tensor
!     │
!     3
! ─ 1 T 2 ─
!     4
!     │ 
	USE MODEL
	TYPE(TENSOR) :: T ! MPO tensor output
	! local variables
	TYPE(TENSOR) :: X, Y, U, S
	COMPLEX, ALLOCATABLE :: WA(:,:)
	COMPLEX :: Q, B
	
! ++++++++ set the vertex tensor here ++++++++
	Q = THETA * ZI
	B = BETA * Z1
	X =  TENSOR([2,2,2,2],[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15],[EXP(2*B),EXP(Q/4.),EXP(Q/4.),EXP(Z0),EXP(Q/4.),CROSS*EXP(-2*B - Q/2.),EXP(Z0),EXP(-Q/4.),EXP(Q/4.),EXP(Z0),CROSS*EXP(-2*B + Q/2.),EXP(-Q/4.),EXP(Z0),EXP(-Q/4.),EXP(-Q/4.),EXP(2*B)])
! ++++++++++++++++++++++++++++++++++++++++++++
	! SVD of X to U14 U32 S
	CALL SYSVD(X,[1,4],[3,2],U,S)
	S%VALS = SQRT(S%VALS) ! split S to half
	U = TEN_PROD(S,U,[2],[3]) ! attach the half S to U
	! set Y tensor
	Y = EYE_TEN([2,2,2])
	! contract U from both sides with Y
	U = TEN_PROD(U,Y,[3],[3])
	T = TEN_TRANS(TEN_PROD(U,U,[2,4],[4,2]),[1,3,2,4])
END FUNCTION GET_MPO0
! ---------- grow column tensor -------------
! grow column tensor
FUNCTION GROW(T, L) RESULT (M)
! input: T - MPO tensor, L - layers to grow
! output: M - L layers of MPO tensors
! staking 3 M 4 =  ─ 3 T 4 ─ 3 T 4 ─ ...
	TYPE(TENSOR), INTENT(IN) :: T
	INTEGER :: L
	TYPE(HUGE_TENSOR) :: M
	! local variables
	INTEGER :: I
	
	M%LEV = 0.
	M%TEN = T ! put the first layer down
	CALL RESCALE(M,[1,3],[2,4])
	DO I = 2, L
		!       3
		!       │
		!       3
		!   ╭ 1 T 2 ╮
		!   │   4   │
		! 1 ┤   │   ├ 2
		!   │   3   │
		!   ╰ 1 M 2 ╯
		!       4
		!       │
		!       4
		M%TEN = REDUCE(TEN_FLATTEN(TEN_PROD(T,M%TEN,[4],[3]),[1,4,0,2,5,0,3,0,6]))
		CALL RESCALE(M,[1,3],[2,4])
		WRITE (*,'(I3,A,F5.1,A)') I,':',SVD_ERR*100,'%'
!		WRITE (*,'(F8.4)') M%LEV/LOG(2.)
	END DO
END FUNCTION GROW
! reduce tensor by SVD
FUNCTION REDUCE(T1) RESULT (T2)
! input: T1 to be reduct, output: T2 as the result
! reduce leg 1,2 by HOSVD
	USE MODEL
	TYPE(TENSOR), INTENT(IN) :: T1
	TYPE(TENSOR) :: T2
	! local tensors
	TYPE(TENSOR), ALLOCATABLE :: US(:)
	
	T2 = T1 ! copy T1 to T2
	! perform HOSVD
	CALL HOSVD(T2, US, [1,2], [DCUT,DCUT])
	!                           │
	!                           3
	! 1 ─ 2 US2 1 ━ 1 US1 2 ─ 1 T2 2 ─ 2
	!                           4
	!                           │
	T2 = TEN_PROD(TEN_PROD(US(2),US(1),[1],[1]),T2,[2],[1])
END FUNCTION REDUCE
! ----------- measurement -------------------
! correlation function
FUNCTION CORR(M, O, N) RESULT (C)
! input: M - column tensor, O - observable, N - data length (lattice 2N)
! output: C - correlation function
	TYPE(HUGE_TENSOR), INTENT(IN) :: M
	TYPE(TENSOR), INTENT(IN) :: O
	INTEGER, INTENT(IN) :: N
	COMPLEX, ALLOCATABLE :: C(:)
	! local tensor
	TYPE(HUGE_TENSOR) :: P, PO, Z
	TYPE(HUGE_TENSOR), ALLOCATABLE :: PS(:)
	INTEGER :: J
	
	REAL :: A
	COMPLEX :: B
	
	! close M 3,4 leg to form transfer matrix
	! uniform
	P%LEV = M%LEV
	P%TEN = TEN_TRACE(M%TEN,[3],[4]) ! self contraction
	! impurity
	PO%LEV = M%LEV
	PO%TEN = TEN_PROD(M%TEN, O, [3,4], [1,2]) ! contract with operator O
	! P and PO is not normalized
	! allocate for space
	ALLOCATE(PS(0:2*N), C(N))
	! accumulate PS
	PS(0)%LEV = 0.
	PS(0)%TEN = EYE_TEN(P%TEN%DIMS) ! initialize with identity mat
	CALL RESCALE(PS(0),[1],[2]) ! normalize
	! construct transfer matrices of all lengths
	DO J = 1, 2*N ! length from 1 to 2*N, 2*N is the lattice length
		PS(J)%LEV = PS(J-1)%LEV + P%LEV
		PS(J)%TEN = TEN_PROD(PS(J-1)%TEN,P%TEN,[2],[1])
		CALL RESCALE(PS(J),[1],[2]) ! normalize
!		WRITE (*,'(F8.4)') PS(J)%LEV/LOG(2.)
	END DO
	! collect correlation
	C = (0.,0.) ! initialize by clearing
	DO J = 1, N ! correlation only cal up to half latt
		Z%LEV = PO%LEV*2 + PS(J-1)%LEV + PS(2*N-J-1)%LEV
		Z%TEN = TEN_PROD(TEN_PROD(TEN_PROD(PO%TEN,PS(J-1)%TEN,[2],[1]),PO%TEN,[2],[1]),PS(2*N-J-1)%TEN,[2,1],[1,2])
		! take the value from Z, also restore the level factor
		! note that the transfer mat of whole latt is in PS(2*N) not PS(N)
		! because PS(2*N) is normalized, so its trace = 1., no need divided
		IF (SIZE(Z%TEN%VALS) == 1) C(J) = Z%TEN%VALS(1)*EXP(Z%LEV - PS(2*N)%LEV)
	END DO
END FUNCTION CORR
! partition function
FUNCTION FREE(M, N) RESULT (F)
	TYPE(HUGE_TENSOR), INTENT(IN) :: M
	INTEGER, INTENT(IN) :: N
	REAL :: F
	! local tensor
	TYPE(HUGE_TENSOR) :: A0, A
	INTEGER :: I
	
	A0%LEV = M%LEV
	A0%TEN = TEN_TRACE(M%TEN,[3],[4])
	A = A0
	DO I = 2, N
		A%LEV = A%LEV + A0%LEV
		A%TEN = TEN_PROD(A%TEN, A0%TEN, [2], [1])
		CALL RESCALE(A, [1], [2])
	END DO
	F = A%LEV + LOG(TR(A%TEN,[1],[2]))
END FUNCTION FREE
! ------------ misc --------------
! evaluate tensor trace
FUNCTION TR(T, LEGS1, LEGS2) RESULT (Z)
	TYPE(TENSOR), INTENT(IN) :: T
	INTEGER, INTENT(IN) :: LEGS1(:), LEGS2(:)
	COMPLEX :: Z
	! local tensor
	TYPE(TENSOR) :: A
	
	Z = (0.,0.) ! initialize
	A = TEN_TRACE(T, LEGS1, LEGS2) ! tensor trace
	IF (SIZE(A%VALS) == 1) Z = A%VALS(1) ! put value
END FUNCTION TR
! rescale huge tensor
SUBROUTINE RESCALE(M, LEGS1, LEGS2)
	TYPE(HUGE_TENSOR), INTENT(INOUT) :: M
	INTEGER, INTENT(IN) :: LEGS1(:), LEGS2(:)
	! local variables
	REAL :: R
	
	R = ABS(TR(M%TEN, LEGS1, LEGS2)) ! get ratio
	IF (R >= 1.E-12) THEN
		M%TEN%VALS = M%TEN%VALS/R ! rescale
		M%LEV = M%LEV + LOG(R) ! update level
	END IF
END SUBROUTINE RESCALE
! end of module PHYSICS
END MODULE PHYSICS
! ################ TASK ####################
MODULE TASK
	USE PHYSICS
	IMPLICIT NONE
CONTAINS
! ------------ Data --------------
SUBROUTINE COLLECT_CORR()
	USE MODEL
	USE MATHIO
	TYPE(TENSOR) :: T, O
	TYPE(HUGE_TENSOR) :: M
	COMPLEX, ALLOCATABLE :: C(:)
	CHARACTER(20) :: FILENAME
	
	T = GET_MPO()
	M = GROW(T, L1)
	O = PAULI_MAT([3])
	C = CORR(M, O, L2)
	! construct file name
	WRITE (FILENAME,'("D",I3.3,"L",I2.2,"X",I2.2,"Q",I1,"C",I1,"B",SP,I5.4)') DCUT,L1,L2,NINT(THETA/PI),NINT(CROSS),NINT(BETA*10000)
	CALL EXPORT('./data center/'//FILENAME,C)
END SUBROUTINE COLLECT_CORR
! ------------ Tests -------------
! test routine
SUBROUTINE TEST()
	TYPE(TENSOR) :: T
	TYPE(HUGE_TENSOR) :: M
	REAL :: F
	
	T = GET_MPO()
!	T = TEN_FLATTEN(TEN_PROD(T,T,[3,4],[4,3]),[1,3,0,2,4])
!	T = TEN_PROD(T,T,[2],[1])
!	T = TEN_PROD(T,T,[1,2],[2,1])
!	CALL TEN_PRINT(T)
	M = GROW(T, 4)
	F = FREE(M,8)
	PRINT *, EXP(F)
END SUBROUTINE TEST
! test correlation
SUBROUTINE TEST_CORR()
	USE MATHIO
	TYPE(TENSOR) :: T, O
	TYPE(HUGE_TENSOR) :: M
	COMPLEX, ALLOCATABLE :: C(:)
	
	T = GET_MPO()
	M = GROW(T, 13)
!	PRINT *, M%LEV, TR(M%TEN,[1,3],[2,4])
	O = PAULI_MAT([3])
	C = CORR(M, O, 20)
	CALL EXPORT('C',C)
END SUBROUTINE TEST_CORR
! end of module TASK
END MODULE TASK
! ############### PROGRAM ##################
PROGRAM MAIN
	USE TASK
	PRINT *, ' ------- ATRG -------- '
	
	CALL COLLECT_CORR()
!	CALL TEST()
!	CALL TEST_CORR()
END PROGRAM MAIN