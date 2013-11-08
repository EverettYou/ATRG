! ############### MODEL ###################
MODULE MODEL
	INTEGER, PARAMETER :: DCUT = 8
END MODULE MODEL
! ############## PHYSICS ###################
MODULE PHYSICS
	USE TENSORIAL
	IMPLICIT NONE
CONTAINS
! end of module PHYSICS
END MODULE PHYSICS
! ################ TASK ####################
MODULE TASK
	USE PHYSICS
	IMPLICIT NONE
CONTAINS
! ------------ Data --------------
! ------------ Tests -------------
! test routine
SUBROUTINE TEST()
	TYPE(TENSOR) :: T
	TYPE(TENSOR), ALLOCATABLE :: US(:)
	
	CALL TEN_LOAD('T',T)
	CALL HOSVD(T,US,[1,2])
	CALL TEN_PRINT(US(1))
	CALL TEN_PRINT(US(2))
	CALL TEN_PRINT(T)
	CALL TEN_SAVE('A',T)
END SUBROUTINE TEST
! end of module TASK
END MODULE TASK
! ############### PROGRAM ##################
PROGRAM MAIN
	USE TASK
	
	CALL TEST()
END PROGRAM MAIN