!*******************************( SORT )********************************
! This subroutine sorts the corresponding arrays "LOCATE" and "VALUE"
! so that the elements of "VALUE" are in order highest to lowest.
! It uses the QUICKSORT procedure by C.A.R. HOARE  and described in many
! programming textbooks.
!***********************************************************************
      SUBROUTINE SORT (NCR, LOCATE, VALUE)

      DIMENSION LOCATE(4000), VALUE(4000), MARRAY(4000), NARRAY(4000)

      M = 1
      N = NCR
      KOUNT = 0

 1111 IF (M .LT. NCR) THEN
         IF (M .LT. N) THEN
! SORT PART OF ARRAY BOUNDED BY M AND N
            I = M - 1
		    J = N
			REF = VALUE (N)
 1112       CONTINUE
               I = I + 1
            IF (VALUE(I).GT.REF) GO TO 1112
 1113       CONTINUE
            J = J - 1
            IF (J .GT.1 .AND. VALUE(J) .LT. REF) GO TO 1113
			IF (I .LT. J) THEN
			   D11 = VALUE(I)
			   VALUE(I) = VALUE(J)
			   VALUE(J) = D11
			   M11 = LOCATE(I)
			   LOCATE(I) = LOCATE(J)
			   LOCATE(J) = M11
			ELSE
			   GO TO 1114
			END IF
			GO TO 1112
 1114       D12 = VALUE(N)
			VALUE(N) = VALUE(I)
			VALUE(I) = D12
			M11 = LOCATE(N)
			LOCATE(N) = LOCATE(I)
			LOCATE(I) = M11
			
! MANIPULATE STACK OF VALUES OF M AND N

		    KOUNT = KOUNT + 1
			MARRAY(KOUNT) = I + 1
			NARRAY(KOUNT) = N
			N = I - 1
		 ELSE
		    M = MARRAY(KOUNT)
			N = NARRAY(KOUNT)
			KOUNT = KOUNT - 1
		 END IF
		 GO TO 1111
	  END IF

! FURTHER SORT THE CELLS SO THAT WITHIN GROUPS HAVING THE SAME
! VALUE, THOSE NEARER THE PLANT STEM HAVE PRIORITY IN THE LIST

      i = 1
      DO WHILE (i .LT. NCR)
		 IF (VALUE(i).GT.VALUE(i+1)) THEN
		    i = i + 1
		 ELSE
		    IF (LOCATE(i) .GT. LOCATE(i+1)) THEN
			   M12 = LOCATE(i)
			   LOCATE(i) = LOCATE(i+1)
			   LOCATE(i+1) = M12
			   i = i - 1
			   IF (i .LT. 1) i = 1
			ELSE
			   i = i + 1
			END IF
		 END IF
	  END DO

      RETURN
      END
