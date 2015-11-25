!********************************************************************************
!** FICHE F.23.  ROUTINE TO SET UP ALPHA FCC LATTICE OF LINEAR MOLECULES       **
!** This FORTRAN code is intended to illustrate points made in the text.       **
!** To our knowledge it works correctly.  However it is the responsibility of  **
!** the user to test it, if it is to be used in a research application.        **
!********************************************************************************
! http://physics.ujep.cz/~mlisal/mol_simul/4_tyden_results.html
! https://github.com/itamblyn/FORTRAN/blob/master/md/mdx/crystal.f90
! https://github.com/libAtoms/QUIP
!http://users-phys.au.dk/santi/codes/doc/modules/crystal/crystal_mod_f90.html


        SUBROUTINE FCC

        COMMON / BLOCK1 / RX, RY, RZ, EX, EY, EZ

!    *******************************************************************
!    ** SETS UP THE ALPHA FCC LATTICE FOR N LINEAR MOLECULES.         **
!    **                                                               **
!    ** THE SIMULATION BOX IS A UNIT CUBE CENTRED AT THE ORIGIN.      **
!    ** N SHOULD BE AN INTEGER OF THE FORM ( 4 * ( NC ** 3 ) ),       **
!    ** WHERE NC IS THE NUMBER OF FCC UNIT CELLS IN EACH DIRECTION.   **
!    ** SEE FIGURE 5.10 FOR A DIAGRAM OF THE LATTICE AND A            **
!    ** DEFINITION OF THE FOUR ORIENTATIONAL SUBLATTICES.             **
!    **                                                               **
!    ** PRINCIPAL VARIABLES:                                          **
!    **                                                               **
!    ** INTEGER N                    NUMBER OF MOLECULES              **
!    ** REAL    RX(N),RY(N),RZ(N)    MOLECULAR POSITIONS              **
!    ** REAL    EX(N),EY(N),EZ(N)    UNIT VECTORS GIVING ORIENTATIONS **
!    ** REAL    RROOT3               1.0 / SQRT ( 3.0 )               **
!    *******************************************************************

        INTEGER     N, NC
        REAL        RROOT3

        PARAMETER ( NC = 3, N = 4 * NC ** 3 )
        PARAMETER ( RROOT3 = 0.5773503 )

        REAL        RX(N), RY(N), RZ(N), EX(N), EY(N), EZ(N)
        REAL        CELL, CELL2
        INTEGER     I, IX, IY, IZ, IREF, M

!    *******************************************************************

!    ** CALCULATE THE SIDE OF THE UNIT CELL **

        CELL  = 1.0 / REAL ( NC )
        CELL2 = 0.5 * CELL

!    ** BUILD THE UNIT CELL **

!    ** SUBLATTICE A **

        RX(1) =  0.0
        RY(1) =  0.0
        RZ(1) =  0.0
        EX(1) =  RROOT3
        EY(1) =  RROOT3
        EZ(1) =  RROOT3

!    ** SUBLATTICE B **

        RX(2) =  CELL2
        RY(2) =  CELL2
        RZ(2) =  0.0
        EX(2) =  RROOT3
        EY(2) = -RROOT3
        EZ(2) = -RROOT3

!    ** SUBLATTICE C **

        RX(3) =  0.0
        RY(3) =  CELL2
        RZ(3) =  CELL2
        EX(3) = -RROOT3
        EY(3) =  RROOT3
        EZ(3) = -RROOT3

!    ** SUBLATTICE D **

        RX(4) =  CELL2
        RY(4) =  0.0
        RZ(4) =  CELL2
        EX(4) = -RROOT3
        EY(4) = -RROOT3
        EZ(4) =  RROOT3

!    ** CONSTRUCT THE LATTICE FROM THE UNIT CELL **

        M = 0

        DO 99 IZ = 1, NC

           DO 98 IY = 1, NC

              DO 97 IX = 1, NC

                 DO 96 IREF = 1, 4

                    RX(IREF+M) = RX(IREF) + CELL * REAL ( IX - 1 )
                    RY(IREF+M) = RY(IREF) + CELL * REAL ( IY - 1 )
                    RZ(IREF+M) = RZ(IREF) + CELL * REAL ( IZ - 1 )

                    EX(IREF+M) = EX(IREF)
                    EY(IREF+M) = EY(IREF)
                    EZ(IREF+M) = EZ(IREF)

96               CONTINUE

                 M = M + 4

97            CONTINUE

98         CONTINUE

99      CONTINUE

!    ** SHIFT CENTRE OF BOX TO THE ORIGIN **

        DO 100 I = 1, N

           RX(I) = RX(I) - 0.5
           RY(I) = RY(I) - 0.5
           RZ(I) = RZ(I) - 0.5

100     CONTINUE

        RETURN
        END




