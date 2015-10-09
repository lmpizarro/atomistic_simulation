module constants

     use types, only :dp
     real(dp), parameter ::                &
        PI           = 3.1415926535898_8, & ! pi
        N_AVOGADRO   = 0.602214129_8,     & ! Avogadro's number in 10^24/mol
        K_BOLTZMANN  = 8.6173324e-11_8,   &
        kb           = 1.0

end module constants
