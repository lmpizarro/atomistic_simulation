module constants

     use types, only :dp
     real(dp), parameter ::                &
        PI           = 3.1415926535898_dp, & ! pi
        N_AVOGADRO   = 0.602214129_dp,     & ! Avogadro's number in 10^24/mol
        K_BOLTZMANN  = 8.6173324e-5_dp,   & ! ev K-1
        kb           = 1.0_dp,            &
        K_B_KJ       = 1.3806488e-23_dp      ! J K-1

end module constants
