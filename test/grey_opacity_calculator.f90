module grey_opacity_calculator
   implicit none
   
   ! Physical constants
   real(8), parameter :: h_planck = 6.626070040d-27  ! erg*s
   real(8), parameter :: c_light = 2.99792458d10     ! cm/s
   real(8), parameter :: k_boltz = 1.38064852d-16    ! erg/K
   real(8), parameter :: eV_to_erg = 1.602176634d-12 ! erg/eV
   real(8), parameter :: eV_to_K = 11604.5218d0      ! K/eV
   
contains

   !================================================================
   ! Calculate Grey Opacity from 30 energy groups
   ! This is the main function you need!
   !================================================================
   subroutine calculate_grey_opacity(n_groups, energy_grid, &
                                    opacity_planck_groups, &
                                    opacity_rosseland_groups, &
                                    temperature_eV, &
                                    grey_planck, grey_rosseland)
      implicit none
      
      integer, intent(in) :: n_groups
      real(8), intent(in) :: energy_grid(n_groups+1)        ! Energy boundaries [eV]
      real(8), intent(in) :: opacity_planck_groups(n_groups)    ! κ_P for each group
      real(8), intent(in) :: opacity_rosseland_groups(n_groups) ! κ_R for each group
      real(8), intent(in) :: temperature_eV                     ! Temperature [eV]
      real(8), intent(out) :: grey_planck                       ! Grey κ_P [cm²/g]
      real(8), intent(out) :: grey_rosseland                    ! Grey κ_R [cm²/g]
      
      real(8) :: weight_planck(n_groups)
      real(8) :: weight_rosseland(n_groups)
      real(8) :: total_weight_p, total_weight_r
      real(8) :: E_mid, dE, x
      real(8) :: B_planck, dB_dT
      integer :: i
      
      ! Temperature in energy units
      real(8) :: kT
      kT = temperature_eV
      
      ! Initialize
      grey_planck = 0.0d0
      grey_rosseland = 0.0d0
      total_weight_p = 0.0d0
      total_weight_r = 0.0d0
      
      ! Calculate weights for each energy group
      do i = 1, n_groups
         ! Energy group center and width
         E_mid = 0.5d0 * (energy_grid(i) + energy_grid(i+1))
         dE = energy_grid(i+1) - energy_grid(i)
         
         ! Dimensionless energy
         x = E_mid / kT
         
         ! Calculate Planck function weight
         ! B_ν ∝ ν³/(exp(hν/kT) - 1) ∝ E³/(exp(E/kT) - 1)
         if (x < 50.0d0) then
            B_planck = E_mid**3 / (exp(x) - 1.0d0)
         else
            B_planck = E_mid**3 * exp(-x)  ! Avoid overflow
         endif
         weight_planck(i) = B_planck * dE
         
         ! Calculate Rosseland (derivative) weight
         ! ∂B_ν/∂T ∝ ν⁴ exp(hν/kT)/(kT² (exp(hν/kT)-1)²)
         if (x < 50.0d0) then
            dB_dT = E_mid**4 * exp(x) / (kT**2 * (exp(x) - 1.0d0)**2)
         else
            dB_dT = E_mid**4 * exp(-x) / kT**2
         endif
         weight_rosseland(i) = dB_dT * dE
         
         total_weight_p = total_weight_p + weight_planck(i)
         total_weight_r = total_weight_r + weight_rosseland(i)
      end do
      
      ! Normalize weights
      if (total_weight_p > 0.0d0) then
         weight_planck = weight_planck / total_weight_p
      else
         weight_planck = 1.0d0 / n_groups  ! Equal weights fallback
      endif
      
      if (total_weight_r > 0.0d0) then
         weight_rosseland = weight_rosseland / total_weight_r
      else
         weight_rosseland = 1.0d0 / n_groups
      endif
      
      ! Calculate grey opacities
      ! Planck mean (arithmetic average)
      do i = 1, n_groups
         grey_planck = grey_planck + opacity_planck_groups(i) * weight_planck(i)
      end do
      
      ! Rosseland mean (harmonic average)
      do i = 1, n_groups
         grey_rosseland = grey_rosseland + &
                         weight_rosseland(i) / opacity_rosseland_groups(i)
      end do
      grey_rosseland = 1.0d0 / grey_rosseland
      
   end subroutine calculate_grey_opacity
   
   !================================================================
   ! Alternative: Simple averaging (less accurate but faster)
   !================================================================
   subroutine simple_grey_opacity(n_groups, &
                                 opacity_planck_groups, &
                                 opacity_rosseland_groups, &
                                 grey_planck, grey_rosseland)
      implicit none
      
      integer, intent(in) :: n_groups
      real(8), intent(in) :: opacity_planck_groups(n_groups)
      real(8), intent(in) :: opacity_rosseland_groups(n_groups)
      real(8), intent(out) :: grey_planck, grey_rosseland
      
      integer :: i
      
      ! Simple arithmetic mean for Planck
      grey_planck = sum(opacity_planck_groups) / n_groups
      
      ! Simple harmonic mean for Rosseland
      grey_rosseland = 0.0d0
      do i = 1, n_groups
         grey_rosseland = grey_rosseland + 1.0d0 / opacity_rosseland_groups(i)
      end do
      grey_rosseland = n_groups / grey_rosseland
      
   end subroutine simple_grey_opacity
   
   !================================================================
   ! Print comparison table
   !================================================================
   subroutine print_opacity_comparison(rho, T, &
                                      grey_planck, grey_rosseland, &
                                      multigroup)
      implicit none
      
      real(8), intent(in) :: rho, T
      real(8), intent(in) :: grey_planck, grey_rosseland
      logical, intent(in) :: multigroup
      
      print *, ''
      print *, '===== Opacity Comparison ====='
      print '(A,F8.3,A,F8.1,A)', ' Density: ', rho, ' g/cm³,  Temperature: ', T, ' eV'
      print *, '------------------------------'
      
      if (multigroup) then
         print *, 'Method: 30-group → Grey conversion'
      else
         print *, 'Method: Direct grey calculation'
      endif
      
      print '(A,E12.4,A)', ' Planck mean:    ', grey_planck, ' cm²/g'
      print '(A,E12.4,A)', ' Rosseland mean: ', grey_rosseland, ' cm²/g'
      print '(A,F8.3)', ' Ratio (P/R):    ', grey_planck/grey_rosseland
      print *, ''
      
   end subroutine print_opacity_comparison
   
end module grey_opacity_calculator