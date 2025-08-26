program test_grey_opacity
   use crash_opacity_reader
   use grey_opacity_calculator
   implicit none
   
   type(crash_opacity_table) :: au_table
   real(8) :: rho_cgs, T_ev
   real(8) :: kappa_p_groups(30), kappa_r_groups(30)
   real(8) :: grey_planck, grey_rosseland
   real(8) :: grey_planck_simple, grey_rosseland_simple
   real(8) :: energy_grid(31)
   integer :: i
   
   print *, '================================================'
   print *, ' Grey Opacity Calculator Test'
   print *, '================================================'
   print *, ''
   print *, 'This program demonstrates how to convert'
   print *, '30-group opacity data to grey opacity κ(ρ,T)'
   print *, ''
   
   ! Load CRASH opacity table for Gold
   call load_crash_opacity('Au', au_table)
   
   ! Copy energy grid
   energy_grid = au_table%energy_grid
   
   print *, ''
   print *, '=== Test 1: Different Temperature Conditions ==='
   print *, ''
   
   rho_cgs = 1.0d0  ! Fixed density
   
   ! Test at different temperatures
   do i = 1, 5
      select case(i)
      case(1)
         T_ev = 1.0d0
      case(2)
         T_ev = 10.0d0
      case(3)
         T_ev = 50.0d0
      case(4)
         T_ev = 100.0d0
      case(5)
         T_ev = 200.0d0
      end select
      
      ! Get 30-group opacities
      call get_crash_opacity(au_table, rho_cgs, T_ev, &
                            grey_planck, grey_rosseland, kappa_p_groups)
      kappa_r_groups = kappa_p_groups / 1.665d0  ! Approximate for test
      
      ! Method 1: Accurate calculation with Planck weighting
      call calculate_grey_opacity(30, energy_grid, &
                                 kappa_p_groups, kappa_r_groups, &
                                 T_ev, grey_planck, grey_rosseland)
      
      print '(A,F6.1,A)', '--- Temperature = ', T_ev, ' eV ---'
      print '(A,E12.4,A)', 'Grey Planck:    ', grey_planck, ' cm²/g'
      print '(A,E12.4,A)', 'Grey Rosseland: ', grey_rosseland, ' cm²/g'
      print '(A,F8.3)', 'Ratio (P/R):    ', grey_planck/grey_rosseland
      
      ! Method 2: Simple averaging
      call simple_grey_opacity(30, kappa_p_groups, kappa_r_groups, &
                              grey_planck_simple, grey_rosseland_simple)
      
      print '(A,E12.4,A)', 'Simple Planck:  ', grey_planck_simple, ' cm²/g'
      print '(A,E12.4,A)', 'Simple Rossel.: ', grey_rosseland_simple, ' cm²/g'
      print '(A,F8.3)', 'Accuracy ratio: ', grey_planck/grey_planck_simple
      print *, ''
   end do
   
   print *, '=== Test 2: Energy Group Weight Analysis ==='
   print *, ''
   print *, 'Showing how different energy groups contribute'
   print *, 'to the grey opacity at T = 10 eV:'
   print *, ''
   
   T_ev = 10.0d0
   call get_crash_opacity(au_table, rho_cgs, T_ev, &
                         grey_planck, grey_rosseland, kappa_p_groups)
   kappa_r_groups = kappa_p_groups / 1.665d0
   
   ! Calculate weights
   print '(A5,A12,A12,A12,A12)', 'Group', 'E_mid[eV]', 'κ_P[cm²/g]', &
          'Weight_P', 'Contribution'
   print *, '--------------------------------------------------------'
   
   do i = 1, 30
      block
         real(8) :: E_mid, x, weight
         E_mid = 0.5d0 * (energy_grid(i) + energy_grid(i+1))
         x = E_mid / T_ev
         
         if (x < 50.0d0) then
            weight = E_mid**3 / (exp(x) - 1.0d0) * &
                    (energy_grid(i+1) - energy_grid(i))
         else
            weight = E_mid**3 * exp(-x) * &
                    (energy_grid(i+1) - energy_grid(i))
         endif
         
         ! Print only significant contributors
         if (weight > 0.001d0) then
            print '(I5,F12.2,E12.4,F12.6,E12.4)', i, E_mid, &
                   kappa_p_groups(i), weight, kappa_p_groups(i)*weight
         endif
      end block
   end do
   
   print *, ''
   print *, '=== Test 3: Grey vs Multigroup Comparison ==='
   print *, ''
   print *, 'Conditions: ρ = 1.0 g/cm³, T = 50 eV'
   
   rho_cgs = 1.0d0
   T_ev = 50.0d0
   
   ! Get multigroup opacities
   call get_crash_opacity(au_table, rho_cgs, T_ev, &
                         grey_planck, grey_rosseland, kappa_p_groups)
   kappa_r_groups = kappa_p_groups / 1.665d0
   
   ! Calculate grey opacities
   call calculate_grey_opacity(30, energy_grid, &
                              kappa_p_groups, kappa_r_groups, &
                              T_ev, grey_planck, grey_rosseland)
   
   print *, ''
   print *, 'Multigroup (30 groups):'
   print '(A,A)', '  Energy range: 0.1 - 20,000 eV'
   print '(A,A)', '  Full spectral resolution'
   print '(A,30F7.1)', '  Groups [eV]: ', (energy_grid(i), i=1,10)
   
   print *, ''
   print *, 'Grey approximation (single value):'
   print '(A,E12.4,A)', '  κ_P(ρ,T) = ', grey_planck, ' cm²/g'
   print '(A,E12.4,A)', '  κ_R(ρ,T) = ', grey_rosseland, ' cm²/g'
   
   print *, ''
   print *, '================================================'
   print *, ' Summary: Grey Opacity Usage Guidelines'
   print *, '================================================'
   print *, ''
   print *, '1. Use Planck mean for:'
   print *, '   - Emission and absorption calculations'
   print *, '   - Optically thin regions'
   print *, '   - Energy balance equations'
   print *, ''
   print *, '2. Use Rosseland mean for:'
   print *, '   - Radiation diffusion'
   print *, '   - Optically thick regions'
   print *, '   - Flux calculations'
   print *, ''
   print *, '3. Typical ratio: κ_P/κ_R ≈ 1.5-2.0'
   print *, '   (Planck > Rosseland always)'
   print *, ''
   print *, '================================================'
   
end program test_grey_opacity