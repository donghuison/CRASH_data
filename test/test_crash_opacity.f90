program test_crash_opacity
   use crash_opacity_reader
   implicit none
   
   type(crash_opacity_table) :: au_table, xe_table, be_table
   real(8) :: rho_cgs, T_ev, kappa_p, kappa_r
   real(8) :: kappa_p_groups(30)
   real(8) :: rho_values(5), T_values(5)
   integer :: i, j, k
   character(len=20) :: materials(3)
   
   print *, '================================================'
   print *, ' CRASH Opacity Table Test Program'
   print *, '================================================'
   print *, ''
   
   ! Test materials
   materials = (/'Au', 'Xe', 'Be'/)
   
   ! Test conditions
   rho_values = (/0.001d0, 0.01d0, 0.1d0, 1.0d0, 10.0d0/)  ! g/cm3
   T_values = (/1.0d0, 5.0d0, 10.0d0, 50.0d0, 100.0d0/)    ! eV
   
   ! Test 1: Load opacity tables for different materials
   print *, '=== Test 1: Loading Opacity Tables ==='
   print *, ''
   
   call load_crash_opacity('Au', au_table)
   print *, ''
   call load_crash_opacity('Xe', xe_table)
   print *, ''
   call load_crash_opacity('Be', be_table)
   print *, ''
   
   ! Test 2: Calculate opacity at various conditions for Gold
   print *, '=== Test 2: Gold Opacity at Various Conditions ==='
   print *, ''
   print *, 'Material: Au (Gold)'
   print *, '-------------------------------------------'
   print '(A8,A10,A12,A12)', 'rho[g/cc]', 'T[eV]', 'kappa_P', 'kappa_R'
   print *, '-------------------------------------------'
   
   do i = 1, 5
      do j = 1, 5
         rho_cgs = rho_values(i)
         T_ev = T_values(j)
         
         call get_crash_opacity(au_table, rho_cgs, T_ev, kappa_p, kappa_r)
         
         print '(F8.3,F10.1,2E12.4)', rho_cgs, T_ev, kappa_p, kappa_r
      end do
   end do
   print *, ''
   
   ! Test 3: Compare opacity between materials at same conditions
   print *, '=== Test 3: Material Comparison ==='
   print *, ''
   print *, 'Conditions: rho = 1.0 g/cm3, T = 10.0 eV'
   print *, '-------------------------------------------'
   print '(A10,A12,A12,A8)', 'Material', 'kappa_P', 'kappa_R', 'P/R ratio'
   print *, '-------------------------------------------'
   
   rho_cgs = 1.0d0
   T_ev = 10.0d0
   
   ! Gold
   call get_crash_opacity(au_table, rho_cgs, T_ev, kappa_p, kappa_r)
   print '(A10,2E12.4,F8.3)', 'Gold', kappa_p, kappa_r, kappa_p/kappa_r
   
   ! Xenon
   call get_crash_opacity(xe_table, rho_cgs, T_ev, kappa_p, kappa_r)
   print '(A10,2E12.4,F8.3)', 'Xenon', kappa_p, kappa_r, kappa_p/kappa_r
   
   ! Beryllium
   call get_crash_opacity(be_table, rho_cgs, T_ev, kappa_p, kappa_r)
   print '(A10,2E12.4,F8.3)', 'Beryllium', kappa_p, kappa_r, kappa_p/kappa_r
   print *, ''
   
   ! Test 4: Energy group analysis
   print *, '=== Test 4: Energy Group Analysis for Au ==='
   print *, ''
   print *, 'Conditions: rho = 1.0 g/cm3, T = 10.0 eV'
   print *, '-------------------------------------------'
   
   rho_cgs = 1.0d0
   T_ev = 10.0d0
   
   call get_crash_opacity(au_table, rho_cgs, T_ev, kappa_p, kappa_r, kappa_p_groups)
   
   print '(A6,A12,A12,A12)', 'Group', 'E_low[eV]', 'E_high[eV]', 'kappa_P'
   print *, '-------------------------------------------'
   
   do k = 1, 30
      print '(I6,2F12.2,E12.4)', k, au_table%energy_grid(k), &
             au_table%energy_grid(k+1), kappa_p_groups(k)
   end do
   
   print *, ''
   print *, 'Weighted average kappa_P = ', kappa_p, ' cm2/g'
   print *, ''
   
   ! Test 5: Temperature dependence
   print *, '=== Test 5: Temperature Dependence ==='
   print *, ''
   print *, 'Material: Au, rho = 1.0 g/cm3'
   print *, '-------------------------------------------'
   print '(A10,A12,A12,A8)', 'T[eV]', 'kappa_P', 'kappa_R', 'P/R ratio'
   print *, '-------------------------------------------'
   
   rho_cgs = 1.0d0
   do i = 1, 10
      T_ev = 0.5d0 * 2.0d0**i  ! 1, 2, 4, 8, 16, 32, 64, 128, 256, 512 eV
      if (T_ev > 200.0d0) exit  ! Stay within table range
      
      call get_crash_opacity(au_table, rho_cgs, T_ev, kappa_p, kappa_r)
      print '(F10.1,2E12.4,F8.3)', T_ev, kappa_p, kappa_r, kappa_p/kappa_r
   end do
   
   print *, ''
   print *, '================================================'
   print *, ' Test completed successfully!'
   print *, '================================================'
   
end program test_crash_opacity