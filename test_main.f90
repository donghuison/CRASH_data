program test_crash_opacity
   use crash_opacity_reader
   implicit none

   type(crash_opacity_table) :: au_table
   real(8) :: rho_cgs, T_ev, kappa_p, kappa_r

   ! Gold opacity 테이블 로드
   call load_crash_opacity('Au', au_table)

   ! 테스트 조건
   rho_cgs = 1.0    ! g/cm3
   T_ev = 10.0      ! eV

   ! Opacity 계산 (Planck & Rosseland 모두 얻음!)
   call get_crash_opacity(au_table, rho_cgs, T_ev, kappa_p, kappa_r)

   print *, 'Density:', rho_cgs, ' g/cm3'
   print *, 'Temperature:', T_ev, ' eV'
   print *, 'Planck mean opacity:', kappa_p, ' cm2/g'
   print *, 'Rosseland mean opacity:', kappa_r, ' cm2/g'

end program

