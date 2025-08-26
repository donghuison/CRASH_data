module crash_opacity_reader
   implicit none
   
   type :: crash_opacity_table
      ! Grid information
      integer :: n_rho = 201
      integer :: n_temp = 201
      integer :: n_groups = 30
      
      ! Range (log scale)
      real(8) :: log_rho_min, log_rho_max  ! log10(kg/m3)
      real(8) :: log_temp_min, log_temp_max ! log10(eV)
      
      ! Energy group boundaries [eV]
      real(8) :: energy_grid(31)
      
      ! Opacity data - store as log10 values
      real(8), allocatable :: planck_opacity(:,:,:)   ! (201,201,30)
      real(8), allocatable :: rosseland_opacity(:,:,:) ! (201,201,30)
      
      ! Grid arrays
      real(8), allocatable :: log_rho_grid(:)
      real(8), allocatable :: log_temp_grid(:)
      
      character(len=100) :: material_name
   end type
   
contains

   !================================================================
   ! Load CRASH opacity table from CRASH_data/LookupTables
   !================================================================
   subroutine load_crash_opacity(material, table)
      implicit none
      character(len=*), intent(in) :: material ! 'Au', 'Xe', 'Be', 'Pl', 'Ay'
      type(crash_opacity_table), intent(out) :: table
      
      character(len=200) :: param_file, energy_file
      integer :: unit, i, j, k
      logical :: file_exists
      character(len=100) :: line
      
      ! Set file paths
      param_file = '../LookupTables/../Param/CRASH/TABLE_'// &
                   trim(material)//'_OPACITY.in'
      energy_file = '../LookupTables/EnergyGrid30.dat'
      
      table%material_name = material
      
      ! Read parameter ranges from TABLE file
      ! For simplicity, use hardcoded values based on actual files
      select case(trim(material))
      case('Au')
         table%log_rho_min = log10(0.3270688818768d0)
         table%log_rho_max = log10(3.9248265825216d4)
         table%log_temp_min = log10(0.2d0)
         table%log_temp_max = log10(1.0d2)
      case('Xe')
         table%log_rho_min = log10(0.21801430476d0)
         table%log_rho_max = log10(2.1801430476d4)
         table%log_temp_min = log10(0.03d0)
         table%log_temp_max = log10(1.0d3)
      case('Be')
         table%log_rho_min = log10(0.1849d0)
         table%log_rho_max = log10(1.849d4)
         table%log_temp_min = log10(0.2d0)
         table%log_temp_max = log10(1.0d2)
      case('Pl')
         table%log_rho_min = log10(0.125d0)
         table%log_rho_max = log10(1.25d4)
         table%log_temp_min = log10(0.2d0)
         table%log_temp_max = log10(1.0d2)
      case('Ay')
         table%log_rho_min = log10(0.12d0)
         table%log_rho_max = log10(1.2d4)
         table%log_temp_min = log10(0.2d0)
         table%log_temp_max = log10(1.0d2)
      case default
         print *, 'ERROR: Unknown material ', material
         stop
      end select
      
      ! Read energy grid
      inquire(file=energy_file, exist=file_exists)
      if (.not. file_exists) then
         print *, 'ERROR: Cannot find ', trim(energy_file)
         stop
      endif
      
      open(newunit=unit, file=energy_file, status='old')
      do i = 1, 31
         read(unit,*) table%energy_grid(i)
      end do
      close(unit)
      
      ! Allocate memory
      allocate(table%planck_opacity(201, 201, 30))
      allocate(table%rosseland_opacity(201, 201, 30))
      allocate(table%log_rho_grid(201))
      allocate(table%log_temp_grid(201))
      
      ! Create grid arrays
      do i = 1, 201
         table%log_rho_grid(i) = table%log_rho_min + &
            (table%log_rho_max - table%log_rho_min) * (i-1) / 200.0
         table%log_temp_grid(i) = table%log_temp_min + &
            (table%log_temp_max - table%log_temp_min) * (i-1) / 200.0
      end do
      
      ! For testing, generate synthetic opacity data
      ! In real case, would read from binary .dat files
      print *, 'Generating synthetic opacity data for testing...'
      do i = 1, 201
         do j = 1, 201
            do k = 1, 30
               ! Simple model: opacity depends on density, temperature, and energy
               table%planck_opacity(i,j,k) = &
                  3.0 + 0.5*table%log_rho_grid(i) - 2.0*table%log_temp_grid(j) &
                  + 0.1*log10(table%energy_grid(k))
               
               table%rosseland_opacity(i,j,k) = &
                  2.8 + 0.5*table%log_rho_grid(i) - 2.0*table%log_temp_grid(j) &
                  + 0.08*log10(table%energy_grid(k))
            end do
         end do
      end do
      
      print *, 'Loaded opacity table for ', trim(material)
      print *, '  Density range: 10^', table%log_rho_min, ' - 10^', &
               table%log_rho_max, ' kg/m3'
      print *, '  Temperature range: 10^', table%log_temp_min, ' - 10^', &
               table%log_temp_max, ' eV'
      
   end subroutine load_crash_opacity
   
   !================================================================
   ! Interpolate opacity (Planck & Rosseland)
   !================================================================
   subroutine get_crash_opacity(table, rho_cgs, T_ev, &
                                kappa_planck, kappa_rosseland, &
                                kappa_planck_groups)
      implicit none
      type(crash_opacity_table), intent(in) :: table
      real(8), intent(in) :: rho_cgs  ! [g/cm3]
      real(8), intent(in) :: T_ev     ! [eV]
      real(8), intent(out) :: kappa_planck    ! [cm2/g]
      real(8), intent(out) :: kappa_rosseland ! [cm2/g]
      real(8), intent(out), optional :: kappa_planck_groups(30)
      
      real(8) :: log_rho, log_T
      real(8) :: x, y, dx, dy
      integer :: i, j, k
      real(8) :: planck_groups(30), ross_groups(30)
      real(8) :: planck_weight, total_weight
      real(8) :: E_mid, kT, dE
      real(8) :: f00, f10, f01, f11
      
      ! Unit conversion and log transform
      log_rho = log10(rho_cgs * 1000.0d0)  ! g/cm3 → kg/m3
      log_T = log10(T_ev)
      
      ! Limit to valid range
      log_rho = max(table%log_rho_min, min(table%log_rho_max, log_rho))
      log_T = max(table%log_temp_min, min(table%log_temp_max, log_T))
      
      ! Find grid position
      i = 1
      do while (i < 201 .and. table%log_rho_grid(i) < log_rho)
         i = i + 1
      end do
      i = max(1, min(200, i-1))
      
      j = 1
      do while (j < 201 .and. table%log_temp_grid(j) < log_T)
         j = j + 1
      end do
      j = max(1, min(200, j-1))
      
      ! Calculate interpolation weights
      dx = (log_rho - table%log_rho_grid(i)) / &
           (table%log_rho_grid(i+1) - table%log_rho_grid(i))
      dy = (log_T - table%log_temp_grid(j)) / &
           (table%log_temp_grid(j+1) - table%log_temp_grid(j))
      
      dx = max(0.0d0, min(1.0d0, dx))
      dy = max(0.0d0, min(1.0d0, dy))
      
      ! Bilinear interpolation weights
      f00 = (1.0d0 - dx) * (1.0d0 - dy)
      f10 = dx * (1.0d0 - dy)
      f01 = (1.0d0 - dx) * dy
      f11 = dx * dy
      
      ! Interpolate for each energy group
      do k = 1, 30
         planck_groups(k) = &
            f00 * table%planck_opacity(i,j,k) + &
            f10 * table%planck_opacity(i+1,j,k) + &
            f01 * table%planck_opacity(i,j+1,k) + &
            f11 * table%planck_opacity(i+1,j+1,k)
         
         ross_groups(k) = &
            f00 * table%rosseland_opacity(i,j,k) + &
            f10 * table%rosseland_opacity(i+1,j,k) + &
            f01 * table%rosseland_opacity(i,j+1,k) + &
            f11 * table%rosseland_opacity(i+1,j+1,k)
         
         ! Convert from log to linear
         planck_groups(k) = 10.0d0**planck_groups(k)
         ross_groups(k) = 10.0d0**ross_groups(k)
      end do
      
      ! Calculate Planck-weighted average
      kT = T_ev
      kappa_planck = 0.0d0
      total_weight = 0.0d0
      
      do k = 1, 30
         E_mid = sqrt(table%energy_grid(k) * table%energy_grid(k+1))
         dE = table%energy_grid(k+1) - table%energy_grid(k)
         
         ! Planck weight function
         if (E_mid/kT < 50.0d0) then
            planck_weight = E_mid**3 / (exp(E_mid/kT) - 1.0d0)
         else
            planck_weight = E_mid**3 * exp(-E_mid/kT)
         endif
         
         planck_weight = planck_weight * dE
         
         kappa_planck = kappa_planck + planck_groups(k) * planck_weight
         total_weight = total_weight + planck_weight
      end do
      
      if (total_weight > 0.0d0) then
         kappa_planck = kappa_planck / total_weight
      else
         kappa_planck = sum(planck_groups) / 30.0d0
      endif
      
      ! Rosseland mean (harmonic mean with derivative weighting)
      kappa_rosseland = 0.0d0
      total_weight = 0.0d0
      
      do k = 1, 30
         E_mid = sqrt(table%energy_grid(k) * table%energy_grid(k+1))
         dE = table%energy_grid(k+1) - table%energy_grid(k)
         
         ! Rosseland weight: derivative of Planck function
         if (E_mid/kT < 50.0d0) then
            planck_weight = E_mid**4 * exp(E_mid/kT) / &
                           (kT**2 * (exp(E_mid/kT) - 1.0d0)**2)
         else
            planck_weight = E_mid**4 * exp(-E_mid/kT) / kT**2
         endif
         
         planck_weight = planck_weight * dE
         
         kappa_rosseland = kappa_rosseland + planck_weight / ross_groups(k)
         total_weight = total_weight + planck_weight
      end do
      
      if (total_weight > 0.0d0) then
         kappa_rosseland = total_weight / kappa_rosseland
      else
         kappa_rosseland = 30.0d0 / sum(1.0d0/ross_groups)
      endif
      
      ! Return group values if requested
      if (present(kappa_planck_groups)) then
         kappa_planck_groups = planck_groups
      endif
      
   end subroutine get_crash_opacity
   
end module crash_opacity_reader