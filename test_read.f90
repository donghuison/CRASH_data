module crash_opacity_reader
   implicit none

   type :: crash_opacity_table
      ! 그리드 정보
      integer :: n_rho = 201
      integer :: n_temp = 201
      integer :: n_groups = 30

      ! 범위 (BATSRUS TABLE_*.in 파일 기준)
      real(8) :: log_rho_min, log_rho_max  ! log10(kg/m3)
      real(8) :: log_temp_min, log_temp_max ! log10(eV)

      ! 에너지 그룹 경계 [eV]
      real(8) :: energy_grid(31)

      ! Opacity 데이터 [cm2/g]
      real(8), allocatable :: planck_opacity(:,:,:)   ! (201,201,30)
      real(8), allocatable :: rosseland_opacity(:,:,:) ! (201,201,30)

      character(len=100) :: material_name
   end type

contains

   !================================================================
   ! CRASH 테이블 로드 (~/CRASH_data/LookupTables 사용)
   !================================================================
   subroutine load_crash_opacity(material, table)
      implicit none
      character(len=*), intent(in) :: material ! 'Au', 'Xe', 'Be', 'Pl', 'Ay'
      type(crash_opacity_table), intent(out) :: table

      character(len=200) :: filename, energy_file
      integer :: unit, i, j, k
      real(8) :: dummy
      logical :: file_exists

      ! 파일 경로 설정
      filename = trim(adjustl(getenv('HOME')))//'/CRASH_data/LookupTables/'// &
         trim(material)//'_opac_CRASH.dat'
      energy_file =
      trim(adjustl(getenv('HOME')))//'/CRASH_data/LookupTables/EnergyGrid30.dat'

      ! 파일 존재 확인
      inquire(file=filename, exist=file_exists)
      if (.not. file_exists) then
         print *, 'ERROR: Cannot find ', trim(filename)
         stop
      endif

      table%material_name = material

      ! 물질별 범위 설정 (TABLE_*_OPACITY.in 기준)
      select case(material)
       case('Au')
         table%log_rho_min = log10(0.327)
         table%log_rho_max = log10(39248.0)
         table%log_temp_min = log10(0.2)
         table%log_temp_max = log10(100.0)
       case('Xe')
         table%log_rho_min = log10(0.218)
         table%log_rho_max = log10(21801.0)
         table%log_temp_min = log10(0.03)
         table%log_temp_max = log10(1000.0)
       case('Be')
         table%log_rho_min = log10(0.185)
         table%log_rho_max = log10(18500.0)
         table%log_temp_min = log10(0.2)
         table%log_temp_max = log10(100.0)
       case default
         ! 기본값
         table%log_rho_min = -1.0
         table%log_rho_max = 5.0
         table%log_temp_min = -1.0
         table%log_temp_max = 3.0
      end select

      ! 에너지 그리드 읽기
      open(newunit=unit, file=energy_file, status='old')
      do i = 1, 31
         read(unit,*) table%energy_grid(i)
      end do
      close(unit)

      ! 메모리 할당
      allocate(table%planck_opacity(201, 201, 30))
      allocate(table%rosseland_opacity(201, 201, 30))

      ! 바이너리 파일 읽기
      open(newunit=unit, file=filename, status='old', &
         form='unformatted', access='stream')

      ! 헤더 건너뛰기 (있을 경우)
      ! CRASH 형식은 직접 데이터부터 시작

      do i = 1, 201
         do j = 1, 201
            ! 30개 Planck mean opacity 읽기
            read(unit) (table%planck_opacity(i,j,k), k=1,30)
            ! 30개 Rosseland mean opacity 읽기
            read(unit) (table%rosseland_opacity(i,j,k), k=1,30)
         end do
      end do

      close(unit)

      print *, 'Loaded CRASH opacity table for ', trim(material)
      print *, '  Density range: 10^', table%log_rho_min, ' - 10^', &
         table%log_rho_max, ' kg/m3'
      print *, '  Temperature range: 10^', table%log_temp_min, ' - 10^', &
         table%log_temp_max, ' eV'

   end subroutine load_crash_opacity

   !================================================================
   ! Opacity 보간 (Planck & Rosseland)
   !================================================================
   subroutine get_crash_opacity(table, rho_cgs, T_ev, &
      kappa_planck, kappa_rosseland)
      implicit none
      type(crash_opacity_table), intent(in) :: table
      real(8), intent(in) :: rho_cgs  ! [g/cm3]
      real(8), intent(in) :: T_ev     ! [eV]
      real(8), intent(out) :: kappa_planck    ! [cm2/g]
      real(8), intent(out) :: kappa_rosseland ! [cm2/g]

      real(8) :: log_rho, log_T
      real(8) :: x, y, dx, dy
      integer :: i, j, k
      real(8) :: planck_groups(30), ross_groups(30)
      real(8) :: planck_weight, total_weight
      real(8) :: E_mid, kT

      ! 단위 변환 및 로그 변환
      log_rho = log10(rho_cgs * 1000.0)  ! g/cm3 → kg/m3
      log_T = log10(T_ev)

      ! 범위 제한
      log_rho = max(table%log_rho_min, min(table%log_rho_max, log_rho))
      log_T = max(table%log_temp_min, min(table%log_temp_max, log_T))

      ! 그리드 위치 계산
      x = (log_rho - table%log_rho_min) / &
         (table%log_rho_max - table%log_rho_min) * 200.0
      y = (log_T - table%log_temp_min) / &
         (table%log_temp_max - table%log_temp_min) * 200.0

      i = int(x) + 1
      j = int(y) + 1
      dx = x - (i-1)
      dy = y - (j-1)

      ! 경계 처리
      i = max(1, min(200, i))
      j = max(1, min(200, j))

      ! 각 에너지 그룹에 대해 이중선형 보간
      do k = 1, 30
         planck_groups(k) = &
            (1-dx)*(1-dy) * table%planck_opacity(i,j,k) + &
            dx*(1-dy) * table%planck_opacity(i+1,j,k) + &
            (1-dx)*dy * table%planck_opacity(i,j+1,k) + &
            dx*dy * table%planck_opacity(i+1,j+1,k)

         ross_groups(k) = &
            (1-dx)*(1-dy) * table%rosseland_opacity(i,j,k) + &
            dx*(1-dy) * table%rosseland_opacity(i+1,j,k) + &
            (1-dx)*dy * table%rosseland_opacity(i,j+1,k) + &
            dx*dy * table%rosseland_opacity(i+1,j+1,k)
      end do

      ! Planck 가중 평균 계산
      kT = T_ev
      kappa_planck = 0.0
      total_weight = 0.0

      do k = 1, 30
         E_mid = 0.5 * (table%energy_grid(k) + table%energy_grid(k+1))

         ! Planck 가중 함수
         if (E_mid/kT < 50.0) then
            planck_weight = E_mid**3 / (exp(E_mid/kT) - 1.0)
         else
            planck_weight = E_mid**3 * exp(-E_mid/kT)
         endif

         planck_weight = planck_weight * &
            (table%energy_grid(k+1) - table%energy_grid(k))

         kappa_planck = kappa_planck + planck_groups(k) * planck_weight
         total_weight = total_weight + planck_weight
      end do

      if (total_weight > 0.0) then
         kappa_planck = kappa_planck / total_weight
      else
         kappa_planck = sum(planck_groups) / 30.0
      endif

      ! Rosseland mean (조화평균)
      kappa_rosseland = 0.0
      do k = 1, 30
         kappa_rosseland = kappa_rosseland + &
            (table%energy_grid(k+1) - table%energy_grid(k)) / ross_groups(k)
      end do
      kappa_rosseland = (table%energy_grid(31) - table%energy_grid(1)) / &
         kappa_rosseland

   end subroutine get_crash_opacity

end module crash_opacity_reader
