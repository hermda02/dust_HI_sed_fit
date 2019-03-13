program dust_hi_fit
  use healpix_types
  use pix_tools
  use fitstools
  use head_fits
  implicit none

  integer(i4b)        :: i, j, l, m, n, niter, npix, nside, nlheader, nmaps, ordering, pics, times, bands
  real(dp)            :: nullval, h, c, k, T, temp, fre, x, z, p, test,T_sum
  real(dp)            :: dust_T_init, dust_T_sigma
  real(dp)            :: missval = -1.6375d30
  logical(lgt)        :: anynull
  logical(lgt)        :: double_precision  

  character(len=128)              :: data, map1, map2, map3, map4, map5, map6, mapHI, mask, output
  character(len=128)              :: arg1, arg2, arg3, fitsname, filename, file1, file2
  character(len=80)               :: line
  character(len=80), dimension(1) :: line2
  character(len=3)                :: number

  real(dp), allocatable, dimension(:,:)      :: HI_temp, HI_mask, masked_HI, T_map
  real(dp), allocatable, dimension(:,:,:)    :: masked, dummy
  real(dp), allocatable, dimension(:)        :: amps, clamps, dummy_T, new_T, amp_std
  real(dp), allocatable, dimension(:)        :: y, freq, sum1, sum2
  character(len=80),  dimension(180)         :: header
  character(len=8),allocatable, dimension(:) :: freqs


  !------------------------------------------------------------------------------------------------------
  ! Daniel Herman 2019                                                                                  |
  !                                                                                                     |  
  ! This program will fit a temperature and an amplitude to the intensity of dust from the HI template. |
  ! We are only looking at high latitudes, where HI column density is < 4e20, because this is where the |
  ! dust/gas ratio is appoximately linear.                                                              |
  !                                                                                                     |  
  ! Start by looking at just 353, 545, 857, 240um, and 100um.                                           |
  !                                                                                                     |  
  !-----------------------------------------------------------------------------------------------------|  

10 format (I1)
11 format (I2)
12 format (I3)
  
  if (iargc() < 3) then
     write(*,*) 'Usage:'
     write(*,*) '   dust_hi_fit [temp estimate (full sky)] [std of temperature (for fitting)] [# of iterations]'
     write(*,*) ''
     stop
  end if

  call getarg(1,arg1)
  read(arg1,*) dust_T_init
  call getarg(2,arg2)
  read(arg2,*) dust_T_sigma
  call getarg(3,arg3)
  read(arg3,*) times

  data  = '../dust_data/sed_data/'
  output= 'results/6bands/'

  map1  = trim(data) // 'npipe6v20_353-5_bmap_QUADCOR_n0064_60arcmin_MJy_calibrated.fits'
  map2  = trim(data) // 'npipe6v20_545-1_bmap_QUADCOR_n0064_60arcmin_MJy_calibrated.fits'
  map3  = trim(data) // 'npipe6v20_857-1_bmap_QUADCOR_n0064_60arcmin_MJy_calibrated.fits'
  map4  = trim(data) // 'DIRBE_240micron_1deg_h064_calibrated.fits'
  map5  = trim(data) // 'DIRBE_140micron_1deg_h064_calibrated.fits'
  map6  = trim(data) // 'DIRBE_100micron_1deg_h064_calibrated.fits'
  mapHI = trim(data) // 'HI_vel_filter_60arcmin_0064.fits'
  mask  = trim(data) // 'HI_mask.fits'

  i     = getsize_fits(mapHI,nside=nside,ordering=ordering,nmaps=nmaps)
  npix  = nside2npix(nside)
  bands = 6

  allocate(masked(0:npix-1,nmaps,bands),dummy(0:npix-1,nmaps,bands),masked_HI(0:npix-1,nmaps),dummy_T(0:npix-1))
  allocate(new_T(0:npix-1),HI_temp(0:npix-1,nmaps),HI_mask(0:npix-1,nmaps),T_map(0:npix-1,nmaps))
  allocate(amps(bands),clamps(bands),freqs(bands),amp_std(bands))
  allocate(y(bands),freq(bands),sum1(bands),sum2(bands))
  
  niter = 1000

  freq(1) = 353.d0
  freq(2) = 545.d0
  freq(3) = 857.d0
  freq(4) = 1249.d0
  freq(5) = 2141.d0
  freq(6) = 2998.d0

  freqs(1) = '353_GHz'
  freqs(2) = '545_GHz'
  freqs(3) = '857_GHz'
  freqs(4) = '1249_GHz'
  freqs(5) = '2141_GHz'
  freqs(6) = '2998_GHz'

  call read_bintab(mapHI, masked_HI, npix, nmaps, nullval, anynull, header=header)
  call read_bintab(mask, HI_mask, npix, nmaps, nullval, anynull, header=header)

  call read_bintab(map1, masked(:,:,1), npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map2, masked(:,:,2), npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map3, masked(:,:,3), npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map4, masked(:,:,4), npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map5, masked(:,:,5), npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map6, masked(:,:,6), npix, nmaps, nullval, anynull, header=header)


  ! Initialize header for writing maps
  nlheader = size(header)
  do i=1,nlheader
     line  = header(i)
     line2 = ''
     if (line(1:8) == 'ORDERING') then
        if (ordering == 1) then
           call add_card(line2, "ORDERING","RING", "Pixel ordering scheme, either RING or NESTED")
        else
           call add_card(line2, "ORDERING","NESTED", "Pixel ordering scheme, either RING or NESTED")
        end if
     end if
     if (line(1:5) == 'NSIDE') then
        call add_card(line2,"NSIDE", nside, "Resolution parameter for HEALPIX")
     end if
     if (line(1:7) == 'LASTPIX') then
        call add_card(line2,"LASTPIX",npix-1,"Last pixel # (0 based)")
     end if
     if (trim(line2(1)) /= '') header(i) = line2(1)
  end do

  ! Initialize dust temperature map
  do i=0,npix-1
     new_T(i) = dust_T_init
  end do

  ! Mask maps
  do j = 1, nmaps
    do i = 0, npix-1
      if (HI_mask(i,j) < 0.5d0) then
        masked(i,:,:)  = missval
        masked_HI(i,:) = missval
        new_T(i)       = missval
      end if
    end do
  end do

  !-----------------------------------------------------------------------------------------------------|  
  ! Here we calculate what the amplitude per map, and temperature per pixel should be for the best fit
  !
  ! with data = I_nu = A_nu * NHI * B_nu(T)
  ! Following the physical model:
  !             I_nu = kappa_nu * r * m_p * NHI * B_nu(T)
  ! So the amplitude we solve for is A_nu = kappa_nu * r * m_p = tau_nu which encapsulates the dust 
  ! emissivity cross section per dust grain per NHI, in units of [cm^2].
  !-----------------------------------------------------------------------------------------------------|  
  
  ! Initializing
  !--------------
  sum1 = 0.d0
  sum2 = 0.d0

  write(*,*) 'Initial amplitudes: '
  write(*,*) '-------------------'
  
  do i=0,npix-1
    if (abs((new_T(i)-missval)/missval) < 1.d-8 .or. abs((masked_HI(i,1)-missval)/missval) < 1.d-8) then
       cycle
    else
       do j=1,bands
          dummy(i,1,j)  = masked_HI(i,1)*planck(freq(j)*1.d9,new_T(i))
          sum1(j) = sum1(j) + (masked(i,1,j)*dummy(i,1,j))
          sum2(j) = sum2(j) + (dummy(i,1,j)**2.d0)
       end do
    endif
  end do

  do j=1,bands
    amps(j) = sum1(j)/sum2(j)
  end do

  amps(1) = 1.0d-7
  amps(2) = 1.0d-6
  amps(3) = 1.0d-6
  amps(4) = 1.0d-5
  amps(5) = 1.0d-5

  clamps = amps

  do n=1,bands
    write(*,*) freqs(n) // ': ', clamps(n)
  end do

  write(*,*) '-------------------------------------'
  write(*,*) ''

  !------------------------------------------------------------
  !------------------------------------------------------------
  filename   = trim(output) // 'amplitudes.dat'
  open(35,file=trim(filename))

  do m=1,times

     write(*,*) 'Iteration', m
     write(*,*) '-----------------------'

     ! Here is where all of the calculations happen
     ! --------------------------------------------------
     do n=1,bands
        amp_std(n) = abs(0.9d0*clamps(n))
     end do

     dummy_T    = create_T(new_T,npix)
     amps       = sample_A(new_T,npix,clamps,amp_std)

     new_T      = dummy_T
     clamps     = amps

     write(*,*) 'Amplitudes: '
     do n=1,bands
        write(*,*) freqs(n) // ': ', clamps(n)
     end do

     ! --------------------------------------------------
     pics  = 0
     T_sum = 0.d0

     do n=0,npix-1
        if (new_T(n) .gt. 10.d0 .and. new_T(n) .lt. 30.d0 ) then
           T_map(n,1) = new_T(n)
           T_sum      = T_sum + T_map(n,1)
           pics       = pics + 1
        else
           T_map(n,1) = missval
        end if
     end do

     write(*,*) 'T = ', T_sum/pics

     ! Write result maps
     if ( mod(m,100) .EQ. 0) then
        if (m .lt. 10) then
           write(number,10) m
        else if (m .gt. 9 .and. m .lt. 100) then
           write(number,11) m
        else if (m .gt. 99) then
           write(number,12) m
        endif

        fitsname   = trim(output) // 'dust_Td_' // trim(number) // '.fits'

        call write_maps(npix,nmaps,header)
     end if
    
     write(35,'(6(E17.8))') clamps(1), clamps(2), clamps(3), clamps(4), clamps(5), clamps(6)
     write(*,*) '-------------------------------------'
     write(*,*) ''

  end do

  deallocate(masked,masked_HI)
  deallocate(dummy)
  close(35)
  

contains

  function planck(fre,T)
        implicit none
        real(dp), intent(in)  :: fre
        real(dp), intent(in)  :: T
        real(dp)              :: h = 6.626d-34
        real(dp)              :: c = 3.0d8 
        real(dp)              :: k = 1.38d-23
        real(dp)              :: planck
        ! Output in units of [W sr^-1 m^-2 Hz^-1]
        planck  = ((2.d0*h*fre**3.d0)/(c**2.d0))*(1.d0/(exp((h*fre)/(k*T))-1))
  end function planck


  function rand_normal(mean,stdev) result(c)
       double precision :: mean,stdev,c,temp(2),theta,r
       if (stdev <= 0.0d0) then
          write(*,*) "Standard Deviation must be positive."
       else
          call RANDOM_NUMBER(temp)
          r=(-2.0d0*log(temp(1)))**0.5
          theta = 2.0d0*PI*temp(2)
          c= mean+stdev*r*sin(theta)
    end if
  end function


  function sample_A(T,npix,amp,std)
    implicit none

    integer(i4b), intent(in)                  :: npix
    real(dp), dimension(0:npix-1), intent(in) :: T
    real(dp), dimension(bands), intent(in)    :: amp, std
    real(dp), dimension(bands)                :: chisq1, sample_A, tau, r, sigma, mu
    real(dp), dimension(npix,bands)           :: vals
    real(dp), dimension(2)                    :: a
    real(dp)                                  :: chisq2, p, num
    integer(i4b)                              :: b

    tau  = amp
    a(1) = 1.d0

    chisq1 = 0.d0
    chisq2 = 0.d0

    b = 0
    if (b .EQ. 0) then
       do i=0,npix-1
          if (abs((T(i)-missval)/missval) < 1.d-8 .or. abs((masked_HI(i,1)-missval)/missval) < 1.d-8) then
             cycle
          else
             b = b + 1
             do j=1,bands
                dummy(i,1,j) = masked_HI(i,1)*planck(freq(j)*1.d9,T(i))
                chisq1(j)    = chisq1(j) + (masked(i,1,j) - tau(j)*dummy(i,1,j))**2.d0
                vals(b,j)    = masked(i,1,j)
             end do
          endif
       end do
    end if

    ! Compute standard deviation for the chisq
    do j=1,bands
       mu(j)     = sum(vals(1:b,j))/b
       sigma(j)  = sqrt(sum((vals(1:b,j)-mu(j))**2.d0)/(b-1))
       chisq1(j) = chisq1(j)/(sigma(j)**2.d0)
    end do


    do j=1,bands
       ! write(*,*) freqs(j) // ':'
       do n=1,niter
          chisq2 = 0.d0
          r(j) = rand_normal(tau(j),std(j))
          do i=0,npix-1
             if (abs((T(i)-missval)/missval) < 1.d-8 .or. abs((masked_HI(i,1)-missval)/missval) < 1.d-8) then
                cycle
             else
                dummy(i,1,j) = masked_HI(i,1)*planck(freq(j)*1.d9,T(i))
                chisq2       = chisq2 + ((masked(i,1,j) - r(j)*dummy(i,1,j))**2.d0)/(sigma(j)**2.d0)
             end if
          end do

          a(2) = chisq2/chisq1(j)
          p    = minval(a)
          call RANDOM_NUMBER(num)

          if (num > p) then
             ! write(*,*) p
             if (r(j) .gt. 0.d0 .and. r(j) .lt. 1.d-3) then
                tau(j)    = r(j)
                chisq1(j) = chisq2
             end if
          end if
       end do
    end do

    write(*,*) 'Chi-square values :'
    do j=1,bands
      write(*,*) freqs(j) // ': ', chisq1(j)
    end do
    write(*,*) ''

    sample_A = tau

  end function sample_A


  function sample_T(x,y,z,T,sigma,amp)
    implicit none

    real(dp), intent(in)                   :: x
    real(dp), dimension(bands), intent(in) :: y
    real(dp), intent(in)                   :: z
    real(dp), dimension(bands), intent(in) :: amp
    real(dp), intent(in)                   :: T
    real(dp), intent(in)                   :: sigma
    real(dp), dimension(2)                 :: a
    real(dp), dimension(bands)             :: test
    real(dp)                               :: sample_T, temp, r, p, b, c, num

    temp = T
    a(1) = 1.d0
    c    = x

    do i=1,niter
       r    = rand_normal(temp,sigma)
       test = 0.d0

       do j=1,bands
          test(j) = amp(j)*z*planck(freq(j)*1.d9,r)
       end do

       b    = sum((test-y)**2.d0)
       a(2) = b/c
       p    = minval(a)

       call RANDOM_NUMBER(num)

       if (num > p) then
          if (r .gt. 15.d0 .and. r .lt. 30.d0 .and. abs(temp - r) .gt. 0.001d0) then
             temp = r
             c    = b
          end if
       end if
    end do

    sample_T = temp

  end function sample_T


  function create_T(T,npix)
    implicit none

    integer(i4b), intent(in)                     :: npix
    real(dp), dimension(0:npix-1), intent(in)    :: T
    real(dp), dimension(0:npix-1)                :: create_T
    real(dp), dimension(bands)                   :: y
    real(dp)                                     :: z
    real(dp)                                     :: x

    do l=0,npix-1
       if (abs((T(l)-missval)/missval) < 1.d-8 .or. abs((masked_HI(l,1)-missval)/missval) < 1.d-8) then
          cycle
       else
          do j=1,bands
             y(j) = masked(l,1,j)
          end do
          z    = masked_HI(l,1)
          x    = sum((clamps(:)*dummy(l,1,:) - y(:))**2.d0)
          create_T(l) = sample_T(x,y,z,T(l),dust_T_sigma,clamps)
       endif
    end do
  end function create_T

  subroutine write_maps(npix,nmaps,header)
    ! Outputs the temperature map, the modeled map, and the residual map
    integer(i4b), intent(in)                        :: npix, nmaps
    integer(i4b)                                    :: y, b
    character(len=80),  dimension(180), intent(in)  :: header
    character(len=80)                               :: file1, file2, file3
    real(dp), dimension(0:npix-1,nmaps,bands)           :: model, resid

    file1 = trim(output) // 'dust_Td_' // trim(number) // '.fits'

    do y=0,npix-1
       model(y,1,:) = amps(1)*dummy(y,1,:)
       resid(y,1,:) = masked(y,1,:) - model(y,1,:)
    end do

    call write_bintab(T_map, npix, nmaps, header, nlheader, file1)
    do y=1,bands
       file2 = trim(output) // 'model_'// trim(freqs(y)) // '_' // trim(number) // '.fits'
       file3 = trim(output) // 'resid_'// trim(freqs(y)) // '_' // trim(number) // '.fits'
       call write_bintab(model(:,:,y), npix, nmaps, header, nlheader, file2)
       call write_bintab(resid(:,:,y), npix, nmaps, header, nlheader, file3)
    end do

  end subroutine write_maps

end program dust_hi_fit
