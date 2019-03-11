program dust_hi_fit
  use healpix_types
  use pix_tools
  use fitstools
  use head_fits
  implicit none

  integer(i4b)        :: i, j, l, m, n, niter, npix, nside, nlheader, nmaps, ordering, pics, times
  real(dp)            :: nullval, h, c, k, T, temp, fre, x, z, p, test
  real(dp)            :: dust_T_init, dust_T_sigma, dummy
  real(dp)            :: missval = -1.6375d30
  logical(lgt)        :: anynull
  logical(lgt)        :: double_precision  

  character(len=128)              :: data, map1, map2, map3, map4, map5, mapHI, mask
  character(len=128)              :: arg1, arg2, arg3, fitsname, filename
  character(len=80)               :: line
  character(len=80), dimension(1) :: line2
  character(len=2)                :: number

  real(dp), allocatable, dimension(:,:)  :: HI_temp, HI_mask, masked_HI, T_map
  real(dp), allocatable, dimension(:,:)  :: masked_353, masked_545, masked_857, masked_240, masked_100
  real(dp), allocatable, dimension(:,:)  :: dummy_353, dummy_545, dummy_857, dummy_240, dummy_100
  real(dp), allocatable, dimension(:)    :: sum1, sum2, amps, y, freq, dummy_T, new_T
  character(len=80),  dimension(180)     :: header


  !------------------------------------------------------------------------------------------------------
  ! Daniel Herman 2019                                                                                  |
  !                                                                                                     |  
  ! This program will fit a temperature and an amplitude to the intensity of dust from the HI template. |
  ! We are only looking at high latitudes, where HI column density is < 4e20, because this is where     |
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

  data  = '../dust_data/'

  map1  = '../dust_data/npipe6v20_353-5_bmap_QUADCOR_n0064_60arcmin_MJy_no_cmb.fits'
  map2  = '../dust_data/npipe6v20_545-1_bmap_QUADCOR_n0064_60arcmin_MJy_no_cmb.fits'
  map3  = '../dust_data/npipe6v20_857-1_bmap_QUADCOR_n0064_60arcmin_MJy_no_cmb.fits'
  map4  = '../dust_data/DIRBE_240micron_1deg_h064_v2.fits'
  map5  = '../dust_data/DIRBE_100micron_Nside064_60a.fits'
  mapHI = '../dust_data/HI_vel_filter_60arcmin_0064.fits'
  mask  = '../dust_data/HI_mask.fits'

  i     = getsize_fits(mapHI,nside=nside,ordering=ordering,nmaps=nmaps)
  npix  = nside2npix(nside)

  allocate(masked_353(0:npix-1,nmaps),masked_545(0:npix-1,nmaps),masked_857(0:npix-1,nmaps))
  allocate(masked_240(0:npix-1,nmaps),masked_100(0:npix-1,nmaps),masked_HI(0:npix-1,nmaps))
  allocate(dummy_353(0:npix-1,nmaps),dummy_545(0:npix-1,nmaps),dummy_857(0:npix-1,nmaps))
  allocate(dummy_240(0:npix-1,nmaps),dummy_100(0:npix-1,nmaps))
  allocate(new_T(0:npix-1),HI_temp(0:npix-1,nmaps),HI_mask(0:npix-1,nmaps),T_map(0:npix-1,nmaps))
  allocate(sum1(5),sum2(5),amps(5),y(5),freq(5),dummy_T(0:npix-1))

  pics = 0

  niter = 1000

  freq(1) = 353.d0
  freq(2) = 545.d0
  freq(3) = 857.d0
  freq(4) = 1240.d0
  freq(5) = 2998.d0

  call read_bintab(mapHI, masked_HI, npix, nmaps, nullval, anynull, header=header)
  call read_bintab(mask, HI_mask, npix, nmaps, nullval, anynull, header=header)

  call read_bintab(map1, masked_353, npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map2, masked_545, npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map3, masked_857, npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map4, masked_240, npix, nmaps, nullval, anynull, header=header)
  call read_bintab(map5, masked_100, npix, nmaps, nullval, anynull, header=header)

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
        masked_353(i,:) = missval
        masked_545(i,:) = missval
        masked_857(i,:) = missval
        masked_240(i,:) = missval
        masked_100(i,:) = missval
        masked_HI(i,:)  = missval
        new_T(i)        = missval
      end if
    end do
  end do


  ! Here we calculate what the amplitude per map, and temperature per pixel should be for the best fit
  !
  ! with data = I_nu = A_nu * NHI * B_nu(T)
  ! Following the physical model:
  !             I_nu = kappa_nu * r * m_p * NHI * B_nu(T)
  ! So the amplitude we solve for is A_nu = kappa_nu * r * m_p which encapsulates the dust emissivity
  ! cross section per dust grain per NHI, in units of [cm^2]

  write(*,*) 'Initial amplitudes: '
  
  amps  = sample_A(new_T,npix)

  do n=1,5
    write(*,*) amps(n)
  end do
  write(*,*) '---------------'

  do m=1,times

     write(*,*) 'Iteration: ', m
     
     if (m .lt. 10) then
       write(number,10) m
     else if (m .gt. 9 .and. m .lt. 100) then
       write(number,11) m
     else if (m .gt. 99) then
       write(number,12) m
     endif

     filename   = 'amplitudes_' // trim(number) // '.dat'
     fitsname   = 'dust_Td_' // trim(number) // '.fits'

     dummy_T    = create_T(new_T,npix)
     amps       = sample_A(new_T,npix)

     new_T      = dummy_T

     do n=0,npix-1
        if (new_T(n) .gt. 10.d0 .and. new_T(n) .lt. 30.d0 ) then
           T_map(n,1) = new_T(n)
        else
           T_map(n,1) = missval
        end if
     end do

     open(35,file=trim(filename))

     do n=1,5
       write(*,*) amps(n)
       write(35, '(F20.8)') amps(n)
     end do

     write(*,*) '---------------'

     call write_bintab(T_map, npix, nmaps, header, nlheader, trim(fitsname))
     close(35)

   end do

   deallocate(masked_353,masked_545,masked_857,masked_100,masked_240,masked_HI)
   deallocate(dummy_353,dummy_545,dummy_857,dummy_240,dummy_100)
  

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


  function sample_A(T,npix)
     implicit none

     integer(i4b), intent(in)                     :: npix
     real(dp), dimension(0:npix-1), intent(in)    :: T
     real(dp), dimension(5)                       :: sample_A

     sum1 = 0.d0
     sum2 = 0.d0

     ! Calculated using chi^2 minimization for a single variable across all pixels:
     !
     ! chi^2 = sum((y_i - model_i)^2/sigma^2) =  sum( (data - amp*NHI*B_nu)^2/sigma)
     !       = sum( data^2 - 2*amp*NHI*B_nu*data + (amp*NHI*B_nu)^2)/sigma^2
     !
     ! Minimize:
     !
     ! d(chi^2)/d(amp) = sum(-2*NHI*B_nu*data + 2*amp*NHI^2*B_nu^2)/simga^2 = 0
     !
     ! amp = sum(NHI*B_nu*data / NHI^2*B_nu*2)
     !
     ! Summed over all pixels and all frequencies.
     
     do i=0,npix-1
        if (abs((T(i)-missval)/missval) < 1.d-8 .or. abs((masked_HI(i,1)-missval)/missval) < 1.d-8) then
           cycle
        else
           dummy_353(i,1)  = masked_HI(i,1)*planck(353.d0*1.d9,T(i))
           dummy_545(i,1)  = masked_HI(i,1)*planck(545.d0*1.d9,T(i))
           dummy_857(i,1)  = masked_HI(i,1)*planck(857.d0*1.d9,T(i))
           dummy_240(i,1)  = masked_HI(i,1)*planck(1240.d0*1.d9,T(i))
           dummy_100(i,1)  = masked_HI(i,1)*planck(2998.d0*1.d9,T(i))

           sum1(1) = sum1(1) + (masked_353(i,1)*dummy_353(i,1))
           sum1(2) = sum1(2) + (masked_545(i,1)*dummy_545(i,1))
           sum1(3) = sum1(3) + (masked_857(i,1)*dummy_857(i,1))
           sum1(4) = sum1(4) + (masked_240(i,1)*dummy_240(i,1))
           sum1(5) = sum1(5) + (masked_100(i,1)*dummy_100(i,1))

           sum2(1) = sum2(1) + (dummy_353(i,1)**2.d0)
           sum2(2) = sum2(2) + (dummy_545(i,1)**2.d0)
           sum2(3) = sum2(3) + (dummy_857(i,1)**2.d0)
           sum2(4) = sum2(4) + (dummy_240(i,1)**2.d0)
           sum2(5) = sum2(5) + (dummy_100(i,1)**2.d0)

        endif
      end do

      ! Amplitude in units of [cm^2]
      sample_A(1) = sum1(1)/sum2(1)
      sample_A(2) = sum1(2)/sum2(2)
      sample_A(3) = sum1(3)/sum2(3)
      sample_A(4) = sum1(4)/sum2(4)
      sample_A(5) = sum1(5)/sum2(5)

  end function sample_A


  function sample_T(x,y,z,T,mu,amp)
    implicit none

    real(dp), intent(in)                  :: x
    real(dp), dimension(5), intent(in)    :: y
    real(dp), intent(in)                  :: z
    real(dp), dimension(5), intent(in)    :: amp
    real(dp), intent(in)                  :: T
    real(dp), intent(in)                  :: mu
    real(dp), dimension(2)                :: a
    real(dp), dimension(5)                :: test
    real(dp)                              :: sample_T
    real(dp)                              :: temp
    real(dp)                              :: r
    real(dp)                              :: p
    real(dp)                              :: b
    real(dp)                              :: num

    temp = T
    a(1) = 1.d0

    do i=1,niter
      r    = rand_normal(temp,mu)
      test = 0.d0

      do j=1,5
        test(j) = amp(j)*z*planck(freq(j)*1.d9,r)
      end do

      b    = sum((abs(test-y)/y))
      a(2) = b/x
      p    = minval(a)
      
      call RANDOM_NUMBER(num)

      if (num > p) then
        if (r .gt. 10.d0 .and. r .lt. 30.d0 ) then
          temp = r
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
    real(dp), dimension(5)                       :: y
    real(dp)                                     :: z
    real(dp)                                     :: x

    do l=0,npix-1
      if (abs((T(l)-missval)/missval) < 1.d-8 .or. abs((masked_HI(l,1)-missval)/missval) < 1.d-8) then
        cycle
      else
        y(1) = masked_353(l,1)
        y(2) = masked_545(l,1)
        y(3) = masked_857(l,1)
        y(4) = masked_240(l,1)
        y(5) = masked_100(l,1)
        z    = masked_HI(l,1)
        x    = (abs(amps(1)*dummy_353(l,1) - y(1))/y(1) + &
                abs(amps(2)*dummy_545(l,1) - y(2))/y(2) + &
                abs(amps(3)*dummy_857(l,1) - y(3))/y(3) + &
                abs(amps(4)*dummy_240(l,1) - y(4))/y(4) + &
                abs(amps(5)*dummy_100(l,1) - y(5))/y(5))   

        create_T(l) = sample_T(x,y,z,T(l),dust_T_sigma,amps)
      endif
    end do
  end function create_T

end program dust_hi_fit
