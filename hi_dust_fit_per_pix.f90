program dust_hi_fit
  use healpix_types
  use pix_tools
  use fitstools
  use head_fits
  implicit none

  integer(i4b)        :: i, j, l, m, n, niter, npix, nside, nlheader, nmaps, ordering, pix, times, bands
  real(dp)            :: nullval, h, c, k, T, temp, fre, x, z, p, test, T_sum
  real(dp)            :: dust_T_init, dust_T_sigma, chisq, thresh
  real(dp)            :: missval = -1.6375d30
  logical(lgt)        :: anynull
  logical(lgt)        :: double_precision  

  character(len=128)              :: mapHI, output, mapfile, rmsfile, freqfile, bandfile
  character(len=128)              :: arg1, arg2, arg3, arg4, arg5, arg6
  character(len=128)              :: chi_file, amp_dist, fitsname, filename, file1, file2
  character(len=80)               :: line, data, version 
  character(len=80), dimension(1) :: line2
  character(len=4)                :: number

  real(dp), allocatable, dimension(:)          :: model_T, new_T, y, freq, norm, s
  real(dp), allocatable, dimension(:,:)        :: HI, T_map,  sum1, sum2, amps, clamps
  real(dp), allocatable, dimension(:,:,:)      :: maps, model, rmss, cov
  character(len=80),  dimension(180)           :: header
  character(len=100),allocatable, dimension(:) :: freqs, map, rms


  !------------------------------------------------------------------------------------------------------
  ! Daniel Herman 2019                                                                                  |
  !                                                                                                     |  
  ! This program will fit a temperature and an amplitude to the intensity of dust from the HI template. |
  ! We are only looking at high latitudes, where HI column density is < 4e20, because this is where the |
  ! dust/gas ratio is appoximately linear.                                                              |
  !                                                                                                     |  
  !-----------------------------------------------------------------------------------------------------|  

10 format (I1)
11 format (I2)
12 format (I3)
13 format (I4)
  
  if (iargc() < 6) then
     write(*,*) 'Usage:'
     write(*,*) '   dust_hi_fit [temp estimate (full sky)] [std of temperature (for fitting)]'
     write(*,*) '               [# of iterations] [HI threshold] [# of bands] [version # (v10)]'
     write(*,*) ''
     stop
  end if

  call getarg(1,arg1)
  read(arg1,*) dust_T_init
  call getarg(2,arg2)
  read(arg2,*) dust_T_sigma
  call getarg(3,arg3)
  read(arg3,*) times
  call getarg(4,arg4)
  read(arg4,*) thresh
  call getarg(5,arg5)
  read(arg5,*) bands
  call getarg(6,arg6)
  read(arg6,*) version

  data  = '../dust_data/sed_data/'
  output= 'results/' // trim(version)// '/'

  write(*,*) bands

  call system('mkdir -p ./' // trim(output))

  mapfile  = 'sed_maps_v3.txt'
  rmsfile  = 'sed_rms_v3.txt'
  freqfile = 'sed_freqs_v3.txt'
  bandfile = 'sed_bands_v3.txt'
  
  mapHI = trim(data) // 'HI_vel_filter_60arcmin_0064.fits'

  i     = getsize_fits(mapHI,nside=nside,ordering=ordering,nmaps=nmaps)
  npix  = nside2npix(nside)

  allocate(maps(0:npix-1,nmaps,bands),model(0:npix-1,nmaps,bands),HI(0:npix-1,nmaps),model_T(0:npix-1))
  allocate(new_T(0:npix-1),T_map(0:npix-1,nmaps),map(bands),rms(bands))
  allocate(rmss(0:npix-1,nmaps,bands), cov(0:npix-1,nmaps,bands),norm(bands),s(bands))
  allocate(y(bands),freq(bands),sum1(bands),sum2(bands),amps(bands),clamps(bands),freqs(bands))

  ! Open and read frequencies, map, and rms map names

  open(27,file=trim(bandfile))
  open(28,file=trim(rmsfile))
  open(29,file=trim(freqfile))
  open(30,file=trim(mapfile))
  read(29,*) freq

  do i=1,bands
     read(28,fmt='(a)') rms(i)
     read(27,fmt='(a)') freqs(i)
     read(30,fmt='(a)') map(i)
     map(i) = trim(data) // trim(map(i))
     rms(i) = trim(data) // trim(rms(i))
     freqs(i) = trim(freqs(i))
  end do

  close(27)
  close(28)
  close(29)
  close(30)

  ! Read maps
  call read_bintab(mapHI, HI, npix, nmaps, nullval, anynull, header=header)

  do i=1,bands
     call read_bintab(map(i), maps(:,:,i), npix, nmaps, nullval, anynull, header=header)
     call read_bintab(rms(i), rmss(:,:,i), npix, nmaps, nullval, anynull, header=header)
  end do
     

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

  niter = 5000
  pix  = 0

  ! Mask maps
  do j = 1, nmaps
    do i = 0, npix-1
      if (HI(i,j) .gt. thresh) then
        maps(i,:,:) = missval
        rmss(i,:,:) = missval
        new_T(i)    = missval
        HI(i,:)     = missval
      else
        pix        = pix + 1
      end if
    end do
  end do

  !-----------------------------------------------------------------------------------------------------|  
  ! Here we calculate what the amplitude per map, and temperature per pixel should be for the best fit
  ! with data = I_nu = A_nu * NHI * B_nu(T)
  !
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
    if (abs((HI(i,1)-missval)/missval) < 1.d-8) then
       cycle
    else
       do j=1,bands
!          cov(i,1,j)   = rmss(i,1,j)**2.d0
!          model(i,1,j) = HI(i,1)*planck(freq(j)*1.d9,new_T(i))
!          sum1(j)      = sum1(j) + (maps(i,1,j)*cov(i,1,j)*model(i,1,j))
!          sum2(j)      = sum2(j) + (model(i,1,j)**2.d0)*cov(i,1,j)
!          norm(j)      = norm(j) + cov(i,1,j)*model(i,1,j)**2.d0
       end do
    endif
  end do

  do j=1,bands
     amps(j) = sum1(j)/sum2(j)
     norm(j) = sqrt(pix*norm(j))
  end do

  clamps = amps

  do n=1,bands
    write(*,*) freqs(n) // ': ', clamps(n)
  end do

  write(*,*) '----------------------------------'
  write(*,*) ''

  filename      = trim(output) // 'amplitudes.dat'
  chi_file      = trim(output) // 'chi_sq.dat'
  amp_dist      = trim(output) // 'amplitude_distribution.dat'

  open(35, file = trim(filename))
  open(36, file = trim(chi_file))

  !------------------------------------------------------------
  !------------------------------------------------------------

  do m=1,times

     write(*,*) 'Iteration', m
     write(*,*) '-----------------------'

     ! Here is where all of the calculations happen
     ! --------------------------------------------------

     model_T = create_T(new_T,npix)
     new_T   = model_T
     amps    = sample_A(new_T,npix)

     clamps  = amps

     chisq   = compute_chisq(new_T,clamps)

     T_sum = 0.d0

     do n=0,npix-1
        if (new_T(n) .gt. 10.d0 .and. new_T(n) .lt. 30.d0 ) then
           T_map(n,1) = new_T(n)
           T_sum      = T_sum + T_map(n,1)
        else
           T_map(n,1) = missval
        end if
     end do

     write(*,*) 'Total Chi-square: '
     write(*,*) chisq
     write(*,*) ''

     write(*,*) 'Mean dust temperature: '
     write(*,*) T_sum/pix
     write(*,*) ''

     write(*,*) 'Amplitudes: '
     do n=1,bands
        write(*,*) freqs(n) // ': ', clamps(n)
     end do

     ! --------------------------------------------------

     ! Write result maps
     if ( mod(m,100) .EQ. 0) then
        if (m .lt. 10) then
           write(number,10) m
        else if (m .gt. 9 .and. m .lt. 100) then
           write(number,11) m
        else if (m .gt. 99) then
           write(number,12) m
        else if (m .gt. 999) then
           write(number, 13) m
        endif

        fitsname   = trim(output) // 'dust_Td_' // trim(number) // '.fits'
        write(*,*) '----------------------------------'
        call write_maps(npix,nmaps,header)
     end if
    
     write(35,*) clamps(:)
     write(36,'(1(E17.8))') chisq
     write(*,*) '----------------------------------'
     write(*,*) ''

  end do

!  open(37, file = trim(amp_dist))
!  do n=1,1000
!     do j=1,bands
!        s(j) = clamps(j) + rand_normal(0.d0,1.d0)/norm(j)
!     end do
!     write(37,'(6(E17.8))') s(1), s(2), s(3), s(4), s(5), s(6)
!  end do

  deallocate(maps,HI,norm)
  deallocate(model,rmss,cov)
  close(35)
  close(36)
  close(37)

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

    integer(i4b), intent(in)                  :: npix
    real(dp), dimension(0:npix-1), intent(in) :: T
    real(dp), dimension(bands)                :: sample_A, tau, r
    real(dp), dimension(2)                    :: a
    real(dp)                                  :: chisq2, p, num
    integer(i4b)                              :: b

    a(1) = 1.d0

    sum1 = 0.d0
    sum2 = 0.d0

    do i=0,npix-1
       if (abs((HI(i,1)-missval)/missval) < 1.d-8) then
          cycle
       else
          do j=1,bands
             model(i,1,j) = HI(i,1)*planck(freq(j)*1.d9,T(i))
             sum1(j)      = sum1(j) + (maps(i,1,j)*model(i,1,j)*cov(i,1,j))
             sum2(j)      = sum2(j) + (model(i,1,j)**2.d0*cov(i,1,j))
          end do
       endif
    end do

    tau = sum1/sum2

    sample_A = tau

  end function sample_A

  function create_T(T,npix)
    implicit none

    integer(i4b), intent(in)                     :: npix
    real(dp), dimension(0:npix-1), intent(in)    :: T
    real(dp), dimension(0:npix-1)                :: create_T
    ! real(dp), dimension(niter+1)                 :: chi, tump, accept, prob
    real(dp), dimension(bands)                   :: y,covs,test
    real(dp), dimension(2)                       :: a
    real(dp)                                     :: x,temp,r,b,c,num!,naccept

    do i=0,npix-1
       if (abs((HI(i,1)-missval)/missval) < 1.d-8) then
          cycle
       else
          x = 0.d0
          do j=1,bands
             y(j) = maps(i,1,j)
!             covs(j) = cov(i,1,j)
             x = x + (clamps(j)*model(i,1,j) - y(j))**2.d0
          end do

          temp = T(i)
          a(1) = 1.d0
          c    = x

          ! prob(1)   = 1.d0
          ! chi(1)    = x
          ! tump(1)   = temp
          ! accept(1) = 0.d0
          ! naccept   = 0

          do l=1,niter
             r    = rand_normal(temp,dust_T_sigma)
             test = 0.d0
             do j=1,bands
                test(j) = clamps(j)*HI(i,1)*planck(freq(j)*1.d9,r)
             end do
             b    = sum((test-y)**2.d0)
             a(2) = exp(-b+c)
             p    = minval(a)

             call RANDOM_NUMBER(num)

             if (num < p) then
                if (r .lt. 35 .and. r .gt. 10) then
                   temp  = r
                   c     = b
                   ! naccept = naccept + 1
                end if
             end if

             ! if (i == 600) then
             !  prob(l+1)   = p 
             !  chi(l+1)    = c
             !  tump(l+1)   = temp
             !  accept(l+1) = naccept/l
             ! end if
          end do

          ! if (i == 600) then

          !   open(40,file = trim(output) // 'prob.dat')
          !   open(41,file = trim(output) // 'chi.dat')
          !   open(42,file = trim(output) // 'temps.dat')
          !   open(43,file = trim(output) // 'accept.dat')

          !   do l = 1,niter+1
          !     write(40,*) prob(l)
          !     write(41,*) chi(l)
          !     write(42,*) tump(l)
          !     write(43,*) accept(l)
          !   end do

          ! end if
          create_T(i)  = temp
       endif
    end do
  end function create_T

  function compute_chisq(T,amp)
    implicit none
    real(dp), dimension(0:npix-1), intent(in) :: T
    real(dp), dimension(bands),    intent(in) :: amp
    real(dp)                                  :: chi
    real(dp)                                  :: compute_chisq

    chi = 0.d0
    do j=1,bands
       do i=0,npix-1
          if (abs((HI(i,1)-missval)/missval) < 1.d-8) then
             cycle
          else
             chi = chi + (maps(i,1,j) - amp(j)*model(i,1,j))**2.d0/cov(i,1,j)
          end if
       end do
    end do

    compute_chisq = (chi/(pix*bands))

  end function compute_chisq

  subroutine write_maps(npix,nmaps,header)
    ! Outputs the temperature map, the modeled map, and the residual map
    implicit none
    integer(i4b), intent(in)                        :: npix, nmaps
    integer(i4b)                                    :: y, b
    character(len=80),  dimension(180), intent(in)  :: header
    character(len=80)                               :: file1, file2, file3
    real(dp), dimension(0:npix-1,nmaps,bands)       :: modl, resid

    write(*,*) "Writing maps!"

    file1 = trim(output) // 'dust_Td_' // trim(number) // '.fits'

    do y=0,npix-1
       modl(y,1,:)  = amps(:)*model(y,1,:)
       resid(y,1,:) = maps(y,1,:) - modl(y,1,:)
    end do

    call write_bintab(T_map, npix, nmaps, header, nlheader, file1)
    do y=1,bands
       file2 = trim(output) // 'model_'// trim(freqs(y)) // '_' // trim(number) // '.fits'
       file3 = trim(output) // 'resid_'// trim(freqs(y)) // '_' // trim(number) // '.fits'
       call write_bintab(modl(:,:,y), npix, nmaps, header, nlheader, file2)
       call write_bintab(resid(:,:,y), npix, nmaps, header, nlheader, file3)
    end do

  end subroutine write_maps

end program dust_hi_fit
