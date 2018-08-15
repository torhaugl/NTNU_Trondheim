!  Author:           NTNU Trondheim, iGEM 18
!  Written by:       Haugland, Tor Strømsem
!  Acknowledgements: Pettersen, Jakob
!  Github:           torhaugl/NTNU_Trondheim
!
!  Program which simulates time-steps of a 3D grid of voxels to
!  calculate biofilm (EPS) and growth of bacteria/cell depending
!  on quorom sensing molecules (QSM). QSM activates bacteria so
!  that they produce much more EPS.
!
!
!  input --
!     Corresponds to an input file. Is used by many
!     subroutines to reduce the amount of arguments
!     passed. Using this is also faster than passing.
!
module input
   implicit none
   real,    parameter :: dt            = 0.1/60.0                 ! time step (min)
   real,    parameter :: final_time    = 14.0*60.0                ! minutes
   integer, parameter :: NumTrials     = floor(final_time / dt)   ! #steps to finish calculation
   real,    parameter :: c_bulk        = 0.2                      ! Concentration bulk substrate
   real,    parameter :: v_width       = 17.0                     ! Voxel length
   real,    parameter :: diff_s        = 40680.0
   real,    parameter :: diff_q        = 33300.0                  ! Diffusion substrate/QSM
   integer, parameter :: v_size(3)     = (/10,10,50/)             ! Sife of 3D grid of voxels in x,y,z direction TODO 10x10x100
   integer, parameter :: v_count       = v_size(1) * v_size(2) * v_size(3) ! #voxels
   real,    parameter :: Vmax          = 0.046                    ! Maximum substrate uptake
   real,    parameter :: Ks            = 0.00234                  ! Half-saturaton const (substrate uptake)
   real,    parameter :: Zqd           = 8.3
   real,    parameter :: Zqu           = 1230.0                   ! QSM production
   real,    parameter :: Kq            = 10.0                     ! Half-saturation QSM
   real,    parameter :: Ymax          = 0.444                    ! Bacteria yield percentage
   real,    parameter :: maintenance   = 0.006                    ! Bacteria eating
   real,    parameter :: avg_mass_cell = 420.0                    ! Average mass of bacteria
   real,    parameter :: density_cell  = 290.0                    ! Density of bacteria
   real,    parameter :: max_vol_frac  = 0.52                     ! Volume fraction that can be occupied in voxel
   real,    parameter :: Mmax          = 14700.0                  ! Maximum mass per particle
   integer, parameter :: pmax          = floor(density_cell * max_vol_frac * v_width**3 / Mmax) ! Max #particle per particle
   integer, parameter :: Nmax          = pmax                     ! Same as above, but used as length of arrays
   real,    parameter :: alpha         = 1.33                     ! Stochastic (down 2 up)
   real,    parameter :: beta          = 10.0                     ! Stochastic (up 2 down)
   real,    parameter :: gamma         = 0.1                      ! Stochastic (importance of QSM)
   real,    parameter :: eps_mass      = avg_mass_cell            ! Mass per EPS particle
   real,    parameter :: Zed           = 0.0
   real,    parameter :: Zeu           = 0.001                    ! Production of EPS, (d)own/(u)p
   real,    parameter :: mu            = 0.001                    ! Transfer coefficient out of voxel
   integer, parameter :: S_q           = 10
   integer, parameter :: S_s           = 10                       ! #Substeps for calculation of concentration
end
!
!  variable --
!     This module contains all variables that are changed
!     during the timestep. This makes it easier to visualize
!     which variables are available, and reduces arguments
!     passed into subroutines. Also is a speed-up.
!
!  TODO
!     Implement use variable in more subroutines
module variable
   use input
   implicit none
   integer, dimension(v_count,2)       :: eps_count
   real,    dimension(v_count,2)       :: eps_amount
   real,    dimension(v_count,2)       :: pressure
   real,    dimension(v_count,2+S_s)   :: c_s
   real,    dimension(v_count,2+S_q)   :: c_q
   real,    dimension(Nmax,v_count,2)  :: biomass
   integer, dimension(Nmax,v_count,2)  :: up
   real, dimension(9)                  :: timer ! Times total time of functions
   real :: curr_time
end
!
!  simulate --
!     This is the main program of the modeling project.
!     It calculates N time-steps and writes them to file
!     in the data/ folder.
program simulate
   use input      ! v_count, Nmax, diff_s, diff_q, v_size(3)
   use Variable   ! all variables
   implicit none

   ! Functions
   real, external :: avg

   ! Variables
   integer                    :: i, n, m ! Indexes for loop
   integer                    :: i_0, i_1 ! Indexes for bulk
   real                       :: r ! Random number
   real                       :: start, finish ! Total time taken
   real                       :: start_update, finish_update
   character(10)              :: time ! Output time started in string
   character(20)              :: filename

   ! Initialize arrays
   curr_time = 0.0
   c_s(:,:) = c_bulk
   c_q(:,:) = 0.0
   up(:,:,:) = 0
   biomass(:,:,:) = 0.0
   eps_count(:,:) = 0
   eps_amount(:,:) = 0.0
   pressure(:,:) = 0.0
   timer(:) = 0.0
   i_0 = v_count - v_size(1)*v_size(2)
   i_1 = v_count
   m = 0

   ! Insert particles into biomass: 1-2 bacteria, 400-800 mass, inactive
   do i = 1, v_size(1)*v_size(2)
      call random_number(r)
      call biomass_append(i, 400.0 + r*400.0, 0, biomass, up)
      biomass(:,i,1) = biomass(:,i,2)
   end do

   !!! CODE
   print*,"CODE START"
   call date_and_time(TIME=time)
   print*,"Time:",time(1:4)
   print*,"Number of steps:", NumTrials
   print*,

   call cpu_time(start)

! TODO Updates independent of neighbours (all except concentration)
!      can be calculated outside voxel loop
   do n = 1,NumTrials
      ! Print
      if (mod(n,floor(NumTrials/100.0)) == 0 .OR. n == 1) then
         print*, floor((real(n)/real(NumTrials)*100.0))
         write (filename,"(A5,I0.3,A4)") "data/", m, ".csv"
         filename = trim(filename)
         call write_all(filename)
         m = m + 1
      endif

      ! Cheap calculation
      call cpu_time(start_update)
      do i = 1, v_count
         call update_eps(i, biomass, up, eps_amount, eps_count)
         call update_mass(i)
         call update_division(i)
         call update_stochastics(i)
         call update_pressure(i, biomass, eps_count, pressure)
         call update_displacement(i, pressure, biomass, eps_count,up)
      enddo
      call cpu_time(finish_update)
      timer(1) = (finish_update - start_update) + timer(1)

      ! Concentration (Expensive calculation)
      call cpu_time(start_update)
      call update_concentration_s(biomass, diff_s, c_s) !Substrate
      call update_concentration_q(biomass,up, diff_q, c_q) !QSM
      call cpu_time(finish_update)
      timer(9) = (finish_update - start_update) + timer(9)

      ! Insert the newly calculated concentrations
      call cpu_time(start_update)
      c_s(:,1) = c_s(:,2)
      c_q(:,1) = c_q(:,2)
      biomass(:,:,1) = biomass(:,:,2)
      eps_count(:,1) = eps_count(:,2)
      eps_amount(:,1) = eps_amount(:,2)
      up(:,:,1) = up(:,:,2)
      pressure(:,1) = pressure(:,2)

      ! bulk
      c_s(i_0:i_1,:) = c_bulk
      c_q(i_0:i_1,:) = 0.0
      biomass(:,i_0:i_1,:) = 0.0
      eps_count(i_0:i_1,:) = 0.0
      eps_amount(i_0:i_1,:) = 0.0
      up(:,i_0:i_1,:) = 0.0
      pressure(i_0:i_1,:) = 0.0
      call cpu_time(finish_update)
      timer(7) = (finish_update - start_update) + timer(7)

      curr_time = curr_time + dt
   end do

   call cpu_time(finish)

   call date_and_time(TIME=time)
   print*,"Finish Time:",time(1:4)
   print*,"CPU time(s):", finish-start
   print*,"Model time(min):", final_time
   print*, timer
   print*,
   print*,"Top"   ,c_s(i_0-1,1)  ,c_q(i_0-1,1)  ,eps_amount(i_0-1,1)
   print*,"Mid"   ,c_s(i_1/2,1)   ,c_q(i_1/2,1)   ,eps_amount(i_1/2,1)
   print*,"Bot"   ,c_s(1,1)      ,c_q(1,1)      ,eps_amount(1,1)
   print*,"AVG"   ,avg(c_s(:,1)) ,avg(c_q(:,1)) ,avg(eps_amount(:,1))
   print*,"Bio"   ,biomass(:4,1,1)
   print*,"upp"   ,up(:4,1,1)
   call get_count_up(1,i)
   print*,"Up:"   ,i
   call get_count_down(1,i)
   print*,"Dow"   ,i
   call get_count_particles(biomass(:,1,1),eps_count(1,1),i)
   print*,"Par"   ,i


   ! Write to file
   !filename = "data/conc_count2.dat"
   !call write_count(filename)
   !filename = "data/conc_cq2.dat"
   !call write_concentration(c_q, filename)
   !filename = "data/conc_cs4.dat"
   !call write_concentration(c_s, filename)
   !print*, "Finished writing to *.dat files"


end program simulate

!!!!!!!!!!!!!!!
! SUBROUTINES !
!!!!!!!!!!!!!!!

!  write_all --
!     Outputs all important data to comma-separated values file (.csv)
!     Each row corresponds to a unique voxel
!     cs, cq, biomass, up, eps_count, eps_amount, pressure
!
!  Arguments:
!     filename    A new file to write to
!
subroutine write_all(filename)
   use input
   use variable
   use csv_file ! Module from FLIBS
      ! csv_write( lun, value, advance)
   implicit none

   character(len=20) :: filename
   integer :: pos(3), i
   logical :: exist

   open(1,file=filename,status="new",action="write")
   write(1,*) "time,x,y,z,cs,cq,biomass,up,eps_count,eps_amount"

   do i =1, v_count
      call index2xyz(i,pos)

      call csv_write(1,curr_time             ,.FALSE.)
      call csv_write(1,pos(1)                ,.FALSE.)
      call csv_write(1,pos(2)                ,.FALSE.)
      call csv_write(1,pos(3)                ,.FALSE.)
      call csv_write(1,c_s(i,1)              ,.FALSE.)
      call csv_write(1,c_q(i,1)              ,.FALSE.)
      call csv_write(1,sum(biomass(:,i,1))/avg_mass_cell   ,.FALSE.)
      call csv_write(1,sum(up(:,i,1))        ,.FALSE.)
      call csv_write(1,eps_count(i,1)        ,.FALSE.)
      call csv_write(1,eps_amount(i,1)       ,.TRUE.)
   enddo
   close(1)
end

!  write_count --
!     Outputs the particle count to a filetype
!     which is easy to 3D-scatter plot. Each particle
!     is placed at a random position in the voxel.
!
!  Arguments:
!     filename    A new file to write to
!
subroutine write_count(filename)
   use input
   use variable
   implicit none
   integer :: i, count, pos(3), j
   character(len=20) :: filename
   real r(3)

   open(unit=1, file=filename, status="new")
   do i = 1,v_count
      call index2xyz(i,pos)
      call get_count_particles(biomass(:,i,1), eps_count(i,1), count)

      do j = 1,count
         call random_number(r(1))
         call random_number(r(2))
         call random_number(r(3))
         write(1,*) pos(1)+r(1), pos(2)+r(2), pos(3)+r(3)
      enddo
   enddo

   close(1)
end

subroutine write_cellcount(biomass, filename)
   use input
   implicit none

   real, intent(in) :: biomass(Nmax, v_count, 2)
   character(len=20) :: filename
   integer :: sum=0, tsum=0, i,j, count

   open(unit=1,file=filename,status='new')
   do j=1,v_count
      tsum=0
      do i=1,Nmax
         call mass2cell_count(biomass(i,j,1), count)
         tsum = tsum + count
      enddo
      write(1,*) tsum
      sum = sum + tsum
   enddo
   write(1,*) "total", sum
end

subroutine write_concentration(c, filename)
   ! Easy scatter
   use input
   implicit none
   real, intent(in) :: c(v_count, 2)
   real c_count(v_count)
   character(len=20) :: filename
   integer :: i, j, count, pos(3)
   real r(3)

   c_count = 10 * c(:,1) !0.2 mol/L -> 20 particles in scatter

   open(unit=1, file=filename, status="new")
   do i = 1,v_count
      count = nint(c_count(i))
      call index2xyz(i,pos)

      do j = 1,count
         call random_number(r(1))
         call random_number(r(2))
         call random_number(r(3))
         write(1,*) pos(1)+r(1), pos(2)+r(2), pos(3)+r(3)
      enddo
   enddo

   close(1)
end

subroutine update_division(i)
   use input
   use variable !biomass
   implicit none
   integer, intent(in) :: i
   real :: randf, r
   real :: newmass
   integer :: n, j, count, newup

   do j = 1,Nmax
      if (biomass(j,i,1) > Mmax) then
         call random_number(randf)
         randf = 0.4 + 0.2*randf ! 0.4-0.6

         ! Mass
         newmass = randf * biomass(j,i,1)
         biomass(j,i,1) = (1-randf) * biomass(j,i,1)


         ! Up cells hypergeometric distribution
         newup = 0
         call mass2cell_count(biomass(j,i,1), count)
         do n=1,up(j,i,1)
            call random_number(r)
            if ( r < real(up(j,i,1)) / real(count)) then
               up(j,i,1) = up(j,i,1) - 1
               newup = newup + 1
            endif
            count = count - 1
         enddo
         ! ...

         ! call biomass_append(....)
         call biomass_append(i, newmass, newup, biomass, up)
      endif
   enddo
end

subroutine update_pressure(i, biomass, eps_count, pressure)
   ! If the number of particles is above max, pressure = 1
   use input
   implicit none

   integer, intent(in) :: i, eps_count(v_count,2)
   real, intent(in) :: biomass(Nmax,v_count,2)
   real, intent(out) :: pressure(v_count,2)
   integer :: count

   call get_count_particles(biomass(:,i,1), eps_count(i,1), count)

   if (count >= pmax) then
      pressure(i,2) = real(pmax)
   else
      pressure(i,2) = real(count) / (real(pmax) - real(count))
   end if

end

subroutine update_displacement(i, pressure, biomass, eps_count, up)
   ! Depends on neighbours
   use input
   implicit none

   integer, intent(in) :: i
   real, intent(in) :: pressure(v_count,2)
   integer, intent(out) :: eps_count(v_count,2), up(Nmax,v_count,2)
   real, intent(out) :: biomass(Nmax,v_count,2)

   integer :: up_temp,j,k, j_list_neigh(6), count, count_neigh, chosen, rand_int, count_displaced
   real :: mass,tot_pressure, rand, P(6) ! Cumulative probability
   logical :: eps_displaced
   integer :: pos(3),pos2(3)


   call get_count_particles(biomass(:,i,1), eps_count(i,1), count)
   if (mu*pressure(i,1)*count < 1 ) return ! Save lots of time

   ! Calculate how many particles are displaced

   count_displaced = 0
   call get_index_neighbours(i, j_list_neigh)

   if (count < Nmax) then !Update transfer depending on neighbours
      do j =1,6
         k = j_list_neigh(j)
         if (k < 1) cycle
         call get_count_particles(biomass(:,k,1), eps_count(k,1), count_neigh)
         if (pressure(i,1) > pressure(k,1) .AND. count > count_neigh) then
            count_displaced = count_displaced + floor(mu * (pressure(i,1) - pressure(k,1)) *(count - count_neigh) )
         endif
      end do
   else if (count >= Nmax) then ! Irrelevant of neighbours
      count_displaced = count - Nmax
   end if

   if (count_displaced < 0) print*, "Error, particle_transfer<0", count_displaced


   ! If paritcles are displaced
   if (count_displaced > 0) then
      ! Calculate probability for each neighbour
      tot_pressure = 0
      do j =1,6
         k = j_list_neigh(j)
         if (k < 1) cycle
         if (pressure(i,1) > pressure(k,1)) then
            tot_pressure = tot_pressure + (pressure(i,1) - pressure(k,1))
         end if
      end do


      do j =1,6
         k = j_list_neigh(j)
         if (pressure(i,1) > pressure(k,1) .AND. k > 0) then
            P(j) = (pressure(i,1) - pressure(k,1)) / tot_pressure
         else
            P(j) = 0
         end if

         if (j > 1) then ! Cumulative sum
            P(j) = P(j-1) + P(j)
         end if
      end do

      if (all(P == 0.0)) then
         count_displaced = 0
         !print*, "Voxel overloaded", i
         return
      endif

      ! Calculate which neighbour(index) and particle is chosen
      ! neighbour
      do while (count_displaced > 0)
         call random_number(rand)
         do j =1,6
            if ( rand < P(j) ) then
               chosen = j_list_neigh(j)
               exit
            endif
         enddo
         call index2xyz(i,pos)
         call index2xyz(chosen,pos2)
         !print*, "Displace (pos):", pos, "to", pos2

         ! Choose particle type
         call random_number(rand)
         rand_int = 1 + floor(count * rand) ! 1 - count
         eps_displaced = (rand_int <= eps_count(i,1))

         ! Remove and Append particles to neighbour
         if (eps_displaced) then
            eps_count(i,2) = eps_count(i,2) - 1
            eps_count(chosen,2) = eps_count(chosen,2) + 1
         else !biomass displaced
            call biomass_remove_random(i,biomass,mass,up,up_temp)
            call biomass_append(chosen,mass,up_temp,biomass,up)
         endif
         count_displaced = count_displaced - 1
      enddo
   endif
end

real function avg(arr)
   ! Calculates the avg value of an array of v_count values
   use input
   implicit none
   real, intent(in) :: arr(v_count)
   avg = sum(arr) / size(arr)
   return
end

subroutine update_eps(i, biomass, up, eps_amount, eps_count)
   ! Calculate new eps. Create particle if eps > eps_mass
   use input !v_count Nmax Zed Zeu eps_mass
   implicit none

   integer, intent(in) :: i
   real, dimension(Nmax, v_count,2), intent(in)  :: biomass
   integer, dimension(Nmax,v_count,2),intent(in) :: up
   real, dimension(v_count,2), intent(out) :: eps_amount
   integer, intent(out) :: eps_count(v_count,2)
   integer :: count_up, count_down

   call get_count_up(i,count_up)
   call get_count_down(i,count_down)

   eps_amount(i,2) = eps_amount(i,2) + Zed*count_down + Zeu*count_up

   if (eps_amount(i,2) >= eps_mass) then
      eps_amount(i,2) = eps_amount(i,2) - eps_mass
      eps_count(i,2) = eps_count(i,2) + 1
   end if
end

subroutine update_mass(i)
   use input !v_count Nmax Ymax, Ks, Vmax, m, dt
   use variable !c_s, biomass
   implicit none

   integer, intent(in) :: i

   biomass(:,i,2) = biomass(:,i,2) + dt*Ymax* ( Vmax*c_s(i,2)/(Ks + c_s(i,2) ) - maintenance)*biomass(:,i,2)
end

subroutine get_count_up(i,count_up)
   ! TODO Worth it?
   use input
   use variable
   implicit none

   integer, intent(in) :: i
   integer, intent(out) :: count_up
   integer :: count,j
   real :: M

   count_up = sum(up(:,i,1))
   M = sum(biomass(:,i,1))
   call mass2cell_count(M,count)
   if (count_up > count) count_up = count

end

subroutine get_count_down(i,count_down)
   ! TODO Worth it?
   use input
   use variable
   implicit none

   integer, intent(in) :: i
   integer, intent(out) :: count_down
   integer :: count,j
   real :: M=0 !Tot mass in voxel

   M = sum(biomass(:,i,1))
   call mass2cell_count(M,count)
   count_down = count - sum(up(:,i,1))
   if (count_down < 0) count_down = 0

end

subroutine get_count_particles(biomass_particle,eps_count_particle, count)
   ! TODO Worth it?
   ! Particle count in voxel!
   ! Nonzero entries in biomass_particle are a particle
   use input
   implicit none

   integer, intent(in)  :: eps_count_particle
   real,    intent(in)  :: biomass_particle(Nmax)
   integer, intent(out) :: count

   call get_count_nonzero(biomass_particle, count)
   count = count + eps_count_particle
end

subroutine get_count_nonzero(arr, n)
   ! Counts
   use input !Nmax
   implicit none

   real,    intent(in)  :: arr(Nmax)
   integer, intent(out) :: n
   integer :: j

   n = count(arr .NE. 0)
end

subroutine update_stochastics(i)
   !TODO error due to up(j,i,1), j is not correct
   use input !v_count, dt
   use variable !c_q, biomass, up
   implicit none

   integer, intent(in)  :: i
   real :: d2u, u2d, rand
   integer :: count_up, count_down,count, count_d2u, count_u2d,n,j

   do j=1,Nmax
      count_up = up(j,i,1)
      call mass2cell_count(biomass(j,i,1), count)
      count_down = count - count_up

      call probability_down2up(i,d2u)
      call probability_up2down(i,u2d)

      ! Down -> Up
      count_d2u = 0
      do n=1,count_down
         call random_number(rand)
         if (rand < dt * d2u) then
            count_d2u = count_d2u + 1
         end if
      end do
      ! Up -> Down
      count_u2d = 0
      do n=1,count_up
         call random_number(rand)
         if (rand < dt * u2d) then
            count_u2d = count_u2d + 1
         end if
      end do

      up(j,i,2) = up(j,i,2) + (count_d2u-count_u2d)
   enddo

end

subroutine probability_down2up(i, prob)
   ! Calculates probability down to up
   use input !v_count, alpha, gamma
   use variable !c_q
   implicit none
   integer, intent(in) :: i
   real, intent(out) :: prob

   prob = alpha*c_q(i,1) / (1 + gamma *c_q(i,1))
end

subroutine probability_up2down(i, prob)
   ! Calculates probability down to up
   use input !v_count, beta, gamma
   use variable
   implicit none
   integer, intent(in) :: i
   real, intent(out) :: prob

   prob = beta / (1 + gamma *c_q(i,1))
end


subroutine update_concentration_s(biomass, diff, c)
   ! Updates conentration (Forward Euler step / IMEX)
   ! Input index, diffusion constant, voxel width, c before/after
   ! Output new concentration
   use input
   implicit none
   integer i
   real, dimension(Nmax,v_count,2) :: biomass
   real, dimension(v_count,2+S_s), intent(out) :: c
   real, intent(in) :: diff
   integer :: j_list_neigh(6), j, k, s
   real, dimension(v_count) :: M, prod, prod_der

   M = sum(biomass(:,:,1), dim=1)

   do s = 1, S_s ! Substeps
      c(:,2+s) = c(:,1+s)
      prod(:) = -Vmax * M(:) * c(:,1+s) / (Ks + c(:,1+s)) !Update producion_s
      prod_der(:) = -Vmax*M(:) * Ks / ((Ks + c(:,1+s))**2)

      do i = 1, v_count
         call get_index_neighbours(i, j_list_neigh)
         do j = 1,6 !For every neighbour
            k = j_list_neigh(j)
            if ( k > 0 ) then
               if (s == 1) then
                  c(i,3) = c(i,3) + 1/(1-dt/S_s * prod_der(i)/(v_width**3)) &
                  *dt/S_s * (c(k,1) - c(i,1)) *  diff/(v_width*v_width) + dt/S_s*prod(i)/(v_width**3)
               else
                  c(i,2+s) = c(i,2+s) + 1/(1-dt/S_s * prod_der(i)/(v_width**3)) &
                  *dt/S_s*((c(k,1+s) - c(i,1+s)) *  diff/(v_width*v_width) + prod(i)/(v_width**3))
               endif
            end if
         end do
         if(c(i,2+s) < 0.0) c(i,2+s) = 0.0
      enddo
   enddo
   c(:,2) = c(:,2+S_s)
end

subroutine update_concentration_q(biomass, up, diff, c)
   ! Updates conentration (Forward Euler step)
   ! Input index, diffusion constant, voxel width, c before/after
   ! Output new concentration
   use input
   implicit none
   real,    dimension(Nmax,v_count,2), intent(in) :: biomass
   integer, dimension(Nmax,v_count,2), intent(in) :: up
   real,    intent(in) :: diff
   real,    dimension(v_count,2+S_q), intent(out) :: c
   real,    dimension(v_count) :: prod
   integer :: j_list_neigh(6), i, j, k, s
   integer :: count_up, count_down
   !TODO Move count function before s-loop

   do s = 1, S_q ! Substeps
      c(:,2+s) = c(:,1+s)

      do i = 1,v_count
         call get_count_up(i,count_up)
         call get_count_down(i,count_down)
         if (Kq == 0) then
            prod(i) = Zqd*count_down + Zqu*count_up
         else
            prod(i) = Zqd*count_down + Zqu*count_up * c(i,1+s)/(Kq + c(i,1+s) )
         endif

         call get_index_neighbours(i, j_list_neigh)
         do j = 1,6 !For every neighbour
            k = j_list_neigh(j)
            if ( k > 0 ) then
               if (s == 1) then
                  c(i,3) = c(i,3) + dt/S_q * (c(k,1) - c(i,1)) *  diff/(v_width*v_width) + dt/S_q*prod(i)/(v_width**3)
               else
                  c(i,2+s) = c(i,2+s) + dt/S_q*(c(k,1+s) - c(i,1+s)) *  diff/(v_width*v_width) + dt/S_q*prod(i)/(v_width**3)
               endif
            endif
         enddo
         if(c(i,2+s) < 0.0) c(i,2+s) = 0.0
      enddo
   enddo
   c(:,2) = c(:,2+S_q)
end

subroutine biomass_remove_random(i, biomass, mass, up, up_temp)
   use input
   implicit none

   integer, intent(in) :: i
   real, intent(out) :: biomass(Nmax, v_count, 2)
   integer, intent(out) :: up(Nmax,v_count,2)
   real, intent(out) :: mass
   integer, intent(out) :: up_temp
   real :: rand
   integer :: rand_int, count_nonzero, j, n

   call get_count_nonzero(biomass(:,i,2), count_nonzero)
   call random_number(rand)
   rand_int = 1 + floor(count_nonzero * rand)

   ! Only check nonzeo values. 'n' variable keeps count.
   n = 1
   do j=1,Nmax
      if (biomass(j,i,2) == 0.0) cycle

      if (rand_int == n) then
         mass = biomass(j,i,2)
         biomass(j,i,2) = 0.0

         up_temp = up(j,i,2)
         up(j,i,2) = 0
         exit
      else
         n = n + 1
      endif
   enddo

end

subroutine biomass_append(i, mass,up_temp,biomass,up)
   ! TODO check total particle count
   ! Append mass at first non-zero biomass
   use input !Nmax
   implicit none


   integer, intent(in) :: i
   real, intent(in) :: mass
   real, intent(out) :: biomass(Nmax,v_count,2)
   integer, intent(out) :: up(Nmax,v_count,2)
   integer :: j
   integer :: up_temp

   do j = 1,Nmax
      if (biomass(j,i,2) == 0.0) then
         biomass(j,i,2) = mass
         up(j,i,2) = up_temp
         exit
      endif

      !if (j == Nmax) print*, "Couldn't append!", biomass(:,i,2)
   enddo

end

subroutine index2xyz(i, pos)
   ! Input index i
   ! Returns x y z in pos list
   use input
   implicit none
   integer, intent(in) :: i
   integer :: pos(3)
   integer :: x, y, z

   x = mod(i-1, v_size(1))
   y = mod( (i-1)/v_size(1), v_size(2))
   z = (i-1)/(v_size(1) * v_size(2))

   pos = (/x, y, z/)
end

subroutine xyz2index(pos, i)
   ! Input xyz in pos list
   ! Returns integer
   use input
   implicit none
   integer, intent(in) :: pos(3)
   integer, intent(out) :: i
   integer :: x,y,z
   x = pos(1)
   y = pos(2)
   z = pos(3)

   i = x + y * v_size(1) + z * v_size(1) * v_size(2) + 1

   if ((i < 1) .OR. (i > v_size(1)*v_size(2)*v_size(3))) then
      print*, "Error, i<=0 or i > v_count",i,x,y,z
      i = 0
   end if
end

subroutine get_index_neighbours(i, i_list_neigh)
   ! Input index
   ! Return list of indexes to neighbours
   use input
   implicit none
   integer, intent(in):: i
   integer, intent(out) :: i_list_neigh(6)
   integer :: pos(3), npos(3)
   i_list_neigh = (/0,0,0,0,0,0/)

   call index2xyz(i,pos)

   ! x +- 1 (Boundary condition)
   npos = pos + (/1,0,0/)
   call continuous_boundary_condition(npos(1), v_size(1))
   call xyz2index(npos, i_list_neigh(1))

   npos = pos - (/1,0,0/)
   call continuous_boundary_condition(npos(1), v_size(1))
   call xyz2index(npos, i_list_neigh(2))

   ! y +- 1
   npos = pos + (/0,1,0/)
   call continuous_boundary_condition(npos(2), v_size(2))
   call xyz2index(npos, i_list_neigh(3))

   npos = pos - (/0,1,0/)
   call continuous_boundary_condition(npos(2), v_size(2))
   call xyz2index(npos, i_list_neigh(4))

   ! z +- 1
   npos = pos + (/0,0,1/)
   if (.NOT. (npos(3) >= v_size(3))) then
      call xyz2index(npos, i_list_neigh(5))
   end if

   npos = pos - (/0,0,1/)
   if (.NOT. (npos(3) < 0)) then
      call xyz2index(npos, i_list_neigh(6))
   end if
end

subroutine continuous_boundary_condition(pos_q, v_size_q)
   ! Input position in one direction with voxel size,
   ! Returns new position according to CBC
   implicit none
   integer, intent(out)::pos_q
   integer, intent(in)::v_size_q

   if (pos_q < 0) then
      pos_q = v_size_q-1
   else if (pos_q >= v_size_q) then
      pos_q = 0
   end if
end



subroutine mass2cell_count(mass, count)
   use input
   implicit none

   real, intent(in):: mass
   integer, intent(out) :: count

   count = ceiling(mass/avg_mass_cell)
end
