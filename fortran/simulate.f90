module input
implicit none

! Input file / Parameters
real, parameter :: dt = 0.07/60.0 !time step (min)
real, parameter :: final_time = 1.0*60.0 !minutes
integer, parameter :: NumTrials = final_time / dt
real, parameter :: c_bulk = 0.2 ! Concentration bulk substrate
real, parameter :: v_width = 17.0 !Voxel length
real, parameter :: diff_s = 40680.0, diff_q = 33300.0 !Diffusion substrate/QSM
integer, parameter :: v_size(3) = (/10,10,100/)
integer, parameter :: v_count = v_size(1) * v_size(2) * v_size(3)
end

program simulate
! Simulates 1 time-step
use input
implicit none

! Variables
!Concentration substrate (before and after timestep)
real, dimension(v_count, 2) :: c_s, c_q
real :: prod_s, prod_q
integer :: pos(3) !xyz position&


! Loops
integer :: i, n ! Standard indexes for loop
integer :: i_list_neigh(6) !Index for neighbours of a voxel
real :: start, finish




! Initialize array
c_s(:,:) = c_bulk
call xyz2index((/5,5,50/), i)
c_s(i,:) = 0.8*v_count

c_q(:,:) = 0.0

!!! CODE
print*,"CODE START"
print*,"Number of steps:", NumTrials

call cpu_time(start)
do n = 1, NumTrials

   do i = 1, v_count
      call update_concentration(i, prod_s, diff_s, c_s) !Substrate
      call update_concentration(i, prod_q, diff_q, c_q) !QSM
   end do
   c_s(:,1) = c_s(:,2) ! Insert the newly calculated concentrations
   c_q(:,1) = c_q(:,2)
end do
call cpu_time(finish)

print*,"CPU time(s):", finish-start
print*,"Model time(min):", final_time

print*,c_s(20,1)
print*,c_s(10,1)
print*,c_s(1,1)
print*,sum(c_s(:,1))

end program simulate










subroutine update_concentration(i, prod, diff, c)
   ! Updates conentration
   ! Input index, diffusion constant, vortex width, c before/after
   ! Output new concentration
   use input
   implicit none
   integer, intent(in) :: i
   real, dimension(v_count,2), intent(out) :: c
   real, intent(in) :: prod, diff
   integer :: j_list_neigh(6), j
   
   call get_index_neighbours(i, j_list_neigh)
   do j = 1,6 !For every neighbour
      if (.NOT. (j_list_neigh(j)  < 1 ) ) then
         c(i,2) = c(i,2) + dt * (c(j_list_neigh(j),1) - c(i,1)) *  diff/(v_width*v_width)
      end if
   end do

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
   ! Input position in one direction with vortex size,
   ! Returns new position according to CBC
   integer, intent(out)::pos_q
   integer, intent(in)::v_size_q

   if (pos_q < 0) then
      pos_q = v_size_q-1
   else if (pos_q >= v_size_q) then
      pos_q = 0
   end if
end
