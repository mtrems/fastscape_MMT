program fastscape_MMT

implicit none !this is good fortran practice

integer i,j,ij,ii,jj,ij_max_slope,ijk,istep,nfreq,nprint
integer nstep,nx,ny,nn
integer nstack,ncolor
double precision xl,yl,dt,dx,dy,totaltime,k,m,n,u,drop,max_drop,length_max_slope
double precision slope,max_slope,length_drop,c

integer,dimension(:),allocatable:: rec,ndon,stack,color,stack2id,longestriver
integer,dimension(:,:),allocatable:: donors
double precision,dimension(:),allocatable:: h,x,y,length,area

character cs*4

nstep = 1000 !number of time steps
nx = 101 !number of nodes in the x-direction
ny = 101 !number of nodes in the y-direction
dt = 1000. !years

nn = nx*ny !total number of nodes is the product of nx and ny
totaltime = nstep*dt !calculate the total time to be modeled

xl = 100.e3 !x-dimension of the model, in meters
yl = 100.e3 !y-dimension of the model, in meters
dx = xl/(nx-1) !grid size in x-direciton, in meters
dy = yl/(ny-1) !grid size in y-direction, in meters

k = 1.e-4 !erodibility/susceptibility constant, in units of...?
n = 1. !slope exponent, a constant
m = 0.4 !area exponent, a constant
u = 2.e-3 !uplift rate, in units of meters per year

nfreq = 20 !frequency with which to output results to a text file
nprint = 0 !counter for printing

allocate (h(nn),x(nn),y(nn),rec(nn),length(nn),ndon(nn))
allocate (donors(8,nn))
allocate (stack(nn),color(nn),stack2id(nn))
allocate (area(nn))

!initialize variables
do i = 1,nn
   h(i) = 0
   x(i) = 0
   y(i) = 0
enddo

!define x and y positions
do j = 1,ny
   do i = 1,nx
      ij = i+(j-1)*nx
      x(ij) = dx*i - dx
      y(ij) = dy*j - dy
   enddo
enddo

!here, give each node in the model which is not on the edge of the model space an elevation using the random number generator 
do j = 2,ny-1
   do i = 2,nx-1
      ij = i+(j-1)*nx !we are indexing this using only one number, ij. 
      call random_number(h(ij)) 
      !fortran random number generator will not return negative values
   enddo
enddo

do istep = 1,nstep

!initialize rec and length arrays
do ij = 1,nn
   rec(ij) = ij 
enddo 
length=0

!Calculate flow directions and build receiver array
do j = 2,ny-1
   do i = 2,nx-1
      ij = i+(j-1)*nx 
      !initialize slope, length, and index corresponding to the maximum slope
      max_slope=0.d0 
      ij_max_slope=ij 
      length_max_slope=0 
      do jj=j-1,j+1
         do ii=i-1,i+1 !look at the eight nodes adjacent to ij
            !calculate the change in elevation with all adjacent nodes
            drop=h(ij)-h(ii+(jj-1)*nx) 
            !calculate the distance between all adjacent nodes
            length_drop=((dx*(ii-i))**2 + (dy*(jj-j))**2)**0.5             
            slope=drop/length_drop !calculate the slope with all adjacent nodes
            if (slope.gt.max_slope)then 
!if the slope between node ij and an adjacent node is greater than the previously defined maximum slope, change slope, index, and length
               ij_max_slope=ii+(jj-1)*nx 
               max_slope=slope 
               length_max_slope=length_drop
            endif
         enddo
      enddo
      !loop through adjacent nodes, define the receiver as the node with the maximum slope
      rec(ij) = ij_max_slope 
      length(ij) = length_max_slope !define the length corresponding to the receiver array
   enddo
enddo

!initialize ndon and donors arrays
ndon=0 
donors=0

!Calculate donor array and number of donors for each cell
do ij = 1,nn
   if (rec(ij).ne.ij) then 
      ijk = rec(ij)
      !adds 1 to ndon, the number of donors, for each receiver of a given donor
      ndon(ijk) = ndon(ijk) + 1 
      donors(ndon(ijk),ijk) = ij !say that for node ijk, one of my donors is ij
   endif
enddo

!Build the stack, a list of all nodes in an order that makes computing discharge easy
nstack = 0
ncolor = 0

do ij = 1,nn
   stack(ij) = 0
   color(ij) = 0
   stack2id(ij) = 0 
enddo

do ij = 1,nn
   if(rec(ij).eq.ij) then
     nstack = nstack + 1
     ncolor = ij
     stack(nstack) = ij
     color(nstack) = ncolor
     stack2id(ij) = nstack
     !call recursive function
     call add_to_stack (ij,nstack,stack,ndon,donors,nn,ncolor,color,stack2id)
   endif
enddo

!set drainage area of each node to the size of a single node initially
do ij = 1,nn
   area(ij) = dx*dy
enddo

!calculate upstream drainage area using the reverse stack
do ij = nn,1,-1
   ijk = stack(ij)
   if (rec(ijk).ne.ijk) then
      area(rec(ijk)) = area(rec(ijk)) + area(ijk)
   endif
enddo

!evolve the landscape
do j = 2,ny-1
   do i = 2,nx-1
      ijk = i+(j-1)*nx
      !do the uplift component to first. Do not apply to nodes on the boundary.
      h(ijk) = h(ijk) + (u*dt) 
      enddo
enddo
!then do the erosion component. Only apply this to nodes who are not their own receivers (i.e., not boundaries and not sinks)
do i = 1,nn
   ijk = stack(i)
   if (rec(ijk).ne.ijk) then
         c = k*(area(ijk)**m)*dt*(1./length(ijk))
         h(ijk) = (h(ijk) + (c*h(rec(ijk))))*(1./(1.+c))
   endif
enddo

!for writing results to file. Gives each time step saved a unique file name
write(cs,'(i4)') istep
if (istep.lt.10) cs(1:3)='000'
if (istep.lt.100) cs(1:2)='00'
if (istep.lt.1000) cs(1:1)='0'

!write x, y, height, color ID, drainage area, and time in the model to file
if (mod(istep,nfreq).eq.0) then
   open(unit=8,status="new",file="out"//cs//".txt")
   do i = 1,nn
      write(8,*) x(i),y(i),h(i),color(i),area(i),istep*dt
   enddo
endif
close(8)

enddo

!extract the longest river profile from the landscape
open (7,file='profile.txt',status='unknown')
i = maxloc(area,1)
do while (ndon(i).ne.0)
   i = donors(maxloc(area(donors(1:ndon(i),i)),1),i)
   write(7,*) length(i),h(i),area(i),(h(i)-h(rec(i)))/length(i)
enddo
close(7)

end program fastscape_MMT

!Recursive subroutine that will go to progressively upstream nodes and build the stack
recursive subroutine add_to_stack (ij,nstack,stack,ndon,donors,nn,ncolor,color,stack2id)

implicit none

integer nn,stack(nn),ndon(nn),donors(8,nn),color(nn),stack2id(nn)
integer ij,j,nstack,ncolor



do j = 1,ndon(ij)
      nstack = nstack + 1
      stack(nstack) = donors(j,ij)
      color(nstack) = ncolor
      stack2id(ij) = nstack
      call add_to_stack(donors(j,ij),nstack,stack,ndon,donors,nn,ncolor,color,stack2id)
enddo

return

end subroutine add_to_stack
