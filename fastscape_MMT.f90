program fastscape_MMT

implicit none !this is good fortran practice

integer i,j,ij,ii,jj,ij_max_slope,ijk
integer nstep,nx,ny,nn
double precision xl,yl,dt,dx,dy,totaltime,k,m,n,u,drop,max_drop,length_max_slope
double precision slope,max_slope,length_drop

integer,dimension(:),allocatable:: rec,fixed,ndon
integer,dimension(:,:),allocatable:: donors
double precision,dimension(:),allocatable:: h,length

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
u = 2.e-3 !uplift rate, in units of...?

allocate (h(nn),rec(nn),length(nn),fixed(nn),ndon(nn))
allocate (donors(8,nn))

call random_number(h) !fortran random number generator will not return negative values

!here, give each node in the model which is not on the edge of the model space an elevation using the random number generator 
do j=2,ny-1
   do i=2,nx-1
      ij=i+(j-1)*nx !we are indexing this using only one number, ij. 
      h(ij)=0.+h(ij)
   enddo
enddo
!print*,'height',h !print the height to test if working

!initialize rec and length arrays
!fixed=1 !initialize fixed. we haven't done anything with this yet
do ij = 1,nn
   rec(ij) = ij 
enddo 
length=0

!Calculate flow directions and build receiver array
print*,"Calculating flow directions and building receiver array"
do j = 2,ny-1
   do i = 2,nx-1
      ij = i+(j-1)*nx 
      !initialize slope, length, and index corresponding to the maximum slope
      max_slope=0.d0 
      ij_max_slope=ij 
      length_max_slope=0 
      do jj=j-1,j+1
         do ii=i-1,i+1 !look at the eight nodes adjacent to ij
            drop=h(ij)-h(ii+(jj-1)*nx) !calculate the change in elevation with all adjacent nodes
            length_drop=((dx*(ii-i))**2 + (dy*(jj-j))**2)**0.5 !calculate the distance between all adjacent nodes
            slope=drop/length_drop !calculate the slope with all adjacent nodes
            if (slope.gt.max_slope)then 
!if the slope between node ij and an adjacent node is greater than the previously defined maximum slope, change slope, index, and length
               ij_max_slope=ii+(jj-1)*nx 
               max_slope=slope 
               length_max_slope=length_drop
            endif
         enddo
      enddo
      rec(ij) = ij_max_slope !loop through all adjacent nodes, define the receiver as the node with the maximum slope
      length(ij) = length_max_slope !define the length corresponding to the receiver array
   enddo
enddo
print*,"Receivers calculated!"

!initialize ndon and donors arrays
ndon=0 
donors=0

!Calculate donor array and number of donors for each cell
print*,"Calculating donor array"
do ij = 1,nn 
   ijk = rec(ij)
   ndon(ijk) = ndon(ijk) + 1 !adds 1 to ndon, the number of donors, for each receiver of a given donor
   donors(ndon(ijk),ijk) = ij !say that for node ijk, one of my donors is ij
enddo
print*,"Donors calculated!"

end
