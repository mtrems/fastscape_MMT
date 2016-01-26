program fastscape_MMT

implicit none !this is good fortran practice

integer i,j,ij,ii,jj,ij_max_slope
integer nstep,nx,ny,nn
double precision xl,yl,dt,dx,dy,totaltime,k,m,n,u,drop,max_drop,length_max_slope
double precision slope,max_slope,length_drop

integer,dimension(:),allocatable:: rec,fixed
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

allocate (h(nn),rec(nn),length(nn),fixed(nn)) !allocate some variables

call random_number(h) !call the random number generator in fortran. this random number generator will not return negative values

!here, give each node in the model which is not on the edge of the model space an elevation using the random number generator 
do j=2,ny-1
   do i=2,nx-1
      ij=i+(j-1)*nx !we are indexing this using only one number, ij. starts in the lower left corner
      h(ij)=0.+h(ij)
   enddo
enddo

fixed=1 !initialize fixed. we haven't done anything with this yet
rec=0 !initialize receiver array with zeros. these zeros will be replaced with the index of the node which receives the flow of water from the current cell
length=0 !initialize length array with zeros. these zeros will be replaced with the length between the current node and its receiver node.

do j = 2,ny-1
   do i = 2,nx-1
      ij = i+(j-1)*nx
      max_slope=0.d0
      ij_max_slope=0
      length_max_slope=0
      do jj=j-1,j+1
         do ii=i-1,i+1
            drop=h(ij)-h(ii+(jj-1)*nx)
            length_drop=((dx*(ii-i))**2 + (dy*(jj-j))**2)**0.5
            slope=drop/length_drop
            if (slope.gt.max_slope)then
               ij_max_slope=ii+(jj-1)*nx
               max_slope=slope
               length_max_slope=length_drop
            endif
         enddo
      enddo
      rec(ij) = ij_max_slope
      length(ij) = length_max_slope
   enddo
enddo 

end
