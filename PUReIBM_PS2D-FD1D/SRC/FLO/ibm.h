!    PUReIBM-PS2D-FD1D is a three-dimensional psudeo-spectral particle-resolved
!    direct numerical simulation solver for detailed analysis of homogeneous
!    fixed and freely evolving fluid-particle suspensions. PUReIBM-PS2D-FD1D
!    is a continuum Navier-Stokes and scalar solvers based on Cartesian grid that utilizes
!    Immeresed Boundary method to represent particle surfuces. The details about the solvers
!    can be found in the below papers in SUBRAMANIAM's group. 
!    Copyright (C) 2015, Shankar Subramaniam, Rahul Garg, Sudheer Tenneti, Bo Sun, Mohammad Mehrabadi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    For acknowledgement, please refer to the following publications:
!     For hydrodynamic solver :
!     (1) TENNETI, S. and SUBRAMANIAM, S., 2014, Particle-resolved direct numerical
!         simulation for gas-solid flow model development. Annu. Rev. Fluid Mech.
!         46 (1) 199-230.
!     (2) M. Mehrabadi, S. Tenneti, R. Garg, and S. Subramaniam, 2015, Pseudo-turbulent 
!         gas-phase velocity fluctuations in homogeneous gas-solid flow: fixed particle
!         assemblies and freely evolving suspensions. J. Fluid Mech. 770 210-246.
!
!     For scalar solver :
!     (3) S. Tenneti, B. Sun, R. Garg, S. Subramaniam, 2013, Role of fluid heating in dense
!         gas-solid flow as revealed by particle-resolved direct numerical simulation.
!         International Journal of Heat and Mass Transfer 58 471-479.

#define SQR(a) ((a)*(a))

#if PARALLEL

#define PARALLEL_START() CALL MPI_INIT(err_code)
#define PARALLEL_FINISH()CALL MPI_FINALIZE(err_code)  
#define GET_NPROCS(comm,psize) CALL MPI_COMM_SIZE(comm,psize,err_code)
#define GET_PROCESSOR_RANK(comm,rank) CALL MPI_COMM_RANK(comm, rank, err_code)
#define BROADCAST_INT(item,nitems,from,comm)   CALL MPI_BCAST(item,nitems,MPI_INTEGER,from,comm,err_code)
#define BROADCAST_DOUBLE(item,nitems,from,comm)   CALL MPI_BCAST(item,nitems,MPI_DOUBLE_PRECISION,from,comm,err_code)
#define BROADCAST_REAL(item,nitems,from,comm)   CALL MPI_BCAST(item,nitems,MPI_DOUBLE_REAL,from,comm,err_code)
#define BROADCAST_CHARARR(str,strlen,from,comm)   CALL MPI_BCAST(str,strlen,MPI_CHARACTER,from,comm,err_code)
#define BROADCAST_LOGICAL(item,nitems,from,comm)   CALL MPI_BCAST(item,nitems,MPI_LOGICAL,from,comm,err_code)
#define BROADCAST_STRING(str,strlen,from,comm)				\
  if(I_AM_NODE_ZERO) strlen = LEN_TRIM(str); BROADCAST_INT(strlen,1,from,comm); BROADCAST_CHARARR(str,strlen,from,comm)
#define SEND_STRING(str,strlen,to,sitag,sctag,comm) strlen = LEN_TRIM(str); CALL MPI_SEND(strlen,1,MPI_INTEGER,to,sitag,comm,err_code); CALL MPI_SEND(str,strlen,MPI_CHARACTER,to,sctag,comm,err_code)
#define RECV_STRING(str,strlen,from,ritag,rctag,comm,status) CALL MPI_RECV(strlen,1,MPI_INTEGER,from,ritag,comm,status,err_code); CALL MPI_RECV(str,strlen,MPI_CHARACTER,from,rctag,comm,status,err_code)

#define CREATE_CART_TOPOLOGY(commold,ndims,psize,isperiodic,reorder,newcomm) CALL MPI_CART_CREATE(commold,ndims,psize,isperiodic,reorder,newcomm,err_code)
#define GET_SHIFT_PROCS(comm,dir,shift,src,dest) CALL MPI_CART_SHIFT(comm,dir,shift,src,dest,err_code)
!#define DOMAIN_1D_DECOMP(nitems,psize,rank,start,end) CALL MPE_DECOMP1D(nitems,psize,rank,start,end)
#define GLOBAL_INT_SUM(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_INTEGER,MPI_SUM,comm,err_code)
#define GLOBAL_LONGINT_SUM(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_INTEGER8,MPI_SUM,comm,err_code)
#define GLOBAL_DOUBLE_SUM(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_DOUBLE_PRECISION,MPI_SUM,comm,err_code)
#define GLOBAL_COMPLEX_SUM(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_DOUBLE_COMPLEX,MPI_SUM,comm,err_code)
#define GLOBAL_DOUBLE_MAX(sendbuf,recvbuf,nitems,comm) CALL MPI_ALLREDUCE(sendbuf,recvbuf,nitems,MPI_DOUBLE_PRECISION,MPI_MAX,comm,err_code)
#define I_AM_IN_THIS_PROC(ii) ((xstart-1.le.ii).and.(ii.le.xend+1)) 
#define LOCAL_INDEX(ii) ii = (ii-xstart+1)
#define GLOBAL_INDEX(ii)(ii+xstart-1)

#define WEST_PERIODIC_POINT_IMAGE(ii,iimage) if(((ii).lt.1))iimage = mxf+(ii)-1
#define EAST_PERIODIC_POINT_IMAGE(ii,iimage) if(((ii).ge.mxf))iimage = (ii)-(mxf-1)

#define WEST_PERIODIC_IMAGE(ii,iimage,x,ximage) if(((ii).lt.1))then ;iimage = mxf+(ii)-1;ximage = mxf+(x)-1;endif
#define EAST_PERIODIC_IMAGE(ii,iimage,x,ximage) if(((ii).ge.mxf))then ;iimage = (ii)-(mxf-1);ximage = (x)-(mxf-1);endif

#define WEST_PERIODIC_IMAGE_MOD(ii,iimage,x,ximage) if(((ii).le.1))then ;iimage = mxf+(ii)-1;ximage = mxf+(x)-1;endif
#define EAST_PERIODIC_IMAGE_MOD(ii,iimage,x,ximage) if(((ii).ge.mxf-1))then ;iimage = (ii)-(mxf-1);ximage = (x)-(mxf-1);endif

#define WEST_PERIODIC_IMAGE_PRES(ii,iimage,x,ximage) if(((ii).lt.0))then ;iimage = mxf+(ii)-1;ximage = mxf+(x)-1;endif
#define EAST_PERIODIC_IMAGE_PRES(ii,iimage,x,ximage) if(((ii).ge.mxf-1))then ;iimage = (ii)-(mxf-1);ximage = (x)-(mxf-1);endif

#define RSENDRECV(sitem,nsitems,to,stag,ritem,nritems,from,rtag,comm,status) CALL MPI_SENDRECV(sitem,nsitems,MPI_DOUBLE_PRECISION,to,stag,ritem,nritems,MPI_DOUBLE_PRECISION,from,rtag,comm,status,err_code)
#define CSENDRECV(sitem,nsitems,to,stag,ritem,nritems,from,rtag,comm,status) CALL MPI_SENDRECV(sitem,nsitems,MPI_DOUBLE_COMPLEX,to,stag,ritem,nritems,MPI_DOUBLE_COMPLEX,from,rtag,comm,status,err_code)

#define POINT_IN_VEL_GRID(x) ((xstart.le.x).and.(x.lt.xend+1))
#define CELL_IN_VEL_GRID(x) ((xstart.le.x).and.(x.le.xend))

#define POINT_IN_PRESS_GRID(x) (((xstart-half).le.x).and.(x.lt.(xend+half)))
#define CELL_IN_PRESS_GRID(x) (((xstart-1).le.x).and.(x.lt.(xend)))

#define POINT_IN_PROC(x) ((xstart.le.x).and.(x.le.xend))
#define CELL_IN_PROC(ii) ((xstart-1.le.ii).and.(ii.le.xend))

#define RPR_POINT_IN_PROC(x) ((xstart-2.le.x).and.(x.lt.xend+2))
#define RPR_CELL_IN_PROC(x) ((xstart-2.le.x).and.(x.lt.xend+2))

#define EAST_NO_MANS_LAND(x) (x).eq.(xend)
#define WEST_NO_MANS_LAND(x) (x).eq.(xstart-1)

#define CONCAVE(x,n,m) (x(n).gt.xc(m,n))

#define CREATE_2D_LSLICE(count,blocklength,stride,newtype) CALL MPI_TYPE_VECTOR(count,blocklength,stride,MPI_LOGICAL,newtype,err_code)

#define CREATE_2D_RSLICE(count,blocklength,stride,newtype) CALL MPI_TYPE_VECTOR(count,blocklength,stride,MPI_DOUBLE_PRECISION,newtype,err_code)
#define CREATE_2D_CSLICE(count,blocklength,stride,newtype) CALL MPI_TYPE_VECTOR(count,blocklength,stride,MPI_DOUBLE_COMPLEX,newtype,err_code)
#define COMMIT(dtype) CALL MPI_TYPE_COMMIT(dtype,err_code)
#define VECSENDRECV(sitem,nsitems,vectype,to,stag,ritem,nritems,from,rtag,comm,status) CALL MPI_SENDRECV(sitem,nsitems,vectype,to,stag,ritem,nritems,vectype,from,rtag,comm,status,err_code)
#define VECSENDRECV2(sitem,nsitems,vectype1,to,stag,ritem,nritems,vectype2,from,rtag,comm,status) CALL MPI_SENDRECV(sitem,nsitems,vectype1,to,stag,ritem,nritems,vectype2,from,rtag,comm,status,err_code)
#define BARRIER(comm) CALL MPI_BARRIER(comm,err_code)
!#define FILL_WEST_SEND_RBUF4(buff,arr) buff => arr(1,1:my,1:mz)
#else

#define PARALLEL_START() 
#define PARALLEL_FINISH() 
#define GET_NPROCS(comm,psize) nproc = 1
#define GET_PROCESSOR_RANK(comm,rank) rank = 0
#define BROADCAST_INT(item,nitems,from,comm)
#define BROADCAST_DOUBLE(item,nitems,from,comm)
#define BROADCAST_REAL(item,nitems,from,comm)
#define BROADCAST_CHARARR(str,strlen,from,comm)
#define BROADCAST_STRING(str,strlen,from,comm)
#define SEND_STRING(str,strlen,to,sitag,sctag,comm)		
#define RECV_STRING(str,strlen,from,ritag,rctag,comm,status)		
#define BROADCAST_LOGICAL(item,nitems,from,comm)
#define CREATE_CART_TOPOLOGY(commold,ndims,psize,isperiodic,reorder,newcomm)
#define GET_SHIFT_PROCS(comm,dir,shift,src,dest)
#define DOMAIN_1D_DECOMPOSE(nitems,psize,rank,start,end)
!#define DOMAIN_1D_DECOMP(nitems,psize,rank,start,end)
#define GLOBAL_INT_SUM(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_LONGINT_SUM(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_DOUBLE_MAX(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_DOUBLE_SUM(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define GLOBAL_COMPLEX_SUM(sendbuf,recvbuf,nitems,comm) recvbuf = sendbuf
#define LOCAL_INDEX(ii)
#define GLOBAL_INDEX(ii)  (ii+xstart-1) 
#define I_AM_IN_THIS_NODE(ii) .TRUE.
#define CELL_IN_PROC(ii) ((xstart.le.ii).and.(ii.le.xend))
#define RSENDRECV(sitem,nsitems,to,stag,ritem,nritems,from,rtag,comm,status)
#define CSENDRECV(sitem,nsitems,to,stag,ritem,nritems,from,rtag,comm,status)
#define VECSENDRECV(sitem,nsitems,vectype,to,stag,ritem,nritems,from,rtag,comm,status)
#define VECSENDRECV2(sitem,nsitems,vectype1,to,stag,ritem,nritems,vectype2,from,rtag,comm,status)
#define CREATE_2D_RSLICE(count,blocklength,stride,newtype)
#define CREATE_2D_CSLICE(count,blocklength,stride,newtype)
#define CREATE_2D_LSLICE(count,blocklength,stride,newtype)
#define BARRIER(comm)
#define COMMIT(newtype)
#endif

#define I_AM_NODE_ZERO (myid.eq.0)
#define I_AM_LAST_PROC (xend.eq.mx1)
#define CELL_IN_BOX(ii) ((1.le.ii).and.(ii.le.mx1))
#define WEST_PERIODIC_IMAGE(ii,iimage,x,ximage) if(((ii).lt.1))then ;iimage = mxf+(ii)-1;ximage = mxf+(x)-1;endif
#define EAST_PERIODIC_IMAGE(ii,iimage,x,ximage) if(((ii).ge.mxf))then ;iimage = (ii)-(mxf-1);ximage = (x)-(mxf-1);endif
#define EAST_PERIODIC_IMAGE_MOD(ii,iimage,x,ximage) if(((ii).ge.mxf-1))then ;iimage = (ii)-(mxf-1);ximage = (x)-(mxf-1);endif
