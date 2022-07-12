module GridCellConnectionSetType
  !This module is for creating the 2D grid connections for lateral GW flow.
  !Han Qiu 2021.12
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use shr_infnan_mod  , only : isnan => shr_infnan_isnan,nan => shr_infnan_nan, assignment(=)
  use decompMod      , only : bounds_type
  implicit none
  save
  public

  type, public :: gridcell_connection_set_type
     Integer, pointer  :: nconn         => null() ! number of connections
     Integer, pointer  :: grid_id_up(:) => null() ! list of ids of upwind cells
     Integer, pointer  :: grid_id_dn(:) => null() ! list of ids of downwind cells
     Real(r8), pointer :: dist(:)       => null() ! list of distance vectors
     Real(r8), pointer :: face_length(:)=> null() ! list of edge of faces normal to distance vectors
     Real(r8), pointer :: uparea(:)     => null() ! list of up cell areas of horizaontal faces 
     Real(r8), pointer :: downarea(:)   => null() ! list of down cell areas of horizaontal faces 
     Real(r8), pointer :: dzg(:)        => null() ! list of areas of dz between downwind and upwind cells
     Real(r8), pointer :: facecos(:)    => null() ! dot product of the cell face normal vector and cell centroid vector
     Real(r8), pointer :: vertcos(:)    => null() ! dot product of the cell face normal vector and cell centroid vector for vertical flux, the rank for vertcos
     ! is from 1 to column size which is different from rank of lateral faces
   contains
     procedure, public :: Init => col_connect_init
  end type gridcell_connection_set_type

  public :: get_natveg_column_id 

  type (gridcell_connection_set_type), public, target :: conn   ! connection type

contains

  ! ************************************************************************** !
  subroutine col_connect_init(this, bounds)

    implicit none

    type(bounds_type), intent(in)       :: bounds
    class(gridcell_connection_set_type) :: this 
    Integer                             :: n, iconn,begc,endc,begg, endg,ii,jj 
    Integer                             :: g,nx,ny
    Real(r8)                            :: x(11,11),y(11,11),zh(11,11),dx(10,9),dy(9,10),dz(10,10),slopex(10,9),slopey(9,10),slopexx(10,10),slopeyy(10,10)

    begc = bounds%begc;  endc = bounds%endc
    begg = bounds%begg;  endg = bounds%endg

    nx = 10
    ny = 10
    iconn=0
    n = (nx-1)*ny+(ny-1)*nx

    ! number of connections in each layer
    allocate(this%nconn)          ;  this%nconn = n                   
    allocate(this%grid_id_up(n))  ;  this%grid_id_up(:) = 0
    allocate(this%grid_id_dn(n))  ;  this%grid_id_dn(:) = 0
    allocate(this%face_length(n)) ;  this%face_length(:) = 0
    allocate(this%uparea(n))      ;  this%uparea(:) = 0
    allocate(this%downarea(n))    ;  this%downarea(:) = 0
    allocate(this%dist(n))        ;  this%dist(:) = 0
    allocate(this%dzg(n))         ;  this%dzg(:) = 0
    allocate(this%facecos(n))     ;  this%facecos(:) = 0
    allocate(this%vertcos(nx*ny)) ;  this%vertcos(:) = 0

    dz = 1.0_r8
    do ii = 1, nx+1
       do jj = 1,ny+1
          x (jj , ii) = ii*10._r8 - 10._r8
          y (jj , 11-ii+1) = (0.2_r8*(ii-1)+(1._r8 - 0.02_r8*(ii-1._r8)*2._r8)*(jj-1._r8))*10._r8 ! converge
       end do
    end do


    do ii = 1, nx+1
       do jj = 1,ny+1
          zh(ii,jj) = 10*SIN(x(ii,jj)*DACOS(-1._r8)/100._r8-DACOS(-1._r8)/2.0_r8)+20._r8  !no y slope
       end do
    end do

    do ii = 1, ny           
       do jj = 1,nx-1
          slopex(ii,jj) = (zh(ii,jj+2)+zh(ii+1,jj+2)-zh(ii,jj)-zh(ii+1,jj))/(x(ii,jj+2)+x(ii+1,jj+2)-x(ii,jj)-x(ii+1,jj))
          dx(ii,jj) = abs((x(ii,jj+2)+x(ii+1,jj+2)-x(ii,jj)-x(ii+1,jj))/4._r8)
       end do
    end do

    if(ny > 1) then
       do ii = 1, ny-1
          do jj = 1,nx
             slopey(ii,jj) = (zh(ii+2,jj)+zh(ii+2,jj+1)-zh(ii,jj)-zh(ii,jj+1))/(y(ii+2,jj)+y(ii+2,jj+1)-y(ii,jj)-y(ii,jj+1))
             dy(ii,jj) = abs((y(ii+2,jj)+y(ii+2,jj+1)-y(ii,jj)-y(ii,jj+1))/4._r8)
          end do
       end do
    endif

    slopey=0._r8


    if(ny>1) then
       do ii= 1, ny
          do jj = 1, nx-1
             ! assume grid is counted along nx, cross-sectional area is dy*dz, top and bottom area is , dy-bar*dx
             iconn                  = iconn+1
             this%grid_id_up(iconn) = (ii-1)*nx+jj                                                            !  Step-2: Eventually will need to read from surface dataset
             this%grid_id_dn(iconn) = (ii-1)*nx+jj+1                                                          !  There is already some code that we will be able to
             this%face_length(iconn)= abs((y(ii+1,jj+1)-y(ii,jj+1))*dz(ii,jj))                                ! cross-sectional area, rectangular
             this%uparea(iconn)     = abs((y(ii+1,jj)-y(ii,jj)+y(ii+1,jj+1)-y(ii,jj+1))*dx(ii,jj))/2.0_r8     ! up cell vertical area, trapzoid
             this%downarea(iconn)   = abs((y(ii+1,jj+1)-y(ii,jj+1)+y(ii+1,jj+2)-y(ii,jj+2))*dx(ii,jj))/2.0_r8 ! down cell vertical area
             this%dist(iconn)       = sqrt((slopex(ii,jj)*dx(ii,jj))**2._r8+dx(ii,jj)**2._r8)                 ! triangle law
             this%dzg(iconn)        = slopex(ii,jj)*dx(ii,jj)                                                 ! down cell elevation - up cell elevation, positive go up hill, negative go down hill
             this%facecos(iconn)    = 1/sqrt(1 + slopex(ii,jj)**2._r8)
          enddo
       enddo

       do ii= 1, ny-1
          do jj = 1, nx
             ! assume grid is counted along ny, cross-sectional area is dx*dz
             iconn                  = iconn+1
             this%grid_id_up(iconn) = jj+ nx*(ii-1)                                                                         !  Step-2: Eventually will need to read from surface dataset
             this%grid_id_dn(iconn) = jj+ nx*ii                                                                             !  There is already some code that we will be able to
             this%face_length(iconn)= abs((x(ii+1,jj+1)-x(ii+1,jj))*dz(ii,jj))                                              !  cross-sectional area
             this%uparea(iconn)     = abs((y(ii+1,jj)-y(ii,jj)+y(ii+1,jj+1)-y(ii,jj+1))*(x(ii+1,jj+1)-x(ii+1,jj)))/2._r8;     ! up cell vertical area
             this%downarea(iconn)   = abs((y(ii+2,jj)-y(ii+1,jj)+y(ii+2,jj+1)-y(ii+1,jj+1))*(x(ii+1,jj+1)-x(ii+1,jj)))/2._r8; ! down cell vertical area
             this%dist(iconn)       = sqrt((slopey(ii,jj)*dy(ii,jj))**2._r8 + dy(ii, jj)**2._r8)                              ! triangle law
             this%dzg(iconn)        = slopey(ii,jj)*dy(ii,jj)
             this%facecos(iconn)    = 1._r8/sqrt(1 + slopey(ii,jj)**2._r8)
          enddo
       enddo
    else
       do ii= 1, ny
          do jj = 1, nx-1
             ! assume grid is counted along nx
             iconn                  = iconn+1
             this%grid_id_up(iconn) = jj                                                                      !  Step-2: Eventually will need to read from surface dataset
             this%grid_id_dn(iconn) = jj+1                                                                    !  There is already some code that we will be able to
             this%face_length(iconn)= abs(dx(jj,1)*dz(jj,1))                                                  !  use to fill this data structure
             this%uparea(iconn)     = abs((y(ii+1,jj)-y(ii,jj)+y(ii+1,jj+1)-y(ii,jj+1))*dx(ii,jj))/2.0_r8     ! up cell vertical area, trapzoid
             this%downarea(iconn)   = abs((y(ii+1,jj+1)-y(ii,jj+1)+y(ii+1,jj+2)-y(ii,jj+2))*dx(ii,jj))/2.0_r8 ! down cell vertical area
             this%facecos(iconn)    = 1._r8/sqrt(1._r8 + slopex(ii,jj)**2_r8)
             this%dist(iconn)       = sqrt((slopex(jj,1)*dx(jj,1))**2._r8+dx(ii,jj)**2_r8)                    ! triangle law !may not used for now
             this%dzg(iconn)        = slopex(ii,jj)*dx(ii,jj) 
          enddo
       enddo
    endif
    iconn = 0
    if (ny > 1) then
       do ii = 1, ny
          do jj = 1, nx
             iconn = iconn+1
             slopexx(ii,jj) = (zh(ii,jj+1)+zh(ii+1,jj+1)-zh(ii,jj)-zh(ii+1,jj))/(x(ii,jj+1)+x(ii+1,jj+1)-x(ii,jj)-x(ii+1,jj))
             slopeyy(ii,jj) = (zh(ii+1,jj)+zh(ii+1,jj+1)-zh(ii,jj)-zh(ii,jj+1))/(y(ii+1,jj)+y(ii+1,jj+1)-y(ii,jj)-y(ii,jj+1))
             this%vertcos(iconn) = 1._r8/sqrt(1._r8 + slopeyy(ii,jj)**2._r8) * 1._r8/sqrt(1._r8 + slopexx(ii,jj)**2._r8)
          enddo
       enddo
    else
       do jj = 1, nx
          iconn = iconn+1
          slopexx(1,jj) = (zh(1,jj+1)+zh(2,jj+1)-zh(1,jj)-zh(2,jj))/(x(1,jj+1)+x(2,jj+1)-x(1,jj)-x(2,jj))
          this%vertcos(iconn) = 1._r8/sqrt(1 + slopexx(1,jj)**2._r8)
       enddo
    endif
  end subroutine col_connect_init

  function get_natveg_column_id(id, bounds) result(id_out)

    implicit none

    type(bounds_type), intent(in) :: bounds
    integer, intent(in)           :: id
    integer                       :: id_out,begc,endc,begg,endg

    begc = bounds%begc;  endc = bounds%endc
    begg = bounds%begg;  endg = bounds%endg


    id_out=begg*16-16+id ! for 2D transect 1 processor for 3d need beg or use filter_hydrologyc


  end function get_natveg_column_id
end module GridCellConnectionSetType

