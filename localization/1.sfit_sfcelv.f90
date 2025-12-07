program semivari_sfcelv
! calculates the semivariances along the river
! 
implicit none
character*128                         :: fname,buf,expdir,camadir,outdir,folddir,allocdir
character*20                          :: line
character*4                           :: mmdd
integer,parameter                     :: latpx=1320,lonpx=1500
real,dimension(lonpx,latpx)           :: rivwth,rivlen,nextdst
real*4,allocatable                    :: globaltrue(:,:,:),globalDr(:,:,:)
integer,allocatable                   :: ocean(:,:),days(:)
integer*4                             :: ind_cent,lon_cent,lat_cent,patch_size,patch_side,patch_nums
integer,parameter                     :: ind=1456 !wlv
!integer,parameter                     :: ind=1011 !discharge
integer*4,dimension(ind)              :: lon_ind,lat_ind
integer*4                             :: i,j,i_m,j_m,pixel
integer                               :: day,ios,countnum
real,allocatable                      :: xf(:)
integer*4,dimension(lonpx,latpx)      :: nextX,nextY
integer*4,dimension(lonpx,latpx,2)    :: nextXY
real                                  :: cov,corr,semivar,std
real,allocatable                      :: rlen(:)
integer                               :: N
integer,allocatable                   :: ux(:),uy(:),xt(:),yt(:)
real,allocatable                      :: svar(:),sd(:)
integer                               :: un,k
character*8                           :: lon,lat,u
logical(4)                            :: results 
integer                               :: ocean_int(lonpx*latpx)
!----------------
call getarg(1,buf)
read(buf,*) patch_size ! radius

call getarg(2,buf)
read(buf,*) N  ! days

call getarg(3,buf)
read(buf,"(A)") camadir
write(*,*) camadir

call getarg(4,buf)
read(buf,"(A)") outdir
write(*,*) outdir

call getarg(5,buf)
read(buf,"(A)") allocdir
write(*,*) allocdir

patch_side=patch_size*2+1
patch_nums=patch_side**2

! read river width
fname=trim(adjustl(camadir))//"map/tej_01min/rivwth.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivwth
    ! ocean is -9999
else
    write(*,*) "no file rivwth",fname
end if
close(34)

! read names of days in a particular year
allocate(days(N))
fname="yr2019.txt"
open(34,file=fname,form="formatted",iostat=ios)
if(ios==0)then
    read(34,*) days
else
    write(*,*) "no days",fname
end if
close(34)

! read river length
fname=trim(adjustl(camadir))//"map/tej_01min/rivlen.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) rivlen
    ! ocean is -9999
else
    write(*,*) "no file rivlen",fname
end if
close(34)

! read next grid information
! read nextX and nextY
fname=trim(adjustl(camadir))//"map/tej_01min/nextxy.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextXY
else
    write(*,*) "no file nextXY at:",fname
end if
nextX=nextXY(:,:,1)
nextY=nextXY(:,:,2)
close(34)

! make ocean mask from storage data (1 is ocean; 0 is not ocean)
ocean_int=transfer((nextX==-9999),ocean_int)
ocean=RESHAPE(ocean_int,[lonpx,latpx])

! read nextdst
fname=trim(adjustl(camadir))//"map/tej_01min/nxtdst.bin"
open(34,file=fname,form="unformatted",access="direct",recl=4*latpx*lonpx,status="old",iostat=ios)
if(ios==0)then
    read(34,rec=1) nextdst
    ! ocean is -9999
else
    write(*,*) "no file nextdst",fname
end if
close(34)
write(*,*) nextX(517,768),nextY(517,768)
write(*,*) ocean(517,768)

!--
! read true WSE
allocate(globaltrue(lonpx,latpx,N),globalDr(lonpx,latpx,N))
do day=1,N
    write(mmdd,'(i4.4)') days(day)
    fname=trim(adjustl(allocdir))//"/sim_wlv2019/E00_"//mmdd//".bin"
    open(34,file=fname,form="unformatted",access="stream",status="old",iostat=ios)
    if(ios==0)then
        read(34) globaltrue(:,:,day)
    else
        write(*,*) "no true",fname
        stop
    end if
    close(34)
end do

! read site info
fname=trim(adjustl(allocdir))//'/lon_index.txt'
open(34,file=fname,status='old',action='read',iostat=ios)
j=0
if(ios==0)then
    do i=1,ind
        read(34,*) line
        read(line,*) lon_ind(j+1)
        j=j+1
    end do
else
    write(*,*) "no file lon_ind:",fname
end if
close(34)

fname=trim(adjustl(allocdir))//'/lat_index.txt'
open(34,file=fname,form='formatted',iostat=ios)
j=0
if(ios==0)then
    do i=1,ind
        read(34,*) line
        read(line,*) lat_ind(j+1)
        j=j+1
    end do
else
    write(*,*) "no file lon_ind:",fname
end if
close(34)

allocate(xt(patch_nums),yt(patch_nums),xf(patch_nums))
allocate(ux(patch_nums),uy(patch_nums),rlen(patch_nums),svar(patch_nums),sd(patch_nums))

fname=trim(adjustl(outdir))//"/semivar/lonlat_list.txt"
open(79,file=fname,form="formatted",status='replace',iostat=ios)
write(79,'(a4,4x,a4,4x,a3,4x,a3)')"lon","lat","up","dn"

21  format(i4.4,2x,i4.4,2x,f8.2,2x,e14.5,2x,e14.5)
22  format(a4,2x,a4,2x,a8,2x,a14,2x,a14)
23  format(a4,4x,a4,4x,i5.5,4x,i5.5)
write(*,*) "start calculation"

!$OMP PARALLEL DO PRIVATE(lon_cent, lat_cent, countnum, lon, lat, un, ux, uy, u,folddir, fname, ios, xt, yt, rlen, k, semivar, std, svar, sd) SHARED(lon_ind,lat_ind, ocean, rivwth, nextX, nextY, lonpx, latpx, patch_size, outdir,globaltrue, N) 
do ind_cent = 1,ind
    lon_cent = lon_ind(ind_cent)
    lat_cent = lat_ind(ind_cent)
         !remove ocean
         if (ocean(lon_cent,lat_cent) /= 0) then
             cycle
         end if
         ! remove rivwth <= 0m
         if (rivwth(lon_cent,lat_cent) <=0.0) then
             cycle
         end if
         !--traget pixel
         countnum=1
         write(lon,'(i4.4)')lon_cent
         write(lat,'(i4.4)')lat_cent
         !$OMP CRITICAL
         write(*,*)"================================="
         write(*,*) lon,lat
         write(*,*)"*********************************"
         !$OMP END CRITICAL

         call up_river_pixel(lon_cent,lat_cent,lonpx,latpx,patch_size,ocean,rivwth,nextX,nextY,ux,uy,un)
         !$OMP CRITICAL
         write(*,*) 'check',un
         !$OMP END CRITICAL
         do i=1,un
             write(u,'(a2, i5.5)') "up",i
             folddir=trim(adjustl(outdir))//"/semivar/"//trim(lon)//trim(lat)//"/"
             if (i == 1) then
                call system('mkdir -p '//trim(folddir))
             end if
             fname=trim(adjustl(outdir))//"/semivar/"//trim(lon)//trim(lat)//"/"//trim(u)//".svg"
             open(34,file=fname,form="formatted",status='replace',iostat=ios)
             if (ios /= 0) then 
                 write(*,*) "Cannot make file"
             endif
             write(34,22)"lon","lat","dis","gamma","sig"
             xt=-9999
             yt=-9999

! 计算上游网格与目标网格间的距离
             call river_up(ux(i),uy(i),lon_cent,lat_cent,patch_size,lonpx,latpx,nextX,nextY,nextdst,xt,yt,rlen,k)
             do j=1,k
                 call semi_var(globaltrue(lon_cent,lat_cent,:),globaltrue(xt(j),yt(j),:),N,semivar,std)
                 svar(j)=semivar
                 sd(j)=std
                 write(34,21)xt(j),yt(j),rlen(j),semivar,std
             end do
             close(34)
         end do
         xt=-9999
         yt=-9999
! 计算下游网格与目标网格间的距离
         call river_dn(lon_cent,lat_cent,patch_size,lonpx,latpx,ocean,rivwth,nextX,nextY,nextdst,xt,yt,rlen,k) 
         write(u,'(a2, i5.5)') "dn",0
         fname=trim(adjustl(outdir))//"/semivar/"//trim(lon)//trim(lat)//"/"//trim(u)//".svg"
         open(34,file=fname,form="formatted",status='replace',iostat=ios)
         if (ios /= 0) then 
             !$OMP CRITICAL
             write(*,*) "Cannot make file"
             !$OMP END CRITICAL
             cycle
         endif
         write(34,22)"lon","lat","dis","gamma","sig"
         do j=1,k
             call semi_var(globaltrue(lon_cent,lat_cent,:),globaltrue(xt(j),yt(j),:),N,semivar,std)
             svar(j)=semivar
             sd(j)=std 
             write(34,21)xt(j),yt(j),rlen(j),semivar,std
         end do
         close(34)
         if (k>0) then 
             k=1
         else
             k=0
         end if
         !$OMP CRITICAL
         write(79,23)lon,lat,un,k
         write(*,*)"##################################"
         !$OMP END CRITICAL
end do
!$OMP END PARALLEL DO
deallocate(ocean,days,globaltrue,globalDr,xt,yt,xf,ux,uy,rlen,svar,sd) 
close(35)   
close(79) 
end program semivari_sfcelv
!*****************************************************************
function roundx(ix, nx)
implicit none
!-- for input -----------
integer                     :: ix, nx
!-- for output ----------
integer                     :: roundx
!------------------------
if (ix .ge. 1) then
  roundx = ix - int((ix -1)/nx)*nx
else
  roundx = nx - abs(mod(ix,nx))
end if 
return
end function roundx
!*****************************************************************
subroutine ixy2iixy(ix,iy, nx, ny, iix, iiy)
implicit none
!- for input -----------------
integer                   :: ix, iy, nx, ny
!- for output ----------------
integer                   :: iix, iiy,roundx
!-----------------------------
if (iy .lt. 1) then
  iiy = 2 - iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else if (iy .gt. ny) then
  iiy = 2*ny -iy
  iix = ix + int(nx/2.0)
  iix = roundx(iix, nx)
else
  iiy = iy
  iix = roundx(ix, nx)
end if
return
end subroutine ixy2iixy
!*****************************************************************
! uord: 从目标网格所在窗口中的网格出发，寻找下游网格，直到到达海洋、河嘴、陆地或者目标网格为止
subroutine uord(i,j,x,y,nx,ny,nextX,nextY,ud)
implicit none 
integer                     :: i,j,x,y,nx,ny
integer,dimension(nx,ny)    :: nextX,nextY
real,dimension(nx,ny)       :: rivlen
integer                     :: ix,iy,iix,iiy,tx,ty,pixel,ud
!--
tx=x
ty=y
ix=i
iy=j
ud=-1
! 从目标网格的周围网格出发，寻找下游网格
! 由于nextxy文件只有下游信息，所以我们需要找出上游信息
do while (ix/=tx .or. iy/=ty) 
  iix=ix
  iiy=iy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  if (ix==-9 .or. iy==-9) then
    ud=+1 
    exit
  end if
  if (ix == -10 .or. iy == -10) then  ! inland termination
    ud=+1 
    exit
  end if
  if (ix == -9999 .or. iy == -9999) then
    ud=+1
    exit
  end if
end do
!---
! 从目标网格出发，寻找下游网格，直到遇到窗口网格
! -9999：窗口网格和目标网格不在同一条河流，没有相遇
if (ud==+1) then
  !-
  tx=i
  ty=j
  ix=x
  iy=y
  do while (ix/=tx .or. iy/=ty) 
    iix=ix
    iiy=iy
    ix=nextX(iix,iiy)
    iy=nextY(iix,iiy)
    if (ix ==-9 .or. iy==-9) then
      ud=-9999
      exit
    end if 
    if (ix == -10 .or. iy == -10) then ! inland termination
      ud=-9999
      exit
    end if 
    if (ix == -9999 .or. iy == -9999) then ! ocean
      ud=-9999
    exit
    end if
  end do
end if
return
! ud=-1,窗口网格在目标网格上游
! ud=+1,窗口网格在目标网格下游
!---
end subroutine uord 
!*****************************************************************
! 计算上游网格到目标网格的距离
subroutine river_up(i,j,x,y,patch,nx,ny,nextX,nextY,nextdst,xpixel,ypixel,rlen,k)
implicit none 
!--
integer                             :: i,j,x,y,patch,nx,ny
integer,dimension(nx,ny)            :: nextX,nextY
real,dimension(nx,ny)               :: nextdst
!--
integer,dimension((2*patch+1)**2) :: xpixel,ypixel,lx,ly
real,dimension((2*patch+1)**2)    :: rlen,rlen1
integer                             :: ix,iy,iix,iiy,tx,ty,k,l
real                                :: length,rl
!--
xpixel = -9999
ypixel = -9999
!---
!write(*,*)i,i,x,y
rlen=-9999.0
rlen1=-9999.0
rl=0.0
length=0.0
!--
xpixel(1) =i
ypixel(1) =j
!---
rl=anint((nextdst(i,j)/1000.0)*100)/100.0
rlen1(1)=0.0
! k:窗口网格到目标网格共有多少个网格
k=2
!-
tx=x
ty=y
ix=i
iy=j
do while (ix/=tx .or. iy/=ty) 
  iix=ix
  iiy=iy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  if (ix==-9 .or. iy==-9) then
    exit
  end if
  if (ix==-10 .or. iy==-10) then
    exit
  end if
  if (ix==-9999 .or. iy==-9999) then
    exit
  end if  
  xpixel(k) =ix
  ypixel(k) =iy 
  rl=anint((nextdst(iix,iiy)/1000.0)*100)/100.0
  length=length+rl
  rlen1(k)=length
  k=k+1
end do
!---- 
k=k-1
rl=0.0
do l=1,k 
  rlen(l)=rlen1(k)-rlen1(k-l+1)
  lx(l)=xpixel(k-l+1)
  ly(l)=ypixel(k-l+1)
  end do 
xpixel=lx 
ypixel=ly
return
!---
end subroutine river_up
!*****************************************************************
! 计算下游网格和目标网格间的距离
subroutine river_dn(x,y,patch,nx,ny,ocean,rivwth,nextX,nextY,nextdst,xpixel,ypixel,rlen,k)
implicit none 
integer                             :: i,j,x,y,patch,nx,ny
integer,dimension(nx,ny)            :: ocean,nextX,nextY
real,dimension(nx,ny)               :: rivwth,nextdst
integer,dimension((2*patch+1)**2) :: xpixel,ypixel
real,dimension((2*patch+1)**2)    :: rlen
integer                             :: ix,iy,iix,iiy,tx,ty,k,p,q,ud,c,l
real                                :: length,rl
!--
xpixel = -9999
ypixel = -9999
xpixel(1) =x
ypixel(1) =y
rlen= 0.0
rl=0.0
length=0.0
!--
rl=anint((nextdst(x,y)/1000.0)*100)/100.0
rlen(1)=0.0
!-
k=2
ix=x
iy=y
do while (ix/=tx .or. iy/=ty) 
  iix=ix
  iiy=iy
  ix=nextX(iix,iiy)
  iy=nextY(iix,iiy)
  if (ix==-9 .or. iy==-9) then
    exit
  end if
  if (ix==-10 .or. iy==-10) then ! inland termination
    exit
  end if
  if (ix==-9999 .or. iy==-9999) then
    exit
  end if 
  xpixel(k) =ix
  ypixel(k) =iy 
  ! 计算目标网格与周围网格间的距离
  rl=anint((nextdst(iix,iiy)/1000.0)*100)/100.0
  length=length+rl!/2.0
  rlen(k)=length
  k=k+1
  if (k>patch) then
    exit
  end if  
end do
!-- 
k=k-1
return
!---
end subroutine river_dn
!*****************************************************************
subroutine up_river_pixel(x,y,nx,ny,patch,ocean,rivwth,nextX,nextY,ux,uy,un)
implicit none 
!--
integer                             :: i,j,x,y,nx,ny,patch
integer,dimension(nx,ny)            :: ocean,nextX,nextY
real,dimension(nx,ny)               :: rivwth
integer,dimension((2*patch+1)**2) :: ux,uy,ux_m,uy_m
integer                             :: i_m,j_m,ix,iy,iix,iiy,tx,ty,ucount,m
integer                             :: k,p,q,r,c,l,ud,g,un
integer,allocatable                 :: lx(:),ly(:) 

ux=-9999
uy=-9999
ux_m=-9999
uy_m=-9999
c=0
ucount=1
do k=1,patch
  p=patch - k + 1
  q=2*(2*p+1)+2*(2*(p-1)+1)
  allocate(lx(q),ly(q))

! outer_lxly: 将x,y,p,q转化为可读且合理的网格点
  call outer_lxly(x,y,p,q,nx,ny,lx,ly,c)
  do l=1,c
    i=lx(l)
    j=ly(l)
    if (ocean(i,j)/=0 .and. rivwth(i,j)<=0.) then
      cycle
    end if

!-- find up or down of x,y  
    call uord(i,j,x,y,nx,ny,nextX,nextY,ud)
    ! ud==-1,窗口网格在目标网格上游
    if (ud/=-1) then
      cycle
    end if 
    ux(ucount)=i
    uy(ucount)=j
    ucount=ucount+1
  end do
  deallocate(lx,ly)
end do

do m=1,ucount
  ix=ux(m)
  iy=uy(m)
  ! 排除目标网格
  if (ix==i .and. iy==j) then
    cycle
  end if 
  if (ix==-9999 .or. iy==-9999) then
    cycle
  end if 
    !-- if downstream of a selected pixel  
    call uord(ix,iy,i,j,nx,ny,nextX,nextY,ud)
    ! ud==1,窗口网格在目标网格下游,排除目标网格下游网格，只保留目标网格上游网格
  if (ud==1) then
    ux(m)=-9999
    uy(m)=-9999
  end if
end do 

un=0
k=1
do m=1,(2*patch+1)**2
  ! 窗口四条边上的网格
  if (ux(m)>0 .and. uy(m)>0) then
    ux_m(k)=ux(m)
    uy_m(k)=uy(m)
    un=un+1
    k=k+1
  end if
end do
ux=ux_m 
uy=uy_m  

! 统计窗口中有多少个上游网格
! un:目标网格的上游网格数量
return
!---
end subroutine up_river_pixel
!*****************************************************************
! 给出目标网格周围网格的索引
subroutine outer_lxly(x,y,p,q,nx,ny,lx,ly,c)
implicit none
! return the outer boundry of patch
integer                   :: x,y,p,q,nx,ny
integer,dimension(q)      :: lx,ly
integer                   :: c
integer                   :: i,j,i_m,j_m,c1 
!--
c=1
do i=x-p,x+p,2*p
  do j=y-p,y+p,1
    call ixy2iixy(i,j,nx,ny,i_m,j_m)
    lx(c)=i_m
    ly(c)=j_m
    c=c+1
  end do
end do
!--
do j=y-p,y+p,2*p
  do i=x-p+1,x+p-1,1
    call ixy2iixy(i,j,nx,ny,i_m,j_m)
    lx(c)=i_m
    ly(c)=j_m
    c=c+1
  end do
end do
c=c-1
return 
end subroutine outer_lxly 
!*****************************************************************
! 根据dis计算周围网格与目标网格之间的相关关系
subroutine semi_var(t,h,n,semivar,std)
implicit none
!--
integer                   :: n,i,pixel
real,dimension(n)         :: t,h
!--
real                      :: semivar,std
!-
real                      :: p,v
!------------------------------------
p=0.0
v=0.0
semivar=0.0
do i=1,n
  p=p+(h(i)-t(i))**2.0
  v=v+(((h(i)-t(i))**2.0)/2.0)**2.0
end do
semivar=p/(2.0*real(n))
std=sqrt((v/real(n-1))-((semivar**2.0)*(real(n)/real(n-1))))
return
end subroutine semi_var
! ***********************
