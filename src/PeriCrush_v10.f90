module bondBasedPeridynamics3D
use omp_lib
implicit none
! 单个块体的几何信息
type geometry3d
    real*8 :: length  ! 长度
    real*8 :: width   ! 宽度
    real*8 :: height  ! 高度
    integer*4 :: ndim ! 维度数
    integer*4 :: xn  ! x方向划分质点数目
    integer*4 :: yn  ! y方向划分质点数目
    integer*4 :: zn  ! z方向划分质点数目
    integer*4 :: xon ! x方向外部边界质点数目
    integer*4 :: yon ! y方向外部边界质点数目
    integer*4 :: zon ! z方向外部边界质点数目
    real*8, allocatable :: coor(:,:)   ! 物质点坐标
    real*8, allocatable :: vol(:,:)    ! 物质点体积
    real*8, allocatable :: rprad(:,:)  ! 物质点代表性大小
    real*8, allocatable :: rpside(:,:) ! 物质点的代表性边长
    real*8 :: dx ! 物质点x方向划分的长度
    real*8 :: dy ! 物质点y方向划分的长度
    real*8 :: dz ! 物质点z方向划分的长度
    real*8 :: delta_x ! 网格尺寸
    real*8 :: delta   ! 邻域半径
    real*8 :: m       ! Mvalue-m值
    integer*4 :: pn   ! 总的物质点数目
    integer*4 :: hpn  ! 总的邻域物质点数目
    integer*4, allocatable :: h(:)  ! 邻域
    integer*4, allocatable :: hn(:) ! 每个物质点的邻域物质点数目
    integer*4, allocatable :: hs(:) ! 开始位置
    integer*4, allocatable :: he(:) ! 结束位置
    integer*4, allocatable :: mark(:) ! 标记整数
    real*8, allocatable :: cpos(:,:) ! 当前物质点位置
    real*8, allocatable :: ivel(:,:) ! 初始速度
    real*8, allocatable :: cvel(:,:) ! 当前速度
    real*8, allocatable :: pvel(:,:) ! 上一步速度
    real*8, allocatable :: idis(:,:) ! 初始位移
    real*8, allocatable :: cdis(:,:) ! 当前位移
    real*8, allocatable :: pdis(:,:) ! 上一步位移
    real*8, allocatable :: iacc(:,:) ! 初始加速度
    real*8, allocatable :: cacc(:,:) ! 当前加速度
    real*8, allocatable :: cfdens(:,:) ! 当前力密度矢量
    real*8, allocatable :: pfdens(:,:) ! 上一步力密度矢量
    real*8, allocatable :: bfdens(:,:) ! 体积力密度矢量
    real*8 :: emod      !杨氏模量
    real*8 :: smod      !剪切模量
    real*8 :: bmod      !体积模量
    real*8 :: lameConst !拉梅常数
    real*8 :: pr        !泊松比
    real*8 :: c         !弹簧系数-理论值
    real*8 :: c1        !弹簧系数-（理论值-主方向）
    real*8 :: c2        !弹簧系数-（理论值-切方向）
    real*8 :: massdens  !质量密度
    real*8, allocatable :: cc(:,:)    ! 弹簧系数-（计算值-1主方向-2切方向）（第1列第2列）
    logical*1, allocatable :: fail(:) ! （定义键是否完好-完好为真-断裂为假）
    real*8, allocatable :: dmg(:,:)   ! 损伤（局部损伤表征）
    real*8 :: scr                     ! 临界伸长率
    real*8 :: energyReleaseRate       ! 能量释放率
    real*8 :: dt                      ! 时间步长
end type geometry3d
type localContact
    integer*4, allocatable :: cp(:, :) ! 接触对
    real*8, allocatable :: rs(:)       ! 接触距离
    real*8 :: c                        ! 接触刚度
    integer*4 :: ncp                   ! 接触对个数
    real*8, allocatable :: cf(:,:)     ! 接触力
    real*8, allocatable :: cft(:,:)    ! 切向接触力
    real*8 :: mu                       ! 静摩擦系数
    real*8 :: tcf(3)                   ! 总接触力
    real*8, allocatable :: overlap(:)  ! 初始重叠量
end type localContact

contains
!################################
!# SUBROUTINE [01] 方形网格生成 #
!################################
subroutine tetrahedronGrids(a)
implicit none
type(geometry3d), intent(inout) :: a
real*8 :: length ! x方向的长度
real*8 :: width  ! y方向的长度
real*8 :: height ! z方向的长度
integer*4 :: lenNum    ! x方向划分数目
integer*4 :: widNum     ! y方向划分数目
integer*4 :: heightNum    ! z方向划分数目
integer*4 :: xPN ! 单边x方排布点数目
integer*4 :: yPN ! 单边y方排布点数目
integer*4 :: zPN ! 单边z方排布点数目
integer*4 :: particleNumber
real*8 :: evol
real*8 :: tvol
integer*4, allocatable :: xNumberArray(:) ! x方向序列数
integer*4, allocatable :: yNumberArray(:) ! y方向序列数
integer*4, allocatable :: zNumberArray(:) ! z方向序列数
real*8, allocatable :: xcoordinate(:) ! x方向坐标序列
real*8, allocatable :: ycoordinate(:) ! y方向坐标序列
real*8, allocatable :: zcoordinate(:) ! z方向坐标序列
integer*4 :: i, j, k, m
length = a%length
width  = a%width
height = a%height
lenNum = a%xn
widNum  = a%yn
heightNum = a%zn
xPN = a%xn + 2*a%xon ! x总方向划分物质点数目（包含边界质点数目）
yPN = a%yn + 2*a%yon
zPN = a%zn + 2*a%zon
particleNumber = xPN * yPN * zPN! 总的物质点数
evol = length * width * height / ( (lenNum-1) * (widNum-1) * (heightNum-1) )
tvol = evol * particleNumber
allocate( xNumberArray( xPN ) )
allocate( yNumberArray( yPN ) )
allocate( zNumberArray( zPN ) )
allocate( xcoordinate( xPN ) )
allocate( ycoordinate( yPN ) )
allocate( zcoordinate( zPN ) )
if(mod(xPN,2).eq.1) then
    xNumberArray = (/(i, i = - (xPN-1)/2, (xPN-1)/2, 1)/)
    xcoordinate = length * xNumberArray/(lenNum-1)
else
    xNumberArray = (/(i, i = - xPN+1, xPN-1, 2)/)
    xcoordinate = length * xNumberArray/2/(lenNum-1)
endif
if(mod(yPN,2).eq.1) then
    yNumberArray = (/(i, i = -(yPN-1)/2, (yPN-1)/2, 1)/)
    ycoordinate = width * yNumberArray / ( widNum-1 )
else
    yNumberArray = (/(i, i = -yPN+1, yPN-1, 2)/)
    ycoordinate = width * yNumberArray/2/(widNum-1)
endif
if(mod(zPN,2).eq.1) then
    zNumberArray = (/(i, i = -(zPN-1)/2, (zPN-1)/2, 1)/)
    zcoordinate = height * zNumberArray / ( heightNum - 1 )
else
    zNumberArray = (/(i, i = -zPN+1, zPN-1, 2)/)
    zcoordinate = height * zNumberArray /2/(heightNum-1)
endif
allocate( a%coor( particleNumber, a%ndim ) )
allocate( a%vol( particleNumber, 1 ) )
do i = 1, xPN, 1
    do j = 1, yPN, 1
        do k = 1, zPN, 1
            m = (i-1)*yPN*zPN + (j-1)*zPN + k
            a%coor(m, 1) = xcoordinate(i)
            a%coor(m, 2) = ycoordinate(j)
            a%coor(m, 3) = zcoordinate(k)
            a%vol (m, 1) = evol
        enddo
    enddo
enddo
a%pn = particleNumber
a%delta_x = evol**(1d0/3d0)
a%delta = a%m * a%delta_x
allocate ( a%rprad ( a%pn, 1 ) )
allocate ( a%rpside ( a%pn, 1 ) )
a%rprad = (evol*3d0/ ( 4d0 * dacos(-1d0)) )**(1d0/3d0)
a%dx = a%length / (a%xn - 1)
a%dy = a%width  / (a%yn - 1)
a%dz = a%height / (a%zn - 1)
a%rpside = a%delta_x
! 释放不必要的内存
deallocate( xNumberArray )
deallocate( yNumberArray )
deallocate( zNumberArray )
deallocate( xcoordinate )
deallocate( ycoordinate )
deallocate( zcoordinate )
end subroutine tetrahedronGrids

subroutine regHorizonSearch(a)
implicit none
type(geometry3d), intent(inout) :: a
integer*4 :: xPN ! 单边x方排布点数目
integer*4 :: yPN ! 单边y方排布点数目
integer*4 :: zPN ! 单边z方排布点数目
integer*4 :: floorMvalue     ! m值的取整
integer*4, allocatable :: xNumberArray(:) ! x方向序列数
integer*4, allocatable :: yNumberArray(:) ! y方向序列数
integer*4, allocatable :: zNumberArray(:) ! z方向序列数
real*8, allocatable :: xcoordinate(:)     ! x方向坐标序列
real*8, allocatable :: ycoordinate(:)     ! y方向坐标序列
real*8, allocatable :: zcoordinate(:)     ! z方向坐标序列
real*8 :: totalParticleVolume             ! 所有物质点的体积
real*8 :: singleParticleVolume            ! 单个物质点的体积
! 标准Horizon的几何信息
real*8 :: elength  ! x方向的长度
real*8 :: ewidth   ! y方向的长度
real*8 :: eheight  ! z方向的长度
integer*4 :: eparticleNumber  ! 标准horizon的物质点个数
integer*4 :: elenNum    ! x方向划分数目
integer*4 :: ewidNum     ! y方向划分数目
integer*4 :: eheightNum    ! z方向划分数目
integer*4 :: exPN ! 单边x方排布点数目
integer*4 :: eyPN ! 单边y方排布点数目
integer*4 :: ezPN ! 单边z方排布点数目
integer*4 :: efloorMvalue       ! m值的取整
real*8 :: etotalParticleVolume  ! 所有物质点的体积
real*8 :: esingleParticleVolume ! 单个物质点的体积
real*8 :: dx, dy, dz !
integer*4 :: kx0, ky0, kz0, n0 !
integer*4 :: i1, j1, k1
integer*4 :: i, j, k
integer*4 :: m, n
integer*4, allocatable :: eboxH(:,:), eSpHH(:,:)
integer*4 :: eXdN, eYdN, eZdN
integer*4 :: nStart, nEnd
integer*4, allocatable :: horizon(:)
integer*4, allocatable :: horPN(:)
real*8 :: bondLength
real*8 :: delta
integer*4 :: maxSinPNN
integer*4 :: particleNumber
integer*4 :: SinPNN
integer*4 :: totNPN
dx = a%dx
dy = a%dy
dz = a%dz
particleNumber = a%pn
xPN = a%xn + 2*a%xon
yPN = a%yn + 2*a%yon
zPN = a%zn + 2*a%zon
delta = a%delta
eXdN = ceiling(delta/dx)
eYdN = ceiling(delta/dy)
eZdN = ceiling(delta/dz)
kx0 = xPN/2
ky0 = yPN/2
kz0 = zPN/2
allocate(eboxH((2*eXdN+1)*(2*eYdN+1)*(2*eZdN+1),3))
maxSinPNN = 0
do i = - eXdN, eXdN, 1
    do j = - eYdN, eYdN, 1
        do k = - eZdN, eZdN, 1
            i1 = kx0 + i 
            j1 = ky0 + j
            k1 = kz0 + k
            n = (i1-1)*yPN*zPN + (j1-1)*zPN + k1
            n0 = (kx0-1)*yPN*zPN + (ky0-1)*zPN + kz0
            bondLength = dsqrt( &
            & (a%coor(n,1) - a%coor(n0,1))**2 + &
            & (a%coor(n,2) - a%coor(n0,2))**2 + &
            & (a%coor(n,3) - a%coor(n0,3))**2 )
            if((bondLength.le.delta).and.(n.ne.n0)) then
                maxSinPNN = maxSinPNN + 1
                eboxH(maxSinPNN, 1) = i
                eboxH(maxSinPNN, 2) = j
                eboxH(maxSinPNN, 3) = k
            endif
        enddo
    enddo
enddo
allocate(eSpHH(maxSinPNN, 3))
do i = 1, maxSinPNN, 1
    eSpHH(i,1) = eboxH(i, 1)
    eSpHH(i,2) = eboxH(i, 2)
    eSpHH(i,3) = eboxH(i, 3)
enddo
allocate(horizon(maxSinPNN * particleNumber))
allocate(horPN(particleNumber))
horizon = 0
horPN = 0
totNPN = 0
do i = 1, xPN, 1
    do j = 1, yPN, 1
        do k = 1, zPN, 1
            n0 = (i-1)*yPN*zPN + (j-1)*zPN + k
            SinPNN = 0
            do m = 1, maxSinPNN, 1
                i1 = i + eSpHH(m,1)
                j1 = j + eSpHH(m,2)
                k1 = k + eSpHH(m,3)
                if(((i1.le.xPN).and.(i1.ge.1)) .and. &
                & ((j1.le.yPN).and.(j1.ge.1)) .and.&
                & ((k1.le.zPN).and.(k1.ge.1))) then
                    n = (i1-1)*yPN*zPN + (j1-1)*zPN + k1
                    totNPN = totNPN + 1
                    SinPNN = SinPNN + 1
                    horizon(totNPN) = n
                endif
            enddo
            horPN(n0) =  SinPNN
        enddo
    enddo
enddo
a%hpn = totNPN
allocate( a%h ( a%hpn ))
allocate( a%hn ( a%pn) )
allocate( a%hs ( a%pn) )
allocate( a%he ( a%pn) )
a%h(1:a%hpn:1) = horizon(1:a%hpn:1)
a%hn( 1:a%pn:1 ) = horPN ( 1:a%pn:1 )
nEnd = 0
nStart = 1
do i = 1, particleNumber, 1
    nEnd = nEnd + horPN(i)
    a%hs(i) = nStart
    a%he(i) = nEnd
    nStart  = nStart + horPN(i)
enddo
! 释放不必要的内存
deallocate(eboxH)
deallocate(eSpHH)
deallocate(horizon)
deallocate(horPN)
end subroutine regHorizonSearch

function particleDelete(a, particle2bdeleted)result(b)
implicit none
! 输入变量
type(geometry3d), intent(inout) :: a
logical :: particle2bDeleted(a%pn)
type(geometry3d) :: b
integer*4, allocatable :: halfHorizon(:)
integer*4, allocatable :: halfStart(:)
integer*4, allocatable :: halfEnd(:)
integer*4, allocatable :: deleteParticleOrder(:)
integer*4, allocatable :: restParticleOrder(:)
integer*4, allocatable :: particleNewOrder(:)
integer*4 :: deleteParticleNumber
integer*4 :: restParticleNumber
integer*4 :: restHorizonLength
integer*4 :: i, j, k, i1, j1, k1, m, n, nCount, nStart, nEnd
integer*4 :: nCount1
deleteParticleNumber = count(particle2bDeleted)
restParticleNumber = a%pn - deleteParticleNumber
allocate(deleteParticleOrder(deleteParticleNumber))
allocate(restParticleOrder(restParticleNumber))
deleteParticleOrder = pack((/(i,i=1, a%pn, 1)/), particle2bDeleted)
restParticleOrder = pack((/(i,i=1, a%pn, 1)/), .not.particle2bDeleted)
restHorizonLength = sum(a%hn(restParticleOrder))
allocate( halfHorizon( restHorizonLength ))
allocate( halfStart( restParticleNumber ) )
allocate( halfEnd( restParticleNumber ) )
allocate( particleNewOrder(a%pn) )
! 定义物质点新编号
nCount = 0
do i = 1, a%pn, 1
    if(particle2bDeleted(i))then
        particleNewOrder(i) = 0
    else
        nCount = nCount + 1
        particleNewOrder(i) = nCount
    endif
enddo
! 去掉所去除点的所有邻域物质点
halfStart = 0
halfEnd = 0
nStart = 1
nEnd = 0
do i = 1, restParticleNumber, 1
    ! 将第i个剩余物质点j的邻域物质点找出来赋值给halfHorizon
    j = restParticleOrder(i) 
    nEnd = nEnd + ( a%he(j) - a%hs(j) + 1)
    halfStart(i) = nStart
    halfEnd(i) = nEnd
    halfHorizon( nStart:nEnd:1 ) = a%h( a%hs(j):a%he(j):1 )
    nStart = nStart + ( a%he(j) - a%hs(j) + 1)
enddo
! 将剩余物质点的邻域质点是可以删除物质点的编号设为0
! 重设物质点顺序
do i = 1, restHorizonLength, 1
    j = halfHorizon(i)
    if(particle2bDeleted(j))then
        halfHorizon(i) = 0
    else
        halfHorizon(i) = particleNewOrder(j)
    endif
enddo
! 重新选取物质点的邻域（将邻域密实）
b%hpn  = count(halfHorizon.ne.0)
b%pn   = restParticleNumber
b%ndim = a%ndim
! 李代桃僵
allocate(b%coor(b%pn,b%ndim))
b%coor(:,:) = a%coor(restParticleOrder(:),:)
deallocate(a%coor)
allocate(b%vol(b%pn,1))
b%vol(:,:) = a%vol(restParticleOrder(:),:)
deallocate(a%vol)
b%m = a%m
b%delta_x = a%delta_x
b%delta = a%delta
allocate ( b%rprad ( b%pn, 1 ) )
allocate ( b%rpside ( b%pn, 1 ) )
b%rprad(:,:) = a%rprad(restParticleOrder(:),:)
deallocate(a%rprad)
b%rpside(:,:) = a%rpside(restParticleOrder(:),:)
deallocate(a%rpside)
b%length = a%length
b%width = a%width
b%height = a%height
b%dx = a%dx
b%dy = a%dy
b%dz = a%dz
b%xn = a%xn 
b%xon = a%xon
b%yn = a%yn 
b%yon = a%yon
b%zn = a%zn 
b%zon = a%zon
deallocate(a%h)
deallocate(a%hn)
deallocate(a%hs)
deallocate(a%he)
allocate(b%h(b%hpn))
allocate(b%hn(b%pn))
allocate(b%hs(b%pn))
allocate(b%he(b%pn))
nCount = 0
nCount1 = 0
do i = 1, b%pn, 1
    nStart = halfStart(i)
    nEnd = halfEnd(i)
    nCount1 = 0
    do j = nStart, nEnd, 1
        if(halfHorizon(j).ne.0)then
            nCount = nCount + 1
            b%h(nCount) = halfHorizon(j)
            nCount1 = nCount1 + 1
        endif
    enddo
    b%hn(i) = nCount1
enddo
nEnd = 0
nStart = 1
do i = 1, b%pn, 1
    nEnd = nEnd + b%hn(i)
    b%hs(i) = nStart
    b%he(i) = nEnd
    nStart = nStart + b%hn(i)
enddo
! 去掉不必要的内存
deallocate( deleteParticleOrder )
deallocate( restParticleOrder )
deallocate( halfHorizon )
deallocate( halfStart )
deallocate( halfEnd )
deallocate( particleNewOrder )
end function particleDelete
! 将模型完全清空
subroutine modelDelete(a)
implicit none
type(geometry3d), intent(inout) :: a
a%length = 0d0    ! 长度
a%width  = 0d0    ! 宽度
a%height = 0d0    ! 高度
a%ndim   = 3      ! 维度数
a%xn     = 0      ! x方向划分质点数目
a%yn     = 0      ! y方向划分质点数目
a%zn     = 0      ! z方向划分质点数目
a%xon    = 0      ! x方向外部边界质点数目
a%yon    = 0      ! y方向外部边界质点数目
a%zon    = 0      ! z方向外部边界质点数目
if(allocated( a%coor ))    deallocate( a%coor )   ! 物质点坐标
if(allocated( a%vol ))     deallocate( a%vol )    ! 物质点体积
if(allocated( a%rprad ))   deallocate( a%rprad )  ! 物质点代表性大小
if(allocated( a%rpside ))  deallocate( a%rpside ) ! 物质点的代表性边长
a%dx      = 0d0                                   ! 物质点x方向划分的长度
a%dy      = 0d0                                   ! 物质点y方向划分的长度
a%dz      = 0d0                                   ! 物质点z方向划分的长度
a%delta_x = 0d0                                   ! 网格尺寸
a%delta   = 0d0                                   ! 邻域半径
a%m       = 3.015d0                               ! m值
a%pn      = 0                                     ! 总的物质点数目
a%hpn     = 0                                     ! 总的邻域物质点数目
if(allocated( a%h ))    deallocate( a%h )     ! 邻域
if(allocated( a%hn ))   deallocate( a%hn )    ! 每个物质点的邻域物质点数目
if(allocated( a%hs ))   deallocate( a%hs )    ! 开始位置
if(allocated( a%he ))   deallocate( a%he )    ! 结束位置
if(allocated( a%mark )) deallocate( a%mark )  ! 标记整数
if(allocated( a%cpos )) deallocate( a%cpos )
if(allocated( a%ivel )) deallocate( a%ivel )
if(allocated( a%cvel )) deallocate( a%cvel )
if(allocated( a%pvel )) deallocate( a%pvel )
if(allocated( a%idis )) deallocate( a%idis )
if(allocated( a%cdis )) deallocate( a%cdis )
if(allocated( a%pdis )) deallocate( a%pdis )
if(allocated( a%iacc )) deallocate( a%iacc )
if(allocated( a%cacc )) deallocate( a%cacc )
if(allocated( a%cfdens )) deallocate( a%cfdens )
if(allocated( a%pfdens )) deallocate( a%pfdens )
if(allocated( a%bfdens )) deallocate( a%bfdens )
a%emod      = 0d0
a%smod      = 0d0
a%bmod      = 0d0
a%lameConst = 0d0
a%pr        = 0d0
a%c         = 0d0
a%massdens  = 0d0
if(allocated( a%cc ))    deallocate( a%cc )
if(allocated( a%fail ))  deallocate( a%fail )
if(allocated( a%dmg ))   deallocate( a%dmg )
a%scr               = 0d0
a%energyReleaseRate = 0d0
a%dt                = 0d0
end subroutine modelDelete

subroutine showModelGeometryInfo(a)
type(geometry3d), intent(in) :: a
write(*, *) '模型[未修改前]长度', a%length
write(*, *) '模型[未修改前]宽度', a%width
write(*, *) '模型[未修改前]高度', a%height
write(*, *) '维度数', a%ndim
write(*, *) 'x方向划分质点数目', a%xn
write(*, *) 'y方向划分质点数目', a%yn
write(*, *) 'z方向划分质点数目', a%zn
write(*, *) 'x方向外部边界质点数目', a%xon
write(*, *) 'y方向外部边界质点数目', a%yon
write(*, *) 'z方向外部边界质点数目', a%zon
write(*, *) 'xmax = ', maxval(a%coor(:,1)) 
write(*, *) 'xmin = ', minval(a%coor(:,1))
write(*, *) 'ymax = ', maxval(a%coor(:,2)) 
write(*, *) 'ymin = ', minval(a%coor(:,2)) 
write(*, *) 'zmax = ', maxval(a%coor(:,3)) 
write(*, *) 'zmin = ', minval(a%coor(:,3)) 
write(*, *) '物质点总体积：', sum(a%vol(:,1))
write(*, *) '物质点平均体积：', sum(a%vol(:,1))/a%pn
write(*, *) '最大代表半径：', maxval(a%rprad)
write(*, *) '最小代表半径：', minval(a%rprad)
write(*, *) '最大代表边长：', maxval(a%rpside)
write(*, *) '最小代表边长：', minval(a%rpside)
write(*, *) '物质点x方向划分的长度', a%dx
write(*, *) '物质点y方向划分的长度', a%dy
write(*, *) '物质点z方向划分的长度', a%dz
write(*, *) '网格尺寸', a%delta_x
write(*, *) '邻域半径', a%delta
write(*, *) 'm值', a%m
write(*, *) '总的物质点数目', a%pn
write(*, *) '总的邻域物质点数目', a%hpn
write(*, *) '单个邻域物质点最大数目', maxval(a%hn)
write(*, *) '单个邻域物质点最小数目', minval(a%hn)
end subroutine showModelGeometryInfo

subroutine kineticVariableAllocate(a)
type(geometry3d), intent(inout) :: a
allocate(a%dmg(a%pn,1))
allocate(a%bfdens(a%pn,a%ndim))
allocate(a%cfdens(a%pn,a%ndim))
allocate(a%iacc(a%pn,a%ndim))
allocate(a%cacc(a%pn,a%ndim))
allocate(a%pvel(a%pn,a%ndim))
allocate(a%cvel(a%pn,a%ndim))
allocate(a%ivel(a%pn,a%ndim))
allocate(a%pdis(a%pn,a%ndim))
allocate(a%cdis(a%pn,a%ndim))
allocate(a%idis(a%pn,a%ndim))
allocate(a%cpos(a%pn,a%ndim))
allocate(a%fail(a%hpn))
a%dmg    = 0.0d0
a%bfdens = 0.0d0
a%cfdens = 0.0d0
a%iacc   = 0.0d0
a%cacc   = 0.0d0
a%pvel   = 0.0d0
a%cvel   = 0.0d0
a%ivel   = 0.0d0
a%pdis   = 0.0d0
a%cdis   = 0.0d0
a%idis   = 0.0d0
a%cpos   = 0.0d0
a%fail   = .true.
end subroutine kineticVariableAllocate

subroutine modelElasticConstant(a)
implicit none
type(geometry3d), intent(inout) :: a
a%smod = a%emod/( 2d0*(a%pr+1d0) )
a%bmod = a%emod/( 3d0*(-2d0*a%pr+1d0) )
a%c = 18d0*a%bmod/(dacos(-1d0)*a%delta**4)
a%c1 = 6d0 * a%emod /(dacos(-1d0) * a%delta**4 * (1d0 - 2d0*a%pr))
a%c2 = 6d0 * a%emod * (1d0-4d0*a%pr) &
& / (dacos(-1d0) * a%delta**4 * (1d0 - 2d0*a%pr) * (1d0+a%pr))
a%lameConst = a%pr * a%emod / (1d0 + a%pr) / (1d0 - 2d0*a%pr)
end subroutine modelElasticConstant

subroutine surfaceCorrection(a)
implicit none
type(geometry3d), intent(inout) :: a
real*8 :: bulkStrainEnergy
real*8 :: totalBondEnergy
real*8 :: ibl
real*8 :: cbl
integer*4 :: particleNumber
real*8 :: delta
real*8 :: delta_x
real*8 :: vcf
real*8 :: bondStretch
integer*4 :: i, j, k
allocate(a%cc(a%pn,1))
particleNumber = a%pn
delta    = a%delta
delta_x  = a%delta_x
a%cdis   = 1d-4*a%coor
bulkStrainEnergy = 0.5d0*a%bmod*(3d0*1d-4)**2
a%cpos   = a%coor + a%cdis
a%cc     = 0d0
do i = 1, a%pn, 1
    totalBondEnergy = 0d0
    do j = a%hs(i), a%he(i), 1
        k = a%h(j)
        ibl = dsqrt( &
        & ( a%coor(k,1) -  a%coor(i,1))**2 + &
        & ( a%coor(k,2) -  a%coor(i,2))**2 + &
        & ( a%coor(k,3) -  a%coor(i,3))**2 )
        cbl =  dsqrt( &
        & ( a%cpos(k,1) - a%cpos(i,1))**2 + &
        & ( a%cpos(k,2) - a%cpos(i,2))**2 + &
        & ( a%cpos(k,3) - a%cpos(i,3))**2 )
        if(ibl.le.delta-delta_x*0.5d0)then
            vcf = 1d0
        else
            vcf = (delta - ibl + delta_x*0.5d0)/delta_x
        endif
        bondStretch = (cbl-ibl)/ibl
        totalBondEnergy = totalBondEnergy + &
        & 0.25d0*ibl*bondStretch*bondStretch * &
        & vcf*a%vol(k,1)
    enddo
    a%cc(i,1) = bulkStrainEnergy/totalBondEnergy
enddo
a%cdis = 0.0d0
a%cpos = 0.0d0
write(*,*)'Average Spring Constant:', sum(a%cc)/particleNumber
write(*,*)'Maximun Spring Constant:', maxval(a%cc(:,1))
write(*,*)'Minimun Spring Constant:', minval(a%cc(:,1))
write(*,*)'theoretical Spring Constant:', a%c
end subroutine surfaceCorrection

subroutine bondBasedForceDensityVector(a)
implicit none
type(geometry3d) :: a
real*8 :: bulkStrainEnergy
real*8 :: totalBondEnergy
real*8 :: ibl
real*8 :: cbl
integer*4 :: particleNumber
real*8 :: delta
real*8 :: delta_x
real*8 :: vcf
real*8 :: bondStretch
integer*4 :: i, j, k
real*8 :: dmg, dmg0
particleNumber = a%pn
delta = a%delta
delta_x = a%delta_x
a%cpos = a%coor + a%cdis
a%cfdens = 0d0
!$ call omp_set_num_threads(16)
!$omp parallel do default(none) private(i, j, k, vcf, bondStretch, &
!$omp&dmg, dmg0, ibl, cbl)   &
!$omp&shared(particleNumber, a, delta, delta_x)
do i = 1, particleNumber, 1
    dmg = 0d0
    dmg0 = 0d0
    do j = a%hs(i), a%he(i), 1
        k = a%h(j)
        ibl = dsqrt( &
        & ( a%coor(k,1) -  a%coor(i,1))**2 + &
        & ( a%coor(k,2) -  a%coor(i,2))**2 + &
        & ( a%coor(k,3) -  a%coor(i,3))**2 )
        cbl = dsqrt( &
        & ( a%cpos(k,1) - a%cpos(i,1))**2 + &
        & ( a%cpos(k,2) - a%cpos(i,2))**2 + &
        & ( a%cpos(k,3) - a%cpos(i,3))**2 )
        if(ibl.le.delta-delta_x*0.5d0)then
            vcf = 1d0
        else
            vcf = (delta - ibl + delta_x*0.5d0)/delta_x
        endif
        bondStretch = (cbl-ibl)/ibl
        if(a%fail(j).and.bondStretch.gt.a%scr)then
            a%fail(j) = .false.
        endif
        dmg0 = dmg0 + a%vol(k,1) * vcf
        if(a%fail(j))then
            a%cfdens(i, 1) = a%cfdens(i, 1) + 0.5d0*(a%cc(i,1) + a%cc(k,1)) * &
            & bondStretch * (a%cpos(k,1) - a%cpos(i,1))/cbl * a%vol(k,1) * &
            & vcf
            a%cfdens(i, 2) = a%cfdens(i, 2) + 0.5d0*(a%cc(i,1) + a%cc(k,1)) * &
            & bondStretch * (a%cpos(k,2) - a%cpos(i,2))/cbl * a%vol(k,1) * &
            & vcf
            a%cfdens(i, 3) = a%cfdens(i, 3) + 0.5d0*(a%cc(i,1) + a%cc(k,1)) * &
            & bondStretch * (a%cpos(k,3) - a%cpos(i,3))/cbl * a%vol(k,1) * &
            & vcf
            dmg = dmg + a%vol(k,1) * vcf
        endif
    enddo
    a%dmg(i,1) = 1d0 - dmg/dmg0
enddo
!$omp end parallel do
a%cfdens = a%cfdens + a%bfdens
end subroutine bondbasedForceDensityVector

subroutine contactForceDensity(a, b, c)
implicit none
type(geometry3d), intent(inout) :: a, b
type(localContact), intent(inout) :: c
integer*4 :: i, j, k
real*8 :: rx, ry, rz, rs, r, r2
real*8 :: fabx, faby, fabz, fab
real*8 :: rdx, rdy, rdz, rd ! 相对位移矢量
real*8 :: rdtx, rdty, rdtz ! 切向方向相对位移矢量
real*8 :: nx, ny, nz, rdn
real*8 :: rdt2 ! 切向方向相对位移矢量大小
integer*4 :: ia, ib
c%cf = 0d0
!!$ call omp_set_num_threads(28)
!!$omp parallel do default(none) private(i, ia, ib, rx, ry, rz, &
!!$omp & rs, r2, r, fab, fabx, faby, fabz)   &
!!$omp & shared(c, b, a)
do i = 1, c%ncp, 1
    ia = c%cp(i, 1)
    ib = c%cp(i, 2)
    rx = b%cpos(ib,1) - a%cpos(ia,1)
    ry = b%cpos(ib,2) - a%cpos(ia,2)
    rz = b%cpos(ib,3) - a%cpos(ia,3)
    rs = c%rs(i)
    r2 = rx*rx + ry*ry + rz*rz
    r = dsqrt(r2)
    if(c%overlap(i) + r .lt. 2d0*rs) then
        ! 接触刚度新的定义方法
        c%c = a%cc(ia,1)*b%cc(ib,1)/( a%cc(ia,1)+b%cc(ib,1) )*1.0d0 
        fab = c%c * (0.5* ( r+c%overlap(i) )/rs - 1d0)
        fabx = rx/r * fab
        faby = ry/r * fab
        fabz = rz/r * fab
        !!$omp critical ! 写入操作怕出问题
        a%cfdens(ia, 1) = a%cfdens(ia, 1) + fabx * b%vol(ib,1)
        a%cfdens(ia, 2) = a%cfdens(ia, 2) + faby * b%vol(ib,1)
        a%cfdens(ia, 3) = a%cfdens(ia, 3) + fabz * b%vol(ib,1)
        ! 相对位置x方向的力密度矢量分量
        b%cfdens(ib, 1) = b%cfdens(ib, 1) - fabx * a%vol(ia,1) 
        b%cfdens(ib, 2) = b%cfdens(ib, 2) - faby * a%vol(ia,1)
        b%cfdens(ib, 3) = b%cfdens(ib, 3) - fabz * a%vol(ia,1)
        c%cf(i, 1) = fabx * b%vol(ib,1) * a%vol(ia,1)
        c%cf(i, 2) = faby * b%vol(ib,1) * a%vol(ia,1)
        c%cf(i, 3) = fabz * b%vol(ib,1) * a%vol(ia,1)
        !!$omp end critical
    endif
enddo
!!$omp end parallel do
!!$ call omp_set_num_threads(3)
!$omp parallel do default(none) private(i) shared(c)
do i = 1, 3, 1 
    c%tcf(i) = sum(c%cf(:,i)+c%cft(:,i))
enddo
!$omp end parallel do
end subroutine contactForceDensity

subroutine timeIntegration(a)
! 该函数负责时间积分部分
implicit none
type(geometry3d), intent(inout) :: a
a%cacc = a%cfdens/a%massDens
a%cdis = 2d0*a%idis - a%pdis + a%cacc * a%dt * a%dt
a%pdis = a%idis
a%idis = a%cdis
end subroutine timeIntegration

subroutine modelTranslation(a, vec)
! 该函数可以将模型平移
implicit none
type(geometry3d), intent(inout) :: a
real*8, intent(in) :: vec(a%ndim)
integer*4 :: i
do i = 1, a%ndim, 1
    a%coor(:,i) = a%coor(:,i) + vec(i)
enddo
end subroutine modelTranslation

subroutine contactPairs(a, b, pa, pb, r, c)
! 该子程序通过两个模型a和b以及可能的接触对pa,pb
! 对识别半径内的物质点进行接触识别并得到接触的重叠量
implicit none
type(geometry3d), intent(inout) :: a, b ! 两个物体
type(localContact), intent(inout) :: c ! 一种接触对
real*8 :: r ! 识别接触对半径
logical :: pa(a%pn), pb(b%pn) ! 选出接触可能的质点
integer*4, allocatable :: mcp(:, :) ! mayContactPair
logical, allocatable :: isWithinRadius(:)
integer*4 :: mxContactNumber, contactPairNumber
integer*4 :: i, j, k, nCount, minpos, minpos0(1)
real*8 :: ztrans,ztrans1,ztrans2
real*8, allocatable :: pdis(:)
mxContactNumber = count(pa)*count(pb)
allocate(mcp(mxContactNumber,2))
allocate(isWithinRadius(mxContactNumber))
nCount = 0
do i = 1, a%pn, 1
    do j = 1, b%pn, 1
        if(pa(i).and.pb(j)) then
            nCount = nCount + 1
            mcp(nCount,1:2:1) = (/i, j/)
        endif
    enddo
enddo
isWithinRadius = &
& (a%coor(mcp(:,1),1) - b%coor(mcp(:,2),1))**2 + &
& (a%coor(mcp(:,1),2) - b%coor(mcp(:,2),2))**2 + &
& (a%coor(mcp(:,1),3) - b%coor(mcp(:,2),3))**2 .le. &
& r**2
c%ncp = count(isWithinRadius)
allocate(c%cp(c%ncp,2))
allocate(c%rs(c%ncp))
nCount = 0
do i = 1, mxContactNumber, 1
    if(isWithinRadius(i))then
        nCount = nCount + 1
        c%cp(nCount,:) = mcp(i,:)
    endif
enddo
c%rs = 0.5d0* (a%rprad(c%cp(:,1),1) + b%rprad(c%cp(:,2),1))
deallocate(mcp)
deallocate(isWithinRadius)
allocate(c%cf(c%ncp, 3))          ! 接触力
allocate(c%cft(c%ncp, 3))         ! 法向接触力
c%c = (a%c * b%c) / (a%c+b%c)*5d0 ! 接触刚度 ! 取5倍
allocate(c%overlap(c%ncp))
allocate(pdis(c%ncp))
pdis = dsqrt ( &
& (a%coor(c%cp(:,1),1) - b%coor(c%cp(:,2),1))**2 + &
& (a%coor(c%cp(:,1),2) - b%coor(c%cp(:,2),2))**2 + &
& (a%coor(c%cp(:,1),3) - b%coor(c%cp(:,2),3))**2 ) - 2d0*c%rs
minpos0 = minloc(pdis)
write(*,*)'初始最大重叠量：', minval(pdis)
minpos = minpos0(1)
ztrans1 = dsqrt( 4d0 * c%rs(minpos)**2 - &
& (a%coor(c%cp(minpos,1),1) - b%coor(c%cp(minpos,2),1))**2 - &
& (a%coor(c%cp(minpos,1),2) - b%coor(c%cp(minpos,2),2))**2 ) + &
& a%coor(c%cp(minpos,1), 3) - b%coor(c%cp(minpos,2), 3)
ztrans2 = -dsqrt( 4d0 * c%rs(minpos)**2 - &
& (a%coor(c%cp(minpos,1),1) - b%coor(c%cp(minpos,2),1))**2 - &
& (a%coor(c%cp(minpos,1),2) - b%coor(c%cp(minpos,2),2))**2 ) + &
& a%coor(c%cp(minpos,1), 3) - b%coor(c%cp(minpos,2), 3)
if ((b%coor(c%cp(minpos,2),3)-a%coor(c%cp(minpos,1),3))*ztrans1.ge.0) then
    ztrans = ztrans1
endif
if ((b%coor(c%cp(minpos,2),3)-a%coor(c%cp(minpos,1),3))*ztrans2.ge.0) then
    ztrans = ztrans2
endif
!call modelTranslation(b, (/0d0,0d0,ztrans/)) ! 移到最恰当的地方
pdis = &
& (a%coor(c%cp(:,1),1) - b%coor(c%cp(:,2),1))**2 + &
& (a%coor(c%cp(:,1),2) - b%coor(c%cp(:,2),2))**2 + &
& (a%coor(c%cp(:,1),3) - b%coor(c%cp(:,2),3))**2 
where (pdis .lt. 4d0*c%rs**2)
    c%overlap = 2d0*c%rs - dsqrt(pdis) ! 计算初始重叠量
endwhere
write(*,*)'最大重叠量',maxval(c%overlap)
end subroutine contactPairs

!##################################
!############## 创建文件夹 ########
!##################################
subroutine creatfolder(foldername)
! 该子程序负责创建模型和计算结果输出的文件夹
implicit none
logical :: istat
real :: a
integer :: b = 0
character(len=100), intent(inout) :: foldername
write(*,*)'为了保存输出结果，请输入一个文件夹名：'
read(*,*)foldername
inquire(file='./'//trim(adjustl(foldername)),exist = istat)
if(istat)then
    do while(istat)
        write(*,*) '此文件夹已经存在，', &
        & '你确定把结果文件输入到这个文件夹吗？', &
        & '（是输入1，不是输入0）'
        read(*,*)b
        if(b>0)then
            exit
        else
            write(*,*)'请再次输入一个新的文件夹名：'
            read(*,*)foldername
            inquire(file='./'//trim(adjustl(foldername)),exist = istat)
            if(.not.istat)then
                call system('mkdir '//trim(adjustl(foldername)))
                exit
            endif
        endif
    end do
else
    call system('mkdir '//trim(adjustl(foldername)))
endif
end subroutine creatfolder
!#################################################
!################  完成文件夹的创建 #############
!#################################################
subroutine modelInformationOutput(a, fileDir, modelName)
! 该子程序将模型基本的信息输出
implicit none
type(geometry3d) :: a
character(len=100), intent(in) :: fileDir
character(len=100), intent(in) :: modelName
character(len=100) :: fileName
character(len=300) :: prefile
integer*4 :: i
prefile = trim(adjustl(filedir))//'./'//trim(adjustl(modelName))
! 输出坐标
fileName = 'Coordinate.bin'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)), &
& access='direct',form='unformatted',recl=a%pn*8)
do i = 1, a%ndim, 1
    write(10,rec=i) a%coor(:,i)
enddo
close(10)
! 输出体积
fileName = 'ParticleVolume.bin'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)), &
& access='direct',form='unformatted',recl=a%pn*8)
write(10,rec=1) a%vol(:,1)
close(10)
! 输出代表性边长
fileName = 'RepresentativeSideLength.bin'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)), &
& access='direct',form='unformatted',recl=a%pn*8)
write(10,rec=1) a%rpside(:,1)
close(10)
! 输出代表性半径
fileName = 'RepresentativeRadius.bin'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)), &
& access='direct',form='unformatted',recl=a%pn*8)
write(10,rec=1) a%rprad(:,1)
close(10)
! 输出hn
fileName = 'horizonParticleNumber.bin'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)), &
& access='direct',form='unformatted',recl=a%pn*4)
write(10,rec=1) a%hn(:)
close(10)
! 输出hs
fileName = 'HorizonParticleStartNumber.bin'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)), &
& access='direct',form='unformatted',recl=a%pn*4)
write(10,rec=1) a%hs(:)
close(10)
! 输出he
fileName = 'HorizonParticleEndNumber.bin'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)), &
& access='direct',form='unformatted',recl=a%pn*4)
write(10,rec=1) a%he(:)
close(10)
! 输出Horizon
fileName = 'Horizon.bin'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)), &
& access='direct',form='unformatted',recl=a%hpn*4)
write(10,rec=1) a%h(:)
close(10)
! 输出常数
fileName = 'GeometryMaterialConstant.dat'
open (10,file=trim(adjustl(prefile))//trim(adjustl(fileName)))
write(10,*)'length', a%length  ! 长度
write(10,*)'width', a%width    ! 宽度
write(10,*)'height', a%height  ! 高度
write(10,*)'ndim', a%ndim ! 维度数
write(10,*)'xn', a%xn     ! x方向划分质点数目
write(10,*)'yn', a%yn     ! y方向划分质点数目
write(10,*)'zn', a%zn     ! z方向划分质点数目
write(10,*)'xon', a%xon   ! x方向外部边界质点数目
write(10,*)'yon', a%yon   ! y方向外部边界质点数目
write(10,*)'zon', a%zon   ! z方向外部边界质点数目
write(10,*)'dx', a%dx     ! 物质点x方向划分的长度
write(10,*)'dy', a%dy     ! 物质点y方向划分的长度
write(10,*)'dz', a%dz     ! 物质点z方向划分的长度
write(10,*)'delta_x', a%delta_x     ! 网格尺寸
write(10,*)'delta', a%delta         ! 邻域半径
write(10,*)'mvalue', a%m            ! Mvalue
write(10,*)'particleNumber', a%pn   ! 总的物质点数目
write(10,*)'horPN', a%hpn  ! 总的邻域物质点数目
write(10,*)'YoungsModulus', a%emod
write(10,*)'ShearModulus', a%smod
write(10,*)'BulkModulus', a%bmod
write(10,*)'LameConstant', a%lameConst
write(10,*)'PoissonsRatio', a%pr
write(10,*)'SpringConstant', a%c
write(10,*)'MassDensity', a%massdens
write(10,*)'CriticalStretch', a%scr
write(10,*)'energyReleaseRate', a%energyReleaseRate
write(10,*)'dt', a%dt
close(10)
end subroutine modelInformationOutput

subroutine modelKineticTimeStepOutPut(a, fileDir, modelName, timeStep)
! 该子程序输出模型a的各个时程变量到fileDir里面去, modelName为自定义输出模型名字
! timeStep为输出时间步
implicit none
type(geometry3d), intent(in) :: a
character(len=100),intent(in) :: fileDir, ModelName
integer*4, intent(in) :: timestep
integer*4 :: i
character(len=100) :: ch, fileName
character(len=300) :: prefile
prefile = trim(adjustl(filedir))//'./'//trim(adjustl(modelName))
write(ch, '(I10)')timestep
! 输出位移
fileName = trim(adjustl(prefile))//'Displacement'// trim(adjustl(ch))//'.bin'
! EXAMPLE :: mydir/mymodelDisplacement100.txt
open (10,file=fileName, access='direct',form='unformatted',recl=a%pn*8)
do i = 1, a%ndim, 1
    write(10,rec=i) a%cdis(:,i)
enddo
close(10)
! 输出损伤
fileName = trim(adjustl(prefile))//'damage'// trim(adjustl(ch))//'.bin'
open (10,file=fileName, access='direct',form='unformatted',recl=a%pn*8)
write(10,rec=1) a%dmg(:,1)
close(10)
! 输出fail
fileName = trim(adjustl(prefile))//'fail'// trim(adjustl(ch))//'.bin'
open (10,file=fileName, access='direct',form='unformatted',recl=a%hpn*1)
write(10,rec=1) a%fail(:)
close(10)
end subroutine modelKineticTimeStepOutPut

subroutine modelRotate(a, rc, r) 
! 该子程序将模型a绕着旋转中心rc旋转角度r
implicit none
type(geometry3d), intent(inout) :: a
real*8, intent(in) :: rc(3), r(3)
real*8 :: rt(3)
real*8, allocatable :: rot3d(:,:), tempcoor(:,:)
allocate(rot3d(a%pn,a%ndim), tempcoor(a%pn,a%ndim))
rt = r * dacos(-1d0)/180d0
rot3d = a%coor
tempcoor = rot3d
rot3d(:,1) = (tempcoor(:,1)-rc(1))*cos(rt(2)) - &
& (tempcoor(:,3)-rc(3))*sin(rt(2)) + rc(1) ! y
rot3d(:,3) = (tempcoor(:,1)-rc(1))*sin(rt(2)) + &
& (tempcoor(:,3)-rc(3))*cos(rt(2)) + rc(3) ! y
tempcoor = rot3d
rot3d(:,2) = (tempcoor(:,2)-rc(2))*cos(rt(3)) - &
& (tempcoor(:,1)-rc(1))*sin(rt(3)) + rc(2) ! z
rot3d(:,1) = (tempcoor(:,2)-rc(2))*sin(rt(3)) + &
& (tempcoor(:,1)-rc(1))*cos(rt(3)) + rc(1) ! z
tempcoor = rot3d
rot3d(:,3) = (tempcoor(:,3)-rc(3))*cos(rt(1)) - &
& (tempcoor(:,2)-rc(2))*sin(rt(1)) + rc(3) ! x
rot3d(:,2) = (tempcoor(:,3)-rc(3))*sin(rt(1)) + &
& (tempcoor(:,2)-rc(2))*cos(rt(1)) + rc(2) ! x
a%coor = rot3d ! 对坐标进行更新
deallocate(tempcoor,rot3d)
end subroutine modelRotate

subroutine fileInputMaterialPoint(a, fileName, aspectRatio)
! 该函数通过文件名读入模型信息
! a -> geometry3d 的派生类
! fileName -> 文件名字, 文件为文本文件, 分为4列
! 文本中应有第一行标注物质点个数和列数4
! 以后每一行有坐标x,y,z和体积vol
! aspectRatio -> 几何缩放因子
implicit none
type(geometry3d), intent(inout) :: a
character(len=100), intent(in) :: fileName
real*8, intent(in) :: aspectRatio
real*8 :: evol
integer*4 :: i, j, k, np, ncol
open(11, file = fileName)
read(11,*) np, ncol
a%pn = np
allocate(a%coor(a%pn, a%ndim))
allocate(a%vol(a%pn, 1))
do i = 1, np, 1
    read(11,*)a%coor(i,1), a%coor(i,2), a%coor(i,3), a%vol(i, 1)
enddo
close(11)
a%coor = a%coor * aspectRatio
a%vol = a%vol * aspectRatio**(3d0)
evol = sum(a%vol(:,1))/a%pn
a%delta_x = evol**(1d0/3d0)
a%delta = a%m * a%delta_x
allocate ( a%rprad ( a%pn, 1 ) )
allocate ( a%rpside ( a%pn, 1 ) )
a%rprad = (a%vol*3d0/ ( 4d0 * dacos(-1d0)) )**(1d0/3d0)
a%rpside = a%vol**(1d0/3d0)
a%dx = a%delta_x
a%dy = a%delta_x
a%dz = a%delta_x
end subroutine fileInputMaterialPoint

subroutine irregSearchHorizon3d(a)
implicit none
type(geometry3d), intent(inout) :: a
integer*4,allocatable :: hrec(:), hrep(:)
integer*4 :: pcount
integer*4 :: maxhn, maxhl, hn, ct0, ct1
logical*1,allocatable :: ishor(:)
integer*4,allocatable :: nparr(:)
integer*4,allocatable :: ih(:), jh(:), hpi(:)
real*8, allocatable :: coor(:,:)
real*8 :: delta, dx, pi
integer*4 :: pNum, ndim, nh
integer*4 :: i
integer*4 :: nEnd, nStart
pi = dacos(-1d0)
delta = a%delta
dx = a%delta_x
pNum = a%pn
ndim = a%ndim
maxhn = nint(4.0d0/3.0d0*pi*(delta/dx+1d0)**3d0)
maxhl = maxhn*pNum
write(*,*)maxhn, maxhl
allocate(hrec(maxhl), ishor(pNum), nparr(pNum))
allocate(coor(pNum, ndim))
coor = a%coor
allocate( a%hn ( a%pn) )
allocate( a%hs ( a%pn) )
allocate( a%he ( a%pn) )
nparr = (/(i, i= 1, pNum, 1)/)
ct0=1; ct1=0
do i = 1, pNum, 1
    ishor = .false.
    where( dabs(coor(:,1)-coor(i,1))<=delta+0.5d0*dx &
    & .and.dabs(coor(:,2)-coor(i,2))<=delta+0.5d0*dx &
    & .and.dabs(coor(:,3)-coor(i,3))<=delta+0.5d0*dx )
        ishor = ( & 
        & (coor(:,1)-coor(i,1))*(coor(:,1)-coor(i,1)) + &
        & (coor(:,2)-coor(i,2))*(coor(:,2)-coor(i,2)) + &
        & (coor(:,3)-coor(i,3))*(coor(:,3)-coor(i,3))) <= &
        & (delta+0.5d0*dx)*(delta+0.5d0*dx)
    endwhere
    ishor(i) = .false.
    hn = count(ishor)
    a%hn(i) = hn ! 每个物质点邻域物质点数目
    ct1 =  ct1 + hn
    hrec(ct0:ct1:1) = pack(nparr, ishor)
    ct0 = ct1+1
enddo
nh = ct1
a%hpn = nh
allocate(a%h(a%hpn))
a%h(1:a%hpn:1) = hrec(1:a%hpn:1)
deallocate(hrec,ishor,nparr)
nEnd = 0
nStart = 1
do i = 1, a%pn, 1
    nEnd = nEnd + a%hn(i)
    a%hs(i) = nStart
    a%he(i) = nEnd
    nStart = nStart + a%hn(i)
enddo
end subroutine irregSearchHorizon3d
end module bondBasedPeridynamics3D

program PeriCrush_v10
use bondBasedPeridynamics3D
implicit none
integer*4 :: totalTimeStep
integer*4 :: i,j,k,m,n,p,tt ! 循环专用标号
real*8 :: t1, t2
real*8, allocatable :: histContactForce(:,:)
! 物质点之间的连接状态! 逻辑变量与逻辑数组
logical*1, allocatable :: isParticleConnected(:) 
character*100 :: inputFileName ! 输入文件名! 符号变量
character*100 :: outputFileDirectory ! 输出文件目录
character*100 :: outPutFileName ! 输出文件名
type(geometry3d) :: a ! 基本body
type(geometry3d) :: b(3) ! 实际计算body
! 接触结构体[包含接触对接触的物理信息、接触质点信息]
type(localContact) :: conta(2)
real*8 :: appliedVelocity
character(len=100) :: fileDir
character(len=100) :: modelName
character(len=100) :: fileName
real*8, allocatable :: arrTemp(:), alphaRatio(:)
logical*1, allocatable :: lgcTemp(:)
real*8 :: roty
integer*4 :: loadPlatenDivNumber, SphereLengthDivNumber
!############## 模型信息 ###########
roty = 45d0  ! 
SphereLengthDivNumber = 43
loadPlatenDivNumber   = 13
a%m = 3.015d0
!############## 模型信息 ###########
! 几何信息 [Object b(1)]
a%length = 1.5d-3
a%width  = 1.5d-3
a%height = 1.5d-3
a%xn = SphereLengthDivNumber
a%yn = SphereLengthDivNumber
a%zn = SphereLengthDivNumber
a%xon = floor(a%m)
a%yon = floor(a%m)
a%zon = floor(a%m)
a%ndim = 3
call tetrahedronGrids(a)
call regHorizonSearch(a)
b(1) = particleDelete(a, sum(a%coor**2, dim=2) .ge. 0.25d0*(a%length + a%delta_x)**2) ! 得到第一个物体
allocate(arrTemp(b(1)%pn), lgcTemp(b(1)%pn), alphaRatio(b(1)%pn))
arrTemp = dsqrt(sum(b(1)%coor**2, dim=2))
lgcTemp =  arrTemp.ge.0.5d0*(b(1)%length - b(1)%delta_x)
alphaRatio = ( 0.5d0*b(1)%length - arrTemp + 0.5d0*b(1)%delta_x)/b(1)%delta_x
where (lgcTemp)
    b(1)%vol(:,1)  = b(1)%vol(:,1) * alphaRatio
    b(1)%coor(:,1) = b(1)%coor(:,1)/arrTemp * ( 0.5d0*b(1)%length - 0.5d0 * alphaRatio * b(1)%delta_x)
    b(1)%coor(:,2) = b(1)%coor(:,2)/arrTemp * ( 0.5d0*b(1)%length - 0.5d0 * alphaRatio * b(1)%delta_x)
    b(1)%coor(:,3) = b(1)%coor(:,3)/arrTemp * ( 0.5d0*b(1)%length - 0.5d0 * alphaRatio * b(1)%delta_x)
endwhere
b(1)%rpSide = b(1)%vol**(1d0/3d0)                        ! 代表性长度
b(1)%rpRad = (b(1)%vol/(4d0/3d0*dacos(-1d0)))**(1d0/3d0) ! 代表性长度
deallocate(arrTemp, lgcTemp, alphaRatio)
call modelRotate(b(1), (/0d0,0d0,0d0/),(/0d0,roty,0d0/)) ! 将模型绕y轴旋转一个角度
call modelDelete(a)                                      ! 删除模型所有信息
call showModelGeometryInfo(b(1))                         ! 显示模型信息
call kineticVariableAllocate(b(1))                       ! 给运动量分配空间
! 给出模型参数
b(1)%emod = 1.00d11
b(1)%pr   = 0.25d0
b(1)%massdens = 2650d0
b(1)%dt       = 0.1d0 * b(1)%delta_x/dsqrt(b(1)%emod/b(1)%massDens)
call modelElasticConstant(b(1))
call surfaceCorrection(b(1))
b(1)%idis = 0d0
b(1)%iacc = 0d0
b(1)%ivel = 0d0
b(1)%cdis = 0d0
b(1)%cfdens = 0d0
b(1)%energyReleaseRate = 10d0
b(1)%scr = dsqrt(5d0*b(1)%energyReleaseRate/(9d0*b(1)%bmod*b(1)%delta))
! 几何信息[Object b(2)]
a%length = 5d0/15d0*b(1)%length
a%width  = 5d0/15d0*b(1)%width
a%height = 5d0/15d0*b(1)%height
a%xn = loadPlatenDivNumber
a%yn = loadPlatenDivNumber
a%zn = loadPlatenDivNumber
a%m = 3.015d0
a%xon = 0
a%yon = 0
a%zon = 0
a%ndim = 3
call tetrahedronGrids(a) ! 生成规则网格
call regHorizonSearch(a) ! 得到规则邻域
b(2) = particleDelete(a, abs(a%coor(:,3)).ge.100d0)     ! 得到第一个物体
b(2)%rpSide = b(2)%vol**(1d0/3d0)                        ! 代表性长度
b(2)%rpRad = (b(2)%vol/(4d0/3d0*dacos(-1d0)))**(1d0/3d0) ! 代表性长度
call modelDelete(a)                ! 删除模型所有信息
call showModelGeometryInfo(b(2))   ! 显示模型信息
call kineticVariableAllocate(b(2)) ! 给运动量分配空间
! 给出模型参数
b(2)%emod = 380d9
b(2)%pr   = 0.25d0
b(2)%massdens = 3850d0
b(2)%dt = b(1)%dt
call modelElasticConstant(b(2))
call surfaceCorrection(b(2))
b(2)%idis = 0d0
b(2)%iacc = 0d0
b(2)%ivel = 0d0
b(2)%cdis = 0d0
b(2)%cfdens = 0d0
b(2)%energyReleaseRate = 1000000d0
b(2)%scr = dsqrt(5d0*b(2)%energyReleaseRate/(9d0*b(2)%bmod*b(2)%delta))
call modelTranslation(b(2), (/0d0, 0d0, 0.5d0 * (b(1)%height + &
& b(2)%height ) + 0.5d0*sum(b(2)%rpRad(:,1))/b(2)%pn /))
! 几何信息[Object b(3)]
a%length = 5d0/15d0*b(1)%length
a%width  = 5d0/15d0*b(1)%width
a%height = 5d0/15d0*b(1)%height
a%xn = loadPlatenDivNumber
a%yn = loadPlatenDivNumber
a%zn = loadPlatenDivNumber
a%m = 3.015d0
a%xon = 0
a%yon = 0
a%zon = 0
a%ndim = 3
call tetrahedronGrids(a) ! 生成规则网格
call regHorizonSearch(a) ! 得到规则邻域
b(3) = particleDelete(a, abs(a%coor(:,3)) .ge. 100d0)    ! 得到第一个物体
b(3)%rpSide = b(3)%vol**(1d0/3d0) ! 代表性长度
b(3)%rpRad = (b(3)%vol/(4d0/3d0*dacos(-1d0)))**(1d0/3d0) ! 代表性长度
call modelDelete(a)                ! 删除模型所有信息
call showModelGeometryInfo(b(3))   ! 显示模型信息
call kineticVariableAllocate(b(3)) ! 给运动量分配空间
! 给出模型参数
b(3)%emod = 380d9
b(3)%pr = 0.25d0
b(3)%massdens = 3850d0
b(3)%dt = b(1)%dt
call modelElasticConstant(b(3))
call surfaceCorrection(b(3))
write(*, *)b(3)%pn
write(*, *)sum(b(3)%cc(:,1))/b(3)%pn
write(*, *)minval(b(3)%cc(:,1))
write(*, *)minval(b(3)%vol(:,1))
b(3)%idis = 0d0
b(3)%iacc = 0d0
b(3)%ivel = 0d0
b(3)%cdis = 0d0
b(3)%cfdens = 0d0
b(3)%energyReleaseRate = 1000000d0
b(3)%scr = dsqrt(5d0*b(3)%energyReleaseRate/(9d0*b(3)%bmod*b(3)%delta))
call modelTranslation(b(3), (/0d0, 0d0, - 0.5d0 * (b(1)%height + &
& b(3)%height ) - 0.5d0*sum(b(3)%rpRad(:,1))/b(3)%pn /))
call contactPairs(b(1), b(2), &
& b(1)%coor(:,3).ge. b(1)%height/2d0 - b(1)%height/5d0, &
& b(2)%coor(:,3).le. (b(1)%height/2d0 + b(2)%height/3d0), &
& b(1)%height/5d0, conta(1))
! 获取接触对
call contactPairs(b(1), b(3), &
& b(1)%coor(:,3).le. -b(1)%height/2d0 + b(1)%height/5d0, &
& b(3)%coor(:,3).ge. -(b(1)%height/2d0 + b(3)%height/3d0), &
& b(1)%height/5d0, conta(2))
conta(1)%mu = 0.0d0 ! 接触静摩擦系数
conta(2)%mu = 0.0d0 ! 接触静摩擦系数
write(*, *)'Contact 1 ContactPairNumber:',conta(1)%ncp
write(*, *)'Contact 2 ContactPairNumber:',conta(2)%ncp
totalTimeStep = 2000000
allocate(histContactForce(totalTimeStep, 4))
! 重新校订临界时间步长
b(1)%dt = 0.5d0 * minval((/ & 
& minval(b(1)%rpSide(:,1))/dsqrt(b(1)%emod/b(1)%massdens),&
& minval(b(2)%rpSide(:,1))/dsqrt(b(2)%emod/b(2)%massdens),&
& minval(b(3)%rpSide(:,1))/dsqrt(b(3)%emod/b(3)%massdens),&
& dacos(-1d0)*minval(b(1)%rpRad(:,1))* &
& dsqrt(b(1)%massdens/b(1)%smod)/(0.887d0 + 0.163d0 * b(1)%pr),&
& dacos(-1d0)*minval(b(2)%rpRad(:,1))* &
& dsqrt(b(2)%massdens/b(2)%smod)/(0.887d0 + 0.163d0 * b(2)%pr),&
& dacos(-1d0)*minval(b(3)%rpRad(:,1))* &
& dsqrt(b(3)%massdens/b(3)%smod)/(0.887d0 + 0.163d0 * b(3)%pr)/))
b(2)%dt = b(1)%dt
b(3)%dt = b(1)%dt
! 初始化-1步的位移
b(1)%pdis = b(1)%idis + 0.5d0*b(1)%iacc*b(1)%dt*b(1)%dt - b(1)%ivel * b(1)%dt
b(2)%pdis = b(2)%idis + 0.5d0*b(2)%iacc*b(2)%dt*b(2)%dt - b(2)%ivel * b(2)%dt
b(3)%pdis = b(3)%idis + 0.5d0*b(3)%iacc*b(3)%dt*b(3)%dt - b(3)%ivel * b(3)%dt
write(*,*)'临界时间步长=',b(1)%dt,'秒'
! 模型数据输出
call creatfolder(fileDir)
modelName = 'S1'
call modelInformationOutput(b(1), fileDir, modelName)
modelName = 'C1'
call modelInformationOutput(b(2), fileDir, modelName)
modelName = 'C2'
call modelInformationOutput(b(3), fileDir, modelName)
appliedVelocity = 0.1d0
do tt = 1, totalTimeStep, 1
    write(*,*)'Current Time Step:', tt
    write(*,*)'Current Time:',tt*b(1)%dt*1d6,' us'
    b(2)%cdis(:,3) = - appliedVelocity * tt* b(1)%dt
    b(3)%cdis(:,3) = appliedVelocity * tt* b(1)%dt
    histContactForce(tt,4) = appliedVelocity * tt* b(1)%dt
    if(mod(tt, 200) .eq. 0) then
        ! 模型数据输出
        modelName = 'S1'
        call modelKineticTimeStepOutPut(b(1), fileDir, modelName, tt)
        modelName = 'C1'
        call modelKineticTimeStepOutPut(b(2), fileDir, modelName, tt)
        modelName = 'C2'
        call modelKineticTimeStepOutPut(b(3), fileDir, modelName, tt)
    endif
    call bondBasedForceDensityVector(b(1)) ! 第一个模型的力
    call bondBasedForceDensityVector(b(2)) ! 第二个模型的力
    call bondBasedForceDensityVector(b(3)) ! 第三个模型的力
    call contactForceDensity(b(1), b(2), conta(1)) ! 加入接触力密度部分
    call contactForceDensity(b(1), b(3), conta(2)) ! 加入接触力密度部分
    call timeIntegration(b(1))																				
    call timeIntegration(b(2))
    call timeIntegration(b(3))
    histContactForce(tt,1:3:1) = conta(1)%tcf(1:3:1)
    if(mod(tt,1000).eq.0)then
        fileName = trim(adjustl(filedir))//'./'//'ContactForce.txt'
        open(10, file = fileName)
        do i = 1, tt, 100
            write(10,*) histContactForce(i,:)
        enddo
        close(10)
        ! 假如力已经很小了，就结束了
        if(dabs(histContactForce(tt,3))/maxval(dabs(histContactForce(:,3))).le.1d-3 & 
        &.and.maxval(dabs(histContactForce(:,3))).gt.10d0.and.tt.gt.20000) exit 
    endif
enddo
write(*,*)'Program Finished'
read(*,*) ! 如果任务停止可以看到任务输出详情
end program PeriCrush_v10