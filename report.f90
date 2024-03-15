
!!!!!! 简单反推，但是不确定的
! Fuelrod 存疑。直接del了，实例化指针rod，但是没有指向内容
! Get_time_data 存疑，似乎只影响功率取定，直接删除了
! liquid_props 存疑，疑似是输出(h = havg_ch,p=pij,T=tfluid)，设为常值即可，但单位含义未知
! material_props 存疑，疑似是输出(imat=matfuel,t=tref,cp=cp,k=kc)，设为常值即可，但单位含义未知
! unitf 存疑，对这个库已知的只有psfrel2psia函数，但是有非only的use
module type_add
! rods 反推并实例化
type :: rods_type
    integer :: jmax,kmax
    real :: linear_power
    real,allocatable,dimension (:) :: x
    type(surf_type),allocatable,dimension(:) :: surf
end type
type :: surf_type
    real,allocatable,dimension(:) :: htcl
end type

type(rods_type) :: rods(1)
end module type_add
!!!!!! 有很大问题的
! pin_sc_conn 完全不知道作用，注释掉了,疑似是得到一个换热相关参数chnum
! flmesh 流体网格类的实例，只出现了一次，反推不出来定义
! ch 疑似是从chnum得到换热系数，返回值似乎是类

!!!!!! 大概没什么问题的
! ihfout 读写通道号，大概取定常数>10即可
! Psfrel2psia 猜测【从磅力每平方英尺psf转换为磅力每平方英寸psia】
program report
    implicit none
    real,external :: Psfrel2psia
end program report
function Psfrel2psia(psf)
    implicit none
    real,intent(in) :: psf
    real :: Psfrel2psia
    Psfrel2psia = psf/144.0
end function

! Rcold 输出RFuel棒密度，设为常数
! nrrod 棒数量 设为1
integer :: nrrod = 1
! dt 时间步长，单位未知，计算中/3600，可能为小时
real :: dt = 1e-3

