!#######################################################################

      module Hf_Geom

!#######################################################################
!> @param dhc     翼片圆弧直径  [in]/[m]  (HCF)
!> @param phc     螺距 [in]/[m] (HCF)
!> @param lhc     翼片直段长度 [in]/[m] (HCF)
!> @param rhc     翼片内弧半径 [in]/[m] (HCF)
!> @param dmhc    截面最大直径 [in]/[m] (HCF)
!> @param thick   clad thick [in]/[m] (HCF)
!> @paarm facnonu Maximum circumferential heat flux factor
        implicit none
        real :: pitch,roughness
        real :: dhc,phc,lhc,rhc,dmhc
        real :: thick
        real :: facnonu
        integer :: hf_c,hcf
        !integer :: TotMeshNum
        !integer :: FuelMeshNum,CladMeshNum
        
        !!!!!Helical fuel heat conduction
        integer :: matfuel,matclad                 !material index of fuel and clad
        integer :: itermax=30000
        integer :: N1,N2,N3,N4,N5,N6                      !网格划分参数
        real :: ntcp                                      !Number of thermal conduction path
        logical :: HF_first_iter
        
        real,allocatable,dimension (:,:) :: nodeCoordX1_clad,nodeCoordY1_clad  !block 1_clad
        real,allocatable,dimension (:,:) :: nodeCoordX2_clad,nodeCoordY2_clad  !block 2_clad
        real,allocatable,dimension (:,:) :: nodeCoordX3_clad,nodeCoordY3_clad  !block 3_clad
        real,allocatable,dimension (:,:) :: nodeCoordX4_clad,nodeCoordY4_clad  !block 4_clad
        real,allocatable,dimension (:,:) :: nodeCoordX1_fuel,nodeCoordY1_fuel  !block 1_fuel
        real,allocatable,dimension (:,:) :: nodeCoordX2_fuel,nodeCoordY2_fuel  !block 2_fuel
        real,allocatable,dimension (:,:) :: nodeCoordX3_fuel,nodeCoordY3_fuel  !block 3_fuel
        real,allocatable,dimension (:,:) :: nodeCoordX4_fuel,nodeCoordY4_fuel  !block 4_fuel
        real,allocatable,dimension (:,:) :: nodeCoordX5_fuel,nodeCoordY5_fuel  !block 5_fuel
        real,allocatable,dimension (:,:) :: nodeCoordX6_fuel,nodeCoordY6_fuel  !block 6_fuel
        
        real,allocatable,dimension (:,:) :: nodeCoordX1_2clad,nodeCoordY1_2clad  !block 1_2clad
        real,allocatable,dimension (:,:) :: nodeCoordX2_2clad,nodeCoordY2_2clad  !block 2_2clad
        real,allocatable,dimension (:,:) :: nodeCoordX3_2clad,nodeCoordY3_2clad  !block 3_2clad
        real,allocatable,dimension (:,:) :: nodeCoordX4_2clad,nodeCoordY4_2clad  !block 4_2clad
        real,allocatable,dimension (:,:) :: nodeCoordX1_2fuel,nodeCoordY1_2fuel  !block 1_2fuel
        real,allocatable,dimension (:,:) :: nodeCoordX2_2fuel,nodeCoordY2_2fuel  !block 2_2fuel
        real,allocatable,dimension (:,:) :: nodeCoordX3_2fuel,nodeCoordY3_2fuel  !block 3_2fuel
        real,allocatable,dimension (:,:) :: nodeCoordX4_2fuel,nodeCoordY4_2fuel  !block 4_2fuel
        real,allocatable,dimension (:,:) :: nodeCoordX5_2fuel,nodeCoordY5_2fuel  !block 5_2fuel
        real,allocatable,dimension (:,:) :: nodeCoordX6_2fuel,nodeCoordY6_2fuel  !block 6_2fuel


        real,allocatable,dimension (:) :: x1_clad,y1_clad,x1_2clad,y1_2clad
        real,allocatable,dimension (:) :: x2_clad,y2_clad,x2_2clad,y2_2clad
        real,allocatable,dimension (:) :: x3_clad,y3_clad,x3_2clad,y3_2clad
        real,allocatable,dimension (:) :: x4_clad,y4_clad,x4_2clad,y4_2clad
        real,allocatable,dimension (:) :: x5_clad,y5_clad,x5_2clad,y5_2clad
        real,allocatable,dimension (:) :: x6_clad,y6_clad,x6_2clad,y6_2clad
        real,allocatable,dimension (:) :: x7_clad,y7_clad,x7_2clad,y7_2clad
        real,allocatable,dimension (:) :: x8_clad,y8_clad,x8_2clad,y8_2clad
      
        real,allocatable,dimension (:) :: x1_fuel,y1_fuel,x1_2fuel,y1_2fuel
        real,allocatable,dimension (:) :: x2_fuel,y2_fuel,x2_2fuel,y2_2fuel
        real,allocatable,dimension (:) :: x3_fuel,y3_fuel,x3_2fuel,y3_2fuel
        real,allocatable,dimension (:) :: x4_fuel,y4_fuel,x4_2fuel,y4_2fuel
        real,allocatable,dimension (:) :: x5_fuel,y5_fuel,x5_2fuel,y5_2fuel
        real,allocatable,dimension (:) :: x6_fuel,y6_fuel,x6_2fuel,y6_2fuel
        real,allocatable,dimension (:) :: x7_fuel,y7_fuel,x7_2fuel,y7_2fuel
        real,allocatable,dimension (:) :: x8_fuel,y8_fuel,x8_2fuel,y8_2fuel
        real,allocatable,dimension (:) :: x9_fuel,y9_fuel,x9_2fuel,y9_2fuel
        real,allocatable,dimension (:) :: x10_fuel,y10_fuel,x10_2fuel,y10_2fuel
        
       
        type :: block
            !>控制体左上角坐标
            real,allocatable :: wn_x(:,:),wn_y(:,:)
            !>控制体左下角坐标
            real,allocatable :: ws_x(:,:),ws_y(:,:)            
            !>控制体右上角坐标
            real,allocatable :: en_x(:,:),en_y(:,:)             
            !>控制体右下角坐标
            real,allocatable :: es_x(:,:),es_y(:,:)
            !>控制体中心坐标
            real,allocatable :: c_x(:,:),c_y(:,:)
            !>控制体上面坐标
            real,allocatable :: n_x(:,:),n_y(:,:)
            !>控制体下面坐标
            real,allocatable :: s_x(:,:),s_y(:,:)
            !>控制体左面坐标
            real,allocatable :: w_x(:,:),w_y(:,:)
            !>控制体右面坐标
            real,allocatable :: e_x(:,:),e_y(:,:)
            !>控制体上、下、左、右边长度
            real,allocatable :: LsN(:,:),LsS(:,:),LsW(:,:),LsE(:,:)
            !> 控制体中心到上、下、左、右边距离
            real,allocatable :: DsN(:,:),DsS(:,:),DsW(:,:),DsE(:,:)
            !> 控制体中心到相邻控制体中心距离
            real,allocatable :: Ds2N(:,:),Ds2S(:,:),Ds2W(:,:),Ds2E(:,:)
            !> 控制体dx,dy
            real,allocatable :: deltaX(:,:),deltaY(:,:)
         !> Angle of opposite side center line
            real,allocatable :: in_ang_sin(:,:)
            !>Area of control volume
            real,allocatable :: area(:,:)
            ! >控制体编号
            integer,allocatable :: mesh_num(:,:)
            !>控制体相邻控制体编号
            integer,allocatable :: neighb_n(:,:),neighb_s(:,:),neighb_w(:,:),neighb_e(:,:)        
        end type
                 
        real,allocatable,dimension (:,:) :: res           !residual error
        
     type, public :: HFCond
        real,allocatable,dimension (:,:,:) :: fblk1_T,fblk2_T,fblk3_T,fblk4_T,fblk5_T,fblk6_T                !Fuel temperature,region 1
        real,allocatable,dimension (:,:,:) :: fblk1_2T,fblk2_2T,fblk3_2T,fblk4_2T,fblk5_2T,fblk6_2T          !Fuel temperature,region 2       
        real,allocatable,dimension (:,:,:) :: fblk1_Tn,fblk2_Tn,fblk3_Tn,fblk4_Tn,fblk5_Tn,fblk6_Tn          !Fuel temperature in old time step, region 1
        real,allocatable,dimension (:,:,:) :: fblk1_2Tn,fblk2_2Tn,fblk3_2Tn,fblk4_2Tn,fblk5_2Tn,fblk6_2Tn    !Fuel temperature in old time step, region 2        
        real,allocatable,dimension (:,:,:) :: cblk1_T,cblk2_T,cblk3_T,cblk4_T                                !Cladding temperature,region 1
        real,allocatable,dimension (:,:,:) :: cblk1_2T,cblk2_2T,cblk3_2T,cblk4_2T                            !Cladding temperature,region 2        
        real,allocatable,dimension (:,:,:) :: cblk1_Tn,cblk2_Tn,cblk3_Tn,cblk4_Tn                            !Cladding temperature in old time step.region 1
        real,allocatable,dimension (:,:,:) :: cblk1_2Tn,cblk2_2Tn,cblk3_2Tn,cblk4_2Tn                            !Cladding temperature in old time step.region 2        
        real,allocatable,dimension (:,:)   :: twall                                                          !Cladding surface temperature
        real,allocatable,dimension (:,:)   :: heatFlux                                                       !Cladding surface heat flux
        real :: fuel_source                                                                            !Fuel source term
        real :: cladding_source                                                                        !Cladding source term
     end type HFCond 
     
     type(block) :: blk1_fuel,blk2_fuel,blk3_fuel,blk4_fuel,blk5_fuel,blk6_fuel
     type(block) :: blk1_clad,blk2_clad,blk3_clad,blk4_clad  
     type(block) :: blk1_2fuel,blk2_2fuel,blk3_2fuel,blk4_2fuel,blk5_2fuel,blk6_2fuel
     type(block) :: blk1_2clad,blk2_2clad,blk3_2clad,blk4_2clad                

     end module Hf_Geom
