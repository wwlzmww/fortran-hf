!#######################################################################

     module HF_init
     use Hf_Geom,  only: HFcond
     !use Rodtab,     only: nrrod
          
     implicit none

     type(HFCond    ), allocatable, public, target :: HFrods(:) !没有实例化

     integer :: nrrod = 1
     real :: dt = 1e-3

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


     end module HF_init
