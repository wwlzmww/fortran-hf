!#######################################################################

     module HF_init
     use Hf_Geom,  only: HFcond
     use Rodtab,     only: nrrod
          
     implicit none

     type(HFCond    ), allocatable, public, target :: HFrods(:)
 
     end module HF_init
