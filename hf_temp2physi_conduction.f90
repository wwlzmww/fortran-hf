!!!! in:温度场：HFrod%fblk1_T(:,:,:)-->HFrod%fblk6_T(:,:,:)
!!!! in:温度场：HFrod%cblk1_T(:,:,:)-->HFrod%cblk4_T(:,:,:)
!!!! in:网格
!!!! out:位移场
module temp2physi_conduction
    use hf_heat_conduction
    implicit none

    !Zr弹性模量关于温度的函数值
    real function Elasticmo(T)  !T单位K
        implicit none
        real,intent(in) :: T
        Elasticmo = 921 - 4.05e-2*T  !弹性模量，单位GPa
        return
    end

    !Zr泊松比关于温度的函数值
    real function poission(T)  !单位K
    use Elasticmo, only: Elasticmo
        implicit none
        real,intent(in) :: T
        G = 349 - 1.66e-2*T  !剪切模量，单位GPa
        E = Elasticmo(T)
        poission = E/2/G - 1
        return
    end

    !只需要知道单元编号，对应结点坐标与关联单元，对应温度，结点即可，形式可以转换
    !现假设以获得的单元控制体的标号为ele,对应的温度为ele_t
    !ele的形式应该为:1 1 2 3 4
    !               2 1 2 5 6
    !其中第一列表示单元编号，后四位数字为单元对应角结点编号，ele_t同理，第一列表示单元编号，后一位数字表示单元温度

    !先对每个单元的力学量进行计算
    real dimension(ele_row) :: Elmo
    real dimension(ele_row) :: poi
    forall(i = 1:ele_row)
        Elmo(i) = Elasticmo(ele_t(i,2))
        poi(i) = poission(ele_t(i,2)) 
    end forall

    !进行单元刚度矩阵计算，采取四结点四边形计算
    subroutine stiffness(E,nu,ele,node,p) !E为弹性模量，nu为泊松比
    !node为结点编号以及坐标，第一列为结点编号，第二列为x坐标，第三列为y坐标，p为压力
        implicit none
        real, intent(in) :: E(:),nu(:),ele(:,:),node(:,:),p !还未指定p的大小
        real dimension(ele_row), intent(out) :: k_stiffness
        real, parameter :: const = 1/sqrt(3), pi = 3.141592654
        do i = 1:ele_row
            !定义结点真实坐标，设单元结点编号右上角为1，右下为2，左下为3，左上为4
            !同时刚度矩阵积分采取高斯近似，选取正负根号3分之1的4个样本点，权重为1
            real dimension(4) :: xreal, yreal,
            real dimension(4,2) :: sample
                forall(j = 1:4)
                xreal(j) = node(ele(i,j+1),2)
                yreal(j) = node(ele(i,j+1),3)
                sample(j,1) = sin((j/2-0.25)*pi) / abs(sin((j/2-0.25)*pi)) * const
                sample(j,2) = sin((-j/2+0.75)*pi) / abs(sin((-j/2+0.75)*pi)) * const
                end forall
            
            do j = 1:4
            
            a = 0.25*( (1+sample(j,2))*xreal(1) + (1-sample(j,2))*xreal(2) + (sample(j,2)-1)*xreal(3) - (sample(j,2)+1)*xreal(4))

        end do

    end subroutine stiffness

end module temp2physi_conduction