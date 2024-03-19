!!!! in:温度场：HFrod%fblk1_T(:,:,:)-->HFrod%fblk6_T(:,:,:)
!!!! in:温度场：HFrod%cblk1_T(:,:,:)-->HFrod%cblk4_T(:,:,:)
!!!! in:网格
!!!! out:位移场
module temp2physi_conduction
    use hf_heat_conduction
    implicit none

    !Zr与UO2弹性模量关于温度的函数值
    real function Elasticmo(T,type)  !Elasticmo单位GPa , T单位K , type表示材料类型
        implicit none
        real, intent(in) :: T
        integer, intent(in) :: type
            select case (type)
            !type = 1 ：Zr
            !type = 2 : UO2
            case (1)
                Elasticmo = 921 - 4.05e-2*T  
            case (2)
                Elasticmo = 233.4 - 2.5476e-2*T 
            case default
                write(*,*) "未定义材料"
            end select 
        return
    end

    !Zr与UO2泊松比关于温度的函数值
    real function poission(T,type)  !T单位K , type表示材料类型
        use Elasticmo, only: Elasticmo
        implicit none
        real, intent(in) :: T
        integer, intent(in) :: type
            select case (type)
            !type = 1 ：Zr
            !type = 2 : UO2
            case (1)
                G = 349 - 1.66e-2*T  !剪切模量，单位GPa
                E = Elasticmo(T,1)
                poission = E/2/G - 1
            case (2)
                poission = 0.328
            case default
                write(*,*) "未定义材料"
            end select 
        return
    end

    !只需要知道单元编号，对应结点坐标与关联单元，对应温度，结点即可，形式可以转换
    !现假设以获得的单元控制体的标号矩阵为ele,对应的温度为ele_t
    !ele的形式应该为:1 1 2 3 4
    !               2 1 2 5 6
    !其中第一列表示单元编号，后四位数字为单元对应角结点编号，ele_t同理，第一列表示单元编号，后一位数字表示单元温度
    !先对每个单元的力学量进行计算

    real dimension(ele_row) :: Elmo
    real dimension(ele_row) :: poi
    forall(i = 1:ele_uo2)
    !编号1~ele_uo2表示材料为UO2的单元编号，ele_uo2~ele_row表示材料为zr的单元编号
        Elmo(i) = Elasticmo(ele_t(i,2),2)
        Elmo(ele_row + 1 - i) = Elasticmo(ele_t(ele_row + 1 - i,2),1)
        poi(i) = poission(ele_t(i,2),2) 
        poi(ele_row + 1 - i) = poission(ele_t(ele_row + 1 - i,2),1) 
    end forall

    !进行单元刚度矩阵计算，采取四结点四边形计算
    subroutine stiffness(E,nu,ele,node,type)
    !E为弹性模量，nu为泊松比.node为结点编号以及坐标，第一列为结点编号，第二列为x坐标，第三列为y坐标，type为类型
        implicit none
        real, intent(in) :: E(:),nu(:),ele(:,:),node(:,:)
        real dimension(8,8,ele_row), intent(out) :: k_stiffness
        real, parameter :: con = 1/sqrt(3)
        real dimension(4,1) :: xreal, yreal
        real dimension(4,2) :: sample
        real dimension(2,4) :: temp
        real dimension(8,8) :: tempk = 0
        real dimension(3,8) :: strainmatrix
        real dimension(3,3) :: stressmatrix
        real dimension(8) :: formcoff
        real :: a, b, c, d, coffA, J
        do i = 1:ele_row
            !定义结点真实坐标，设单元结点编号右上角为1，右下为2，左下为3，左上为4
            !同时刚度矩阵积分采取高斯近似，选取正负根号3分之1的4个样本点，权重为1
            sample = reshape([con,con,-con,-con,con,-con,-con,con],[4,2])
            forall(j = 1:4)
            xreal(j,1) = node(ele(i,j+1),2)
            yreal(j,1) = node(ele(i,j+1),3)
            end forall

            !以下计算第i个单元的刚度矩阵
            do j = 1:4
                temp = reshape([ 1+sample(j,2) ,  1+sample(j,1) , -1+sample(j,2) , -1-sample(j,1) ,&
                              & -1+sample(j,2) , -1+sample(j,1) , -1-sample(j,2) ,  1-sample(j,1)],[2,4])
                a = 0.25*matmul(temp,xreal)(1,1)
                b = 0.25*matmul(temp,xreal)(2,1)
                c = 0.25*matmul(temp,yreal)(1,1)
                d = 0.25*matmul(temp,yreal)(2,1)
                !计算雅各布矩阵行列式
                J = a*d - b*C
                !计算形函数系数矩阵
                formcoff = 0.25*[-1+sample(j,2) , -1+sample(j,1) ,  1-sample(j,2) , -1-sample(j,1) ,&!是否可以直接用数乘一个数组？
                           & 1+sample(j,2) ,  1+sample(j,1) , -1-sample(j,2) ,  1-sample(j,1)]
                    !计算应变矩阵strainmatrix
                    do w = 1:4
                        strainmatrix(1,2*w-1) = 0.25*(d*formcoff(2*w-1) - c*formcoff(2*w))/J
                        strainmatrix(2,2*w)   = 0.25*(a*formcoff(2*w) - b*formcoff(2*w-1))/J
                        strainmatrix(3,2*w-1) = 0.25*(a*formcoff(2*w) - b*formcoff(2*w-1))/J
                        strainmatrix(3,2*w)   = 0.25*(d*formcoff(2*w-1) - c*formcoff(2*w))/J
                    end do
                    !计算应力矩阵stressmatrix
                    select case(type)
                    !type = 1:平面应力
                    !type = 2:平面应变
                    case(1)
                        coffA = E(i)/(1-nu(i)**2)
                        stressmatrix = reshape([coffA,coffA*nu(i),0,coffA*nu(i),coffA,0,0,0,coffA*(1-nu(i))/2],[3,3])
                    case(2)
                        coffA = E(i)/(1+nu(i))/(1-2*nu(i))
                        stressmatrix = reshape([coffA*(1-nu(i)),coffA*nu(i),0,coffA*nu(i),coffA*(1-nu(i)),0,0,0,coffA*(1-2*nu(i))/2],[3,3])
                    case default
                        write(*,*) "unknown"
                    end select
                    
                    tempk = tempk + J* matmul(matmul(transpose(strainmatrix),stressmatrix),strainmatrix)
            end do 
        k_stiffness(:,:,i) = tempk
        end do
    end subroutine stiffness

end module temp2physi_conduction