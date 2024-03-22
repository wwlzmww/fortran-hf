!!!! in:温度场：HFrod%fblk1_T(:,:,:)-->HFrod%fblk6_T(:,:,:)
!!!! in:温度场：HFrod%cblk1_T(:,:,:)-->HFrod%cblk4_T(:,:,:)
!!!! in:网格
!!!! out:位移场
module mathcal
    implicit none

    !定义数与一维数组的数乘
    real function multnV(a,V)
        implicit none
        integer :: i
        real, intent(in) :: a, V(:)
        real, dimension(size(V)), intent(out) :: multnV 
        forall(i = 1 : size(V))
            multnV(i) = a*V(i)
        end forall
        return
    end
        
    !定义数与矩阵的数乘
    real function multnM(a,M)
        implicit none
        integer :: i, j
        real, intent(in) :: a,M(:,:)
        real, dimension(size(M,1) , size(M,2)), intent(out) :: multnM 
        forall(i = 1 : size(M,1) , j = 1 : size(M,2))
            multnM(i,j) = a*M(i,j)
        end forall
        return
    end

    !返回曲线侧的x与y方向分量,方向由对应点指向曲率圆心
    real function cirlencom(p,x,y,cirx,cirx)
        !p为需要转化的量的长度，x，y为对应点X,Y坐标，cirx,ciry为曲率圆心坐标
        implicit none
        real, intent(in) :: p, x, y, cirx, ciry
        real, dimension(2), intent(out) :: cirlencom !第一列表示x分量，第二列表示y分量
        real :: radius
            radius = sqrt( (cirx - x)**2 + (ciry - y)**2 )
            cirlencom(1) = p * (cirx - x)/radius
            cirlencom(2) = p * (ciry - y)/radius
        return 
    end

    !删除矩阵对应第i行并拼接
    real function revector(i,vector)
        implicit none
        integer :: j, k 
        integer, intent(in) :: i
        real, intent(in) :: vector(:,:) 
        integer :: row = size(vector,1), col = size(vector,2)
        real, dimension(row - 1,col), intent(out) :: revector
            forall(j = 1 : i - 1 , k = 1 : col)
                revector(j,k) = vector(j,k)
            end forall

            forall(j = 1 : row - i , k = 1 : col)
                revector(i+j-1,k) = vector(i+j,k)
            end forall
        return
    end

    !删除矩阵对应第i行，第j列并拼接矩阵
    real function remat(i,j,mat)
        use revector
        implicit none
        integer, intent(in) :: i, j 
        real, intent(in) :: mat(:,:) 
        integer :: row = size(mat,1), col = size(mat,2)
        real, dimension(row - 1,col - 1), intent(out) :: remat
            remat = revector(i,mat)
            remat = transpose(revector(j,transpose(mat)))
        return
    end

    !利用LU分解对矩阵求逆的方法
    real function inv(A)
        implicit none
        integer :: N = size(A,1)
        real, intent(in) :: A
        real, intent(out) :: inv
        real, dimension(N,N) :: L, U, invL, invU
        integer :: i, j, k

        !LU分解
        do i = 1 , N 
            do j = 1 , N 
                if(i <= j) then 
                    U(i,j) = A(i,j)
                        do k = 1 , i-1
                            U(i,j) = U(i,j) - L(i,k)*U(k,j)
                        end do
                    L(j,i) = 0.0d0
                else
                    L(j,i) = A(j,i)
                        do k = 1 , j-1
                            L(j,i) = L(j,i) - L(j,k)*U(k,i)
                        end do
                    L(j,i) = L(j,i) / U(i,i)
                    U(j,i) = 0.0d0
                end if
            end do
        end do

        !求L,U的逆
        do i = 1 , N 
            invL(i,i) = 1.0d0 / L(i,i)
            do j = 1 , i-1
                invL(i,j) = -L(i,j) / L(i,i)
            end do
            do j = i+1 , N 
                invL(i,j) = 0.0d0
            end do
        end do

        do i = 1 , N 
            invU(i,i) = 1.0d0 / U(i,i)
            do j = 1 , i-1
                invU(i,j) = 0.0d0
            end do
            do j = i+1 , N 
                invU(i,j) = -U(i,j) / U(i,i)
            end do
        end do

        !求A的逆
        inv = matmul(invU, invL)

    end

end module mathcal

module temp2physi_conduction
    use hf_heat_conduction
    use mathcal

    implicit none

    !输入相关格式
     !只需要知道单元编号，对应结点坐标与关联单元，对应温度，结点即可，形式可以转换如下
     !现假设以获得的单元控制体的编号矩阵为ele
     !ele的形式应该为:1 1 2 3 4 300
     !               2 1 2 5 6 301
     !其中第一列表示单元编号，后四位数字为单元对应角结点编号，且依次是：左下，右下，右上，左上,第6列表示单元温度
    integer ele_row = size(ele,1) !获取总单元数量
    real, dimension(ele_row) :: Elmo !定义单元弹性模量，ele_row表示单元总数
    real, dimension(ele_row) :: poi !定义单元泊松比
    real, dimension(2*node_num,2*node_num) :: k_stiffness = 0 !定义整体刚度矩阵，node_num为结点总数
    real, dimension(2*node_num,3) :: force !定义结点载荷列向量，第一列为结点载荷，第二列为结点x编号，第三列为结点y编号
    real, dimension(2*node_num,3) :: displacement
        !定义结点位移列向量,第一列为结点位移，第二列为结点x编号，第三列为结点y编号
        !example: force = [15.5 1 0] ————X1
        !                 [  0  0 1] ————Y1
        !表示结点1在x方向受载荷15.5,在y方向不受载荷,即p1x=15.5,p1y=0

    real, dimension(2,2) :: cir 
        !保存曲率圆心坐标
        !十字螺旋1/4模型中，其实只有3个圆弧侧，2个对称的正曲率圆弧，1个负曲率圆弧，故cir()可以设为2×2矩阵
        !cir = [cirx1,ciry1] ——————负曲率圆弧的圆心坐标
        !      [ mcir , 0  ] ——————圆心在x轴上的正曲率圆弧圆心坐标，圆心在y轴上的为[0, mcir]
    
    !初始化 
        forall(i = 1 : node_num)
            force(2*i-1,1) = 0
            force(2*i-1,2) = i
            force(2*i-1,3) = 0
            force(2*i,1) = 0
            force(2*i,2) = 0
            force(2*i,3) = i
        end forall
        displacement = force
    !求解
    call solutionD() !得到结点位移列向量

    !Zr与UO2弹性模量关于温度的函数
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

    !Zr与UO2泊松比关于温度的函数
    real function poission(T,type)  !T单位K , type表示材料类型
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
    
    !对单元力学量进行赋值
    subroutine material(ele_uo2,ele_row,Elmo,poi)
        implicit none
        integer :: i
        real, intent(in,out) :: Elmo(:), poi(:)
        integer, intent(in) :: ele_uo2, ele_row 
        forall(i = 1 : ele_uo2)
        !编号1~ele_uo2表示材料为UO2的单元编号，ele_uo2~ele_row表示材料为zr的单元编号,分别赋值
            Elmo(i) = Elasticmo(ele_t(i,2),2)
            poi(i) = poission(ele_t(i,2),2) 
        end forall

        forall(i = ele_uo2 + 1 : ele_row)
            Elmo(i) = Elasticmo(ele_t(i,1),1)
            poi(i) = poission(ele_t(i,1),1)
        end forall
    end subroutine material

    !进行刚度矩阵计算，采取四结点四边形有限元计算
    subroutine stiffness(E,nu,ele,node,type)
        !E为弹性模量，nu为泊松比.node为结点编号、坐标以及边界限制，第一列为结点编号，第二列为x坐标，第三列为y坐标
        !第四列为x方向限制，第五列为y方向限制，取值为0表示固定，1表示自由, 不为0，1的自然数表示对应区域的受外载荷边界
        !node的形式应该为:1 0.01 0.01 0 0 
        !                2 0.03 0.03 1 1
        !type为选择平面应力或平面应变的类型

        call material(ele_uo2,ele_row,Elmo,poi)
        implicit none
        integer, intent(in) :: type
        real, intent(in) :: E(:),nu(:),ele(:,:),node(:,:)
        real, dimension(8,8) :: ele_k_stiffness !单元刚度矩阵
        real, dimension(4,1) :: xreal, yreal !储存单元的四个结点坐标
        real, dimension(4,2) :: sample!样本矩阵
        real, dimension(2,4) :: temp  !过渡用矩阵
        real, dimension(3,8) :: strainmatrix !应变矩阵
        real, dimension(3,3) :: stressmatrix !应力矩阵
        real, dimension(8) :: formcoff !形函数系数矩阵
        real, dimension(8,1) :: thermal_stress !单元热应力向量
        real, parameter :: con = 1/sqrt(3) !定义样本需要的常数
        integer :: i, j, w, k
        real :: a, b, c, d, s, t, coffA, J !过渡系数与雅各比行列式
        sample = reshape([con,con,-con,-con,con,-con,-con,con],[4,2])

        do i = 1 , ele_row
            !定义结点真实坐标，设单元结点编号左下为1，右下为2，右上为3，左上为4
            !同时刚度矩阵积分采取高斯积分近似，选取正负根号3分之1的4个样本点，权重为1
            forall(j = 1 : 4)
            xreal(j,1) = node(ele(i,j+1),2)
            yreal(j,1) = node(ele(i,j+1),3)
            end forall

            ele_k_stiffness = 0.0 !初始化单元刚度矩阵
            thermal_stress = 0.0 !初始化单元热应力矩阵

            !计算第i个单元的刚度矩阵
            do j = 1 , 4
                s = sample(j,1) !基坐标下的样本点横坐标
                t = sample(j,2) !基坐标下的样本点纵坐标
                temp = reshape([ -1+s , -1+t , -1-s ,  1-t ,&
                              &   1+s ,  1+t ,  1-s , -1-t],[2,4])
                temp = multnM(0.25,temp)
                a = matmul(temp,yreal)(1,1)
                b = matmul(temp,yreal)(2,1)
                c = matmul(temp,xreal)(2,1)
                d = matmul(temp,xreal)(1,1)

                !计算雅各布矩阵行列式
                J = a*c - b*d

                !计算形函数系数矩阵
                formcoff = [-1+t , -1+s ,  1-t , -1-s ,&
                           & 1+t ,  1+s , -1-t ,  1-s]
                formcoff = multnV(0.25,formcoff)

                    !计算应变矩阵strainmatrix
                    do w = 1 , 4
                        strainmatrix(1,2*w-1) = 0.25*(a*formcoff(2*w-1) - b*formcoff(2*w))/J
                        strainmatrix(2,2*w)   = 0.25*(c*formcoff(2*w) - d*formcoff(2*w-1))/J
                        strainmatrix(3,2*w-1) = 0.25*(c*formcoff(2*w) - d*formcoff(2*w-1))/J
                        strainmatrix(3,2*w)   = 0.25*(a*formcoff(2*w-1) - b*formcoff(2*w))/J
                    end do

                !计算应力矩阵stressmatrix
                select case(type)
                !type = 1:平面应力
                !type = 2:平面应变
                    case(1)
                        coffA = E(i)/(1-nu(i)**2)
                        stressmatrix = reshape([1,nu(i),0,nu(i),1,0,0,0,(1-nu(i))/2],[3,3])
                        stressmatrix = multnM(coffA,stressmatrix)
                    case(2)
                        coffA = E(i)/(1+nu(i))/(1-2*nu(i))
                        stressmatrix = reshape([1-nu(i),nu(i),0,nu(i),1-nu(i),0,0,0,(1-2*nu(i))/2],[3,3])
                         stressmatrix = multnM(coffA,stressmatrix)
                    case default
                        write(*,*) "unknown input"
                end select

                !利用高斯积分近似，采取四个样本点计算积分
                ele_k_stiffness = ele_k_stiffness + J* matmul(matmul(transpose(strainmatrix),stressmatrix),strainmatrix) 
                call thermal_load(strainmatrix,stressmatrix,thermal_stress,J,i,type) 

            end do 
        
            !将该单元刚度矩阵整合到整体刚度矩阵中
            call assemblestiffness(k_stiffness,ele_k_stiffness,ele(i,:))

            !将热应力整合到载荷向量中
            do k = 1 , 4 
                force( 2*ele(ele_num,k+1) - 1 , 1 ) = force( 2*ele(ele_num,k+1) - 1 , 1 ) + thermal_stress( 2*k - 1 , 1 )
                force( 2*ele(ele_num,k+1) , 1 ) = force( 2*ele(ele_num,k+1) , 1 ) + thermal_stress( 2*k , 1 )
            end do

        end do
    end subroutine stiffness

    !组装刚度矩阵的子程序
    subroutine assemblestiffness(k_stiffness,ele_k_stiffness,ele_num)
     !k_stiffness为整体刚度矩阵，ele_k_stiffness为单元刚度矩阵, ele_num为对应单元编号
        implicit none
        integer :: i, j
        integer, intent(in) :: ele_num(:)
        real, intent(in) :: ele_k_stiffnesss(:,:)
        real, intent(in,out) :: k_stiffnesss(:,:)
        integer :: disx1, disx2, disy1, disy2 !代替坐标,用于确定单元刚度矩阵元素在整体刚度矩阵中的位置

            do  i = 1 , 4
                disx1 = 2*ele_num(i+1) - 1
                disy1 = 2*ele_num(i+1)
                do j = 1 , 4
                    disx2 = 2*ele_num(j+1) - 1
                    disy2 = 2*ele_num(j+1)

                    !将单元刚度矩阵对应结点对应到整体刚度矩阵的对应位置
                    !example:      u4   v4  u7  v7                        7  8     13  14
                    !         u4  [ 1   2   3   4  ]                 7  [ 1  2  ... 3  4  ]  
                    !         v4  [ 5   6   7   8  ]   则在整体刚度   8  [ 5  6  ... 7  8  ]
                    !         u7  [ 9   10  11  12 ]   矩阵的对应为              ...
                    !         v7  [ 13  14  15  16 ]                 13 [ 9  10 ... 11 12 ]
                    !                                                14 [ 13 14 ... 15 16 ]

                    k_stiffnesss(disx2,disx1) = k_stiffnesss(disx2,disx1) + ele_k_stiffness(2*j-1 , 2*i-1)
                    k_stiffnesss(disy2,disx1) = k_stiffnesss(disy2,disx1) + ele_k_stiffness(2*j , 2*i-1)
                    k_stiffnesss(disx2,disy1) = k_stiffnesss(disx2,disy1) + ele_k_stiffness(2*j-1 , 2*i)
                    k_stiffnesss(disy2,disy1) = k_stiffnesss(disy2,disy1) + ele_k_stiffness(2*j , 2*i)
                
                end do
            end do

    end subroutine assemblestiffness

    !创建结点载荷,位移列向量与确定边界条件的子程序
    subroutine load(p,node,force,displacement,cir)
        implicit none
        integer :: i, nodextype, nodeytype !定义结点x与y方向的边界类型
        real, intent(in) :: node(:,:), cir(:,:), p !cir为曲率圆心坐标矩阵，p为外界恒定压力
        real, intent(in,out) :: force(:,:), displacement(:,:)
        real, dimension(2*node_num,3), intent(out) :: tempdisplacement = displacement !用于计算非固定边界的过渡位移向量
        call stiffness(Elmo,poi,ele,node,1)

        !对结点载荷进行赋值
        do i = 1 , node_num
            nodextype = node(i,4)
            nodeytype = node(i,5)
            !nodetype = 0 : 方向固定，同时删去刚度矩阵对应第行与列、过渡载荷和位移列向量的对应行
            !nodetype = 1 : 方向自由，载荷列向量对应行设为0
            !nodetype = 2 : 受x方向载荷的结点类型
            !nodetype = 3 : 受y方向载荷的结点类型
            !nodetype = 4 : 圆心在x轴上的正曲率圆周侧载荷的结点类型
            !nodetype = 5 : 圆心在y轴上的正曲率圆周侧载荷的结点类型
            !nodetype = 6 : 负曲率圆周侧载荷的结点类型

            !判断x方向载荷类型
            select case(nodextype)
                case(0)
                    k_stiffness = remat(2*i-1, 2*i-1, k_stiffness)
                    tempdisplacement = revector(2*i-1, tempdisplacement)
                    force = revector(2*i-1, force)
                case(1,3)
                    force(2*i-1, 1) = 0
                case(2)
                    force(2*i-1, 1) = force(2*i-1, 1) - p
                case(4)
                    force(2*i-1, 1) = force(2*i-1, 1) + cirlencom(p,node(i,2),node(i,3),cir(2,1),cir(2,2))(1)
                case(5)
                    force(2*i-1, 1) = force(2*i-1, 1) + cirlencom(p,node(i,2),node(i,3),cir(2,2),cir(2,1))(1)
                case(6)
                    force(2*i-1, 1) = force(2*i-1, 1) + cirlencom(-p,node(i,2),node(i,3),cir(1,1),cir(1,2))(1)
                case default
                    write(*,*) "未定义载荷类型"
            end select 

            !判断y方向载荷类型
            select case(nodeytype)
                case(0)
                    k_stiffness = remat(2*i, 2*i, k_stiffness)
                    tempdisplacement = revector(2*i, tempdisplacement)
                    force = revector(2*i, force)
                case(1,2)
                    force(2*i , 1) = 0
                case(3)
                    force(2*i , 1) = force(2*i , 1) - p
                case(4)
                    force(2*i, 1) = force(2*i , 1) + cirlencom(p,node(i,2),node(i,3),cir(2,1),cir(2,2))(2)
                case(5)
                    force(2*i , 1) = force(2*i , 1) + cirlencom(p,node(i,2),node(i,3),cir(2,2),cir(2,1))(2)
                case(6)
                    force(2*i , 1) = force(2*i , 1) + cirlencom(-p,node(i,2),node(i,3),cir(1,1),cir(1,2))(2)
                case default
                    write(*,*) "未定义载荷类型"
            end select      

        end do

    end subroutine load

    !热应变可以视为在每个结点上施加了热应力载荷，求热应力载荷的子程序
    subroutine thermal_load(B,D,J,thermal_stress,ele_num,type)
        implicit none
        real, intent(in) :: B, D, J
        integer, intent(in) :: ele_num, type
        real, intent(in,out) :: thermal_stress(:,:)
        real, dimension(3,1) :: tempvt = reshape([ele(ele_num,6) - 300 , ele(ele_num,6) - 300 , 0],[3,1]) !以300K为初始温度
        real :: alpha !热膨胀系数
            !uo2:1e-5
            !zr:6e-6

            if (ele_num <= ele_uo2) then
                alpha = 1e-5
            else
                alpha = 6e-6
            end if

            select case(type)
            !type = 1:平面应力
            !type = 2:平面应变
                case(1)
                    tempvt = multnM(alpha,tempvt)
                case(2)
                    tempvt = multnM(alpha*(1+poi(ele_num)),tempvt)
                case default
                    write(*,*) "unknown input"
            end select

            thermal_stress = thermal_stress + J* matmul(matmul(transpose(B),D),tempvt)

    end subroutine thermal_load

    !求解子程序
    subroutine solutionD()
        !先解出结点非固定坐标的位移向量，最后再把固定位移的结点坐标整合在一起
        call load(15.5,node,force,displacement,cir)
        implicit none
        integer :: i, tempsize
        real, intent(in,out) :: force(:,:), tempdisplacement(:,:), displacement(:,:)
        real, intent(in) :: k_stiffness(:,:)
        tempsize = size(tempdisplacement,1)
        tempdisplacement(:,1) = matmul(inv(k_stiffness),force(:,1)) 
        do i = 1 , tempsize

            if ( tempdisplacement(i,2) /= 0 ) then
                displacement(2*tempdisplacement(i,2) - 1,1) = tempdisplacement(i,1) !将过渡位移对应编号的x方向的位移赋值给位移向量
            end if

            if ( tempdisplacement(i,3) /= 0 ) then
                displacement(2*tempdisplacement(i,2),1) = tempdisplacement(i,1) !将过渡位移对应编号的y方向的位移赋值给位移向量
            end if

        end do
    
    
    end subroutine solutionD


end module temp2physi_conduction