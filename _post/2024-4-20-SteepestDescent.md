---
layout:     post   				    # 使用的布局（不需要改）
title:      Steepest Descent 				# 标题 
subtitle:     #副标题
date:       2024-04-20 				# 时间
author:     BlueArt 						# 作者
header-img: img/post-bg-2015.jpg 	#这篇文章标题背景图片
catalog: true 						# 是否归档
tags:								#标签
    - Convex Optimization
    - Steepest Descent
---


# 问题引入
本质是求二次规划问题(Quadratic Form)的最小值,形式如下
$$x^* \in \min_{x} \, f(x) = \min_{x}\,\frac{1}{2}x^T A x -b^T x+c$$
如果$A$是对称正定，那么该问题其实就相当于求解方程$Ax^* = b$，因为根据求极值公式，
$$f^{'}(\bar{x})=0$$
即可求出$\bar{x}$对应的最小值，又因为$f^{'}(\bar{x})=\frac{1}{2}(A+A^T)\bar{x}-b$，所以当$A=A^T$，$f^{'}(\bar{x})=0$等价为$A\bar{x}-b=0$，所以说求解Quadratic Form形式的最优解就相当于求解一个线性系统(linear system)。***那么$A\in \mathbb{R}_{++}^n$又有什么用呢？*** 也就是说他能给解决问题带来什么好处呢？\
All in all, $A$是正定矩阵，严格保证了满足$Ax=b$的就是问题的全局最优解(global optimal)。假定在最优点$x^*$附近的任意扰动$e$，也就是说$x^*$满足$Ax^8=b$，那么根据$f(x)$表达式，可以知道，
$$
\begin{equation}
\begin{aligned}
 f(x^*+e)&=\frac{1}{2}(x^*+e)^T A(x^*+e)-b^T(x^*+e)+c \\
         &=\frac{1}{2}(x^{*})^TAx^* - b^Tx^* +c + \frac{1}{2}e^TAe-b^Te+e^TAx^* \\
         &= f(x^*)+e^Tb+\frac{1}{2}e^TAe-b^Te \quad \bm{(使用等式Ax^*=b证明该式)}\\
         &= f(x^*) + \frac{1}{2}e^TAe

\end{aligned}
\end{equation}$$
由于$A$是正定矩阵，所以$e^TAe >0$，也就是说任何异于$x^*$点都不可能取到最优值，也就是说最优值有且只有一个$x^*$。从几何直观上理解，$A$为正定矩阵，有且只有一个最优解，$A$为半正定矩阵，可能有解，但是一旦有解就可能会有多个对应相同最优值的最优解。***注意$f(x)$中的系数各有作用，$A$控制着形状，特征值$\lambda(A)$控制着轴的长度，特征向量规范化所决定的正交矩阵决定了旋转的方向；而$A$和$b$共同控制着最优解的取值，而他们一起控制最优值。***

# Steepest Descent
下降法都有一个一般的更新公式,
$$x(k+1)=x(k)+\alpha(k)d(k)$$
其中$\alpha(k)$为第k代的步长(step size)，$d(k)$为第k代的搜索方向(search direction)。一般这种算法都会规定一个搜索方向，然后根据精确的(exact)或者不精确的(inexact)的搜索来决定步长，这两个确定步长的差异在于，exact方法确定的步长是精确的，并且可以较少代数收敛，但是每次计算exact步长的代价较大，inexact方法的步长是不精确的，需要较多代数收敛，但是每次计算代价小，可以得到工程应用的近似最优解。  
## Search Direction 

同样的，令$d(k)$为待优化问题的负梯度方向便可以得到steepest descent的搜索方向，即$d(k)=-\nabla f(x(k))=b-Ax(k)$，$b-Ax(k)$的物理意义是残差(residual)，于是定义$r(k)=b-Ax(k)=-f^{'}(x(k))$，误差$e(k)=x(k)-x^*$为偏离最优解的程度。那么$r(k)$可以表示为$r(k)=Ax^* -Ax(k) = A(x^*-x(k))=-Ae(k)$，于是有以下关系，
$$
\begin{equation}
\begin{aligned}
 d(k)&=b-Ax(k)=-f^{'}(x(k))=-Ae(k) \\
 r(k)&=b-Ax(k)=-f^{'}(x(k))=-Ae(k)\\
 e(k)&=x(k)-x^*=-A^T f^{'}(x(k))= -A^T r(k)\\
 r(k)&=d(k)
\end{aligned}
\end{equation}$$
搜索方向即为残差方向。
## Step Size 
### Exact Linesearch
如果当前在第$0$代，那么这个时候要选择一个第$0$代的步长使得前进到$x(1)$时对应的$f(x(1))$最小，同样其他代数也可以类似的描述，如果当前在第$k$代，那么这个时候要选择一个第$k$代的步长使得前进到$x(k+1)$时对应的$f(x(k+1))$最小，数学描述如下，
$$
\begin{equation}
\alpha(k) = \mathop{\arg\min}_{\alpha}f(x(k+1))
\end{equation}
$$
由于$f(x(k+1))$是可微的，所以对应的$\alpha(k)$可以用求极值条件求得，
$$
\begin{equation}
\begin{aligned}
\frac{d}{d \alpha} f(x(k+1)) &= \frac{d}{d \alpha} f(x(k)+\alpha d(k)) \\
          &= [f^{'}(x(k)+\alpha d(k))]^T \, d(k) \\
          &= [f^{'}(x(k)+\alpha d(k))]^T \,r(k) \\
          &= 0
\end{aligned}
\end{equation}
$$
因此有一个重要条件$\bm{\langle f^{'}(x(k+1)), d(k)\rangle=0}$，等价于$\bm{\langle d(k+1), d(k)\rangle=\langle r(k+1), r(k)\rangle=0}$，这个条件只是精确的线搜索的一个by-product，下面确定$\alpha$取值，
$$\begin{equation}
\begin{aligned}
    &r^T(k+1)r(k)=0 \\
    &\Leftrightarrow -f^{'}(x(k+1))r(k)=0 \\
    &\Leftrightarrow  (b-Ax(k+1))^T r(k) = 0\\
    &\Leftrightarrow  \{b-A[x(k)+\alpha(k)r(k)] \}^T r(k) = 0\\
    &\Leftrightarrow  \alpha(k)r^T(k)A^T r(k) = (b-Ax(k))^Tr(k) \\
    &\Leftrightarrow  \alpha(k) = \frac{r(k)^T r(k)}{r^T(k)A^T r(k)}=\frac{d(k)^T d(k)}{d^T(k)A^T d(k)} \quad \bm{(A=A^T,所以这里A和A^T可以互相转化)}
\end{aligned}
\end{equation}$$
于是可以得到下面两套算法流程，\

> $\underline{Algorithm 1: Steepest Descent with Exact Line Search}$\
  $Initializing:$ $arbitrarily$ $select$ $x(0)$ \
  $for$ $k = 0,1,2,3...$ \
    $\quad \quad d(k) = b - Ax(k)$\
    $\quad \quad \alpha(k)=\frac{d(k)^T d(k)}{d^T(k)A^T d(k)}$\
    $\quad \quad x(k+1)=x(k)+\alpha(k)d(k)$ \

算法二则是将算法一的$x$的更新表达式两边同时乘矩阵$A$,然后同时加上$b$，可以得到新的$d(k+1)$更新表达式为$d(k+1) = d(k)-A\alpha(k)d(k)$，于是便有算法二\
> $\underline{Algorithm 2: Another Steepest Descent with Exact Line Search}$\
  $Initializing:$ $arbitrarily$ $select$ $x(0),d(0)=b-Ax(0)$ \
  $for$ $k = 0,1,2,3...$ \
    $\quad \quad \alpha(k)=\frac{d(k)^T d(k)}{d^T(k)A^T d(k)}$\
    $\quad \quad x(k+1)=x(k)+\alpha(k)d(k)$\
     $\quad \quad d(k+1)=d(k)-A\alpha(k)d(k)$

算法二方向(也可以说残差)更新，没有用到$x(k)$，因此会造成没有反馈(feedback)，因此积累的浮点数误差会使得算法二中$x$只能收敛到$x^*$的邻域内。

***一些重要发现***
1. 商业求解器SCS求解Quadratic form形式问题不能到达最优解，仍然有一定误差
2. steepest descent exact line search形式在第n代（n为feature dim，也就是待优化变量的维数）就可以收敛，如果算法继续运行，特别是当前解接近最优解时算法速度非常慢。

### Inexact Linesearch
