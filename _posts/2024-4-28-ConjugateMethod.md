---
layout:     post   				    # 使用的布局（不需要改）
title:      Conjugate Method 				# 标题 
subtitle:     #副标题
date:       2024-04-28 				# 时间
author:     BlueArt 						# 作者
header-img: img/post-bg-2015.jpg 	#这篇文章标题背景图片
catalog: true 						# 是否归档
tags:								#标签
    - Convex Optimization
    - Conjugate Method
---

# 共轭方法由来

最速下降法只选取正交的方向，难免出现后面迭代的方向与之前迭代的方向是相同的，这多个相同的方向的迭代可以通过选取恰当的步长一次迭代替代即可，可以想象二维平面有两点$A,B$，从$A$到$B$是走阶梯状的曲线，你可以上下上下上下多次重复，也可以上三步下三步走到B点。
自然想到令当前代$k$的方向和$k+1$代的残差正交，可以表示为$\langle d(k),e(k+1) \rangle=0$，可以继续转化$d(k)^T(e(k)+\alpha(k)d(k))=0$，于是
$$\alpha(k)=-\frac{d(k)^Te(k)}{d(k)^Td(k)}$$
这里$e(k)$是未知的，因为$e(k)=x(k)-x^*$，最优解此时还未求出。之前利用已知搜索方向，求搜索步长使得目标函数最小，然后中间得到了一个副产品(by-product)，也就是搜索方向和残差的关系，
$$
\begin{split} 
\frac{d}{d\alpha}f(x(k+1)) &=\frac{d}{d\alpha}f(x(k)+\alpha d(k)) \\
&=f^{'}(x(k+1))^T d(k)=0
\end{split}
$$
残差的定义还是和最速下降法一样，$r(k)=b-Ax(k)=-f^{'}(x(k))=-Ae(k)$，于是上式可以写成$r(k+1)^Td(k)=0$，最速下降法直接令$d(k)=r(k)$，然而共轭方法则是做了一步代换$r(k+1)=-Ae(k+1)$，之后便得到$\bm{d(k)^TAe(k+1)=0}$，这么做的理由是：***为了使得算法能够避免掉不必要的搜索，之前的想法是找满足$d(k)^Te(k+1)=0$的搜索方向，但是这种方法无法求解，于是延续这个思想，找$d(k)$和$e(k+1)$的关系式就可以不必要垂直，这放松了条件使得$d(k)$可能通过$e(k+1)$表示。***
于是共轭方法各个变量有如下关系，
$$
\begin{equation}
\begin{split}
    &r(k)=b-Ax(k)=-f^{'}(x(k))=-Ae(k) \quad\bm{(这个定义和最速下降法一致)} \\
    &e(k)=x(k)-x^*=-A^T f^{'}(x(k))= -A^T r(k)\\
    & d(k)^TAe(k+1)=0
\end{split}
\end{equation}
$$
这里只知道$d(k)$的隐式表达，而不知道显示表达，其实这个隐式表达就已经决定了$d(k)$的方向，而不必具体知道$d(k)$的数值，因为具体的数值可以使用$\alpha(k)$修正。

# Step Size
假设$d(k)$已知，
$$
\begin{equation}
\begin{split}
    &d(k)^TAe(k+1)=0 \\
    \Leftrightarrow \quad &d(k)^TA(e(k)+\alpha(k)d(k))=0 \\
    \Leftrightarrow \quad &\alpha(k)=-\frac{d(k)^TAe(k)}{d(k)^TAd(k)} \\
    \Leftrightarrow \quad &\alpha(k)=\frac{d(k)^Tr(k)}{d(k)^TAd(k)}=\frac{d(k)^T(b-Ax(k)))}{d(k)^TAd(k)}
\end{split}
\end{equation}
$$
如果搜索方向就是残差($d(k)=r(k)$)，那么共轭方法就会蜕化成最速下降法。

# The Selection of Conjugate Direction

使用$Gram-Schmidt Procedure$将一组线性独立的向量转化成一组$orthogonal$或者$conjugate$的向量，方法描述如下，
1. 给定一组$linearly indepentdent$向量$u_0,u_1,...,u_n$，如果这些向量是$A-conjugate$，那么他们是必然线性独立的。
2. 令$d_0=u_0$,其余的共轭方向可以选择为$d_k=u_k+\sum_{i=0}^{k-1}\beta_{k,i}u_i$
3. $\forall j \in {0,1,2,...,k-1},d_k^TAd_j=u_k^TAd_j+\sum_{i=0}^{k-1}\beta_{k,i}u_i^TAd_j=u_k^TAd_j+\beta_{k,j}u_j^TAd_j=0$\
$\Rightarrow \beta_{k,j}=-\frac{u_k^TAd_j}{u_j^TAd_j}$
因此,
$$
d_k= \begin{cases} u_0  \qquad\qquad \qquad\qquad,k=0 \\ u_k - \sum_{i=0}^{k-1}\frac{u_k^TAd_i}{d_i^TAd_i} d_i\quad ,else\end{cases}
$$

# 一群正交和A-conjugate性质

1. $d_i^TAe_k=0,\forall i \in \{0,1,2,...,k-1\}$
2. $d_i^Tr_k=0,\forall i \in \{0,1,2,...,k-1\}$
3. $u_i^T r_k =0,\forall i \in \{0,1,2,...,k-1\}$

*证明见附录。*
根据上一节最速下降法已经证明，优化$f(x_k)$等价于优化$||e_k||_A$，类似的共轭方法也可以从两个角度分析，
- 从$||e_k||_A$角度分析，由于$e_k=e_0+\sum_{j=0}^{k-1}\alpha_jd_j$，每次迭代就是选择一组$d_0,d_1,...,d_{k-1}$的线性组合，如果记$\mathcal{D}_i=span\{d_0,d_1,...,d_{i-1}\}$，那么每次迭代就是选择$e_k\in e_0+\mathcal{D}_{k}$，使得$||e_k||_A$最小。
- 从$f(x_k)$角度分析，由于$x_k=x_0+\sum_{j=0}^{k-1}\alpha_jd_j$，和上面一样每次迭代就是选择一组$d_0,d_1,...,d_{k-1}$的线性组合，如果记$\mathcal{D}_i=span\{d_0,d_1,...,d_{i-1}\}$，那么每次迭代就是选择$x_k\in x_0+\mathcal{D}_{k}$，使得$f(x_k)$最小。

直观上看，如果把$f(x)$画成空间中一簇簇等高线，第一次迭代就是找一条从$x_0$出发的线，使得沿着这条线$f(x_k)$最小；第二次迭代就是找一个面，从$x_0$出发，使得$f(x_k)$最小，同样也可以等价理解成找一个从$x_1$出发的线，并且这条线和之前的是$A-conjugate$的，使得沿着这条线$f(x_k)$最小。
# 重要性质
## 1. Conjugate Method至多$n$步到达最优解，这里$n$为优化变量的特征维度
$Definition1:$ $d_i,d_j,$  are $A-conjugate$ $\Leftrightarrow d_i^TAd_j =0，\forall i \ne j$$.
通过这个定义可以推出两个$A-conjugate$的向量不能为$0$，如果两个向量至少其中一个为0时，那么$i=j$时，$d_i^TAd_i=0$依然成立，和假设矛盾。

$Definition2:$ A set ${d_0,d_1,...,d_{n-1}}$ is A-conjugate $\Leftrightarrow d_i^TAd_j =0,\forall (i\ne j) \cup(i,j=0,1,...,n-1)$.


$Property1:$ If a set ${d_0,d_1,...,d_{n-1}}$ is $A-conjugate$, then all components i.e.$d_0,d_1,...,d_{n-1}$ are linearly independent，i.e. these vectors can span whole $\mathbb{R}^n$ space($n$ is feature-dim). If $d_i \in \mathbb{R}^n$, then there exists at most $n$ components that are $A-conjugate.$
$\underline{proof:}$ 
>反证法: 假设满足上述条件，set ${d_0,d_1,...,d_{n-1}}$是线性相关的，那么必然存在某个向量能够用其他向量线性表示。
也就是说集合$\Lambda=\{d_0,d_1,...,d_{n-1}\},\Omega=\{0,1,2,...,n-1\},\exists k \in \Omega,s.t.d_k\in\Lambda$ $and$ $ d_k=\sum_{j\in\Omega/k}\alpha_j d_j$。
两边同时乘$d_k^TA$，有，
$$
d_k^TAd_k= \sum_{j\in\Omega/k}\alpha_j d_k^TAd_j=0
$$
由于$d_k\ne \bm{0}$，因此和假设矛盾，得证。

由于$d(0),d(1),...,d(n-1)$是线性无关的，同时由于$d_i \in \mathbb{R}^n$，于是$d(0),d(1),...,d(n-1)$可以$span$整个$\mathbb{R}^n$空间，于是任意$\mathbb{R}^n$向量都可以用$d(0),d(1),...,d(n-1)$线性表示。于是，
$$e(0)=\sum_{j=0}^{n-1} \delta(j) d(j)$$
两边同时乘$d(k)^TA$，于是，
$$ 
\begin{split}
d(k)^TAe(0) &= \sum_{j=0}^{n-1} \delta(j) d(k)^TAd(j) \\
&= \delta(k) d(k)^TAd(k)
\end{split}
$$
因此，
$$
\begin{split}
\delta(k) &= \frac{d(k)^TAe(0)}{d(k)^TAd(k)} \\
&= \frac{d(k)^TAe(0)+d(k)^TA(\sum_{j=0}^{k-1} \alpha(j)d(j))}{d(k)^TAd(k)}\\
&= \frac{d(k)^TA[e(0)+\sum_{j=0}^{k-1} \alpha(j)d(j)]}{d(k)^TAd(k)} \\
&= \frac{d(k)^TAe(k)}{d(k)^TAd(k)} 
\end{split}
$$
根据上面求得
$$\alpha(k)=-\frac{d(k)^TAe(k)}{d(k)^TAd(k)}$$
于是$\delta(k)=-\alpha(k)$。于是
$$
\begin{split}
e(k) &= e(0) + \sum_{j=0}^{k-1} \alpha(j)d(j) \\
&= \sum_{j=0}^{n-1} \delta(j) d(j) - \sum_{j=0}^{k-1} \delta(j)d(j) \\
&= - \sum_{j=k}^{n-1} \delta(j)d(j) \\
&= \sum_{j=k}^{n-1} \alpha(j)d(j)
\end{split}
$$
对比两个式子,
$$
\begin{split}
e(k) &= \sum_{j=k}^{n-1} \alpha(j)d(j) \\
x(k) &= x(0) + \sum_{j=0}^{k-1} \alpha(j)d(j)
\end{split}
$$
随着迭代的进行，$k$逐渐增大，误差函数所加的项数逐渐减少，优化变量所加项数逐渐增加，当刚好到第$n$代时误差函数值为0，达到最优解。

## 2. 如果$u_0,u_1,...,u_n$向量中下个向量$u_{n+1}$能够其余向量线性表示，那么由这个向量产生的$d_{n+1}$一定是0
通过这个结论可知，在产生一组A-conjugate的向量时，没有必要指定初始向量都是线性独立的，只需要在计算$d_k$时候舍弃掉为0的那些向量即可，剩余的向量同样会span成和$u_0,u_1,...,u_n$的相同空间。

$Lemma1:$ 如果$d_0,d_1,...,d_n$是$A-conjugate$ $\Rightarrow$ $d_i\ne 0,\forall i \in \{0,1,2,...,n+1\}$

> $\underline{Proof:}$ $d_0,d_1,...,d_n$是$A-conjugate$，那么其中任意一项都不是$0$是显然的，下面证明$d_{n+1}$的情况，
如果$d_{n+1}=0$，那么，
$$d_{n+1}= u_{n+1} - \sum_{i=0}^{n}\frac{u_{n+1}^TAd_i}{d_i^TAd_i} d_i=0$$
可以得到，
$$u_{n+1} = \sum_{i=0}^{n}\frac{u_{n+1}^TAd_i}{d_i^TAd_i} d_i$$
可以看出$u_{n+1}$可以由$d_0,d_1,...,d_n$线性表示，由于$d_0,d_1,...,d_n$和$u_0,u_1,...,u_n$可以span相同的空间，因此说明$u_{n+1}$可以由$u_0,u_1,...,u_n$线性表示，也就是说$u_0,u_1,...,u_n,u_{n+1}$是线性相关的，这与假设$u_0,u_1,...,u_n,u_{n+1}$线性独立相矛盾。$\quad \blacksquare$


$Thm1:$ 如果$u_0,u_1,...,u_n$向量中下个向量$u_{n+1}$能够其余向量线性表示，那么由这个向量产生的$d_{n+1}$一定是0

> $\underline{Proof:}$ 
根据条件，
$$u_{n+1} = \sum_{i=0}^{n} \lambda_i u_i$$
于是，
$$
\begin{split}
d_{n+1} &= u_{n+1} - \sum_{i=0}^{n}\frac{u_{n+1}^TAd_i}{d_i^TAd_i} d_i \\
&=\sum_{i=0}^{n} \lambda_i d_i-\sum_{i=0}^{n}\mu_i d_i \quad \bm{(u_i和d_i可以span出相同空间,\forall i\in \{ 0,1,2,...,n\})} \\
&= \sum_{i=0}^{n} \gamma_i d_i
\end{split}$$
上面等式两边同时乘$d_k^TA,\forall k \in \{ 0,1,2,...,n\}$，
$$
\begin{split} 
&d_k^TAd_{n+1} = \sum_{i=0}^{n} \gamma_i d_k^TAd_i \\
\Leftrightarrow &d_k^TAd_{n+1} = 0=\gamma_k d_k^TAd_k \\
&\bm{由于d_k^TAd_k\ne 0,\forall k\in \{0,1,2,...,n\}}\\
\Leftrightarrow &\gamma_k=0 \\
\Leftrightarrow &d_{n+1}=0 \qquad \blacksquare\\
\end{split}$$ 


##  3.1 $d_i^TAe_k=0,\forall i \in \{0,1,2,...,k-1\}$ 
$e_k=\sum_{j=k}^n \alpha_j d_j$同时左乘$d_i^TA,\forall i\in \{0,1,2,...,k-1\}$
##  3.2 $d_i^Tr_k=0,\forall i \in \{0,1,2,...,k-1\}$ 
$r_k=-Ae_k$，因此该式和上面一个性质是等价的
##  3.3 $u_i^T r_k =0,\forall i \in \{0,1,2,...,k-1\}$
$d_i=u_i+\sum_{j=0}^{i-1}\beta_{i,j}u_j$转置之后右乘$r_k$，根据上面那些性质即得证。
有个副产品，由于
$d_i^Tr_k=u_i^Tr_k+\sum_{j=0}^{i-1}\beta_{i,j}u_j^Tr_k$，当$k=i$，$\sum_{j=0}^{i-1}\beta_{i,j}u_j^Tr_k$仍然为0，因为残差和他之前的所有搜索方向 **(转化前转化后都成立)** 正交，于是有$\bm{d_i^Tr_i=u_i^Tr_i}$.