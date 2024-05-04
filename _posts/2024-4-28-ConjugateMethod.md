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

# 重要性质
## Conjugate Method至多$n$步到达最优解，这里$n$为优化变量的特征维度
$Definition1:$ $d_i,d_j, \forall i \ne j$ are A-conjugate $\Leftrightarrow d_i^TAd_j =0$.


$Definition2:$ A set ${d_0,d_1,...,d_{n-1}}$ is A-conjugate $\Leftrightarrow d_i^TAd_j =0,\forall (i\ne j) \cup(i,j=0,1,...,n-1)$.


$Property1:$ If a set ${d_0,d_1,...,d_{n-1}}$ is A-conjugate, then all components i.e.$d_0,d_1,...,d_{n-1}$ are linearly independent，i.e. these vectors can span whole $\mathbb{R}^n$ space($n$ is feature-dim).
$\underline{proof:}$ 
>反证法: 假设满足上述条件，set ${d_0,d_1,...,d_{n-1}}$是线性相关的，那么必然存在某个向量能够用其他向量线性表示。
也就是说集合$\Lambda=\{d_0,d_1,...,d_{n-1}\},\Omega=\{0,1,2,...,n-1\},\exists k \in \Omega,s.t.d_k\in\Lambda$ $and$ $ d_k=\sum_{j\in\Omega/k}\alpha_j d_j$。
两边同时乘$d_k^TA$，有，
$$
d_k^TAd_k= \sum_{j\in\Omega/k}\alpha_j d_k^TAd_j=0
$$
由于$d_k\ne \bm{0}$，因此和假设矛盾，得证。

