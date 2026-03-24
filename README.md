# 20260324
运行 *benchmark* 文件夹中的  `EffMainPlot_15/90/165.m`  可以得到文章中的Fig.1。 15/90/165的意思是极向角。其余各种等离子体参数都在程序里进行了注释。这个程序是主要是实现计算：对每一个给定的$n_{\parallel}$， 遍历所有给定的 $y=\frac{\ell \omega_{c}}{\omega}$ 值，在每一个数据点上结合其他等离子体参数计算ECCD和SCD的效率。


*Support* 文件夹中是各种所需函数的显示表达式。文章中的各种表达式为了方便调用被写为函数。另外，为了提高计算效率，我把它们用Mathematica写成了方便数值计算的形式。

该代码仅作学术用途，不可用于商用等牟利，使用此代码时应该引用文章：
```
@article{Chen_2026,
doi = {10.1088/1741-4326/ae1fa0},
url = {https://doi.org/10.1088/1741-4326/ae1fa0},
year = {2025},
month = {nov},
publisher = {IOP Publishing},
volume = {66},
number = {1},
pages = {016035},
author = {Chen, Shao-Yong and Chen, Jia-Ning and Zhang, Ye-Min and Mou, Mao-Lin and Tang, Chang-Jian},
title = {A linear model of synergetic current drive with lower-hybrid wave and electron cyclotron wave},
journal = {Nuclear Fusion},
abstract = {A linear model of synergetic current drive (SCD) with lower-hybrid wave (LHW) and electron cyclotron wave (ECW) is developed to quickly calculate the quantitative SCD efficiency and reveal the conditions for the occurrence of the synergy effect. In this model, the response function dominated by collisions in the presence of LHW is derived from the adjoint equations by using perturbation and Green-function techniques, where the relativistic and trapping effects are considered. The SCD results compared to LinLiu’s ECW current drive (ECCD) efficiency show two features of the synergy effect, one is that it is inclined to occur at smaller  with the fixed ECW parallel refractive index , and the other is that the threshold values of y, at which the synergy effect becomes sufficiently significant, shifts towards higher values with a decreasing ; the SCD results without the trapping effect and the response function without relativistic effect have been provided. The quasilinear simulation on ECCD and SCD efficiency with a two-dimensional Fokker–Planck code is consistent with the results of the linear model in trends. Given considerable parameters improvement nowadays, the SCD efficiencies with parameters of current mainstream reactors are calculated in this work. We provide the current profile of ECCD and SCD by combining with GENRAY, showing the potential application of integrating SCD into ray-tracing code. Criteria for the occurrence and the sufficient significance of the synergy effect are suggested, which indicate that the synergy effect is dependent on the power factor that quantifies the degree of the overlap of the two waves’ quasilinear domains, the LHW power, and synergy electrons.}
}
```
