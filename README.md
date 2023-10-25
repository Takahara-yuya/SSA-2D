# SSA-2D
This is a program to employ SSA on the geophysics data(gravity or magnetic)

Author:
Zuwei Huang
hzw1498218560@tongji.edu.cn
School of Ocean and Earth Science, Tongji University
Integrated Geophysics Group

MSSA.cpp: main program

par.dat: startup file.the format below told you how to use it.
input_filename： data/model_dt90inc.grd -->input file
output_reginal: data/reginal -->output file root(reginal field)
output_local: data/local -->output file root(local field)
output_eigenvalue: data/eigen_value.dat -->output eigenvalue analysis
Lr: 20 -->windows size(less then the 0.5 of the rows of the inputdata)
Lc: 20 -->windows size(less then the 0.5 of the cols of the inputdata)
number_of_k: 5 -->how many order of the eigenvalue you want to output
order: 5 7 9 15 20 -->order
reconstruct: 1 -->do you want to reconstruct the signal?(always 1)

reference:
[1] Golyandina N E , Usevich K D , Florinsky I V .Filtering of digital terrain models by 2D Singular Spectrum Analysis[J].International Journal of Ecology and Development, 2007, 8(F07):81-94.
[2] Zhu D, Liu T Y, Li H W. 2018. Separation of potential field based on singular spectrum analysis. Chinese J. Geophys. (in Chinese), 61(9): 3800-3811, doi: 10.6038/cjg2018L0240.
