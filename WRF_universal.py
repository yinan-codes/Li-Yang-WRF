# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 16:51:11 2025
3.26封装出三个函数：截断速度，筛矩形框与计算wave ray flux (WRF)
5.7完善筛矩形框的函数；加入波数限制条件的选择；为避免调整坐标轴的不便，加入两轮地图（0-720°）
5.10加入计数到达目标区的波射线数量
5.29测试完好
6.12在WRF_calu时加入每个网格穿越的波射线数量 ray_count
7.23在region_threshold增加参数判断经过筛选区域是否进行截断
9.9在WRF_calu新增波列传播到对应网格的平均时间wave_propagation_time（单位：day）；推荐使用pcolormesh进行马赛克填色，形式会比较漂亮
   同时新增波列在对应网格的平均传播速度wave_propagation_u, wave_propagation_v
9.11为加强使用者对计算结果形式的理解同时方便将计算结果绘图，在WRF_calu中一并输出对应的lon, lat。感谢王琳洁师姐的建议
9.16诊断过程中发现存在东半球的波向西传到西半球，因此在西侧再加入一轮坐标（-360°-0°）。感谢朱勉师兄的发现与建议
9.18由于后续WRF计算涉及到累加，将0值输出为nan，改为输出0，使用者可以绘图时自行设置
9.23对Fun2_region_threshold的输出改进为如果没有波射线经过目标区，则输出空数组，而不是报错；
    对WRF_calu输出追溯波源的变量source_count，适用于全球波源
10.8对Fun1_threshold将经向波数连续不稳定的时间长度L设置为外部参数。感谢张公俊师兄提供的实例
12.5对Fun1_threshold中截断速度的计算进行调整，感谢刘一凡同学的提示

@author: 杨艺楠1; supervisor: 李建平 教授1,2
1 Frontiers Science Center for Deep Ocean Multi-spheres and Earth System (DOMES)/Key Laboratory of Physical Oceanography/
  Academy of Future Ocean/College of Oceanic and Atmospheric Sciences/Center for Ocean Carbon Neutrality, 
  Ocean University of China, Qingdao 266100, China.
2 Laboratory for Ocean Dynamics and Climate, Qingdao Marine Science and Technology Center, Qingdao 266237, China.
Email: yyn2064@stu.ouc.edu.cn
$References$
Yang, Y. N. and J. P. Li, 2025: Novel monsoon indices based on vector projection and directed angle 
for measuring the East Asian summer monsoon. Clim. Dyn., 63, 210, https://doi.org/10.1007/s00382-025-07696-7
"""

import numpy as np
import sys
sys.path.append('D:/Function/to/your/path')

from Fun1_threshold import load_wave_ray_data
from Fun1_threshold import threshold
from Fun2_region_threshold import region_threshold
from Fun3_WRF_calculate import WRF_calu


# 初始参数======================================================================
a = 6371e3; # 地球半径
space = 4; # WRF计算的空间分辨率（网格大小）space° x space°
time_step = 3600 # 时间步长 单位：s
velo_threshold = 150 # 截断速度阈值（可选）
conver = np.full((1), np.deg2rad(1) * a) # 经向换算因子（单位：m/°） 
wn_min = 3; wn_max = 5 # 所需波数的最小，最大值；默认跑全了1-7波
target_region = (90 + 360, 140 + 360, 15, 40) # 设置目标区域参数 格式为lon_min; lon_max; lat_min; lat_max 注意需要按照-360至720的经度来（三圈Mercator地球坐标）；
lon_min, lon_max, lat_min, lat_max = target_region
# =============================================================================

path1 = 'D:/data/to/your/path' # 波射线计算结果（nc）
rlon, rlat, rmwn = load_wave_ray_data(path1)

# 设置速度阈值并剔除异常值=======================================================
rlon1, rlat1 = threshold(rlon, rlat, wn_min, wn_max, conver, time_step, velo_threshold, 
                          rzwn=None, rmwn=rmwn, check_wn=False, wn_threshold=7, L=5) # check_wn 是否进行经向波数制约


# cut = True表示离开目标区时截断；False保留完整轨迹 
rlon_selected, rlat_selected, entered_count, multi_entry_count = region_threshold(lon_min, lon_max, lat_min, lat_max, rlon1, rlat1, wn_min, wn_max, cut=True);

# 计算wave ray flux (WRF)矢量，波射线数量，波列传播到对应网格的平均时间（单位：day），波列在对应网格的平均传播速度
lon, lat, WRF_u, WRF_v, ray_count, wave_propagation_time, wave_propagation_u, wave_propagation_v, source_count = WRF_calu(rlon_selected, rlat_selected, a, space, conver, time_step)
