# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 21:44:42 2025
设置速度阈值并剔除异常值
@author: 杨艺楠
"""
import netCDF4 as nc
import numpy as np

def load_wave_ray_data(path):
    """
    从指定路径的 NetCDF 文件中读取波射线变量：
    rlon、rlat、rzwn、rmwn，并进行 squeeze 处理。

    参数：
        path (str): NetCDF 文件路径

    返回：
        rlon, rlat, rzwn, rmwn: 全部为 ndarray
    """
    with nc.Dataset(path) as ds:
        rlon = ds.variables['rlon'][:].squeeze()
        rlat = ds.variables['rlat'][:].squeeze()
        rzwn = ds.variables['rzwn'][:].squeeze()
        rmwn = ds.variables['rmwn'][:].squeeze()
    
    return rlon, rlat, rzwn, rmwn
#==============================================================================
def threshold(rlon0, rlat0, wn_min, wn_max, conver, time_step, velo_threshold, rzwn=None, rmwn=None, wn_threshold=None, check_wn=False):
    
    # 计算后差 
    rlon_diff = rlon0[1:] - rlon0[:-1]
    rlat_diff = rlat0[1:] - rlat0[:-1]
    
    # 计算每对的中间纬度用于计算每一点是否超速时转换
    lat_center = (rlat0[1:] + rlat0[:-1])/2
    # 先转为弧度，才能计算余弦值
    cos_lat = np.cos(np.deg2rad(lat_center))
    
    
    # 计算差分的平方和
    diff_square = np.sqrt((rlon_diff*conver/time_step)**2 + (rlat_diff*cos_lat*conver/time_step)**2) # 考虑纬度变化

    # 找到满足条件的索引
    exceed_threshold = diff_square > velo_threshold  # 设置速度阈值
    # 遍历每个组合，找到首次超过阈值的位置
    for j in range(3): # rlon80.shape[1] 3个解
        for loc in range(rlon0.shape[2]): # 波源
            for wn in range(wn_min-1, wn_max): # 波数 wave number
                # 找到第一次超过阈值的位置
                exceed_indices = np.where(exceed_threshold[:, j, loc, wn])[0]
                if len(exceed_indices) > 0:
                    first_exceed_idx = exceed_indices[0]
                    # 将超过阈值后面的时间步全部设为 nan
                    rlon0[first_exceed_idx:, j, loc, wn] = np.nan
                    rlat0[first_exceed_idx:, j, loc, wn] = np.nan
                    
                # ✅ 可选：波数异常判断
                if check_wn and rzwn is not None and rmwn is not None and wn_threshold is not None:
                    for var in [rzwn, rmwn]:
                        series = var[:, j, loc, wn]  # shape (T,)
                        if np.isnan(series).all():
                            continue
                        d2 = series[2:] - 2 * series[1:-1] + series[:-2]
                        mask = np.abs(d2) > wn_threshold
                        for t in range(len(mask) - 2):
                            if mask[t] and mask[t + 1] and mask[t + 2]:
                                cut_idx = t + 2
                                rlon0[cut_idx:, j, loc, wn] = np.nan
                                rlat0[cut_idx:, j, loc, wn] = np.nan
                                break
    
    return(rlon0, rlat0)