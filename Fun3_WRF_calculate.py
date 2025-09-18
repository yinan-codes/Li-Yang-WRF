# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 09:07:53 2025

@author: 杨艺楠
"""
import numpy as np

def WRF_calu(rlon_selected, rlat_selected, a, space, conver, time_step):
    # 经度和纬度格点数
    num_lon = int(1080/space) + 1           # 经度（-360到720度）
    num_lat = int(180/space) + 1            # 纬度（-90到90度）
    
    # 生成经度和纬度格点（中心点）
    lon = np.linspace(-360, 720, num_lon, endpoint=True)   # -360 ~ 720
    lat = np.linspace(-90, 90, num_lat, endpoint=True)
    
    # 设置网格分箱边界（闭区间端点+step）
    lon_bins = np.arange(-360, 721, space)  # -360 ~ 720（包含 720） # 经度每 space° 一格
    lat_bins = np.arange(-90, 91, space)                            # 纬度每 space° 一格
    num_lon_bins = len(lon_bins) 
    num_lat_bins = len(lat_bins) 
    # num_lon = int(720/space)+1
    # num_lat = int(180/space)+1
    # # 生成经度和纬度格点（中心点）
    # lon = np.linspace(0, 720, num_lon, endpoint=True)  # 经度（0到720度）
    # lat = np.linspace(-90, 90, num_lat, endpoint=True)  # 纬度（-90到90度）
    # # 设置网格
    # lon_bins = np.arange(0, 721, space)  # 经度每 space° 一格
    # lat_bins = np.arange(-90, 91, space)  # 纬度每 space° 一格
    # num_lon_bins = len(lon_bins) 
    # num_lat_bins = len(lat_bins) 
    
    # 预分配 NaN 数组，最后一维增加1来记录时间索引
    shape = (num_lon_bins, num_lat_bins, rlon_selected.shape[1])
    first_entry_lon = np.full(shape, np.nan)
    first_entry_lat = np.full(shape, np.nan)
    last_exit_lon = np.full(shape, np.nan)
    last_exit_lat = np.full(shape, np.nan)
    ray_count = np.zeros((num_lon_bins, num_lat_bins))  # 每个网格穿越的波数

    
    # 新增数组记录首次进入和最后离开的时间步索引
    first_entry_time = np.full(shape, np.nan)
    last_exit_time = np.full(shape, np.nan)
    # 遍历波射线
    # 每条波射线（3,n,7）条都有一个切片（lon, lat）记录这条波射线在全球网格的进、出
    for n in range(rlon_selected.shape[1]):  # 遍历波的索引
        # 将位置索引映射到网格
        lon_indices = np.clip(np.digitize(rlon_selected[:, n], lon_bins) - 1, 0, num_lon_bins - 1)
        lat_indices = np.clip(np.digitize(rlat_selected[:, n], lat_bins) - 1, 0, num_lat_bins - 1)
    
        visited_grid = {}  # 记录访问的网格
    
        # 找到非 nan 的最后一个索引
        valid_idx = np.where(~np.isnan(rlon_selected[:, n]))[0]
        if len(valid_idx) > 1:
            last_valid_idx = valid_idx[-1]  # 最后一个非nan索引
    
            for i in range(last_valid_idx - 1):
                grid = (lon_indices[i], lat_indices[i])  # 当前时间步的位置在网格中的索引
                next_grid = (lon_indices[i + 1], lat_indices[i + 1])  # 下一个时间步的位置在网格中的索引
    
                # 🔹 记录首次进入
                if grid not in visited_grid:
                    first_entry_lon[grid[0], grid[1], n] = rlon_selected[i, n]
                    first_entry_lat[grid[0], grid[1], n] = rlat_selected[i, n]
                    # ➡️ 记录时间步
                    first_entry_time[grid[0], grid[1], n] = i
                    visited_grid[grid] = True
                    
                    ray_count[grid[0], grid[1]] += 1
    
                # 🔹 记录真正离开
                if grid != next_grid:
                    last_exit_lon[grid[0], grid[1], n] = rlon_selected[i + 1, n]
                    last_exit_lat[grid[0], grid[1], n] = rlat_selected[i + 1, n]
                    # ➡️ 记录时间步
                    last_exit_time[grid[0], grid[1], n] = i + 1
                    
    delta_lon = last_exit_lon - first_entry_lon #要转换成实际的米！
    delta_lat = last_exit_lat - first_entry_lat
    delta_time = last_exit_time - first_entry_time  #得到的是时间步长(单位: h)，要再转化为秒    
    
    # 将纬度转换为弧度（纬度要广播到72列）
    lat_rad = np.deg2rad(lat)  # (36,)
    lat_rad = np.tile(lat_rad, (num_lon, 1), )  # 复制为 (72, 36)

    # ✅ 纬度方向换算因子（单位：m/°）不变的那个
    # 不变
    deg_to_m_lat = np.full((num_lon, num_lat), np.deg2rad(1) * a)

    # ✅ 经度方向换算因子（单位：m/°）—— 与纬度相关
    # 变
    deg_to_m_lon = conver * np.cos(lat_rad)
    deg_to_m_lon = deg_to_m_lon[:, :, np.newaxis]         #拓展为3维 shape: (181, 46, 1)
    
    
    velocity_u = delta_lon/(delta_time*time_step)*deg_to_m_lon; velocity_v = delta_lat/(delta_time*time_step)*conver # 单位：°/s
    WRF_u = np.nansum(velocity_u, axis=2); WRF_v = np.nansum(velocity_v, axis=2)
    # WRF_u[WRF_u==0] = np.nan; WRF_v[WRF_v==0] = np.nan
    
    
    wave_propagation_time = np.nanmean(first_entry_time, axis=2) * time_step/3600.0/24 # 计算波列传播到对应网格的平均时间（单位：day）
    wave_propagation_u = np.nanmean(velocity_u, axis=2); wave_propagation_v = np.nanmean(velocity_v, axis=2) # 波列在对应网格的平均传播速度
    return(lon, lat, WRF_u, WRF_v, ray_count, wave_propagation_time, wave_propagation_u, wave_propagation_v)

# def WRF_calu(rlon_selected, rlat_selected, a, space, conver, time_step):
#     # 经度和纬度格点数
#     num_lon = int(360/space)+1
#     num_lat = int(180/space)+1
#     # 生成经度和纬度格点（中心点）
#     lon = np.linspace(0, 360, num_lon, endpoint=True)  # 经度（0到360度）
#     lat = np.linspace(-90, 90, num_lat, endpoint=True)  # 纬度（-90到90度）
#     # 设置网格
#     lon_bins = np.arange(0, 361, space)  # 经度每 space° 一格
#     lat_bins = np.arange(-90, 91, space)  # 纬度每 space° 一格
#     num_lon_bins = len(lon_bins) 
#     num_lat_bins = len(lat_bins) 
    
#     # 预分配 NaN 数组，最后一维增加1来记录时间索引
#     shape = (num_lon_bins, num_lat_bins, rlon_selected.shape[1])
#     first_entry_lon = np.full(shape, np.nan)
#     first_entry_lat = np.full(shape, np.nan)
#     last_exit_lon = np.full(shape, np.nan)
#     last_exit_lat = np.full(shape, np.nan)
    
#     # 新增数组记录首次进入和最后离开的时间步索引
#     first_entry_time = np.full(shape, np.nan)
#     last_exit_time = np.full(shape, np.nan)
#     # 遍历波射线
#     # 每条波射线（3,294,7）条都有一个切片（lon, lat）记录这条波射线在全球网格的进、出
#     for n in range(rlon_selected.shape[1]):  # 遍历波的索引
#         # 将位置索引映射到网格
#         lon_indices = np.clip(np.digitize(rlon_selected[:, n], lon_bins) - 1, 0, num_lon_bins - 1)
#         lat_indices = np.clip(np.digitize(rlat_selected[:, n], lat_bins) - 1, 0, num_lat_bins - 1)
    
#         visited_grid = {}  # 记录访问的网格
    
#         # 找到非 nan 的最后一个索引
#         valid_idx = np.where(~np.isnan(rlon_selected[:, n]))[0]
#         if len(valid_idx) > 1:
#             last_valid_idx = valid_idx[-1]  # 最后一个非nan索引
    
#             for i in range(last_valid_idx - 1):
#                 grid = (lon_indices[i], lat_indices[i])  # 当前时间步的位置在网格中的索引
#                 next_grid = (lon_indices[i + 1], lat_indices[i + 1])  # 下一个时间步的位置在网格中的索引
    
#                 # 🔹 记录首次进入
#                 if grid not in visited_grid:
#                     first_entry_lon[grid[0], grid[1], n] = rlon_selected[i, n]
#                     first_entry_lat[grid[0], grid[1], n] = rlat_selected[i, n]
#                     # ➡️ 记录时间步
#                     first_entry_time[grid[0], grid[1], n] = i
#                     visited_grid[grid] = True
    
#                 # 🔹 记录真正离开
#                 if grid != next_grid:
#                     last_exit_lon[grid[0], grid[1], n] = rlon_selected[i + 1, n]
#                     last_exit_lat[grid[0], grid[1], n] = rlat_selected[i + 1, n]
#                     # ➡️ 记录时间步
#                     last_exit_time[grid[0], grid[1], n] = i + 1
                    
#     delta_lon = last_exit_lon - first_entry_lon #要转换成实际的米！
#     delta_lat = last_exit_lat - first_entry_lat
#     delta_time = last_exit_time - first_entry_time  #得到的是时间步长(2h)，要再转化为秒    
#     '''
#     # 将纬度转换为弧度（纬度要广播到72列）
#     lat_rad = np.deg2rad(lat)  # (36,)
#     lat_rad = np.tile(lat_rad, (num_lon, 1))  # 复制为 (72, 36)

#     # ✅ 纬度方向换算因子（单位：m/°）不变的那个
#     # 不变
#     deg_to_m_lat = np.full((num_lon, num_lat), np.deg2rad(1) * a)

#     # ✅ 经度方向换算因子（单位：m/°）—— 与纬度相关
#     # 变
#     deg_to_m_lon = np.deg2rad(1) * a * np.cos(lat_rad)
#     '''
#     velocity_u = delta_lon/(delta_time*time_step)*conver; velocity_v = delta_lat/(delta_time*time_step)*conver # 单位：°/s
#     WRF_u = np.nansum(velocity_u, axis=2); WRF_v = np.nansum(velocity_v, axis=2)
#     WRF_u[WRF_u==0] = np.nan; WRF_v[WRF_v==0] = np.nan
#     return(WRF_u, WRF_v)



