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
    
    # 3D 记录（保留，不再用它们做“首入-末出”求差）
    shape = (num_lon_bins, num_lat_bins, rlon_selected.shape[1])
    first_entry_lon = np.full(shape, np.nan)
    first_entry_lat = np.full(shape, np.nan)
    last_exit_lon   = np.full(shape, np.nan)
    last_exit_lat   = np.full(shape, np.nan)
    first_entry_time = np.full(shape, np.nan)
    last_exit_time   = np.full(shape, np.nan)

    # —— 轻量“每次经过”的容器（只存索引与首末信息）——
    pass_ix, pass_iy = [], []
    pass_fe_lon, pass_fe_lat, pass_fe_time = [], [], []
    pass_le_lon, pass_le_lat, pass_le_time = [], [], []

    # 通过次数（按“每次经过”计数）
    ray_count = np.zeros((num_lon_bins, num_lat_bins), dtype=np.int32)

    # “该格点各射线的第一次进入时间”（用于 wave_propagation_time）
    first_entry_time_first = np.full((num_lon_bins, num_lat_bins, rlon_selected.shape[1]), np.nan)

    # —— 遍历每条射线 —— 
    for n in range(rlon_selected.shape[1]):
        # 将位置索引映射到网格
        lon_indices = np.clip(np.digitize(rlon_selected[:, n], lon_bins) - 1, 0, num_lon_bins - 1)
        lat_indices = np.clip(np.digitize(rlat_selected[:, n], lat_bins) - 1, 0, num_lat_bins - 1)

        # 当前射线：各格子的“本次进入”暂存
        cur_entry_lon2  = np.full((num_lon_bins, num_lat_bins), np.nan)
        cur_entry_lat2  = np.full((num_lon_bins, num_lat_bins), np.nan)
        cur_entry_time2 = np.full((num_lon_bins, num_lat_bins), np.nan)

        # 有效索引范围
        valid = ~np.isnan(rlon_selected[:, n]) & ~np.isnan(rlat_selected[:, n])
        idx_valid = np.where(valid)[0]
        if idx_valid.size < 2:
            continue
        last_valid_idx = int(idx_valid[-1])

        for i in range(last_valid_idx - 1):
            if not (valid[i] and valid[i+1]):
                # 断点：丢弃未闭合的一次通过
                continue

            grid      = (lon_indices[i],   lat_indices[i])
            next_grid = (lon_indices[i+1], lat_indices[i+1])

            # —— 首次进入（旧逻辑保留，用于记录“该射线首次进此格”的信息）——
            if np.isnan(first_entry_time[grid[0], grid[1], n]):
                first_entry_lon[grid[0], grid[1], n]  = rlon_selected[i, n]
                first_entry_lat[grid[0], grid[1], n]  = rlat_selected[i, n]
                first_entry_time[grid[0], grid[1], n] = i

            # —— 本次进入（用于“每次经过”的切片记录）——
            if np.isnan(cur_entry_time2[grid[0], grid[1]]):
                cur_entry_lon2[grid[0], grid[1]]  = rlon_selected[i, n]
                cur_entry_lat2[grid[0], grid[1]]  = rlat_selected[i, n]
                cur_entry_time2[grid[0], grid[1]] = i

            # —— 记录“该射线对该格”的第一次进入时间（用于平均到达时间）——
            if np.isnan(first_entry_time_first[grid[0], grid[1], n]):
                first_entry_time_first[grid[0], grid[1], n] = i

            # —— 检测离开：一旦离开当前格，立刻把这一次通过追加到轻量列表 —— 
            if grid != next_grid:
                # 更新“最后离开”的 3D 记录（保留旧接口）
                last_exit_lon[grid[0], grid[1], n]  = rlon_selected[i + 1, n]
                last_exit_lat[grid[0], grid[1], n]  = rlat_selected[i + 1, n]
                last_exit_time[grid[0], grid[1], n] = i + 1

                # 若本次进入有效，则形成一段通过
                if not np.isnan(cur_entry_time2[grid[0], grid[1]]):
                    gx, gy = grid
                    pass_ix.append(gx); pass_iy.append(gy)
                    pass_fe_lon.append(cur_entry_lon2[gx, gy])
                    pass_fe_lat.append(cur_entry_lat2[gx, gy])
                    pass_fe_time.append(cur_entry_time2[gx, gy])
                    pass_le_lon.append(rlon_selected[i + 1, n])
                    pass_le_lat.append(rlat_selected[i + 1, n])
                    pass_le_time.append(i + 1)

                    # 按“每次经过”计数
                    ray_count[gx, gy] += 1

                    # 清空该格的“本次进入”暂存，等待后续可能的重入
                    cur_entry_lon2[gx, gy]  = np.nan
                    cur_entry_lat2[gx, gy]  = np.nan
                    cur_entry_time2[gx, gy] = np.nan

                # 不预设 next_grid 的进入；下一循环会自然处理

    # —— 用轻量“每次经过”记录生成 WRF —— 
    if len(pass_ix) == 0:
        WRF_u = np.zeros((num_lon_bins, num_lat_bins), dtype=np.float64)
        WRF_v = np.zeros_like(WRF_u)
        wave_propagation_u = np.full_like(WRF_u, np.nan)
        wave_propagation_v = np.full_like(WRF_v, np.nan)
    else:
        ix = np.asarray(pass_ix, dtype=np.int32)
        iy = np.asarray(pass_iy, dtype=np.int32)
        fe_lon = np.asarray(pass_fe_lon, dtype=np.float64)
        fe_lat = np.asarray(pass_fe_lat, dtype=np.float64)
        fe_t   = np.asarray(pass_fe_time, dtype=np.float64)
        le_lon = np.asarray(pass_le_lon, dtype=np.float64)
        le_lat = np.asarray(pass_le_lat, dtype=np.float64)
        le_t   = np.asarray(pass_le_time, dtype=np.float64)

        dt   = (le_t - fe_t) * time_step               # 秒
        dlam = (le_lon - fe_lon)                       # °（经度差）
        dphi = (le_lat - fe_lat)                       # °（纬度差）
        phiM = np.deg2rad(0.5*(le_lat + fe_lat))       # 段中点纬度（弧度）

        with np.errstate(invalid='ignore', divide='ignore'):
            # conver = (π/180)*a，单位 m/°
            u_seg = (dlam/dt) * (conver * np.cos(phiM))  # m/s
            v_seg = (dphi/dt) * (conver)                 # m/s

        # 散点累加到网格（与 np.nansum(axis=2) 的结果等价）
        WRF_u = np.zeros((num_lon_bins, num_lat_bins), dtype=np.float64)
        WRF_v = np.zeros_like(WRF_u)
        np.add.at(WRF_u, (ix, iy), u_seg)
        np.add.at(WRF_v, (ix, iy), v_seg)

        # 每格“每次经过的平均速度”（可视化/诊断用）
        with np.errstate(invalid='ignore'):
            wave_propagation_u = np.where(ray_count>0, WRF_u/np.maximum(ray_count,1), np.nan)
            wave_propagation_v = np.where(ray_count>0, WRF_v/np.maximum(ray_count,1), np.nan)

    # 平均到达时间（天）：按“各射线第一次进入该格”的平均
    wave_propagation_time = np.nanmean(first_entry_time_first, axis=2) * (time_step/3600.0/24)

    return (lon, lat, WRF_u, WRF_v, ray_count, wave_propagation_time, wave_propagation_u, wave_propagation_v)


