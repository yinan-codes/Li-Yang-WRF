# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 21:44:42 2025
设置速度阈值并剔除异常值
@author: 杨艺楠
"""
import netCDF4 as nc
import numpy as np

# def load_wave_ray_data(path):
#     """
#     从指定路径的 NetCDF 文件中读取波射线变量：
#     rlon、rlat、rzwn、rmwn，并进行 squeeze 处理。

#     参数：
#         path (str): NetCDF 文件路径

#     返回：
#         rlon, rlat, rzwn, rmwn: 全部为 ndarray
#     """
#     with nc.Dataset(path) as ds:
#         rlon = ds.variables['rlon'][:].squeeze()
#         rlat = ds.variables['rlat'][:].squeeze()
#         rzwn = ds.variables['rzwn'][:].squeeze()
#         rmwn = ds.variables['rmwn'][:].squeeze()
    
#     return rlon, rlat, rzwn, rmwn
def load_wave_ray_data(path):
    with nc.Dataset(path) as ds:
        rlon = ds.variables['rlon'][:].squeeze()
        rlat = ds.variables['rlat'][:].squeeze()
        rmwn = ds.variables['rmwn'][:].squeeze()
    
    return rlon, rlat, rmwn
#==============================================================================
def threshold(rlon0, rlat0, wn_min, wn_max, conver, time_step, velo_threshold, *, # *表示从此处开始只能关键字传参，不得位置传参
              rzwn=None, rmwn=None, check_wn=False, wn_threshold=None, L=None):
              # wn_mode="abs",               # "abs" 或 "ratio"
              # wn_threshold=None, N=None,   # 仅 abs 模式使用（|Δ²x_t| > wn_threshold）
              # wn_ratio_threshold=None):    # 仅 ratio 模式使用（|Δ²x_t|/均幅 > wn_ratio_threshold）
    """
    N = L/h: 连续不稳定的时间步数
    # wn_mode:
    #   - "abs"   : 数值制约 -> |Δ²x_t| > wn_threshold
    #   - "ratio" : 比率制约 -> |Δ²x_t| / mean(|x_{t-1}|,|x_t|,|x_{t+1}|) > wn_ratio_threshold

    # 注意：check_wn=False 时不做任何波数制约。
    """
    # # —— 预先准备：有关经向波数制约的参数校验 ——
    # k_steps = 3 if (N is None) else int(N)
    # if k_steps < 1:
    #     raise ValueError("N(=连续步数) 必须是 >=1 的整数")
    # if wn_mode not in ("abs", "ratio"):
    #     raise ValueError("wn_mode 只能是 'abs' 或 'ratio'")
    # if check_wn:
    #     if (rzwn is None) or (rmwn is None):
    #         raise ValueError("check_wn=True 需提供 rzwn 与 rmwn")
    #     if wn_mode == "abs" and (wn_threshold is None):
    #         raise ValueError("wn_mode='abs' 需提供 wn_threshold（绝对阈值）")
    #     if wn_mode == "ratio" and (wn_ratio_threshold is None):
    #         raise ValueError("wn_mode='ratio' 需提供 wn_ratio_threshold（比率阈值）")
            
    
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
                    # —— 波数异常（可选） ——
                if not check_wn:
                    continue
    
                cut_done = False
                for var in (rzwn, rmwn):
                    if cut_done:
                        break
                    if var is None:
                        continue
                    series = var[:, j, loc, wn]  # shape (T,)
                    if series.size < 3 or np.isnan(series).all():
                        continue
    
                    d2 = series[2:] - 2*series[1:-1] + series[:-2]   # 长度 T-2
    
                    # if wn_mode == "abs":
                    #     metric = np.abs(d2)  # |Δ²x_t|
                    #     thr = wn_threshold
                    # else:  # "ratio"
                    #     den = (np.abs(series[2:]) + np.abs(series[1:-1]) + np.abs(series[:-2])) / 3.0
                    #     with np.errstate(divide='ignore', invalid='ignore'):
                    #         # metric = np.abs(d2) / den   # 不加 τ，den=0 时可能得到 inf
                    #         ratio_tau_alpha = 0.1; ax = np.abs(series[np.isfinite(series)]) # —— ratio：|Δ²x_t| / (三点均幅 + τ)，τ = α * median(|x|)
                    #         tau = ratio_tau_alpha * (np.median(ax) if ax.size else 0.0)
                    #         metric = np.abs(d2) / (den + tau)
                    #     thr = wn_ratio_threshold
    
                    # mask = (metric > thr).astype(np.int8)
                    # if mask.size >= k_steps:
                    #     hit = np.convolve(mask, np.ones(k_steps, dtype=np.int8), mode='valid')
                    #     idx = np.flatnonzero(hit == k_steps)
                    #     if idx.size > 0:
                    #         cut_idx = int(idx[0] + (k_steps - 1) + 2)  # d2→原时轴 +2
                    #         rlon0[cut_idx:, j, loc, wn] = np.nan
                    #         rlat0[cut_idx:, j, loc, wn] = np.nan
                    #         cut_done = True

                    
                # ✅ 可选：波数异常判断
                if check_wn and rzwn is not None and rmwn is not None and wn_threshold is not None:
                    N = int(np.ceil((L*3600) / time_step))
                    for var in [rzwn, rmwn]:
                        series = var[:, j, loc, wn]  # shape (T,)
                        if np.isnan(series).all():
                            continue

                        d2 = series[2:] - 2 * series[1:-1] + series[:-2] # 进行二阶差分
                        mask = (np.abs(d2) > wn_threshold).astype(np.int8)
                        hit = np.convolve(mask, np.ones(N, dtype=np.int8), mode='valid')  # 连续 True 计数
                        idx = np.flatnonzero(hit == N)
                        if idx.size > 0:
                            cut_idx = idx[0] + (N - 1) + 2   # +2 把 d2 轴映射回原时间轴
                            rlon0[cut_idx:, j, loc, wn] = np.nan
                            rlat0[cut_idx:, j, loc, wn] = np.nan
                            # rzwn[cut_idx:, j, loc, wn] = np.nan # for tests
                            # rmwn[cut_idx:, j, loc, wn] = np.nan # for tests
                            break

    
    return(rlon0, rlat0)