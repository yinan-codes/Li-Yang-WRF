# -*- coding: utf-8 -*-
"""
Created on Wed Mar 26 22:35:11 2025
2025.5.7更新筛选区域函数，加入波数制约条件
2025.7.23增加参数判断经过筛选区域是否进行截断
@author: 杨艺楠
"""

import numpy as np
def region_threshold(lon_rec1, lon_rec2, lat_rec1, lat_rec2,
                      rlon1, rlat1, wn_min, wn_max, cut=True):
    """
    保留满足以下条件的轨迹：
    - 如果进入过矩形框 2 次及以上，则在最后一次离开时截断，保留从起点到最后一次离开；
    - 如果只进入过一次，则采用 mask_out_of_box_after_entry() 的逻辑。
    同时返回：
    - 至少进入一次的轨迹数 entered_count
    - 重复进入（两次及以上）的轨迹数 multi_entry_count
    """
    xall, yall = [], []
    multi_entry_count = 0

    for loc in range(rlat1.shape[2]):
        for wn in range(wn_min - 1, wn_max):
            for j in range(3):
                x = rlon1[:, j, loc, wn]
                y = rlat1[:, j, loc, wn]

                if np.isnan(x[1:]).all() or np.isnan(y[1:]).all():
                    continue

                valid = ~np.isnan(x) & ~np.isnan(y)
                if not np.any(valid):
                    continue

                x_valid = x[valid]
                y_valid = y[valid]

                inside = (lon_rec1 <= x_valid) & (x_valid <= lon_rec2) & \
                         (lat_rec1 <= y_valid) & (y_valid <= lat_rec2)
                if not np.any(inside):
                    continue  # 从未进入框


                original_idx = np.where(valid)[0]
                shifted = np.roll(inside, 1)
                shifted[0] = False
                entry_flags = inside & (~shifted)
                exit_flags = (~inside) & shifted
                entry_indices = original_idx[entry_flags]
                exit_indices = original_idx[exit_flags]

                if len(entry_indices) >= 2:
                    multi_entry_count += 1  # ✅ 重复进入
                    if len(exit_indices) > 0 and exit_indices[-1] > entry_indices[-1]:
                        cut_idx = exit_indices[-1]
                        x_masked = x.copy()
                        y_masked = y.copy()
                        if cut:  # ✅ 只有当 cut=True 时才截断
                            x_masked[cut_idx+1:] = np.nan
                            y_masked[cut_idx+1:] = np.nan
                        # x_masked[cut_idx+1:] = np.nan
                        # y_masked[cut_idx+1:] = np.nan
                        xall.append(x_masked)
                        yall.append(y_masked)
                        continue
                    else:
                        xall.append(x.copy())
                        yall.append(y.copy())
                        continue

                # 只进入一次的情况
                x_masked, y_masked = mask_out_of_box_after_entry(
                    x.copy(), y.copy(), lon_rec1, lon_rec2, lat_rec1, lat_rec2, cut=cut # 输入到函数中cut的判断
                )
                valid_masked = ~np.isnan(x_masked) & ~np.isnan(y_masked)
                if not np.any(valid_masked):
                    continue
                
                if cut:
                    xi_last, yi_last = x_masked[valid_masked][-1], y_masked[valid_masked][-1]
                    inside_last = (lon_rec1 <= xi_last <= lon_rec2 and lat_rec1 <= yi_last <= lat_rec2)
                    if inside_last or np.count_nonzero(valid_masked) < np.count_nonzero(valid):
                        xall.append(x_masked)
                        yall.append(y_masked)
                else:
                    # ✅ 不截断的情况下，只要进入过一次就保留
                    xall.append(x_masked)
                    yall.append(y_masked)
          
                    

    xout = np.array(xall).T; entered_count = xout.shape[1]
    return np.array(xall).T, np.array(yall).T, entered_count, multi_entry_count

def mask_out_of_box_after_entry(x, y, lon_min, lon_max, lat_min, lat_max, cut=True):
    """
    轨迹一旦进入框内后，如果后续有出框的点，则从该点开始将轨迹设为 NaN。
    若轨迹未离开框，不进行截断。
    返回处理后的 x, y。
    """
    if not cut:
        return x, y

    valid = ~np.isnan(x) & ~np.isnan(y)
    if not np.any(valid):
        return x, y

    x_valid = x[valid]
    y_valid = y[valid]
    inside = (lon_min <= x_valid) & (x_valid <= lon_max) & \
              (lat_min <= y_valid) & (y_valid <= lat_max)

    if not np.any(inside):
        return x, y  # 从未进入框

    first_in = np.where(valid)[0][np.argmax(inside)]

    for idx in range(first_in, len(x)):
        xi, yi = x[idx], y[idx]
        if np.isnan(xi) or np.isnan(yi):
            continue
        if not (lon_min <= xi <= lon_max and lat_min <= yi <= lat_max):
            x[idx:] = np.nan
            y[idx:] = np.nan
            break  # 截断一次后即退出

    return x, y
# def region_threshold(lon_rec1, lon_rec2, lat_rec1, lat_rec2,
#                       rlon1, rlat1, wn_min, wn_max):
#     """
        # 只考虑进入一次的情况
#     保留满足以下任一条件的轨迹：
#     1. 曾进入框，后来离开（被截断后仍有有效点）；
#     2. 曾进入框，未离开（最后一点仍在框内）。
#     """
#     xall, yall = [], []

#     for loc in range(rlat1.shape[2]):
#         for wn in range(wn_min - 1, wn_max):
#             for j in range(3):
#                 x = rlon1[:, j, loc, wn]
#                 y = rlat1[:, j, loc, wn]

#                 if np.isnan(x[1:]).all() or np.isnan(y[1:]).all():
#                     continue

#                 valid = ~np.isnan(x) & ~np.isnan(y)
#                 if not np.any(valid):
#                     continue

#                 x_masked, y_masked = mask_out_of_box_after_entry(
#                     x.copy(), y.copy(),
#                     lon_rec1, lon_rec2, lat_rec1, lat_rec2
#                 )

#                 valid_masked = ~np.isnan(x_masked) & ~np.isnan(y_masked)
#                 if not np.any(valid_masked):
#                     continue

#                 # 判断是否是“进入但未离开”，需要检查最后一点是否仍在框内
#                 xi_last, yi_last = x_masked[valid_masked][-1], y_masked[valid_masked][-1]
#                 inside_last = (lon_rec1 <= xi_last <= lon_rec2 and lat_rec1 <= yi_last <= lat_rec2)

#                 # 只要最后还有点 & 曾进入过，就保留
#                 if inside_last or np.count_nonzero(valid_masked) < np.count_nonzero(valid):
#                     xall.append(x_masked)
#                     yall.append(y_masked)

#     return np.array(xall).T, np.array(yall).T




# def region_threshold(lon_rec1, lon_rec2, lat_rec1, lat_rec2, rlon1, rlat1, wn_min, wn_max):
#     xall, yall = [], []
#     for loc in range(rlat1.shape[2]):  # 22*13 个波源点
#         for wn in range(wn_min-1, wn_max):  # 纬向波数
#             for j in range(3):  # 每个解
#                 x, y = rlon1[:, j, loc, wn], rlat1[:, j, loc, wn]
#                 if np.isnan(x[1:]).all() or np.isnan(y[1:]).all():
#                     continue

#                 x_masked, y_masked = mask_out_of_box_after_entry(
#                     x.copy(), y.copy(), lon_rec1, lon_rec2, lat_rec1, lat_rec2
#                 )

#                 valid = ~np.isnan(x_masked) & ~np.isnan(y_masked)
#                 if not np.any(valid):
#                     continue

#                 # 检查最后一个有效点是否在框内
#                 xi_last, yi_last = x_masked[valid][-1], y_masked[valid][-1]
#                 if lon_rec1 <= xi_last <= lon_rec2 and lat_rec1 <= yi_last <= lat_rec2:
#                     xall.append(x_masked)
#                     yall.append(y_masked)

#     return np.array(xall).T, np.array(yall).T


# def mask_out_of_box_after_entry(x, y, lon_min, lon_max, lat_min, lat_max):
#     """
#     轨迹一旦进入框内后，如果后续有出框的点，则从该点开始将轨迹设为 NaN。
#     若轨迹最后一个有效点仍在框内，则不设 NaN。
#     返回处理后的 x, y
#     """
#     valid = ~np.isnan(x) & ~np.isnan(y)
#     if not np.any(valid):
#         return x, y  # 全是 NaN，直接返回

#     x_valid = x[valid]
#     y_valid = y[valid]
#     inside = (lon_min <= x_valid) & (x_valid <= lon_max) & \
#               (lat_min <= y_valid) & (y_valid <= lat_max)
    
#     if not np.any(inside):
#         return x, y  # 从未进入框，保留原轨迹（或你可选择丢弃）

#     # 进入框的第一个有效索引（对应原始 x/y 的索引）
#     first_in = np.where(valid)[0][np.argmax(inside)]

#     # 原始 x/y 中 first_in 之后的全部点（含）
#     for idx in range(first_in, len(x)):
#         xi, yi = x[idx], y[idx]
#         if np.isnan(xi) or np.isnan(yi):
#             continue
#         if not (lon_min <= xi <= lon_max and lat_min <= yi <= lat_max):
#             # 出框了，开始设为 nan
#             x[idx:] = np.nan
#             y[idx:] = np.nan
#             break  # 只截断一次

#     return x, y