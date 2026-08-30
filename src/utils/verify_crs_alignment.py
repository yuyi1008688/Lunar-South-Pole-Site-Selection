# ──────────────────────────────────────────────────────────────────
# 章节  : Utils · 跨章节工具：CRS/网格对齐验证（逐参数比对 transform，判断能否逐像元比较）
# 来源  : 竞赛提交包 3工程文件/Utils_*/verify_crs_alignment.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
CRS对齐验证脚本
验证gailv_240m.tif与AVGVISIB_EARTH的transform是否完全一致
"""
import os
import rasterio

PATH_LIGHT = os.environ.get('LUNAR_LIGHT_TIF', os.path.join('data', 'rasters', 'AVGVISIB_probability.tif'))
PATH_EARTH = os.environ.get('LUNAR_EARTH_TIF', 'avgvisib_65s_240m_earth.jp2')

print("=" * 60)
print("CRS对齐验证")
print("=" * 60)

with rasterio.open(PATH_LIGHT) as src1:
    with rasterio.open(PATH_EARTH) as src2:
        print(f"\n光照文件 (gailv_240m.tif):")
        print(f"  形状: {src1.shape}")
        print(f"  CRS: {src1.crs}")
        print(f"  Transform: {src1.transform}")
        print(f"  Bounds: {src1.bounds}")
        print(f"  Resolution: {src1.res}")
        
        print(f"\n通信文件 (AVGVISIB_EARTH):")
        print(f"  形状: {src2.shape}")
        print(f"  CRS: {src2.crs}")
        print(f"  Transform: {src2.transform}")
        print(f"  Bounds: {src2.bounds}")
        print(f"  Resolution: {src2.res}")
        
        print(f"\n{'='*60}")
        print(f"对齐检查:")
        print(f"  分辨率相同: {src1.res == src2.res}")
        print(f"  Transform完全一致: {src1.transform == src2.transform}")
        
        # 检查transform的每个参数
        t1 = src1.transform
        t2 = src2.transform
        print(f"\n  Transform参数对比:")
        print(f"    a (x分辨率): {t1.a} vs {t2.a} -> {'OK' if t1.a == t2.a else 'FAIL'}")
        print(f"    b (x旋转):   {t1.b} vs {t2.b} -> {'OK' if t1.b == t2.b else 'FAIL'}")
        print(f"    c (x原点):   {t1.c} vs {t2.c} -> {'OK' if t1.c == t2.c else ''}")
        print(f"    d (y旋转):   {t1.d} vs {t2.d} -> {'OK' if t1.d == t2.d else 'FAIL'}")
        print(f"    e (y分辨率): {t1.e} vs {t2.e} -> {'OK' if t1.e == t2.e else 'FAIL'}")
        print(f"    f (y原点):   {t1.f} vs {t2.f} -> {'OK' if t1.f == t2.f else 'FAIL'}")
        
        # 判断是否可以直接逐像元比较
        if src1.transform == src2.transform:
            print(f"\n  >>> 结论: Transform完全一致，可以直接逐像元比较")
        else:
            print(f"\n  >>> 结论: Transform不一致，不能直接逐像元比较！")
            print(f"  >>> 需要使用rasterio.warp.reproject或窗口读取进行空间对齐")
            
            # 检查是否只是原点不同但网格对齐
            if src1.res == src2.res:
                print(f"\n  分辨率相同但原点不同，检查网格是否对齐...")
                # 计算偏移量是否为分辨率的整数倍
                dx = abs(t1.c - t2.c)
                dy = abs(t1.f - t2.f)
                res_x = abs(t1.a)
                res_y = abs(t1.e)
                print(f"    X方向偏移: {dx:.2f}m = {dx/res_x:.4f}个像元")
                print(f"    Y方向偏移: {dy:.2f}m = {dy/res_y:.4f}个像元")
                
                if abs(dx/res_x - round(dx/res_x)) < 0.01 and abs(dy/res_y - round(dy/res_y)) < 0.01:
                    print(f"    偏移是像元的整数倍，网格对齐，可以用窗口读取")
                else:
                    print(f"    偏移不是像元的整数倍，网格未对齐，必须重采样！")
