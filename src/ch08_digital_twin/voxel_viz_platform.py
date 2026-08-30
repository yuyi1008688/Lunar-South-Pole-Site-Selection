# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch08 · 数字孪生：PyQt5 + PyVista 交互式水冰体元可视化平台（阈值/配色/图层可调，--test 自检）
# 来源  : 竞赛提交包 3工程文件/Ch08_*/体元可视化平台.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
19_voxel_qt_app.py  水冰密度三维体元可视化平台 v6 (PyQt5 + pyvista)

基于 18_voxel_qt_app.py, 保留 v5 的磨砂玻璃胶囊图例方案, 微调:
  - 冰密度配色增加 viridis: YlOrRd / Blues / Greens / coolwarm(默认) / RdYlBu / viridis
  - 冰密度胶囊与面板腰线固定显示完整数据范围 [vmin, vmax], 不随阈值变化, 仅作展示范围用
  - 3D 视图中的冰体网格仍随阈值和配色实时变化
  - DEM 胶囊显示原始高程范围 [nanmin, nanmax](约 -4242 ~ +1954 m), 不再是平移后的 0 起
  - 18 保持原样不动

正常运行: python 19_voxel_qt_app.py
测试(不弹窗): python 19_voxel_qt_app.py --test
"""
import os
import sys
import numpy as np
import pyvista as pv
import vtk
import matplotlib
from PyQt5.QtWidgets import QWidget, QVBoxLayout
from PyQt5.QtGui import (QPainter, QColor, QPen, QFont, QLinearGradient,
                         QFontMetrics)
from PyQt5.QtCore import Qt, QRect

# 路径自适应：基于本脚本所在目录，数据放在同级 data/ 文件夹
_HERE = os.path.dirname(os.path.abspath(__file__))
VTK_PATH = os.path.join(_HERE, 'data', 'ice_voxel.vtk')
DEM_NPY = os.path.join(_HERE, 'data', 'dem_384.npy')
OUT_DIR = os.path.join(_HERE, 'outputs')
os.makedirs(OUT_DIR, exist_ok=True)

CMAP_TERRAIN = 'gist_earth'   # 地形配色(独立, 不抢冰的色)
# 冰密度配色选项(coolwarm 默认) — 已移除 jet/turbo, 新增 viridis
CMAP_OPTIONS = ['YlOrRd', 'Blues', 'Greens', 'coolwarm', 'RdYlBu', 'viridis']

# 统一基准(与master一致): 384x384 @ 240m, ±46080m
NX, NY = 384, 384
DX, DY = 240.0, 240.0
X0, Y0 = -46080.0, -46080.0

MODES = ['富集核(阈值)', '体渲染(半透明)', '等值面(立体)',
         '三向切片', '地形+冰壳', '埋藏深度(冰壳)']
MODE_CORE, MODE_VOLUME, MODE_ISO, MODE_SLICES, MODE_TERRAIN, MODE_BURIAL = 0, 1, 2, 3, 4, 5

BACKGROUNDS = {'暖白': ('#faf3e8', 'black'), '暖阳渐变': (('#f0c784', '#faf0dd'), 'black'),
               '纯白': ('white', 'black'), '深空黑': ('#0a0e1a', 'white'),
               '渐变蓝': (('#1a2a4a', '#0a0e1a'), 'white')}

# 光照参数(低环境光->明暗对比强->立体感强)
LIGHT = dict(smooth_shading=True, ambient=0.18, diffuse=0.85, specular=0.35,
             specular_power=25)


def load_grid():
    grid = pv.read(VTK_PATH)
    d = grid.cell_data['ice_density']
    valid = d[~np.isnan(d)]
    stats = {
        'vmin': float(valid.min()), 'vmax': float(valid.max()),
        'p50': float(np.percentile(valid, 50)),
        'p75': float(np.percentile(valid, 75)),
        'p90': float(np.percentile(valid, 90)),
        'n_valid': int(valid.size),
    }
    grid_vol = grid.copy()
    grid_vol.cell_data['ice_density'] = np.nan_to_num(d, nan=0.0)
    return grid, grid_vol, stats


def build_terrain():
    """
    返回:
      dem: StructuredGrid, 几何 Z 为平移后的正高程, 颜色标量 elevation 为原始高程
      dem_max_shift: 平移后的最大高程(用于定位冰壳悬浮高度)
      h_disp: 平移后的高程矩阵(2D)
      dem_min_orig, dem_max_orig: 原始 DEM 高程范围(用于图例)
    """
    h = np.load(DEM_NPY)
    dem_min_orig = float(np.nanmin(h))
    dem_max_orig = float(np.nanmax(h))
    # 几何需要非负, 因此平移; 但颜色保留原始高程, 使图例显示真实范围
    h_filled = np.nan_to_num(h, nan=dem_min_orig)
    h_disp = h_filled - dem_min_orig
    xc = X0 + (np.arange(NX) + 0.5) * DX
    yc = Y0 + (np.arange(NY) + 0.5) * DY
    gx, gy = np.meshgrid(xc, yc, indexing='ij')
    dem = pv.StructuredGrid()
    dem.points = np.column_stack([gx.ravel(), gy.ravel(), h_disp.ravel()])
    dem.dimensions = (NX, NY, 1)
    dem.point_data['elevation'] = h_filled.ravel()
    return dem, float(h_disp.max()), h_disp, dem_min_orig, dem_max_orig


def set_bg(plotter, bg_name):
    bg, txt_color = BACKGROUNDS[bg_name]
    if isinstance(bg, tuple):
        plotter.set_background(bg[0], top=bg[1])
    else:
        plotter.set_background(bg)
    return txt_color


# =================== VTK 色条 (仅 --test 预览用, 保留) ===================
def sb_ice(color, title='Ice density'):
    return {
        'title': title, 'color': color, 'n_labels': 6, 'fmt': '%.2f',
        'position_x': 0.50, 'position_y': 0.92, 'width': 0.44, 'height': 0.06,
        'vertical': False, 'title_font_size': 12, 'label_font_size': 10,
    }


def sb_ice_vertical(color, title='Ice density'):
    return {
        'title': title, 'color': color, 'n_labels': 6, 'fmt': '%.2f',
        'position_x': 0.88, 'position_y': 0.22, 'width': 0.045, 'height': 0.50,
        'vertical': True, 'title_font_size': 13, 'label_font_size': 11,
    }


def sb_dem_vertical(color):
    return {
        'title': 'DEM (m)', 'color': color, 'n_labels': 5, 'fmt': '%.0f',
        'position_x': 0.78, 'position_y': 0.22, 'width': 0.045, 'height': 0.50,
        'vertical': True, 'title_font_size': 13, 'label_font_size': 11,
    }


def style_legend(plotter, title, txt):
    try:
        actor = plotter.scalar_bars[title]
    except Exception:
        return
    try:
        is_dark = (txt == 'white')
        card = (0.05, 0.06, 0.09) if is_dark else (1.0, 1.0, 1.0)
        fg = (1.0, 1.0, 1.0) if is_dark else (0.08, 0.08, 0.10)
        op = 0.72 if is_dark else 0.90
        actor.DrawBackgroundOn()
        bp = actor.GetBackgroundProperty()
        bp.SetColor(*card)
        bp.SetOpacity(op)
        actor.DrawFrameOn()
        fp = actor.GetFrameProperty()
        fp.SetColor(*fg)
        fp.SetLineWidth(1.2)
        tp = actor.GetTitleTextProperty()
        tp.SetColor(*fg); tp.SetFontSize(15); tp.SetBold(1); tp.SetShadow(0)
        lp = actor.GetLabelTextProperty()
        lp.SetColor(*fg); lp.SetFontSize(11); lp.SetBold(0); lp.SetShadow(0)
    except Exception as e:
        print('style_legend skip:', e)


# ============================ 测试模式 ============================
def run_test():
    grid, grid_vol, stats = load_grid()
    dem, dem_max_shift, h_disp, dem_min_orig, dem_max_orig = build_terrain()
    vmax = stats['vmax']
    th = stats['p75']
    cases = [
        ('v19_core.png', MODE_CORE, th, 'coolwarm', '暖白', True),
        ('v19_volume.png', MODE_VOLUME, 0.0, 'viridis', '暖白', False),
        ('v19_iso.png', MODE_ISO, th, 'RdYlBu', '暖白', False),
        ('v19_slices.png', MODE_SLICES, 0.0, 'RdYlBu', '暖白', True),
        ('v19_terrain.png', MODE_TERRAIN, th, 'Greens', '暖白', False),
        ('v19_burial.png', MODE_BURIAL, th, 'Blues', '暖白', False),
    ]
    for fname, mode, t, cmap, bg, floor in cases:
        pl = pv.Plotter(off_screen=True, window_size=(1400, 1000))
        txt = set_bg(pl, bg)
        if mode == MODE_TERRAIN or mode == MODE_BURIAL:
            pl.add_mesh(dem.outline(), color=txt, line_width=1.2)
        else:
            pl.add_mesh(grid.outline(), color=txt, line_width=1.2)
        if mode == MODE_CORE:
            thm = grid.threshold(t)
            pl.add_mesh(thm, scalars='ice_density', cmap=cmap, clim=[t, vmax],
                        show_scalar_bar=True, scalar_bar_args=sb_ice(txt), **LIGHT)
            style_legend(pl, 'Ice density', txt)
        elif mode == MODE_VOLUME:
            pl.add_volume(grid_vol, scalars='ice_density', cmap=cmap,
                          clim=[0, stats['p90']], opacity='linear',
                          show_scalar_bar=True, scalar_bar_args=sb_ice(txt))
            style_legend(pl, 'Ice density', txt)
        elif mode == MODE_ISO:
            grid_pt = grid.cell_data_to_point_data()
            iso = grid_pt.contour(isosurfaces=[t], scalars='ice_density')
            pl.add_mesh(iso, scalars='ice_density', cmap=cmap, clim=[t, vmax],
                        show_scalar_bar=True, scalar_bar_args=sb_ice(txt), **LIGHT)
            style_legend(pl, 'Ice density', txt)
        elif mode == MODE_SLICES:
            sl = grid.slice_orthogonal()
            pl.add_mesh(sl, cmap=cmap, clim=[stats['vmin'], stats['p90']],
                        show_scalar_bar=True, scalar_bar_args=sb_ice(txt), **LIGHT)
            style_legend(pl, 'Ice density', txt)
            if floor:
                fl = grid.slice(normal='z', origin=(0, 0, 0.251))
                pl.add_mesh(fl, cmap=cmap, clim=[stats['vmin'], stats['p90']],
                            opacity=0.45, show_scalar_bar=False)
        elif mode == MODE_TERRAIN:
            pl.add_mesh(dem, scalars='elevation', cmap=CMAP_TERRAIN,
                        show_scalar_bar=True, scalar_bar_args=sb_dem_vertical(txt),
                        smooth_shading=True)
            style_legend(pl, 'DEM (m)', txt)
            ice = grid.threshold(t)
            pts = ice.points.copy()
            pts[:, 2] = pts[:, 2] * 1500 + (dem_max_shift + 2500)
            ice.points = pts
            pl.add_mesh(ice, scalars='ice_density', cmap=cmap, clim=[t, vmax],
                        show_scalar_bar=True, scalar_bar_args=sb_ice_vertical(txt), **LIGHT)
            style_legend(pl, 'Ice density', txt)
        else:  # 埋藏深度(冰壳)
            pl.add_mesh(dem, scalars='elevation', cmap=CMAP_TERRAIN, opacity=0.45,
                        show_scalar_bar=True, scalar_bar_args=sb_dem_vertical(txt),
                        smooth_shading=True)
            style_legend(pl, 'DEM (m)', txt)
            ice = grid.threshold(t)
            pts = ice.points.copy()
            ix = np.clip(((pts[:, 0] - X0) / DX + 0.5).astype(int), 0, NX - 1)
            iy = np.clip(((pts[:, 1] - Y0) / DY + 0.5).astype(int), 0, NY - 1)
            dem_val = h_disp[ix, iy]
            pts[:, 2] = dem_val - pts[:, 2] * 1500
            ice.points = pts
            pl.add_mesh(ice, scalars='ice_density', cmap=cmap, clim=[t, vmax],
                        show_scalar_bar=True, scalar_bar_args=sb_ice_vertical(txt), **LIGHT)
            style_legend(pl, 'Ice density', txt)
        if mode in (MODE_TERRAIN, MODE_BURIAL):
            pl.set_scale(zscale=1)
        else:
            pl.set_scale(zscale=1500)
        pl.camera_position = 'iso'
        pl.camera.azimuth = 35
        pl.camera.elevation = 22
        pl.screenshot(os.path.join(OUT_DIR, fname))
        pl.close()
        print('saved', fname)
    print('TEST OK')


# ============================ Qt 应用 ============================
QSS = """
QMainWindow, QWidget { background: #0d1117; color: #e6edf3;
    font-family: "Microsoft YaHei", "Segoe UI"; font-size: 13px; }
QFrame#panel { background: #161b22; border-right: 1px solid #30363d; }
QLabel#title { color: #58d6ff; font-size: 16px; font-weight: bold; }
QLabel#section { color: #8b949e; font-size: 11px; letter-spacing: 2px; }
QLabel#value { color: #58d6ff; font-weight: bold; }
QPushButton { background: #21262d; border: 1px solid #30363d; border-radius: 6px;
    padding: 9px; color: #e6edf3; }
QPushButton:hover { background: #2d333b; border-color: #58a6ff; }
QPushButton:pressed { background: #1f6feb; }
QComboBox { background: #21262d; border: 1px solid #30363d; border-radius: 6px;
    padding: 6px; color: #e6edf3; }
QComboBox QAbstractItemView { background: #21262d; color: #e6edf3;
    selection-background-color: #1f6feb; }
QSlider::groove:horizontal { height: 6px; background: #30363d; border-radius: 3px; }
QSlider::handle:horizontal { width: 16px; height: 16px; margin: -5px 0;
    border-radius: 8px; background: #58d6ff; }
QSlider::sub-page:horizontal { background: #1f6feb; border-radius: 3px; }
QCheckBox { color: #e6edf3; spacing: 8px; }
QStatusBar { background: #161b22; color: #8b949e; }
"""


def _cmap_gradient(name):
    try:
        cmap = matplotlib.colormaps[name]
    except Exception:
        import matplotlib.cm as mcm
        cmap = mcm.get_cmap(name)
    g = QLinearGradient(0, 0, 1, 0)
    g.setCoordinateMode(QLinearGradient.ObjectBoundingMode)
    for i in range(12):
        t = i / 11.0
        r, gg, b, _ = cmap(t)
        g.setColorAt(t, QColor(int(r * 255), int(gg * 255), int(b * 255)))
    return g


class LegendCapsule(QWidget):
    """磨砂玻璃圆角胶囊图例 (3D 视图左下角悬浮)"""

    def __init__(self, title, cmap, vmin, vmax, fmt='%.2f', parent=None):
        super().__init__(parent)
        self.title = title
        self._cmap = cmap
        self.vmin, self.vmax = vmin, vmax
        self.fmt = fmt
        self.setFixedSize(196, 46)
        self.setAttribute(Qt.WA_TranslucentBackground, True)

    def setCmap(self, c):
        self._cmap = c
        self.update()

    def setRange(self, a, b):
        self.vmin, self.vmax = a, b
        self.update()

    def paintEvent(self, e):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        rect = self.rect()
        p.setPen(Qt.NoPen)
        p.setBrush(QColor(16, 20, 28, 150))
        p.drawRoundedRect(rect, 13, 13)
        p.setBrush(Qt.NoBrush)
        p.setPen(QPen(QColor(255, 255, 255, 45), 1))
        p.drawRoundedRect(rect.adjusted(0, 0, -1, -1), 13, 13)
        p.setBrush(QColor(255, 255, 255, 18))
        p.setPen(Qt.NoPen)
        p.drawRoundedRect(rect.adjusted(2, 2, -2, -(rect.height() - 12)), 11, 11)
        p.setPen(QColor(232, 237, 242, 235))
        p.setFont(QFont('Microsoft YaHei', 9, QFont.Bold))
        p.drawText(13, 17, self.title)
        bx, by, bw, bh = 13, 21, self.width() - 26, 11
        p.setBrush(_cmap_gradient(self._cmap))
        p.setPen(Qt.NoPen)
        p.drawRoundedRect(QRect(bx, by, bw, bh), 5, 5)
        p.setFont(QFont('Microsoft YaHei', 8))
        p.setPen(QColor(220, 225, 230, 215))
        p.drawText(bx, self.height() - 4, self.fmt % self.vmin)
        fm = QFontMetrics(p.font())
        s = self.fmt % self.vmax
        p.drawText(self.width() - 13 - fm.horizontalAdvance(s), self.height() - 4, s)
        p.end()


class LegendWaist(QWidget):
    """右侧参数面板底部的装饰腰线(与面板等宽)"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._cmap = 'coolwarm'
        self.vmin, self.vmax = 0, 1
        self.setFixedHeight(48)
        self.setAttribute(Qt.WA_TranslucentBackground, True)

    def set_cmap(self, c):
        self._cmap = c
        self.update()

    def set_range(self, a, b):
        self.vmin, self.vmax = a, b
        self.update()

    def paintEvent(self, e):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        rect = self.rect()
        inner = rect.adjusted(6, 3, -6, -3)
        p.setPen(Qt.NoPen)
        p.setBrush(QColor(13, 17, 23, 165))
        p.drawRoundedRect(inner, 10, 10)
        p.setBrush(Qt.NoBrush)
        p.setPen(QPen(QColor(255, 255, 255, 28), 1))
        p.drawRoundedRect(inner, 10, 10)
        p.setPen(QColor(232, 237, 242, 215))
        p.setFont(QFont('Microsoft YaHei', 9, QFont.Bold))
        p.drawText(14, 18, 'Ice density')
        bx, by, bw, bh = 14, 24, self.width() - 28, 11
        p.setBrush(_cmap_gradient(self._cmap))
        p.setPen(Qt.NoPen)
        p.drawRoundedRect(QRect(bx, by, bw, bh), 5, 5)
        p.setFont(QFont('Microsoft YaHei', 8))
        p.setPen(QColor(220, 225, 230, 200))
        p.drawText(bx, self.height() - 4, '%.2f' % self.vmin)
        fm = QFontMetrics(p.font())
        s = '%.2f' % self.vmax
        p.drawText(self.width() - 14 - fm.horizontalAdvance(s), self.height() - 4, s)
        p.end()


class ViewArea(QWidget):
    """3D 视图容器: 承载 QtInteractor, 并可在其上叠加悬浮图例"""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._lay = QVBoxLayout(self)
        self._lay.setContentsMargins(0, 0, 0, 0)
        self.plotter = None
        self.overlays = []

    def setPlotter(self, p):
        self.plotter = p
        self._lay.addWidget(p)

    def addOverlay(self, w):
        w.setParent(self)
        w.setAttribute(Qt.WA_TransparentForMouseEvents, True)
        w.show()
        self.overlays.append(w)
        self._place()

    def _place(self):
        if not self.overlays:
            return
        m = 14
        w = self.width()
        h = self.height()
        x = m
        y = h - m
        for ov in reversed(self.overlays):
            ow, oh = ov.width(), ov.height()
            ov.setGeometry(x, y - oh, ow, oh)
            y -= (oh + 8)

    def resizeEvent(self, e):
        super().resizeEvent(e)
        self._place()


def run_app():
    from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QFrame,
                                 QHBoxLayout, QVBoxLayout, QLabel, QPushButton,
                                 QComboBox, QSlider, QCheckBox, QFileDialog)
    from PyQt5.QtCore import Qt, QTimer
    from pyvistaqt import QtInteractor

    grid, grid_vol, stats = load_grid()
    dem, dem_max_shift, h_disp, dem_min_orig, dem_max_orig = build_terrain()
    vmin, vmax = stats['vmin'], stats['vmax']

    class VoxelApp(QMainWindow):
        def __init__(self):
            super().__init__()
            self.setWindowTitle('月球南极水冰密度三维体元可视化平台 v6')
            self.resize(1500, 950)
            self._actors = []
            self._camera_init = False
            self._updating = False
            self._core_filter = None
            self._core_actor = None
            self._last = {}
            self._debounce = QTimer(self)
            self._debounce.setSingleShot(True)
            self._debounce.setInterval(100)
            self._debounce.timeout.connect(lambda: self.update_scene(from_slider=True))

            central = QWidget()
            self.setCentralWidget(central)
            lay = QHBoxLayout(central)
            lay.setContentsMargins(0, 0, 0, 0)
            lay.setSpacing(0)

            panel = QFrame()
            panel.setObjectName('panel')
            panel.setFixedWidth(290)
            pl_ = QVBoxLayout(panel)
            pl_.setContentsMargins(18, 18, 18, 18)
            pl_.setSpacing(10)

            t = QLabel('水冰密度体元\n三维可视化 v6')
            t.setObjectName('title')
            pl_.addWidget(t)

            pl_.addWidget(self._sec('显示模式'))
            self.cb_mode = QComboBox()
            self.cb_mode.addItems(MODES)
            self.cb_mode.currentIndexChanged.connect(self.update_scene)
            pl_.addWidget(self.cb_mode)

            pl_.addWidget(self._sec('密度阈值'))
            row = QHBoxLayout()
            self.sl_th = QSlider(Qt.Horizontal)
            self.sl_th.setRange(0, 1000)
            self.sl_th.setValue(int(1000 * (stats['p50'] - vmin) / (vmax - vmin)))
            self.sl_th.valueChanged.connect(self._on_th_slider)
            self.lb_th = QLabel()
            self.lb_th.setObjectName('value')
            row.addWidget(self.sl_th, 1)
            row.addWidget(self.lb_th)
            pl_.addLayout(row)

            pl_.addWidget(self._sec('Z轴夸张 / 埋藏夸张'))
            row2 = QHBoxLayout()
            self.sl_z = QSlider(Qt.Horizontal)
            self.sl_z.setRange(500, 5000)
            self.sl_z.setSingleStep(100)
            self.sl_z.setValue(1500)
            self.sl_z.valueChanged.connect(self._on_z_slider)
            self.lb_z = QLabel('x1500')
            self.lb_z.setObjectName('value')
            row2.addWidget(self.sl_z, 1)
            row2.addWidget(self.lb_z)
            pl_.addLayout(row2)

            pl_.addWidget(self._sec('冰密度配色'))
            self.cb_cmap = QComboBox()
            self.cb_cmap.addItems(CMAP_OPTIONS)
            self.cb_cmap.setCurrentText('coolwarm')
            self.cb_cmap.currentIndexChanged.connect(self.update_scene)
            pl_.addWidget(self.cb_cmap)

            pl_.addWidget(self._sec('背景'))
            self.cb_bg = QComboBox()
            self.cb_bg.addItems(list(BACKGROUNDS.keys()))
            self.cb_bg.setCurrentText('暖白')
            self.cb_bg.currentIndexChanged.connect(self.update_scene)
            pl_.addWidget(self.cb_bg)

            self.ck_floor = QCheckBox('显示底部密度平面')
            self.ck_floor.setChecked(True)
            self.ck_floor.stateChanged.connect(self.update_scene)
            pl_.addWidget(self.ck_floor)

            self.lb_hint = QLabel()
            self.lb_hint.setWordWrap(True)
            self.lb_hint.setObjectName('section')
            pl_.addWidget(self.lb_hint)

            btn = QPushButton('保存当前视角截图 (PNG)')
            btn.clicked.connect(self.save_shot)
            pl_.addWidget(btn)

            pl_.addStretch(1)
            self.lb_info = QLabel()
            self.lb_info.setWordWrap(True)
            pl_.addWidget(self.lb_info)

            self.waist = LegendWaist()
            pl_.addWidget(self.waist)

            lay.addWidget(panel)

            self.view = ViewArea()
            lay.addWidget(self.view, 1)
            self.plotter = QtInteractor(self.view)
            self.view.setPlotter(self.plotter)

            self.cap_ice = LegendCapsule('Ice density', 'coolwarm', vmin, vmax)
            self.cap_dem = LegendCapsule('DEM (m)', 'gist_earth',
                                         dem_min_orig, dem_max_orig, fmt='%.0f')
            self.view.addOverlay(self.cap_ice)
            self.view.addOverlay(self.cap_dem)
            self.cap_dem.hide()

            self.statusBar().showMessage('就绪')
            self.update_scene()

        def _sec(self, text):
            lb = QLabel(text)
            lb.setObjectName('section')
            return lb

        @property
        def threshold(self):
            return vmin + (vmax - vmin) * self.sl_th.value() / 1000.0

        def _on_th_slider(self, _):
            self.lb_th.setText(f'{self.threshold:.3f}')
            self._debounce.start()

        def _on_z_slider(self, _):
            self.lb_z.setText(f'x{self.sl_z.value()}')
            self._debounce.start()

        def _mesh_ice_range(self, mode, th):
            """返回 3D 视图里冰体网格实际使用的 scalar_range"""
            if mode == MODE_VOLUME:
                return 0.0, stats['p90']
            if mode == MODE_SLICES:
                return vmin, stats['p90']
            return th, vmax

        def _refresh_labels(self, th, n_shown, cmap):
            self.lb_info.setText(
                f"冰密度配色: {cmap}\n"
                f"体元总量: {grid.n_cells:,}\n"
                f"有效体元: {stats['n_valid']:,}\n"
                f"P50 / P75 / P90:\n{stats['p50']:.3f} / {stats['p75']:.3f} / {stats['p90']:.3f}\n"
                f"当前显示: {n_shown:,} 体元")
            self.statusBar().showMessage(
                f"{self.cb_mode.currentText()} | 阈值 {th:.3f} | 配色 {cmap} | "
                f"显示 {n_shown:,} 体元 | 左键旋转 滚轮缩放")

        def update_scene(self, from_slider=False):
            if self._updating:
                return
            self._updating = True
            try:
                mode = self.cb_mode.currentIndex()
                bg_name = self.cb_bg.currentText()
                floor_on = self.ck_floor.isChecked()
                z = self.sl_z.value()
                th = self.threshold
                cmap = self.cb_cmap.currentText()
                self.lb_th.setText(f'{th:.3f}')
                self.lb_z.setText(f'x{z}')
                self.sl_th.setEnabled(mode in (MODE_CORE, MODE_ISO, MODE_TERRAIN, MODE_BURIAL))
                self.ck_floor.setEnabled(mode not in (MODE_TERRAIN, MODE_BURIAL))
                if mode == MODE_BURIAL:
                    self.lb_hint.setText('半透明地表之下即冰体, 其与地表的缝隙=夸张后的埋藏深度')
                else:
                    self.lb_hint.setText('')

                params = dict(mode=mode, bg=bg_name, floor=floor_on, cmap=cmap)
                last = self._last

                if (from_slider and mode == MODE_CORE
                        and self._core_actor is not None
                        and last.get('params') == params
                        and last.get('z') == z):
                    self._core_filter.SetThresholdFunction(
                        vtk.vtkThreshold.THRESHOLD_UPPER)
                    self._core_filter.SetUpperThreshold(th)
                    self._core_filter.Update()
                    self._core_actor.mapper.scalar_range = (th, vmax)
                    n_shown = int(self._core_filter.GetOutput().GetNumberOfCells())
                    last['th'] = th
                    # 图例始终展示完整范围
                    self.cap_ice.setCmap(cmap)
                    self.cap_ice.setRange(vmin, vmax)
                    self.waist.set_cmap(cmap)
                    self.waist.set_range(vmin, vmax)
                    self._refresh_labels(th, n_shown, cmap)
                    return

                for a in self._actors:
                    self.plotter.remove_actor(a)
                self._actors = []
                self._core_actor = None

                txt = set_bg(self.plotter, bg_name)
                if mode in (MODE_TERRAIN, MODE_BURIAL):
                    a_out = self.plotter.add_mesh(dem.outline(), color=txt,
                                                  line_width=1.5)
                else:
                    a_out = self.plotter.add_mesh(grid.outline(), color=txt,
                                                  line_width=1.5)
                self._actors.append(a_out)

                ice_clim = self._mesh_ice_range(mode, th)

                if mode == MODE_CORE:
                    if self._core_filter is None:
                        f = vtk.vtkThreshold()
                        f.SetInputData(grid)
                        f.SetInputArrayToProcess(
                            0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS,
                            'ice_density')
                        self._core_filter = f
                    self._core_filter.SetThresholdFunction(
                        vtk.vtkThreshold.THRESHOLD_UPPER)
                    self._core_filter.SetUpperThreshold(th)
                    self._core_actor = self.plotter.add_mesh(
                        self._core_filter, scalars='ice_density', cmap=cmap,
                        clim=ice_clim, show_scalar_bar=False, **LIGHT)
                    self._actors.append(self._core_actor)
                    n_shown = int(self._core_filter.GetOutput().GetNumberOfCells())
                    if floor_on:
                        fl = grid.slice(normal='z', origin=(0, 0, 0.251))
                        self._actors.append(self.plotter.add_mesh(
                            fl, cmap=cmap, clim=[vmin, stats['p90']],
                            opacity=0.45, show_scalar_bar=False))

                elif mode == MODE_VOLUME:
                    self._actors.append(self.plotter.add_volume(
                        grid_vol, scalars='ice_density', cmap=cmap,
                        clim=ice_clim, opacity='linear', show_scalar_bar=False))
                    n_shown = stats['n_valid']

                elif mode == MODE_ISO:
                    grid_pt = grid.cell_data_to_point_data()
                    iso = grid_pt.contour(isosurfaces=[th], scalars='ice_density')
                    self._actors.append(self.plotter.add_mesh(
                        iso, scalars='ice_density', cmap=cmap, clim=ice_clim,
                        show_scalar_bar=False, **LIGHT))
                    n_shown = iso.n_cells

                elif mode == MODE_SLICES:
                    sl = grid.slice_orthogonal()
                    self._actors.append(self.plotter.add_mesh(
                        sl, cmap=cmap, clim=ice_clim,
                        show_scalar_bar=False, **LIGHT))
                    n_shown = stats['n_valid']
                    if floor_on:
                        fl = grid.slice(normal='z', origin=(0, 0, 0.251))
                        self._actors.append(self.plotter.add_mesh(
                            fl, cmap=cmap, clim=[vmin, stats['p90']],
                            opacity=0.45, show_scalar_bar=False))

                elif mode == MODE_TERRAIN:
                    self._actors.append(self.plotter.add_mesh(
                        dem, scalars='elevation', cmap=CMAP_TERRAIN,
                        show_scalar_bar=False, smooth_shading=True))
                    ice = grid.threshold(th)
                    pts = ice.points.copy()
                    pts[:, 2] = pts[:, 2] * z + (dem_max_shift + 2500)
                    ice.points = pts
                    self._actors.append(self.plotter.add_mesh(
                        ice, scalars='ice_density', cmap=cmap, clim=ice_clim,
                        show_scalar_bar=False, **LIGHT))
                    n_shown = ice.n_cells

                else:  # 埋藏深度(冰壳)
                    self._actors.append(self.plotter.add_mesh(
                        dem, scalars='elevation', cmap=CMAP_TERRAIN, opacity=0.45,
                        show_scalar_bar=False, smooth_shading=True))
                    ice = grid.threshold(th)
                    pts = ice.points.copy()
                    ix = np.clip(((pts[:, 0] - X0) / DX + 0.5).astype(int), 0, NX - 1)
                    iy = np.clip(((pts[:, 1] - Y0) / DY + 0.5).astype(int), 0, NY - 1)
                    dem_val = h_disp[ix, iy]
                    pts[:, 2] = dem_val - pts[:, 2] * z
                    ice.points = pts
                    self._actors.append(self.plotter.add_mesh(
                        ice, scalars='ice_density', cmap=cmap, clim=ice_clim,
                        show_scalar_bar=False, **LIGHT))
                    n_shown = ice.n_cells

                if mode in (MODE_TERRAIN, MODE_BURIAL):
                    self.plotter.set_scale(zscale=1)
                else:
                    self.plotter.set_scale(zscale=z)
                if not self._camera_init:
                    self.plotter.camera_position = 'iso'
                    self.plotter.camera.azimuth = 35
                    self.plotter.camera.elevation = 22
                    self._camera_init = True

                # ===== 自绘图例更新: 冰密度始终展示完整范围; DEM 展示原始高程范围 =====
                self.cap_ice.setCmap(cmap)
                self.cap_ice.setRange(vmin, vmax)
                self.cap_ice.show()
                self.waist.set_cmap(cmap)
                self.waist.set_range(vmin, vmax)
                if mode in (MODE_TERRAIN, MODE_BURIAL):
                    self.cap_dem.setRange(dem_min_orig, dem_max_orig)
                    self.cap_dem.setCmap(CMAP_TERRAIN)
                    self.cap_dem.show()
                else:
                    self.cap_dem.hide()
                self.view._place()

                self._last = {'params': params, 'th': th, 'z': z}
                self._refresh_labels(th, n_shown, cmap)
            finally:
                self._updating = False

        def save_shot(self):
            path, _ = QFileDialog.getSaveFileName(
                self, '保存截图', os.path.join(OUT_DIR, 'voxel_shot.png'),
                'PNG (*.png)')
            if path:
                self.plotter.screenshot(path)
                self.statusBar().showMessage(f'已保存: {path}')

    app = QApplication(sys.argv)
    app.setStyleSheet(QSS)
    win = VoxelApp()
    win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    if '--test' in sys.argv:
        run_test()
    else:
        run_app()
