# coding=utf-8
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch08 · 数字孪生：体元可视化平台启动器（VOXEL_PYTHON 可覆盖解释器）
# 状态  : 【可选工具留档】启动的是本地 PyVista GUI（非商业软件）；主链不依赖。
# 来源  : 竞赛提交包 3工程文件/Ch08_*/VoxelVizLauncher.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
# =====================================================================
# VoxelVizLauncher.py
# 启动器 —— 一键启动 PyQt 体元三维可视化平台
# ---------------------------------------------------------------------
# 路径自适应：
#   - 体元可视化主程序默认取本文件上一级目录的「体元可视化平台.py」
#   - Python 解释器默认当前解释器（sys.executable），可用环境变量 VOXEL_PYTHON 覆盖
#   - 主程序路径可用环境变量 VOXEL_SCRIPT / UI 参数 scriptpath 覆盖
# 调试日志写到系统临时目录 voxel_launcher_debug.log，运行异常可查。
# =====================================================================

import os
import subprocess
import sys
import traceback
import tempfile
from datetime import datetime

_HERE = os.path.dirname(os.path.abspath(__file__))

# 默认 Python 解释器（本机）；别人机器可设环境变量 VOXEL_PYTHON
DEFAULT_PYTHON = sys.executable
# 默认主程序：本启动器所在目录的上一级「体元可视化平台.py」
DEFAULT_SCRIPT = os.path.normpath(os.path.join(_HERE, "..", "体元可视化平台.py"))

DEBUG_LOG = os.path.join(tempfile.gettempdir(), "voxel_launcher_debug.log")
CREATE_NO_WINDOW = 0x08000000


def _log(msg):
    try:
        with open(DEBUG_LOG, "a", encoding="utf-8") as f:
            f.write("[%s] %s\n" % (datetime.now().strftime("%Y-%m-%d %H:%M:%S"), msg))
    except Exception:
        pass


class VoxelVizLauncher(object):
    """一键打开独立的体元可视化窗口。"""

    def execute(self, keyargs):
        _log("========== execute 开始 ==========")
        _log("python version: %s" % sys.version)
        _log("keyargs: %s" % repr(keyargs))

        try:
            # 主程序路径优先级：UI 参数 > 环境变量 > 默认相对路径
            script = (str(keyargs.get("scriptpath") or "")).strip() \
                or os.environ.get("VOXEL_SCRIPT", "").strip() \
                or DEFAULT_SCRIPT
            # Python 解释器优先级：环境变量 > 默认
            python = os.environ.get("VOXEL_PYTHON", "").strip() or DEFAULT_PYTHON

            _log("script path: %s" % script)
            _log("python path: %s" % python)
            _log("script exists? %s" % os.path.isfile(script))
            _log("python exists? %s" % os.path.isfile(python))

            if not os.path.isfile(script):
                return {"status": "ERROR: 主程序不存在 -> %s（请设置环境变量 VOXEL_SCRIPT 或检查路径）" % script}
            if not os.path.isfile(python):
                return {"status": "ERROR: Python 解释器不存在 -> %s（请安装依赖环境或设置环境变量 VOXEL_PYTHON）" % python}

            workdir = os.path.dirname(script)
            _log("workdir: %s" % workdir)

            cmd = [python, script]
            _log("cmd: %s" % repr(cmd))

            proc = subprocess.Popen(
                cmd,
                cwd=workdir,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                creationflags=CREATE_NO_WINDOW,
            )
            _log("subprocess started, pid=%d" % proc.pid)
            return {"status": "OK: 已启动体元可视化平台 (PID=%d)，调试日志: %s" % (proc.pid, DEBUG_LOG)}
        except Exception as e:
            tb = traceback.format_exc()
            _log("EXCEPTION: %s" % str(e))
            _log(tb)
            return {"status": "ERROR: 启动失败 -> %s" % str(e)}
