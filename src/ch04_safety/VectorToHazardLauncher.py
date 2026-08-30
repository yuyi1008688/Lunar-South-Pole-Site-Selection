# -*- coding: utf-8 -*-
# ──────────────────────────────────────────────────────────────────
# 章节  : Ch04 · 安全势场：VectorToHazard 外部计算进程启动器
# 状态  : 【历史留档，主链已不使用】主分析链经函数级/直接脚本调用即可，
#         本启动器仅作为当年'外部进程隔离'工程化方案的留档保留。
# 来源  : 竞赛提交包 3工程文件/Ch04_*/VectorToHazardLauncher.py（算法逻辑保持原样，仅整理路径配置）
# 路径  : 已改为环境变量可覆盖；复现时请按文件内 docstring 说明准备输入数据
# ──────────────────────────────────────────────────────────────────
"""
VectorToHazard_launcher.py —— 外部计算进程启动器
"""
import os
import subprocess
import traceback

# ===== 配置（改你的实际路径） =====
PYTHON_EXE = os.environ.get("LUNAR_PYTHON", sys.executable)
MAIN_SCRIPT = os.environ.get("VECTORTOHAZARD_MAIN",
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "VectorToHazard_main.py"))
# ===================================

# ===== 调用参数定义 =====
# input, vector_path, string, 矢量数据路径,
# input, diam_field, string, 直径字段名,
# input, output_dir, string, 输出目录,
# output, status, string, 执行状态,


class VectorToHazardLauncher(object):
    def execute(self, keyargs):
        try:
            vector_path = keyargs.get("vector_path", "").strip()
            diam_field = keyargs.get("diam_field", "diameter_k").strip()
            output_dir = keyargs.get("output_dir", "").strip()

            if not vector_path:
                return {"status": "ERROR: 缺少 vector_path 参数"}
            if not os.path.isfile(vector_path):
                return {"status": f"ERROR: 矢量不存在: {vector_path}"}
            if not output_dir:
                output_dir = os.path.dirname(vector_path)
            if not os.path.isdir(output_dir):
                os.makedirs(output_dir, exist_ok=True)

            cmd = [PYTHON_EXE, MAIN_SCRIPT,
                   "--vector_path", vector_path,
                   "--output_dir", output_dir,
                   "--diam_field", diam_field]

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            if result.returncode != 0:
                return {"status": f"ERROR: 执行失败\n{result.stderr}"}

            return {"status": f"OK: {result.stdout.strip()}"}

        except subprocess.TimeoutExpired:
            return {"status": "ERROR: 执行超时(120秒)"}
        except Exception as e:
            return {"status": f"ERROR: {str(e)}\n{traceback.format_exc()}"}


if __name__ == "__main__":
    tool = VectorToHazardLauncher()
    print(tool.execute({
        "vector_path": "./craters_clean_study_PS.shp",
        "diam_field": "diameter_k",
        "output_dir": "./"
    }))