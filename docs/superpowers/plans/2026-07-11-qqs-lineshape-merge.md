# qqs/lineshape/ 光致发光算法统一 实施计划

> **给 agentic workers：** 必选子技能：使用 superpowers:subagent-driven-development（推荐）或 superpowers:executing-plans 按任务逐步实施。步骤使用复选框（`- [ ]`）语法跟踪。

**目标：** 把 `qqs/lineshape/src/photonics2/photoluminescence.py`、`qqs/lineshape/src/photonics2/plott.py`、`qqs/lineshape/src/pl.py` 三个差异文件统一到单一版本（来源是 `lineshape_new_ref/` + 主版决策），通过 `run_compare.py` 验证数值与 `BASELINE.md` 一致，通过 `python src/pl.py` 冒烟测试产生 6 张 PNG。

**Spec：** `docs/superpowers/specs/2026-07-11-qqs-lineshape-merge-design.md`

**架构：** 三个文件的备份与替换 + `run_compare.py` 回归 + `pl.py` 冒烟测试 + 文档更新。无 API 变更，不删 `lineshape_new_ref/`，不删旧版独占模块。

**技术栈：** shell（`cp`、`sed`、`rm`）、Python（运行 `run_compare.py`、`python src/pl.py`），无需安装新包。

---

## 全局约束

- `lineshape_new_ref/` 子包的 `photoluminescence.py`、`plott.py` 是对应文件的替换来源；`pl.py` 的合并源是 spec 决策表
- 三个被替换文件在替换前各自留一份 `<file>.bak`，验证通过后删除
- `src/INCAR` 在 `pl.py` 冒烟测试期间临时改 `resolution` 为 `1000`，测试结束用 `git checkout -- src/INCAR` 恢复
- 不改 `phonon_struct_op.py`、`cases/`、`output/`、`run_compare.py`、`BASELINE.md`
- 不删任何旧版独占模块（`jt.py`、`jtsoc.py`、`ht.py`、`schrodinger.py`、`embedding.py`、`ebedding.py`、`hermite.py`、`nu.py`、`un.py`、`constants.py`）

---

### Task 1：验证三个替换源都到位

**涉及文件：**
- 读取：`qqs/lineshape/src/lineshape_new_ref/photonics2/photoluminescence.py`
- 读取：`qqs/lineshape/src/lineshape_new_ref/photonics2/plott.py`
- 读取：`qqs/lineshape/src/pl.py`（主版，作为合并基线参考）

**接口：**
- 进入：进入实施前的最后一步只读检查
- 产出：确认三个源文件存在且非空

- [ ] **Step 1：确认两个 `lineshape_new_ref` 源存在并记录行数**

```bash
test -f qqs/lineshape/src/lineshape_new_ref/photonics2/photoluminescence.py \
  && wc -l qqs/lineshape/src/lineshape_new_ref/photonics2/photoluminescence.py \
  && test -f qqs/lineshape/src/lineshape_new_ref/photonics2/plott.py \
  && wc -l qqs/lineshape/src/lineshape_new_ref/photonics2/plott.py
```

期望：`photoluminescence.py` 约 359 行；`plott.py` 约 235 行。

---

### Task 2：备份三个被替换文件

**涉及文件：**
- 创建：`qqs/lineshape/src/photonics2/photoluminescence.py.bak`
- 创建：`qqs/lineshape/src/photonics2/plott.py.bak`
- 创建：`qqs/lineshape/src/pl.py.bak`

**接口：**
- 进入：Task 1 已确认源存在
- 产出：三个回滚点

- [ ] **Step 1：复制三个原文件到备份**

```bash
cp qqs/lineshape/src/photonics2/photoluminescence.py \
   qqs/lineshape/src/photonics2/photoluminescence.py.bak
cp qqs/lineshape/src/photonics2/plott.py \
   qqs/lineshape/src/photonics2/plott.py.bak
cp qqs/lineshape/src/pl.py \
   qqs/lineshape/src/pl.py.bak
ls -l qqs/lineshape/src/photonics2/*.bak qqs/lineshape/src/pl.py.bak
```

期望：三个备份文件存在且大小非零。

---

### Task 3：替换 `photonics2/photoluminescence.py`

**涉及文件：**
- 修改：`qqs/lineshape/src/photonics2/photoluminescence.py`（覆盖为 `lineshape_new_ref` 的向量化版本）

- [ ] **Step 1：用 `lineshape_new_ref` 版本覆盖**

```bash
cp qqs/lineshape/src/lineshape_new_ref/photonics2/photoluminescence.py \
   qqs/lineshape/src/photonics2/photoluminescence.py
```

- [ ] **Step 2：校验 29 条决策（grep）**

```bash
F=qqs/lineshape/src/photonics2/photoluminescence.py
# 行数 / 导入
wc -l $F
grep -c "^import multiprocessing" $F
grep -c "^import yaml" $F
grep -c "^import yaml" $F
# #1 / #2 sys.exit → raise in phonopy_read_*
grep -c "raise$" $F
grep -c "sys.exit" $F
# #3 fold 注释删除（允许 0）
grep -c "print(nml)" $F
# #4 plotdata 删除（允许 0）
grep -c "def plotdata" $F
# #5 get_phonon 手写流（无 yaml.load）
grep -c "yaml.load" $F
grep -c "lattice_inform" $F
grep -c "phonon_band_in_one_kpoint" $F
# #6 / #7 get_phonon 返回 modes / frequencies + 归一化
grep -c "norms\[norms == 0\] = 1" $F
# #8 read_grd_ex_pos 多 path 参数
grep -c "def read_grd_ex_pos(self, str_g, str_e, shiftVector, n_defect, defect_range, path)" $F
# #9 D_R 局部变量
grep -c "^        D_R = \[\]" $F
# #10 n_defect 质量加权
grep -c "sumd += self.D_R\[i\] \* weight" $F
# #11 / #12 get_S_omega / get_C_omega / occupy 实例方法删除
grep -c "def get_S_omega" $F
grep -c "def get_C_omega" $F
grep -c "def occupy" $F
# #13 el_ph 调模块级向量化
grep -c "S_omega_vectorized(" $F
grep -c "C_omega_vectorized(" $F
# #14 el_ph_ 删除（多进程版）
grep -c "def el_ph_" $F
# #15 write_S 删除
grep -c "def write_S" $F
# #16 C_omega.data 落盘保留
grep -c "C_omega.data" $F
# #17 S_omega_vectorized 低频 0.3×sigma
grep -c "np.where(freqs < 0.005, sigma, 0.3 \* sigma)" $F
# #18 / #19 PL / PLA 错误处理
grep -c 'raise ValueError(f"Unknown process' $F
# #20 PL 中 Gt 向量化
grep -c "Gt = np.exp(St + Ct + Ct\[::-1\]" $F
# #21 PL 中 A 移位
grep -c "idx = (shift - np.arange(n)) % n" $F
# #22 PLA 列表累加 → numpy
grep -c "self.I = self.A \* (t \* r)\*\*3" $F
# #23 __init__ "large" 分支删除
grep -c '"large" in method' $F
# #24 defect_range 默认 -1.0
grep -c 'defect_range", -1.0' $F
# #25 skipmodes 阈值 ≤ 0.00
grep -c "self.frequencies\[i\] <= 0.00" $F
# #26 HuangRhyes "large" 分支删除
grep -c '"large" in self.method' $F
# #27 HuangRhyes IPR/q 向量化
grep -c "participation = np.sum(self.Modes\*\*2, axis=2)" $F
# #28 print_table with-open 风格
grep -c "with open(\"main_modes.data\"" $F
# #29 print_table 触发条件无括号
grep -c "if self.S\[i\] > k and self.frequencies\[i\] > 0.04 or i in visualmode:" $F
```

期望清单：

| 检查 | 期望 |
|------|------|
| 行数 | ≈ 359 |
| `^import multiprocessing` | 0 |
| `^import yaml` | 0 |
| `raise$` | ≥ 2（#1, #2） |
| `sys.exit` | 0 |
| `print(nml)` | 0（#3 注释清理） |
| `def plotdata` | 0（#4 删除） |
| `yaml.load` | 0（#5） |
| `lattice_inform` | 0（#5） |
| `phonon_band_in_one_kpoint` | 0（#5） |
| `norms[norms == 0] = 1` | ≥ 1（#7 归一化） |
| `def read_grd_ex_pos(self, str_g, str_e, shiftVector, n_defect, defect_range, path)` | ≥ 1（#8） |
| `^        D_R = \[\]` | ≥ 1（#9 局部变量） |
| `sumd += self.D_R[i] \* weight` | ≥ 1（#10 质量加权） |
| `def get_S_omega` | 0（#11 删除） |
| `def get_C_omega` | 0（#12 删除） |
| `def occupy` | 0（#12 删除） |
| `S_omega_vectorized(` | ≥ 1（#13） |
| `C_omega_vectorized(` | ≥ 1（#13） |
| `def el_ph_` | 0（#14 删除） |
| `def write_S` | 0（#15 删除） |
| `C_omega.data` | ≥ 1（#16 保留） |
| `np.where(freqs < 0.005, sigma, 0.3 * sigma)` | ≥ 1（#17） |
| `raise ValueError(f"Unknown process` | ≥ 2（#18, #19） |
| `Gt = np.exp(St + Ct + Ct[::-1]` | ≥ 1（#20） |
| `idx = (shift - np.arange(n)) % n` | ≥ 1（#21） |
| `self.I = self.A * (t * r)**3` | ≥ 1（#22） |
| `"large" in method` | 0（#23 删除） |
| `defect_range", -1.0` | ≥ 1（#24） |
| `self.frequencies[i] <= 0.00` | ≥ 1（#25） |
| `"large" in self.method` | 0（#26 删除） |
| `participation = np.sum(self.Modes**2, axis=2)` | ≥ 1（#27） |
| `with open("main_modes.data"` | ≥ 1（#28） |
| `if self.S[i] > k and self.frequencies[i] > 0.04 or i in visualmode:` | ≥ 1（#29） |

任一项不符即停止，回滚 `.bak`。

---

### Task 4：合并 `photonics2/plott.py`

**涉及文件：**
- 修改：`qqs/lineshape/src/photonics2/plott.py`（以 `lineshape_new_ref` 为基线，按 spec 决策表删除 4 行）

**接口：**
- 进入：Task 3 已完成（`photoluminescence.py` 已替换为向量化版）
- 产出：合并后的 `plott.py`；行为满足 spec 决策表 10 条

- [ ] **Step 1：用 `lineshape_new_ref` 版本覆盖（基线）**

```bash
cp qqs/lineshape/src/lineshape_new_ref/photonics2/plott.py \
   qqs/lineshape/src/photonics2/plott.py
```

- [ ] **Step 2：删除决策 2/5/6 的 3 处 `np.savetxt`（落盘）**

每处都是单行 `np.savetxt(...)`，整行删除。用 sed：

```bash
sed -i "/np\.savetxt('C_omega,S_omega.data'/d" qqs/lineshape/src/photonics2/plott.py
sed -i "/np\.savetxt('lineshape_eV.csv'/d" qqs/lineshape/src/photonics2/plott.py
sed -i "/np\.savetxt('lineshape_nm.csv'/d" qqs/lineshape/src/photonics2/plott.py
grep -c "np.savetxt" qqs/lineshape/src/photonics2/plott.py
```

期望：`np.savetxt` 计数为 `0`。

- [ ] **Step 3：删除决策 7 的 `split` 行**

```bash
sed -i '/split = para\.get("split", 0)/d' qqs/lineshape/src/photonics2/plott.py
grep -c 'split = para\.get' qqs/lineshape/src/photonics2/plott.py
```

期望：计数为 `0`。

- [ ] **Step 3.5：删除决策 3 的 `use_line_collection=True`（matplotlib ≥ 3.10 不兼容）**

new_ref 的 `plott.py` 在 `Shw` 分支末尾用 `ax2.stem(..., use_line_collection=True)`。`matplotlib.Axes.stem` 自 3.10 起移除该参数。`run_compare.py` 不受影响（不走 `plott.py`），但 `python src/pl.py` 冒烟测试会抛 `TypeError`。手动 sed 删除：

```bash
sed -i 's/, use_line_collection=True//' qqs/lineshape/src/photonics2/plott.py
grep -c "use_line_collection" qqs/lineshape/src/photonics2/plott.py
```

期望：`use_line_collection` 计数为 `0`。

- [ ] **Step 4：校验 10 条决策**

```bash
# 决策 1 — Shw 含 C_omega 红线、无 fill_between
grep -c "C_omega, color=\"tomato\"" qqs/lineshape/src/photonics2/plott.py
grep -c "fill_between" qqs/lineshape/src/photonics2/plott.py
# 决策 4 — Shw.png 写入：1 处 Shw 分支 + 1 处 Shw+Acm 分支（后者合法，非重复 bug）
grep -c 'plt.savefig("Shw.png"' qqs/lineshape/src/photonics2/plott.py
# 决策 8 — Iabs_new
grep -c "Iabs_new" qqs/lineshape/src/photonics2/plott.py
# 决策 9 — Shw+Acm 内无 xcm
awk '/Shw\+Acm/,/^$/' qqs/lineshape/src/photonics2/plott.py | grep -c "xcm"
```

期望（按顺序）：
- `C_omega, color="tomato"` ≥ 1
- `fill_between` = 0
- `plt.savefig("Shw.png"` = 2（不是 1：1 处 Shw 分支末尾 + 1 处 Shw+Acm 组合图合法写入）
- `Iabs_new` ≥ 1
- `xcm` in Shw+Acm block = 0

注：决策 3 `use_line_collection=True` 的 grep 已移到 Step 3.5 删除，校验列表中不再列。

---

### Task 5：合并并替换 `src/pl.py`

**涉及文件：**
- 修改：`qqs/lineshape/src/pl.py`（按 spec 决策表合并，覆盖现有内容）

**接口：**
- 进入：Task 3、4 已完成（`photoluminescence.py`、`plott.py` 已替换为向量化版）
- 产出：合并后的 `src/pl.py`，行为满足 spec 决策表 9 条决策

- [ ] **Step 1：写合并后的 `src/pl.py` 内容（直接覆盖）**

按下面整段内容写入 `qqs/lineshape/src/pl.py`：

```python
#!/usr/bin/env python3
"""Driver for the photoluminescence pipeline.

Reads band.yaml + ground/excited POSCARs and produces Huang-Rhys data + 6 plots.

Usage:
    python src/pl.py [path/to/INCAR]
"""
import matplotlib
matplotlib.use('Agg')                       # decision #2 — headless runs
import os
import sys
import numpy as np
import matplotlib.pyplot as plt

from photonics2.photoluminescence import Photoluminescence
from photonics2.plott import plot_S_I


def get_pl(path_phonopy, ground_struct_path, excited_struct_path, **parameter):
    # decision #3 — read resolution from kwargs (INCAR) so it can be overridden
    res = parameter.get("resolution", 1000)
    proc = parameter.get("process", "emission")
    plot_width = parameter.get("plot_width", 1.0)
    if "emi" in proc:
        # decision #4
        Amin = -0.8 * plot_width
        Amax = 0.2 * plot_width
    elif "abs" in proc:
        # decision #5
        Amin = -0.2 * plot_width
        Amax = 0.8 * plot_width

    print("################# reading band.yaml Start! ###################")
    mass = parameter.get("mass", [])
    if len(mass) > 0:
        p = Photoluminescence(path_phonopy, "phonopy", m=mass,
                              POSCAR_GRD=ground_struct_path, POSCAR_EX=excited_struct_path,
                              n_defect=1, resolution=res)
    else:
        p = Photoluminescence(path_phonopy, "phonopy",
                              POSCAR_GRD=ground_struct_path, POSCAR_EX=excited_struct_path,
                              n_defect=1, resolution=res)
    print("################# band.yaml OK! ###################")

    print("################# HuangRhyes Start! ###################")
    HR = p.HuangRhyes()
    print("HuangRhyes=", HR)

    f = open("partial.HuangRhyes.data", "w")
    f.write("mode No.\t partial.HuangRhyes \n")
    i = 0
    for s in p.S:
        i += 1
        f.write(str(i) + "\t" + str(s) + "\n")
    f.close()
    p.print_table(0.1, [])
    print("################# HuangRhyes OK! ###################")

    print("################# el_ph Start! ###################")
    gw = parameter.get("gw", 1e-3)
    p.el_ph(delta_width=gw, temperature=parameter.get("T", 0), jtmodes=[])
    print("################# el_ph OK! ###################")

    print("################# HuangRhyes lineshape Start! ###################")
    zpl = parameter.get("zpl", 2.5)
    A = p.PL(gamma=parameter.get("gamma", 10e-3) * res, SHR=0, EZPL=zpl, process=proc)
    p.PLA()
    print("################# HuangRhyes lineshape OK! ###################")
    # decision #6 — 6 plots
    plot_S_I(p, parameter.get("title", "any"),
             [p.EZPL + Amin * 0.5, p.EZPL + Amax * 0.5],
             "Shw Sk PLeV  PLnm A Phon")
    return p


def read_parameters_from_incar(incar_path):
    params = {}
    with open(incar_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, value = line.split("=", 1)
                key = key.strip()
                value = value.strip()
                try:
                    if "." in value:
                        value = float(value)
                    else:
                        value = int(value)
                except ValueError:
                    pass
                params[key] = value
    return params


# decision #7 — INCAR resolved relative to this script's directory
incar_path = os.path.join(os.path.dirname(__file__), "INCAR")
if len(sys.argv) > 1:
    incar_path = sys.argv[1]

incar_params = read_parameters_from_incar(incar_path)
get_pl(
    # decision #8 — post-restructure cases/ layout
    incar_params.get("BAND_YAML", "../cases/zlq/band.yaml"),
    incar_params.get("GRD_POSCAR", "../cases/zlq/GS"),
    incar_params.get("EXC_POSCAR", "../cases/zlq/ES"),
    gamma=incar_params.get("gamma", 10e-3),
    zpl=incar_params.get("zpl", 1.0),
    gw=incar_params.get("gw", 1e-3),
    resolution=incar_params.get("resolution", 1000),
    T=incar_params.get("T", 0),
    plot_width=incar_params.get("plot_width", 2.5),
    title=incar_params.get("title", "any"),
    process=incar_params.get("process", "emission")
)
```

写入命令（保持单次 Write 即可，使用 Write 工具直接覆写 `qqs/lineshape/src/pl.py`）。

- [ ] **Step 2：校验 9 条决策在文件中都体现**

```bash
grep -n "^#!/usr/bin/env python3" qqs/lineshape/src/pl.py
grep -n "matplotlib.use('Agg')" qqs/lineshape/src/pl.py
grep -n 'parameter.get("resolution"' qqs/lineshape/src/pl.py
grep -n -- "-0.8 \* plot_width" qqs/lineshape/src/pl.py
grep -n -- "-0.2 \* plot_width" qqs/lineshape/src/pl.py
grep -n 'Shw Sk PLeV  PLnm A Phon' qqs/lineshape/src/pl.py
grep -n 'os.path.dirname(__file__)' qqs/lineshape/src/pl.py
grep -n '../cases/zlq' qqs/lineshape/src/pl.py
```

期望：每个 grep 都至少命中 1 行。

---

### Task 6：跑 baseline 回归（`run_compare.py`）

**涉及文件：**
- 只读：`qqs/lineshape/run_compare.py`、`qqs/lineshape/BASELINE.md`、`qqs/lineshape/cases/`

- [ ] **Step 1：跑 `run_compare.py`**

```bash
cd qqs/lineshape && python run_compare.py
```

期望：退出码 0；8 个案例全部 `✅ ✅ ✓`。

- [ ] **Step 2：交叉对比 `BASELINE.md`**

```bash
cd qqs/lineshape && python - <<'PY'
import re, pathlib
baseline = pathlib.Path("BASELINE.md").read_text()
rows = {}
for line in baseline.splitlines():
    m = re.match(r"\|\s*([\w-]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|\s*(-?\d+)\s*\|\s*(\d+)\s*\|", line)
    if m:
        rows[m.group(1)] = tuple(m.group(i) for i in (2,3,4,5,6))
print("baseline cases:", sorted(rows))
PY
```

对比每案例数值与 `run_compare.py` 输出。期望：A 列、B 列都与 baseline 对应行一致。

若失败：停止。回滚三个文件：
```bash
cp qqs/lineshape/src/photonics2/photoluminescence.py.bak qqs/lineshape/src/photonics2/photoluminescence.py
cp qqs/lineshape/src/photonics2/plott.py.bak qqs/lineshape/src/photonics2/plott.py
cp qqs/lineshape/src/pl.py.bak qqs/lineshape/src/pl.py
```
然后回到 spec 复核。

---

### Task 7：冒烟测试 `python src/pl.py`

**涉及文件：**
- 临时修改：`qqs/lineshape/src/INCAR`（`resolution` 改为 `1000`）
- 生成：`Shw.png`、`Sk.png`、`PLev.png`、`PLnm.png`、`AeV.png`、`Phon.png`、`partial.HuangRhyes.data` 等产物到 cwd

- [ ] **Step 1：临时把 `src/INCAR` 的 `resolution` 改成 `1000`**

```bash
sed -i 's/^resolution = .*/resolution = 1000/' qqs/lineshape/src/INCAR
grep "^resolution" qqs/lineshape/src/INCAR
```

期望：`resolution = 1000`。

- [ ] **Step 2：从 `qqs/lineshape/src/` 跑 `python pl.py`**

```bash
cd qqs/lineshape/src && python pl.py
```

期望：脚本打印各阶段 OK；进程退出码 0。

- [ ] **Step 3：检查 6 张 PNG 都生成**

```bash
cd qqs/lineshape/src && ls -1 Shw.png Sk.png PLev.png PLnm.png AeV.png Phon.png 2>&1
```

期望：6 个文件全部存在，无 `No such file or directory`。

- [ ] **Step 4：确认无意外产物（spec 验收 #6）**

```bash
cd qqs/lineshape/src && ls C_omega,S_omega.data lineshape_eV.csv lineshape_nm.csv 2>&1
```

期望：`No such file or directory`（3 个被禁的产物都不存在）。

- [ ] **Step 5：清理生成产物，恢复 `src/INCAR`**

```bash
cd qqs/lineshape/src && rm -f Shw.png Sk.png PLev.png PLnm.png AeV.png Phon.png partial.HuangRhyes.data main_modes.data modes.data
git checkout -- qqs/lineshape/src/INCAR
grep "^resolution" qqs/lineshape/src/INCAR
```

期望：6 个 PNG 与 `partial.HuangRhyes.data` 被删除；`src/INCAR` 恢复为 `resolution = 4000`。

---

### Task 8：删除三个备份

**涉及文件：**
- 删除：`qqs/lineshape/src/photonics2/photoluminescence.py.bak`
- 删除：`qqs/lineshape/src/photonics2/plott.py.bak`
- 删除：`qqs/lineshape/src/pl.py.bak`

- [ ] **Step 1：删除备份**

```bash
rm qqs/lineshape/src/photonics2/photoluminescence.py.bak \
   qqs/lineshape/src/photonics2/plott.py.bak \
   qqs/lineshape/src/pl.py.bak
ls qqs/lineshape/src/photonics2/*.bak qqs/lineshape/src/pl.py.bak 2>&1
```

期望：`No such file or directory`。

---

### Task 9：更新 AGENTS.md

**涉及文件：**
- 修改：`AGENTS.md` —— 描述 `qqs/lineshape/` 的条目

- [ ] **Step 1：定位现有描述**

```bash
grep -n "qqs/lineshape" AGENTS.md
```

- [ ] **Step 2：替换该条目**

把现有那段：

> `qqs/lineshape/` (legacy, unmaintained) — earlier version of the photoluminescence and Jahn-Teller code from the original repository author. Contains 17 Python source files (`src/photonics2/`), two standalone scripts (`src/pl.py`, `src/phonon_struct_op.py`), 8 case study directories with VASP/phonopy input data (`cases/`), and regeneratable computation outputs (`output/`). Includes an alternative earlier snapshot with numpy-vectorized implementations (`src/lineshape_new_ref/`). Distinct code conventions: no type hints, no `from __future__ import annotations`, bare `except:` with `print()+sys.exit()`, absolute intra-package imports, `multiprocessing`/`concurrent.futures` imports.

替换为：

> `qqs/lineshape/` (legacy, unmaintained) — earlier version of the photoluminescence and Jahn-Teller code from the original repository author. Contains 16 Python source files (`src/photonics2/`: 1 vectorized PL driver + 10 unmaintained JT/embedding/Schrödinger helpers + 5 shared files; `src/pl.py`, `src/phonon_struct_op.py` standalone), 8 case study directories with VASP/phonopy input data (`cases/`), and regeneratable computation outputs (`output/`). 自 2026-07-11 起，`src/photonics2/photoluminescence.py` 和 `src/photonics2/plott.py` 与 `src/lineshape_new_ref/` 内的同名文件算法等价；`src/pl.py` 按 spec 决策表合并（详见 `docs/superpowers/specs/2026-07-11-qqs-lineshape-merge-design.md`）。保留 `lineshape_new_ref/` 是为了 `run_compare.py` 的 A/B 回归。代码规约：no type hints, no `from __future__ import annotations`, bare `except:` with `print()+sys.exit()`, absolute intra-package imports, `multiprocessing`/`concurrent.futures` imports.

- [ ] **Step 3：校验**

```bash
grep -A1 "2026-07-11" AGENTS.md
```

期望：含 `算法等价` 与 spec 链接。

---

### Task 10：更新 CLAUDE.md

**涉及文件：**
- 修改：`CLAUDE.md` —— "Legacy `qqs/lineshape/` directory" 段

- [ ] **Step 1：在该段目录树后追加一行**

在以 `src/lineshape_new_ref/` 开头那行后加：

> 自 2026-07-11 起，`src/photonics2/{photoluminescence.py, plott.py}` 与 `src/lineshape_new_ref/` 内的同名文件算法等价；`src/pl.py` 按 spec 决策表合并。`lineshape_new_ref/` 保留是为了 `run_compare.py` A/B 回归；详见 `docs/superpowers/specs/2026-07-11-qqs-lineshape-merge-design.md`。

- [ ] **Step 2：校验**

```bash
grep -A1 "2026-07-11" CLAUDE.md
```

期望：在 legacy 段出现 `算法等价`。

---

### Task 11：终验

- [ ] **Step 1：确认 `photoluminescence.py` 与来源字节相同**

```bash
diff -q qqs/lineshape/src/photonics2/photoluminescence.py \
        qqs/lineshape/src/lineshape_new_ref/photonics2/photoluminescence.py
```

期望：无输出。

- [ ] **Step 2：确认 `plott.py` 已按决策表合并（4 处删除 + 1 处 use_line_collection 删除、其余取 new_ref）**

```bash
grep -c "np.savetxt" qqs/lineshape/src/photonics2/plott.py
grep -c 'split = para\.get' qqs/lineshape/src/photonics2/plott.py
grep -c "use_line_collection" qqs/lineshape/src/photonics2/plott.py
grep -c "C_omega, color=\"tomato\"" qqs/lineshape/src/photonics2/plott.py
grep -c 'plt.savefig("Shw.png"' qqs/lineshape/src/photonics2/plott.py
```

期望（按顺序）：`np.savetxt` = 0，`split` = 0，`use_line_collection` = 0，`C_omega, color="tomato"` ≥ 1，`plt.savefig("Shw.png"` = 2。

- [ ] **Step 3：确认 `src/pl.py` 不再含旧版 shebang**

```bash
head -1 qqs/lineshape/src/pl.py
```

期望：`#!/usr/bin/env python3`。

- [ ] **Step 4：确认无残留备份**

```bash
ls qqs/lineshape/src/photonics2/*.bak qqs/lineshape/src/pl.py.bak 2>&1
```

期望：`No such file or directory`。

- [ ] **Step 5：最终再跑一次 `run_compare.py`**

```bash
cd qqs/lineshape && python run_compare.py
```

期望：8 个案例 `✅ ✅ ✓`。

- [ ] **Step 6：确认 `src/INCAR` 仍是 `resolution = 4000`**

```bash
grep "^resolution" qqs/lineshape/src/INCAR
```

期望：`resolution = 4000`。

- [ ] **Step 7：记录合并后 A 列 skipmodes（用于更新 `BASELINE.md`）**

```bash
cd qqs/lineshape && python run_compare.py 2>&1 | tee /tmp/post_merge_compare.log
```

期望：所有 8 case A/B HR 一致；A 列 skipmodes 现在应等于 B 列（合并后两边都用 `freq ≤ 0.0` 阈值）。记下新数字，**不要**直接改 `BASELINE.md`，等 Task 12 一并改。

- [ ] **Step 8：更新 `AGENTS.md` `qqs/lineshape/` 段**

line 57 现有描述 "Contains 17 Python source files (`src/photonics2/`)" 在合并后过期（应是 16）。把整段：

> `qqs/lineshape/` (legacy, unmaintained) — earlier version of the photoluminescence and Jahn-Teller code from the original repository author. Contains 17 Python source files (`src/photonics2/`), two standalone scripts (`src/pl.py`, `src/phonon_struct_op.py`), 8 case study directories with VASP/phonopy input data (`cases/`), and regeneratable computation outputs (`output/`). Includes an alternative earlier snapshot with numpy-vectorized implementations (`src/lineshape_new_ref/`). Distinct code conventions: no type hints, no `from __future__ import annotations`, bare `except:` with `print()+sys.exit()`, absolute intra-package imports, `multiprocessing`/`concurrent.futures` imports.

替换为：

> `qqs/lineshape/` (legacy, unmaintained) — earlier version of the photoluminescence and Jahn-Teller code from the original repository author. 自 2026-07-11 合并后：`src/photonics2/photoluminescence.py` 和 `src/photonics2/plott.py` 与 `src/lineshape_new_ref/` 内的同名文件算法等价（cp 覆盖 + spec 决策表合并）；`src/pl.py` 按 spec 决策表合并。`src/photonics2/` 共 16 个 Python 源文件（含 10 个无人维护 JT/embedding/Schrödinger 旧模块）；`src/pl.py`、`src/phonon_struct_op.py` 顶层脚本；`src/lineshape_new_ref/` 保留供 `run_compare.py` A/B 回归；8 个 VASP/phonopy 案例在 `cases/`；regeneratable 计算输出在 `output/`。详见 `docs/superpowers/specs/2026-07-11-qqs-lineshape-merge-design.md`。代码规约：no type hints, no `from __future__ import annotations`, bare `except:` with `print()+sys.exit()`, absolute intra-package imports, `multiprocessing`/`concurrent.futures` imports.

- [ ] **Step 9：更新 `BASELINE.md` Notes**

合并前 `BASELINE.md` line 32 写："A uses threshold `freq ≤ 0.005 eV` while B uses `freq ≤ 0.0 eV`。合并后 A 和 B 都用 `freq ≤ 0.0 eV`。把：

> - **skipmodes differs** because A uses threshold `freq ≤ 0.005 eV` while B uses `freq ≤ 0.0 eV`. This does not affect the final HR since skipped modes contribute negligibly.

替换为：

> - **skipmodes identical** between A and B after the 2026-07-11 merge (both use threshold `freq ≤ 0.0 eV`).
> - **Recorded commit**: see git log; **pre-merge commit** was `d43f491`.

如果 Step 7 跑出的 A 列 skipmodes 数字与 B 列不同（合并后理论上应该相同），记录差异到 spec/plan 附录，但**不要**擅自改 `BASELINE.md` 的数字列；只更新 Notes 段。

---

## 回滚

若任何 Task 6 或 Task 7 失败且未进入 Task 8：

```bash
cp qqs/lineshape/src/photonics2/photoluminescence.py.bak qqs/lineshape/src/photonics2/photoluminescence.py
cp qqs/lineshape/src/photonics2/plott.py.bak qqs/lineshape/src/photonics2/plott.py
cp qqs/lineshape/src/pl.py.bak qqs/lineshape/src/pl.py
git checkout -- qqs/lineshape/src/INCAR
```

若 Task 9/10 文档已改但代码回滚：

```bash
git checkout -- AGENTS.md CLAUDE.md
```