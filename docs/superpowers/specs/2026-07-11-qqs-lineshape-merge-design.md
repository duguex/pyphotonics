# qqs/lineshape/ 光致发光算法统一 设计文档

**日期**：2026-07-11
**状态**：已批准
**相关文档**：`2026-07-10-qqs-lineshape-reorganization-design.md`（目录重整，已完成）

## 目标

把 `qqs/lineshape/src/` 与 `qqs/lineshape/src/lineshape_new_ref/` 下三个存在差异的文件统一为单一版本，同时**保留两份代码库在磁盘上**，以便 `run_compare.py` 继续做 A/B 回归对比。

三个存在差异的文件：

1. `photonics2/photoluminescence.py` —— PL 算法核心
2. `photonics2/plott.py` —— 纯格式调整
3. `pl.py`（顶层 driver 脚本）—— 存在真实的行为差异（详见"Differences"）

## 背景

`qqs/lineshape/` 携带了同一光致发光计算的两种实现：

| 位置 | 作用 |
|------|------|
| `qqs/lineshape/src/photonics2/` | 原始实现：逐模 Python 循环、`multiprocessing.Pool`、PyYAML 解析。仍由 `pl.py` 调用 |
| `qqs/lineshape/src/lineshape_new_ref/` | 早期的 numpy 向量化快照：S(ω) 沿模轴广播、手写 YAML 流解析、单线程 |

之前的重整 spec（`2026-07-10`）只是把文件移到新的 `src/cases/output/` 目录结构，没有触及算法层面的分叉。

`run_compare.py` 会用两种版本跑全部 8 个案例并比较 HR / ΔQ。`BASELINE.md` 记录了用来校验两版的数值基线。

## 两版文件清单对比

对比 `src/photonics2/`（16 个 `.py`）与 `src/lineshape_new_ref/`（6 个 `.py`，含顶层 1 个 `pl.py` + 子包 5 个）。

### 同名文件：6 对

| 文件 | 状态 | 备注 |
|------|------|------|
| `pl.py`（顶层） | **不同** | driver 脚本 —— 详见 "Differences" |
| `photonics2/photoluminescence.py` | **不同** | PL 核心类 —— 详见 "Differences" |
| `photonics2/plott.py` | **不同** | 纯格式调整，无功能影响 |
| `photonics2/configuration_coordinate.py` | 字节级相同 | 84 行 |
| `photonics2/xyz.py` | 字节级相同 | 20 行 |
| `photonics2/__init__.py` | 字节级相同（均空） | |

同名对中**3 对不同、3 对相同**。

### 仅在 `photonics2/` 中存在（旧版独占，10 个文件）

`constants.py`、`ebedding.py`、`embedding.py`、`hermite.py`、`ht.py`、`jt.py`、`jtsoc.py`、`nu.py`、`schrodinger.py`、`un.py`。

**本次不动。** 它们是无人维护的 Jahn-Teller / embedding / Schrödigner 辅助模块；当前 PL 流水线无人 import（已用 `grep -rE "from photonics2\.(jt|jtsoc|ht|schrodinger|embedding|ebedding|hermite|nu|un)|^import photonics2\.(jt|jtsoc|ht|schrodinger|embedding|ebedding|hermite|nu|un)"` 验证）。未来可单独清理 spec 处理，本次超出范围。

### 仅在 `lineshape_new_ref/` 中存在

无。

## 详细差异

### `photonics2/photoluminescence.py`

差异细分 29 条。决策表（依据 `2026-07-11` 复核）：

| # | 维度 | `photonics2/`（当前） | `lineshape_new_ref/`（目标） | **决策** |
|---|------|------------------------|------------------------------|----------|
| 1 | `phonopy_read_modes` 错误处理 | `sys.exit()` | `raise` | 取 new_ref |
| 2 | `phonopy_read_frequencies` 错误处理 | `sys.exit()` | `raise` | 取 new_ref |
| 3 | `fold` 注释 | 注释掉的 `print(nml)` | 删除 | 取 new_ref（清理） |
| 4 | `plotdata` 实例方法 | 存在（debug 用 `plt.show`） | **删除** | 取 new_ref（删除；外部无人调用） |
| 5 | `get_phonon` 实现 | `yaml.load(Loader=yaml.CLoader)` + `lattice_inform` + `phonon_band_in_one_kpoint` | 手写流解析 + 直接 numpy | 取 new_ref |
| 6 | `get_phonon` 返回值 | `(m, modes, fre)` 来自包装类 | `(m, modes, frequencies)` 直接 numpy | 取 new_ref（API 不变） |
| 7 | `get_phonon` 内部归一化 | 包装类内 `np.linalg.norm(reshape(-1,1))` | `get_phonon` 末尾 `modes = modes / sqrt_m; modes = modes / norms` | 取 new_ref |
| 8 | `read_grd_ex_pos` 签名 | 无 `path` 参数 | **新增 `path` 参数**（`__init__` 内传） | 取 new_ref（API 变化） |
| 9 | `D_R` 列表构建 | `self.D_R = []` + `.append()` | 局部 `D_R = []` 后 `self.D_R = np.array(D_R)` | 取 new_ref（行为等价） |
| 10 | `n_defect > 0` 时 D_R 偏移加权 | 等权：`sumd += D_R[i]; sumd /= count` | **质量加权**：`sumd += D_R[i] * m[i]; sumd /= Σm` | **取 new_ref（质量加权）**。`run_compare.py` 在 8 个 case 上传 `n_defect=0`，本分支不执行，故 baseline 不受影响；若以后跑 `n_defect>0` 需补 baseline |
| 11 | `get_S_omega` 实例方法 | 存在，调用 `self.S[k] * gaussian(...)` | **删除**，被模块级 `S_omega_vectorized` 替代 | 取 new_ref（删除） |
| 12 | `get_C_omega` + `occupy` 实例方法 | 存在 | **删除**，被模块级 `C_omega_vectorized` 替代 | 取 new_ref（删除） |
| 13 | `el_ph` 实现 | 调 `self.get_S_omega` / `self.get_C_omega`（实例方法 + 列表推导） | 调模块级 `S_omega_vectorized` / `C_omega_vectorized`（numpy） | 取 new_ref |
| 14 | `el_ph_` 方法（多进程版） | 存在（`multiprocessing.Pool`，8 线程） | **删除** | 取 new_ref（删除；`run_compare.py` 走 `el_ph` 不走 `el_ph_`） |
| 15 | `write_S` 方法 | 存在（写 S_omega 到文件） | **删除** | 取 new_ref（删除；外部无人调用） |
| 16 | `el_ph` 中 `C_omega.data` 落盘 | `with open("C_omega.data", "w")` 写每一项 | 同 | 保留（两边一致） |
| 17 | `S_omega_vectorized` 低频 sigma | 旧版 `get_S_omega` 内部 `if frequencies[k]<0.035: sigma else: sigma` —— **两个分支完全相同**（bug 残留） | `sigma_k = np.where(freqs < 0.005, sigma, 0.3*sigma)` —— 低频用 0.3×sigma | 取 new_ref（修 bug） |
| 18 | `PL()` 错误处理 | `sys.exit()` | `raise ValueError(f"Unknown process: {self.process}")` | 取 new_ref |
| 19 | `PLA()` 错误处理 | `sys.exit()` | `raise ValueError(f"Unknown process: {self.process}")` | 取 new_ref |
| 20 | `PL()` `Gt` 构造 | Python for 循环 `Gt += [np.exp(...)]` | numpy 向量化 `Gt = np.exp(St + Ct + Ct[::-1] - St[0] - 2*Ct[0] - SHR - self.gamma*np.abs(t))` | 取 new_ref |
| 21 | `PL()` A 循环移位 | `tA.copy()` + 显式 for 循环赋索引 | `shift = ...; idx = (shift - np.arange(n)) % n; self.A = A[idx]` | 取 new_ref |
| 22 | `PLA()` 列表累加 | `I += [self.A[i]*((t)*r)**3]` | `self.I = self.A * (t * r)**3` | 取 new_ref |
| 23 | `__init__` `"large"` 方法分支 | 存在（`if "large" in method: ... read_force ...`） | **删除**，只剩 `else: self.read_grd_ex_pos(...)` | 取 new_ref（删除；`run_compare.py` 走 phonopy 路径不触发） |
| 24 | `__init__` 调用 `read_grd_ex_pos` 时 `defect_range` 默认值 | `4.0` | `-1.0` | 取 new_ref |
| 25 | `__init__` skipmodes 触发条件 | `if self.frequencies[i] <= 0.005:` | `if self.frequencies[i] <= 0.00:` | 取 new_ref。阈值从 ≤0.005 eV 降到 ≤0.0 eV。8 个 case baseline 一致即证明实际影响为零 |
| 26 | `HuangRhyes` `"large"` 方法分支 | 存在（`if "large" in self.method: ...` 用 force_grd/force_ex） | **删除** | 取 new_ref（删除；`run_compare.py` 不触发） |
| 27 | `HuangRhyes` IPR / q 计算 | 嵌套 for 循环按模式 i | numpy 向量化 `participation = np.sum(Modes**2, axis=2); IPR = 1/Σ; q = sqrt(m)·D_R·Modes` | 取 new_ref |
| 28 | `print_table` 风格 | 多个 open + 手动 close | `with open(...) as f, open(...) as f2:` + f-string | 取 new_ref（行为等价） |
| 29 | `print_table` 触发条件 | `if(self.S[i]>k and self.frequencies[i]>0.04 or i in visualmode):` | `if self.S[i] > k and self.frequencies[i] > 0.04 or i in visualmode:` | 取 new_ref（括号差异不影响 Python `and` 优先级高于 `or`） |

合并动作：**整文件以 `lineshape_new_ref/photonics2/photoluminescence.py` 为基线 cp 覆盖**。29 条差异全部取 new_ref；其中 #4/#11/#12/#14/#15/#23/#26 是删除若干方法/函数（`plotdata`/`get_S_omega`/`get_C_omega`/`occupy`/`el_ph_`/`write_S`/`__init__ "large"` 分支/`HuangRhyes "large"` 分支），不在 `Photoluminescence` 公开 API 上（公开 API 为 `__init__`/`HuangRhyes`/`el_ph`/`PL`/`PLA`/`print_table`，不变）。

新增 API 风险（spec 风险表需补一行）：`read_grd_ex_pos` 多一个 `path` 位置参数（#8），外部脚本若按旧签名调用会 TypeError。已知 `run_compare.py` 通过 `__init__` 间接调用，不受影响；其他外部调用方需在 spec 文档化。

不破坏 baseline 的关键点（spec 风险表需补一行）：`run_compare.py` 与 `BASELINE.md` 都建立在 `n_defect=0` 路径上；决策 #10（质量加权）只影响 `n_defect>0` 分支，不影响现有 8 个 case 的数值对齐。

### `photonics2/plott.py`

`plot_S_I(p, name, plrange, image, **para)` 签名相同，但分支内部存在真实差异。决策表（依据 `2026-07-11` 复核）：

| # | 维度 | `photonics2/`（当前） | `lineshape_new_ref/`（目标） | **决策** |
|---|------|------------------------|------------------------------|----------|
| 1 | `Shw` 分支视觉 | `ax.fill_between(..., facecolor='lightgrey')` 灰填充带 | 删 fill_between；新增 `ax.plot(... C_omega, color="tomato")` 红色 C(ħω,T) 线 | **取 new_ref**（红线 C_ω） |
| 2 | `Shw` 分支落盘 | 无 | 新增 `np.savetxt('C_omega,S_omega.data', ...)` | **删除该行**（不要落盘） |
| 3 | `Shw` 分支 stem | `ax2.stem(...)` 不带 `use_line_collection` | 加 `use_line_collection=True`（matplotlib 3.1+ 兼容；**3.10+ 已移除**，执行时需 sed 删除该参数） | **取 new_ref 语义意图，但参数 matplotlib ≥ 3.10 不可用**；执行 plan 需 `sed -i 's/, use_line_collection=True//' plott.py` |
| 4 | `Shw` 分支 savefig | `plt.savefig("Shw.png", dpi=500)` **写两次**（line 80-81） | 只写一次 | **取 new_ref**（修 bug） |
| 5 | `PLeV` 分支落盘 | 无 | 新增 `np.savetxt('lineshape_eV.csv', ...)` | **删除该行**（不要落盘） |
| 6 | `PLnm` 分支落盘 | 无 | 新增 `np.savetxt('lineshape_nm.csv', ...)` | **删除该行**（不要落盘） |
| 7 | `PLcm` 分支 `split` 参数 | `split = para.get("split", 0)`（接受外部参数） | 整段删除该行 | **删除**（旧版未实际引用） |
| 8 | `PLnm` 分支变量复用 | 直接 `Iabs[int(...)]` 切片 | 新建 `Iabs_new` 再切片 | **取 new_ref**（避免污染后续） |
| 9 | `Shw+Acm` 内死代码 `xcm` | 计算 `xcm = x/0.000124` 但未使用 | 同 | **取 new_ref**（删除死代码） |
| 10 | 风格 / 空行 / `=` 两侧空格 / 末尾注释 | 旧版风格 | 清理版 | **取 new_ref** |

合并产出：以 `lineshape_new_ref/photonics2/plott.py` 为基线，按决策 2/5/6 删除三处 `np.savetxt`，按决策 7 删除 `split = para.get(...)` 行；其余 7 处直接沿用 new_ref。

合并后文件**与 `lineshape_new_ref/photonics2/plott.py` 字节不同**（少了 4 处）。验收不再用 byte-equal，改用 grep 校验所有决策点。

### `pl.py`（顶层，两个目录）

这一对**确实不同**。决策表（依据 `2026-07-11` 复核）：

| # | 维度 | `src/pl.py` | `src/lineshape_new_ref/pl.py` | **决策** |
|---|------|-------------|--------------------------------|----------|
| 1 | shebang | `#!/home/qqs/miniconda3/envs/pymatgen/bin/python` | `#!/share/home/ckduan/anaconda3/envs/my_pydefect/bin/python` | 替换为 `#!/usr/bin/env python3`（可移植） |
| 2 | `import matplotlib; matplotlib.use('Agg')` | 无 | 有 | **取 new_ref** —— 无显示器运行必需 |
| 3 | `res` 默认值 | 硬编码 `1000`（覆盖 INCAR） | `parameter.get("resolution", 1000)` | **取 new_ref** —— 尊重 INCAR 字段 |
| 4 | emission 的窗口 `Amin/Amax` | `-0.9 * pw / 0.1 * pw` | `-0.8 * pw / 0.2 * pw` | **取 new_ref：-0.8 / 0.2** |
| 5 | absorption 的窗口 `Amin/Amax` | `-0.1 * pw / 0.9 * pw` | `-0.2 * pw / 0.8 * pw` | **取 new_ref：-0.2 / 0.8** |
| 6 | `plot_S_I` flag 字符串 | `"Shw Sk PLeV  PLnm A Phon"`（6 图） | `"Shw PLeV  PLnm "`（3 图） | **取主版：6 图**（`Shw`、`Sk`、`PLeV`、`PLnm`、`A`、`Phon`） |
| 7 | INCAR 默认路径 | `os.path.join(os.path.dirname(__file__), "INCAR")` | `"./INCAR"` | **取主版** —— 目录重整后脚本相对路径仍能解析 |
| 8 | 案例默认路径 | `"../cases/zlq/{band.yaml,GS,ES}"` | `"./zlq/{band.yaml,GS,ES}"` | **取主版** —— 与重整后的 `cases/` 布局一致 |
| 9 | 注释 / 末尾 docstring | 长；若干注释掉的调用模板 + 模块 docstring | 清理后 | **取主版** —— 保留 `masses_kg` 用法提示与 PL 相关注释；删掉长模块 docstring（README 才是正式文档） |

合并产出的 `src/pl.py` 行为按上表第 3–9 行确定；两边各取部分行在 `get_pl` 内混合，不是简单整段替换。

`run_compare.py` **不**调用 `pl.py` —— 它自己构造 driver。所以这里的改动不会被 A/B 回归覆盖；plan 里的冒烟测试会直接 `python src/pl.py` 跑 `cases/zlq/`。

## 目标状态

```
qqs/lineshape/src/
├── photonics2/
│   ├── photoluminescence.py      ← 替换为 lineshape_new_ref 的向量化版本
│   ├── plott.py                  ← 按 "Differences → plott.py" 决策表合并（基线为 lineshape_new_ref 的版本）
│   ├── configuration_coordinate.py
│   ├── xyz.py
│   ├── __init__.py
│   └── {constants,ebedding,embedding,hermite,ht,jt,jtsoc,nu,schrodinger,un}.py  ← 不动
├── lineshape_new_ref/
│   └── {photonics2/, pl.py}      ← 不动（保留作 A/B 回归）
├── pl.py                         ← 按 "Differences → pl.py" 决策表合并后替换
├── phonon_struct_op.py           ← 不动
└── INCAR                         ← 不动
```

被替换的三个文件在替换前各自留一份 `<file>.bak`，验证通过后删除。

## 约束

- **算法等价，非 API 变更。** 公开 `Photoluminescence` API（`__init__`、`HuangRhyes`、`el_ph`、`PL`、`PLA`、`print_table`，属性 `frequencies`、`Modes`、`m`、`S`、`skipmodes`、`numModes`、`numAtoms`、`Delta_R`、`Delta_Q`）不变。
- **不删** `lineshape_new_ref/` —— 它仍作为 `run_compare.py` 的参照实现保留。
- **不删**旧版独占模块（10 个无人维护文件）。超出范围。
- **不改** `phonon_struct_op.py`、`INCAR`、`cases/`。
- **`pl.py` 决策已定**，实现阶段不要回头改。

## 验收标准

1. `run_compare.py` 在全部 8 个案例上退出码为 0。
2. 每一案例："A" 列（现在跑向量化版本）与 "B" 列在 `HuangRhyes`、`Delta_R`、`Delta_Q` 上匹配到 4 位小数。
3. 每一案例：匹配值与 `BASELINE.md` 对应行精度一致。
4. `python src/pl.py` 对 `cases/zlq/` 跑通。**运行前**临时把 `src/INCAR` 里的 `resolution = 4000` 改成 `resolution = 1000`（验证后用 `git checkout -- src/INCAR` 恢复成 4000）。然后跑 `python src/pl.py`，确认产出 6 张 PNG（`Shw.png`、`Sk.png`、`PLev.png`、`PLnm.png`、`AeV.png`、`Phon.png`）。
5. 替换后 `src/plott.py` 通过决策表 grep 校验：决策 1 命中 C_omega 红线、决策 3 `use_line_collection` 计数 = 0（matplotlib ≥ 3.10 已删除）、决策 4 `plt.savefig("Shw.png"` 计数 = 2、决策 2/5/6 三处 `np.savetxt` 不存在、决策 7 `split` 行不存在、决策 8 `Iabs_new` 存在、决策 9 `Shw+Acm` 分支内无 `xcm`。
5b. 替换后 `src/photonics2/photoluminescence.py` 通过决策表 29 条 grep 校验（按"决策—grep 表达式"清单逐项断言）。
6. 跑 `python src/pl.py` 期间，cwd 除已知的 `*.png` 与 `partial.HuangRhyes.data`/`main_modes.data`/`modes.data` 外，**不**出现 `C_omega,S_omega.data`、`lineshape_eV.csv`、`lineshape_nm.csv`。

## 风险

| 风险 | 可能性 | 缓解 |
|------|--------|------|
| 逐模循环与 numpy 广播的浮点累计顺序差异超过 `1e-4` | 低 | numpy 行为确定；广播顺序固定；`run_compare.py` 会暴露漂移 |
| 移除 `multiprocessing.Pool` 启动开销后 wall-clock 反而变长 | 中 | 可接受；用户后续可自行 profile。记为已知行为变化 |
| `get_phonon` 不再 `import yaml` —— 其他脚本若依赖该 yaml 导入可能断 | 未发现 | `grep -r yaml qqs/lineshape/src/` 仅命中旧版 `photonics2.photoluminescence` |
| `read_grd_ex_pos` 新增 `path` 位置参数（决策 #8）—— 外部脚本若按旧签名调用会 TypeError | 低 | 已知 `run_compare.py` 走 `__init__` 间接调用，不受影响；其他外部调用方需在 spec 文档化 |
| 决策 #10（`n_defect` 质量加权）不影响现有 baseline（8 个 case 走 `n_defect=0` 分支） | 低 | 验收 #2/#3 在 baseline 对齐上验证；以后跑 `n_defect>0` 需补 baseline |
| `.bak` 备份残留 | 低 | Task 5 验证通过后删除 |
| `src/INCAR` 临时改 `resolution` 后忘了恢复 | 低 | 验收步骤里显式 `git checkout -- src/INCAR` 收尾 |

## 非目标

- 清理旧版独占模块（`jt.py`、`jtsoc.py` 等）
- `Photoluminescence` API 重设计
- 加重类型提示 / 改写到 `pyphotonics/` 包规约
- 修改 `cases/`

## 后续

暂无。