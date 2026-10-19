# CERPS 机制开关登记表（Mechanism Registry）

> **用途**：把代码中的每一个"机制开关"变成可受控对比的"科学变量"。
> 每次做机制对比实验前，先在此确认：开关位置、默认、依赖/互斥、文献依据、
> 预期效应方向、正交性风险、验证状态。
> **约定**：`[推断]` = 由代码/文献推导的占位结论，需领域判断确认；`[确证]` = 有原文/数值支持。
> 本表与 `cerps/config.py`、`cerps/batch.py` 的 `validate()` 交叉核对。

---

## 轴 1 — 细胞内表观连续动力学 + 四种遗传核  (`EpiConfig`)

| 参数 | 位置 | 默认 | 取值 | 依赖/互斥 | 文献依据 | 预期效应方向 `[推断]` | 正交性风险 | 验证 |
|---|---|---|---|---|---|---|---|---|
| `continuous` | `EpiConfig.continuous` | `False` | bool | 关→`u` 静态继承；开→`u` 随 `act` 演化 | HSL2024 §2.2（epi 仅继承）| 开：表型更具可塑性、稳态 `u` 偏离继承值 | 与 `kernel` 交互：`u` 被动力学改写 | ✅ |
| `w0` / `w1` | `EpiConfig` | `0.01 / 0.06` | ≥0 | 仅 `continuous` 生效 | Berry 2017（擦除∝转录率）、Miyamoto 2015（`τ dθ/dt=(Θ−αx)−θ`）| `w1↑` → 活性基因 `u↑` 更快 | `w1` 与 `grn.lam` 共驱表达 | ✅ |
| `e0` / `e1` | `EpiConfig` | `0.01 / 0.06` | ≥0 | 仅 `continuous` 生效 | 同上 | `e1↑` → 低活性基因 `u↓` 更快 | 与 `w1` 方向相反 | ✅ |
| `kernel` | `EpiConfig.kernel` | `"beta"` | beta\|gaussian\|cpg\|deterministic | 独立 | β 核：Lei 2019 / Huang&Lei 2017-19；高斯：Miyamoto 2015；cpg：Haerter 2013；确定：Furusawa&Kaneko 2001 | — | **⚠️ 见下「内核正交性」** | ✅ |
| `m0,m1,m2,m3,m4` | `EpiConfig` | `60,0.08,0.8,2,2` | — | β 核专用；`m2=0.8` | HSL2024 表 2（注：表 2 列 1.0，论文正文用 0.8，见协议矛盾）| 调参改变 `φ` 形状与不动点 | `m2` 与论文内部矛盾，需记录 | ✅ |
| `sigma_g` | `EpiConfig` | `0.045` | ≥0 | gaussian 专用 | Miyamoto 2015 | `σ_g↑` → 分配噪声↑ | **⚠️ 折叠截断在 u≈0/1 偏差**（已修：反解前置均值）| ✅ |
| `cpg_sites` | `EpiConfig` | `60` | int | cpg 专用 | Haerter 2013（NAR 42:2235）| 位点越多越近连续 | **⚠️ 需 `cpg_sites==m0+1` 才与 β 核方差一致** | ✅ |
| `cpg_mode` | `EpiConfig` | `"moment"` | moment\|mechanistic | cpg 专用 | — | `moment`：位点概率=`φ(v)`，与 β 核矩一致；`mechanistic`：仿射均值映射（固定点 0.1429≠0.1325）| **⚠️ 旧默认 mechanistic 改写了漂移项，使四核非受控对比**（已改为 moment）| ✅ |
| `cpg_p_maintain` / `cpg_p_denovo` | `EpiConfig` | `0.88 / 0.02` | [0,1] | cpg `mechanistic` 专用 | Haerter 2013 规则 | 决定维持/新写概率 | 仅 `mechanistic` 模式，勿与 `moment` 混用 | ✅ |

**内核正交性红线**（跨内核对比前必读）：四种内核只有在 `cpg_mode="moment"` 且 `cpg_sites==m0+1` 时才有相同的均值/方差映射；否则换内核=换确定性漂移，对比不是受控实验。

---

## 轴 2 — 分裂时基因表达分子级分配  (`DivisionConfig`)

| 参数 | 位置 | 默认 | 取值 | 依赖/互斥 | 文献依据 | 预期效应方向 `[推断]` | 正交性风险 | 验证 |
|---|---|---|---|---|---|---|---|---|
| `partition` | `DivisionConfig.partition` | `"exact"` | exact\|binomial | — | HSL2024 §2.3（自陈对称分裂简化）；Berg 1978；Huh&Paulsson 2011 | `binomial`：引入分配噪声，稳态方差↑ | 与 `alpha` 互斥（见下）| ✅ |
| `volume` | `DivisionConfig.volume` | `200.0` | >0 | — | — | 分子数尺度；CV∝1/√n | `n=round(V·x)` 小 x 离散偏差 | ✅ |
| `asymmetry` | `DivisionConfig.asymmetry` | `0.5` | (0,1) | — | Bertaux 2018（体积比例浮动）| 固定 0.5 = 假设子体等大 | 文献标准按体积比例，`p=0.5` 是理想化 | ✅ |
| `alpha` | `DivisionConfig.alpha` | `1.0` | ≥0 | 与 `partition` 协同（exact≡α=0，binomial≡α=1）| Soltani 2016 式(5)：`Var[x⁺|x]=αx/4` | `α=0` 确定性；`α=1` 二项；`α>1` 成簇；`0<α<1` 主动有序 | 二档开关丢失 `α>1` 与 `0<α<1` 区间 | ✅ |

**注意**：`batch.validate` 会警告"exact 模式下 `alpha` 无效"。奇数拷贝数下 exact 无法严格实现，作理想化基准。

---

## 轴 3 — 多基因 GRN（乘积/加和交叉抑制）  (`GRNConfig`)

| 参数 | 位置 | 默认 | 取值 | 依赖/互斥 | 文献依据 | 预期效应方向 `[推断]` | 正交性风险 | 验证 |
|---|---|---|---|---|---|---|---|---|
| `regulation` | `GRNConfig.regulation` | `"product"` | product\|additive | — | product：Duddu 2026（J R Soc Interface 23:20251065 式4.5）逐项同构；additive：HSL2024 | G=2 时两式**严格全等**；G≥3 时拓扑不同（product→**单高主导 `[确证]`**，additive→对称）| **⚠️ G≥3 时两种组合给出不同稳态，不可混用对比**；product 单高主导已为 3-基因扫描实证（docs/multigene_attractors.md：单高~93%）| ✅ |
| `b_rep` | `GRNConfig.b_rep` | `()` = 全互抑 | 矩阵 | 默认全互抑；自定义边 | — | 抑制强度↑ → 靶基因压低 | 与 `s_self` 协同 | ✅ |
| `s_self` | `GRNConfig.s_self` | `(0.5,0.5)` | 每基因阈值 | — | HSL2024 式(1) | 阈值↑ → 更难触发表达 | **⚠️ 交叉项阈值用靶基因 `s[j]`（已修）；s1≠s2 时 legacy 精确退化式(1)，可复现论文图10** | ✅ |
| `c_act` | `GRNConfig.c_act` | `()` | 矩阵 | 与 `b_rep` 并存 | — | 激活边↑ → 靶基因升 | 与 `b_rep` 反向 | ✅ |
| `a_tf` | `GRNConfig.a_tf` | `(0.0,0.0)` | 每基因 | — | HSL2024 式(9) | 外源 TF↑ → 对应基因表达↑ | 仅作用于被指定基因 | ✅ |
| `sigma` / `zeta` | `GRNConfig` | `0.05 / 1.0` | — | 驱动 OU 外噪声 | — | `σ↑` 或 `ζ↓` → 表达噪声↑ | OU 是**马氏**过程（见轴 P6）| ✅ |

**四稳态张力**：3 基因全互抑下"全开+3 单高"格局在同拓扑文献缺支持（Duddu 2020：三阳~4%且被排除）。CERPS 自身 3-基因扫描（docs/multigene_attractors.md）实证：全开仅 ~0.5-1%、单高主导 ~93%、双高尾 ~7%——与 Duddu 统计方向相反，印证"纯互抑网络不能稳健复现全开=干性"隐喻。多能核心实为**互激活**（Boyer 2005），建议分层结构。

---

## 轴 4 — 连续 stemness 依赖的增殖/死亡速率  (`RatesConfig`)

| 参数 | 位置 | 默认 | 取值 | 依赖/互斥 | 文献依据 | 预期效应方向 `[推断]` | 正交性风险 | 验证 |
|---|---|---|---|---|---|---|---|---|
| `rate_mode` | `RatesConfig.rate_mode` | `"legacy"` | legacy\|continuous | — | HSL2024 表 1 | — | legacy=离散分类；continuous=连续 | ✅ |
| `p_switch` | `RatesConfig.p_switch` | `12.0` | >0 | continuous | HSL2024 式(10) 极限 | `p↑` → 越接近 legacy 硬阈值 | `p` 小→过渡带宽，分类边界外移 | ✅ |
| `K_switch` | `RatesConfig.K_switch` | `None` | ≥0 | continuous | — | `None`→自动对齐 legacy 边界（双基因 0.25）| **⚠️ 旧默认 0.5 掩盖 TA/SC 互补阈值 bug（已修 `1-K`）** | ✅ |
| `theta0` / `beta0` | `RatesConfig` | `200 / 0.04` | — | niche 反馈 | HSL2024 §2.2 | `θ0↑` → 反馈更弱；`β0↑` → 增殖更快 | — | ✅ |
| `max_div_ta` | `RatesConfig.max_div_ta` | `15` | int | 衰老 | C++ 代码 | 阈值↑ → TA 可分裂更多次 | — | ✅ |
| `delay_steps` | `RatesConfig.delay_steps` | `0` | ≥0 | **P8 反馈时滞** | —（新增，见轴 P8）| `>0` → niche 信号延迟 τ 步生效 | 默认 0 退化为瞬时，向后兼容 | ⬜ 实现中 |
| `theta_rho` / `rho0` | `RatesConfig` | `12 / 6` | — | 接触抑制 | Van Liedekerke 2019/20 | `ρ0↑` → 抑制更晚 onset | 仅空间模式生效（经 `rho_local`）| ✅ |

**尖锐极限声明收紧**：`p_switch→∞` 只恢复速率**幅值**（表 1）；分类边界 `S=K_switch` 与 legacy 矩形域在 (x1,x2) 平面有 ~8% 不一致（结构必然）。

---

## 轴 5 — niche 反馈 / 局部密度 / 旁分泌场  (`RatesConfig` + `SpaceConfig`)

| 参数 | 位置 | 默认 | 取值 | 依赖/互斥 | 文献依据 | 预期效应方向 `[推断]` | 正交性风险 | 验证 |
|---|---|---|---|---|---|---|---|---|
| `β_SC=β0·θ0/(θ0+Q)` | `rates.legacy_rates` / `rates` | — | — | Q=SC 数（瞬时或 `delay_steps` 延迟）| HSL2024 §2.2 | Q↑ → SC 增殖↓（负反馈）| 与 `delay_steps` 耦合（P8）| ✅ |
| `space.paracrine` | `SpaceConfig.paracrine` | `False` | bool | 需 `space.enabled` | Du,Nie&Holmes 2015（PLoS Comput Biol 11:e1004285，逐点同构）| 开：场驱动 `response` 加到基因表达 | 场与 GRN 的剂量-响应标定缺文献先例 | ✅ |
| `diff_rate` | `SpaceConfig.diff_rate` | `25.0` | >0 | 旁分泌 | Du 2015 | `D↑` → 场扩散更快、更平滑 | **⚠️ 旧 CFL 静默 clamp（已改自适应子步，名义 D 生效）** | ✅ |
| `decay_rate` / `secretion` | `SpaceConfig` | `0.02 / 0.05` | — | 旁分泌 | Du 2015 | 衰减↑→场更局部；分泌↑→场更强 | — | ✅ |
| `response` | `SpaceConfig.response` | `(0,0)` | 每基因 | 旁分泌 | — | 场对基因表达的激活量 | 需剂量-响应标定 | ✅ |
| `normalize_source` | `SpaceConfig.normalize_source` | `False` | bool | 旁分泌 | — | 开：注入除以 h³ 得真浓度 | 影响 `response` 量纲 | ✅ |
| `contact_inhibition` / `theta_rho` / `rho0` | `SpaceConfig` + `RatesConfig` | `True / 12 / 6` | — | 空间 | Van Liedekerke 2019/20 | 邻居超 `ρ0` → 增殖↓ | 与 `beta` 反馈相乘 | ✅ |
| `inhibition_radius` | `SpaceConfig.inhibition_radius` | `1.6` | — | 接触抑制 | — | 半径↑ → 更多邻居计入 | — | ✅ |

**稳定性警示**：Du 2015 原文"Paneth 来源 Wnt 单独既不能维持 niche 也不能促进稳态"提示场–增殖正反馈可能失稳；开启 paracrine 前做稳定性扫描。

---

## 轴 6 — 三维空间离格力学  (`SpaceConfig`)

| 参数 | 位置 | 默认 | 取值 | 依赖/互斥 | 文献依据 | 预期效应方向 `[推断]` | 正交性风险 | 验证 |
|---|---|---|---|---|---|---|---|---|
| `enabled` | `SpaceConfig.enabled` | `False` | bool | — | Drasdo&Höhme 2005；Van Liedekerke 2019/20 | 开：引入位置/力学 | — | ✅ |
| `relax` | `SpaceConfig.relax` | `0.15` | (0,1) | 每步重叠松弛比例（`dt=0.25` 基准）| CBM 标准 | `relax↑` → 力学更快收敛 | 旧 `relax(dt)` 不用 dt（已改 `relax_time`）| ✅ |
| `relax_time` | `SpaceConfig.relax_time` | `None` | >0 | 显式松弛时间 | 过阻尼 `ṙ=F/γ` | 设显式值使力学 `dt` 无关 | `None`→由 `relax` 推导（dt=0.25 行为不变）| ✅ |
| `motility` | `SpaceConfig.motility` | `0.0` | ≥0 | — | — | `↑` → 随机游走更强 | — | ✅ |
| `cell_radius` / `domain_radius` | `SpaceConfig` | `1.0 / 60.0` | — | 反射球边界 | — | 域相对越大→聚集越松 | 反射边界（已改，幽灵信号 8%→1e-6）| ✅ |
| `grid_h` | `SpaceConfig.grid_h` | `4.0` | >0 | 旁分泌场网格 | Du 2015 | 网格粗→CFL 更易超 | 与 `diff_rate/dt` 共同决定 Fo | ✅ |

**CBM 已知缺陷**：Van Liedekerke 2020 批评纯成对力的中心基模型在密聚集块失真；ν=0.5 不可压缩细胞接触面积仅粗略近似。

---

## 机制基元覆盖表（决定可承载的假说）

| 基元 | 原生支持 | 对应开关 | 状态 |
|---|---|---|---|
| **P1 niche 反馈** | ✅ | `theta0`, `β_SC` 公式, 接触抑制, paracrine | 完备 |
| **P2 周期驱动** | ✅ | `tau0`, SC/TA 周期, `max_div_ta` 衰老 | 完备 |
| **P3 组分稀释** | ✅ | `partition`/`alpha`（二项≡α=1）| 完备 |
| **P4 quorum sensing** | ✅ | `paracrine` 扩散场（+CFL 子步）| 完备 |
| **P5 状态继承** | ✅ | 四种核 + `div_count` 继承 | 完备 |
| **P6 非马氏记忆** | ✅ | 新增 `MemoryConfig`（泄漏积分记忆核，见下）| ✅ 已实现+测试 |
| **P7 空间耦合** | ✅ | `space.*` 离格 + 接触抑制 | 完备 |
| **P8 反馈时滞** | ✅ | `RatesConfig.delay_steps` | ✅ 已实现+测试 |

---

## 新增基础组件（本次补全）

### P8 — 反馈时滞  (`RatesConfig.delay_steps`, 默认 0)
- **位置**：`cerps/config.py` `RatesConfig`；`cerps/simulator.py` 维护 `self._q_hist`；`cerps/rates.py` 经 `effective_niche_signal()` 取延迟值。
- **语义**：niche 信号 Q（连续模式 `Q_soft` / legacy 模式 SC 计数）延迟 `delay_steps` 步生效。
- **预期效应方向 `[推断-方向]`**：正延迟 → 反馈环路引入相位滞后，可能放大振荡或改变稳态 SC 占比（类比 Huh&Paulsson 2011）。**机制已落地并测试**（delay_steps=0 退化为基线；>0 时 `self._q_hist` 长度=步数且动力学可测量变化），方向性结论待生物学标定。
- **正交性**：仅改 Q 的时间取值，不改分类/核；与 `rate_mode`、`paracrine` 正交。
- **验证**：`delay_steps=0` 与现状逐位一致；`>0` 时 `self._q_hist` 长度 = 步数且动力学可测量变化。

### P6 — 非马氏记忆  (`MemoryConfig`, 默认 `enabled=False`)
- **位置**：`cerps/config.py` `MemoryConfig`；`cerps/simulator.py` 维护 `self.mem`（每细胞泄漏积分器）；分裂时继承给两子细胞。
- **语义**：每细胞维护驱动信号（默认 stemness `S∈[0,1]`）的指数遗忘平均 `mem`；可选经 `mem_couple` 耦合进 `beta`。
- **预期效应方向 `[推断-方向]`**：`mem_couple>0` 且 mem 偏高（细胞近期多处于干性态）→ 增殖倾向增强；提供跨当前时刻的历史依赖（非马氏成分）。**机制已落地并测试**（enabled=False 退化为基线；True 时 mem 演化并被子细胞继承）。
- **正交性**：`enabled=False` 时 `mem≡0`、速率不变（向后兼容）；与 P8 可叠加（先记忆、后时滞）。
- **验证**：`enabled=False` 结果与基线一致；`enabled=True` 时 `mem` 随驱动演化并被子细胞继承；`mem_couple≠0` 时动力学改变。

---

## 使用铁律（实验前必读）
1. **永远保留 `legacy_config()` 作校准零点**；任何新机制的效应表述为"相对该基准的偏离量"。
2. **跨内核/跨 `regulation` 对比前**，先确认"内核正交性红线"与"G≥3 topological 差异"两条。
3. **每次批量前跑 `batch.health_check`**：dt 减半一致性 + CFL 合规（dt 0.25→0.125 可使稳态总数变 ~8.5%）。
4. **记录本表未覆盖的自定义开关**到此处，避免"代码开关≠科学变量"。
