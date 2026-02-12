# MCO Mn(II)-oxidation PU learning toolkit

用于从多铜氧化酶（MCO）序列中，对 **Mn(II)-oxidizing 倾向**进行排序与不确定性评估。该项目将问题定义为 **Positive-Unlabeled (PU) 学习**，不会把全部 unlabeled 序列硬当作负例。

## 1) 安装与运行示例

### Linux / macOS

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

### Windows (PowerShell)

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Windows (CMD)

```bat
py -3 -m venv .venv
.venv\Scripts\activate.bat
pip install -r requirements.txt
```


### Windows 快速安装（推荐，不需要激活脚本）

```powershell
./scripts/windows_setup.ps1
# 或（仓库根目录包装脚本）
./windows_setup.ps1
```

```bat
REM CMD
scripts\windows_setup.bat
REM 或（仓库根目录包装脚本）
windows_setup.bat
```

安装后可直接使用 venv 里的 Python（无需 `Activate.ps1`）：

```powershell
.\.venv\Scripts\python.exe -m mco_mnox.cli --help
.\.venv\Scripts\python.exe -m streamlit run mco_mnox/gui.py
```

训练：

```bash
python -m mco_mnox.cli train \
  --pos positives.fasta \
  --unl unlabeled.fasta \
  --neg negatives.fasta \
  --outdir runs/run1 \
  --config configs/default.yaml
```

预测：

```bash
python -m mco_mnox.cli predict \
  --model runs/run1/model.pkl \
  --fasta unlabeled.fasta \
  --out runs/run1/results.tsv
```

快速 toy 运行：

```bash
python -m mco_mnox.cli train --pos data/toy/positives.fasta --unl data/toy/unlabeled.fasta --neg data/toy/negatives.fasta --outdir runs/toy --config configs/default.yaml
python -m mco_mnox.cli predict --model runs/toy/model.pkl --fasta data/toy/unlabeled.fasta --out runs/toy/results.tsv
```


## GUI 使用（Streamlit）

安装依赖后可直接启动本地图形界面：

```bash
python scripts/run_gui.py
# 或
python -m streamlit run mco_mnox/gui.py
```

Windows 下也可以直接运行：

```powershell
python scripts/run_gui.py
# 或 PowerShell 脚本
./scripts/run_gui.ps1
```

```bat
REM CMD
scripts\run_gui.bat
```

GUI 提供两个页面：
- **训练**：填写 positives/unlabeled/(可选)negatives、策略、embedder、输出目录后训练；
- **预测**：加载 `model.pkl` + FASTA，输出 `results.tsv` 与 `results.json`，并在界面中展示表格与候选优先级。


## Windows 使用建议（重点）

- 路径可使用 `C:\data\positives.fasta` 或 `C:/data/positives.fasta`，程序内部使用 `pathlib` 处理。
- 如果 PowerShell 执行策略拦截激活脚本，可先执行：

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

- 若无法安装 `streamlit`（公司代理/离线环境），CLI 仍可正常使用；GUI 仅在安装 `streamlit` 后可用。

- 若 `py -3 -m venv .venv` 报 `Could not find platform independent libraries <prefix>`：
  1) 先试 `python -m venv .venv`；
  2) 执行 `py -0p` 检查 Python 安装路径是否异常；
  3) 建议重装官方 Python（勾选 *Add python.exe to PATH*），或用 `winget install Python.Python.3.11`。
- 若 PowerShell 提示 `禁止运行脚本`（ExecutionPolicy）：
  - 临时仅当前窗口放行：

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
```

  - 或长期当前用户：

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

  - 也可以完全绕过激活脚本，直接使用 `.\.venv\Scripts\python.exe ...`（推荐）。


- 若提示 `./scripts/windows_setup.ps1`（或 `./windows_setup.ps1`）无法识别，通常是**当前目录不是仓库根目录**或文件不存在。先执行：

```powershell
Get-Location
Get-ChildItem .
Get-ChildItem .\scripts\windows_setup*
```

  如果找不到脚本，请切换到仓库根目录后重试，或直接运行根目录包装脚本：

```powershell
.\windows_setup.ps1
```

- 在 PowerShell 中调用可执行文件建议使用调用运算符 `&`，例如：

```powershell
& .\.venv\Scripts\python.exe -m mco_mnox.cli --help
& .\.venv\Scripts\python.exe -m streamlit run mco_mnox/gui.py
```

  如果提示 `\.venv\Scripts\python.exe 无法识别`，说明 venv 尚未成功创建，先运行 `windows_setup.ps1/.bat`。


### 真实 ESM2 使用（离线/在线）

当前代码已支持真实 `esm2_t12_35M` / `esm2_t33_650M` 推理。建议先在 `configs/default.yaml` 中设置：

```yaml
features:
  embedder: esm2_t33_650M
  esm_local_dir: "D:/models/esm2_t33_650M_UR50D"  # 离线目录，留空则在线拉取
  esm_force_local: true
  esm_device: auto
  esm_batch_size: 2
  esm_max_length: 1024
```

- 离线模式：把模型目录完整下载到本地后，设置 `esm_local_dir` + `esm_force_local: true`。
- 在线模式：`esm_local_dir` 留空，`esm_force_local: false`，程序会按模型名在线下载并缓存。

中国大陆建议：
- 优先离线模式（最稳定）；
- 在线模式可配代理或镜像；
- 如需 HuggingFace 镜像，可尝试在环境变量中设置 `HF_ENDPOINT`。

## 2) 如何准备 positives/unlabeled 数据

- `positives.fasta`：实验确认具有 Mn(II)-oxidizing 活性的 MCO（通常约 20 条）。
- `unlabeled.fasta`：尚未验证活性的大量 MCO。
- `negatives.fasta`（可选）：明确不具备 Mn 氧化活性的强负例；可为空。

建议：
- 统一氨基酸大写、去除非法字符。
- 用非冗余策略（如 mmseqs/cd-hit）减少高度同源重复。
- 保留 metadata（物种、定位、实验条件）用于后续分析。

## 3) 特征与 PU 建模说明

### 特征构成
最终特征 = `PLM embedding` + `可解释手工特征`。

- Embedding：CLI 保留 `esm2_t33_650M | esm2_t12_35M | protT5 | none` 选项。
  - 当前默认/离线稳定模式使用 `none`（k-mer 向量 fallback）以保证 CPU 可运行。
  - 同一序列 embedding 会缓存到 JSON，避免重复计算。
- 可解释特征：
  - 四个铜结合 motif 的命中、位置、间距、完整度；
  - T1 位点邻域代理（疏水性、芳香族比例、MLF比例）；
  - 信号肽/Tat (`RR..FLK`) 与跨膜螺旋启发式；
  - 酸性富集（D/E 比例、酸性簇）、估计 pI；
  - His/Cys 密度、低复杂度比例、长度等。

### PU 策略（可切换）
配置文件 `pu.strategy`：
- `bagging`：Bootstrap PU，多轮抽样 unlabeled 子集作为临时负例，输出 `score_mean` + `score_std`。
- `rn`：Reliable Negative 两阶段，先弱监督选 RN，再用 P vs RN 训练。

## 4) 输出字段解释与不确定性

`results.tsv` 核心字段：
- `score_mean`：PU 集成后 Mn-oxidizing 倾向分数（用于排序）。
- `score_std`：集成标准差，不确定性估计（越高越不确定）。
- `rank`：基于 `score_adjusted`（对非典型 MCO motif 完整度不足者降权）排名。
- `nearest_positive` + `nearest_positive_cosine`：embedding 空间最近正例证据。
- `motif*_hit`、`mco_motif_completeness`：MCO motif 证据。
- `signal_peptide_flag`、`tat_rr_flag`、`tm_helix_flag`：定位/膜锚定代理。
- `acidic_ratio`、`estimated_pI`：酸性与理化代理。
- `noncanonical_mco_flag`：可能非典型 MCO/注释错误。
- `explanation`：线性模型系数 top features（全局解释）。

同时输出 `results.json` 三类实验优先级列表：
1. 高分低不确定性（优先验证命中）；
2. 高分高不确定性（优先验证以提升模型）；
3. 中分但 motif/定位证据强（边界候选）。

## 5) 防泄漏评估建议

项目提供 k-mer Jaccard 近似聚类与 GroupKFold 切分接口，避免同簇泄漏；可替换为 mmseqs/cd-hit 聚类标签。评估应以排序指标为主（PR-AUC、top-k hit rate、不同 seed 稳定性），并明确这是 PU/弱监督场景。

## 6) 如何在获得新标签后增量更新模型（Active Learning）

1. 每轮从 `results.tsv` 中挑选：
   - 高分低不确定（确认命中）；
   - 高分高不确定（最大信息增益）；
   - 边界案例（校正决策边界）。
2. 将新实验结果加入：
   - 真阳性追加到 `positives.fasta`；
   - 明确阴性追加到 `negatives.fasta`。
3. 重新训练并对比：
   - top-N 命中、排名稳定性、score_std 收敛。
4. 周期迭代，逐步提高模型区分能力。

## 7) 离线模型权重准备（ESM2/ProtT5）

当前代码已内置 ESM2 (`esm2_t12_35M` / `esm2_t33_650M`) 推理入口。
- 离线：下载模型到本地后配置 `features.esm_local_dir`，并设置 `features.esm_force_local: true`；
- 在线：`esm_local_dir` 留空，程序自动下载并缓存；
- 若缺少依赖（`torch`/`transformers`）会提示安装，或可回退 `embedder: none`。

> `protT5` 目前仍为占位接口，默认回退到 k-mer embedding；如需我也可继续补齐 ProtT5 本地推理。

## 8) 可复现性

- 固定随机种子（配置项 `seed`）。
- 保存运行配置到 `config.used.json`。
- 训练日志输出到 `train.log`。

## 9) 测试

```bash
pytest -q
```
