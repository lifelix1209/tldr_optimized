# De Novo 评估逻辑说明（当前实现）

本文档描述 `scripts/parent_support.py` 在 `--denovo` 模式下的实际评估逻辑（以当前代码为准）。

## 1. 整体流程

1. 读取候选位点（`--candidates`，支持 `tldr` 表格或 JSON）。
2. 对每个候选位点在父母 BAM 中统计：
   - 覆盖深度（`mom_depth_mean` / `dad_depth_mean`）
   - alt 证据（soft-clip + CIGAR insertion）
   - alt 的左右断点分布与位置离散度
3. 在 `--denovo` 开启时，调用 `evaluate_denovo(...)` 给出分类：
   - `PASS_DENOVO`
   - `UNCERTAIN`
   - `FAIL`

## 2. 输入与标准化

候选记录由 `_normalize_candidate` 归一化，关键字段：

- `uuid`（来自 `UUID`）
- `chrom`（`Chrom` 或 `Chromosome`）
- `bp_left` / `bp_right`（优先 `bp_left` / `bp_right`，否则回退 `Start` / `End`）
- `wiggle`（若候选自带则优先用于父母侧搜索窗口）
- `te_family`（`TE_family` 或 `Family`）
- `child_support`（原样保留）

坐标校验规则：

- `bp_left`、`bp_right` 必须可解析为整数
- 坐标不能为负
- 必须满足 `bp_right >= bp_left`（允许单点断点）

不满足则候选被丢弃。

候选预过滤（`analyze_parent_support`）：

- 默认会丢弃 `TE_family` 为 `NA` 的候选（可通过 `--allow_na_te_family` 关闭）
- 可选 `--require_filter_pass`，仅分析 `Filter=PASS` 的候选

## 3. 父母证据统计（每个位点、每个父母）

### 3.1 深度统计

函数：`calculate_window_depth(...)`

- 在左右断点各取 `window_size` 窗口做 pileup。
- 左右窗口深度取中位数，最终 `depth_mean` 为两侧合并后的中位数。

### 3.2 alt 证据统计

函数：`calculate_alt_support_parent(...)`

alt 定义由两部分组成：

- Soft-clip 读段（`find_softclip_reads`）
- CIGAR insertion（`find_insertion_cigar_reads`，`op == 1`）

统计字段包括：

- `softclip_count` / `insertion_count` / `total_alt`
- `total_left` / `total_right`（靠近左/右断点的 alt 数）
- `has_both_sides`（`total_left > 0 and total_right > 0`）
- `bp_position_std` / `bp_position_iqr`（alt 位置离散度）

说明：

- Soft-clip 与 insertion 都会被归类到更近的左或右断点。
- `wiggle` 控制“靠近断点”的窗口范围。
- 若候选包含 `wiggle` 字段，则父母侧优先使用候选 `wiggle`；否则使用 CLI `--wiggle`（默认 80）。
- 若 `--te_enhance` 且提供 `--te_fasta`，soft-clip 需通过轻量 k-mer 匹配才计数。
- v1 增加质量门控：`min_softclip_len` / `min_insertion_len` / `min_mapq`。
- v1 默认按 read 去重统计 `total_alt`（同一 read 的多种证据不会重复计入总 alt）。

## 4. de novo 判定逻辑（`evaluate_denovo`）

## 4.1 自适应父母深度阈值

先从所有候选的 `mom_depth_mean` / `dad_depth_mean`（>0）计算总体中位数 `median_depth`，再计算：

`min_parent_depth = max(min_parent_depth_default, min_parent_depth_frac * median_depth)`

默认参数：

- `min_parent_depth_default = 8`
- `min_parent_depth_frac = 0.25`

## 4.2 child support 解析

从 `child_support` 解析 `child_support_count`，优先识别：

- `useable:<n>` / `usable:<n>`

解析不到则为 `None`（不会触发 `m_child` 阈值判定）。

## 4.3 决策顺序（重要）

### 步骤 A：先看父母 alt 是否应直接判 `FAIL`

当某个父母 `alt > max_parent_alt`（默认 `0`）时，进一步判断其是否“足够一致”：

- `side_ok`：
  - 若未开启 `--require_both_sides`，恒为真
  - 若开启，则要求 `has_both_sides == True`
- `pos_ok`：
  - 固定阈值：`bp_position_std_threshold`（默认 15）、`bp_position_iqr_threshold`（默认 20）
  - v1 默认启用 wiggle 自适应下限：
    - `std_threshold = max(bp_position_std_threshold, wiggle * bp_std_wiggle_factor)`（默认系数 0.25）
    - `iqr_threshold = max(bp_position_iqr_threshold, wiggle * bp_iqr_wiggle_factor)`（默认系数 0.5）

若 `side_ok && pos_ok`，该父母触发 `FAIL` 原因（如 `Mom alt too high (x)`）。

若不满足一致性（单侧/离散过大），不直接 `FAIL`，而是追加 `UNCERTAIN` 原因（如 `alt present but not supported on both breakpoints`、`position-diffuse`）。

### 步骤 B：若没有 `FAIL` 原因，再看 `UNCERTAIN` 条件

会累计以下 `UNCERTAIN` 原因：

- 父母覆盖不足（任一或双方 `< min_parent_depth`）
- `child_support_count < m_child`（仅在可解析时生效，默认 `m_child=1`）
- 步骤 A 中产生的一致性不足原因（单侧/离散）

### 步骤 C：最终分类

- 若存在 `FAIL` 原因：`evaluation = FAIL`
- 否则若存在任意 `UNCERTAIN` 原因：`evaluation = UNCERTAIN`
- 否则：`evaluation = PASS_DENOVO`

## 5. CLI 参数（与 de novo 评估直接相关）

在 `scripts/parent_support.py --denovo` 下可用：

- `--wiggle`（默认 80，候选 `wiggle` 优先）
- `--min_softclip_len`（默认 30）
- `--min_insertion_len`（默认 50）
- `--min_mapq`（默认 20）
- `--no_dedup_by_read`（关闭 read 去重计数）
- `--allow_na_te_family`
- `--require_filter_pass`
- `--m_child`（默认 1）
- `--min_parent_depth_default`（默认 8）
- `--min_parent_depth_frac`（默认 0.25）
- `--max_parent_alt`（默认 0）
- `--require_both_sides`（默认关闭）
- `--bp_position_std_threshold`（默认 15）
- `--bp_position_iqr_threshold`（默认 20）
- `--no_adaptive_bp_thresholds`
- `--bp_std_wiggle_factor`（默认 0.25）
- `--bp_iqr_wiggle_factor`（默认 0.5）

TE 增强相关：

- `--te_enhance`
- `--te_fasta`

## 6. 输出关键字段（de novo 结果）

每条记录会包含：

- 基础位点：`uuid`, `chrom`, `bp_left`, `bp_right`, `te_family`
- 统计上下文：`analysis_wiggle`
- child：`child_support`, `child_support_count`
- 父母深度与 alt：`mom_depth`, `dad_depth`, `mom_alt`, `dad_alt`
- 父母左右/一致性：`mom_total_left/right`, `dad_total_left/right`, `mom_has_both_sides`, `dad_has_both_sides`
- 父母位置离散度：`mom_bp_position_std/iqr`, `dad_bp_position_std/iqr`
- 判定阈值：`bp_std_threshold`, `bp_iqr_threshold`
- 判定：`evaluation`, `reasons`

## 7. 现实现行为注意点

- `m_child` 仅在 `child_support` 可解析为数值时生效；解析失败不会自动判 `UNCERTAIN`。
- `FAIL` 优先级高于 `UNCERTAIN`：只要有任一父母形成“coherent alt”即直接 `FAIL`。
- `--require_both_sides` 仅影响“alt 是否足以判 FAIL”，不影响 alt 计数本身。
