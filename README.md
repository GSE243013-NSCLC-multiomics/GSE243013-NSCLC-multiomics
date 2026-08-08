# GSE243013 非小细胞肺癌多组学分析
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21816175.svg)](https://doi.org/10.5281/zenodo.21816175)
## Archived Release

The exact version of the code accompanying the manuscript is archived
in Zenodo:

- Version: v1.0.0
- DOI: https://doi.org/10.5281/zenodo.21816175
- Release date: 6 August 2026

 1. 研究概述

本代码库包含一项针对接受新辅助抗PD-1免疫治疗的非小细胞肺癌（NSCLC）患者免疫异质性单细胞RNA测序（scRNA-seq）研究的分析代码。该研究利用GSE243013数据（243名患者，1,254,749个细胞）来识别与病理学无应答相关的免疫组分表达程序。

**Key findings**
The all-immune glycolysis program was enriched toward pathological nonresponse in the primary analysis (NES = −2.3589, FDR = 3.00×10^−11) and in the nested strict chemo-immunotherapy sensitivity cohort (NES = −2.4126, FDR = 7.51×10^−12).
Across eight prespecified immune-cell strata, all eight NES values were negative, indicating concordant directionality, whereas only one of eight strata met FDR < 0.05.
In TCGA, the glycolysis score was associated with overall survival in LUAD (HR = 1.46, 95% CI 1.24–1.73, P < 0.001) but not in LUSC (HR = 1.02, 95% CI 0.89–1.18, P = 0.76).
Across 145 program-level cross-cohort tests, four met BH FDR < 0.05.
Because heterogeneity for the all-immune glycolysis association was high (I2 = 90.4%), the pooled fixed-effect estimate is considered descriptive rather than confirmatory.

## 2. 仓库范围

此仓库包含：
- 计算流程的所有分析脚本（步骤 00-09A）
- 手稿生成与审计脚本 (M1-M7C)
- 配置模板和参数定义
 支持手稿的小型最终结果表
- 分析工作流和数据字典文档

此仓库不**包含：
- 原始患者水平单细胞RNA测序数据（必须从GEO下载）
- TCGA 多组学数据（通过 curatedTCGAData 下载）
- 大型中间结果（BPCells 矩阵、edgeR 对象、评分矩阵）
- 患者级别清单（从 GEO 元数据本地重建）
- 完成 renv.lock (参见 [依赖限制](docs/DEPENDENCY_LIMITATIONS.md)

## 3. 主要发现

 **糖酵解程序：**免疫区糖酵解相关表达程序满足预设的内部整合证据标准（meta-FDR < 0.05）
2. **八种细胞类型：**所有8种评分细胞类型均呈负NES（FDR < 0.05）
3. **TCGA-LUAD 外部关联：** HR = 1.46 (CI: 1.24-1.73, P = 5.68e-06)
4. **TCGA-LUSC:** 无显著关联 (HR = 1.02, P = 0.750)
5. **异质性：** I-squared 约为 90%；固定效应汇总估计具有描述性
6. **** LUAD 甲基化（30 个 CpG 位点，最高 cg02952918），RPPA（86 种蛋白质，最高 Cyclin B1）；未生成程序水平的 CNV 关联结果
7. **CollecTRI TF 推理：** 尚未完成任何细胞类型的分析

## 4. 目录结构

```
GSE243013-NSCLC-multiomics/
├── config/                    # 配置文件模板
│   ├── paths.example.yml      # 项目路径模板
│   ├── analysis_parameters.yml # 所有分析参数
│   └── *.txt                  # 特定步骤的分析定义
├── 脚本/
│   ├── analysis/              # 步骤 00-09A 分析脚本
│   ├── manuscript/            # M1-M7C 手稿脚本
│   └── run_pipeline.R         # 管道运行器
├── 数据/
│   ├── README.md              # 数据下载说明

├── 资源/
│   └── resource_versions.tsv  # 外部资源版本
├── 结果/
│   └── final_tables/          # 小型最终结果表
├── manuscript/                # 手稿输出文件
├── docs/                      # 文档
├── tests/                     # 流水线验证测试
├── reproduce.sh               # Main reproduction script
├── README.md                  # This file
├── LICENSE                    # MIT License (code)
├── CITATION.cff               # Citation metadata
└── .gitignore                 # Git ignore rules
```

## 5. Data Availability

- **GSE243013:** Available from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE243013). Four supplementary files are required (see [data/README.md](data/README.md)).
- **TCGA:** Available via curatedTCGAData R package (version 1.34.0). LUAD and LUSC are analyzed separately.
- **MSigDB:** MSigDB Hallmark gene sets (v2024.1.Hs) were used for the pathway–enrichment results reported in the manuscript.
Reactome results are not reported in the manuscript.
- **Patient manifest:** Reconstructed locally from downloaded GEO metadata; not redistributed in this repository.

## 6. Hardware and Storage Requirements

- **RAM:** 16 GB minimum (32 GB recommended)
- **Disk:** ~110 GB free space (including raw data)
- **CPU:** Multi-core recommended for edgeR and fgsea
- **OS:** macOS (ARM64) or Linux (x86_64/ARM64)

See [docs/HARDWARE_REQUIREMENTS.md](docs/HARDWARE_REQUIREMENTS.md) for details.

## 7. Software Requirements

Tested with R 4.6.0. Key packages:

| Package | Purpose | Source |
|---------|---------|--------|
| edgeR | Patient-level pseudobulk DE | Bioconductor |
| fgsea | Preranked pathway enrichment | Bioconductor |
| BPCells | Column-major sparse matrix | GitHub |
| metafor | Fixed-effect meta-analysis | CRAN |
| survival | Cox PH models | CRAN |
| curatedTCGAData | TCGA data access | Bioconductor |

See [docs/DEPENDENCY_LIMITATIONS.md](docs/DEPENDENCY_LIMITATIONS.md) for full list and renv status.

## 8. Environment Restoration

```bash
git clone https://github.com/GSE243013-NSCLC-multiomics/GSE243013-NSCLC-multiomics.git
cd GSE243013-NSCLC-multiomics
export GSE243013_PROJECT_ROOT="$(pwd)"
cp config/paths.example.yml config/paths.yml
# Edit config/paths.yml for your environment
```

## 9. Data Download

Before running the pipeline, download 4 files from GEO:
1. `GSE243013_NSCLC_immune_scRNA_counts.mtx.gz`
2. `GSE243013_NSCLC_immune_scRNA_metadata.csv.gz`
3. `GSE243013_genes.csv.gz`
4. `GSE243013_barcodes.csv.gz`

See [data/README.md](data/README.md) for detailed instructions.

TCGA data is downloaded automatically by Step 08A via `curatedTCGAData`.

## 10. Full Execution Order

**Optional prerequisites** (run only if needed):
```bash
bash reproduce.sh --step 04A    # BPCells binary install
bash reproduce.sh --step 07A   # MSigDB extraction install
# Or include both:
bash reproduce.sh --include-install-steps
```

**Main pipeline:**
```bash
bash reproduce.sh              # Run all steps (00-09A)
bash reproduce.sh --from-step 05 --to-step 07  # Partial run
bash reproduce.sh --dry-run    # Show what would run
```

| Step | Description | Runtime |
|------|-------------|---------|
| 00 | R environment check | Short |
| 01 | Metadata download and inspection | Short |
| 02 | Patient manifest construction | Long |
| 03 | Cohort definition repair | Long |
| 04 | Count matrix import (BPCells) | Long |
| 05 | Patient-level pseudobulk | Long |
| 06 | edgeR differential expression | Long |
| 07 | Pathway enrichment (fgsea) | Long |
| 08A | TCGA data download | Long |
| 08B1 | TCGA clinical modeling (Cox) | Long |
| 08B2 | TCGA multi-omics integration | Long |
| 09 | Final audit and evidence report | Long |
| 09A | Gap closure and revision | Long |

## 11. Expected Outputs

Key outputs in `results/final_tables/`:
- `Table_1_cohort_characteristics.csv` - Patient cohort summary (aggregate counts only)
- `Table_3_primary_edgeR_results.csv` - Primary DE results
- `Table_4_pathway_TF_programs.csv` - Pathway and TF programs
- `Table_5_TCGA_clinical_validation.csv` - TCGA clinical models (external association)
- `GSE243013_final_evidence_tiers_revised.csv` - Evidence tier classification
- `GSE243013_cohort_counts_public.csv` - Aggregate cohort counts

## 12. Reproducing Manuscript Tables and Figures

```bash
# Run manuscript generation scripts (M1-M7C)
bash reproduce.sh --from-step M1 --to-step M7C
```

## 13. Known Limitations

- **CollecTRI:** TF inference not completed for any cell type
- **PROGENy:** Pathway activity inference not completed
- **renv:** Full renv.lock not generated due to binary package dependencies
- **BPCells:** Binary installation required for ARM64 macOS
- **TCGA:** Not an immunotherapy cohort; provides expression-clinical correlation data only
- **CNV:** No program-level association result was generated
- **Heterogeneity:** I-squared approximately 90%; fixed-effect pooled estimate is descriptive

## 14. Citation

```bibtex
@software{Xie2026GSE243013,
  title = {GSE243013 NSCLC Multi-Omics Analysis},
  author = {Xie, Xiaowei and Guo, Ying and Zhang, Yang},
  year = {2026},
  url = {https://github.com/GSE243013-NSCLC-multiomics/GSE243013-NSCLC-multiomics}
}
```

See [CITATION.cff](CITATION.cff) for full citation metadata.

## 15. License

- **Code:** MIT License (see [LICENSE](LICENSE))
- **Data:** Original data providers retain all rights. See GEO and TCGA data use agreements.
- **Third-party resources:** MSigDB, curatedTCGAData, and other packages retain their original licenses.

## 16. Contact

- **Corresponding author:** Xiaowei Xie
- **Email:** xiexiaowei@fjmu.edu.cn
- **ORCID:** 0009-0009-4221-2777
