# Hardware Requirements

## Minimum Requirements

| Resource | Minimum | Recommended |
|----------|---------|-------------|
| RAM | 16 GB | 32 GB |
| Disk | 110 GB free | 150 GB free |
| CPU | 2 cores | 8+ cores |
| OS | macOS ARM64 | macOS ARM64 or Linux x86_64 |

## Detailed Requirements by Step

### Step 00-03: Metadata Processing
- **RAM:** 4 GB
- **Disk:** 1 GB
- **CPU:** 1 core
- **Runtime:** < 5 minutes

### Step 04: BPCells Import
- **RAM:** 8 GB
- **Disk:** 50 GB (temporary for count matrix download)
- **CPU:** 2 cores
- **Runtime:** 30-60 minutes

### Step 05: Pseudobulk
- **RAM:** 8 GB
- **Disk:** 5 GB
- **CPU:** 2 cores
- **Runtime:** 10-30 minutes

### Step 06: edgeR DE
- **RAM:** 16 GB
- **Disk:** 5 GB
- **CPU:** 4+ cores
- **Runtime:** 1-2 hours

### Step 07: Pathway Enrichment
- **RAM:** 8 GB
- **Disk:** 2 GB
- **CPU:** 2 cores
- **Runtime:** 30-60 minutes

### Step 08A: TCGA Download
- **RAM:** 4 GB
- **Disk:** 5 GB
- **CPU:** 1 core
- **Runtime:** 10-30 minutes

### Step 08B1-B2: TCGA Analysis
- **RAM:** 8 GB
- **Disk:** 2 GB
- **CPU:** 2 cores
- **Runtime:** 30-60 minutes

### Step 09-09A: Final Audit
- **RAM:** 4 GB
- **Disk:** 1 GB
- **CPU:** 1 core
- **Runtime:** 5-10 minutes

## Storage Breakdown

| Component | Size |
|-----------|------|
| GSE243013 raw data | 7 GB |
| BPCells format | 20 GB |
| TCGA data | 3 GB |
| Intermediate results | 30 GB |
| Final results | 1 GB |
| Manuscript | 1 GB |
| **Total** | **~62 GB** |

## Platform Notes

### macOS ARM64
- BPCells requires binary installation (see Step 04A)
- R libraries installed in `~/Library/R/arm64/4.6/library`

### Linux x86_64
- BPCells can be compiled from source
- R libraries installed in standard library path

### Windows
- Not tested; may require WSL2
- Path handling may differ
