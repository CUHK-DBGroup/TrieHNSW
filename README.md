# TrieHNSW

A label-filtered Approximate Nearest Neighbor (ANN) search method that combines a **Trie** structure with **HNSW** (Hierarchical Navigable Small World) graphs. The Trie organizes vectors by their label sets, while per-node HNSW graphs enable efficient nearest neighbor search under label constraints.

Supported query types: **containment**, **overlap**, and **equality**.

## Prerequisites

- C++17 compiler (GCC 7+ or Clang 5+)
- CMake 3.22+
- OpenMP
- CPU with AVX/AVX-512 support

## Build

```bash
mkdir build && cd build
cmake ..
make
```

The compiled binary `LabelFilteredANNS` will be generated in the `build/` directory.

## Data Format

- **Vector files** (`.fvecs`): Binary format. Each vector is stored as a 4-byte integer (dimension) followed by `dim` floats.
- **Label files** (`.csv`): CSV text format. Each line contains comma-separated integer label IDs for the corresponding vector. Example:
  ```
  1,2,3,4
  2,3,4,5
  4
  ```
- **Ground truth files** (`.bin`): Binary format. Each entry is a `(uint32_t id, float distance)` pair. The file contains `num_queries * K` entries.

## Quick Start

A sample dataset (`data/siftsmall`) is included for testing. It contains 10,000 base vectors, 100 queries, 128 dimensions, and 5 labels.

### 1. Build Index

```bash
./build/LabelFilteredANNS \
    --mode build \
    --base_path data/siftsmall_base.fvecs \
    --query_path data/siftsmall_query.fvecs \
    --base_label_path data/base_labels.csv \
    --query_label_path data/query_labels.csv \
    --num_base 10000 \
    --nq 100 \
    --dim 128 \
    --m 32 \
    --ef 100 \
    --label_num 6 \
    --alpha 1.0 \
    --dataset siftsmall
```

Expected output:

```
./index/siftsmall/
Index built and saved to: ./index/siftsmall/siftsmall_nb10000_dim128_m32_ef100_label6_alpha1.000000.ind2
```

### 2. Search

```bash
./build/LabelFilteredANNS \
    --mode search \
    --query_path data/siftsmall_query.fvecs \
    --query_label_path data/query_labels.csv \
    --gt_path data/siftsmall_gt.bin \
    --query_type containment \
    --result_file result.csv \
    --num_base 10000 \
    --nq 100 \
    --dim 128 \
    --m 32 \
    --ef 100 \
    --label_num 6 \
    --alpha 1.0 \
    --dataset siftsmall
```

Expected output:

```
data/siftsmall_gt.bin,containment,10,30070.9,0.964
data/siftsmall_gt.bin,containment,15,32613.4,0.985
data/siftsmall_gt.bin,containment,20,27983.6,0.995
Results saved to: result.csv
```

### 3. Insert

Build from 80% of the data, then incrementally insert the remaining 20%:

```bash
./build/LabelFilteredANNS \
    --mode insert \
    --base_path data/siftsmall_base.fvecs \
    --query_path data/siftsmall_query.fvecs \
    --base_label_path data/base_labels.csv \
    --query_label_path data/query_labels.csv \
    --gt_path data/siftsmall_gt.bin \
    --query_type containment \
    --result_file result_insert.csv \
    --num_base 10000 \
    --nq 100 \
    --dim 128 \
    --m 32 \
    --ef 100 \
    --label_num 6 \
    --alpha 1.0 \
    --build_ratio 0.8 \
    --dataset siftsmall_insert
```

Expected output:

```
Insert done. build_ratio=0.8 insert_count=2000 insert_time=0.689 insert_qps=2903.49
data/siftsmall_gt.bin,containment,80,2000,0.689,2903.49,10,26943.9,0.963
data/siftsmall_gt.bin,containment,80,2000,0.689,2903.49,15,26568.6,0.987
data/siftsmall_gt.bin,containment,80,2000,0.689,2903.49,20,24054.5,0.995
Results saved to: result_insert.csv
```

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `--mode` | `search` | Operation mode: `build`, `search`, or `insert` |
| `--dataset` | `sift` | Dataset name (used for index file naming) |
| `--num_base` | `1000000` | Number of base vectors |
| `--nq` | `10000` | Number of query vectors |
| `--dim` | `128` | Vector dimensionality |
| `--m` | `32` | HNSW max connections per node |
| `--ef` | `100` | HNSW construction ef parameter |
| `--label_num` | `52` | Number of distinct labels |
| `--alpha` | `1.0` | Trie construction alpha parameter |
| `--build_ratio` | `0.8` | Fraction of data for initial build (insert mode only) |
| `--query_type` | - | Query filter type: `containment`, `overlap`, or `equality` |
| `--result_file` | - | Output CSV file path for search results |
| `--base_path` | - | Path to base vector file (`.fvecs`) |
| `--query_path` | - | Path to query vector file (`.fvecs`) |
| `--base_label_path` | - | Path to base label file (`.csv`) |
| `--query_label_path` | - | Path to query label file (`.csv`) |
| `--gt_path` | - | Path to ground truth file (`.bin`) |

## Output Format

**Search mode:**

```
gt_path, query_type, ef_search, QPS, recall
```

**Insert mode:**

```
gt_path, query_type, build_percent, insert_count, insert_time, insert_qps, ef_search, QPS, recall
```
