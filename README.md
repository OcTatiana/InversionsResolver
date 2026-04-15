# InversionResolver

## Short Description

**InversionResolver** is a command-line tool for detecting and resolving inversions in genomic data.  
It supports two modes of operation:

- **`sort`** — runs the inversion resolving algorithm on preprocessed data  
- **`genome`** — runs a full pipeline including preprocessing and inversion resolution

---

## Installation

Install directly from GitHub:

```bash
pip install git+https://github.com/OcTatiana/InversionsResolver.git
```

After installation, the command-line tool becomes available:

```bash
InversionResolver --help
```

---

## Commands Overview

### 1. `sort` — Algorithm Only

Runs the inversion resolving algorithm on an already prepared input file with signed permutation.

**Input:**

* File with space-separated signed permutation

**Output:**

* Permutation sequence (.perm file)
* Changed IDs (.id file)
* Visualization (.png file)


---

### 2. `genome` — Full Pipeline

Runs the complete workflow:

```
Raw data → Filtering → Chaining → Inversion resolving → Output
```

**Input:**

* Parsed psl file (tab-separated) with synteny blocks coordinates for query and target

**Output:**

* Log files with filtered synteny 
* Filtered .xlsx  and .csv files with chains' strand and coordinates
* Number-coded permutation for all chromosomes
* Resolved permutation sequence and IDs
* Visualization for every chromosome in 2 modes: unscaled and scaled bt syntheny chain length

---

## Quick Start (Test Data)

Example using test data from the repository:

### Download data
```bash 
wget https://raw.githubusercontent.com/OcTatiana/InversionsResolver/refs/heads/main/test/test_permutation.txt
```

### Run algorithm only

```bash
InversionResolver sort \
  -i test_permutation.txt \
  -o single_perm_output 
```

---

### Run full pipeline

TBA

---

## Parameters Description

### `sort` command

| Option              | Description                     |
| ------------------- | ------------------------------- |
| `-i`, `--input`     | Path to input file (required)   |
| `-o`, `--output`    | Output file name (required)     |
| `-v`, `--visualize` | Enable visualization (optional) |
| `-s`, `--seed`      | Random seed (optional)          |

---

### `genome` command

#### Required parameters

| Option               | Description                        |
| -------------------- | ---------------------------------- |
| `-i`, `--input`      | Path to raw input file             |
| `-o`, `--output-dir` | Output directory                   |
| `-sp`, `--species`   | Species name (used in file naming) |

---

#### Optional parameters

| Option                     | Default | Description                |
| -------------------------- |---------| -------------------------- |
| `-v`, `--visualize`        | False   | Enable visualization       |
| `-s`, `--seed`             | 30      | Random seed                |
| `-ol`, `--overlap-length`  | 100000  | Minimum overlap length     |
| `-op`, `--overlap-percent` | 0.1     | Minimum overlap percentage |
| `-ro`, `--remove-overlaps` | False   | Remove all overlaps        |
| `-cs`, `--chain-step`      | 1000000 | Chaining step size         |
| `-nc`, `--no-chaining`     | False   | Disable chaining           |
