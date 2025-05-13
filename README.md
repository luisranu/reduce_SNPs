# 🧬 reduce_SNPs.sh

A simple Bash script to filter biallelic SNPs from a VCF file by removing nearby variants that carry redundant information. Optionally, the user can define a distance (in base pairs) and a minimum 

## 📄 Description

The script performs the following steps:

1. **SNP Selection**: Retains only biallelic SNPs (those with a single ALT allele).
2. **Redundancy Filtering**: Compares the genotype (GT) string of each SNP with the three previously selected SNPs. If an exact match is found, the SNP is excluded.
3. **Distance & Frequency Filtering**:
   - If the distance to the last selected SNP is **greater than** the specified threshold (`-d`), the SNP is included.
   - If the distance is **less than or equal to** the threshold, the SNP is included **only if** its allele frequency differs from each of the three previous SNPs by at least the specified frequency threshold (`-f`).

> **Note**: By default, both `-d` and `-f` are set to `0`, meaning all SNPs are compared regardless of distance, and **any** allele frequency difference is considered.


## 💻 Implementation

This script is written in **Bash**, designed to work efficiently in Unix-based environments.

## ⚙️ Requirements

- **Operating System**: Linux or macOS
  (Windows users can run it via WSL or a Unix-like shell environment)
- **Bash shell**: Bash (tested with version 5+)
- **Dependencies**:
  - [`bcftools`](https://samtools.github.io/bcftools/): Must be installed and available in your system's `PATH`.


## 📥 Input

The input file must be a compressed VCF (`.vcf.gz`).

### Parameters:

* `-I` : Input VCF file (`.vcf.gz`, required)
* `-O` : Output VCF file (`.vcf.gz`, required)
* `-d` : Distance (in bp) to keep between SNPs (integer, default: `0`)
* `-f` : Minimum allele frequency difference (decimal, default: `0.00`)
* `-h` : Help message


## 🚀 Usage

```bash
bash reduce_SNPs.sh -i input.vcf.gz -o output.vcf.gz [-d distance] [-f frequency]
```

### 🧪 Example:

```bash
bash reduce_SNPs.sh -i input.vcf.gz -o filtered.vcf.gz -d 1000 -f 0.05
```

### Example logic:
With a distance threshold of 1000 nucleotides (`-d 1000`) and a frequency threshold of 0.05 (`-f 0.05`), a new **non redundant** SNP will be included if:
- It is **more than 1000 nucleotides** away from the last selected SNP, **or**
- It is within 1000 nucleotides but its allele frequency differs by at least 0.05 from **each** of the last three selected SNPs.

## 🎯 Why Reduce SNPs?

Working with a smaller, cleaner set of SNPs has many practical advantages. Especially when you're dealing with large-scale genetic data.

### ✅ Advantages of Reducing the SNP Set

- **Lower Computational Cost**: Fewer SNPs reduce memory usage and runtime, especially in large-scale GWAS or population genomics.
- **Improved Statistical Power Through Reduced Redundancy**: In genome-wide studies, nearby SNPs are often strongly correlated due to linkage disequilibrium (LD), meaning they carry overlapping information. By retaining only representative SNPs, the dataset becomes less redundant, resulting in cleaner inputs and fewer statistical tests. This, in turn, improves the power to detect true associations after multiple testing correction (e.g., Bonferroni, FDR).
- **More Reliable Models**: Too many similar SNPs can confuse statistical models and make results unstable. Removing correlated SNPs helps prevent multicollinearity—a condition where predictors (SNPs) provide overlapping information—leading to more stable coefficient estimates, smaller standard errors, and better detection of true associations in regression or mixed models.

### 📌 When to Filter Nearby SNPs with Similar Allele Frequencies

Filtering out close-together SNPs with very similar allele frequencies is especially useful when:

- **You're Working with Low-Coverage Data**: In low-coverage sequencing, genotyping errors are more frequent due to insufficient read depth. As a result, small fluctuations in read counts can lead to inconsistent genotype calls and artificially inflate differences in allele frequencies, making nearby non-informative SNPs appear informative. Filtering out such SNPs with similar frequencies helps eliminate redundant variants or those in high linkage disequilibrium, improving the robustness and interpretability of downstream analyses.
- **You’re Working with Poor-Quality Samples**: With ancient DNA, environmental samples, or degraded material, filtering can help minimize false positives.
- **You're Running Pilot Analyses**: A reduced SNP set makes exploratory analyses faster and easier.

In contrast, if you're performing fine-mapping, haplotype reconstruction, or focusing on regions where every variant matters, it's best **not** to filter out nearby SNPs—even if their frequencies are similar.


