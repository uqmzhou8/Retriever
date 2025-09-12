# Tutorial: Reference-Free Imputation with Retriever and Beagle

*Based on Zhou et al. (2025) "Chimeric Reference Panels for Genomic Imputation"*  
*Developed by the Ortiz-Barrientos Lab: www.ortizbarrientoslab.org*

This tutorial walks you through using Retriever, a method that creates chimeric reference panels from your own samples, combined with Beagle for genotype imputation when external reference panels are unavailable.

## Table of Contents

1. [Overview and Background](#1-overview-and-background)
2. [Installation](#2-installation)
3. [Input File Requirements](#3-input-file-requirements)
4. [Parameter Selection](#4-parameter-selection)
5. [Running Retriever](#5-running-retriever)
6. [Beagle Imputation](#6-beagle-imputation)
7. [Troubleshooting](#7-troubleshooting)
8. [Performance Expectations](#8-performance-expectations)
9. [Advanced Usage](#9-advanced-usage)

## 1. Overview and Background

### The Reference Panel Problem

Genotype imputation has been useful for genomics research, allowing studies to expand from thousands of directly genotyped variants to millions of imputed variants. However, conventional imputation methods require external reference panels with hundreds of fully genotyped individuals, creating practical barriers for most non-model organisms.

### How Retriever Works

Retriever constructs a chimeric reference panel directly from target samples using a sliding-window approach to identify and retrieve genomic partitions with complete data. By exploiting the complementary distribution of missing data across samples, Retriever assembles a panel that preserves local patterns of linkage disequilibrium and captures novel variants.

The approach assumes that while individual samples may have missing genotypes, the collective dataset often contains complete information at each position across different subsets of individuals.

### Biological Foundation

Imputation relies on linkage disequilibrium (LD) - the correlation between nearby genetic variants due to limited recombination. This correlation follows the relationship:

P(linked after t generations) = (1 - r)^t

Where r ≈ d × 10^-8 for distance d in base pairs. LD patterns vary between populations, with African populations typically showing faster decay and inbred lines showing extensive blocks.

## 2. Installation

Retriever makes use of external packages, so installation in a virtual environment is recommended:

```bash
# Method 1: Using provided environment file
conda env create -f Retriever.yml
conda activate Retriever

# Method 2: Manual installation
conda create -n Retriever
conda activate Retriever
conda install numpy pandas cyvcf2 scikit-learn
```

These packages are prerequisites for running Retriever and are available in both Anaconda and pip.

### Testing Installation

Test your installation using the provided sample script:

```bash
python Sample_script.py
```

This should finish execution within 2 seconds using the dummy dataset. Ensure both `dummy_dataset.vcf` and `Sample_script.py` are in the same folder.

### Download Beagle

Download Beagle4.1 (recommended over Beagle5 for small reference panels):

```bash
wget https://faculty.washington.edu/browning/beagle/beagle.27Jan18.7e1.jar
```

**Note:** Validation studies found that Beagle5 performs poorly with small reference panels (<100 individuals), while Beagle4.1 maintains stable performance with panels of 25-50 individuals.

## 3. Input File Requirements

### File Format

Retriever is designed for VCF version 4.0 onwards with all samples aligned to a single reference genome. Both compressed (`file_name.vcf.gz`) and uncompressed (`file_name.vcf`) files are compatible.

### Required Preprocessing

**Population Structure:**
- Split files so each contains only a single population
- This reduces population structure bias and improves accuracy

**Chromosome Organization:**
- Each chromosome should be in its own VCF file
- Enables parallel processing during chimeric reference panel assembly

**Variant Filtering:**
- Retain only biallelic variants
- Exclude variants with more than 20% missing genotypes
- Remove duplicated positions
- Ensure positions are sorted in ascending order
- Use diploid individuals only

### Filtering with bcftools

Example command to extract one chromosome, remove multiallelic sites, and retain sites with less than 20% missing data:

```bash
bcftools view -m2 -M2 -v snps \
    --regions 'chr1' \
    input_samples.vcf.gz \
    -i 'F_MISSING<0.2' \
    -o chr1_filtered.vcf
```

Replace `chr1`, `input_samples.vcf.gz`, and `chr1_filtered.vcf` with your specific file names.

### Post-Imputation File Merging

After imputation, merge chromosomes back into a single file if needed:

```bash
bcftools merge chr1_imputed.vcf.gz chr2_imputed.vcf.gz -o merged_files.vcf.gz
```

### Sample Size Considerations

While imputation can be performed with approximately 200 samples, larger sample sizes generally lead to higher-quality reference panels and more accurate imputation results. Validation testing used datasets ranging from 165-1135 individuals.

## 4. Parameter Selection

### Evidence-Based Recommendations

Optimization testing across different masking levels (1-30% missing data) and multiple species found:

**Optimal Panel Size: 25-50 individuals**
- Approximately 25% of total sample size or 50 individuals
- Performance peaked around 25-50 individuals before declining
- Larger panels (>75) sometimes reduced accuracy due to small partition effects

**Window Size: 1000 bp (default)**
- Tested range: 100-10,000 bp  
- Window size had minimal impact on accuracy
- Default setting performed adequately across genomic contexts

**Recombination Rate Sensitivity:**
- Tested 0.5-2.0 cM/Mb (Beagle default: 1.0 cM/Mb)
- Changes had negligible effect on accuracy
- Advantage for organisms lacking detailed genetic maps

## 5. Running Retriever

### Import and Basic Usage

Prior to using Retriever, import the package:

```python
from Retriever import *
```

### Assembly of Chimeric Reference Panel

```python
vcf2ref('file_name.vcf', num_ref, window_size=1000, outfile='chimeric_ref_gts.vcf')
```

#### Required Arguments

**file_name:** Name of the input VCF file containing genomic data for chimeric reference genome creation. Must be in VCF format.

**num_ref:** Number of chimeric reference genomes to create. Recommended: approximately 25% of total sample size or 50 individuals.

#### Optional Arguments

**window_size:** Window length in base pairs for chimeric reference genome assembly. Larger values result in larger complete partitions but may increase computational time. Default = 1000

**outfile:** Output file name for chimeric reference genomic data. Default = 'chimeric_ref_gts.vcf'

### Single Chromosome Example

```python
from Retriever import *

# Basic implementation with validated parameters
vcf2ref('chr22_filtered.vcf', 
        num_ref=50, 
        window_size=1000, 
        outfile='chr22_chimeric_ref.vcf')

print("Chimeric reference panel assembly complete")
```

### Parallelization Across Chromosomes

The multiprocessing tool allows several chromosomes to be processed concurrently:

```python
from Retriever import *
from multiprocessing import Pool

# Define chromosome names
chromosome_numbers = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5']

def par_Ret(chr_name):
    vcf2ref(f'filtered_{chr_name}.vcf', 
            num_ref=50, 
            window_size=1000, 
            outfile=f'chimeric_ref_{chr_name}.vcf')

# Run in parallel
with Pool(4) as p:  # Adjust number based on available CPU cores
    p.map(par_Ret, chromosome_numbers)
```

For more information on multiprocessing, visit: https://docs.python.org/3/library/multiprocessing.html

## 6. Beagle Imputation

### Single Chromosome Imputation

```bash
java -Xmx8g -jar beagle.27Jan18.7e1.jar \
    gt=chr22_filtered.vcf \
    ref=chr22_chimeric_ref.vcf \
    out=chr22_imputed \
    nthreads=4
```

### Batch Processing Multiple Chromosomes

```bash
for chr in {1..22}; do
    echo "Processing chromosome ${chr}..."
    java -Xmx8g -jar beagle.27Jan18.7e1.jar \
        gt=chr${chr}_filtered.vcf \
        ref=chr${chr}_chimeric_ref.vcf \
        out=chr${chr}_imputed \
        nthreads=4
done
```



## 7. Troubleshooting

### Common Errors and Solutions

**"Not enough individuals in window" Error:**
- **Cause:** Window contains fewer complete genotypes than required panel size
- **Solutions:** 
  - Reduce `num_ref` parameter (try 25 instead of 50)
  - Increase `window_size` to 2000-5000 bp
  - Check input data preprocessing

**Memory Issues:**
- **Cause:** Buckets stored in local memory exceed available RAM
- **Solutions:**
  - Allocate at least 3 times expected data size in memory
  - Process chromosomes separately
  - Reduce panel size if necessary

**Low Imputation Accuracy (<90%):**
- **Causes:** Population structure, insufficient sample size, poor preprocessing
- **Solutions:**
  - Verify samples are from single population
  - Check for multiallelic sites or excessive missingness
  - Ensure adequate sample size (≥200 individuals recommended)

### Special Considerations

- Ensure only 1 chromosome per VCF file input to Retriever
- Only 1 CPU core required for computing; additional cores do not expedite the program
- CPU availability allows parallelization across multiple chromosome files
- If insufficient individuals in any window, program will terminate with error
- Position must include at least as many genotypes as chimeric panel size

## 8. Performance Expectations

### Computational Requirements

Based on scaling analysis with human chromosome 1:

**Runtime Expectations:**
- 1% missing data: ~2-4 hours
- 10% missing data: ~6-8 hours  
- 20% missing data: ~15-20 hours
- 30% missing data: ~25-30 hours

**Memory Requirements:**
- Allocate 3x expected data size in available memory
- Temporary data structures stored during assembly

**CPU Usage:**
- Retriever uses single CPU core by design (preserves LD patterns)
- Beagle can utilize multiple threads (tested up to 96 threads)
- Parallel processing recommended across chromosomes

### Dataset Size Effects

Validation showed performance scales with:
- **Sample Size:** Minimum ~200 individuals, accuracy improves with larger samples
- **Missing Data:** Method tolerates up to 30% overall missing data
- **Population Structure:** Single populations perform better than mixed populations

## 9. Advanced Usage

### Novel Variant Handling

Testing showed that when genomic positions were removed from external reference panels:
- 10% positions removed: External panel accuracy dropped to ~85%
- 20% positions removed: External panel accuracy dropped to ~75%  
- Retriever maintained >98% accuracy regardless

This suggests utility for datasets containing population-specific variants not in existing reference resources.

### Comparison with Other Methods

Validation benchmarking against machine learning approaches:
- **Retriever + Beagle4:** >95% accuracy consistently
- **missForest:** 93-97% accuracy (best reference-free method)
- **k-Nearest Neighbors:** 90-94% accuracy

### Case Study: Downstream Analysis Effects

Analysis of *Helianthus annuus* showed that imputation affected recombination rate estimates. Missing data appeared to artificially inflate estimates, while imputed data produced rates more consistent with biological expectations.

### Integration Considerations

- Imputed data compatible with standard population genetic analysis tools
- Consider how imputation might influence specific downstream analyses
- Quality filtering based on R² values recommended for critical applications

### Future Directions

The method provides a foundation for extending imputation to non-model organisms. Areas for development include adaptation for polyploid species, integration with long-read sequencing data, and optimization for ancient DNA applications.

---

For detailed methodology and validation results, see Zhou et al. (2025) "Chimeric Reference Panels for Genomic Imputation" and visit the Retriever GitHub repository (https://github.com/uqmzhou8/Retriever) and the Ortiz-Barrientos Lab website (www.ortizbarrientoslab.org).
