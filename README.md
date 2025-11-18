# Brain Transcriptomic Drug Response Predictor for ADHD

Genotype-based transcriptomic imputation and drug response prediction system for personalized ADHD treatment optimization.

## Overview

This repository contains a computational framework that predicts individual-specific transcriptomic profiles across brain tissues from genotype data, then uses these profiles to forecast response to psychiatric medications and prioritize the most effective drugs for ADHD treatment.

## Key Features

- **Multi-Tissue Transcriptome Imputation**: Predicts gene expression across 10 brain regions from genotype data
- **Drug Response Prediction**: Models individual response to stimulant and non-stimulant ADHD medications
- **Personalized Drug Prioritization**: Ranks medications by predicted effectiveness for each individual
- **Polygenic Score Integration**: Incorporates cognitive and ADHD genetic risk profiles
- **Clinical Validation**: Tested on SPARK and ABCD cohorts

## Methodology

### Transcriptomic Imputation
- **Input**: Genotype data (SNP arrays or WGS)
- **Model**: Tissue-specific prediction models using GTEx reference
- **Output**: Predicted expression levels for 12,000+ genes across brain regions

### Drug Response Prediction
- **Mechanism**: Maps transcriptomic profiles to drug perturbation signatures (LINCS L1000)
- **Drugs Covered**: Methylphenidate, Amphetamines, Atomoxetine, Guanfacine, Clonidine
- **Response Metric**: Continuous score from -1:1 indicating predicted effectiveness

### Genetic Risk Integration
- **PGS Calculation**: ADHD polygenic risk scores from latest GWAS
- **Cognitive gFactor**: General cognitive ability polygenic scores
- **Interaction Modeling**: How genetic risk modifies drug response

### Ethical Considerations

**Clinical Use**: Research tool only - not for direct clinical decision-making

**Genetic Privacy**: Local processing option to avoid data transfer

**Bias Mitigation**: Population-specific models for diverse ancestry groups

**Interpretability**: Transparent feature importance and confidence scores
