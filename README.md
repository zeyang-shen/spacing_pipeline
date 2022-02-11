[![python-version](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![DOI](https://zenodo.org/badge/DOI/10.1101/2020.04.02.021535.svg)](https://www.biorxiv.org/content/10.1101/2020.04.02.021535v1.full)

# Spacing pipeline
The scripts provided in this repository are used to compute and characterize the spacing relationships of transcription factors. 

Here is the overview of the method:

<p align="center">
<img src="https://github.com/zeyang-shen/spacing_pipeline/blob/main/ENCODE_processing_pipeline.png" width="850" height="163">
</p>

## Dependencies
* Python ([https://www.python.org/downloads/](https://www.python.org/downloads/))
* NumPy 1.15.4 ([https://numpy.org/install/](https://numpy.org/install/))
* pandas 1.1.4 ([https://pandas.pydata.org/](https://pandas.pydata.org/))
* Biopython 1.70 ([https://biopython.org/wiki/Download](https://biopython.org/wiki/Download))
* SciPy 1.1.0 ([https://www.scipy.org](https://www.scipy.org))
* Matplotlib 3.3.3 ([https://matplotlib.org/](https://matplotlib.org/))
* Seaborn 0.11.0 ([https://seaborn.pydata.org/installing.html](https://seaborn.pydata.org/installing.html))
* HOMER ([http://homer.ucsd.edu/homer/download.html](http://homer.ucsd.edu/homer/download.html))

## Quick Usage
[identify_motif.py](https://github.com/zeyang-shen/spacing_pipeline/blob/main/scripts/identify_motif.py) can find motifs given a peak file, a FASTA file for peak sequences, and a motif file. The recommended parameters are as below to filter for motifs passing a false positive rate <0.1% (--cutoff) and a location <50 bp from peak centers (-d 50): 
```bash
python identify_motif.py ../ENCODE_processed_files/CTCF_idr.fa CTCF --motif_path ../motifs/ --cutoff -d 50
```

To identify motifs and simultaneously separate peaks into those falling at **repetitive** and nonrepetitive DNA regions, please download the repeats annotations first and run `identify_motif.py` script by specifying `--repeat`:

```bash
wget https://homer.ucsd.edu/zeyang/hg38_repeats.tar.gz
tar -zxvf hg38_repeats.tar.gz
python identify_motif.py ../ENCODE_processed_files/CTCF_idr.fa CTCF --motif_path ../motifs/ --cutoff -d 50 --repeat hg38_repeats/hg38_repeats_merged.nodup.all.txt
```

[characterize_spacing.py](https://github.com/zeyang-shen/spacing_pipeline/blob/main/scripts/characterize_spacing.py) can take in two processed files from `identify_motif.py` for a pair of transcription factors and output results of spacing relationships. The basic usage is as below:
```bash
python characterize_spacing.py ../ENCODE_processed_files/ GATA1 TAL1 --motif_path ../motifs/
```

## Citation
If you use our findings or scripts, please cite our paper: [https://doi.org/10.7554/eLife.70878](https://doi.org/10.7554/eLife.70878).

## Data
`motifs/` folder stores the PWM files in the JASPAR format used in the paper.

`ENCODE_processed_files/` folder includes the processed data of this paper based on ENCODE ChIP-seq data:
* _idr.tsv -- ChIP-seq peaks in HOMER peak file format after running IDR
* _idr.fa -- sequences of ChIP-seq peaks in _idr.tsv
* _idr_cutoff.tsv -- ChIP-seq peaks that have been identified to have valid motifs
* _idr_cutoff_inmask.tsv -- Peaks in _idr_cutoff.tsv that fall into repetitive regions
* _idr_cutoff_masked.tsv -- Peaks in _idr_cutoff.tsv that fall into nonrepetitive regions

## Contact
If you enconter a problem when using the scripts, you can
1. post an issue on [Issue](https://github.com/zeyang-shen/spacing_pipeline/issues) section
2. or email Zeyang Shen by zes017@ucsd.edu

## License
[This project is licensed under GNU GPL v3](https://github.com/zeyang-shen/spacing_pipeline/blob/main/LICENSE)

## Contributors
The scripts were developed primarily by Zeyang Shen and Rick Zhenzhi Li. Supervision for the project was provided by Christopher K. Glass. 
