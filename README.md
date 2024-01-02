# Integration-Sites-IFO0880

### THE PAPER!
This repository accompanies the work ["XYZ"](https://www.google.com).

Add ifo0880 folders to code and data folder of ["CRISPR-COPIES"] (https://github.com/Zhao-Group/COPIES) repository.

First, run the script to incorporate putative essential genes of *Rhodotorula toruloides*:
```
python code/ifo0880/deg_update.py
```

Add the genome, feature table, and protein sequences of *Rhodotorula toruloides* strain NBRC 0880 from NCBI to data/ifo0880 folder. Utilize the updated database of essential genes to run the CRISPR-COPIES command line tool and obtain genome-wide integration sites- 
```
python code/main_rt.py -g ../data/ifo0880/GCA_000988875.2_ASM98887v2_genomic.fna -t ../data/ifo0880/GCA_000988875.2_ASM98887v2_feature_table.txt -p NGG -o 3prime -l 20 --edit_dist 5 --intspace 400 --distal_end_len 20000 -hr_l 1000 --protein_file ../data/ifo0880/GCA_000988875.2_ASM98887v2_protein.faa --blast_org 'Rhodosporidium toruloides IFO0880' --GC_grna 20,80 --polyG_grna 6 --polyT_grna 6 --RE_grna CCTGCAGG,GGCGCGCC,GTTTAAAC --RE_hr CCTGCAGG,GGCGCGCC,GTTTAAAC -sl 10 -out ../data/ifo0880/output_10_5.csv
```

Next, run the scripts on the generated output file to integrate transcriptomics information:
```
python code/ifo0880/transcript.py
```

To reproduce the results for zero-shot predictions of on-target efficiency for gRNAs characterized in this study (Supplementary Figure S1), run:
```
python code/ifo0880/on-target.py
python code/Uni-deepSG/PredictionTools/predict.py (Note: Requires CUDA)
```

### Reference
<details>
<summary>If you using these scripts, please cite us:</summary>

```bibtex

```
</details>
