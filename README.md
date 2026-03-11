# Monte Carlo mutation rate analysis for missense variants in exons

Requires:
  - BED file containing exons (https://github.com/MJLW/extract_exons_from_gff)
  - DNMs in TSV format (CHROM\tPOS\tREF\tALT)
  - VCF/TSV containing mutation rates per base

Use extract_consequence.sh to convert VCF containing an INFO/MR (mutation rate) to a TSV (CHROM\tPOS\tREF\tALT\tMR):
```
bash extract_consequence.sh <input.vcf> <output.tsv>
```

Calculate p-values:
```
python3 run.py --exons <exons.bed> --dnms <dnms.tsv> --roulette <mr.tsv> --out <out.tsv>
```
