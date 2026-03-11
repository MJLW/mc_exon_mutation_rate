bcftools +split-vep -i 'Consequence="missense_variant"' -f "chr%CHROM\t%POS\t%REF\t%ALT\t%INFO/MR" $1 > $2
