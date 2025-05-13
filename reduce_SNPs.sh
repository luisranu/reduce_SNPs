#!/bin/bash


show_help() {
    echo "Usage: bash $0 -I <input_file> -O <output_file> -d <distance> -f <frequency>"
    echo
    echo "Options:"
    echo "  -I    Input file in .vcf or .vcf.gz format (required)"
    echo "  -O    Output file to be generated (required)"
    echo "  -d    Distance in base pairs (integer, default: 0)"
    echo "  -f    Minimum allele frequency difference (decimal, default: 0.00)"
    echo "  -h    Show this help message"
    echo
    echo "Example:"
    echo "    bash $0 -I input.vcf.gz -O output.vcf.gz -d 1000 -f 0.05"
    exit 1
}



# Show help if "help" is passed as the first argument
if [[ "$1" == "help" || "$1" == "--help" || "$1" == "-h" ]]
then
    show_help
fi

# Default values
nb=0
daf=0.00

# Read options
while getopts ":I:O:d:f:h" opt; do
  case $opt in
    I) input_file="$OPTARG" ;;
    O) output_file="$OPTARG" ;;
    d) nb="$OPTARG" ;;
    f) daf="$OPTARG" ;;
    h) show_help ;;
    \?) echo "Invalid option: -$OPTARG" >&2; show_help ;;
    :) echo "Option -$OPTARG requires an argument." >&2; show_help ;;
  esac
done

# Validate that input and output files are provided
if [[ -z "$input_file" || -z "$output_file" ]]; then
  echo "❌ Error: You must specify an input file (-I) and an output file (-O)"
  show_help
fi

# Validate that input file exists and is readable
if [[ ! -f "$input_file" || ! -r "$input_file" ]]; then
  echo "❌ Error: Input file '$input_file' does not exist or is not readable"
  exit 1
fi

# Validate that nb is a positive integer
if ! [[ "$nb" =~ ^[0-9]+$ ]]; then
  echo "❌ Error: The value of -d (distance) must be a positive integer"
  exit 1
fi

# Validate that daf is a valid decimal between 0 and 1 using a dot
if ! [[ "$daf" =~ ^0\.[0-9]+$|^1\.0+$ ]]; then
  echo "❌ Error: The value of -f (frequency) must be a decimal number between 0 and 1 (use dot as decimal separator)"
  exit 1
fi

# Show parameters for confirmation
echo "✅ Input file         : $input_file"
echo "✅ Output file        : $output_file"
echo "✅ Maximum distance   : $nb"
echo "✅ Minimum difference : $daf"


# Make temporal directory
# -----------------------

dir=temp_reduce_SNPs
mkdir -p $dir




# Filtering input 
# ---------------
# -m2 -M2: Filters to include only variants with exactly 2 alleles.
# -v snps: Filters to include only SNPs.
bcftools view -m2 -M2 -v snps $input_file -Ou | \
# get Information of CHR, POS, FREC y GT
bcftools query -f '%CHROM\t%POS\tAF=%AF[\t%GT]\n' | sed 's/AF=//g' > $dir/temp.GT.info.tsv

# Get header and add line with info about the run
# -----------------------------------------------
bcftools view -h $input_file > $dir/temp.preheader.txt 
current_date=$(date +"%Y-%m-%d %H:%M:%S")
cat $dir/temp.preheader.txt | grep "^##" > $dir/temp.header.txt 
echo "##reduced_SNPs=reduced_SNPs -I $input_file -O $output_file -d $nb -f $daf; Date= $current_date" >> $dir/temp.header.txt
tail -n 1 $dir/temp.preheader.txt >> $dir/temp.header.txt 

# Functions
# ---------

# Refreshes the variables

refresh() {
    old_cr=$cr
    old_pos2=$old_pos1
    old_pos1=$pos
    old_af2=$old_af1
    old_af1=$af
    old_gt2=$old_gt1
    old_gt1=$gt
}

# Calculate the differences in allele frequencies
dif_freqs() {
    local a="$1"
    local b="$2"
    if [[ -n "$a" && -n "$b" ]]; then
        awk -v d1="$a" -v d2="$b" 'BEGIN {print (d1>d2) ? d1-d2 : d2-d1 }'
    else
        echo ""  # Do not calculate if a value is missing
    fi
}


# SCRIPT 

# If daf=0.00, then comparations are senseless. So nb=0
if [ "$(echo "$daf == 0" | bc)" -eq 1 ]
then
    nb=0
fi



while read -r line
do 
    # Get variables
    cr=$(echo $line | cut -d ' ' -f1)
    pos=$(echo $line | cut -d ' ' -f2)
    af=$(echo $line | cut -d ' ' -f3)
    gt=$(echo $line | cut -d ' ' -f4-)
    # If SNP changes of chromosome, write and refresh
    if [[ $cr != $old_cr ]]                                                                                                          # Cond1
    then 
        echo -e "$cr\t$pos" >> $dir/temp.reduced.snps
        refresh
    # If SNP do not change of chomosome
    else      
    # If the genotypes are different from the previous 5                                  
        if [[ "$gt" != "$old_gt1" && "$gt" != "$old_gt2" ]]                                                                          # Cond 2
        then
    # If nb is equal to zero                                                                                                         # Cond 3
            if [ $nb -eq 0 ]
            then
                echo -e "$cr\t$pos" >> $dir/temp.reduced.snps
                refresh
            else 
    # If the new position exceeds the limit of base pairs relative to the previous position
                res=$(echo "$old_pos1 + $nb" | bc )
                if (( res < pos ))                                                                                                   # Cond 5
                then
                    echo -e "$cr\t$pos" >> $dir/temp.reduced.snps
                    refresh
    # If the new position do not exceeds the limit of base pairs relative to the previous position
                else
    # Get differences in alelic frequencies 
                    dif1=$(dif_freqs "$old_af1" "$af")
                    dif2=$(dif_freqs "$old_af2" "$af")
    # Sort and keep the smallest difference                    
                    dif_min=$(printf "%s\n" "$dif1" "$dif2" | grep -v '^$' | sort -n | head -n 1)
    # If the difference is greater than the one selected by the user
                    if (( $(awk -v d1="$dif_min" -v daf="$daf" 'BEGIN {print (d1 > daf) }' ) ))                                      # Cond 5
                    then 
                        echo -e "$cr\t$pos" >> $dir/temp.reduced.snps
                        refresh
                    fi                                                                                                               # FIN Cond 5
                fi                                                                                                                   # FIN Cond 4
            fi                                                                                                                       # FIN Cond 3
        fi                                                                                                                           # FIN Cond 2
    fi                                                                                                                               # FIN Cond 1
done < $dir/temp.GT.info.tsv
   
orig=$(bcftools view -H $input_file | wc -l)   
binary=$(cat $dir/temp.GT.info.tsv | wc -l)
reduc=$(cat $dir/temp.reduced.snps | wc -l)
perc=$(awk -v b=$binary -v r=$reduc 'BEGIN { print  (100*r/b) }')
   
# Select from original VCF the positions obtained by reduction and join with header

bcftools view -H -R $dir/temp.reduced.snps -O z $input_file  > $dir/temp.corpus.vcf.gz

cat $dir/temp.header.txt | gzip -c > $dir/temp.header.gz

zcat $dir/temp.header.gz $dir/temp.corpus.vcf.gz | gzip -c > $output_file


echo "The file $input_file contains $orig variants. Of these, $binary are SNPs with only one ALT allele."  
echo "After running reduce_SNPs.sh with a distance threshold of $nb nucleotides and a frequency threshold of $daf,"
echo "$output_file is obtained with $reduc SNPs — which corresponds to $perc % of the biallelic SNPs."

rm -r $dir
   
   
    
