#!/bin/bash

#############################################
# Bioinformatics
# 2017
# Exercise 1: Sequence Processing
#############################################


#The start codon has the sequence "AUG", and the stop codon has the sequence "UAG", "UAA", or "UGA"
#In DNA is TAC, whose complementary is ATG
#Stop codons in DNA: TAA, TGA, TAG
#

function codon_to_aminoacid ()
{
case $codon_to_string in
     	GCT|GCC|GCA|GCG)
        echo "A"
        ;;
    	CGT|CGC|CGA|CGG|AGA|AGG)
        echo "R"
        ;;
    	AAT|AAC)
        echo "N"
	;;
	GAT|GAC)
        echo "D"
	;;
	TGT|TGC)
	echo "C"
	;;
        CAA|CAG)
        echo "Q"
	;;
        GAA|GAG)
        echo "E"
	;;
	GGT|GGC|GGA|GGG)
        echo "G"
	;;
        CAT|CAC)
        echo "H"
	;;
	ATT|ATC|ATA)
        echo "I"
	;;
        ATG)
        echo "M"
	;;
  	TTA|TTG|CTT|CTC|CTA|CTG)
        echo "L"
	;;
        AAA|AAG)
        echo "K"
	;;
  	TTT|TTC)
        echo "F"
	;;
        CCT|CCC|CCA|CCG)
        echo "P"
	;;
	TCT|TCC|TCA|TCG|AGT|AGC)
        echo "S"
	;;
        ACT|ACC|ACA|ACG)
        echo "T"
	;;
        TGG)
        echo "W"
	;;
        TAT|TAC)
        echo "Y"
	;;
        GTT|GTC|GTA|GTG)
        echo "V"
	;;
	TAA|TGA|TAG)
        echo "Z"


        #exit 1
        ;;
    *)
        echo "Unknown_Codon"
esac
}


function identify_codon () 
{
unset aminoacids_array[@]

counter=${#gene_aux[@]}
while [ $counter -ge 3 ]
    do
		codon=(${gene_aux[0]} ${gene_aux[1]} ${gene_aux[2]}) 
		codon_to_string="${codon[0]}${codon[1]}${codon[2]}"
		aa=$(codon_to_aminoacid $codon_to_string)
		aminoacids_array+=("$aa")
		cant_elements=${#gene_aux[@]}
                gene_aux=(${gene_aux[@]:3:$cant_elements})
            counter=${#gene_aux[@]}
	done


}

function dna_to_complementary () 
{
gene_aux=(${gene[@]})    
for ((i=0;i < ${#gene_aux[@]}; i++))

do
    aux=${gene_aux[i]}
        case $aux in
        A)
        gene_aux[i]=T
        ;;
        T)
        gene_aux[i]=A
        ;;
        C)
        gene_aux[i]=G
        ;;
        G)
        gene_aux[i]=C
        ;;
    *)
        echo "Unknown"
esac
done
}

convert_string_into_array () {

uppercase=$(echo $1 | tr [:lower:] [:upper:])

unset result[@]
for (( x=0; x<${#1}; x++))
do
        result[$x]=${uppercase:x:1}
done
}

function format_output_fasta () {

function join { local IFS="$1"; shift; echo "$*"; }
result_aa=$(join "" ${aminoacids_array[@]})

}

function frame ()
{

aux=$1

gene_aux=(${gene[@]})

#Deletes aux amount of nucleotid (for reading frames)
gene_aux=("${gene_aux[@]:$aux}")
identify_codon
format_output_fasta

}

function frameR ()
{

aux=$1

gene_aux=(${reversed_array[@]})

#Deletes aux amount of nucleotid (for reading frames)
gene_aux=("${gene_aux[@]:$aux}")
identify_codon
format_output_fasta

}

function string_gene_description () {

orf=$1
direction=$2
locus=$(awk '/LOCUS/ {print $2}' $genbak_file)
definition=$(awk '/DEFINITION/ {$1=""; print $0}' $genbak_file)

echo -n ">"$locus"-"$orf$2
echo "" $definition
}

###############################################################
###############################################################
# Body
###############################################################

# Genbank file is received as an argument
genbak_file=$1

# Delete temp files generated thourgh previous runs
rm -f clean_gb_file
rm -f input_file

### Text processing & formatting
#########################################################################################
# Takes the nucleotids sequence from the genbank file by using the word ORIGIN as a pattern matching for regular exp (awk)
awk '$0 ~ /ORIGIN/ { vart = NR }{ arr[NR]=$0 } END { for (i = vart; i<=NR ; i++) print arr[i]  }' $genbak_file > clean_gb_file
awk '{$1=""; print $0}' clean_gb_file > input_file
# Creates a single line sequence (suitable for array generation)
string_sequence=$(cat input_file | tr -d " \t\n\r")
#########################################################################################

# Functions
convert_string_into_array $string_sequence 
gene=(${result[@]})
gene_aux=(${gene[@]})
cant=${#gene[@]}

#5'3' Frame 1
frame 0
frame5_3_1=(${result_aa[@]})
#5'3' Frame 2
frame 1
frame5_3_2=(${result_aa[@]})
#5'3' Frame 3
frame 2
frame5_3_3=(${result_aa[@]})


dna_to_complementary
reversed_array=$(echo ${gene_aux[@]} | rev)

#3'5' Frame 1
frameR 0
frame3_5_1=(${result_aa[@]})
#3'5' Frame 2
frameR 1
frame3_5_2=(${result_aa[@]})
#3'5' Frame 3
frameR 2
frame3_5_3=(${result_aa[@]})

amino_seq_final=$(echo "$frame5_3_1" | tr Z "*")
amino_seq_final_2=$(echo "$frame5_3_2" | tr Z "*")
amino_seq_final_3=$(echo "$frame5_3_3" | tr Z "*")

amino_seq_final_4=$(echo "$frame3_5_1" | tr Z "*")
amino_seq_final_5=$(echo "$frame3_5_2" | tr Z "*")
amino_seq_final_6=$(echo "$frame3_5_3" | tr Z "*")


#Cuts every line at row 70 in order to be FASTA compliant!
# ORF - 5.3
################
string_gene_description 0 F > Ex1_output.fas
printf "$amino_seq_final \n" | fold -w 70 >> Ex1_output.fas
string_gene_description 1 F >> Ex1_output.fas
printf "$amino_seq_final_2 \n" | fold -w 70  >> Ex1_output.fas
string_gene_description 2 F >> Ex1_output.fas
printf "$amino_seq_final_3 \n" | fold -w 70  >> Ex1_output.fas
#################
# ORF 3.5
#################
string_gene_description 0 R >> Ex1_output.fas
printf "$amino_seq_final_4 \n" | fold -w 70 >> Ex1_output.fas
string_gene_description 1 R >> Ex1_output.fas
printf "$amino_seq_final_5 \n" | fold -w 70  >> Ex1_output.fas
string_gene_description 2 R >> Ex1_output.fas
printf "$amino_seq_final_6 \n" | fold -w 70  >> Ex1_output.fas
