# Bioinformatics

DNA Sequencing script that matches every 3-nucleotide-array, denominated codon, to the corresponding aminoacid.
The 6 possible reading frames, according to ORF, are printed in a fasta formated output file.

To run the script you must grant execution permissions to the user by doing

<code>chmod +x Ex1.sh</code><br>

And then, in order to properly execute it, the messenger DNA file (Genbank or .gb) must be passed as an argument. Following our example, using Tourette syndrome, the execution should be <br>

<code>./Ex1.sh sequence_tourette.gb</code><br>

It is essential to note that the script is only able to process messenger DNA sequences, which have already been cleaned from noncoding regions, called instrons. Here, only coding regions, or exons, are considered.
