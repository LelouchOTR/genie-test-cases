## Fastq

| Test #   | Test                                  | Description                                                                                              |
| :------: | :------------------------------------ | :------------------------------------------------------------------------------------------------------- |
| FASTQ_01 | Single End – constant read length     | Three reads in a single fastq file with constant read length.                                            |
| FASTQ_02 | Single End – variable read length     | Three reads in a single fastq file with the second read having a different length than the other two.    |
| FASTQ_03 | Paired End - both mates same length   | Two paired fastq files with three reads in each file. The read length is constant for all reads.         |
| FASTQ_04 | Paired End - mates different length   | Two paired fastq files with three reads in each file. The read length is constant in each file but different between the files. |
| FASTQ_05 | (fastq.gz input)                      | Three reads in a single fastq file with constant read length, but gzipped.                               |
| FASTQ_06 | (fastq.gz output)                     | Three reads in a single fastq file with constant read length, but the decoded output is gzipped.         |
| FASTQ_07 | Paired End - different read names     | Two paired fastq files with three reads in each file. The read names between mates differs.              |
| FASTQ_08 | Alphabet test Quality score range     | A single fastq file with a single read containing the full range of allowed quality values.              |
| FASTQ_09 | Alphabet test nucleotide alphabet ACTG  | A single fastq file with a single read containing the full range of the ACGT nucleotide alphabet         |
| FASTQ_10 | Alphabet test nucleotide alphabet IUPAC | A single fastq file with a single read containing the full range of the IUPAC nucleotide alphabet        |
| FASTQ_11 | Read name alphabet test               | A single fastq file with a single read containing the full range of the allowed read name characters     |
| FASTQ_12 | Paired End - unequal file lengths     | A paired end fastq file, but there are reads at the end of one file with no mate in the other file.      |

## SAM/BAM

| Test #  | Test                                            | Description                                                                                                                               |
| :-----: | :---------------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------------- |
| SAM_01  | Unmapped read – single end                      | Sam file containing one unpaired and unmapped read.                                                                                       |
| SAM_02  | Unmapped pair                                   | Sam file containing three unpaired and unmapped read pairs. The second pair has the „first read" flag in the other order as the other two pairs |
| SAM_03  | Half-mapped read pair                           | Sam file containing three half mapped read pairs. The second pair has the „first read" flag in the other order as the other two pairs      |
| SAM_04  | Mapped read single end                          | Sam file containing three unpaired and mapped reads, all perfect matches to the reference sequence. The second pair has the „first read" flag in the other order as the other two pairs |
| SAM_05  | Mapped read pair – same position + TLEN         | Sam file containing a perfectly mapped pair. The mates are mapped to the same position                                                    |
| SAM_06  | Mapped read pair – enclosed + TLEN              | Sam file containing a perfectly mapped pair. One read encloses the other read, i.a. read A begins before read B begins and ends after the end of read B. |
| SAM_07  | Mapped read pair – overlapping + TLEN           | Sam file containing a perfectly mapped pair. One read overlaps the other read, i.a. read A begins before read B begins and ends before the end of read B. |
| SAM_08  | Mapped read pair – no overlapping + TLEN        | Sam file containing a perfectly mapped pair. One read does not overlap the other read, i.a. read A ends before read B has started         |
| SAM_09  | Mapped read pair – long distance + TLEN         | Sam file containing a perfectly mapped pair. One read does not overlap the other read, i.a. read A ends before read B has started. Distance between the two reads is > 1000000 |
| SAM_10  | Mapped read pair – different reference + TLEN   | Sam file containing a perfectly mapped pair. The mates map to different reference sequences.                                              |
| SAM_11  | Secondary alignment                             | Sam file containing a perfectly mapped single read and a secondary alignment for the same read.                                           |
| SAM_12  | Supplementary / chimeric alignment              | Sam file containing a perfectly mapped chimeric read with at least supplementary alignments.                                              |
| SAM_13  | Base substitution (M, =, X)                     | Sam file containing a read with M, =, X in the CIGAR string                                                                               |
| SAM_14  | Base insertion                                  | Sam file containing a read with I in the CIGAR string                                                                                     |
| SAM_15  | Base deletion                                   | Sam file containing a read with D in the CIGAR string                                                                                     |
| SAM_16  | Softclips                                       | Sam file containing a read with S in the CIGAR string                                                                                     |
| SAM_17  | Padding (P)                                     | Sam file containing a read with P in the CIGAR string                                                                                     |
| SAM_18  | Hardclips                                       | Sam file containing a read with H in the CIGAR string                                                                                     |
| SAM_19  | PCR duplicate flag                              | Sam file with three reads, each containing different values for the flags „PCR duplicate", „Quality checks" and „Properly aligned"        |
| SAM_20  | Paired end – different flags per mate           | A pair where each mate contains a different combination of the flags „PCR duplicate", „Quality checks" and „Properly aligned"             |
| SAM_21  | next read … flags – unmapped                    |                                                                                                                                           |
| SAM_22  | next read … flags – half mapped                 |                                                                                                                                           |
| SAM_23  | next read … flags – short distance              |                                                                                                                                           |
| SAM_24  | next read … flags – long distance               |                                                                                                                                           |
| SAM_25  | Short intron / splice (N)                       |                                                                                                                                           |
| SAM_26  | Long intron / splice (N)                        |                                                                                                                                           |
| SAM_27  | Empty read (all bases deleted)                  |                                                                                                                                           |
| SAM_28  | Empty read (all bases softclipped)              |                                                                                                                                           |
| SAM_29  | Empty read (all bases hardclipped)              |                                                                                                                                           |
| SAM_30  | Empty read (no nucleotides in read / * in sam)  |                                                                                                                                           |
| SAM_31  | Quality scores absent                           |                                                                                                                                           |
| SAM_32  | Optional tags                                   |                                                                                                                                           |
| SAM_33  | Read groups                                     |                                                                                                                                           |
| SAM_34  | Reverse Complement (different + same) – short distance |                                                                                                                                           |
| SAM_35  | Reverse Complement (different + same) – long distance |                                                                                                                                           |
| SAM_36  | Reverse Complement (different + same) – unmapped |                                                                                                                                           |
| SAM_37  | Reverse Complement (different + same) – half mapped |                                                                                                                                           |
| SAM_38  | Circular reference                              | Sam file with a read overlapping the end of a circular reference.                                                                         |
| SAM_39  | (bam input)                                     |                                                                                                                                           |
| SAM_40  | (bam output)                                    |                                                                                                                                           |
| SAM_41  | (cram input)                                    |                                                                                                                                           |
| SAM_42  | (cram output)                                   |                                                                                                                                           |
