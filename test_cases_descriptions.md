## Fastq

| Test                                  | Description                                                                                              |
| :------------------------------------ | :------------------------------------------------------------------------------------------------------- |
| Single End – constant read length     | Three reads in a single fastq file with constant read length.                                            |
| Single End – variable read length     | Three reads in a single fastq file with the second read having a different length than the other two.    |
| Paired End - both mates same length   | Two paired fastq files with three reads in each file. The read length is constant for all reads.         |
| Paired End - mates different length   | Two paired fastq files with three reads in each file. The read length is constant in each file but different between the files. |
| (fastq.gz input)                      | Three reads in a single fastq file with constant read length, but gzipped.                               |
| (fastq.gz output)                     | Three reads in a single fastq file with constant read length, but the decoded output is gzipped.         |
| Paired End - different read names     | Two paired fastq files with three reads in each file. The read names between mates differs.              |
| Alphabet test Quality score range     | A single fastq file with a single read containing the full range of allowed quality values.              |
| Alphabet test nucleotide alphabet ACTG  | A single fastq file with a single read containing the full range of the ACGT nucleotide alphabet         |
| Alphabet test nucleotide alphabet IUPAC | A single fastq file with a single read containing the full range of the IUPAC nucleotide alphabet        |
| Read name alphabet test               | A single fastq file with a single read containing the full range of the allowed read name characters     |
| Paired End - unequal file lengths     | A paired end fastq file, but there are reads at the end of one file with no mate in the other file.      |

## SAM/BAM

| Test                                            | Description                                                                                                                               |
| :---------------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------------- |
| Unmapped read – single end                      | Sam file containing one unpaired and unmapped read.                                                                                       |
| Unmapped pair                                   | Sam file containing three unpaired and unmapped read pairs. The second pair has the „first read" flag in the other order as the other two pairs |
| Half-mapped read pair                           | Sam file containing three half mapped read pairs. The second pair has the „first read" flag in the other order as the other two pairs      |
| Mapped read single end                          | Sam file containing three unpaired and mapped reads, all perfect matches to the reference sequence. The second pair has the „first read" flag in the other order as the other two pairs |
| Mapped read pair – same position + TLEN         | Sam file containing a perfectly mapped pair. The mates are mapped to the same position                                                    |
| Mapped read pair – enclosed + TLEN              | Sam file containing a perfectly mapped pair. One read encloses the other read, i.a. read A begins before read B begins and ends after the end of read B. |
| Mapped read pair – overlapping + TLEN           | Sam file containing a perfectly mapped pair. One read overlaps the other read, i.a. read A begins before read B begins and ends before the end of read B. |
| Mapped read pair – no overlapping + TLEN        | Sam file containing a perfectly mapped pair. One read does not overlap the other read, i.a. read A ends before read B has started         |
| Mapped read pair – long distance + TLEN         | Sam file containing a perfectly mapped pair. One read does not overlap the other read, i.a. read A ends before read B has started. Distance between the two reads is > 1000000 |
| Mapped read pair – different reference + TLEN   | Sam file containing a perfectly mapped pair. The mates map to different reference sequences.                                              |
| Secondary alignment                             | Sam file containing a perfectly mapped single read and a secondary alignment for the same read.                                           |
| Supplementary / chimeric alignment              | Sam file containing a perfectly mapped chimeric read with at least supplementary alignments.                                              |
| Base substitution (M, =, X)                     | Sam file containing a read with M, =, X in the CIGAR string                                                                               |
| Base insertion                                  | Sam file containing a read with I in the CIGAR string                                                                                     |
| Base deletion                                   | Sam file containing a read with D in the CIGAR string                                                                                     |
| Softclips                                       | Sam file containing a read with S in the CIGAR string                                                                                     |
| Padding (P)                                     | Sam file containing a read with P in the CIGAR string                                                                                     |
| Hardclips                                       | Sam file containing a read with H in the CIGAR string                                                                                     |
| PCR duplicate flag                              | Sam file with three reads, each containing different values for the flags „PCR duplicate", „Quality checks" and „Properly aligned"        |
| Paired end – different flags per mate           | A pair where each mate contains a different combination of the flags „PCR duplicate", „Quality checks" and „Properly aligned"             |
| next read … flags – unmapped                    |                                                                                                                                           |
| next read … flags – half mapped                 |                                                                                                                                           |
| next read … flags – short distance              |                                                                                                                                           |
| next read … flags – long distance               |                                                                                                                                           |
| Short intron / splice (N)                       |                                                                                                                                           |
| Long intron / splice (N)                        |                                                                                                                                           |
| Empty read (all bases deleted)                  |                                                                                                                                           |
| Empty read (all bases softclipped)              |                                                                                                                                           |
| Empty read (all bases hardclipped)              |                                                                                                                                           |
| Empty read (no nucleotides in read / * in sam)  |                                                                                                                                           |
| Quality scores absent                           |                                                                                                                                           |
| Optional tags                                   |                                                                                                                                           |
| Read groups                                     |                                                                                                                                           |
| Reverse Complement (different + same) – short distance |                                                                                                                                           |
| Reverse Complement (different + same) – long distance |                                                                                                                                           |
| Reverse Complement (different + same) – unmapped |                                                                                                                                           |
| Reverse Complement (different + same) – half mapped |                                                                                                                                           |
| Circular reference                              | Sam file with a read overlapping the end of a circular reference.                                                                         |
| (bam input)                                     |                                                                                                                                           |
| (bam output)                                    |                                                                                                                                           |
| (cram input)                                    |                                                                                                                                           |
| (cram output)                                   |                                                                                                                                           |
