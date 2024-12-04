# Copyright (C) 2024 Tobias Jakobi
#
# @Author: Tobias Jakobi <tjakobi>
# @Email:  tjakobi@arizona.edu
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import csv
import gzip
import sys
import pybedtools


def read_annotation_file(annotation_file, entity="exon", string=False):
    """Reads a GTF file
    Will halt the program if file not accessible
    Returns a BedTool object only containing gene sections
    """

    try:
        file_handle = open(annotation_file)
    except PermissionError:
        message = ("Input file " + str(
            annotation_file) + " cannot be read, exiting.")
        sys.exit(message)
    else:

        with file_handle:
            line_iterator = iter(file_handle)
            bed_content = ""
            print("Start parsing GTF file")
            for line in line_iterator:
                # we skip any comment lines
                if line.startswith("#"):
                    continue

                # split up the annotation line
                columns = line.split('\t')

                if not (columns[2] == entity):
                    continue

                # we do not want any 0-length intervals -> bedtools segfault
                if int(columns[4]) - int(columns[3]) == 0:
                    continue

                # extract chromosome, start, stop, score(0), name and strand
                # we hopefully have a gene name now and use this one for the entry


                entry = [
                    columns[0],
                    columns[3],
                    columns[4],
                    "name",
                    str(0),
                    columns[6],
                ]
                print(entry)

                # TODO: function to split up additional arguments for gene id

                # concatenate lines to one string
                bed_content += '\t'.join(entry) + "\n"

        # if not bed_content:
        #     exit(-1)
        #
        # if string:
        #     return bed_content
        # else:
        #     # create a "virtual" BED file
        #     virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)
        #     print("Start merging GTF file")
        #
        #     # we trust that bedtools >= 2.27 is installed. Otherwise this merge will probably fail
        #     return virtual_bed_file.sort().merge(s=True,  # strand specific
        #                                          c="4,5,6",
        #                                          # copy columns 5 & 6
        #                                          o="distinct,distinct,distinct")  # group

def process_gtf_file(file_path):
    """
    Reads a refFlat file dumped by UCSC Genome Browser and outputs sorted
    exon output used for the main script

    refFlat schema:
    https://genome.ucsc.edu/cgi-bin/hgTables?hgta_doSchemaDb=hg38&hgta_doSchemaTable=refFlat

    Expected output:
    chr1    11868   12227   LOC102725121_exon_0_0_chr1_11869_f      0       +

    name = genename + exon # + 0 + chr + start+1 + strand (f or r)


    Args:
        file_path (str): Path to the input file.
    """
    output_exons = ""
    output_genes = ""

    try:
        with gzip.open(file_path, mode='rt') as gz_file:
            csv_reader = csv.reader(gz_file, delimiter='\t')
            next(csv_reader, None)

            for row_number, row in enumerate(csv_reader, start=1):
                geneName, name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, exonCount, exonStarts, exonEnds = row

                starts = exonStarts.split(',')
                stops = exonEnds.split(',')

                output_genes += "\t".join(
                  [chrom, txStart, txEnd, geneName, str(0), strand]) + "\n"

                for exon_num in range(int(exonCount)-1):
                    strand_tag  = "f" if strand == "+" else "r"

                    name_tag = "_".join([geneName,"exon",str(exon_num),str(0),chrom,str(int(starts[exon_num])+1),strand_tag])

                    output_genes += "\t".join([chrom,starts[exon_num],stops[exon_num],geneName,str(0),strand]) + "\n"
                    #output_genes += "\t".join([chrom,starts[exon_num],stops[exon_num],name_tag,str(0),strand]) + "\n"
                    output_exons += "\t".join([chrom,starts[exon_num],stops[exon_num],name_tag,str(0),strand]) + "\n"

        virtual_bed_file = pybedtools.BedTool(output_exons, from_string=True)

        virtual_bed_file_genes = pybedtools.BedTool(output_genes, from_string=True).sort()

        # the unique gene level file is different in regard to printing out all genes with comma instead of just choosing one
        # and omit the other co-optimal hits

        virtual_bed_file_genes = virtual_bed_file_genes.merge(s=True, o= ["distinct","count_distinct","distinct"], c=[4,4,6])

        virtual_bed_file_sorted = virtual_bed_file.sort()

        virtual_bed_file_merged = virtual_bed_file_sorted.merge(s=True, o= ["collapse","count","distinct"], c=[4,4,6])

        file_base = file_path.replace(".gz","")

        with open(file_base+".sort.bed", "w") as file:
            file.write(str(virtual_bed_file_sorted))

        with open(file_base+".unique.bed", "w") as file:
            file.write(str(virtual_bed_file_genes))

        with open(file_base+".merged.bed", "w") as file:
            file.write(str(virtual_bed_file_merged ))

    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
    except Exception as e:
        print(f"An error occurred: {e}")


if __name__ == "__main__":
    read_annotation_file(sys.argv[1])