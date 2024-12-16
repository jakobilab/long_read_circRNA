#!/usr/bin/env python3

# Copyright (C) 2024 Tobias Jakobi
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either self.version 3 of the License, or
# (at your option) any later self.version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
import argparse
import csv
import gzip
import os
import sys

import pybedtools
import requests
from tqdm import tqdm

import yaml

def is_writeable(directory):
    try:
        with open(os.path.join(directory, 'testfile'), 'w'):
            pass
    except PermissionError:
        return False
    else:
        return True

def get_id_from_column_9(input_string:str, entity:str):

    splits = input_string.strip().split(';')

    feature_dict = {}

    for split in splits:
        item = split.strip().split(" ")
        if len(item) == 2:
            feature_dict[item[0]] = item[1].replace("\"","")

    if entity in feature_dict:
        return feature_dict[entity]
    else:
        return None


def read_annotation_file(annotation_file, entity="exon"):
    """Reads a GTF file
    Will halt the program if file not accessible
    Returns a BedTool object only containing gene sections
    """
    """
    Reads a GTF file and outputs exons output used for the main script

    Download (needs to be the version__lift__version file as GTF, e.g. v37):
    https://www.gencodegenes.org/human/release_37lift37.html

    Expected output (each line one exon):
    chr1    11868   12227   DDX11L1 100     +       ENST00000456328.2_1     ENSG00000223972.5_4

    Args:
        annotation_file (str): Path to the GTF file.
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

                entry = [
                    columns[0],
                    columns[3],
                    columns[4],
                    get_id_from_column_9(columns[8], "gene_name"),
                    "100",
                    columns[6],
                    get_id_from_column_9(columns[8], "transcript_id"),
                    get_id_from_column_9(columns[8], "gene_id")
                ]

                # concatenate lines to one string
                bed_content += '\t'.join(entry) + "\n"

        if not bed_content:
            exit(-1)


        # create a "virtual" BED file
        virtual_bed_file = pybedtools.BedTool(bed_content, from_string=True)

        return virtual_bed_file.sort()



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


def process_refflat_file(file_path):
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



def gencode():
    bed_file = read_annotation_file(sys.argv[1])

    file_base = sys.argv[1].replace(".gtf", "")

    with open(file_base + ".exon.bed", "w") as file:
        file.write(str(bed_file))

    virtual_bed_file_merged = bed_file.merge(s=True,
                                             o=["collapse",
                                                "count",
                                                "distinct"],
                                             c=[4, 4, 6])

    with open(file_base + ".exon.merge.bed", "w") as file:
        file.write(str(virtual_bed_file_merged))

def process_data(configuration: str, data_path: str):
    with open(configuration, 'r') as f:
        config = (yaml.safe_load(f))

        full_data_path=os.path.join(data_path, config['dataset'])

        if is_writeable(data_path):
            print("Writing data to {}".format(data_path))

            # create folder, e.g. h19
            if not os.path.exists(full_data_path):
                os.makedirs(full_data_path)

            for item in config:

                if 'url' in config[item]:

                    url = config[item]['url']
                    file_name = os.path.join(full_data_path, config[item]['name'])

                    file_type = config[item]['type']

                    file_name_unzipped = file_name.replace("."+file_type,
                                                           "")

                    if not os.path.exists(file_name_unzipped):

                        with requests.get(url, stream=True) as r:
                            r.raise_for_status()
                            total_size = int(r.headers.get('content-length', 0))

                            with open(file_name, 'wb') as f:
                                with tqdm(total=total_size, unit='B',
                                          unit_scale=True, desc=file_name) as pbar:
                                    for chunk in r.iter_content(chunk_size=8192):
                                        f.write(chunk)
                                        pbar.update(len(chunk))

                        # most files need to be unpacked
                        if file_type == 'gz':
                            print("Unpacking.")
                            os.system("gzip -d " + file_name)
                            print("Done.")

                        if 'postprocess' in config[item]:
                            print("doing post")

                    # file already exists
                    else:
                        print("Skipping, " + config[item]['name'] +
                              ", file already exists.")


def main():

    parser = argparse.ArgumentParser(
        description="Download required data from Nanopore circRNA pipeline")

    # REQUIRED ARGUMENTS
    group = parser.add_argument_group("Required options")

    group.add_argument("-c",
                       "--config",
                       dest="config",
                       help="Path to the configuration file, e.g. hg19.yml",
                       required=True
                       )

    group.add_argument("-p",
                       "--path",
                       dest="path",
                       help="Path to store the downloaded data files",
                       required=True
                       )

    args = parser.parse_args(sys.argv[1:])

    if not os.path.isfile(args.config):
        print("Configuration file not accessible.")
    else:
        process_data(configuration=args.config, data_path=args.path)

if __name__ == "__main__":
    main()
