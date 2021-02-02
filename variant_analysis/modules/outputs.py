import os
import random
import hashlib
#from downloads import *


# función donde se generan los identificadores para los archivos
def generate_id(taxonomy_variants, genes_variants, output_files_ids):
    for id in taxonomy_variants:
        variants = taxonomy_variants.get(id)
        variant_genes = genes_variants.get(id)
        id_variants = {}
        for variant in variants:
            gene = id + variant_genes.get(variant[-1])  # TODO: cambiar a GBid
            hashed_id = str(int(hashlib.sha256(
                gene.encode('utf-8')).hexdigest(), 16) % 10**6).zfill(6)
            id_variants[variant[-1]] = hashed_id
        output_files_ids[id] = id_variants


def save_output(taxonomy_variants, output_files_ids, genes_variants, output_formats, taxid_genbankid, output_dir_path):

    is_fasta_generated = False
    for format in output_formats:  # se recorren los formatos seleccionados por el usuario
        taxid = list(taxonomy_variants.keys())[0]
        GBIDs = '_'.join([x.identifier for x in taxid_genbankid[taxid]])
        if format == 'mothur':
            mothur_dir_path = '{}/mothur/{}.taxonomy'.format(
                output_dir_path, GBIDs)
            with open(mothur_dir_path, 'w+') as taxonomy_file:
                for id in taxonomy_variants:
                    variants = taxonomy_variants.get(id)
                    for variant in variants:
                        line = '>{}   '.format(
                            (output_files_ids.get(id)).get(variant[-1]))
                        for elem in variant:
                            line += '{};'.format(elem)
                        taxonomy_file.write('{}\n'.format(line))
                if is_fasta_generated is False:
                    generate_fasta(taxonomy_variants, output_files_ids,
                                   genes_variants, taxid_genbankid, output_dir_path)
                    is_fasta_generated = True

        elif format == 'dada2':
            dada2_dir_path = '{}/dada2/{}.fa'.format(output_dir_path, GBIDs)

            with open(dada2_dir_path, 'w+') as taxonomy_file:
                for id in taxonomy_variants:
                    variants = taxonomy_variants.get(id)
                    for variant in variants:
                        line = ">"
                        for elem in variant:
                            line += "{};".format(elem)
                        taxonomy_file.write("{}\n".format(line))
                        variant_gene = genes_variants.get(id)
                        taxonomy_file.write("{}\n".format(
                            variant_gene.get(variant[-1])))

        elif format == "qiime":
            qiime_categories = ["k__", "p__",
                                "c__", "o__", "f__", "g__", "s__"]

            qiime_dir_path = '{}/qiime/{}.taxonomy'.format(
                output_dir_path, GBIDs)

            with open(qiime_dir_path, 'w+') as taxonomy_file:
                for id in taxonomy_variants:
                    variants = taxonomy_variants.get(id)
                    for variant in variants:
                        line = ">{}   ".format(
                            (output_files_ids.get(id)).get(variant[-1]))
                        i = 0
                        for elem in variant:
                            if i < 6:
                                line += "{0}{1};".format(
                                    qiime_categories[i], elem)
                            elif i == 6:
                                line += " {0}{1}".format(
                                    qiime_categories[i], elem)
                            else:
                                line += "{}".format(elem)
                            i = i + 1
                        taxonomy_file.write("{}\n".format(line))
            if is_fasta_generated is False:
                generate_fasta(taxonomy_variants, output_files_ids,
                               genes_variants, taxid_genbankid, output_dir_path)
                is_fasta_generated = True


# se genera el archivo común a mothur y qiime con los genes
def generate_fasta(taxonomy_variants, output_files_ids, genes_variants, taxid_genbankid, output_dir_path):
    taxid = list(taxonomy_variants.keys())[0]
    GBIDs = '_'.join([x.identifier for x in taxid_genbankid[taxid]])
    fasta_dir_path = '{}/fasta/{}.fa'.format(
        output_dir_path, GBIDs)

    with open(fasta_dir_path, 'w+') as fasta_file:
        for id in taxonomy_variants:
            variants = taxonomy_variants.get(id)
            for variant in variants:
                line = ">{0} | {1} | {2} \n".format((output_files_ids.get(id)).get(variant[-1]), variant[-2],
                                                    variant[-1])
                fasta_file.write(line)
                variant_gene = genes_variants.get(id)
                fasta_file.write("{}\n".format(variant_gene.get(variant[-1])))
