import os
import urllib.request
import pathlib
import modules.config as config
import time
import json
import regex as re
from datetime import datetime
from pathlib import Path
from Bio import Entrez
import xml.etree.ElementTree as ET
from classes.Genome import Genome


url_tax_id = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nucleotide&id={}&retmode=json&api_key={}&email={}&tool={}"

url_tax = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id={}&api_key={}&email={}&tool={}"

Entrez.email = config.EMAIL
Entrez.api_key = config.API_KEY
Entrez.tool = config.TOOL
Entrez.sleep_between_tries = config.SLEEP_BETWEEN_TRIES


def get_genomes(genome_ids):  # se obtiene el genoma del NCBI
    genomes = []
    for genome_id in genome_ids:
        handle = Entrez.efetch(
            db="nucleotide", id=genome_id.strip(), rettype="fasta", retmode='xml')
        genome = __parse_fasta_xml(handle.read())
        genomes.append(genome)

    return genomes


# Método para el guardado de un FASTA de un genoma en un objeto de la clase Genome
def __parse_fasta_xml(fasta_xml):

    root = ET.fromstring(fasta_xml)
    genome_id = root.find('TSeq_accver').text
    tax_id = root.find('TSeq_taxid').text
    description = root.find('TSeq_defline').text
    header = f'>{genome_id} {description}'
    sequence = root.find('TSeq_sequence').text

    genome = Genome(header, genome_id, tax_id, sequence, None)

    return genome

# se obtiene la taxonomía


def get_info_genomes(genomes: list[Genome]):

    for genome in genomes:
        genome = get_info_genome(genome)

    return genomes


def get_info_genome(genome: Genome):

    basic_categories = ["superkingdom", "phylum",
                        "class", "order", "family", "genus"]

    if genome.taxonomy_identifier is None:
        handle = Entrez.esummary(
            db="nucleotide", id=genome.identifier, retmode='json')
        genome_info_json = handle.read()

        info = json.loads(genome_info_json)
        gbid = info['result']['uids'][0]
        genome.taxonomy_identifier = str(info['result'][gbid]['taxid'])
        genome_type = str(info['result'][gbid]['genome'])
        genome.is_plasmid = 'plasmid' in genome_type
        genome.database_name = str(info['result'][gbid]['title'])
        genome.database_length = info['result'][gbid]['slen']

    handle = Entrez.efetch(
        db="taxonomy", id=genome.taxonomy_identifier, retmode='xml')
    genome_taxonomy_xml = handle.read()

    taxonomy_list = []
    root = ET.fromstring(genome_taxonomy_xml)
    global_taxon = root.find('Taxon')
    strain = global_taxon.find('ScientificName')
    taxons = global_taxon.find('LineageEx')
    for taxon in taxons.findall('Taxon'):
        rank = taxon.find('Rank').text
        scientificName = taxon.find('ScientificName').text
        # if rank.text == 'superkingdom' and (name.text != 'Bacteria' and name.text != 'Archaea'):
        #             rejected_ids.append(genome.identifier)
        if rank in basic_categories:
            taxonomy_list.append(scientificName)
    taxonomy_list.append(strain)
    genome.taxonomy = taxonomy_list

    return genome
