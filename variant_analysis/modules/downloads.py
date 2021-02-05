
from typing import Optional
import modules.config as config
import json
from Bio import Entrez
import xml.etree.ElementTree as ET
from classes.Genome import Genome

Entrez.email = config.EMAIL
Entrez.api_key = config.API_KEY
Entrez.tool = config.TOOL
Entrez.sleep_between_tries = config.SLEEP_BETWEEN_TRIES


def get_genomes(genome_ids: list[str]) -> list[Genome]:
    handle = Entrez.efetch(
        db="nucleotide", id=genome_ids, rettype="fasta", retmode='xml')
    genomes = __parse_entrez_list(Entrez.read(handle))
    handle.close()

    return genomes


def __parse_entrez_list(entrez_list) -> list[Genome]:
    genomes = []
    for elem in entrez_list:
        genome_id = str(elem['TSeq_accver'])
        tax_id = str(elem['TSeq_taxid'])
        description = str(elem['TSeq_defline'])
        header = f'>{genome_id} {description}'
        sequence = str(elem['TSeq_sequence'])

        genome = Genome(header, genome_id, tax_id, sequence, None)
        genomes.append(genome)
    return genomes


def get_info_genomes(genomes: list[Genome]) -> tuple[list[Genome], Optional[str]]:
    valid_genomes = []
    rejected_ids = []
    msg = None
    for genome in genomes:
        genome = get_info_genome(genome)
        if genome.taxonomy[0] == 'Bacteria' or genome.taxonomy[0] == 'Archaea':
            valid_genomes.append(genome)
        else:
            rejected_ids.append(f'{genome.identifier}({genome.taxonomy[0]})')
    if rejected_ids:
        msg = f'Genomes {", ".join(rejected_ids)} ignored'
    return (valid_genomes, msg)


def get_info_genome(genome: Genome) -> Genome:

    basic_categories = ["superkingdom", "phylum",
                        "class", "order", "family", "genus"]

    handle = Entrez.esummary(
        db="nucleotide", id=genome.identifier, retmode='json')
    genome_info_json = handle.read()
    handle.close()
    info = json.loads(genome_info_json)
    gbid = info['result']['uids'][0]
    if genome.taxonomy_identifier is None:
        genome.taxonomy_identifier = str(info['result'][gbid]['taxid'])
    genome_type = str(info['result'][gbid]['genome'])
    genome.is_plasmid = 'plasmid' in genome_type
    genome.database_name = str(info['result'][gbid]['title'])
    genome.database_length = info['result'][gbid]['slen']

    handle = Entrez.efetch(
        db="taxonomy", id=genome.taxonomy_identifier, retmode='xml')
    genome_taxonomy_xml = handle.read()
    handle.close()
    taxonomy_list = []
    root = ET.fromstring(genome_taxonomy_xml)
    global_taxon = root.find('Taxon')
    strain = global_taxon.find('ScientificName').text
    taxons = global_taxon.find('LineageEx')
    for taxon in taxons.findall('Taxon'):
        rank = taxon.find('Rank').text
        scientificName = taxon.find('ScientificName').text
        if rank in basic_categories:
            taxonomy_list.append(scientificName)
    taxonomy_list.append(strain)
    genome.taxonomy = taxonomy_list

    return genome
