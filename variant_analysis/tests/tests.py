from downloads import get_genome, get_taxonomy, get_gen16s
from outputs import generate_id, generate_fasta, save_output
from views import genes16s_comparison
import pathlib
import shutil
import os
import uuid


# Create your tests here.

"""
class ExampleUnitTest(SimpleTestCase):

    def setUp(self):
        if not os.path.exists('downloaded_files/'):
            os.makedirs('downloaded_files/')

    def test_get_genome(self):
        current_user = uuid.uuid4()
        downloaded_files_path = pathlib.Path(
            current_dir).joinpath("downloaded_files")
        # se crea la carpeta para el usuario
        downloaded_files_path.joinpath(str(current_user)).mkdir()

        get_genome(str(current_user), ['CP000863'], [])

        ruta_archivo = downloaded_files_path.joinpath(
            str(current_user)).joinpath('NCBI_Download.fa')

        self.assertTrue(os.path.isfile(ruta_archivo))

        self.assertNotEquals(os.stat(ruta_archivo).st_size, 0)

    def test_get_taxonomy(self):
        current_user = uuid.uuid4()
        downloaded_files_path = pathlib.Path(
            current_dir).joinpath("downloaded_files")
        # se crea la carpeta para el usuario
        downloaded_files_path.joinpath(str(current_user)).mkdir()
        ruta_archivo = downloaded_files_path.joinpath(
            str(current_user)).joinpath("listaID.csv")
        with open(ruta_archivo, "w+") as f:
            f.write('CP000863')
        taxonomy = {}
        taxid_genbankid = {}
        rejected_ids = []
        valid_ids = []
        get_taxonomy(taxonomy, taxid_genbankid, str(
            current_user), rejected_ids, valid_ids)

        self.assertNotEquals(taxonomy, {})
        self.assertNotEquals(taxid_genbankid, {})
        self.assertNotEquals(valid_ids, {})

    def test_get_genes16s(self):
        current_user = uuid.uuid4()
        downloaded_files_path = pathlib.Path(
            current_dir).joinpath("downloaded_files")
        # se crea la carpeta para el usuario
        downloaded_files_path.joinpath(str(current_user)).mkdir()
        ruta_archivo = downloaded_files_path.joinpath(
            str(current_user)).joinpath("NCBI_Download.fa")
        ruta_origen = 'aplicacion/test_files/NCBI_Download.fa'
        shutil.copyfile(ruta_origen, ruta_archivo)
        genes16s = {}
        genes16s_strand = {}
        genes16s_begin_end = {}
        training_db = 'HOMD'

        get_gen16s(genes16s, training_db, genes16s_strand,
                   genes16s_begin_end, str(current_user))

        self.assertTrue(genes16s)
        self.assertTrue(genes16s_strand)
        self.assertTrue(genes16s_begin_end)

    def test_generate_output_files(self):
        current_user = uuid.uuid4()
        downloaded_files_path = pathlib.Path(
            current_dir).joinpath("downloaded_files")

        taxonomy_variants = {'405416': [['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'], [
            'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']]}
        output_files_ids = {}

        generate_id(taxonomy_variants, output_files_ids)

        self.assertTrue(output_files_ids)

        # se crea la carpeta para el usuario
        downloaded_files_path.joinpath(str(current_user)).mkdir()
        taxonomy_variants = {'405416': [['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'], [
            'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']]}
        genes_variants = {'405416': {'v1': 'AGAGTGATAG', 'v2': 'AGTAGATAGAT'}}
        ruta_archivo = downloaded_files_path.joinpath(
            str(current_user)).joinpath('fastaTaxonomy.fa')

        generate_fasta(taxonomy_variants, output_files_ids,
                       genes_variants, current_user)

        self.assertTrue(os.path.isfile(ruta_archivo))
        self.assertNotEquals(os.stat(ruta_archivo).st_size, 0)
        taxonomy_variants = {
            '405416': [['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'], ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']]}
        genes_variants = {'405416': {'v1': 'AGAGTGATAG', 'v2': 'AGTAGATAGAT'}}
        output_files_ids = {'405416': {'v1': '123_1233', 'v2': '123_1234'}}
        output_formats = ['mothur', 'dada2', 'qiime']
        save_output(taxonomy_variants, output_files_ids,
                    genes_variants, output_formats, current_user)

        ruta_archivo = downloaded_files_path.joinpath(
            str(current_user)).joinpath('mothurTaxonomy.taxonomy')
        self.assertTrue(os.path.isfile(ruta_archivo))
        self.assertNotEquals(os.stat(ruta_archivo).st_size, 0)
        ruta_archivo = downloaded_files_path.joinpath(
            str(current_user)).joinpath('dada2Taxonomy.fa')
        self.assertTrue(os.path.isfile(ruta_archivo))
        self.assertNotEquals(os.stat(ruta_archivo).st_size, 0)
        ruta_archivo = downloaded_files_path.joinpath(
            str(current_user)).joinpath('qiimeTaxonomy.taxonomy')
        self.assertTrue(os.path.isfile(ruta_archivo))
        self.assertNotEquals(os.stat(ruta_archivo).st_size, 0)

    # def test_obtain_variants(self):
    # 	taxonomy_variants = {'405416': [['a','b','c', 'd', 'e','f','g','h'], ['a','b','c', 'd', 'e','f','g','h']]}
    # 	genes_variants = {'405416': {'v1': 'AGAGTGATAG', 'v2': 'AGTAGATAGAT'}}
    # 	taxid_genbankid = {'CP000863':'405416'}
    # 	rejected_ids = []
    # 	genes16s = {'405416': ['AGAGTGATAG','AGAGTGATAG', 'AGAGTGATAG']}
    # 	variants_gb_id = {}
    # 	taxonomy={'405416': ['a','b','c', 'd', 'e','f','g','h']}
    # 	variants_strand={}
    # 	genes16s_strand ={'405416': ['+','-','+']}
    # 	unique_number_genes16s={}
    # 	total_number_genes16s={}
    # 	genes16s_comparison(taxid_genbankid, rejected_ids, genes16s, variants_gb_id, taxonomy, variants_strand, genes_variants, taxonomy_variants, genes16s_strand, unique_number_genes16s,total_number_genes16s)

    def test_send_email(self):
        current_user = uuid.uuid4()
        email = 'laramaria.vazquez.gonzales@rai.usc.es'
        output_formats = ['dada2']
        downloaded_files_path = pathlib.Path(
            current_dir).joinpath("downloaded_files")
        # se crea la carpeta para el usuario
        downloaded_files_path.joinpath(str(current_user)).mkdir()
        ruta_archivo = downloaded_files_path.joinpath(
            str(current_user)).joinpath("dada2Taxonomy.fa")
        ruta_archivo2 = downloaded_files_path.joinpath(
            str(current_user)).joinpath("statistics2.csv")
        ruta_archivo3 = downloaded_files_path.joinpath(
            str(current_user)).joinpath("statistics1.csv")

        with open(ruta_archivo, "w+") as f, open(ruta_archivo2, "w+") as f1, open(ruta_archivo3, "w+") as f2:
            send_mail(email, output_formats, current_user, False)

    def tearDown(self):
        downloaded_files_path = pathlib.Path(
            current_dir).joinpath("downloaded_files")
        shutil.rmtree(downloaded_files_path, ignore_errors=True)
"""
