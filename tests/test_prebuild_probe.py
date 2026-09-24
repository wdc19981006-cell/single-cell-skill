"""Offline pre-build schema and pooled mapping regressions."""
import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path

import h5py

SCRIPTS = Path(__file__).resolve().parents[1] / '.agents/skills/geo-single-cell-loader/scripts'
sys.path.insert(0, str(SCRIPTS))
from common import digest, read_csv, verify_confirmation, write_csv
from prebuild_probe import run
from text_schema_probe import select_text_reader


def row(gse, sample, path, kind='text'):
    return dict(database=gse, sample=sample, tissue='Synthetic', disease='Synthetic',
                source_type='Synthetic', group='Fixture', local_path=f'data/{gse}/raw/{path}',
                file_type=kind, count_source='counts', count_evidence='fixture raw counts',
                metadata_evidence='fixture source', files_json='[]', delimiter='',
                orientation='', feature_column='', drop_columns='')


class PrebuildProbeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.root = Path(self.temporary.name)

    def tearDown(self):
        self.temporary.cleanup()

    def setup_manifest(self, gse, rows):
        workflow = self.root / 'data' / gse / '.workflow'
        workflow.mkdir(parents=True)
        manifest = workflow / 'sample_manifest.csv'
        write_csv(manifest, rows)
        (workflow / 'sample_manifest.confirmation.json').write_text(json.dumps({
            'confirmed_by': 'user', 'user_statement': 'SIMULATED APPROVAL ONLY',
            'groups': {item['sample']: item['group'] for item in rows},
            'manifest_md5': digest(manifest),
        }), encoding='utf-8')
        return workflow, manifest

    def test_gene_id_column_is_completed_before_build(self):
        gse = 'GSE999999999'
        raw = self.root / 'data' / gse / 'raw'
        raw.mkdir(parents=True)
        (raw / 'counts.csv').write_text('gene_id,AAACCTGAGGCTACGA-1,AAACCTGAGGTTCCTA-1\n'
                                        'ENSG000001,0,2\nENSG000002,1,0\n', encoding='utf-8')
        workflow, manifest = self.setup_manifest(gse, [row(gse, 'GSM1', 'counts.csv')])
        run(manifest, self.root)
        completed = read_csv(manifest)[0]
        self.assertEqual((completed['delimiter'], completed['orientation'], completed['feature_column']),
                         ('comma', 'genes_by_cells', 'gene_id'))
        self.assertEqual(completed['text_reader'], 'fread')
        self.assertEqual((int(completed['matrix_rows']), int(completed['matrix_columns'])), (2, 2))
        self.assertEqual(int(completed['estimated_dense_bytes']), 32)
        self.assertGreater(int(completed['physical_ram_bytes']), 0)
        self.assertIn('below 25%', completed['reader_selection_reason'])
        verify_confirmation(manifest)
        self.assertEqual(json.loads((workflow / 'text_schema_probe.json').read_text())[0]['status'],
                         'VERIFIED_SAMPLE')

    def test_wide_comma_header_skips_wrong_delimiter(self):
        gse = 'GSE999999999'
        raw = self.root / 'data' / gse / 'raw'
        raw.mkdir(parents=True)
        cells = [f'AAACCTGAGGCTACGA-{index:05d}' for index in range(9000)]
        (raw / 'wide.csv').write_text('gene_id,' + ','.join(cells) + '\n' +
                                       'ENSG000001,' + ','.join(['0'] * len(cells)) + '\n',
                                       encoding='utf-8')
        _, manifest = self.setup_manifest(gse, [row(gse, 'GSM1', 'wide.csv')])
        run(manifest, self.root)
        completed = read_csv(manifest)[0]
        self.assertEqual((completed['delimiter'], completed['feature_column']), ('comma', 'gene_id'))

    def test_normalized_text_is_rejected_before_build(self):
        gse = 'GSE999999999'
        raw = self.root / 'data' / gse / 'raw'
        raw.mkdir(parents=True)
        (raw / 'counts.tsv').write_text('Gene\tAAACCTGAGGCTACGA-1\nENSG000001\t0.25\n',
                                        encoding='utf-8')
        workflow, manifest = self.setup_manifest(gse, [row(gse, 'GSM1', 'counts.tsv')])
        with self.assertRaisesRegex(ValueError, 'Normalized/processed'):
            run(manifest, self.root)
        self.assertEqual(read_csv(manifest)[0]['orientation'], '')
        self.assertEqual(json.loads((workflow / 'text_schema_probe.json').read_text())[0]['status'], 'STOP')

    def test_header_without_gene_id_keeps_existing_txt_route(self):
        gse = 'GSE999999999'
        raw = self.root / 'data' / gse / 'raw'
        raw.mkdir(parents=True)
        (raw / 'counts.txt').write_text('AAACCTGAGGCTACGA-1\tAAACCTGAGGTTCCTA-1\n'
                                        'A1BG\t0\t2\nA2M\t1\t0\n', encoding='utf-8')
        workflow, manifest = self.setup_manifest(gse, [row(gse, 'GSM1', 'counts.txt')])
        run(manifest, self.root)
        completed = read_csv(manifest)[0]
        self.assertEqual(completed['feature_column'], '__row_names__')
        self.assertEqual(completed['orientation'], 'genes_by_cells')
        self.assertTrue(json.loads((workflow / 'text_schema_probe.json').read_text())[0]['header_missing_id'])

    def test_blank_csv_id_header_is_not_treated_as_missing_column(self):
        gse = 'GSE999999999'
        raw = self.root / 'data' / gse / 'raw'
        raw.mkdir(parents=True)
        (raw / 'counts.csv').write_text('"",AAACCTGAGGCTACGA-1,AAACCTGAGGTTCCTA-1\n'
                                        'A1BG,0,2\nA2M,1,0\n', encoding='utf-8')
        _, manifest = self.setup_manifest(gse, [row(gse, 'GSM1', 'counts.csv')])
        run(manifest, self.root)
        completed = read_csv(manifest)[0]
        self.assertEqual(completed['feature_column'], '__row_names__')
        self.assertEqual(completed['text_header_missing_id'], 'false')

    def test_explicit_first_column_streaming_route_is_preserved(self):
        gse = 'GSE999999999'
        raw = self.root / 'data' / gse / 'raw'
        raw.mkdir(parents=True)
        (raw / 'raw_counts.csv').write_text('Gene,AAACCTGAGGCTACGA-1\nA1BG,0\nA2M,2\n',
                                            encoding='utf-8')
        explicit = row(gse, 'GSM1', 'raw_counts.csv')
        explicit['feature_column'] = '__row_names__'
        _, manifest = self.setup_manifest(gse, [explicit])
        run(manifest, self.root)
        completed = read_csv(manifest)[0]
        self.assertEqual(completed['feature_column'], '__row_names__')
        self.assertEqual(completed['text_header_missing_id'], 'false')

    def test_gse166504_dimensions_stream_despite_small_compressed_size(self):
        reader, estimated, reason = select_text_reader(
            25127, 82168, 120 * 1024**2, 32 * 1024**3, 'comma', 'genes_by_cells', [])
        self.assertEqual(reader, 'streaming')
        self.assertEqual(estimated, 25127 * 82168 * 8)
        self.assertIn('25% of physical RAM', reason)

    def _tenx(self, gse, features, matrix_dimensions=(2, 2, 2)):
        trio = self.root / 'data' / gse / 'raw' / 'trio'
        trio.mkdir(parents=True)
        rows, columns, nonzero = matrix_dimensions
        (trio / 'matrix.mtx').write_text(
            f'%%MatrixMarket matrix coordinate integer general\n% fixture\n{rows} {columns} {nonzero}\n1 1 1\n2 2 1\n',
            encoding='utf-8')
        (trio / 'features.tsv').write_text('\n'.join(features) + '\n', encoding='utf-8')
        (trio / 'barcodes.tsv').write_text('AAACCTGAGGCTACGA-1\nAAACCTGAGGTTCCTA-1\n', encoding='utf-8')

    def test_10x_missing_gene_name_records_feature_id_fallback(self):
        gse = 'GSE999999999'
        self._tenx(gse, ['ENSG000001\tNA\tGene Expression',
                         'ENSG000002\tA2M\tGene Expression'])
        workflow, manifest = self.setup_manifest(gse, [row(gse, 'GSM1', 'trio', '10x_mtx')])
        run(manifest, self.root)
        completed = read_csv(manifest)[0]
        self.assertEqual(completed['feature_name_fallback'], 'feature_id_for_missing_gene_name')
        report = json.loads((workflow / 'tenx_structure_probe.json').read_text())[0]
        self.assertEqual(report['status'], 'VERIFIED')
        self.assertEqual(report['missing_gene_names'], 1)
        self.assertEqual(report['gene_expression_features'], 2)

    def test_10x_dimension_mismatch_stops_before_build(self):
        gse = 'GSE999999999'
        self._tenx(gse, ['ENSG000001\tA1BG\tGene Expression',
                         'ENSG000002\tA2M\tGene Expression'], matrix_dimensions=(3, 2, 2))
        workflow, manifest = self.setup_manifest(gse, [row(gse, 'GSM1', 'trio', '10x_mtx')])
        with self.assertRaisesRegex(ValueError, 'dimensions mismatch'):
            run(manifest, self.root)
        report = json.loads((workflow / 'tenx_structure_probe.json').read_text())[0]
        self.assertEqual(report['status'], 'STOP')

    def _h5(self, gse, suffixes):
        raw = self.root / 'data' / gse / 'raw'
        raw.mkdir(parents=True)
        with h5py.File(raw / 'pooled.h5', 'w') as handle:
            group = handle.create_group('matrix')
            group.create_dataset('barcodes', data=[f'AAACCTGAGGCTACGA-{suffix}'.encode()
                                                  for suffix in suffixes])

    def _evidence(self, workflow, candidate=None):
        source = 'https://example.org/geo/cell_metadata.csv'
        evidence = {'checks': [
            {'source_type': 'geo_supplementary', 'url': source, 'result': 'candidate found'},
            {'source_type': 'publication_supplement', 'url': 'https://example.org/paper/supplement',
             'result': 'no matching map'},
            {'source_type': 'author_repository', 'reason': 'No public repository link in article',
             'result': 'not_found'},
        ], 'candidate_maps': [candidate] if candidate else []}
        (workflow / 'mapping_evidence.json').write_text(json.dumps(evidence), encoding='utf-8')
        return source

    def test_exhaustive_pooled_mapping_passes(self):
        gse = 'GSE999999999'
        self._h5(gse, [1, 2])
        rows = [row(gse, 'GSM1', 'pooled.h5', '10x_h5'),
                row(gse, 'GSM2', 'pooled.h5', '10x_h5')]
        workflow, manifest = self.setup_manifest(gse, rows)
        candidate = workflow / 'mapping_candidates' / 'cells.csv'
        candidate.parent.mkdir()
        write_csv(candidate, [
            {'cell': 'AAACCTGAGGCTACGA-1', 'sample': 'GSM1'},
            {'cell': 'AAACCTGAGGCTACGA-2', 'sample': 'GSM2'},
        ])
        self._evidence(workflow, dict(source_type='geo_supplementary',
            source_url='https://example.org/geo/cell_metadata.csv',
            local_path=f'data/{gse}/.workflow/mapping_candidates/cells.csv'))
        run(manifest, self.root)
        completed = read_csv(manifest)
        self.assertEqual(len({item['cell_map_path'] for item in completed}), 1)
        self.assertTrue((self.root / completed[0]['cell_map_path']).is_file())
        self.assertEqual(json.loads((workflow / 'pooled_mapping_probe.json').read_text())['status'], 'VERIFIED')
        verify_confirmation(manifest)

    def test_confirmed_subset_keeps_complete_pooled_map(self):
        gse = 'GSE999999999'
        self._h5(gse, [1, 2, 3])
        rows = [row(gse, 'GSM1', 'pooled.h5', '10x_h5'),
                row(gse, 'GSM2', 'pooled.h5', '10x_h5')]
        workflow, manifest = self.setup_manifest(gse, rows)
        write_csv(workflow / 'sample_report.csv',
                  [{'sample': sample} for sample in ('GSM1', 'GSM2', 'GSM3')])
        receipt_path = workflow / 'sample_manifest.confirmation.json'
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        receipt['selected_samples'] = ['GSM1', 'GSM2']
        receipt_path.write_text(json.dumps(receipt), encoding='utf-8')
        candidate = workflow / 'mapping_candidates' / 'cells.csv'
        candidate.parent.mkdir()
        write_csv(candidate, [
            {'cell': f'AAACCTGAGGCTACGA-{index}', 'sample': f'GSM{index}'}
            for index in (1, 2, 3)
        ])
        self._evidence(workflow, dict(source_type='geo_supplementary',
            source_url='https://example.org/geo/cell_metadata.csv',
            local_path=f'data/{gse}/.workflow/mapping_candidates/cells.csv'))
        run(manifest, self.root)
        result = json.loads((workflow / 'pooled_mapping_probe.json').read_text())
        self.assertEqual(result['status'], 'VERIFIED')
        self.assertEqual(result['excluded_samples'], ['GSM3'])
        self.assertEqual(result['excluded_cells'], 1)
        self.assertEqual(len(read_csv(self.root / result['cell_map_path'])), 3)
        verify_confirmation(manifest)

    def test_gse264203_seven_suffixes_six_gsm_stops(self):
        gse = 'GSE264203'
        self._h5(gse, range(1, 8))
        rows = [row(gse, f'GSM{8213939 + index}', 'pooled.h5', '10x_h5')
                for index in range(6)]
        workflow, manifest = self.setup_manifest(gse, rows)
        self._evidence(workflow)
        with self.assertRaisesRegex(ValueError, r'Barcode groups \(7\) and GSM count \(6\) differ'):
            run(manifest, self.root)
        report = json.loads((workflow / 'pooled_mapping_probe.json').read_text())
        self.assertEqual(len(report['barcode_suffix_counts']), 7)
        self.assertEqual(report['status'], 'STOP')
        self.assertFalse((self.root / 'data' / gse / 'seurat_raw.rds').exists())


if __name__ == '__main__':
    unittest.main()
