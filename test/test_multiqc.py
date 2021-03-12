from seqpy.multiqc import MultiQC

INPUT = "assets/test/multiqc_data.json"
SAMPLES = ['Mmus_BAL-TAL_LRTN_FLT_Rep3_F8', 'Mmus_BAL-TAL_RRTN_BSL_Rep3_B9', 'Mmus_BAL-TAL_LRTN_FLT_Rep4_F9',
           'Mmus_BAL-TAL_LRTN_BSL_Rep1_B7', 'Mmus_BAL-TAL_LRTN_FLT_Rep5_F10', 'Mmus_BAL-TAL_LRTN_FLT_Rep2_F7',
           'Mmus_BAL-TAL_LRTN_FLT_Rep1_F6', 'Mmus_BAL-TAL_RRTN_BSL_Rep4_B10', 'Mmus_BAL-TAL_LRTN_GC_Rep3_G9',
           'Mmus_BAL-TAL_LRTN_GC_Rep2_G8', 'Mmus_BAL-TAL_RRTN_BSL_Rep2_B8', 'Mmus_BAL-TAL_LRTN_GC_Rep1_G6',
           'Mmus_BAL-TAL_RRTN_GC_Rep4_G10']

EXPECTED_DATA_KEY = ['forward-percent_gc',
'forward-avg_sequence_length',
'forward-total_sequences',
'forward-percent_duplicates',
'forward-percent_fails',
'reverse-percent_gc',
'reverse-avg_sequence_length',
'reverse-total_sequences',
'reverse-percent_duplicates',
'reverse-percent_fails',
'forward-fastqc_sequence_counts_plot-Unique Reads',
'reverse-fastqc_sequence_counts_plot-Unique Reads',
'forward-fastqc_sequence_counts_plot-Duplicate Reads',
'reverse-fastqc_sequence_counts_plot-Duplicate Reads',
'forward-fastqc_per_base_sequence_quality_plot',
'reverse-fastqc_per_base_sequence_quality_plot',
'forward-fastqc_per_sequence_quality_scores_plot',
'reverse-fastqc_per_sequence_quality_scores_plot',
'forward-fastqc_per_sequence_gc_content_plot-Percentages',
'reverse-fastqc_per_sequence_gc_content_plot-Percentages',
'forward-fastqc_per_sequence_gc_content_plot-Counts',
'reverse-fastqc_per_sequence_gc_content_plot-Counts',
'forward-fastqc_per_base_n_content_plot',
'reverse-fastqc_per_base_n_content_plot',
'forward-fastqc_sequence_duplication_levels_plot',
'reverse-fastqc_sequence_duplication_levels_plot',
'forward-fastqc_overrepresented_sequencesi_plot-Top over-represented sequence',
'reverse-fastqc_overrepresented_sequencesi_plot-Top over-represented sequence',
'forward-fastqc_overrepresented_sequencesi_plot-Sum of remaining over-represented sequences',
'reverse-fastqc_overrepresented_sequencesi_plot-Sum of remaining over-represented sequences',
'forward-fastqc_adapter_content_plot',
'reverse-fastqc_adapter_content_plot',]


import unittest

class TestMultiQCParser(unittest.TestCase):
    mqc = MultiQC(multiQC_json = INPUT,
                  samples = SAMPLES,
                  outlier_comparision_point = "median")

    def test_data_key_loading(self):
        mqc = self.mqc
        self.assertEqual(mqc.sample_wise_data_keys, EXPECTED_DATA_KEY)

    def test_compilation_for_one_value_per_sample(self):
        mqc = self.mqc
        LRTN_subset = [sample for sample in mqc.samples if "LRTN" in sample]

        key = 'reverse-percent_gc'
        mqc.compile_subset(samples_subset = LRTN_subset,
                           subset_name = "Live Return",
                           key = key)

        self.assertEqual(mqc.subsets["Live Return"][key].values, [52.0, 51.0, 52.0, 53.0, 51.0, 52.0, 52.0, 50.0, 53.0, 53.0, 53.0, 53.0, 53.0])

    def test_compilation_for_indexed_values_per_sample(self):
        mqc = self.mqc
        FLT_subset = [sample for sample in mqc.samples if "FLT" in sample]

        key = 'forward-fastqc_per_base_n_content_plot'
        mqc.compile_subset(samples_subset = FLT_subset,
                           subset_name = "Fleet Group",
                           key = key)

        # check length matches expected
        self.assertEqual(len(mqc.subsets["Fleet Group"][key].values), 38)
        # check specific entry value matches expected by sum check
        self.assertEqual(sum(mqc.subsets["Fleet Group"][key].values[2]), 0.748698368731686)

    def test_outlier_detection(self):
        # check with stringent threshold for detection
        deviation = 0.5
        self.assertEqual(len(self.mqc.detect_outliers(key = 'forward-percent_duplicates', deviation = deviation)), 10)

        # check with more lax threshold
        deviation = 1
        self.assertEqual(len(self.mqc.detect_outliers(key = 'forward-percent_duplicates', deviation = deviation)), 5)

        # check with unreasonably high threshold
        deviation = 99999
        self.assertEqual(len(self.mqc.detect_outliers(key = 'forward-percent_duplicates', deviation = deviation)), 0)

    def test_bad_substring_mapping(self):
        with self.assertRaises(ValueError, msg="Bad file substring mapping failed as intended"):        
            mqc = MultiQC(multiQC_json = INPUT,
                          samples = SAMPLES,
                          outlier_comparision_point = "median",
                          file_mapping_substrings = {"UNREAL":"File_X"})
