from seqpy.multiqc import MultiQC

RAW_INPUT = "assets/test/raw_multiqc_data.json"
TRIM_INPUT = "assets/test/trimmed_multiqc_data.json"
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
    mqc = MultiQC(multiQC_json = RAW_INPUT,
                  samples = SAMPLES,
                  outlier_comparision_point = "median")

    def test_data_key_loading(self):
        mqc = self.mqc
        self.assertEqual(sorted(mqc.sample_wise_data_keys), sorted(EXPECTED_DATA_KEY))

    def test_data_access_for_general_stats_data(self):
        mqc = self.mqc
        key = 'reverse-percent_gc'
        test = mqc.data[SAMPLES[0]][key]
        print(test, type(test))
        self.assertEqual(mqc.data[SAMPLES[0]][key].value, 52)
        self.assertEqual(mqc.data[SAMPLES[0]][key].units, 'percent_gc')

    def test_data_access_for_bar_graph_data(self):
        mqc = self.mqc
        key = 'forward-fastqc_sequence_counts_plot-Unique Reads'
        test = mqc.data[SAMPLES[0]][key]
        print(test, type(test))
        self.assertEqual(mqc.data[SAMPLES[0]][key].value, 5420)
        self.assertEqual(mqc.data[SAMPLES[0]][key].units, "Number of reads")

    def test_data_access_for_xy_line_graph_data(self):
        mqc = self.mqc
        key = 'forward-fastqc_per_sequence_quality_scores_plot'
        test = mqc.data[SAMPLES[0]][key]
        print(test, type(test))
        # check count of sequences in forward reads for sample 1
        # with a Mean Sequence Quality (Phred Score) of 40 is equal to 4071
        self.assertEqual(mqc.data[SAMPLES[0]][key].units, 'Count')
        self.assertEqual(mqc.data[SAMPLES[0]][key].bin_units, 'Mean Sequence Quality (Phred Score)')
        self.assertEqual(mqc.data[SAMPLES[0]][key].datakey, key)
        self.assertEqual(mqc.data[SAMPLES[0]][key].values[40.0], 4071.0)

    def test_compilation_for_one_value_per_sample(self):
        mqc = self.mqc
        LRTN_subset = [sample for sample in mqc.samples if "LRTN" in sample]

        key = 'reverse-percent_gc'
        LRTN_reverse_gc_content = mqc.compile_subset(samples_subset = LRTN_subset,
                                                     key = key)

        print(LRTN_subset)
        self.assertEqual(LRTN_reverse_gc_content, [52.0, 52.0, 53.0, 51.0, 52.0, 52.0, 53.0, 53.0, 53.0])

    def test_compilation_for_indexed_values_per_sample(self):
        mqc = self.mqc
        FLT_subset = [sample for sample in mqc.samples if "FLT" in sample]

        key = 'reverse-fastqc_per_base_n_content_plot'
        fleet_reverse_fastqcByBase_n_content = mqc.compile_subset(samples_subset = FLT_subset,
                                                                  key = key)

        # check number of bins is correct
        self.assertEqual(len(fleet_reverse_fastqcByBase_n_content), 38)
        # check number fleet samples is correct
        self.assertEqual(len(fleet_reverse_fastqcByBase_n_content[1]), 5)
        # check specific entry value matches expected by difference of expected sum check
        test_metric = abs(sum(fleet_reverse_fastqcByBase_n_content[1]) - 6.61705047919786)
        self.assertLess(test_metric, 0.0001)

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
            mqc = MultiQC(multiQC_json = RAW_INPUT,
                          samples = SAMPLES,
                          outlier_comparision_point = "median",
                          file_mapping_substrings = {"UNREAL":"File_X"})

class TestTrimmedMultiQCParser(TestMultiQCParser):
    mqc = MultiQC(multiQC_json = TRIM_INPUT,
                  samples = SAMPLES,
                  outlier_comparision_point = "median")

    def test_data_key_loading(self):
        mqc = self.mqc
        expected_data_keys = EXPECTED_DATA_KEY.copy()
        expected_data_keys.remove('reverse-fastqc_adapter_content_plot')
        expected_data_keys.remove('forward-fastqc_adapter_content_plot')
        expected_data_keys.append('reverse-fastqc_sequence_length_distribution_plot')
        expected_data_keys.append('forward-fastqc_sequence_length_distribution_plot')

        [print(key) for key in mqc.sample_wise_data_keys]
        print(set(mqc.sample_wise_data_keys).difference(set(expected_data_keys)))
        self.assertTrue(set(mqc.sample_wise_data_keys) == set(expected_data_keys))

    def test_data_access_for_bar_graph_data(self):
        mqc = self.mqc
        key = 'forward-fastqc_sequence_counts_plot-Unique Reads'
        test = mqc.data[SAMPLES[0]][key]
        print(test, type(test))
        self.assertEqual(mqc.data[SAMPLES[0]][key].value, 5421)
        self.assertEqual(mqc.data[SAMPLES[0]][key].units, "Number of reads")

    def test_outlier_detection(self):
        # check with stringent threshold for detection
        deviation = 0.5
        self.assertEqual(len(self.mqc.detect_outliers(key = 'forward-percent_duplicates', deviation = deviation)), 9)

        # check with more lax threshold
        deviation = 1
        self.assertEqual(len(self.mqc.detect_outliers(key = 'forward-percent_duplicates', deviation = deviation)), 4)

        # check with unreasonably high threshold
        deviation = 99999
        self.assertEqual(len(self.mqc.detect_outliers(key = 'forward-percent_duplicates', deviation = deviation)), 0)

    def test_outlier_detection_on_indexed_values(self):
        key = "forward-fastqc_sequence_length_distribution_plot"
        # check with stringent threshold for detection
        deviation = 0.5
        self.assertEqual(len(self.mqc.detect_outliers(key = key, deviation = deviation)), 196)

        # check with more lax threshold
        deviation = 1
        self.assertEqual(len(self.mqc.detect_outliers(key = key, deviation = deviation)), 109)

        # check with unreasonably high threshold
        deviation = 99999
        self.assertEqual(len(self.mqc.detect_outliers(key = key, deviation = deviation)), 0)

    def test_data_access_for_xy_line_graph_data(self):
        mqc = self.mqc
        key = 'forward-fastqc_per_sequence_quality_scores_plot'
        test = mqc.data[SAMPLES[0]][key]
        print(test, type(test))
        # check count of sequences in forward reads for sample 1
        # with a Mean Sequence Quality (Phred Score) of 40 is equal to 4071
        self.assertEqual(mqc.data[SAMPLES[0]][key].units, 'Count')
        self.assertEqual(mqc.data[SAMPLES[0]][key].bin_units, 'Mean Sequence Quality (Phred Score)')
        self.assertEqual(mqc.data[SAMPLES[0]][key].datakey, key)
        self.assertEqual(mqc.data[SAMPLES[0]][key].values[40.0], 4221.0)

    def test_compilation_for_one_value_per_sample(self):
        mqc = self.mqc
        LRTN_subset = [sample for sample in mqc.samples if "LRTN" in sample]

        key = 'reverse-percent_gc'
        LRTN_reverse_gc_content = mqc.compile_subset(samples_subset = LRTN_subset,
                                                     key = key)

        print(LRTN_subset)
        self.assertEqual(LRTN_reverse_gc_content, [52.0, 52.0, 53.0, 51.0, 52.0, 53.0, 53.0, 54.0, 53.0])

    def test_compilation_for_indexed_values_per_sample(self):
        mqc = self.mqc
        FLT_subset = [sample for sample in mqc.samples if "FLT" in sample]

        key = 'reverse-fastqc_per_base_n_content_plot'
        fleet_reverse_fastqcByBase_n_content = mqc.compile_subset(samples_subset = FLT_subset,
                                                                  key = key)

        # check number of bins is correct
        self.assertEqual(len(fleet_reverse_fastqcByBase_n_content), 38)
        # check number fleet samples is correct
        self.assertEqual(len(fleet_reverse_fastqcByBase_n_content[1]), 5)
        # check specific entry value matches expected by difference of expected sum check
        test_metric = abs(sum(fleet_reverse_fastqcByBase_n_content[1]) - 6.593416959915156)
        self.assertLess(test_metric, 0.0001)
