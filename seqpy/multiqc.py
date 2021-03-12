#! /usr/bin/env python
""" Parser for multiQC data
"""
from __future__ import annotations
from dataclasses import dataclass, field
from collections import namedtuple, defaultdict, namedtuple
from typing import Union
import configparser
import argparse
from pathlib import Path
import gzip
import json
from statistics import stdev, median, mean

@dataclass
class Subset:
    name: str
    datakey: str
    samples: list = field(repr=False)
    values: Union[list, dict] = field(repr=False)

class MultiQC():
    OUTLIER_COMPARISION = {"median":median,
                           "mean":mean}

    def __init__(self, multiQC_json: Path,
                       samples: list[str],
                       file_mapping_substrings: dict[str, str] = {"_R1_":"forward", "_R2":"reverse"},
                       forward_substring: str = "_R1_",
                       reverse_substring: str = "_R2_",
                       outlier_comparision_point: str = "median"):
        try:
            self.outlier_comparision = self.OUTLIER_COMPARISION[outlier_comparision_point]
        except KeyError:
            raise ValueError(f"Outlier comparision point not defined.  Select from {list(self.OUTLIER_COMPARISION.keys())}")
        self.samples = samples
        # substrings used to identify forward and reverse reads by file name
        # Only forward substring is used if single end study is used
        self._file_mapping_substrings = file_mapping_substrings

        # extracts data from multiQC json file.
        self.data = self._extract_multiQC_data(json_file = multiQC_json, samples = samples)
        self.sample_wise_data_keys = list(self.data[samples[0]].keys())

        # holds subsets computed by user
        self.subsets = defaultdict(dict)

    def compile_subset(self, samples_subset, subset_name, key):
        aggregate = self._compile_across_samples(data_mapping = self.data,
                                                 key = key,
                                                 samples_subset = samples_subset,
                                                subset_name = subset_name)
        self.subsets[subset_name][key] = Subset(name = subset_name,
                                                datakey = key,
                                                samples = samples_subset,
                                                values = aggregate)
        print(f"Compiled subset data for '{subset_name}' with samples {samples_subset} for {key}")
        return self.subsets[subset_name][key]

    def _compile_across_samples(self, data_mapping, key, samples_subset, subset_name):
        """ Creates an entry for a subset of samples for value in key.

        For examples, for a set of samples, sample_1 to sample_N.
        A subset of samples can be aggregated.
        Examples:
            # all samples aggregation for key1:
            data_mapping = _compile_across_samples(data_mapping,
                                                   key = "key1",
                                                   samples_subset = [sample_1,... ,sample_N],
                                                   subset_name = "all")
            # a subset of samples aggregation for key1:
            data_mapping = _compile_across_samples(data_mapping,
                                                   key = "key1",
                                                   samples_subset = [sample_1, sample_{i+2}, ...],
                                                   subset_name = "odd")


        Automatically determines if the values are dict like {index: value}
        or a single value per sample and returns the corresponding type.
        """
        aggregate = None

        for sample in self.samples:
            data = data_mapping[sample][key]
            if type(data) == float:
                # initiate aggregate data type to match
                # if already initiated, this does not affect aggregate
                aggregate = list() if not aggregate else aggregate
                aggregate.append(data)
            elif type(data) == dict:
                # these must be {index:value} dicts of length 1

                # initiate aggregate data type to match
                # if already initiated, this does not affect aggregate
                aggregate = defaultdict(list) if not aggregate else aggregate
                for index, value in data.items():
                    aggregate[index].append(value)
            else:
                raise ValueError(f"For {key}, {type(data)} type for data is unexpected.  Aggregation not implemented.")

        return aggregate

    def detect_outliers(self, key: str, deviation: float, subset_name: str = None):
        if subset_name:
            try:
                subset = self.subsets[subset_name][key]
            except KeyError:
                raise KeyError(f"Subset {subset_name} has not been compiled yet for {key}")
        # compile for all samples
        else:
            # compile all subset
            subset = self.compile_subset(samples_subset = self.samples,
                                         subset_name = "All",
                                         key = key)
        print(f"Checking for {key} outliers in {subset}\n"
              f"\tOutlier Criteria: (Standard deviations > {deviation})")

        outliers = list()

        # handle one value per sample style data
        if type(subset.values) == list:
            _stdev = stdev(subset.values)
            _median = median(subset.values)
            if _stdev == 0:
                print("Standard Deviation equals zero for {key}! No outliers detected!")
                return outliers
            for i, value in enumerate(subset.values):
                stdevs_from_median = abs(value - _median) / _stdev
                if stdevs_from_median > deviation:
                    outliers.append((subset.samples[i], stdevs_from_median))

        elif type(subset.values) == defaultdict:
            for index, values in subset.values.items():
                _stdev = stdev(values)
                _median = median(values)
                if _stdev == 0:
                    print(f"Standard Deviation equals zero for {key}:{index}. No outliers values.")
                    continue
                for i, value in enumerate(values):
                    stdevs_from_median = abs(value - _median) / _stdev
                    if stdevs_from_median > deviation:
                        outliers.append((f"{subset.samples[i]}:{index}", stdevs_from_median))
        else:
            print(type(subset.values))
            raise ValueError("Unknown type for outlier detection")
        if outliers:
            print(f"Outliers detected!")
        else:
            print(f"No outliers detected!")
        return outliers


    def _label_file(self, filename: str):
        """ Given a filename.  Return the file label based on a unique substring.

        Substrings to file label mapping should be provided at init.
        """
        matched = False
        for substring, filelabel in self._file_mapping_substrings.items():
            if substring in filename:
                if not matched:
                    matched = filelabel
                else:
                    raise ValueError(f"File name {filename} matched multiple substrings in provided mapping {self._file_mapping_substrings}")
        if matched:
            return matched
        else:
            # no matches
            raise ValueError(f"File name {filename} did not match any substrings in provided mapping {self._file_mapping_substrings}")

    # data extraction functions
    # TODO: These should return a value that will be assigned directly to the data mapping
    def _extract_multiQC_data(self, json_file: Path, samples):
        data_mapping = defaultdict(lambda: defaultdict(dict))
        with open(json_file, "r") as f:
            raw_data = json.load(f)

        ###  extract general stats
        for file_data in raw_data["report_general_stats_data"]:
            for file_name, data in file_data.items():
                cur_sample = [sample for sample in samples if sample in file_name][0]
                filelabel = self._label_file(file_name)
                for key, value in data.items():
                    data_mapping[cur_sample][f"{filelabel}-{key}"] = value

        ### extract plot data

        # extract fastqc_sequence_counts_plot
        for plot_name, data in raw_data["report_plot_data"].items():
            # different kinds of plot data need to be extracted with different methods
            plot_type = data["plot_type"]
            if plot_type == "bar_graph":
                data_mapping = self._extract_from_bar_graph(data, plot_name, data_mapping, samples)
            elif plot_type == "xy_line":
                data_mapping = self._extract_from_xy_line_graph(data, plot_name, data_mapping, samples)
            else:
                raise ValueError(f"Unknown plot type {plot_type}. Data parsing not implemented for multiQC {plot_type}")

        return data_mapping

    def _extract_from_bar_graph(self, data, plot_name, data_mapping, samples):
        # determine data mapping for samples in multiqc (which are files)
        # and samples in a dataset (which may have a forward and reverse read file)
        mqc_samples_to_samples = dict()
        # this should be a list with one entry
        assert len(data["samples"]) == 1
        for i, mqc_sample in enumerate(data["samples"][0]):
            matching_samples = [sample for sample in samples if sample in mqc_sample]
            # only one sample should map
            assert len(matching_samples) == 1

            sample = matching_samples[0]
            mqc_samples_to_samples[i] = (sample, self._label_file(mqc_sample))

        # iterate through data from datasets
        # this should be a list with one entry
        assert len(data["datasets"]) == 1
        for sub_data in data["datasets"][0]:
            name = sub_data["name"]
            values = sub_data["data"]
            for i, value in enumerate(values):
                sample, sample_file = mqc_samples_to_samples[i]
                data_mapping[sample][f"{sample_file}-{plot_name}-{name}"] = values[i]
        return data_mapping

    def _extract_from_xy_line_graph(self, data, plot_name, data_mapping, samples):
        # determine data mapping for samples in multiqc (which are files)
        # and samples in a dataset (which may have a forward and reverse read file)

        # Iterate through datasets (each typically representing one sample)
        # Nested list of list, top level list includes percentage and counts

        # plots with multiple value types are detected
        multiple_value_types = True if len(data["datasets"]) > 1 else False

        # plots with bins are detected
        if "categories" in data["config"].keys():
            bins = [str(bin) for bin in data["config"]["categories"]]
        else:
            bins = False


        # dataset represents an entire plot (i.e. all lines)
        # Note: for xy plots with both percent and raw counts, there will be two datasets
        for i, dataset in enumerate(data["datasets"]):
            if multiple_value_types:
                data_label = f"-{data['config']['data_labels'][i]['name']}"
            else:
                data_label = ""

            # dataset entry represents one sample (i.e. one line from the line plot)
            for dataset_entry in dataset:
                file_name = dataset_entry["name"]

                # values is a list of [index,value]
                values = dataset_entry["data"]
                matching_samples = [sample for sample in samples if sample in file_name]
                # only one sample should map
                assert len(matching_samples) == 1

                sample = matching_samples[0]
                sample_file = self._label_file(file_name)

                # three level nested dict entries for xy graphs
                # {sample: {sample_file-plot_type: {index: value}}}
                data_key = f"{sample_file}-{plot_name}{data_label}"
                data_mapping[sample][data_key] = dict()
                cur_data_mapping_entry = data_mapping[sample][data_key]
                # for plots with bins, add bin string to values iterable
                if bins:
                    values = zip(bins, values)
                # for non-categorical bins, each values should be an [index,value]
                for j, value in values:
                    cur_data_mapping_entry[j] = value
        return data_mapping
