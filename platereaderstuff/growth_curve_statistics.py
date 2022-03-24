#!/usr/bin/env python3
DESCRIPTION = """
Author: Cody Martin
        University of Wisconsin-Madison
        Department of Bacteriology
        Anantharaman lab

Purpose: Process plate reader growth curve data under different 
         strain/media/concentration conditions. Additionally,
         generate summary statistics. Plate reader = 
         Hryckowian plate reader (Epoch 2, using Gen5 v3.10.06 software)  

Usage: ./growth_curve_statistics.py -i <input.csv> -k <plate_setup.csv> -o <output_basename> [OPTIONS]
"""

import argparse
import pandas as pd
import numpy as np
from functools import reduce


class WellPlate:
    def __init__(self, file: str, setup: str, noblank: bool):
        self.file = file
        self.tmpdata = self.read_data(self.file)
        self.update_time()

        self.setup = self.parse_key(setup)

        self.data = self.process()

        self.tech_reps = self.split_technical_reps()

        # TODO: add blanking later
        if not noblank:
            pass

    def read_data(self, file: str) -> pd.DataFrame:
        """Read data file as dataframe

        Args:
            file (str): input .csv file with plate reader data

        Returns:
            pd.DataFrame
        """
        return pd.read_csv(file)

    def read_setup(self, file: str) -> pd.DataFrame:
        """Read plate setup file

        Args:
            file (str): setup file in .csv format
                NOTE: must have the Row and column labels
                and the Row column should be called 'Row'

        Returns:
            pd.DataFrame
        """
        return pd.read_csv(file, index_col="Row")

    def parse_key(self, file: str) -> pd.DataFrame:
        """Parse the plate setup dataframe to return a cleaner
        dataframe that has the well coordinates for each sample

        Args:
            file (str): setup file in .csv format

        Returns:
            pd.DataFrame: data frame with the well coordinates
                of each sample on the plate
                NOTE: this will only parse samples included in
                the setup file, so if there are samples that
                should not be processed from a plate for whatever
                reason, do not include them in the plate setup.
        """
        setup = self.read_setup(file)
        coordinates = list()

        # only parse cells that are not empty/blank
        for col, rowseries in setup.notna().iteritems():
            series = setup.loc[rowseries, col]
            if not series.empty:
                wells = list(map(lambda x: f"{x}{series.name}", series.index.to_list()))
                coordinates.append(list(zip(series.to_list(), wells)))

        sample_coordinates = reduce(lambda x, y: x + y, coordinates)

        return pd.DataFrame(sample_coordinates, columns=["Sample", "Well"])

    def update_time(self):
        """The Gen5 software likes to change the format
        of the time when the run goes for longer than 24 hours.
        The normal format is mm:ss but after 24 h becomes hh:mm:ss.

        This standardizes the times, so that they can be made relative
        to the first time point at t=0.
        """
        time = pd.to_timedelta(
            np.where(
                self.tmpdata.Time.str.count(":") == 1,
                self.tmpdata.Time + ":00",
                self.tmpdata.Time,
            )
        ) / pd.Timedelta(hours=1)

        self.tmpdata.Time = time

    def combine_metadata(
        self, tmpdata: pd.DataFrame, setup: pd.DataFrame
    ) -> pd.DataFrame:
        """Combine the metadata from the parsed setup dataframe with the raw
        data. In other words, for each sample in the parsed setup dataframe,
        add the OD values.

        Args:
            tmpdata (pd.DataFrame): raw data from plate reader in .csv format
            setup (pd.DataFrame): setup dataframe from self.parse_setup()

        Returns:
            pd.DataFrame: data frame with sample name, well coordinates, and OD
                at all time points
        """
        return tmpdata.melt(id_vars="Time", var_name="Well", value_name="OD").merge(
            setup, on="Well", how="right"
        )

    # TODO: add optional blanking
    def process(self) -> pd.DataFrame:
        """For a single plate experiment, combine all the technical replicates
        and summarize them with the mean and standard deviation.

        Returns:
            pd.DataFrame: summarized dataframe collapsing all technical replicates
                for a given sample into only the mean and standard deviation
        """
        data = self.combine_metadata(self.tmpdata, self.setup)
        processed_data = (
            data.groupby(["Time", "Sample"])
            .OD.agg(["mean", "std"])
            .reset_index()
            .rename({"mean": "OD"}, axis=1)
            .sort_values(by=["Sample", "Time"])
            .reset_index(drop=True)
        )

        processed_data[
            ["Strain", "Media", "Concentration"]
        ] = processed_data.Sample.str.split("_", expand=True)

        return processed_data.astype(
            {
                "Sample": "category",
                "Strain": "category",
                "Media": "category",
                "Concentration": "category",
            }
        )

    def split_technical_reps(self) -> list[pd.DataFrame]:
        """Should you want to consider each technical replicate as
        independent replicates, this will split up the combined metadata
        output to split all the technical replicates into separate dataframes.

        This allows the GrowthCurveAnalyzer class to take the same inputs as
        with processing the data in terms of collapsing the technical replicates
        into a single value.

        NOTE: The way the splitting occurs requires that the same number of technical
        replicates for each sample be used. In other words, for a single plate, all
        samples should have 3 technical replicates for example (can be some other number
        as long as all samples have the same number).

        Returns:
            list[pd.DataFrame]: list of dataframes where all technical replicates
                for a given sample are split equally into each dataframe such that
                no single dataframe contains the same sample more than once
        """
        groups = self.setup.groupby("Sample")
        unique_rep_counts = groups.size().unique()

        ## Note: need to have same number of technical replicates per plate
        n_reps = unique_rep_counts[0]
        technical_reps_setup = [groups.nth(i).reset_index() for i in range(n_reps)]

        time = self.tmpdata.Time
        technical_reps = [
            pd.concat([time, self.tmpdata[setup.Well]], axis=1)
            for setup in technical_reps_setup
        ]

        return [
            self.combine_metadata(trep, setup)
            for trep, setup in zip(technical_reps, technical_reps_setup)
        ]

    def combine_technical_reps(self) -> pd.DataFrame:
        """Row stack all the technical replicates back into a single dataframe.

        Returns:
            pd.DataFrame: data frame with all technical replicates collapsed into
                one dataframe. Should be analogous to self.combine_metadata()
        """
        return reduce(lambda x, y: pd.concat([x, y]), self.tech_reps).reset_index(
            drop=True
        )

    def prism_technical_reps(self, droplog=False) -> pd.DataFrame:
        """Pivots the dataframe into a wide format that is more compatible with
        use in GraphPad Prism

        Args:
            droplog (bool, optional): drop the logOD column if it exists.
                Defaults to False.

        Returns:
            pd.DataFrame: wide data frame for exporting to prism
                Each sample is additionally labeled by its well
        """
        combined = (
            self.combine_technical_reps()
            .assign(Sample_well=lambda x: x.Sample + "_" + x.Well)
            .drop(["Well", "Sample"], axis=1)
            .pivot(index="Time", columns="Sample_well")
        )

        if droplog:
            combined = combined.drop(columns=["logOD"])

        # flatten multiindex column naming
        columns = [name[1] for name in combined.columns.to_flat_index().to_list()]
        combined.columns = columns
        return combined[sorted(combined.columns)].reset_index()

    def prism_output(self) -> pd.DataFrame:
        """For the processed data that collapses technical replicates
        into their mean/std, convert dataframe into wide format that
        is more suitable to plotting in GraphPad Prism.

        Returns:
            pd.DataFrame: wide dataframe where each sample occupies two
                columns, one for the mean OD and one for the stdev
                at each time point
        """
        data = self.data.drop(["Strain", "Media", "Concentration"], axis=1).pivot(
            index="Time", columns="Sample"
        )

        columns = [
            f"{sample}_{val}" for val, sample in data.columns.to_flat_index().to_list()
        ]
        data.columns = columns
        return data[sorted(data.columns)].reset_index()


class GrowthCurveAnalyzer:
    def __init__(
        self,
        file: str,
        setup: str,
        output: str,
        noblank: bool,
        prism: bool,
    ):
        self.file = file
        self.setup = setup
        self.baseoutput = output.rsplit(".", 1)[0]

        self.plate = WellPlate(file, setup, noblank)

        # if users want the prism output of the data to plot in GraphPad Prism
        if prism:
            self.plate.prism_output().to_csv(
                f"{self.baseoutput}_processed_data_PRISM.csv", index=False
            )
        else:
            self.plate.data.to_csv(f"{self.baseoutput}_processed_data.csv", index=False)

        # add natural log of OD for subsequent processing
        self.plate.data["logOD"] = np.log(self.plate.data.OD)

    def smooth_data(self, data: pd.DataFrame, column: str, intervals=5) -> pd.DataFrame:
        """Smoothen a given column of data over n time intervals for all samples
        using a rolling average over the n intervals. Each sample has its own data
        smoothed independently of each other.

        Args:
            data (pd.DataFrame): processed data from the WellPlate class
            column (str): column to smooth
            intervals (int, optional): number of time intervals to smooth data over.
                Defaults to 5.

        Returns:
            pd.DataFrame: copy of input data where column values for indicate
                column have been smoothed. Should be otherwise identical to the
                input data.
        """
        return (
            data.groupby("Sample")[column]
            .rolling(intervals, min_periods=1)
            .mean()
            .reset_index()
            .drop("level_1", axis=1)
        )

    def get_max_column_value(
        self, data: pd.DataFrame, column: str, intervals=5
    ) -> pd.DataFrame:
        """For a given column, calculate the max value AFTER smoothing.

        Args:
            data (pd.DataFrame): processed data from the WellPlate class
            column (str): column to grab max value from
            intervals (int, optional): number of time intervals to smooth data over.
                Defaults to 5.

        Returns:
            pd.DataFrame: dataframe of max values of input column for ALL samples
        """
        smoothed = self.smooth_data(data, column, intervals)
        argmax = smoothed.groupby("Sample").idxmax()[column]
        return data.loc[argmax.reset_index(drop=True)].copy()

    def growth_rate(self, data: pd.DataFrame) -> pd.DataFrame:
        """Calculate the growth rate for each sample as the differential
        of the logOD over time for adjacent times. No smoothing applied here.

        The indices in the growth rates vectors are -1 compared to the actual
        sample/time indices, meaning that the index from the growth rate
        corresponds to the time interval of (index, index + 1) states. This
        time interval range will be included in the output.

        NOTE: the units of growth rate are in lnOD/h. You cannot simply
        do exp(lnOD/h) to get the data in OD/h. Best would be to look
        at the time interval and choose the raw data points to do
        that calculation.

        Args:
            data (pd.DataFrame): processed data from the WellPlate class

        Returns:
            pd.DataFrame: dataframe the bounds of the time interval and
                corresponding growth rate
        """

        g_rates = (
            data.groupby("Sample")
            .apply(lambda x: np.diff(x.logOD) / np.diff(x.Time))
            .reset_index()
            .rename({0: "growth_rate"}, axis=1)
        )

        # create a tuple of the time intervals
        time_intervals = list()
        for t in data.Time.unique():
            try:
                time_intervals.append((previous, t))
            except UnboundLocalError:
                pass  # for the first one
            previous: int = t

        time_intervals *= len(data.Sample.unique())
        g_rates_df = pd.DataFrame(
            [
                [name, value]
                for name, growth_rate_vector in zip(g_rates.Sample, g_rates.growth_rate)
                for value in list(growth_rate_vector)
            ],
            columns=["Sample", "growth_rate"],
        )

        return pd.concat(
            [pd.DataFrame(time_intervals, columns=["t1", "t2"]), g_rates_df], axis=1
        )

    def max_growth_rate(self, data: pd.DataFrame, intervals=5) -> pd.DataFrame:
        """Calculate the max growth rate for each sample with the growth
        rate smoothed over n time intervals in units of lnOD/h.

        Args:
            data (pd.DataFrame): processed data from the WellPlate class
            intervals (int, optional): number of time intervals to smooth data over.
                Defaults to 5.

        Returns:
            pd.DataFrame: dataframe with the max growth rate for each sample
                and the time interval for which the max growth rate occurs
        """
        # calculate growth rates and then smooth them
        growth_rates = (
            self.growth_rate(data)
            .groupby("Sample")
            .rolling(intervals, min_periods=1)
            .mean()
            .reset_index()
            .drop("level_1", axis=1)
        )

        # get the index of the max growth rate for each sample
        argmax = growth_rates.groupby("Sample").idxmax().growth_rate

        # return the actual max growth rate values
        return growth_rates.loc[argmax.reset_index(drop=True)]

    def lag_time(self, data: pd.DataFrame, intervals=5) -> pd.DataFrame:
        """Calculate lag time as defined as half time to the max
        growth rate. Since the time at the max growth rate corresponds
        to an interval, the midpoint of the interval is used as the
        time at the max growth rate, then divided by two.

        Data are smoothed over n time intervals.

        Args:
            data (pd.DataFrame): processed data from the WellPlate class
            intervals (int, optional): number of time intervals to smooth data over.
                Defaults to 5.

        Returns:
            pd.DataFrame: dataframe with lagtime in units of hrs for each sample
        """
        # use midpoint of time interval to set a single time as the max
        max_growth_rates = self.max_growth_rate(data, intervals).assign(
            lag_time=lambda x: (x.t1 + (x.t2 - x.t1) / 2) / 2
        )
        return max_growth_rates

    def doubling_time(self, data: pd.DataFrame, intervals=5) -> pd.DataFrame:
        """Calculate doubling time as defined as ln(2) / max_growth_rate for
        each sample. Data are smoothed over n time intervals.

        Args:
            data (pd.DataFrame): processed data from the WellPlate class
            intervals (int, optional): number of time intervals to smooth data over.
                Defaults to 5.

        Returns:
            pd.DataFrame: dataframe with doubling time in units of hours
                for each sample
        """
        log2 = np.log(2)
        return self.max_growth_rate(data, intervals).assign(
            doubling_time=lambda x: log2 / x.growth_rate
        )

    def summary(
        self, data: pd.DataFrame, save: bool, intervals=5, technical_rep=False
    ) -> pd.DataFrame:
        """Calculate the following summary statistics for each sample
        and output them to a dataframe to save:
            - max OD
            - max lnOD
            - final OD
            - final lnOD
            - max growth rate
            - lag time
            - doubling time

        See README.txt file that is generated for more explanation of units and how
        values are calculated.

        Args:
            data (pd.DataFrame): processed data from the WellPlate class
            save (bool): _description_
            intervals (int, optional): number of time intervals to smooth data over.
                Defaults to 5.
            technical_rep (bool, optional): if inputting a dataframe of only a single technical
                replicate, there will be no standard deviation column, so this prevents trying
                to drop that column. Defaults to False.

        Returns:
            pd.DataFrame: _description_
        """
        max_ODs = self.get_max_column_value(data, "OD").rename(
            {"OD": "max_OD", "logOD": "max_logOD"}, axis=1
        )
        final_ODs = (
            self.smooth_data(data, "OD", intervals)
            .groupby("Sample")
            .nth(-1)
            .reset_index()
            .rename({"OD": "final_OD"}, axis=1)
            .assign(final_logOD=lambda x: np.log(x.final_OD))
        )
        lagtime = self.lag_time(data, intervals).rename(
            {"growth_rate": "max_growth_rate"}, axis=1
        )
        doublingtime = self.doubling_time(data, intervals)[["Sample", "doubling_time"]]

        summary_df = reduce(
            lambda left, right: pd.merge(left, right, on="Sample"),
            (max_ODs, final_ODs, lagtime, doublingtime),
        )

        summary_df = summary_df if technical_rep else summary_df.drop("std", axis=1)

        if save:
            summary_df.to_csv(f"{self.baseoutput}_statistical_summary.csv", index=False)
        return summary_df

    def summarize_technical_replicates(self):
        """If wanting to process individual technical replicates instead
        of collapsing them into the mean of technical replicates, this
        will setup passing the list of individual technical replicates
        from WellPlate.tech_reps into the GrowthCurveAnalyzer.summary
        function.
        """
        for trep in self.plate.tech_reps:
            trep["logOD"] = np.log(trep.OD)

        trep_summaries = [
            self.summary(data=group, save=False, technical_rep=True)
            for group in self.plate.tech_reps
        ]
        trep_summary = (
            reduce(lambda x, y: pd.concat([x, y]), trep_summaries)
            .sort_values(by=["Sample", "Well"])
            .reset_index(drop=True)
        )

        trep_summary.to_csv(
            f"{self.baseoutput}_statistical_summary_TECHNICAL_REPS.csv", index=False
        )
        self.plate.prism_technical_reps(droplog=True).to_csv(
            f"{self.baseoutput}_TECHNICAL_REPLICATES.csv", index=False
        )

    def write_readme(self, name: str = "README.txt"):
        """Write a metadata file that has descriptions of the columns in each returned dataset."""

        file_descriptions = [
            (
                "_processed_data.csv",
                "Biological replicate data processed to average technical reps together per sample",
            ),
            (
                "_processed_data_PRISM.csv",
                "Same as _processed_data.csv except formated for use in GraphPad Prism",
            ),
            (
                "_TECHNICAL_REPLICATES.csv",
                "Relabeling columns from input data to the format Sample_Well",
            ),
            (
                "_statistical_summary.csv",
                "Growth curve statistical summary for a single biological replicate with values averaged across technical replicates",
            ),
            (
                "_statistical_summary_TECHNICAL_REPLICATES.csv",
                "Growth curve statistical summary for each technical replicate individually",
            ),
        ]
        stat_summary_metadata = [
            ("Time", "time at max OD or max logOD"),
            ("Sample", "Full sample name"),
            ("max_OD", "max OD value"),
            ("Strain", "Strain, 1st part of Sample"),
            ("Media", "Media, 2nd part of Sample"),
            ("Concentration", "Concentration with units, 3rd part of Sample"),
            ("max_logOD", "max log OD value ** Note: this is natural log"),
            ("final_OD", "final OD value smoothed over the last 5 intervals"),
            ("final_logOD", "final lnOD value smoothed over the last 5 intervals"),
            ("t1", "left side of time interval for max growth rate"),
            ("t2", "right side of time interval for max growth rate"),
            (
                "max_growth_rate",
                "maximum growth rate smoothed over 5 time intervals in units of lnOD/h",
            ),
            ("lag_time", "culture lag time in hours"),
            ("doubling_time", "culture doubling time in units of hours"),
        ]

        with open(name, "w") as readme:
            readme.write("----File descriptions----\n")
            for file, desc in file_descriptions:
                readme.write(f"{file:46s}{desc}\n")

            readme.write("\n----Statistical summary tables----\n")
            for col, desc in stat_summary_metadata:
                readme.write(f"{col:16s}{desc}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-i", "--input", required=True, type=str, help="input file in .csv"
    )
    parser.add_argument("-k", "--key", type=str, required=True, help="plate setup file")
    parser.add_argument(
        "-o", "--output", required=True, type=str, help="output file name"
    )
    parser.add_argument(
        "--noblank",
        default=False,
        action="store_true",
        help="If you don't want to subtract the background. NOTE: must remove the background cells from your plate setup file that is passed to the -k flag. No blanking is actually the default currently, so this flag does nothing. If you want blanking, contact me.",
    )
    parser.add_argument(
        "--prism",
        default=False,
        action="store_true",
        help="use if you want processed data output compatible with GraphPad Prism",
    )
    parser.add_argument(
        "--readme_name",
        default="README.txt",
        help="name of readme file for metadata information, default=%(default)s",
    )
    parser.add_argument(
        "--ignore_treps",
        default=False,
        action="store_true",
        help="use if you want to ignore processing technical reps by themselves",
    )
    args = parser.parse_args()

    analyzer = GrowthCurveAnalyzer(
        args.input,
        args.key,
        args.output,
        args.noblank,
        args.prism,
    )
    analyzer.summary(data=analyzer.plate.data, save=True)

    if not args.ignore_treps:
        analyzer.summarize_technical_replicates()

    analyzer.write_readme(args.readme_name)
