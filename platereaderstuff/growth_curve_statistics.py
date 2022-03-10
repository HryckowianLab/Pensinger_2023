#!/usr/bin/env python3
"""
Author: Cody Martin
        University of Wisconsin-Madison
        Department of Bacteriology
        Anantharaman lab

Purpose: Generate summary statistics for growth curve data from 
         Hryckowian plate reader (Gen5?)

Usage: ./growth_curve_statistics.py -i <input.csv> -k <plate_setup.csv> -o <output_basename> [OPTIONS]
"""

import argparse
import pandas as pd
import numpy as np
from functools import reduce

# TODO: 2) output technical replicates for superplots


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
        return pd.read_csv(file)

    def read_setup(self, file: str) -> pd.DataFrame:
        return pd.read_csv(file, index_col="Row")

    def parse_key(self, file: str) -> pd.DataFrame:
        setup = self.read_setup(file)
        coordinates = list()
        for col, rowseries in setup.notna().iteritems():
            series = setup.loc[rowseries, col]
            if not series.empty:
                wells = list(map(lambda x: f"{x}{series.name}", series.index.to_list()))
                coordinates.append(list(zip(series.to_list(), wells)))

        sample_coordinates = reduce(lambda x, y: x + y, coordinates)

        return pd.DataFrame(sample_coordinates, columns=["Sample", "Well"])

    def update_time(self):
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
        return tmpdata.melt(id_vars="Time", var_name="Well", value_name="OD").merge(
            setup, on="Well", how="right"
        )

    # TODO: add optional blanking
    def process(self) -> pd.DataFrame:
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
        groups = self.setup.groupby("Sample")
        unique_rep_counts = groups.size().unique()
        if len(unique_rep_counts) == 1:
            n_reps = unique_rep_counts[0]
            technical_reps_setup = [groups.nth(i).reset_index() for i in range(n_reps)]
        else:
            print(
                "Not all of your samples have the same number of technical replicates. Contact Cody to handle this. Exiting."
            )
            exit()

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
        return reduce(lambda x, y: pd.concat([x, y]), self.tech_reps).reset_index(
            drop=True
        )

    def prism_technical_reps(self):
        combined = (
            self.combine_technical_reps()
            .assign(Sample_well=lambda x: x.Sample + "_" + x.Well)
            .drop(["Well", "Sample"], axis=1)
            .pivot(index="Time", columns="Sample_well")
        )

        columns = [name[1] for name in combined.columns.to_flat_index().to_list()]
        combined.columns = columns
        return combined[sorted(combined.columns)].reset_index()

    def prism_output(self) -> pd.DataFrame:
        data = self.data.drop(["Strain", "Media", "Concentration"], axis=1).pivot(
            index="Time", columns="Sample"
        )
        # print(data.columns)

        columns = [
            f"{sample}_{val}" for val, sample in data.columns.to_flat_index().to_list()
        ]
        # print(columns)
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

        if prism:
            self.plate.prism_output().to_csv(
                f"{self.baseoutput}_processed_data_PRISM.csv", index=False
            )
        else:
            self.plate.data.to_csv(f"{self.baseoutput}_processed_data.csv", index=False)

        self.plate.data["logOD"] = np.log(self.plate.data.OD)

    def smooth_data(self, data: pd.DataFrame, column: str, intervals=5) -> pd.DataFrame:
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
        smoothed = self.smooth_data(data, column, intervals)
        argmax = smoothed.groupby("Sample").idxmax()[column]
        return data.loc[argmax.reset_index(drop=True)].copy()

    def growth_rate(self, data: pd.DataFrame) -> pd.DataFrame:
        # the indices in the growth rates vectors are -1 compared to the actual sample/time indices
        # this means that the argmax index from the growth rate
        # corresponds to the time interval of
        # (argmax, argmax + 1) states
        g_rates = (
            data.groupby("Sample")
            .apply(lambda x: np.diff(x.logOD) / np.diff(x.Time))
            .reset_index()
            .rename({0: "growth_rate"}, axis=1)
        )

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
        growth_rates = (
            self.growth_rate(data)
            .groupby("Sample")
            .rolling(intervals, min_periods=1)
            .mean()
            .reset_index()
            .drop("level_1", axis=1)
        )
        # max_grates = growth_rates.groupby("Sample").growth_rate.max().reset_index()
        argmax = growth_rates.groupby("Sample").idxmax().growth_rate

        return growth_rates.loc[argmax.reset_index(drop=True)]

    def lag_time(self, data: pd.DataFrame, intervals=5) -> pd.DataFrame:
        # use midpoint of time interval to set a single time as the max
        max_growth_rates = self.max_growth_rate(data, intervals).assign(
            lag_time=lambda x: (x.t1 + (x.t2 - x.t1) / 2) / 2
        )
        return max_growth_rates

    def doubling_time(self, data: pd.DataFrame, intervals=5) -> pd.DataFrame:
        log2 = np.log(2)
        return self.max_growth_rate(data, intervals).assign(
            doubling_time=lambda x: log2 / x.growth_rate
        )

    def summary(
        self, data: pd.DataFrame, save: bool, intervals=5, technical_rep=False
    ) -> pd.DataFrame:
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

        self.plate.prism_technical_reps().to_csv(
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
            ("doubling_time", "culture doubling time in units of 1/hour"),
        ]

        with open(name, "w") as readme:
            readme.write("----File descriptions----\n")
            for file, desc in file_descriptions:
                readme.write(f"{file:46s}{desc}\n")

            readme.write("\n----Statistical summary tables----\n")
            for col, desc in stat_summary_metadata:
                readme.write(f"{col:16s}{desc}\n")

    def plot(self):
        # there's a bug in here...
        import matplotlib.pyplot as plt
        from matplotlib.ticker import ScalarFormatter
        import seaborn as sns

        fig, axes = plt.subplots(1, 2, figsize=(12, 6), sharey="row")
        for ax, (strain, data) in zip(axes, self.plate.data.groupby("Strain")):
            sns.lineplot(
                x="Time",
                y="OD",
                hue="Concentration",
                style="Media",
                linewidth=4,
                ci=None,
                data=data,
                ax=ax,
            )
            ax.tick_params(
                axis="both",
                which="both",
                direction="in",
                top=True,
                right=True,
            )
            ax.set(
                xlabel="Time post inoculation (h)",
                yscale="log",
                ylim=(0.05, 1),
                title=strain,
            )
            ax.xaxis.set_minor_locator(plt.MultipleLocator(1))

        axes[0].yaxis.set_major_formatter(ScalarFormatter())
        axes[0].set_ylabel("OD$_{600}$")
        axes[0].legend(frameon=False, title="Concentration (mM)")
        axes[1].get_legend().remove()

        fig.tight_layout()
        fig.savefig(f"{self.baseoutput}.pdf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process growth curve data and calculate summary statistics"
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
        "--plot",
        default=False,
        action="store_true",
        help="if you want to see a plot example",
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

    if args.plot:
        analyzer.plot()

    analyzer.write_readme(args.readme_name)
