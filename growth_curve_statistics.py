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
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Args:
    input: Path
    key: Path
    output: Path
    prism: bool
    readme: Path
    average_replicates: bool
    window: int
    conditions: list[str]
    condition_delimiter: str


def parse_args() -> Args:
    parser = argparse.ArgumentParser(
        description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    io_group = parser.add_argument_group("I/O -- All required")
    io_group.add_argument(
        "-i",
        "--input",
        metavar="FILE",
        required=True,
        type=Path,
        help="input csv file for well-plate OD measurements",
    )
    io_group.add_argument(
        "-k",
        "--key",
        metavar="FILE",
        type=Path,
        required=True,
        help="plate setup csv file that labels each well with a sample identifier. Blank cells in this file are ignored. (see example)",
    )
    io_group.add_argument(
        "-o",
        "--output",
        metavar="FILE",
        required=True,
        type=Path,
        help="output file name",
    )

    condition_group = parser.add_argument_group(
        "Adjust experiment condition naming in the setup/key file"
    )
    condition_group.add_argument(
        "-c",
        "--conditions",
        nargs="+",
        default=["Strain", "Media", "Concentration"],
        help="space separated names of experimental conditions for each well. Must be equal to the number of divisions in the setup file as specified by the -d delimiter argument. (default: %(default)s)",
    )
    condition_group.add_argument(
        "-d",
        "--delimiter",
        metavar="STR",
        default="_",
        help="delimiter used in the setup file to separate experimental conditions. (default: '%(default)s')",
    )

    misc_args = parser.add_argument_group("OTHER")
    misc_args.add_argument(
        "--prism",
        default=False,
        action="store_true",
        help="use if you want processed data output compatible with GraphPad Prism",
    )
    misc_args.add_argument(
        "--readme",
        metavar="FILE",
        default="README.txt",
        type=Path,
        help="name of readme file for metadata information, default=%(default)s",
    )
    misc_args.add_argument(
        "--average-replicates",
        default=False,
        action="store_true",
        help="use if you want to average results over all technical replicates at the end",
    )
    misc_args.add_argument(
        "-w",
        "--window",
        default=5,
        metavar="INT",
        type=int,
        help="number of time intervals to smoothen data when computing all metrics. Lower this if log phase is not sampled frequently enough. A value of -w=1 corresponds to no smoothening. (default: %(default)s)",
    )
    _args = parser.parse_args()
    args = Args(
        input=_args.input,
        key=_args.key,
        output=_args.output,
        prism=_args.prism,
        readme=_args.readme,
        average_replicates=_args.average_replicates,
        window=_args.window,
        conditions=_args.conditions,
        condition_delimiter=_args.delimiter,
    )
    return args


def split_conditions(data: pd.DataFrame, conditions: list[str], delimiter: str):
    data[conditions] = data["Sample"].str.split(pat=delimiter, expand=True)


class WellPlate:
    def __init__(
        self,
        data_file: Path,
        setup_file: Path,
        conditions: list[str],
        delimiter: str,
        prism: bool,
    ) -> None:
        self._data = self.read_data(data_file)
        self._setup = self.read_setup(setup_file)
        self.data = self.combine()

        self.conditions = conditions
        self.delimiter = delimiter
        split_conditions(data=self.data, conditions=conditions, delimiter=delimiter)
        self.prism = prism

    def read_data(self, file: Path) -> pd.DataFrame:
        data = pd.read_csv(file)
        time = pd.to_timedelta(
            np.where(
                data.Time.str.count(":") == 1,
                data.Time + ":00",
                data.Time,
            )
        ) / pd.Timedelta(hours=1)
        data["Time"] = time
        data = data.melt(id_vars="Time", var_name="Well", value_name="OD").dropna()
        # data["logOD"] = np.log(data["OD"])

        return data

    def read_setup(self, file: Path) -> pd.DataFrame:
        setup = (
            pd.read_csv(file)
            .melt(id_vars="Row", var_name="Column", value_name="Sample")
            .assign(Well=lambda df: df["Row"] + df["Column"])
            .dropna()
        )
        # setup[conditions] = setup["Sample"].str.split(condition_delimiter, expand=True)
        return setup

    def combine(self) -> pd.DataFrame:
        data = self._data.merge(self._setup, how="inner", on="Well")
        return data

    def convert_to_prism_fmt(self, data: pd.DataFrame) -> pd.DataFrame:
        wide = data.pivot(index="Time", columns="Sample", values="OD")
        wide.columns.name = None
        return wide.reset_index()

    def process(self) -> pd.DataFrame:
        data = self.data.groupby(by=["Time", "Sample"])["OD"].mean().reset_index()
        split_conditions(data, self.conditions, self.delimiter)
        return data

    def write(self, output: Path, float_fmt: str):
        data = self.process()
        if self.prism:
            data = self.convert_to_prism_fmt(data)
            outdir = output.parent
            basename = output.stem
            output = outdir.joinpath(f"{basename}_PRISM.csv")

        data.to_csv(output, index=False, float_format=float_fmt)


def smooth_data(data: pd.DataFrame, window: int) -> pd.DataFrame:
    smoothed: pd.DataFrame = (
        data.groupby(["Sample", "Well"])
        .rolling(window, min_periods=1, on="Time")["OD"]
        .mean()
        .reset_index()
    )
    return smoothed


def calculate_growth_rate(data: pd.DataFrame) -> pd.DataFrame:
    data["logOD"] = np.log(data["OD"])
    growth_rates = (
        data.groupby(["Sample", "Well"])
        .apply(lambda x: np.diff(x["logOD"]) / np.diff(x["Time"]))  # type: ignore
        .reset_index()
        .rename({0: "growth_rate"}, axis=1)
        .explode("growth_rate")
        .reset_index(drop=True)
        .astype({"growth_rate": float})
    )

    times = data["Time"].unique()
    time_intervals = list(zip(*np.vstack([times[:-1], times[1:]])))
    time_intervals *= len(growth_rates) // len(time_intervals)

    return pd.concat(
        [pd.DataFrame(time_intervals, columns=["t1", "t2"]), growth_rates], axis=1
    )


def summarize_growth(
    data: pd.DataFrame, window: int = 2, average_replicates: bool = False
) -> pd.DataFrame:
    # smooth once
    data = smooth_data(data, window)
    # calculate growth rates using the differential of logOD
    growth_rates = calculate_growth_rate(data)

    # pos of max growth rate
    argmax = (
        growth_rates.groupby(["Sample", "Well"])["growth_rate"]
        .idxmax()
        .reset_index(drop=True)
    )

    summary: pd.DataFrame = (
        growth_rates.loc[argmax]
        .assign(
            # lag time is half time to max growth rate
            lag_time=lambda x: ((x.t2 + x.t1) / 2) / 2,
            # doubling time is log(2) / max_growth_rate
            doubling_time=lambda x: np.log(2) / x.growth_rate,
        )
        .reset_index(drop=True)
        .rename(columns={"growth_rate": "max_growth_rate"})
    )

    # reorder cols
    summary = summary[
        ["Sample", "Well", "t1", "t2", "max_growth_rate", "lag_time", "doubling_time"]
    ]

    if average_replicates:
        aggfn = ["mean", "std"]
        cols = ["t1", "t2", "max_growth_rate", "lag_time", "doubling_time"]
        aggcols = {
            f"{col}_{agg}": pd.NamedAgg(column=col, aggfunc=agg)
            for col in cols
            for agg in aggfn
        }
        summary = summary.groupby("Sample").agg(**aggcols).reset_index()
    return summary


def write_results(
    data: WellPlate,
    summary: pd.DataFrame,
    output: Path,
):
    outdir = output.parent
    basename = output.stem

    float_fmt = "%.4f"
    data.write(outdir.joinpath(f"{basename}_processed_data.csv"), float_fmt=float_fmt)

    summary.to_csv(
        outdir.joinpath(f"{basename}_statistical_summary.csv"),
        index=False,
        float_format=float_fmt,
    )


def write_readme(file: Path):
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
            "_statistical_summary.csv",
            "Growth curve statistical summary for a single biological replicate with values averaged across technical replicates",
        ),
    ]
    stat_summary_metadata = [
        ("Time", "time at max OD or max logOD"),
        ("Sample", "Full sample name from plate setup file"),
        # ("max_OD", "max OD value"),
        # ("Strain", "Strain, 1st part of Sample"),
        # ("Media", "Media, 2nd part of Sample"),
        # ("Concentration", "Concentration with units, 3rd part of Sample"),
        # ("max_logOD", "max log OD value ** Note: this is natural log"),
        # ("final_OD", "final OD value smoothed over the last 5 intervals"),
        # ("final_logOD", "final lnOD value smoothed over the last 5 intervals"),
        ("t1", "left side of time interval for max growth rate"),
        ("t2", "right side of time interval for max growth rate"),
        (
            "max_growth_rate",
            "maximum growth rate smoothed over n time intervals in units of lnOD/h",
        ),
        ("lag_time", "culture lag time in hours"),
        ("doubling_time", "culture doubling time in units of hours"),
    ]
    if not file.exists():
        with file.open("w") as fp:
            fp.write("----File descriptions----\n")
            for fileext, desc in file_descriptions:
                fp.write(f"{fileext:46s}{desc}\n")

            fp.write("\n----Statistical summary tables----\n")
            for col, desc in stat_summary_metadata:
                fp.write(f"{col:16s}{desc}\n")


def main():
    args = parse_args()
    print(f"1. Processing data at [{args.input}, {args.key}]")
    data = WellPlate(
        data_file=args.input,
        setup_file=args.key,
        conditions=args.conditions,
        delimiter=args.condition_delimiter,
        prism=args.prism,
    )

    print("2. Summarizing growth statistics")
    summary = summarize_growth(
        data=data.data,
        window=args.window,
        average_replicates=args.average_replicates,
    )

    split_conditions(
        data=summary,
        conditions=args.conditions,
        delimiter=args.condition_delimiter,
    )

    print(f"3. Writing results to {args.output.stem}*.csv")
    write_results(data=data, summary=summary, output=args.output)
    write_readme(args.readme)


if __name__ == "__main__":
    main()
