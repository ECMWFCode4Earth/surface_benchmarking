import os
from os import makedirs as mkdir
import sys
from warnings import filterwarnings
import pprint
from xml.sax import make_parser

filterwarnings("ignore")
import argparse
import pandas as pd
import time


def showProgress(iter_ix, num_iters, bar_len=49):
    """
    Show progress bar with percentage.
    Arguments:
        iter_ix: iteration index (zero-based).
        num_iters: total number of iterations.
        bar_len (int): bar length.
    """
    progress = float(iter_ix + 1) / num_iters
    arrow_len = int(round((bar_len * progress)))
    sys.stdout.write(
        "\r[{0}>{1}] {2}%".format(
            "-" * arrow_len,
            " " * (bar_len - arrow_len),
            str(int(round(100 * progress))).zfill(2),
        )
    )
    sys.stdout.flush()
    if iter_ix + 1 == num_iters:
        print("\n")


def ensure_dir(file_path):
    if not os.path.exists(file_path):
        os.makedirs(file_path)

def clean_dir(file_path):
    for f in os.listdir(file_path):
            os.remove(os.path.join(file_path, f))


def get_parser(doc):
    """Get a command-line argument parser for the program."""

    parser = argparse.ArgumentParser(
        description=doc,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("options_file", type=str, help=("Options file"))

    return parser


class LDAS_config(object):
    """
    Class defining default variables for LANDVER
    """

    def __init__(self, dictionary):
        """Constructor"""
        for key in dictionary:
            setattr(self, key, dictionary[key])

        self.check_contents()

    def check_contents(self):
        required_atts = {
            "analysis_period": list,
            "output_dir": str,
            "EXPVER": list,
            "CLASS": list,
            "TYPE": list,
            "ANOFFSET": list,
            "STREAM": list,
            "REPRES": list,
            "LEVTYPE": list,
            "TIME": list,
            #"TIME_RES": list,
            "STEP": list,
            "DOMAIN": list,
            "GRID": list,
            "EXPNAME": list,
        }

        self.config_file = None
        #self.extract_SH = False  # Extract grib files from MARS (required for preprocessing) #not yet implemented for fluxes
        #self.extract_LH = False  # Extract grib files from MARS (required for preprocessing)
        self.pre_process_SH = (
            False  # Preprocess sensible heat flux (required for validation)
        )
        self.pre_process_LH = ( 
            False # Preprocess latent heat flux (required for validation)
        )
        self.lt2utc=False #convert lt to utc: False = no conversion, True = convert insitu data from lt to utc

        self.validate_SH = True  # Validate sensible heat flux
        self.validate_LH = True  # Validate latent heat flux

        self.validation_times=["all_times"] # 00_UTC and 12_UTC are optional
        self.climate_classes=["all_climate"]
        self.land_classes=["all_land"] 
        self.Network = ['ICOS_FLUXNET'] #currently just "ICOS_FLUXNET" for fluxes
        self.in_situ_dir = "/home/lauracma/Documents/ecmwf_proj/data/ISMN_data/in_situ_data/pre_processed_data/"  # adapted
        self.plot_time_series = True  # Plot time series of soil moisture and temperature for analysis and in situ data
        self.Rescale_data = (
            False  # Rescale in situ and analysis SM using max/min values = 1.0/0.0
        )
        self.SH_units = "J"  # units can be "J" for J/m^2 or "W" for W/m^2
        self.LH_units = "J"
        self.fluxes_units = "J"
        self.validation_layer = [
            "Surface"
        ]  # just surface
        self.ST_quality_control = False  # Screen SM where soil temperature less than specified value ST_QC_threshold #not yet adapted
        self.ST_QC_threshold = (
            277.15  # Minimum soil temperature threshold (degrees Kelvin)
        )
        # Screen stations where orography difference between in situ station and analysis grid exceeds specified threshold:
        self.filterOro = False
        self.clim_file = "orog_639l_2"  # orog clim file of the experiments
        self.Oro_diff_threshold = 300.0  # Threshold for difference in height between analysis gridpoint and in situ station
        self.clim_dir = "/perm/rd/dadf/static/"  # Location of orography files
        self.ST_ML10 = [False, False, False, False]
        # Use observations and analysis data on same local time (according to longitude) to compare diurnal cycle.
        self.ST_compare_lt = False
        self.table_format = (
            "ASCII"  # Format for output table scores ('ASCII' or 'PDF').
        )
        self.daily_obs_time_average = False  # If false, the in situ observation is used at 00UTC (+-30 mins). If true, the daily time
        # average is taken (averaged over 24 hours centred on 00 UTC).
        self.stat_sig = True
        self.min_obs = 5
        self.variable = "land_validation"
        self.executable_path = os.path.realpath(sys.argv[0])

        for key in required_atts:
            assert hasattr(
                self, key
            ), f"Required input key {key} missing from input file"
            assert isinstance(
                getattr(self, key), required_atts[key]
            ), f"Required input key {key} must be a {required_atts[key]} but is instead a {type(getattr(self,key))}"

        Networks = ['ICOS_FLUXNET']
        for f in self.Network:
            assert (
                f in Networks
            ), f"Unknown network '{f}'. Must be one of {sorted(Networks)}"

        retrieve_len = len(self.EXPVER)
        retrieve_atts = [
            "EXPVER",
            "CLASS",
            "TYPE",
            "ANOFFSET",
            "STREAM",
            "REPRES",
            "LEVTYPE",
            "TIME",
            "STEP",
            "DOMAIN",
            "GRID",
            "EXPNAME",
        ]
        for key in retrieve_atts:
            assert (
                len(getattr(self, key)) == retrieve_len
            ), f"Retrieval key '{key}' does not have the same length as EXPVER. Check your settings"

    def __repr__(self):
        return str(type(self)) + "\n" + pprint.pformat(self.__dict__)


def get_configuration(options_file):
    config_dict = {}
    try:
        # Execute the specified Python file, storing the result in the
        # dictionary config_dict:
        # NOTE: We could use `execfile()` on Python 2, but opening the
        #       file and calling `exec()` works on both Python 2 and 3.
        with open(options_file, "r") as config_module:
            exec(config_module.read(), globals(), config_dict)
        cfg = LDAS_config(config_dict)
        cfg.config_file = os.path.realpath(options_file)
    except IOError as e:
        raise IOError(f"failed to load configuration from file: {str(e)}")
    return cfg


def web_summary(cfg):
    print(
        """<h4>Summary of experiments</h4>
<table border=1>
<tr><th colspan=3>Being inspected</th></tr>
<tr><th>Exp. ID</th><th>Class</th><th>Description</th></tr>"""
    )
    for expver, class_type, name in zip(cfg.EXPVER, cfg.CLASS, cfg.EXPNAME):
        line = f"<tr><td>{expver}</td><td>{class_type}</td><td>{name}</td></tr>"
        print(line)
    print("</table>")
    return None


def make_clickable(url, name):
    return '<a href="{}" rel="noopener noreferrer" target="_blank">{}</a>'.format(
        url, name
    )


def web_obs(cfg, df):

    indices = [
        "metric",
        "expver",
        "Surface SH",
        "Surface LH",
    ]
    df = (
        df.set_index(["expver", "date", "type", "network", "metric", "variable"])[
            "link"
        ]
        .unstack()
        .reset_index()
    )

    if "Surface SH" in df.columns:
        df["Surface SH"] = df.apply(
            lambda x: make_clickable(x["Surface SH"], "sshf"), axis=1
        )
    #else:
    #    indices.remove("Surface SH")
    if "Surface LH" in df.columns:
        df["Surface LH"] = df.apply(
            lambda x: make_clickable(x["Surface LH"], "slhf"), axis=1
        )
    #else:
    #    indices.remove("Surface LH")

    map_plot = df.loc[df["type"] == "Map"]

    map_plot = map_plot.rename(columns={"metric": "country"})
    indices[0] = "country"

    map_plot = (
        map_plot.drop(["date", "type", "network"], axis=1)
        .set_index(indices)
        .sort_index(axis=0)
    )

    map_plot.sort_index(axis=1)

    #    print(map_plot['link'])
    #    map_plot["link"]=map_plot["link"].values[::-1]

    print("""<h4>Map plots: Pearson R and station locations</h4>""")

    print(
        """<input 
  type="text" 
  id="myInput" 
  onkeyup="myFunction()" 
  placeholder="Filter table..." 
  title="title text">"""
    )

    html1 = map_plot.to_html(
        sparsify=False,
        escape=False,
        header=False,
        table_id="myTable",
        render_links=True,
    )

    print(html1)

    print(
        """<h4>Confidence interval (95%) plots of Pearson R anomalies calculated using a Fisher Z transform with lag-1 autocorrelation. </h4>"""
    )

    ci_plot = df.loc[df["type"] == "CI"]
    indices[0] = "metric"

    ci_plot = (
        ci_plot.drop(["date", "type", "network"], axis=1)  # [df["network"] == "odb"]
        .set_index(indices)
        .sort_index()
    )

    print(
        """<input 
  type="text" 
  id="myInput2" 
  onkeyup="myFunction2()" 
  placeholder="Filter table..." 
  title="title text">"""
    )

    #    df.style.format({'link': make_clickable})

    html2 = ci_plot.to_html(
        sparsify=False,
        escape=False,
        header=False,
        table_id="myTable2",
        render_links=True,
    )

    # # add sorting to columns
    # for i, string in enumerate(indices):
    #     html=html.replace(f'<th>{string}',f'<th onclick="sortTable({i})">{string}')

    print(html2)
    #    print(df.columns())

    table_plot = df.loc[df["type"] == "Table"]

    indices_table = indices
    indices_table.insert(1, "network")

    table_plot = (
        table_plot.drop(["date", "type"], axis=1)  # [df["network"] == "odb"]
        .set_index(indices_table)
        .sort_index()
    )

    #    df.style.format({'link': make_clickable})
    print("""<h4>Tables of scores (R, R_anom, RMSD, Unbiased RMSD)</h4>""")

    print(
        """<input 
    type="text" 
    id="myInput3" 
    onkeyup="myFunction3()" 
    placeholder="Filter table..." 
    title="title text">"""
    )

    html3 = table_plot.to_html(
        sparsify=False,
        escape=False,
        header=False,
        table_id="myTable3",
        render_links=True,
    )

    # # add sorting to columns
    # for i, string in enumerate(indices):
    #     html=html.replace(f'<th>{string}',f'<th onclick="sortTable({i})">{string}')

    print(html3)

    box_plot = df.loc[df["type"] == "Box plot"]
    indices[0] = "metric"

    indices_box_plot = indices

    box_plot = (
        box_plot.drop(["date", "type"], axis=1)  # [df["network"] == "odb"]
        .set_index(indices_box_plot)
        .sort_index()
    )

    #    df.style.format({'link': make_clickable})

    print("""<h4>Box plots of scores (R, R_anom, RMSD, Unbiased RMSD) </h4>""")

    print(
        """<input 
    type="text" 
    id="myInput4" 
    onkeyup="myFunction4()" 
    placeholder="Filter table..." 
    title="title text">"""
    )

    html4 = box_plot.to_html(
        sparsify=False,
        escape=False,
        header=False,
        table_id="myTable4",
        render_links=True,
    )

    # # add sorting to columns
    # for i, string in enumerate(indices):
    #     html=html.replace(f'<th>{string}',f'<th onclick="sortTable({i})">{string}')

    print(html4)

    print("""<h4>Time series plots </h4>""")

    time_series_plot = df.loc[df["type"] == "Time series"]
    indices[0] = "metric"

    time_series_plot = (
        time_series_plot.drop(["date", "type"], axis=1)  # [df["network"] == "odb"]
        .set_index(indices)
        .sort_index()
    )

    print(
        """<input 
  type="text" 
  id="myInput5" 
  onkeyup="myFunction5()" 
  placeholder="Filter table..." 
  title="title text">"""
    )

    #    df.style.format({'link': make_clickable})

    html5 = time_series_plot.to_html(
        sparsify=False,
        escape=False,
        header=False,
        table_id="myTable5",
        render_links=True,
    )

    # # add sorting to columns
    # for i, string in enumerate(indices):
    #     html=html.replace(f'<th>{string}',f'<th onclick="sortTable({i})">{string}')

    print(html5)

    return None


def web_create(cfg, times, land_type, *argv):
    if all(v is None for v in argv):
        return None
    df = pd.concat(argv, ignore_index=True)

    print("Creating webpage...", end="\r")

    start = time.time()

    folder = f"{cfg.output_dir}/{'_'.join(cfg.EXPVER)}"
    mkdir(folder, exist_ok=True)
    index = f"{folder}/{cfg.variable}_{cfg.analysis_period[0]}_{cfg.analysis_period[1]}_{times}_{land_type}_index.html"
    indexhead = f"""<html>
<head>
<style>
h1, h2, h3, h4, h5, p, table, ul {{font-family:sans-serif;}}
table {{border-collapse:collapse;}}
td, th {{padding:3px 7px 2px 7px; }}
th {{background-color: #99ccff;}}
body {{background-color:white;}}
</style>
<title>Land validation - Preview of plots in {cfg.output_dir} </title>
</head>
<body>
<h2><b>LDAS validation for surface heat fluxes</b></h2>
<h3>Generated by {cfg.executable_path} {cfg.config_file}</h3>
<ul>
</ul>
"""

    indexfoot = """
<script>
const myFunction = () => {
  const trs = document.querySelectorAll('#myTable tr:not(.header)')
  const filter = document.querySelector('#myInput').value
  const regex = new RegExp(filter, 'i')
  const isFoundInTds = td => regex.test(td.innerHTML)
  const isFound = childrenArr => childrenArr.some(isFoundInTds)
  const setTrStyleDisplay = ({ style, children }) => {
    style.display = isFound([
      ...children // <-- All columns
    ]) ? '' : 'none' 
  }
  
  trs.forEach(setTrStyleDisplay)
}
const myFunction2 = () => {
  const trs = document.querySelectorAll('#myTable2 tr:not(.header)')
  const filter = document.querySelector('#myInput2').value
  const regex = new RegExp(filter, 'i')
  const isFoundInTds = td => regex.test(td.innerHTML)
  const isFound = childrenArr => childrenArr.some(isFoundInTds)
  const setTrStyleDisplay = ({ style, children }) => {
    style.display = isFound([
      ...children // <-- All columns
    ]) ? '' : 'none' 
  }
  
  trs.forEach(setTrStyleDisplay)
}
const myFunction3 = () => {
  const trs = document.querySelectorAll('#myTable3 tr:not(.header)')
  const filter = document.querySelector('#myInput3').value
  const regex = new RegExp(filter, 'i')
  const isFoundInTds = td => regex.test(td.innerHTML)
  const isFound = childrenArr => childrenArr.some(isFoundInTds)
  const setTrStyleDisplay = ({ style, children }) => {
    style.display = isFound([
      ...children // <-- All columns
    ]) ? '' : 'none' 
  }

  
  trs.forEach(setTrStyleDisplay)
}
const myFunction4 = () => {
  const trs = document.querySelectorAll('#myTable4 tr:not(.header)')
  const filter = document.querySelector('#myInput4').value
  const regex = new RegExp(filter, 'i')
  const isFoundInTds = td => regex.test(td.innerHTML)
  const isFound = childrenArr => childrenArr.some(isFoundInTds)
  const setTrStyleDisplay = ({ style, children }) => {
    style.display = isFound([
      ...children // <-- All columns
    ]) ? '' : 'none' 
  }
  
  trs.forEach(setTrStyleDisplay)
}
const myFunction5 = () => {
  const trs = document.querySelectorAll('#myTable5 tr:not(.header)')
  const filter = document.querySelector('#myInput5').value
  const regex = new RegExp(filter, 'i')
  const isFoundInTds = td => regex.test(td.innerHTML)
  const isFound = childrenArr => childrenArr.some(isFoundInTds)
  const setTrStyleDisplay = ({ style, children }) => {
    style.display = isFound([
      ...children // <-- All columns
    ]) ? '' : 'none' 
  }
  
  trs.forEach(setTrStyleDisplay)
}

</script>
"""

    with open(index, "w") as f:
        orig_stdout = sys.stdout
        sys.stdout = f

        print(indexhead)

        web_summary(cfg)

        web_obs(cfg, df)

        print(indexfoot)

        print("</body></html>")
        # finished writing to file, now revert STDOUT
        sys.stdout = orig_stdout

    print(f"Webpage created at {index} ({(time.time() - start):.2f}s)")
    print(f"Link to webpage: file://{index}")
    return None
