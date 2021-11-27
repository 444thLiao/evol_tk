"""
Merge multiple tabs after the genome assembly


checkM table
seqtk table for assembled
seqtk table for raw reads
"""
import os
from os.path import *

import pandas as pd

from dating_workflow.step_script import process_path


remaining_columns = {
    "checkm": [
        "marker lineage",
        "# markers",
        "# marker sets",
        "Completeness",
        "Contamination",
    ],
    "seqtk_assembly": [
        "metricsContigs_no",
        "metricsContigs_bp",
        "metricsContigs_ok",
        "metricsContigs_GC %",
        "metricsContigs_Max",
        "metricsContigs_avg",
        "metricsContigs_Min",
        "metricsContigs_N50",
        "metricsContigs_lt1K",
        "metricsContigs_lt5K",
    ],
    "seqtk_reads": [
        "metricsReads_avg_len",
        "metricsReads_AvgQual",
        "metricsReads_nReads(single file)",
    ],
}


def main(checkm_path, seqtk_assembly_path, seqtk_reads_path, ofile):
    # checkm_path = "/mnt/ivy/thliao/project/coral_ruegeria/data_processing/checkM_ruegeria_summary.tab"
    # seqtk_assembly_path = "/mnt/ivy/thliao/project/coral_ruegeria/data_processing/20211015/summary_output/seqtk_assembly_accessment.csv"
    # seqtk_reads_path = "/mnt/ivy/thliao/project/coral_ruegeria/data_processing/20211015/summary_output/seqtk_reads_accessment.csv"
    checkm_df = pd.read_csv(checkm_path, sep="\t", index_col=0)
    sa = pd.read_csv(seqtk_assembly_path, sep=",", index_col=0)
    sr = pd.read_csv(seqtk_reads_path, sep=",", index_col=0)

    sub_cm_df = checkm_df.loc[:, remaining_columns["checkm"]]
    sub_sa_df = sa.loc[:, remaining_columns["seqtk_assembly"]]
    sub_sr_df = sr.loc[:, remaining_columns["seqtk_reads"]]

    merged_df = pd.concat([sub_cm_df, sub_sa_df, sub_sr_df], axis=1)
    try:
        merged_df.to_excel(ofile, index=1)
    except:
        merged_df.to_excel(ofile, sep="\t", index=1)


if __name__ == "__main__":
    import sys

    cm, sa, sr, of = sys.argv[1:]
    main(cm, sa, sr, of)

    # python /home-user/thliao/script/evol_tk/bin/merge_tabs.py ./checkM_ruegeria_summary.tab 20211015/summary_output/seqtk_assembly_accessment.csv 20211015/summary_output/seqtk_reads_accessment.csv ./summary_info.xlsx