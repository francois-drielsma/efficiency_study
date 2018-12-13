#!/usr/env python

import sys
import os
import shutil
import subprocess
import time
import utilities.run # check that the path is correct

def run_analysis(jobs, cards_list, logs):
    if utilities.run.Run.is_scarf():
        n_procs = 150
    else:
        n_procs = 3
    config = {
        'jobs':jobs,
        'n_procs':n_procs,
        'extra_args':[],
        'delta_t':60,
        'max_t':60*60*96,
        'logs':logs,
    }
    print "\nStaging", len(jobs), "jobs with cards", cards_list
    analysis_pool = utilities.run.Run(config)
    analysis_pool.prefix = cards_list
    analysis_pool.run_many()

def main_mc_analysis():
    job_list = range(16)
    cards_list = ["scripts/config/config_mc.py",]
    logs = 'logs/mc-logs'
    run_analysis(job_list, cards_list, logs)

def main_reco_analysis():
    job_list = range(16)
    cards_list = ["scripts/config/config_reco.py",]
    logs = 'logs/reco-logs'
    run_analysis(job_list, cards_list, logs)
        
def main_systematics_analysis():
    job_list = range(140)
    cards_list = ["scripts/config/config_mc_systematics.py",] #
    logs = 'logs/systematics-logs'
    run_analysis(job_list, cards_list, logs)

def main_both_analysis():
    job_list = range(16)
    cards_list = ["scripts/config/config_mc.py", "scripts/config/config_reco.py"]
    logs = 'logs/both-logs'
    run_analysis(job_list, cards_list, logs)

def main():
    print "Running\n    ",
    for item in sys.argv:
        print item,
    print
    if len(sys.argv) < 2 or sys.argv[1] not in ["reco", "sys", "mc"]:
        print """Usage:
    python run_many_analyses.py <option>
where option is one of 'reco' 'sys' 'mc'
        """
        return
    if sys.argv[1] == "reco":
        main_reco_analysis()
    elif sys.argv[1] == "sys":
        main_systematics_analysis()
    elif sys.argv[1] == "mc":
        main_mc_analysis()
    return

if __name__ == "__main__":
    main()
