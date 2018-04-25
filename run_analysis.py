#!/usr/env python

import sys
import os
import shutil
import subprocess
import time
import run # check that the path is correct

def run_analysis(jobs, cards_list, logs):
    if run.Run.is_scarf():
        n_procs = 40
    else:
        n_procs = 3
    config = {
        'jobs':jobs,
        'n_procs':n_procs,
        'extra_args':[],
        'delta_t':60,
        'max_t':60*60*24,
        'logs':logs,
    }
    for cards in cards_list:
        print "\nRunning", len(jobs), "with cards", cards_list
        analysis_pool = run.Run(config)
        analysis_pool.prefix = cards
        analysis_pool.run_many()

def main_mc_analysis():
    job_list = range(12)
    cards_list = ["scripts/config_mc.py",] #] # "scripts/config_mc.py", 
    logs = 'logs/mc-logs'
    run_analysis(job_list, cards_list, logs)

def main_reco_analysis():
    job_list = range(12)
    cards_list = ["scripts/config_reco.py",] #] # "scripts/config_mc.py", 
    logs = 'logs/reco-logs'
    run_analysis(job_list, cards_list, logs)
        
def main_systematics_analysis():
    job_list = range(1)
    cards_list = ["scripts/config_mc_systematics.py",] #
    logs = 'logs/systematics-logs'
    run_analysis(job_list, cards_list, logs)

if __name__ == "__main__":
    #main_mc_analysis()
    #main_systematics_analysis()
    main_reco_analysis()
