#!/usr/env python

import sys
import os
import shutil
import subprocess
import time
import xboa.common

def is_scarf():
    uname = subprocess.check_output(['uname', '-a'])
    return 'scarf.rl.ac.uk' in uname

def generate_analysis(unique_id, cards_list):
    unique_id = str(unique_id)
    log_name = "logs/analysis_"+unique_id+".log"
    log_file = open(log_name, 'w')
    for cards in cards_list:
        print "\nRunning", cards, unique_id
        run = ['python', 'calculate_emittance.py', cards, unique_id]
        subproc = subprocess.Popen(run, stdout=log_file, stderr=subprocess.STDOUT)
        yield subproc
        print "\nFinished", cards, unique_id

def poll_local(processes):
    _processes_update = []
    for generator, proc in processes:
        if proc == None or proc.returncode != None:
            #process needs to be started
            try:
                proc = generator.next()
                # running the next set of cards; keep it in the job list
                _processes_update.append((generator, proc))
                sys.stdout.write('x')
            except StopIteration:
                # no more analysis to run now; give up and make way for another job
                sys.stdout.write('x')
                pass
        else:
            proc.poll()
            _processes_update.append((generator, proc))
            sys.stdout.write('.')
    return _processes_update

def poll_bjobs():
    global SCRIPT_NAME
    output = subprocess.check_output(['bjobs', '-prw'])
    count = 0
    for line in output.split('\n'):
        if SCRIPT_NAME in line or 'analysis.py' in line:
            count += 1
    return count

def main(jobs, number_of_concurrent_processes, cards_list):
    print 'Running', len(jobs), 'jobs on scarf?', is_scarf()
    timer = 0
    time_step = 10
    processes = []
    bjobs = 0
    while len(processes) > 0 or len(jobs) > 0:
        # refill the process list with new jobs
        while len(processes) < number_of_concurrent_processes and len(jobs) > 0:
            job = jobs.pop(0)
            processes.append((generate_analysis(job, cards_list), None))
        # remove any processes that have finished
        print timer,
        processes = poll_local(processes)
        print len(jobs)+len(processes)
        timer += time_step
        time.sleep(time_step)
    unique_id = 0
    print "\nFinished ... press <cr> to end"

if __name__ == "__main__":
    jobs = range(12)
    if is_scarf():
        n_procs = min(len(jobs), 100)
    else:
        n_procs = 3
    cards_list = ["scripts/config_mc.py", "scripts/config_reco.py"]
    main(jobs, n_procs, cards_list)

