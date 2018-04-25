#!/usr/env python

import copy
import subprocess
import time
import os
import shutil

class Run(object):
    def __init__(self, config):
        self.prefix = 'simulation'
        self.config = config
        self.jobs = config['jobs']
        self.n_procs = config['n_procs']
        self.extra_args = config['extra_args']
        self.delta_t = config['delta_t']
        self.max_t = config['max_t']
        self.logs = config['logs']
        self.processes = []
        self.make_log_dir()

    def make_log_dir(self):
        try:
            shutil.rmtree(self.logs)
        except OSError:
            pass
        os.makedirs(self.logs+"/bsub")
        

    def run_many(self):
        print 'Running', len(self.jobs), 'jobs with cards', self.prefix,
        if self.is_scarf():
            print "on scarf"
        else:
            print "locally"
        print self.jobs
        timer = 0
        job_list = copy.deepcopy(self.jobs)

        self.processes = []
        while len(self.processes) > 0 or len(job_list) > 0:
            print round(timer/60., 1), '    ',
            self.poll()
            print len(self.processes)
            while len(self.processes) < self.n_procs and len(job_list) > 0:
                job = job_list.pop(0)
                self.processes.append(self.run_one(job))
            timer += self.delta_t
            if timer > self.max_t:
                print 'Out of time ... killing all jobs'
                self.kill_all()
                job_list = []
            time.sleep(self.delta_t)
        try:
            print subprocess.check_output(['bjobs'])
        except OSError:
            pass # not on scarf
        print "Done"

    def run_one(self, unique_id):
        unique_id = str(unique_id)
        log_name = self.logs+"/analysis_"+unique_id+".log"
        cards = self.prefix
        run = ['python', 'calculate_emittance.py', cards, unique_id]+self.extra_args
        bsub_name = ""
        if self.is_scarf():
            bsub = ['bsub',
                    '-n', '1',
                    '-W', '24:00',
                    '-q', 'scarf-ibis',
                    '-o', log_name,
                    '-e', log_name,
                ]
            run = bsub+run
            bsub_name = self.prefix.split(".")[0]
            bsub_name = bsub_name.split("/")[-1]
            bsub_name = self.logs+"/bsub/"+bsub_name+"_"+str(unique_id)+"_bsub.log"
            log_file = open(bsub_name, 'w')
        else:
            log_file = open(log_name, 'w')
        subproc = subprocess.Popen(run, stdout=log_file, stderr=subprocess.STDOUT)
        print "\nRunning", cards, unique_id, "in process id", subproc.pid
        return subproc, bsub_name

    @classmethod
    def is_scarf(cls):
        uname = subprocess.check_output(['uname', '-a'])
        return 'scarf.rl.ac.uk' in uname

    def get_bjob_number(self, bsub_name):
        line = open(bsub_name).readline()
        if line == "":
            return None
        bjob_number = line.split(' ')[1]
        bjob_number = bjob_number.rstrip('>')
        bjob_number = bjob_number.lstrip('<')
        return int(bjob_number)

    def poll(self, verbose = False):
        if verbose:
            print '\nPolling local'
        processes_update = self.poll_local(verbose)
        if self.is_scarf():
            if verbose:
                print 'Polling scarf'
            processes_update += self.poll_scarf(verbose)
        # remove duplicates
        processes_update = list(set(processes_update))
        self.processes = processes_update
            
    def poll_local(self, verbose):
        processes_update = []
        for proc, dir_name in self.processes:
            if proc.returncode == None:
                proc.poll()
                if proc.returncode == None:
                    processes_update.append((proc, dir_name))
            if verbose:
                print '   ', proc.pid, dir_name, proc.returncode
        return processes_update

    def poll_scarf(self, verbose):
        processes_update = []
        for proc, bsub_name in self.processes:
            bjob_number = self.get_bjob_number(bsub_name)
            output = subprocess.check_output(['bjobs', '-prw', str(bjob_number)])
            okay = False
            for line in output.split('\n'):
                is_alive = False
                for alive_key in ['PEND', 'RUN']:
                    is_alive = is_alive or alive_key in line
                if self.prefix in line and str(bjob_number) in line and is_alive:
                    processes_update.append((proc, bsub_name))
                    okay = True
                    break
            if not okay:
                print "Process with id", proc.pid, "and bsub log", bsub_name, "died"
            short_text = subprocess.check_output(['bjobs', str(bjob_number)])
            if verbose:
                print '   ', proc.pid, dir_name, bjob_number, short_text 

        return processes_update

    def kill_all(self):
        self.kill_all_local()
        if self.is_scarf():
            self.kill_all_scarf()
    
    def kill_all_local(self):
        for proc, dir_name in self.processes:
            if proc.returncode != None:
                continue # proc has returned
            pid_str = str(proc.pid)
            subprocess.check_output(['kill', '-9', pid_str])
            
    def kill_all_scarf(self):
        for proc, dir_name in self.processes:
            bjob_number = self.get_bjob_number(dir_name)
            try:
                output = subprocess.check_output(['bkill', str(bjob_number)])
            except Exception:
                pass

