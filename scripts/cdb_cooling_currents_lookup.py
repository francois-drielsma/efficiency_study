import sys
import datetime
import time
import urllib2
import cdb
from cdb import Beamline
from cdb import CoolingChannel


MODULE_LIST = ['SSU', 'FCU', 'SSD']
COIL_LIST = ['SSU-E2', 'SSU-C', 'SSU-E1', 'SSU-M2', 'SSU-M1', 'FCU-C', 'SSD-M1', 'SSD-M2', 'SSD-E1', 'SSD-C', 'SSD-E2']
OPTICS_RUN_DICT = {}

def sort_predicate(module_or_coil):
    module_list = MODULE_LIST+COIL_LIST
    return module_list.index(module_or_coil['name'])

def print_heading(output = sys.stdout):
    print >> output, "run".ljust(8), "date".ljust(12), "beamline".ljust(20),
    for coil in COIL_LIST:
        print >> output, coil.rjust(8),
    print >> output

def print_run(bl_data, channel_data, output = sys.stdout):
    global OPTICS_RUN_DICT
    run = bl_data['run_number']
    optics = bl_data["optics"]
    if optics not in  OPTICS_RUN_DICT:
         OPTICS_RUN_DICT[optics] = []
    OPTICS_RUN_DICT[optics].append(run)
    start_time = bl_data["start_time"]
    print >> output, str(run).ljust(8), str(start_time.date()).ljust(12), ('"'+optics+'"').ljust(20),
    print >> output, channel_data['tag']
    #for module in sorted(channel_data, key=sort_predicate):
    #    for coil in sorted(module['coils'], key=sort_predicate):
    #        print >> output, str(round(coil['iset'], 2)).rjust(8),
    print >> output
    output.flush()

def is_channel_okay(bl_data, channel_data):
    if channel_data == None or bl_data == None:
        print "Channel data:", channel_data == None, "BL Data:", bl_data == None
        return False
    if len(channel_data) < 3:
        print "Len Channel data:", len(channel_data)
        return False
    return True
    channel_data = sorted(channel_data)
    if abs(channel_data[1]['coils'][0]['iset']) < 5.: # FC-C
        print "FC current", channel_data[1]['coils'][0]['iset']
        return False
    return True

def get_settings_by_run(run_start, run_end, run_step, output, do_header = True):
    beamline = Beamline(url='http://cdb.mice.rl.ac.uk')
    channel = CoolingChannel(url='http://cdb.mice.rl.ac.uk')

    if do_header:
        print_heading(output)
    for run in range(run_start, run_end, run_step): #:
        channel_data, bl_data = None, None
        n_tries = 0
        while channel_data == None and n_tries < 5:
            try:
                channel_data = channel.get_coolingchannel_for_run(run)
                print ".",
                sys.stdout.flush()
            except (cdb._exceptions.CdbPermanentError, KeyError, urllib2.URLError):
                print "Failed to get channel data for run", run, "... trying again", n_tries
                n_tries += 1
                time.sleep(0.1)
                continue
        while bl_data == None and n_tries < 5:
            try:
                bl_data = beamline.get_beamline_for_run(run)[run]
            except (cdb._exceptions.CdbPermanentError, KeyError, urllib2.URLError):
                print "Failed to get bl data for run", run, "... trying again", n_tries
                n_tries += 1
                time.sleep(0.1)
                continue
        if not is_channel_okay(bl_data, channel_data):
            print "Channel not okay for run", run
            continue
        print_run(bl_data, channel_data, output)
        sys.stdout.flush()
        time.sleep(1)

    #print >> output, "Finished; note runs with FC < 5 A or incomplete channel magnet data are excluded"

def get_settings_by_date(date_start, date_end):
    beamline = Beamline(url='http://cdb.mice.rl.ac.uk')
    channel = CoolingChannel(url='http://cdb.mice.rl.ac.uk')

    channel_data_dates = channel.get_all_coolingchannels()
    print "Getting cooling channel keys"
    for key in channel_data_dates.keys():
        print type(key), key
    return
    print_heading()
    bl_data_dates = beamline.get_beamlines_for_dates(start_date, end_date)
     
    for run in range(run_start, run_end):
        channel_data, bl_data = None, None
        n_tries = 0
        while channel_data == None and n_tries < 5:
            try:
                channel_data = channel.get_coolingchannel_for_run(run)
            except (cdb._exceptions.CdbPermanentError, KeyError, urllib2.URLError):
                n_tries += 1
                time.sleep(5)
                continue
        while bl_data == None:
            try:
                bl_data = beamline.get_beamline_for_run(run)[run]
            except (cdb._exceptions.CdbPermanentError, KeyError, urllib2.URLError):
                n_tries += 1
                time.sleep(5)
                continue
        if not is_channel_okay(bl_data, channel_data):
            continue
        print_run(bl_data, channel_data)

    #    print "Finished; note runs with FC < 5 A or incomplete channel magnet data are excluded"



#2016-04 8512 to 8948 
#2016-05 8995 to 9212 
def main():
    if len(sys.argv) > 2:
        start = int(sys.argv[1])
        try:
            end = int(sys.argv[2])
        except IndexError:
            end = start+1
        try:
            step = int(sys.argv[3])
        except IndexError:
            step = 1
        try:
            file_name = sys.argv[4]
        except IndexError:
            file_name = "settings.csv"
        output = open(file_name, "w")

        get_settings_by_run(start, end, step, output, True)
    if len(sys.argv) == 2:
        fin = open(sys.argv[1])
        output = open("settings.csv", "w")
        do_header = True
        for line in fin.readlines():
            run = int(line)
            get_settings_by_run(run, run+1, 1, output, do_header)
            do_header = False

if __name__ == "__main__":
    main()

