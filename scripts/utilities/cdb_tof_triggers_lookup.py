import cdb

def get_material(run_number):
    # From 15/12/2017 change to "wedge" 10515 ... 10698
    # From 11/12/2017 change to LiH  10474 ... 10514
    # From 29/11/2017 change to empty 10301 ... 10473
    # From 27/10/2017 Empty lH2   10098 ... 10300
    # From 25/09/2017 Full lH2 9697? .. 10095 
    run_list = [
       ("Wedge",       10515, 10698),
       ("LiH",         10474, 10514),
       ("No absorber", 10301, 10473),
       ("Empty lH2",   10098, 10300),
       ("Full lH2",     9697, 10095),
    ]
    for absorber, run_start, run_end in run_list:
        if run >= run_start and run <= run_end:
            return absorber
    return "Not sure"

def get_beamline(run_number):
    bl = cdb.Beamline()
    beamline_all = bl.get_beamline_for_run(run_number)[run_number]
    return beamline_all

def get_channel(run_number):
    channel = cdb.CoolingChannel()
    channel_all = channel.get_coolingchannel_for_run(run_number)
    return channel_all

def get_tof_triggers_sum(tof_key, run_list):
    trig_list = [get_beamline(run_number)["scalars"][tof_key] for run_number in run_list]
    return sum(trig_list)

def get_beamline_tags(run_list):
    beamline_tags = [get_beamline(run_number)["optics"] for run_number in run_list]
    return beamline_tags

def get_channel_tags(run_list):
    channel_tags = [get_channel(run_number)["tag"] for run_number in run_list]
    return channel_tags

def get_absorber(run_number):
    channel = cdb.CoolingChannel()
    absorber = channel.get_absorber_for_run(run_number)
    return absorber.keys()[0][1]

def get_absorber_list(run_list):
    absorbers = [get_absorber(run_number) for run_number in run_list]
    return absorbers

def get_time_sum(run_list):
    start_time = [get_beamline(run_number)["start_time"] for run_number in run_list]
    end_time = [get_beamline(run_number)["end_time"] for run_number in run_list]
    delta_time = [end_time[i] - start for i, start in enumerate(start_time)]
    total_seconds = sum([delta.total_seconds() for delta in delta_time])
    minutes = int(total_seconds)/60
    hours = minutes/60
    minutes = minutes % 60
    return hours, minutes

def parse_one_setting(run_list):
    run_list = [run for run in run_list if run > 2000]
    tof1_triggers = get_tof_triggers_sum("ToF1 Triggers", run_list)
    tof2_triggers = get_tof_triggers_sum("ToF2 Triggers", run_list)
    lmc1234 = get_tof_triggers_sum("LMC-1234 Count", run_list)
    hours, minutes = get_time_sum(run_list)
    optics_list = list(set(get_beamline_tags(run_list)))
    channel_list = list(set(get_channel_tags(run_list)))
    absorber_list = list(set(get_absorber_list(run_list)))
    return {
        "runs":run_list,
        "tof1_triggers":tof1_triggers,
        "tof2_triggers":tof2_triggers,
        "lmc1234":lmc1234,
        "time":[hours, minutes],
        "bl":optics_list,
        "channel":channel_list,
        "absorber":absorber_list,
    }

if __name__ == "__main__":
    print parse_one_setting([9620] )
