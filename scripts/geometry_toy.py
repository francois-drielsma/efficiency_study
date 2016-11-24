import math


class GeometryToy(object):
    def __init__(self):
        self.tracker_aperture = 150.
        self.absorber_aperture = 160.
        self.tracker_distance = 1818.*2.
        self.tof_distance = 8224.8
        self.c_light = 299.792458
        self.mass = 105.6583715

    def max_path_length(self):
        theta = math.atan((self.tracker_aperture+self.absorber_aperture)/self.tracker_distance/2.)
        max_path_length = math.cos(theta)*self.tof_distance
        return max_path_length

    def p_tot(self, tof12, distance):
        beta = distance/tof12/self.c_light
        gamma = 1./(1-beta**2)**0.5
        p_tot = self.mass*beta*gamma
        return p_tot

    def print_summary(self, tof_low, tof_high, tof_step):
        tof12 = tof_low
        while tof12 < tof_high + tof_step/20.:
            d1 = self.tof_distance
            d2 = self.max_path_length()
            p1 = self.p_tot(tof12, d1)
            p2 = self.p_tot(tof12, d2)
            print str(round(tof12, 5)).ljust(10), \
                  str(round(d1, 3)).ljust(10), str(round(p1, 3)).ljust(10), \
                  str(round(d2, 3)).ljust(10), str(round(p2, 3)).ljust(10)
            tof12 += tof_step

if __name__ == "__main__":
    GeometryToy().print_summary(28.0, 32.0, 0.1)


