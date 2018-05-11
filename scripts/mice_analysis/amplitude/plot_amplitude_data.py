import numpy

import xboa.common

class PlotAmplitudeData(object):
      def __init__(self, amplitude_data, plot_dir, key):
          self.data = amplitude_data
          self.plot_dir = plot_dir
          self.key = key

      def plot(self):
          self.plot_data_1d("emittance_vs_n_events_"+self.key, self.emittance_4d_lambda, "#varepsilon_{4D} [mm]", self.n_events_lambda, "Number of Events")
          self.plot_data_1d("emittance_vs_beta_x_"+self.key, self.emittance_4d_lambda, "#varepsilon_{4D} [mm]", self.beta_x_lambda, "#beta_{x} [mm]")
          self.plot_data_1d("emittance_vs_beta_y_"+self.key, self.emittance_4d_lambda, "#varepsilon_{4D} [mm]", self.beta_y_lambda, "#beta_{y} [mm]")

      def plot_data_1d(self, plot_name, plot_lambda_x, x_label, plot_lambda_y, y_label):
          x_axis = []
          y_axis = []
          if len(self.data.state_list) == 0:
              print "Warning - no data for emittance vs beta plots/etc"
              return
          for state in self.data.state_list:
              x_axis.append(plot_lambda_x(state))
              y_axis.append(plot_lambda_y(state))
          canvas = xboa.common.make_root_canvas("plot")
          hist, graph = xboa.common.make_root_graph(plot_name, x_axis, x_label, y_axis, y_label)
          canvas.Draw()
          hist.Draw()
          graph.SetMarkerStyle(24)
          graph.Draw("SAME P")
          for fmt in ["png", "pdf", "root"]:
              canvas.Print(self.plot_dir+"/"+plot_name+"."+fmt)

      @classmethod
      def emittance_4d_lambda(cls, state):
          return state["emittance"]

      @classmethod
      def n_events_lambda(cls, state):
          return state["n_events"]

      @classmethod
      def beta_x_lambda(cls, state):
          return cls.beta_2d(state, 0)

      @classmethod
      def beta_y_lambda(cls, state):
          return cls.beta_2d(state, 2)

      @classmethod
      def beta_2d(cls, state, axis):
          twod_matrix = [item[axis:axis+2] for item in state["cov"][axis:axis+2]]
          emit = numpy.linalg.det(twod_matrix)**0.5/cls.mu_mass
          beta = twod_matrix[0][0]/emit
          return beta

      mu_mass = xboa.common.pdg_pid_to_mass[13]