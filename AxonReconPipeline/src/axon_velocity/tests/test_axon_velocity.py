import axon_velocity as av
from pathlib import Path
import numpy as np
import unittest

toy_data_folder = Path(__file__).parent / "data" / "toy"


class TestExtractors(unittest.TestCase):
    def setUp(self):
        # load data
        self.template_bifurcation = np.load(toy_data_folder / "bifurcation" / "template.npy")
        self.locations_bifurcation = np.load(toy_data_folder / "bifurcation" / "locations.npy")
        self.template_sinusoidal = np.load(toy_data_folder / "sinusoidal" / "template.npy")
        self.locations_sinusoidal = np.load(toy_data_folder / "sinusoidal" / "locations.npy")
        self.fs = 32000

        self.gtr_bif = av.compute_graph_propagation_velocity(self.template_bifurcation, self.locations_bifurcation,
                                                             self.fs, verbose=True)
        self.gtr_sin = av.compute_graph_propagation_velocity(self.template_sinusoidal, self.locations_sinusoidal,
                                                             self.fs, verbose=True)

    def tearDown(self):
        pass

    def test_axon_velocity(self):
        # output structure looks good
        assert len(self.gtr_bif.branches) > 0
        for branch in self.gtr_bif.branches:
            assert "channels" in branch
            assert "r2" in branch
            assert "velocity" in branch
            assert "peak_times" in branch
            assert "distances" in branch
            assert "offset" in branch

        # detection removed some channels
        assert len(self.locations_bifurcation) > len(self.gtr_bif.selected_channels)

        # output structure looks good
        assert len(self.gtr_sin.branches) > 0
        for branch in self.gtr_sin.branches:
            assert "channels" in branch
            assert "r2" in branch
            assert "velocity" in branch
            assert "peak_times" in branch
            assert "distances" in branch
            assert "offset" in branch

        # detection removed some channels
        assert len(self.locations_bifurcation) > len(self.gtr_sin.selected_channels)

    def test_plotting(self):
        av.plot_template(self.template_bifurcation, self.locations_bifurcation)
        av.plot_peak_latency_map(self.template_bifurcation, self.locations_bifurcation, self.fs)
        av.plot_peak_latency_map(self.template_bifurcation, self.locations_bifurcation, self.fs, plot_image=False)
        av.plot_amplitude_map(self.template_bifurcation, self.locations_bifurcation)
        av.plot_amplitude_map(self.template_bifurcation, self.locations_bifurcation, plot_image=False)
        av.plot_peak_std_map(self.template_bifurcation, self.locations_bifurcation, self.fs)
        av.plot_peak_std_map(self.template_bifurcation, self.locations_bifurcation, self.fs, plot_image=False)
        av.plot_branch_velocities(self.gtr_bif.branches)
        av.plot_axon_summary(self.gtr_bif)
        av.plot_branch_velocities(self.gtr_bif.branches)

        for br in self.gtr_bif.branches:
            av.plot_velocity(br["peak_times"], br["distances"], br["velocity"], br["offset"])
            av.plot_template_propagation(self.template_bifurcation, self.locations_bifurcation, br["channels"])
