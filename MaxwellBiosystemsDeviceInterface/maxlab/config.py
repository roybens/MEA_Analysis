class Config:

    def __init__(self, str):
        config = [m.replace('(', ' ').replace(')', ' ').replace('/', ' ').split() for m in str.split()[0].split(';')[:-1]]
        self.mappings = [self.Mapping(*m) for m in config]

    def get_channels(self):
        return [m.channel for m in self.mappings]

    def get_channels_for_electrodes(self, electrodes):
        return [m.channel for m in self.mappings if m.electrode in electrodes]

    class Mapping:
        def __init__(self, channel, electrode, x, y):
            self.channel = int(channel)
            self.electrode = int(electrode)
            self.x = float(x)
            self.y = float(y)


def electrode_rectangle_indices(xmin, ymin, xmax, ymax):
    return [220 * y + x for y in range(max(ymin, 0), min(ymax + 1, 120)) for x in range(max(xmin, 0), min(xmax + 1, 220))]


def electrode_rectangle_um(xmin, ymin, xmax, ymax):
    return electrode_rectangle_indices(int(xmin / 17.5), int(ymin / 17.5), int(xmax / 17.5), int(ymax / 17.5))
