from distutils.version import StrictVersion
from pathlib import Path
import MEAutility as mu
import numpy as np
import pickle


def import_LFPy_neuron():
    try:
        import LFPy
    except:
        raise ModuleNotFoundError("LFPy is not installed. Install it with 'pip install LFPy'")

    try:
        import neuron
    except:
        raise ModuleNotFoundError("NEURON is not installed. Install it from https://www.neuron.yale.edu/neuron/download")

    if StrictVersion(LFPy.__version__) < StrictVersion('2.1'):
        raise ImportError("LFPy version must be >= 2.1. To use a previous LFPy version, downgrade MEArec to <= 1.4.1")

    return LFPy, neuron



############ CREATE CELL MODELS ##################
def delete_cell(cell):
    for sec in cell.allseclist:
        sec = None
        del sec


def create_morph_with_axon_bifurcations(l_dend=200, d_dend=2, d_soma=10, l_axon=100, d_axon=0.3,
                                        n_bifurcations=1, l_child_axons_1=400, l_child_axons_2=600,
                                        angles_child_axons_1=5 / 4 * np.pi, angles_child_axons_2=7 / 4 * np.pi,
                                        nseg=None, zmod=False, zamp=8, zperiod=100, zseg_len=1,
                                        **cell_kwargs):
    '''
    Creates cell model with axonal bifurcations. The morphology has the z-axiz as the main axis.

    Parameters
    ----------
    l_dend: float
        Length of the dendrite (positive z-direction) in um
    d_dend: float
        Diameter of the dendrite in um
    d_soma: float
        Diameter (and length) of the soma in um
    l_axon: float
        Length of the axon (negative z-direction) in um
    d_axon: float
        Diameter of the axon in um
    n_bifurcations: int
        Number of bifurcations (only 1 is supported for the moment
    l_child_axons_1: float
        Length of the first axonal child  in um
    l_child_axons_2: float
        Length of the second axonal child  in um
    angles_child_axons_1: float
        Angle of the first axonal child in radians (0 is East)
    angles_child_axons_2: float
        Angle of the second axonal child in radians (0 is East)
    cell_kwargs: LFPy.Cell kwargs

    Returns
    -------
    cell: LFPy.Cell

    '''
    LFPy, neuron = import_LFPy_neuron()
    # create dendrite, soma, and axon
    dend, soma, axon = _create_dend_soma_axon(l_dend, d_dend, d_soma, l_axon, d_axon)

    if n_bifurcations == 1:
        assert isinstance(l_child_axons_1, (int, np.integer, float, np.float))
        assert isinstance(l_child_axons_2, (int, np.integer, float, np.float))
        assert isinstance(angles_child_axons_1, (int, np.integer, float, np.float))
        assert isinstance(angles_child_axons_2, (int, np.integer, float, np.float))

        # create bifurcation
        if zmod is False:
            axon1 = neuron.h.Section(name='childaxon1')
            neuron.h.pt3dclear(sec=axon1)
            neuron.h.pt3dadd(0, -l_axon, 0,  d_axon, sec=axon1)
            neuron.h.pt3dadd(np.cos(angles_child_axons_1) * l_child_axons_1, -l_axon +
                             np.sin(angles_child_axons_1) * l_child_axons_1, 0, d_axon, sec=axon1)

            axon2 = neuron.h.Section(name='childaxon2')
            neuron.h.pt3dclear(sec=axon2)
            neuron.h.pt3dadd(0, -l_axon, 0, d_axon, sec=axon2)
            neuron.h.pt3dadd(np.cos(angles_child_axons_2) * l_child_axons_2, -l_axon +
                             np.sin(angles_child_axons_2) * l_child_axons_2, 0, d_axon, sec=axon2)
        else:
            omega = 2 * np.pi * 1 / zperiod

            # child 1
            axon1 = neuron.h.Section(name='childaxon1')
            nsegs1 = int(l_child_axons_1 / zseg_len)
            neuron.h.pt3dclear(sec=axon1)
            neuron.h.pt3dadd(0, -l_axon, 0, d_axon, sec=axon1)
            for l in np.linspace(0, l_child_axons_1, nsegs1):
                xpos = np.cos(angles_child_axons_1) * l
                ypos = -l_axon + np.sin(angles_child_axons_1) * l
                zpos = zamp * np.sin(omega * l)
                neuron.h.pt3dadd(xpos, ypos, zpos, d_axon, sec=axon1)

            # child 2
            axon2 = neuron.h.Section(name='childaxon2')
            nsegs2 = int(l_child_axons_2 / zseg_len)
            neuron.h.pt3dclear(sec=axon2)
            neuron.h.pt3dadd(0, -l_axon, 0, d_axon, sec=axon2)

            for l in np.linspace(0, l_child_axons_2, nsegs2):
                xpos = np.cos(angles_child_axons_2) * l
                ypos = -l_axon + np.sin(angles_child_axons_2) * l
                zpos = -zamp * np.sin(omega * l)  # invert sign for second child
                neuron.h.pt3dadd(xpos, ypos, zpos, d_axon, sec=axon2)

        axon1.connect(axon, 1, 0)
        axon2.connect(axon, 1, 0)
        axon_children = [axon1, axon2]
    else:
        raise NotImplementedError

    sections = [dend, soma, axon] + axon_children

    if nseg is not None:
        for sec in sections:
            sec.nseg = nseg

    morphology = neuron.h.SectionList()
    morphology.wholetree()

    if nseg is not None:
        cell = LFPy.Cell(morphology=morphology, nsegs_method=None, pt3d=True, **cell_kwargs)
    else:
        cell = LFPy.Cell(morphology=morphology, pt3d=True, **cell_kwargs)

    return cell, sections


def create_morph_with_sin_axon(l_dend=200, d_dend=2, d_soma=10, l_axon=100, d_axon=0.3,
                               l_sin_axon=500, sin_period=200, sin_amplitude=50, nseg_sin=100, nseg=None,
                               **cell_kwargs):
    # create dendrite, soma, and axon
    LFPy, neuron = import_LFPy_neuron()
    dend, soma, axon = _create_dend_soma_axon(l_dend, d_dend, d_soma, l_axon, d_axon)

    sinaxon = neuron.h.Section(name='sinaxon')
    neuron.h.pt3dclear(sec=sinaxon)

    seglen_y = l_sin_axon / nseg_sin
    omega = 2 * np.pi * 1 / sin_period

    for i in np.arange(nseg_sin):
        ypos = i * seglen_y
        neuron.h.pt3dadd(sin_amplitude * np.sin(omega * ypos), -l_axon - ypos, 0, d_axon, sec=sinaxon)

    sinaxon.connect(axon, 1, 0)

    sections = [dend, soma, axon, sinaxon]
    if nseg is not None:
        for sec in sections:
            sec.nseg = nseg

    morphology = neuron.h.SectionList()
    morphology.wholetree()

    if nseg is not None:
        cell = LFPy.Cell(morphology=morphology, nsegs_method=None, pt3d=True, **cell_kwargs)
    else:
        cell = LFPy.Cell(morphology=morphology, pt3d=True, **cell_kwargs)

    return cell, sections


def create_morph_with_spiral_axon():
    raise NotImplementedError


def _create_dend_soma_axon(l_dend=200, d_dend=2, d_soma=10, l_axon=100, d_axon=1):
    # create dendrite, soma, and axon
    LFPy, neuron = import_LFPy_neuron()
    dend = neuron.h.Section(name='dend')
    neuron.h.pt3dclear(sec=dend)
    neuron.h.pt3dadd(0, d_soma, 0, d_dend, sec=dend)
    neuron.h.pt3dadd(0, d_soma + l_dend, 0, d_dend, sec=dend)

    soma = neuron.h.Section(name='soma')
    neuron.h.pt3dclear(sec=soma)
    neuron.h.pt3dadd(0, d_soma, 0, d_soma, sec=soma)
    neuron.h.pt3dadd(0, 0, 0, d_soma, sec=soma)

    axon = neuron.h.Section(name='axon')
    neuron.h.pt3dclear(sec=axon)
    neuron.h.pt3dadd(0, 0, 0, d_axon, sec=axon)
    neuron.h.pt3dadd(0, -l_axon, 0, d_axon, sec=axon)

    dend.connect(soma, 0, 0)
    axon.connect(soma, 1, 0)

    return dend, soma, axon


def insert_biophysics(cell, params_dict):
    for sec in cell.allseclist:
        if 'soma' in sec.name():
            sec.insert('na')
            sec.gbar_na = params_dict['na_soma']
            sec.insert('Kv1')
            sec.gbar_Kv1 = params_dict['Kv1_soma']

        elif 'axon' in sec.name():
            sec.insert('nax')
            sec.gbar_nax = params_dict['na_axon']
            sec.insert('Kv1')
            sec.gbar_Kv1 = params_dict['Kv1_axon']

        if 'dend' not in sec.name() and 'apic' not in sec.name():
            sec.ena = params_dict['ena']
            sec.ek = params_dict['ek']

        # passive
        sec.insert('pas')
        sec.g_pas = params_dict['g_pas']
        sec.e_pas = params_dict['e_pas']
        sec.Ra = params_dict['ra']
        sec.cm = params_dict['cm']


def insert_simple_biophysics(cell):
    for sec in cell.allseclist:
        if 'soma' in sec.name() or 'axon' in sec.name():
            sec.insert('hh')
        else:
            sec.insert('pas')


def get_default_biophysics_params(simple_biophysics=False):
    if simple_biophysics:
        v_init = -65  # important for a correct initialization of the Hodgkin-Huxley dynamics
        celsius = None
    else:
        v_init = -85
        celsius = 33

    rm = 15000
    params_dict = {
        'Kv1_soma': 100,
        'Kv1_axon': 400,

        'na_axon': 500,
        'na_soma': 500,

        'rm': rm,
        'ra': 80,
        'cm': 1,
        'celsius': celsius,
        'v_init': v_init,
        'e_pas': v_init,
        'g_pas': 1. / rm,

        'ena': 55,
        'ek': -98,
    }
    return params_dict


############ PROBES and MORPHOLOGY ##############
def create_mea_probe(pitch=17.5, dim=100, elec_size=5, z_offset=10):
    mea_dim = dim  # n rows x n cols
    mea_pitch = pitch  # rows and cols pitch

    mea_info = {'dim': mea_dim,
                'electrode_name': 'hd-mea',
                'pitch': mea_pitch,
                'shape': 'square',
                'size': elec_size,
                'type': 'mea',
                'plane': 'xy'}

    hdmea = mu.return_mea(info=mea_info)
    # Move the MEA out of the neuron plane (yz)
    hdmea.move([0, 0, z_offset])

    return hdmea


def planarize_swc(swc_input_file, swc_output_file=None, planar_dimension="z", span_um=10):
    swc_dtype = np.dtype([('id', np.int32, 1), ('type', np.int16, 1),
                          ('x', np.float32, 1), ('y', np.float32, 1), ('z', np.float32, 1),
                          ('radius', np.float32, 1), ('parent', np.int32, 1)])

    swc_input_file = Path(swc_input_file)
    swc_input = np.loadtxt(swc_input_file, dtype=swc_dtype)

    if span_um == 0:
        swc_input[planar_dimension] = 0
    else:
        span_original = np.ptp(swc_input[planar_dimension])
        swc_input[planar_dimension] = swc_input[planar_dimension] / span_original * span_um

    if swc_output_file is None:
        swc_output_file = swc_input_file.parent / f"{swc_input_file.stem}_planar_{planar_dimension}_span_{span_um}.swc"

    np.savetxt(swc_output_file, swc_input)

    return swc_output_file


def center_cell_xy(cell):
    # center in the xy plane
    xrange = np.ptp(cell.x)
    yrange = np.ptp(cell.y)

    x_shift = - xrange / 2 - np.min(cell.x)
    y_shift = - yrange / 2 - np.min(cell.y)
    cell.set_pos(x_shift, y_shift, 0)


############ IO ##################
def save_cell(cell, cell_name, save_folder):
    '''

    Parameters
    ----------
    cell
    cell_name
    save_path

    Returns
    -------

    '''
    save_folder = Path(save_folder)
    save_folder.mkdir(parents=True, exist_ok=True)

    cell_path = save_folder / f"{cell_name}.pkl"
    sec_path = save_folder / f"{cell_name}_sections.npy"
    sections = []
    for idx in cell.get_idx():
        sections.append(cell.get_idx_name(idx)[1])
    sections = np.array(sections)
    cell.cellpickler(cell_path)
    np.save(sec_path, sections)

    return str(cell_path)


def load_cell(cell_path):
    '''

    Parameters
    ----------
    cell_path

    Returns
    -------

    '''
    cell_path = Path(cell_path)
    sec_path = cell_path.parent / f"{cell_path.stem}_sections.npy"

    assert cell_path.is_file(), "The specified 'cell_path' is not a file"
    assert sec_path.is_file(), "Could not find 'sections.npy' file "

    cell = pickle.load(cell_path.open('rb'))
    sections = np.load(sec_path)

    return cell, sections
