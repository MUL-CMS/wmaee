import numpy as np
from plotly.graph_objects import Figure, Scatter, layout
from pymatgen import Spin, Orbital

def plot_total_dos(output, efermi=False, erange=None):
    """
    Plot the total density of states from a given output.
    :param output: (VASPOutput) the instance of VASPOutput
    :param efermi: (bool) weth to rescale the x-axis such that the Fermi level is zero (default: False)
    :param erange: (float, float) the energy range to clip the plot. If None the whole range will be used (default: None)
    :return: (plotly.graph_objects.Figure) the plot
    """
    # Non spin polarized case
    dos = output.total_dos
    fig = Figure()
    spins = (Spin.up, Spin.down)
    max_val = []
    min_val = []
    for i, spin in enumerate(spins):
        density = dos.densities[spin] * spin.value
        max_val.append(np.amax(density))
        min_val.append(np.amin(density))
        energies = dos.energies
        energies = energies if not efermi else energies - dos.efermi
        scatter = Scatter(x=energies, y=density, name=spin.name)
        fig.add_trace(scatter)
        if i >= len(dos.densities) - 1:
            break

    fermi_energy = dos.efermi if not efermi else 0.0
    min_val, max_val = min(min_val), max(max_val)
    fig.update_layout(shapes=[layout.Shape(type="line", x0=fermi_energy, y0=min_val, x1=fermi_energy, y1=max_val)])
    if erange is not None:
        fig.update_xaxes(range=erange)
    return fig


def plot_projected_dos(output, combine_spins=False, efermi=False, orbitals=None, sum_density=True,
                       combine_orbitals=None, erange=None):
    """
    Plot the projected density of states
    :param output: (VASPOutput) the VASP output object
    :param combine_spins: (bool) wether to sum up the up and down spin density for spin polarized calculations (default: False)
    :param efermi: (bool) wether to rescale the x-axis such that the Fermi level is zero (default: False)
    :param orbitals: (list or tuple of str) which orbitals to display e.g ('s', 'p', 'd') or ('px', 'py'). None means all (default: None)
    :param sum_density: (bool) wether to sum up the densities for the individual orbital
    :param combine_orbitals:
    :param erange: (float, float) the energy range to clip the plot. If None the whole range will be used (default: None)
    :return: (plotly.graph_objects.Figure) the plot
    """
    pdos = output.partial_dos[-1]
    energies = output.total_dos.energies
    energies = energies if not efermi else energies - output.total_dos.efermi
    fig = Figure()
    max_val = []
    min_val = []
    sd = []

    figure_data = []
    flatten = lambda l: [item for sublist in l for item in sublist]
    for orbital, data in pdos.items():
        if orbitals is not None:
            orbital_str = [o.name if isinstance(o, Orbital) else o for o in orbitals]
            if not any([orbital.name.startswith(orb_str) for orb_str in orbital_str]):
                continue

        density = [(Spin.up, np.sum(list(data.values()), axis=0))] if combine_spins else [(spin, spin.value * den) for
                                                                                          spin, den in data.items()]
        max_val.append(np.amax([d for _, d in density]))
        min_val.append(np.amin([d for _, d in density]))

        if combine_spins:
            figure_data.append(dict(x=energies, y=density[0][1], name=orbital.name, orbital=orbital, spin=Spin.up))
            sd.append((Spin.up, density[0][1]))
        else:
            for s, d in density:
                figure_data.append(
                    dict(x=energies, y=d, name='{}-{}'.format(orbital.name, s.name), orbital=orbital, spin=s))
                sd.append((s, d))
    if combine_orbitals is not None:
        combinations = []

        for combination in combine_orbitals:
            current_combination = []
            for dict_ in figure_data:
                orbital_str = [o.name if isinstance(o, Orbital) else o for o in combination]
                orbital = dict_['orbital']
                if any([orbital.name.startswith(orb_str) for orb_str in orbital_str]):
                    current_combination.append(dict_)
            combinations.append(current_combination)
        # summ up all with the same sping
        for dat in flatten(combinations):
            for i, f in enumerate(figure_data):
                if dat['orbital'] == f['orbital'] and dat['spin'] == f['spin']:
                    figure_data.pop(i)
                    break
        combined_data = []
        for combination in combinations:
            sum_up, sum_down = [], []
            name_up, name_down = '', ''
            xaxis = None
            for orbital in [o for o in combination if o['spin'] == Spin.down]:
                sum_down.append(orbital['y'])
                name_down += orbital['orbital'].name + '-'
                xaxis = orbital['x']
            for orbital in [o for o in combination if o['spin'] == Spin.up]:
                sum_up.append(orbital['y'])
                name_up += orbital['orbital'].name + '-'
                xaxis = orbital['x']
            name_up += 'up'
            name_down += 'down'
            # If combine_spins only Spin.up is used
            if len(sum_up) > 0:
                sum_up = np.sum(sum_up, axis=0)
                combined_data.append(dict(
                    x=xaxis,
                    name=name_up,
                    y=sum_up
                ))

            if not combine_spins:
                # Check if it was spin polarized calulation at all
                if len(sum_down) > 0:
                    sum_down = np.sum(sum_down, axis=0)
                    combined_data.append(dict(
                        x=xaxis,
                        name=name_down,
                        y=sum_down
                    ))
    else:
        combined_data = []

    all_data = combined_data + figure_data
    for d in all_data:
        fig.add_trace(Scatter(x=d['x'], y=d['y'], name=d['name']))
        max_val.append(np.amax(d['y']))
        min_val.append(np.amin(d['y']))
    if sum_density:
        if combine_spins:
            sd_up = np.sum([d for s, d in sd if s == Spin.up], axis=0)
            max_val.append(np.amax(sd_up))
            min_val.append(np.amin(sd_up))
            fig.add_trace(Scatter(x=energies, y=sd_up, name='sum'))
        else:
            sd_up = np.sum([d for s, d in sd if s == Spin.up], axis=0)
            max_val.append(np.amax(sd_up))
            min_val.append(np.amin(sd_up))
            fig.add_trace(Scatter(x=energies, y=sd_up, name='sum-up'))
            if len([d for s, d in sd if s == Spin.down]) > 0:
                sd_down = np.sum([d for s, d in sd if s == Spin.down], axis=0)
                max_val.append(np.amax(sd_down))
                min_val.append(np.amin(sd_down))
                fig.add_trace(Scatter(x=energies, y=sd_down, name='sum-down'))
    fermi_energy = output.total_dos.efermi if not efermi else 0.0
    min_val, max_val = min(min_val), max(max_val)
    fig.update_layout(shapes=[layout.Shape(type="line", x0=fermi_energy, y0=min_val, x1=fermi_energy, y1=max_val)])
    if erange is not None:
        fig.update_xaxes(range=erange)
    return fig