"""Reproduce the nine published figures without SSH, Slurm, or local absolute paths."""
from pathlib import Path
import argparse
import json
import math
import os

from published_results import DEFAULT_DATA, load_measurements, verify_assets


def plot(data, output):
    output.parent.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault('MPLCONFIGDIR', str(output.parent / '.mpl_cache'))
    os.environ.setdefault('TEXMFVAR', str(output.parent / '.tex_cache'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    import numpy as np
    from academic_style import apply_academic_style, validate_active_style

    apply_academic_style()
    validate_active_style()
    rows = load_measurements(data)
    output.mkdir(parents=True, exist_ok=False)
    figures = []

    def save(fig, ax, name):
        ax.grid(False, which='both')
        fig.canvas.draw()
        legends = fig.legends + ([ax.get_legend()] if ax.get_legend() else [])
        for legend in legends:
            bounds = legend.get_window_extent(fig.canvas.get_renderer())
            assert bounds.x0 >= 0 and bounds.y0 >= 0
            assert bounds.x1 <= fig.bbox.width and bounds.y1 <= fig.bbox.height
        for suffix in ('png', 'pdf'):
            fig.savefig(output / (name + '.' + suffix), dpi=600)
        plt.close(fig)
        figures.append(name)

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.subplots_adjust(left=.15, right=.96, top=.92, bottom=.41)
    k = np.linspace(0, 80, 500)
    ax.plot(k, np.degrees(np.arctan(np.sqrt(k*.01/2))), color='black')
    ax.plot(k, np.degrees(np.arccos(1-k*.01/4)), color='black', ls='--')
    for row in rows:
        color, marker = ('#D62728', 'o') if row['level'] == 13 else ('#1F77B4', '^')
        ax.errorbar(row['k'], row['theta_deg'],
                    yerr=[[row['theta_deg']-row['directional_min']],
                          [row['directional_max']-row['theta_deg']]],
                    fmt=marker, color=color, ls='none', capsize=3, markeredgewidth=1.6,
                    markerfacecolor=color if row['radius_shape_gate_passed'] else 'white')
    handles = [Line2D([], [], color='black', label='Small-slope theory'),
               Line2D([], [], color='black', ls='--', label='Calibration theory'),
               Line2D([], [], color='#D62728', marker='o', mfc='white', ls='none', label='L13 data'),
               Line2D([], [], color='#1F77B4', marker='^', mfc='white', ls='none', label='L14 data')]
    ax.set(xlabel=r'Interaction strength $k$', ylabel=r'Contact angle $\theta\ ({}^\circ)$',
           xlim=(0, 80), ylim=(0, 40), title=r'$h_*=0.01$, $\gamma=1$')
    ax.legend(handles=handles, loc='upper center', bbox_to_anchor=(.5, -.28), ncol=2)
    fig.text(.5, .025, 'Filled marker: radius/shape gates passed', ha='center', fontsize=20)
    save(fig, ax, 'k_theta')

    colors = {15: '#000000', 25: '#D62728', 35: '#1F77B4'}
    for level in (13, 14):
        fig, ax = plt.subplots(figsize=(8, 6))
        fig.subplots_adjust(left=.16, right=.96, top=.91, bottom=.43)
        for row in rows:
            if row['level'] != level:
                continue
            ax.plot(row['s'], row['hx'], color=colors[row['nominal_angle_deg']],
                    ls='-' if level == 13 else '--',
                    label=r'$k=%.2f$, $t=%g$, $\theta=%.2f^\circ$' %
                          (row['k'], row['time'], row['theta_deg']))
        ax.axhline(.01, color='#777777', lw=1, ls=':')
        ax.set(xlabel=r'Coordinate $x$', ylabel=r'Height $h$', xlim=(-5, 5), ylim=(0, 1.04),
               title='Latest droplet sections, L%d' % level)
        ax.legend(loc='upper center', bbox_to_anchor=(.5, -.26), ncol=1)
        fig.text(.5, .025, 'Angles from bulk fits; foot views use equal scales',
                 ha='center', fontsize=20)
        save(fig, ax, 'sections_L%d' % level)

    for row in rows:
        fit = next(f['parabola'] for f in row['fits'] if
                   f['axis'] == 'x' and f['sign'] == 1 and f['window'] == [.4, .8])
        assert abs(fit['geometric_angle_deg']-row['theta_deg']) < .005
        rc = fit['contact_radius']
        fig, ax = plt.subplots(figsize=(8, 6))
        fig.subplots_adjust(left=.16, right=.96, top=.90, bottom=.35)
        points = [(x-rc, h-.01) for x, h in zip(row['s'], row['hx']) if x > 0]
        ax.plot(*zip(*points), color='black', label='Profile')
        x = np.linspace(-1.1, .6, 500)
        ax.plot(x, fit['intercept']+fit['quadratic_coefficient']*(x+rc)**2-.01,
                color='#D62728', ls='--', label=r'$40$--$80\%$ parabola')
        ax.plot(x, -fit['slope']*x, color='#1F77B4', ls='-.', label='Extrapolated tangent')
        ax.axhline(0, color='#777777', lw=1)
        theta = math.radians(fit['geometric_angle_deg'])
        phi = np.linspace(math.pi-theta, math.pi, 100)
        ax.plot(.55*np.cos(phi), .55*np.sin(phi), color='black', lw=1.6)
        ax.text(.06, .17, r'$\theta=%.2f^\circ$' % row['theta_deg'],
                ha='left', va='center', fontsize=20)
        ax.set(xlabel=r'$r-R_{\mathrm{fit}}$', ylabel=r'$h-h_*$',
               xlim=(-1.1, .68), ylim=(-.06, .78),
               title=r'$k=%.2f$, L%d, $t=%g$' % (row['k'], row['level'], row['time']))
        ax.set_aspect('equal', adjustable='box')
        fig.legend(*ax.get_legend_handles_labels(), loc='lower center',
                   bbox_to_anchor=(.5, .015), ncol=1)
        save(fig, ax, 'angle_foot_'+row['case'])
    (output / 'figure_manifest.json').write_text(json.dumps({'figures': figures}, indent=2))
    return figures


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=DEFAULT_DATA)
    parser.add_argument('--output', type=Path, help='New output directory; existing directories are rejected')
    parser.add_argument('--check', action='store_true', help='Verify hashes and recompute angles without plotting')
    args = parser.parse_args()
    count = verify_assets(args.data)
    rows = load_measurements(args.data)
    print('PASS: %d asset hashes; %d measured profiles' % (count, len(rows)))
    if not args.check:
        output = args.output or args.data / 'regenerated_figures'
        print('Rendered %d figures in %s' % (len(plot(args.data, output)), output))


if __name__ == '__main__':
    main()
