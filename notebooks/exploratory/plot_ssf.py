#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
from numpy import *
from scipy.interpolate import UnivariateSpline
import scipy.odr as odr
import matplotlib
import math
import os.path
from sys import stderr

# define and parse command line arguments
parser = argparse.ArgumentParser(prog='plot_ssf.py')
parser.add_argument('--temperature', type=float, required=True, help='temperature')
parser.add_argument('--cutoff', type=float, required=True, help='potential cutoff')
parser.add_argument('--data', metavar='FILENAME', required=True, help='data file with physical properties of interface system')
parser.add_argument('--particles', type=int, help='filter on particle number')
parser.add_argument('--box-width', type=float, help='filter on box width')
parser.add_argument('--box-height', type=float, help='filter on box height')
parser.add_argument('--dump', metavar='FILENAME', help='dump plot data for γ(q) to file')
parser.add_argument('--qmin', default=5e-2, type=float, help='lower limit on wavenumber range')
parser.add_argument('--qmax', default=4, type=float, help='upper limit on wavenumber range for interface structure factor')
parser.add_argument('--ssf-range', type=float, nargs=2, help='plot range for static structure factors')
parser.add_argument('--num-err', default=5, type=int, help='number of error bars')
parser.add_argument('--fit-limit', default=.5, type=float, help='maximum wavenumber for bulk fits')
parser.add_argument('--periodic-bulk', action='store_true', help='use bulk model with periodic boundaries')
parser.add_argument('--adjust-slab', type=float, help='adjust width of the liquid slab by this amount')
parser.add_argument('--fit-gamma', default=(.1, .5), type=float, nargs=2, help='min and max wavenumbers for fit to γ(q)')
parser.add_argument('--table', action='store_true', help='output fit results in table format')
parser.add_argument('--output', nargs='+', help='output file(s)')
parser.add_argument('-v', '--verbose', action='store_true', help='verbose output')
parser.add_argument('interface', metavar='INTERFACE', help='structure factor of interface')
parser.add_argument('liquid', metavar='LIQUID', help='structure factor of liquid phase')
parser.add_argument('vapour', metavar='VAPOUR', help='structure factor of vapour phase')
args = parser.parse_args()

# support plotting without X11 display
if args.output:
    matplotlib.use('Agg')
#else:
#    matplotlib.use('GtkAgg')

# matplotlib backend must be selected before importing pylab
from pylab import *

dashed_line = [3, 2]
dotted_line = [1, 1]

def ornstein_zernike(params, q, rho, temp):
    """ Static structure factor in Ornstein-Zernike approximation """
    kappa, xi = params
    return rho * temp * kappa / (1 + (q * xi)**2)

def parabola(params, q, gamma_0):
    """ Small-wavenumber form of the effective surface tension up to second order """
    ell, = params
    return gamma_0 * (1 + (q * ell)**2)

def main(args):
    import defaults; defaults.set_plot_defaults()
    from compute_gid_bulk import compute_gid
    from utility import split_blocks

    if args.output:
        axes((.15, .15, .8, .8))
    else:
        rc('figure', figsize=(8, 6), dpi=80)
        rc('legend', labelspacing=0, fontsize=11)

    # read physical properties of interface system
    data = loadtxt(args.data)
    if args.verbose:
        print >>stderr, "{0} data sets found in file '{1}'".format(data.shape[0], args.data)

    # select temperature and cutoff
    data = data[where(data[:,2] == args.cutoff)]
    data = data[where(data[:,3] == args.temperature)]
    if args.verbose:
        print >>stderr, "select cutoff {0} and temperature {1}".format(args.cutoff, args.temperature)

    # apply filters
    idx = arange(data.shape[0])
    if args.particles:
        idx = intersect1d(idx, where(data[:, 1] == args.particles)[0])
    if args.box_width:
        idx = intersect1d(idx, where(data[:, 9] == args.box_width)[0])
    if args.box_height:
        idx = intersect1d(idx, where(data[:, 10] == args.box_height)[0])

    if args.verbose:
        print >>stderr
        print >>stderr, "  particles   box width   box height"
        for i,x in enumerate(data[:]):
            flag = (i in idx) and "*" or " "
            print >>stderr, "{0} {1:<12g}{2:<12g}{3:<12g}".format(flag, x[1], x[9], x[10])
        print >>stderr

    if len(idx) == 0:
        raise SystemExit("no data sets left after filtering");
    elif len(idx) > 1:
        raise SystemExit("too many data sets left, apply more stringent filters")
    data = data[idx[0]]

    # assign names to some columns of the selected data row
    particles = data[1]
    cutoff = data[2]
    temperature = data[3]
    box_width = data[9]
    box_height = data[10]
    density = particles / (pow(box_width, 2) * box_height)
    rho_liquid = data[14]
    rho_liquid_err = data[15]
    rho_vapour = data[16]
    rho_vapour_err = data[17]
    gamma_0 = data[20]
    gamma_0_err = data[21]
    # compute slab width from ρ = φ ρ_l + (1-φ) ρ_v ⇒ w = φ L_z
    # should agree with (data[12] ± data[13])
    # The formula implicitly places the interface at the position of the Gibbs dividing surface
    slab_width = (density - rho_vapour) / (rho_liquid - rho_vapour) * box_height

    adjust_slab = args.adjust_slab or None
    if adjust_slab:
        slab_width += adjust_slab
        print "Adjust slab width by δL = {0:g}".format(args.adjust_slab)

    if args.verbose and not args.table:
        print "Temperature: {0:g}".format(temperature)
        print "Potential cutoff: {0:g}".format(cutoff)
        print "Box size: {0:g} {1:g} {2:g}".format(box_width, box_width, box_height)
        print "Liquid density: {0:g}".format(rho_liquid)
        print "Vapour density: {0:g}".format(rho_vapour)
        print "Interfacial tension: {0:g}".format(gamma_0)
    elif args.table:
        print '# temperature  cutoff  box_width  slab_width  rho_liquid  err  rho_vapour  err  kappa_liquid  err  xi_liquid  err  kappa_vapour  err  xi_vapour  err  gamma_0  err  ell err'
        print '# data files: {0}'.format(args.data)
        print '# structure factors: {0} {1} {2}'.format(args.interface, args.liquid, args.vapour)
        print '# fit limit S(q): {0}, fit limits γ(q): {1} {2}'.format(args.fit_limit, *args.fit_gamma)
        print '{0:g}  {1:g}  {2:g}  {3:.4g}'.format(temperature, cutoff, box_width, slab_width),
        print '{0:.5g} {1:.2g} '.format(rho_liquid, rho_liquid_err),
        print '{0:.5g} {1:.2g} '.format(rho_vapour, rho_vapour_err),

    # load ssf data
    ssf_interface = split_blocks(loadtxt(args.interface))
    ssf_liquid = split_blocks(loadtxt(args.liquid))
    ssf_vapour = split_blocks(loadtxt(args.vapour))

    ssf_interface = array(ssf_interface) # assume equal length of each block
    ssf_interface = concatenate((
        [ssf_interface[0, :, 0]]
      , [mean(ssf_interface[..., 1], axis=0)]
      , [sqrt(sum(pow(ssf_interface[..., 2], 2), axis=0) / ssf_interface.shape[0])]
    ), axis=0).T

    ssf_liquid = array(ssf_liquid) # assume equal length of each block
    ssf_liquid = concatenate((
        [ssf_liquid[0, :, 0]]
      , [mean(ssf_liquid[..., 1], axis=0)]
      , [sqrt(sum(pow(ssf_liquid[..., 2], 2), axis=0) / ssf_liquid.shape[0])]
    ), axis=0).T

    ssf_vapour = array(ssf_vapour) # assume equal length of each block
    ssf_vapour = concatenate((
        [ssf_vapour[0, :, 0]]
      , [mean(ssf_vapour[..., 1], axis=0)]
      , [sqrt(sum(pow(ssf_vapour[..., 2], 2), axis=0) / ssf_vapour.shape[0])]
    ), axis=0).T

    # plot interface structure factor
    x = ssf_interface[:, 0]
    y = ssf_interface[:, 1]
    label = r'$S_\text{tot}(q)$'
    if ssf_interface.shape[1] > 2:
        y_err = ssf_interface[:, 2]
        n = args.num_err
        plot(x, y, '-r', label=label)
        errorbar(x[:n], y[:n], y_err[:n], fmt='or', markeredgecolor='r', markersize=2)
    else:
        plot(x, y, '-or', markeredgecolor='r', markersize=2, label=label)

    # plot CWT singularity
    A_cw = 2 * temperature / density / box_height * (rho_liquid - rho_vapour)**2
    q = linspace(args.qmin, 3)
    h_q = A_cw / (gamma_0 * pow(q, 2))
    l, = plot(q, h_q, ':k', linewidth=.75, zorder=0)
    l.set_dashes(dotted_line)

    # plot structure factor of liquid bulk
    x = ssf_liquid[:, 0]
    y = ssf_liquid[:, 1]
    label=r'$S_\ell(k)$'
    if ssf_liquid.shape[1] > 2:
        y_err = ssf_liquid[:, 2]
        n = args.num_err
        plot(x, y, '-g', label=label, linewidth=.75, zorder=1)
        errorbar(x[:n], y[:n], y_err[:n], fmt='dg', markeredgecolor='g', markersize=2)
    else:
        plot(x, y, '-dg', markeredgecolor='g', markersize=2, label=label, linewidth=.75, zorder=1)

    # plot structure factor of vapour bulk
    x = ssf_vapour[:, 0]
    y = ssf_vapour[:, 1]
    label=r'$S_v(k)$'
    if ssf_vapour.shape[1] > 2:
        y_err = ssf_vapour[:, 2]
        n = args.num_err
        plot(x, y, '-b', label=label, linewidth=.75, zorder=1)
        errorbar(x[:n], y[:n], y_err[:n], fmt='sb', markeredgecolor='b', markersize=2)
    else:
        plot(x, y, '-sb', markeredgecolor='b', markersize=2, label=label, linewidth=.75, zorder=1)

    # fit bulk structure factors to Ornstein-Zernike form
    if args.fit_limit > 0:
        idx, = where(ssf_liquid[:, 0] <= args.fit_limit)
        kappa_liquid = mean(ssf_liquid[:, 1]) / rho_liquid / temperature # initial guess
        # result is a tuple (param, param_err, covariance_matrix)
        param, param_err = odr.odr(
            ornstein_zernike                        # fit model
          , (kappa_liquid, 1)                       # initial parameter values (kappa, xi)
          , ssf_liquid[idx, 1], ssf_liquid[idx, 0]  # data (y, x)
          , extra_args=(rho_liquid, temperature,), full_output=0
        )[:2]
        kappa_liquid, xi_liquid = abs(param)
        kappa_liquid_err, xi_liquid_err = param_err
        if args.verbose and not args.table:
            print 'Compressiblity of liquid: {0:g} ± {1:g}'.format(kappa_liquid, kappa_liquid_err)
            print 'Correlation length of liquid: {0:g} ± {1:g}'.format(xi_liquid, xi_liquid_err)
        elif args.table:
            print '{0:.5g} {1:.2g} '.format(kappa_liquid, kappa_liquid_err),
            print '{0:.5g} {1:.2g} '.format(xi_liquid, xi_liquid_err),

        x = linspace(args.qmin, 3 * args.fit_limit, num=20)
        y = ornstein_zernike((kappa_liquid, xi_liquid), x, rho_liquid, temperature)
        y_err = zeros_like(y)
        y_err[0] = y[0] * kappa_liquid_err / kappa_liquid_err
        l, = plot(x, y, ':g', linewidth=.5, zorder=1)
        l.set_dashes(dotted_line)

        # splice small-wavenumber fit with full SSF at high wavenumbers
        ssf_fit = array((x, y, y_err)).T
        ssf_fit = ssf_fit[where(x < args.fit_limit)]
        ssf_liquid = concatenate((ssf_fit, ssf_liquid[idx[-1] + 1:]), axis=0)

        idx, = where(ssf_vapour[:, 0] <= args.fit_limit)
        kappa_vapour = mean(ssf_vapour[:, 1]) / rho_vapour / temperature # initial guess
        # result is a tuple (param, param_err, covariance_matrix)
        param, param_err = odr.odr(
            ornstein_zernike                        # fit model
          , (kappa_vapour, 1)                       # initial parameter values (kappa, xi)
          , ssf_vapour[idx, 1], ssf_vapour[idx, 0]  # data (y, x)
          , extra_args=(rho_vapour, temperature,), full_output=0
        )[:2]
        kappa_vapour, xi_vapour = abs(param)
        kappa_vapour_err, xi_vapour_err = param_err
        if args.verbose and not args.table:
            print 'Compressiblity of vapour: {0:g} ± {1:g}'.format(kappa_vapour, kappa_vapour_err)
            print 'Correlation length of vapour: {0:g} ± {1:g}'.format(xi_vapour, xi_vapour_err)
        elif args.table:
            print '{0:.5g} {1:.2g} '.format(kappa_vapour, kappa_vapour_err),
            print '{0:.5g} {1:.2g} '.format(xi_vapour, xi_vapour_err),
            print '{0:.5g} {1:.2g} '.format(gamma_0, gamma_0_err),

        x = linspace(args.qmin, 3 * args.fit_limit, num=20)
        y = ornstein_zernike((kappa_vapour, xi_vapour), x, rho_vapour, temperature)
        y_err = zeros_like(y)
        y_err[0] = y[0] * kappa_vapour_err / kappa_vapour_err
        l, = plot(x, y, ':b', linewidth=.5, zorder=1)
        l.set_dashes(dotted_line)

        # splice small-wavenumber fit with full SSF at high wavenumbers
        ssf_fit = array((x, y, y_err)).T
        ssf_fit = ssf_fit[where(x < args.fit_limit)]
        ssf_vapour = concatenate((ssf_fit, ssf_vapour[idx[-1] + 1:]), axis=0)

    # interpolate bulk structure factors
    ssf_liquid_ = UnivariateSpline(ssf_liquid[:, 0], ssf_liquid[:, 1], k=1, s=0)
    ssf_vapour_ = UnivariateSpline(ssf_vapour[:, 0], ssf_vapour[:, 1], k=1, s=0)
    if ssf_liquid.shape[1] > 2:
        ssf_liquid_err_ = UnivariateSpline(ssf_liquid[:, 0], ssf_liquid[:, 2], k=1, s=0)
    if ssf_vapour.shape[1] > 2:
        ssf_vapour_err_ = UnivariateSpline(ssf_vapour[:, 0], ssf_vapour[:, 2], k=1, s=0)

    # compute bulk structure factor of composed system, N/A = ρ L_z
    alpha = rho_liquid * slab_width / (density * box_height)
    beta = rho_vapour * (box_height - slab_width) / (density * box_height)

    x = ssf_interface[:, 0]
    if args.periodic_bulk:
        ssf_bulk = alpha * ssf_liquid_(x) + beta * ssf_vapour_(x)
    else:
        ssf_bulk = rho_liquid * compute_gid(x, ssf_liquid_, 0, slab_width, kmax=max(ssf_liquid[:, 0])) \
                 + rho_vapour * compute_gid(x, ssf_vapour_, 0, box_height - slab_width, kmax=max(ssf_vapour[:, 0]))
        ssf_bulk /= density * box_height

    ssf_bulk_ = UnivariateSpline(x, ssf_bulk, k=1, s=0)

    x = logspace(
        math.log(max(min(ssf_liquid[:, 0]), min(ssf_vapour[:, 0])), 10)
      , math.log(min(max(ssf_liquid[:, 0]), max(ssf_vapour[:, 0])), 10)
      , num=200
    )
    l, = plot(x, ssf_bulk_(x), '--', color='grey', label=r'$S_\text{b}(q)$')
    l.set_dashes(dashed_line)

    if args.verbose and not args.table:
        if args.periodic_bulk:
            R = alpha * ssf_liquid_(x[:2]) / beta / ssf_vapour_(x[:2])
        else:
            R = rho_liquid * compute_gid(x[:1], ssf_liquid_, 0, slab_width, kmax=max(ssf_liquid[:, 0])) \
              / rho_vapour / compute_gid(x[:1], ssf_vapour_, 0, box_height - slab_width, kmax=max(ssf_vapour[:, 0]))
        print 'Contributions to S_b(q→0): {0:.1f}% (liquid) {1:.1f}% (vapour)'.format(
            100 / (1 + 1 / R[0]), 100 / (1 + R[0])
        )

    # plot difference between interfacial and bulk structure factors
    x = ssf_interface[:, 0]
    y = ssf_interface[:, 1] - ssf_bulk
    idx = where(x < args.qmax)  # limit data to args.qmax
    x, y = x[idx], y[idx]
    plot(x, y, '-k', label=args.periodic_bulk and r'$\widetilde H(q)$' or r'$H(q)$')

    # xlabel(r'$q\sigma, k\sigma$')
    text(sqrt(0.1 * 1), -0.05, r'$q\sigma, k\sigma$', ha='center', va='top',
         transform=matplotlib.transforms.blended_transform_factory(gca().transData, gca().transAxes))
    ylabel(r'structure factors')

    # re-order legend items
    handles, labels = gca().get_legend_handles_labels()
    idx = ([0, 4, 3, 2, 1],)
    legend(array(handles)[idx], array(labels)[idx], loc='best', handlelength=1.5, labelspacing=0.3)

    text(.9, .9, r'$T^*={0:.2f}$'.format(temperature), ha='right', va='top', transform=gca().transAxes)
    text(.9, .8, r'$N={0:d}{{,}}{1:03d}$'.format(int(particles / 1000), int(particles % 1000))
      , ha='right', va='top', transform=gca().transAxes
    )

    axis('tight')
    xscale('log')
    yscale('log')
    xlim(args.qmin, xlim()[1])
    if args.ssf_range:
        ylim(args.ssf_range)
    else:
        ylim(.7 * ylim()[0], 1.5 * ylim()[1])
        if args.verbose:
            print >>stderr, 'axis limits: ({0:.3g}, {1:.3g}) ({2:.3g}, {3:.3g})'.format(*(xlim() + ylim()))

    if args.output:
        savefig(args.output[0], bbox_inches='tight', pad_inches=0.05)
    else:
        show()

    # start new figure
    cla()
    if 'alpha' in locals() and 'A_cw' in locals():
        x = ssf_interface[:, 0]
        S_q = ssf_interface[:, 1] - ssf_bulk
        y = A_cw / (S_q * x**2)
        gamma_max = max(y[where(x <= args.qmax)])
        if ssf_interface.shape[-1] > 2 and 'ssf_liquid_err_' in locals() and 'ssf_vapour_err_' in locals():
            S_q_err = sqrt(
                pow(ssf_interface[:, 2], 2)
              + pow(alpha * ssf_liquid_err_(x), 2)
              + pow(beta * ssf_vapour_err_(x), 2)
            )
            y_err = S_q_err / S_q * y
            plot(x, y, '-k', linewidth=.75)
            errorbar(x, y, y_err, fmt='ok', markersize=3)

            if args.verbose:
                num = 6 * args.num_err
                prec = get_printoptions()
                set_printoptions(precision=2)
                print >>stderr, "Error contributions to γ_q"
                print >>stderr, "  q:         {0}".format(ssf_interface[:num:3, 0])
                set_printoptions(precision=0)
                print >>stderr, "  interface: {0}".format(
                    100 * pow(ssf_interface[:num:3, 2] / S_q_err[:num:3], 2))
                print >>stderr, "  liquid:    {0}".format(
                    100 * pow(alpha * ssf_liquid_err_(x[:num:3]) / S_q_err[:num:3], 2))
                print >>stderr, "  vapour:    {0}".format(
                    100 * pow(beta  * ssf_vapour_err_(x[:num:3]) / S_q_err[:num:3], 2))
                set_printoptions(**prec)
        else:
            plot(x, y, '-k')

        axhline(y=gamma_0, xmin=0, xmax=.65, color='k', linewidth=.5)

        # write plot data to file
        if args.dump:
            f = open(args.dump, 'a')
            print >>f, '# temperature: {0:g}, cutoff: {1:g}, gamma_0: {2:.3g}, fit limit: {3:g}'.format(temperature, cutoff, gamma_0, args.fit_limit)
            if adjust_slab:
                print >>f, '# adjust slab width: {0:g}'.format(adjust_slab)
            print >>f, '# interface: {0} (ρ = {1:.3g})'.format(os.path.basename(args.interface), density)
            print >>f, '# liquid: {0} (ρ = {1:.3g})'.format(os.path.basename(args.liquid), rho_liquid)
            print >>f, '# vapour: {0} (ρ = {1:.3g})'.format(os.path.basename(args.vapour), rho_vapour)
            if args.periodic_bulk:
                print >>f, '# use bulk model with periodic boundaries'
            print >>f, '#\n# q   gamma(q)   gamma_err(q)'
            savetxt(f, array((x, y, y_err)).T)
            print >>f, '\n'
            f.close()

        xlabel(r'Wavenumber $q \sigma$')
        ylabel(r'Surface tension $\gamma^*(q) \propto 1 / q^2 H_q$')

        axis('tight')
        xscale('log')
        xlim(args.qmin, args.qmax)
        ylim(0, 1.2 * gamma_max)

        # fit quadratic increase to γ(q)
        # result is a tuple (param, param_err, covariance_matrix)
        idx, = where((x >= args.fit_gamma[0]) & (x <= args.fit_gamma[1]))
        param, param_err = odr.odr(
            parabola                                # fit model
          , (1,)                                    # initial parameter values (ell, )
          , y[idx], x[idx]                          # data (y, x)
          , extra_args=(gamma_0,), full_output=0
        )[:2]
        ell, = abs(param)
        ell_err, = param_err
        if args.verbose and not args.table:
            print 'Length scale ℓ: {0:g} ± {1:g}'.format(ell, ell_err)
        elif args.table:
            print '{0:.5g} {1:.2g} '.format(ell, ell_err),

        # create inset
        axes([.2, .5, .4, .4])
        gca().patch.set_alpha(0)

        plot(x, y - gamma_0, '-b', linewidth=.75)
        if 'y_err' in locals():
            errorbar(x, y - gamma_0, y_err, fmt='ob', markersize=3)

        x = linspace(args.qmin, 3 * args.fit_gamma[1], num=20)
        y = parabola((ell,), x, gamma_0)
        l, = plot(x, y - gamma_0, ':b')
        l.set_dashes(dotted_line)

        axis('tight')
        xscale('log')
        yscale('log')
        xlim(args.qmin, 3 * args.fit_gamma[1])
        ylim(1e-3, 1)

        if not args.output:
            show()
        elif len(args.output) > 1:
            savefig(args.output[1])

if __name__ == '__main__':
    main(args)
