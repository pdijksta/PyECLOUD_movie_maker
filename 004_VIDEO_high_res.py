from __future__ import division, print_function
import re
import os
import argparse

import scipy.io as sio
from scipy.constants import e as qe
import pylab as pl
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

import LHCMeasurementTools.mystyle as ms

pl.close('all')

# Generic config
fontsz=12
min_min=7
plot_delta = 3 # delta for x, y limits of plots
ms.mystyle_arial(fontsz=fontsz, dist_tick_lab=10)
chamber_path = 'LHC_chm_ver.mat'

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('--figtitle', default='PyECLOUD simulation')
parser.add_argument('--folder_sim', default='.')
parser.add_argument('--passlist', type=int, nargs='*')
parser.add_argument('--lower_limit', type=int, help='In logarithmic scale to basis 10, the lower limit of the electron cloud density that is shown in the plot')
parser.add_argument('--dpi', type=int, default=200)
args = parser.parse_args()

figtitle = args.figtitle
folder_sim = args.folder_sim
if args.lower_limit:
    lim_density_manual = True
    lim_dens_min = args.lower_limit
    print('Using manually set lower density limits of %i.' % lim_dens_min)
else:
    lim_dens_min = min_min

pwd = os.getcwd()
try:
    os.chdir(folder_sim)

    if args.passlist:
        passlist = args.passlist
    else:
        all_files = os.listdir('rho_video')
        regex = re.compile('rho_pass(\d+).mat')
        matches = map(regex.match, all_files)
        matches = filter(None, matches)
        passlist = map(lambda x: int(x.group(1)), matches)
        print('Automatic passlist is %s' % passlist)

    x_beam_pos = 0.
    y_beam_pos = 0.
    color_beam = 'r'
    N_dec=1
    outp_filname='pass%d.avi'%passlist[0]


    dict_pyecltest=sio.loadmat('Pyecltest.mat')
    dict_chm = sio.loadmat(chamber_path)
    print('Using chamber %s' % chamber_path)

    Vx = np.squeeze(dict_chm['Vx'])*1e3
    Vy = np.squeeze(dict_chm['Vy'])*1e3

    Vx = list(Vx)
    Vx.append(Vx[0])
    Vx=np.array(Vx)

    Vy = list(Vy)
    Vy.append(Vy[0])
    Vy=np.array(Vy)

    t=np.squeeze(dict_pyecltest['t'])
    lam_b1 = np.squeeze(dict_pyecltest['lam_t_array'])
    dens_1 = dict_pyecltest['el_dens_at_probes'][0,:]
    cendens = np.squeeze(dict_pyecltest['cen_density'])

    i_photog=0
    for pass_ind in passlist:

        filename_rho='rho_video/rho_pass%d.mat'%pass_ind
        filename_efield='efield_video/efield_pass%d.mat'%pass_ind
        dict_ecl_video=sio.loadmat(filename_rho)
        dict_efield=sio.loadmat(filename_efield)

        xg_sc = np.squeeze(dict_ecl_video['xg_sc'])*1e3
        yg_sc = np.squeeze(dict_ecl_video['yg_sc'])*1e3

        xmin=np.min(xg_sc)
        xmax=np.max(xg_sc)
        ymin=np.min(yg_sc)
        ymax=np.max(yg_sc)

        rho_video=-dict_ecl_video['rho_video']/qe
        ex_video=dict_efield['efx_video']
        ey_video=dict_efield['efy_video']
        t_video=np.squeeze(dict_ecl_video['t_video'].real)
        b_spac=np.squeeze(dict_pyecltest['b_spac'].real)
        (nphotog,_,_)=rho_video.shape

        rho_video[rho_video==0] = 1e-20 # suppress warnings while taking the log
        lim_dens_max = int(np.log10(rho_video.max())+1)
        lim_dens = [lim_dens_min, lim_dens_max]
        print("Using automatic upper limit of %i" % lim_dens_max)

        for ii in xrange(0, nphotog, N_dec):
            if ii % 10 == 0:
                print('Pass %d %d/%d' % (pass_ind,ii,nphotog))
            ms.figure(figtitle, figsize=(4.5,4))
            imm = np.squeeze(rho_video[ii,:,:])
            imm = np.log10(imm)
            imm[imm<lim_dens[0]] = lim_dens[0]

            ax = pl.gca()
            im = ax.imshow(imm.T, cmap=None, norm=None, aspect='auto', interpolation=None,
                    alpha=None, vmin=lim_dens[0], vmax=lim_dens[1], origin='lower', extent=[xmin, xmax, ymin, ymax])
            pl.plot(Vx, Vy, 'y', linewidth=1.5)
            pl.axis('equal')
            pl.xlim(np.min(Vx)-plot_delta, np.max(Vx)+plot_delta)
            pl.ylim(np.min(Vy)-plot_delta, np.max(Vy)+plot_delta)
            pl.xlabel('x [mm]')
            pl.ylabel('y [mm]')
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.02)
            clb = pl.colorbar(im, cax=cax)
            clb.set_label('log10(e- dens.)')

            filename = str('Pass%05d_%05d' % (pass_ind,ii)) + '.png'
            pl.subplots_adjust(left=0.2, right=0.7)
            pl.savefig(filename, dpi=args.dpi)
            pl.clf()
            i_photog+=1

        out = 'movie_pass_%i.mp4' % pass_ind
        if os.path.isfile(out):
            os.remove(out)
        command = 'ffmpeg -framerate 1/0.2 -i Pass%05d_%s.png -c:v libx264 -r 30 -pix_fmt yuv420p %s >/dev/null 2>&1' %\
            (pass_ind, '%05d', out)
        exit_status = os.system(command)
        print("%s finished with status %i" % (command, exit_status))

        folderpngname = outp_filname.split('.avi')[0]+'_high_res_pngs'
        os.system('mkdir -p %s'%folderpngname)
        os.system('mv *.png '+folderpngname)
finally:
    os.chdir(pwd)

## Old command:
#command = ('mencoder',
#           'mf://*.png',
#           '-mf',
#           'type=png:w=800:h=600:fps=5',
#           '-ovc',
#           'lavc',
#           '-lavcopts',
#           'vcodec=mpeg4',
#           '-oac',
#           'copy',
#           '-o',
#           outp_filname)

## Old part of for loop
            #t_curr = t_video[ii]
            #cendens_curr = np.interp(t_curr, t, cendens)
            #lam_b1_curr = np.interp(t_curr, t, lam_b1)
            #cendens_curr = np.interp(t_curr, t, cendens)

            #pl.subplot(2,1,2)
            #mask_t_win = np.logical_and(t>t_curr-tbeam_win_length/2., t<t_curr+tbeam_win_length/2.)
            #pl.plot(c*(t[mask_t_win]-t_curr), lam_b1[mask_t_win])
            #pl.plot(-c*(t[mask_t_win]-t_curr), lam_b2[mask_t_win], 'r')
            #pl.axis('tight')
            #pl.grid('on')
            #pl.axvline(0,linestyle='--',color='k')
            #pl.xlabel('(s-s0) [m]')
            #pl.ylabel('beam prof. [p/m]')
            #pl.xlim(np.min(Vx*1e3)+1.5, np.max(Vx*1e3)+1.5)
            #pl.subplots_adjust(top=0.85,right=0.8, left=0.15, hspace=0.3, wspace=0.5)

