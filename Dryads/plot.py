import warnings

#
import numpy as np

#
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as colors
import matplotlib.collections as collections
import styles

#
from astrodendro import Dendrogram
from astrodendro.plot import DendrogramPlotter
import astropy.units as u
import astropy.constants as c

#
from tools import *


class DryadsPlotter(object):

    def __init__(self, dryads):
        self.dryads = dryads

        self.DendroPlotter = DendrogramPlotter(dryads.dendro)

    def plot_tree(self, ax, feature, cmap = 'RdYlBu_r', cscale = 'linear', norm = None, *args, **kargs):

        # feature to plot
        self.feature, self.data_feature, self.unit_feature = feature_builder(self.dryads, feature)


        # color map
        self.cmap = plt.cm.get_cmap(cmap)  ## push the checking to plt.cm
        if isinstance(norm, colors.Normalize):
            self.norm = norm
            self.lognorm = True if type(norm) == colors.LogNorm else False

            warnings.warn('Custom norm; cscale is not used.')
        else:
            if isinstance(cscale, str) and (cscale.lower() == 'linear'):
                self.norm = colors.Normalize(vmin = np.nanpercentile(self.data_feature, 2.5),
                                             vmax = np.nanpercentile(self.data_feature, 97.5))
                self.lognorm = False
            elif isinstance(cscale, str) and (cscale.lower() == 'log'):
                self.norm = colors.LogNorm(vmin = np.nanpercentile(self.data_feature, 2.5),
                                           vmax = np.nanpercentile(self.data_feature, 97.5))
                self.lognorm = True
            else:
                raise ValueError('Choose either "linear" or "log" for normalization.')


        # dendro plot
        self.DendroPlotter.plot_tree(ax, *args, **kargs)

        # dryads plot
        for i in range(len(self.dryads.dendro)):
            vertices = self.DendroPlotter.get_lines(structures = i)\
                                         .get_paths()[0]\
                                         .vertices

            # plot only contours above the lower threshold
            if self.lognorm:
                cmask = (self.dryads.cbins >= vertices[0, 1])&\
                        (np.isfinite(np.log(self.data_feature[i])))
            else:
                cmask = (self.dryads.cbins >= vertices[0, 1])&\
                        (np.isfinite(self.data_feature[i]))
            cbins = self.dryads.cbins[cmask]
            cfeature = self.data_feature[i][cmask]


            # patch list
            #patchl = [patches.Rectangle((vertices[0, 0]-.3, cbins[j]),
            #                            .6, cbins[j+1]-cbins[j],
            #                            facecolor = self.cmap(self.norm(cfeature[j])),
            #                            edgecolor = 'none')\
            #          if j != (len(cbins)-1)\
            #          else\
            #          patches.Rectangle((vertices[0, 0]-.3, cbins[j]),
            #                            .6, vertices[1, 1]-cbins[j],
            #                            facecolor = self.cmap(self.norm(cfeature[j])),
            #                            edgecolor = 'none')\
            #          for j in range(len(cbins))]
            #collection = collections.PatchCollection(patchl)

            # plot
            for j in range(len(cbins)):
                if j != (len(cbins)-1):
                    patch = patches.Rectangle((vertices[0, 0]-.3, cbins[j]),
                                              .6, cbins[j+1]-cbins[j],
                                              facecolor = self.cmap(self.norm(cfeature[j])),
                                              edgecolor = 'none')
                else:
                    patch = patches.Rectangle((vertices[0, 0]-.3, cbins[j]),
                                              .6, vertices[1, 1]-cbins[j],
                                              facecolor = self.cmap(self.norm(cfeature[j])),
                                              edgecolor = 'none')

                ax.add_patch(patch)


    def plot_scatter(self, ax, xfeature, yfeature,\
    cmap = 'RdYlBu_r', cscale = 'linear',  cfeature = 'virial', norm = None, *args, **kargs):

        # feature to plot
        ## x
        self.xfeature, self.data_xfeature, self.unit_xfeature = feature_builder(self.dryads, xfeature)
        ## x
        self.yfeature, self.data_yfeature, self.unit_yfeature = feature_builder(self.dryads, yfeature)
        ## c for color-coding
        self.cfeature, self.data_cfeature, self.unit_cfeature = feature_builder(self.dryads, cfeature)
        self.leading_cfeature = np.array([self.data_cfeature[i][(self.dryads.cbins >= self.dryads.dendro[i].vmin)&\
                                                                np.isfinite(self.data_cfeature[i])][0]\
                                          if sum((self.dryads.cbins >= self.dryads.dendro[i].vmin)&\
                                                 np.isfinite(self.data_cfeature[i])) >= 1.\
                                          else np.nan\
                                          for i in range(len(self.data_cfeature))])

        # color map
        self.cmap = plt.cm.get_cmap(cmap)  ## push the checking to plt.cm
        if isinstance(norm, colors.Normalize):
            self.norm = norm
            self.lognorm = True if type(norm) == colors.LogNorm else False

            warnings.warn('Custom norm; cscale is not used.')
        else:
            if isinstance(cscale, str) and (cscale.lower() == 'linear'):
                self.norm = colors.Normalize(vmin = np.nanpercentile(self.leading_cfeature, 2.5),
                                             vmax = np.nanpercentile(self.leading_cfeature, 97.5))
                self.lognorm = False
            elif isinstance(cscale, str) and (cscale.lower() == 'log'):
                self.norm = colors.LogNorm(vmin = np.nanpercentile(self.leading_cfeature, 2.5),
                                           vmax = np.nanpercentile(self.leading_cfeature, 97.5))
                self.lognorm = True
            else:
                raise ValueError('Choose either "linear" or "log" for normalization.')


        # dryads plot
        for i in range(len(self.dryads.dendro)):
            vmin = self.dryads.dendro[i].vmin

            # plot only contours above the lower threshold
            if self.lognorm:
                mask = (self.dryads.cbins >= vmin)&\
                       (np.isfinite(np.log(self.data_xfeature[i])))&\
                       (np.isfinite(np.log(self.data_yfeature[i])))&\
                       (self.dryads.nobs[i] >= 4.)
            else:
                mask = (self.dryads.cbins >= vmin)&\
                       (np.isfinite(self.data_xfeature[i]))&\
                       (np.isfinite(self.data_yfeature[i]))&\
                       (self.dryads.nobs[i] >= 4.)

            if sum(mask) >= 1.:
                # features to plot
                xfeature = self.data_xfeature[i][mask]
                yfeature = self.data_yfeature[i][mask]
                cfeature = self.leading_cfeature[i]

                # plot
                ax.plot(xfeature, yfeature,
                        linestyle = '-',
                        marker = 'o',
                        color = self.cmap(self.norm(cfeature)),
                        linewidth = 3.,
                        markerfacecolor = self.cmap(self.norm(cfeature)),
                        markeredgecolor = 'none',
                        markersize = 5.)
            else:
                continue
