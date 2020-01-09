from __future__ import print_function
from lenstronomy.Plots.model_plot import ModelPlot
from lenstronomy.Plots.chain_plot import plot_chain_list

import matplotlib.pyplot as plt
import numpy as np
import corner
from lenstronomy.Util.analysis_util import half_light_radius
from lenstronomy.LightModel.light_model import LightModel
import lenstronomy.Util.class_creator as class_creator
import lenstronomy.Util.util as util



def plot_chain(chain_list):
    """
    a fuction to plot chain_list of fitting results
    :param chain_list: chain_list of fitting results
    :return:plot
    """
    for i in range(len(chain_list)):
        f, axes = plot_chain_list(chain_list, index=i)
    f.show()

def plot_mcmc(chain_list_mcmc):
    """
    a fuction to plot chain_list and contour figure of mcmc results
    :param chain_list_mcmc: chain_list of mcmc results
    :return:
    """
    sampler_type, samples_mcmc, param_mcmc, dist_mcmc = chain_list_mcmc[0]
    print("number of non-linear parameters in the MCMC process: ", len(param_mcmc))
    print("parameters in order: ", param_mcmc)
    print("number of evaluations in the MCMC process: ", np.shape(samples_mcmc)[0])
    if not samples_mcmc == []:
        fig=corner.corner(samples_mcmc, labels=param_mcmc, show_titles=True)
    fig.show()

class ClsrAnalysis(object):
    """
    show modeling process, reconstructed source
    """
    def __init__(self, kwargs_data_joint, kwargs_model, multi_band_type='joint-linear'):
        self.kwargs_data_joint = kwargs_data_joint
        self.multi_band_list = self.kwargs_data_joint.get('multi_band_list', [])
        self.multi_band_type = multi_band_type
        self.kwargs_model = kwargs_model


    def plot_modeling(self,kwargs_result,
                      deltaPix_s=0.03,numPix_s=None, text_source='', data_index = 0,
                      text='sys',img_name='sys',font_size=25,scale_size=0.1):
        """
        a function to show modeling process containing data, reconstructed image, residual map,
        and reconstructed source.
        :param kwargs_result: modeling results
        :param deltaPix: pixel scale in the source plane
        :param numPix: pixel numbers in the source plane
        :param multi_band_type:string, e.g., 'joint-linear', 'single-band'
        :param text: string, label of reconstructed image
        :param text_source:string, label of reconstructed source
        :param img_name: string, label of saved images
        :param font_size: font_size
        :return:
        """
        model_plot = ModelPlot(self.multi_band_list, self.kwargs_model, kwargs_result, arrow_size=0.02, cmap_string="gist_heat",
                 multi_band_type=self.multi_band_type)
        num_bands = len(self.kwargs_data_joint['multi_band_list'])
        if num_bands >1:
            f, axes = plt.subplots(num_bands, 3, figsize=(22, 18))
        else:
            f, axes = plt.subplots(num_bands, 3, figsize=(22, 6))
        for band_index in range(num_bands):
            if num_bands >1:
                ax1 = axes[band_index, 0]
                ax2 = axes[band_index, 1]
                ax3 = axes[band_index, 2]
                img_index = band_index
            else:
                ax1 = axes[0]
                ax2 = axes[1]
                ax3 = axes[2]
                img_index = data_index
            model_plot.data_plot(ax=ax1, band_index=band_index, text='Observed'+text+ repr(img_index+1),font_size = font_size)
            model_plot.model_plot(ax=ax2, image_names=True, band_index=band_index,font_size = font_size, text='Modeled'+text+ repr(img_index+1))
            model_plot.normalized_residual_plot(ax=ax3, v_min=-6, v_max=6, band_index=band_index,font_size = font_size)
        f.savefig(img_name+'residual.pdf', bbox_inches='tight')
        if numPix_s is None:
            numPix_s = self.kwargs_data_joint['multi_band_list'][0][0]['image_data'].shape[0]
        f_s, axes_s = plt.subplots(1, 1, figsize=(9, 6) )
        model_plot.source_plot(ax=axes_s, deltaPix_source=deltaPix_s, numPix=numPix_s, band_index=band_index, scale_size=scale_size,
                               font_size=font_size, text ="Source"+text_source,  plot_scale='log', v_min =-3,
                               with_caustics=True)
        #model_plot.error_map_source_plot(ax = axes_s[1],deltaPix_source=deltaPix_s, numPix=numPix_s, font_size=font_size)
        f_s.savefig(img_name +'source.pdf')
        f_s.show()



    def source_flux_rh(self,kwargs_result,deltaPix_s,numPix):
        """
        A function to calculate flux, half light radius of the given modeling result.
        :param kwargs_results: modeling result
        :param deltaPix: pixel scale in the source plane
        :param numPix: pixel numbers in the source plane
        :return: flux (cts/s) and R_e (") in the source plane
        """
        imageModel = class_creator.create_im_sim(self.multi_band_list, multi_band_type='single-band',
                                                 kwargs_model=self.kwargs_model, bands_compute=[True], band_index=0)
        kwargs_source = kwargs_result['kwargs_source']
        _, _, _, _ = imageModel.image_linear_solve(inv_bool=True, **kwargs_result)

        x_center = kwargs_source[0]['center_x']
        y_center = kwargs_source[0]['center_y']
        x_grid_source, y_grid_source = util.make_grid(numPix=numPix, deltapix=deltaPix_s)
        x_grid_source += x_center
        y_grid_source += y_center
        source_light_model = self.kwargs_model['source_light_model_list']
        lightModel = LightModel(light_model_list=source_light_model)
        flux = lightModel.surface_brightness(x_grid_source, y_grid_source, kwargs_source) * deltaPix_s ** 2
        rh = half_light_radius(flux, x_grid_source, y_grid_source, x_center, y_center)
        return flux, rh
