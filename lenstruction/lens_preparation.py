from __future__ import print_function

import numpy as np
import lenstronomy.Util.util as util
from lenstronomy.LensModel.numeric_lens_differentials import NumericLens
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Cosmo.lens_cosmo import LensCosmo


class LensPreparation(object):
    """
    This class is used to do initialize kwargs of the lens models according to the given lens_model_list.
    (It's easy and free to choose lens_model_list, but a little tricky to generate corresponding kwargs.)
    For the situation that you know deflection maps, then you want know shear, convergency, and shift ....
    Or, you know, shear, convergency, shift, then you wanna know deflection maps.

    """
    def __init__(self, zl, zs, xdeflection=None,ydeflection=None, gamma1=None,gamma2=None, kappa=None, lenseq='1'):
       """

       input the initial lens model.
       :param zl: redshift of the lens
       :param zs: redshift of the source
       :param xdeflection: the x deflection map
       :param ydeflection: the y deflection map
       :param gamma1: gamma1 map of shear
       :param gamma2: gamma2 map of shear
       :param kappa: convergence map
       :param lenseq: '1'/'-1', corresponding to lens equation beta = theta +/- alpha
       """
       cosmo = LensCosmo(zl, zs)
       dlsds = cosmo.D_ds / cosmo.D_s
       #alphax
       if xdeflection is None:
           self.alphax = None
       elif lenseq == '1':
           self.alphax = xdeflection* dlsds
       elif lenseq =='-1':
           self.alphax = -xdeflection * dlsds
       #alphay
       if ydeflection is None:
           self.alphay = None
       elif lenseq == '1':
           self.alphay = ydeflection* dlsds
       elif lenseq =='-1':
           self.alphay = -ydeflection * dlsds
       #gamma1
       if gamma1 is None:
           self.gamma1 = None
       elif lenseq == '1':
           self.gamma1 = gamma1 * dlsds
       elif lenseq == '-1':
           self.gamma1 = -gamma1 * dlsds
       #gamma2
       if gamma2 is None:
           self.gamma2 = None
       elif lenseq == '1':
           self.gamma2 = gamma2 * dlsds
       elif lenseq == '-1':
           self.gamma2 = -gamma2 * dlsds
       #kappa
       if kappa is None:
           self.kappa = None
       elif lenseq == '1':
           self.kappa = kappa * dlsds
       elif lenseq == '-1':
           self.kappa = -kappa * dlsds



    def kwargs_lens_configuration(self,ximg_list, yimg_list, kwargs_data_joint, fixed_index=None,
                                  lens_model_list= ['SHIFT','SHEAR','CONVERGENCE','FLEXIONFG'], diff = 0.1):
        """

        :param ximg_list: list, x cooridinate of lensed image
        :param yimg_list: list, y cooridinate of lensed image
        :param kwargs_data_joint: list, image data arguments
        :param fixed_index: int, fixed image index
        :param lens_model_list: string list, name of lens model
        :param diff: float, scale of derivation
        :return:
        """

        if fixed_index is not None:
            img_index = fixed_index
        else:
            img_index = 0
        self.img_index = img_index
        kwagrs_lens_list = []
        kwargs_lens_init = []
        for i in range(len(ximg_list)):
             kwargs_data = kwargs_data_joint['multi_band_list'][i][0]
             kwargs_lens = self.lens_param_initial(x= ximg_list[i], y = yimg_list[i], kwargs_data =kwargs_data,
                                                   lens_model_list=lens_model_list, diff=diff)
             kwagrs_lens_list.append(kwargs_lens)
             kwargs_lens_init += (kwargs_lens)
        kwargs_lens_sigma_tmp = []
        kwargs_fixed_lens_tmp = []
        kwargs_lower_lens_tmp = []
        kwargs_upper_lens_tmp = []
        for i in range(len(ximg_list)):
            for lens_type in lens_model_list:
                if lens_type == 'SHIFT':
                    kwargs_lens_sigma_tmp.append({'alpha_x':0.03,'alpha_y':0.03})
                    if i==img_index:
                        kwargs_fixed_lens_tmp.append(kwagrs_lens_list[img_index][0])
                    else:
                        kwargs_fixed_lens_tmp.append({'alpha_x':kwagrs_lens_list[i][0]['alpha_x'],'alpha_y':kwagrs_lens_list[i][0]['alpha_y']})
                    kwargs_lower_lens_tmp.append({'alpha_x':kwagrs_lens_list[i][0]['alpha_x']-1,'alpha_y':kwagrs_lens_list[i][0]['alpha_y']-1})
                    kwargs_upper_lens_tmp.append({'alpha_x':kwagrs_lens_list[i][0]['alpha_x']+1,'alpha_y':kwagrs_lens_list[i][0]['alpha_y']+1})
                elif lens_type == 'SHEAR':
                    kwargs_lens_sigma_tmp.append({'e1':0.1,'e2':0.1})
                    if i==img_index:
                        kwargs_fixed_lens_tmp.append(kwagrs_lens_list[img_index][1])
                    else:
                        kwargs_fixed_lens_tmp.append({'ra_0': kwagrs_lens_list[i][1]['ra_0'], 'dec_0': kwagrs_lens_list[i][1]['dec_0']})
                    kwargs_lower_lens_tmp.append({'e1':-1,'e2':-1})
                    kwargs_upper_lens_tmp.append({'e1':1,'e2':1})
                elif lens_type == 'CONVERGENCE':
                    kwargs_lens_sigma_tmp.append({'kappa_ext':0.1})
                    if i==img_index:
                        kwargs_fixed_lens_tmp.append(kwagrs_lens_list[img_index][2])
                    else:
                        kwargs_fixed_lens_tmp.append({'ra_0': kwagrs_lens_list[i][2]['ra_0'],'dec_0': kwagrs_lens_list[i][1]['dec_0']})
                    kwargs_lower_lens_tmp.append({'kappa_ext':-1})
                    kwargs_upper_lens_tmp.append({'kappa_ext':1})
                elif lens_type == 'FLEXIONFG':
                    kwargs_lens_sigma_tmp.append({'F1':0.01,'F2':0.01,'G1':0.01,'G2':0.01})
                    kwargs_fixed_lens_tmp.append({'ra_0': kwagrs_lens_list[i][3]['ra_0'],'dec_0': kwagrs_lens_list[i][3]['dec_0'],
                                         'F1':0,'F2':0,'G1':0,'G2':0})
                    kwargs_lower_lens_tmp.append({'F1':-1,'F2':-1,'G1':-1,'G2':-1})
                    kwargs_upper_lens_tmp.append({'F1':1,'F2':1,'G1':1,'G2':1})
        kwargs_lens_sigma=kwargs_lens_sigma_tmp
        kwargs_lower_lens=kwargs_lower_lens_tmp
        kwargs_upper_lens=kwargs_upper_lens_tmp
        kwargs_fixed_lens=kwargs_fixed_lens_tmp
        self.lens_model_list = lens_model_list
        self.num_img = len(ximg_list)
        lens_params= [kwargs_lens_init, kwargs_lens_sigma, kwargs_fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        return lens_params



    def lens_param_initial(self, x, y, kwargs_data, lens_model_list, diff):
        """
        Initial lens parameters, deduced from deflection maps or read directly
        :param x: x coordinate of lensed image
        :param y: y coordinate of lensed image
        :param kwargs_data: arguments of image data
        :param lens_model_list: list of choices of lens model
        :param diff: scale of derivation
        :return: list of lens parameters
        :return:
        """
        if self.alphax is not None and self.alphay is not None:
            kwargs_lens = self._cal_kwargs_lens(x=x, y=y, kwargs_data=kwargs_data, diff=diff, lens_model_list=lens_model_list)
        else:
            kwargs_lens = self._read_kwargs_lens(x=x, y=y, kwargs_data=kwargs_data, lens_model_list=lens_model_list)
        return kwargs_lens



    def _read_kwargs_lens(self,x,y, kwargs_data,lens_model_list=['SHIFT','SHEAR','CONVERGENCE','FLEXIONFG']):
        """
        Initial lens parameters, read from shear, convergence
        :param x: x coordinate of lensed image
        :param y: y coordinate of lensed image
        :param kwargs_data: arguments of image data
        :param lens_model_list: list of choices of lens model
        :return: list of lens parameters
        """
        imageData = ImageData(**kwargs_data)
        r_cut = (np.shape(imageData.data)[0] - 1) / 2
        xaxes, yaxes = imageData.pixel_coordinates
        ra_center = xaxes[int(r_cut + 1), int(r_cut + 1)]
        dec_center = yaxes[int(r_cut + 1), int(r_cut + 1)]
        kwargs_lens = []
        for lens_type in lens_model_list:
            if lens_type == 'SHIFT':
                kwargs_lens.append({'alpha_x': ra_center, 'alpha_y': dec_center})
            elif lens_type=='SHEAR':
                gamma1 = self.gamma1[x,y]
                gamma2 = self.gamma2[x,y]
                kwargs_lens.append({'e1': gamma1, 'e2': gamma2, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type =='CONVERGENCE':
                kappa = self.kappa[x,y]
                kwargs_lens.append({'kappa_ext': kappa, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'FLEXIONFG':
                kwargs_lens.append(
                    {'F1': 0, 'F2': 0, 'G1': 0, 'G2': 0, 'ra_0': ra_center, 'dec_0': dec_center})
        return kwargs_lens

    def _cal_kwargs_lens(self,x, y, kwargs_data, lens_model_list, diff= 0.01,):
        """
        Initial lens parameters, deduced from deflection maps.
        :param x: x coordinate of lensed image
        :param y: y coordinate of lensed image
        :param kwargs_data: arguments of image data
        :param lens_model_list: list of choices of lens model
        :param diff: scale of derivation
        :return: list of lens parameters
        """
        imageData = ImageData(**kwargs_data)
        r_cut = (np.shape(imageData.data)[0] - 1) / 2
        alphax = self.alphax[int(x - r_cut):int(x + r_cut + 1), int(y - r_cut):int(y + r_cut + 1)]
        alphay = self.alphay[int(x - r_cut):int(x + r_cut + 1), int(y - r_cut):int(y + r_cut + 1)]
        xaxes, yaxes = imageData.pixel_coordinates
        ra_center = xaxes[int(r_cut + 1), int(r_cut + 1)]
        dec_center = yaxes[int(r_cut + 1), int(r_cut + 1)]
        kwargs_lens_in = [{'grid_interp_x': xaxes[0], 'grid_interp_y': yaxes[:, 0], 'f_x': alphax,
                           'f_y': alphay}]

        kwargs_lens=[]
        for lens_type in lens_model_list:
            if lens_type=='INTERPOL':
                kwargs_lens.append({'grid_interp_x': xaxes[0], 'grid_interp_y': yaxes[:,0], 'f_x': alphax, 'f_y': alphay})
            elif lens_type=='SHIFT':
                kwargs_lens.append({'alpha_x': ra_center , 'alpha_y': dec_center })
            elif lens_type == 'SHEAR':
                gamma1, gamma2 = NumericLens(['INTERPOL']).gamma(util.image2array(xaxes),
                                                                 util.image2array(yaxes), kwargs=kwargs_lens_in,diff = diff)
                gamma_1_center, gamma_2_center = gamma1.mean(), gamma2.mean()
                kwargs_lens.append({'e1': gamma_1_center, 'e2': gamma_2_center, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'CONVERGENCE':
                kappa = NumericLens(['INTERPOL']).kappa(util.image2array(xaxes),
                                                        util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                kappa_center = kappa.mean()
                kwargs_lens.append({'kappa_ext': kappa_center, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'FLEXION':
                g1, g2, g3, g4 = NumericLens(['INTERPOL']).flexion(util.image2array(xaxes),
                                                                   util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                g1_c, g2_c, g3_c, g4_c = g1.mean(), g2.mean(), g3.mean(), g4.mean()
                kwargs_lens.append({'g1': g1_c,'g2':g2_c,'g3':g3_c,'g4':g4_c, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'FLEXIONFG':
                g1, g2, g3, g4 = NumericLens(['INTERPOL']).flexion(util.image2array(xaxes),
                                                                   util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                g1_c, g2_c, g3_c, g4_c = g1.mean(), g2.mean(), g3.mean(), g4.mean()
                F1_c = (g1_c + g3_c) * 0.5
                F2_c = (g2_c + g4_c) * 0.5
                G1_c = (g1_c - g3_c) * 0.5 - g3_c
                G2_c = (g2_c - g4_c) * 0.5 - g4_c
                kwargs_lens.append({'F1': F1_c, 'F2': F2_c, 'G1': G1_c, 'G2': G2_c, 'ra_0': ra_center, 'dec_0': dec_center})
        return kwargs_lens




    def model_index_list(self,num_img=None, lens_model_list=None):
        """
        :param num_img: number of lensed image
        :param lens_model_list: lens model list
        :return:
        """
        lensmodel_list = []
        index_lens_model_list = []
        if num_img is None:
            num_img = self.num_img
        if lens_model_list is None:
            lens_model_list = self.lens_model_list
        for i in range(num_img):
            lensmodel_list += lens_model_list
            index_lens_model_list.append([i * 4, i * 4 + 1, i * 4 + 2, i * 4 + 3])
            #TODO if the lens model list is not lens=4 anymore
            lens_model_index_list = {'lens_model_list': lensmodel_list,'index_lens_model_list': index_lens_model_list}
        return lens_model_index_list



    #TODO isolate parameters sigma, lower and upper boundary setting
    def _lens_param_fixed(self, lens_param_fixed=None):
        lens_param_fixed = lens_param_fixed
        return lens_param_fixed

    def _lens_param_upper(self, lens_param_upper=None):
        lens_param_upper = lens_param_upper
        return lens_param_upper

    def _lens_param_lower(self, lens_param_lower=None):
        lens_param_lower = lens_param_lower
        return lens_param_lower



