from __future__ import print_function

import numpy as np
import lenstronomy.Util.util as util
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Cosmo.lens_cosmo import LensCosmo
from astropy import wcs


class LensPreparation(object):
    """
    This class is used to do initialize kwargs of the lens models according to the given lens_model_list.
    (It's easy and free to choose lens_model_list, but a little tricky to generate corresponding kwargs.)
    For the situation that you know deflection maps, then you want know shear, convergency, and shift ....
    Or, you know, shear, convergency, shift, then you wanna know deflection maps.

    """
    def __init__(self, zl, zs, xdeflection=None,ydeflection=None, gamma1=None,gamma2=None, kappa=None, hdul = None):
       """

       input the initial lens model.
       :param zl: redshift of the lens
       :param zs: redshift of the source
       :param xdeflection: the x deflection map
       :param ydeflection: the y deflection map
       :param gamma1: gamma1 map of shear
       :param gamma2: gamma2 map of shear
       :param kappa: convergence map
       """
       cosmo = LensCosmo(zl, zs)
       dlsds = cosmo.D_ds / cosmo.D_s
       self.dlsds = dlsds
       self.wcs = wcs.WCS(hdul[0].header)
       #alphax
       if xdeflection is None:
           self.alphax = None
       else:
           self.alphax = xdeflection* dlsds
       #alphay
       if ydeflection is None:
           self.alphay = None
       else:
           self.alphay = ydeflection* dlsds
       #gamma1
       if gamma1 is None:
           self.gamma1 = None
       else:
           self.gamma1 = gamma1 * dlsds
       #gamma2
       if gamma2 is None:
           self.gamma2 = None
       else:
           self.gamma2 = gamma2 * dlsds
       #kappa
       if kappa is None:
           self.kappa = None
       else:
           self.kappa = kappa * dlsds




    def params(self,ximg_list, yimg_list, ra_list, dec_list,
               kwargs_data_joint,
               lens_model_list= ['SHIFT','SHEAR','CONVERGENCE','FLEXIONFG'], diff = 0.03,
               kwargs_sigma = None, kwargs_fixed = None, kwargs_lower = None, kwargs_upper = None):
        """
        lens model parameters configuration in lenstronomy keywords arguments
        :param ximg_list: list, x cooridinate of lensed image
        :param yimg_list: list, y cooridinate of lensed image
        :param kwargs_data_joint: list, image data arguments
        :param lens_model_list: string list, name of lens model
        :param diff: float (arcsec), scale of derivation
        :return: lens model keywords arguments
        """

        kwagrs_lens_list = []
        kwargs_lens_init = []
        magnification_list = []
        for i in range(len(ximg_list)):
             kwargs_data = kwargs_data_joint['multi_band_list'][i][0]
             kwargs_lens, magnification = self.lens_param_initial(x= ximg_list[i], y = yimg_list[i], ra =ra_list[i] , dec= dec_list[i],  kwargs_data =kwargs_data, lens_model_list=lens_model_list, diff=diff)
             kwagrs_lens_list.append(kwargs_lens)
             kwargs_lens_init += (kwargs_lens)
             magnification_list.append(magnification)
        magnification_list = np.abs(magnification_list)
        mag_min_index = np.where(magnification_list == np.min(magnification_list))[0][0]
        kwargs_lens_sigma_tmp = []
        kwargs_fixed_lens_tmp = []
        kwargs_lower_lens_tmp = []
        kwargs_upper_lens_tmp = []
        fixed_index = mag_min_index
        for i in range(len(ximg_list)):
            for lens_type in lens_model_list:
                if lens_type == 'SHIFT':
                    kwargs_lens_sigma_tmp.append({'alpha_x':0.03,'alpha_y':0.03})
                    if i==fixed_index:
                        kwargs_fixed_lens_tmp.append(kwagrs_lens_list[fixed_index][0])
                    else:
                        kwargs_fixed_lens_tmp.append({'alpha_x':kwagrs_lens_list[i][0]['alpha_x'],'alpha_y':kwagrs_lens_list[i][0]['alpha_y']})
                    kwargs_lower_lens_tmp.append({'alpha_x':kwagrs_lens_list[i][0]['alpha_x']-1,'alpha_y':kwagrs_lens_list[i][0]['alpha_y']-1})
                    kwargs_upper_lens_tmp.append({'alpha_x':kwagrs_lens_list[i][0]['alpha_x']+1,'alpha_y':kwagrs_lens_list[i][0]['alpha_y']+1})
                elif lens_type == 'SHEAR':
                    kwargs_lens_sigma_tmp.append({'gamma1':0.1,'gamma2':0.1})
                    if i==fixed_index:
                        kwargs_fixed_lens_tmp.append(kwagrs_lens_list[fixed_index][1])
                    else:
                        kwargs_fixed_lens_tmp.append({'ra_0': kwagrs_lens_list[i][1]['ra_0'], 'dec_0': kwagrs_lens_list[i][1]['dec_0']})
                    kwargs_lower_lens_tmp.append({'gamma1':-1,'gamma2':-1})
                    kwargs_upper_lens_tmp.append({'gamma1':1,'gamma2':1})
                elif lens_type == 'CONVERGENCE':
                    kwargs_lens_sigma_tmp.append({'kappa_ext':0.1})
                    if i==fixed_index:
                        kwargs_fixed_lens_tmp.append(kwagrs_lens_list[fixed_index][2])
                    else:
                        kwargs_fixed_lens_tmp.append({'ra_0': kwagrs_lens_list[i][2]['ra_0'],'dec_0': kwagrs_lens_list[i][1]['dec_0']})
                    kwargs_lower_lens_tmp.append({'kappa_ext':-0.2}) #
                    kwargs_upper_lens_tmp.append({'kappa_ext':3}) #
                elif lens_type == 'FLEXIONFG':
                    kwargs_lens_sigma_tmp.append({'F1':0.01,'F2':0.01,'G1':0.01,'G2':0.01})
                    kwargs_fixed_lens_tmp.append({'ra_0': kwagrs_lens_list[i][3]['ra_0'],'dec_0': kwagrs_lens_list[i][3]['dec_0'],
                                         'F1':0,'F2':0,'G1':0,'G2':0})
                    kwargs_lower_lens_tmp.append({'F1':-1,'F2':-1,'G1':-1,'G2':-1})
                    kwargs_upper_lens_tmp.append({'F1':1,'F2':1,'G1':1,'G2':1})
        if kwargs_sigma is None:
            kwargs_lens_sigma=kwargs_lens_sigma_tmp
        else:
            kwargs_lens_sigma = kwargs_sigma
        if kwargs_lower is None:
            kwargs_lower_lens = kwargs_lower_lens_tmp
        else:
            kwargs_lower_lens = kwargs_lower
        if kwargs_upper is None:
            kwargs_upper_lens = kwargs_upper_lens_tmp
        else:
            kwargs_upper_lens = kwargs_upper
        if kwargs_fixed is None:
            kwargs_fixed_lens=kwargs_fixed_lens_tmp
        else:
            kwargs_fixed_lens = kwargs_fixed
        self.fixed_index = fixed_index
        self.lens_model_list = lens_model_list
        self.num_img = len(ximg_list)
        self.magnification_list = magnification_list
        lens_params= [kwargs_lens_init, kwargs_lens_sigma, kwargs_fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        return lens_params



    def lens_param_initial(self, x, y,ra, dec, kwargs_data, lens_model_list, diff):
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
            kwargs_lens,magnification  = self._cal_kwargs_lens(x=x, y=y, kwargs_data=kwargs_data, diff=diff, lens_model_list=lens_model_list)
        else:
            kwargs_lens, magnification = self._read_kwargs_lens(ra=ra, dec=dec, kwargs_data=kwargs_data, lens_model_list=lens_model_list)
        return kwargs_lens, magnification



    def _read_kwargs_lens(self,ra,dec, kwargs_data,lens_model_list=['SHIFT','SHEAR','CONVERGENCE','FLEXIONFG']):
        """
        Initial lens parameters, read from shear, convergence
        :param ra: x coordinate of lensed image
        :param dec: y coordinate of lensed image
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
                xgamma, ygamma = self.radec2detector(ra, dec)
                gamma1 = self.gamma1[xgamma, ygamma]
                gamma2 = self.gamma2[xgamma, ygamma]
                kwargs_lens.append({'gamma1': gamma1, 'gamma2': gamma2, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type =='CONVERGENCE':
                xkappa, ykappa = self.radec2detector(ra, dec)
                kappa = self.kappa[ xkappa, ykappa]
                kwargs_lens.append({'kappa_ext': kappa, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'FLEXIONFG':
                kwargs_lens.append(
                    {'F1': 0, 'F2': 0, 'G1': 0, 'G2': 0, 'ra_0': ra_center, 'dec_0': dec_center})
        magnification = 1./ ((1 - kappa)**2 - (gamma1**2 + gamma2**2))
        return kwargs_lens, magnification

    def _cal_kwargs_lens(self,x, y, kwargs_data, lens_model_list, diff):
        """
        #TODO: require same pixel size
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
                gamma1, gamma2 = LensModel(['INTERPOL']).gamma(util.image2array(xaxes),
                                                                 util.image2array(yaxes), kwargs=kwargs_lens_in,diff = diff)
                gamma_1_center, gamma_2_center = gamma1.mean(), gamma2.mean()
                kwargs_lens.append({'gamma1': gamma_1_center, 'gamma2': gamma_2_center, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'CONVERGENCE':
                kappa = LensModel(['INTERPOL']).kappa(util.image2array(xaxes),
                                                        util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                kappa_center = kappa.mean()
                kwargs_lens.append({'kappa_ext': kappa_center, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'FLEXION':
                g1, g2, g3, g4 = LensModel(['INTERPOL']).flexion(util.image2array(xaxes),
                                                                   util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                g1_c, g2_c, g3_c, g4_c = g1.mean(), g2.mean(), g3.mean(), g4.mean()
                kwargs_lens.append({'g1': g1_c,'g2':g2_c,'g3':g3_c,'g4':g4_c, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'FLEXIONFG':
                g1, g2, g3, g4 = LensModel(['INTERPOL']).flexion(util.image2array(xaxes),
                                                                   util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                g1_c, g2_c, g3_c, g4_c = g1.mean(), g2.mean(), g3.mean(), g4.mean()
                F1_c = (g1_c + g3_c) * 0.5
                F2_c = (g2_c + g4_c) * 0.5
                G1_c = (g1_c - g3_c) * 0.5 - g3_c
                G2_c = (g2_c - g4_c) * 0.5 - g4_c
                kwargs_lens.append({'F1': F1_c, 'F2': F2_c, 'G1': G1_c, 'G2': G2_c, 'ra_0': ra_center, 'dec_0': dec_center})

        magnification = LensModel(['INTERPOL']).magnification(util.image2array(xaxes),
                                                                  util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
        magnification = np.mean(magnification)
        return kwargs_lens, magnification




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




    def match_pixelsize(self):

        """
        match the pixel size of data and lens map
        :return:
        """
        return 0




    def radec2detector(self,ra,dec):
        """
        WCS transformation from world to pixel coordinates
        :param ra: ra in deg
        :param dec: dec in deg
        :return:
        """

        y_detector= np.int(self.wcs.wcs_world2pix([[ra, dec]], 1)[0][0])
        x_detector = (np.int(self.wcs.wcs_world2pix([[ra, dec]], 1)[0][1]))

        return x_detector, y_detector