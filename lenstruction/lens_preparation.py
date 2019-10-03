import numpy as np
import lenstronomy.Util.util as util
from lenstronomy.LensModel.numeric_lens_differentials import NumericLens
from lenstronomy.Data.imaging_data import ImageData
from lenstronomy.Cosmo.lens_cosmo import LensCosmo
from astropy.io import fits


class LensPreparation(object):
    """
    This class is used to do initialize kwargs of the lens models according to the given lens_model_list.
    (It's easy and free to choose lens_model_list, but a little tricky to generate corresponding kwargs.)
    For the situation that you know deflection maps, then you want know shear, convergency, and shift ....
    Or, you know, shear, convergency, shift, then you wanna know deflection maps.

    """
    def __init__(self,xdeflection,ydeflection,zl,zs,lenseq='1'):
       """
       :param alphax: the x deflection map corresponds to data in dataclass
       :param alphay: the y deflection map corresponds to data in dataclass
       :param zl: redshift of the lens
       :param zs:  redshift of the lens
       :param lenseq:
       """
       alphax_map_hat = -fits.open(xdeflection)[0].data
       alphay_map_hat = -fits.open(ydeflection)[0].data
       if lenseq == '1':
           alphax_map_hat = alphax_map_hat
           alphay_map_hat = alphay_map_hat
       elif lenseq == '-1':
           alphax_map_hat = -alphax_map_hat
           alphay_map_hat = -alphay_map_hat
       cosmo = LensCosmo(zl, zs)
       dlsds = cosmo.D_ds / cosmo.D_s
       self.alphax = alphax_map_hat* dlsds
       self.alphay = alphay_map_hat* dlsds


    def initial_kwargs_lens(self,x,y,kwargs_data,alphax_shift=0, alphay_shift=0, diff= 0.01, lens_model_list=['SHIFT','SHEAR','CONVERGENCE','FLEXIONFG']):
        """
        This function returns list type of kwargs of lens models.
        :param lens_model_list: list of strings with lens model names
        :param alphax_shift: shift the source's x position
        :param alphay_shift: shift the source's y position
        :return:  a list of kwargs of lens models corresponding to the lens models existed in lens_model_list
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
                alpha_x_center = alphax[int(r_cut + 1), int(r_cut + 1)]
                alpha_y_center = alphay[int(r_cut + 1), int(r_cut + 1)]
                kwargs_lens.append({'alpha_x': alpha_x_center - alphax_shift, 'alpha_y': alpha_y_center - alphay_shift})
            elif lens_type == 'SHEAR':
                gamma1, gamma2 = NumericLens(['INTERPOL']).gamma(util.image2array(xaxes),
                                                                 util.image2array(yaxes), kwargs=kwargs_lens_in,diff = diff)
                gamma_1_center, gamma_2_center = gamma1.mean(), gamma2.mean()
                kwargs_lens.append({'e1': gamma_1_center, 'e2': gamma_2_center, 'ra_0': ra_center, 'dec_0': dec_center})
                self.gamma1 = gamma1
                self.gamma2 = gamma2
            elif lens_type == 'CONVERGENCE':
                kappa = NumericLens(['INTERPOL']).kappa(util.image2array(xaxes),
                                                        util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                kappa_center = kappa.mean()
                kwargs_lens.append({'kappa_ext': kappa_center, 'ra_0': ra_center, 'dec_0': dec_center})
                self.kappa = kappa
            elif lens_type == 'FLEXION':
                g1, g2, g3, g4 = NumericLens(['INTERPOL']).flexion(util.image2array(xaxes),
                                                                   util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                self.g1 = g1
                self.g2 = g2
                self.g3 = g3
                self.g4 = g4
                g1_c, g2_c, g3_c, g4_c = g1.mean(), g2.mean(), g3.mean(), g4.mean()
                kwargs_lens.append({'g1': g1_c,'g2':g2_c,'g3':g3_c,'g4':g4_c, 'ra_0': ra_center, 'dec_0': dec_center})
            elif lens_type == 'FLEXIONFG':
                g1, g2, g3, g4 = NumericLens(['INTERPOL']).flexion(util.image2array(xaxes),
                                                                   util.image2array(yaxes), kwargs=kwargs_lens_in, diff = diff)
                self.F1 = (g1 + g3) * 0.5
                self.F2 = (g2 + g4) * 0.5
                self.G1 = (g1 - g3) * 0.5 - g3
                self.G2 = (g2 - g4) * 0.5 - g4
                g1_c, g2_c, g3_c, g4_c = g1.mean(), g2.mean(), g3.mean(), g4.mean()
                F1_c = (g1_c + g3_c) * 0.5
                F2_c = (g2_c + g4_c) * 0.5
                G1_c = (g1_c - g3_c) * 0.5 - g3_c
                G2_c = (g2_c - g4_c) * 0.5 - g4_c
                kwargs_lens.append({'F1': F1_c, 'F2': F2_c, 'G1': G1_c, 'G2': G2_c, 'ra_0': ra_center, 'dec_0': dec_center})
        magnification=NumericLens(['INTERPOL']).magnification(util.image2array(xaxes),
                                                                 util.image2array(yaxes), kwargs=kwargs_lens_in ,diff = diff)
        magnification=np.abs(magnification.mean())
        return kwargs_lens, magnification






    def kwargs_lens_configuration(self,ximg_list, yimg_list, kwargs_data_joint, fixed_index=None,
                                  lens_model_list_in= ['SHIFT','SHEAR','CONVERGENCE','FLEXIONFG'], diff = 0.01):
        self.lens_model_list_in = lens_model_list_in
        kwagrs_lens_list = []
        kwargs_lens_init_list = []
        betax_list = []  # x_center in the source plane
        betay_list = []  # y_center in the source plane
        mag_map_list = []  # magnification
        lens_model_list = []
        kwargs_lens_sigma_tmp = []
        kwargs_fixed_lens_tmp = []
        kwargs_lower_lens_tmp = []
        kwargs_upper_lens_tmp = []
        index_lens_model_list = []
        for i in range(len(ximg_list)):
             kwargs_lens, magnification = self.initial_kwargs_lens(
                        ximg_list[i], yimg_list[i], kwargs_data_joint['multi_band_list'][i][0], diff=diff)
             kwagrs_lens_list.append(kwargs_lens)
             thetax, thetay = kwargs_lens[1]['ra_0'], kwargs_lens[1]['dec_0']
             betax, betay = thetax - kwargs_lens[0]['alpha_x'], thetay - kwargs_lens[0]['alpha_y']
             betax_list.append(betax)
             betay_list.append(betay)
             mag_map_list.append(magnification)
        self.magnification =  mag_map_list
        if fixed_index is not None:
            img_index = fixed_index
        else:
            if len(ximg_list) > 1:
                img_index = np.where(mag_map_list == np.min(mag_map_list))[0][0]
            else:
                img_index = 0
        self.img_index = img_index
        dbetax = betax_list - betax_list[img_index]
        dbetay = betay_list - betay_list[img_index]
        for i in range(len(ximg_list)):
            kwagrs_lens_list[i][0]['alpha_x'] = kwagrs_lens_list[i][0]['alpha_x'] + dbetax[i]
            kwagrs_lens_list[i][0]['alpha_y'] = kwagrs_lens_list[i][0]['alpha_y'] + dbetay[i]
            kwargs_lens_init_list += (kwagrs_lens_list[i])
        deltapix = np.max(kwargs_data_joint['multi_band_list'][0][0]['transform_pix2angle'][0])
        d_p = 3 * deltapix
        for i in range(len(ximg_list)):
            for lens_type in lens_model_list_in:
                if lens_type == 'SHIFT':
                    kwargs_lens_sigma_tmp.append({'alpha_x':d_p,'alpha_y':d_p})
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
            index_lens_model_list.append([i*4,i*4+1,i*4+2,i*4+3])
            lens_model_list+=lens_model_list_in
        kwargs_lens_init=kwargs_lens_init_list
        kwargs_lens_sigma=kwargs_lens_sigma_tmp
        kwargs_lower_lens=kwargs_lower_lens_tmp
        kwargs_upper_lens=kwargs_upper_lens_tmp
        kwargs_fixed_lens=kwargs_fixed_lens_tmp
        self.num_img = len(ximg_list)
        self.betax = betax_list[img_index]
        self.betay = betay_list[img_index]
        self.lens_model_list = lens_model_list
        self.index_lens_model_list = index_lens_model_list
        self.min_mag = np.abs(np.min(mag_map_list))
        lens_params= [kwargs_lens_init, kwargs_lens_sigma, kwargs_fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        return lens_params


    def constrain(self,constrain={}):
        lens_constrian = constrain
        return lens_constrian

    def model_index_list(self):
        lens_model_index_list = {'lens_model_list': self.lens_model_list,'index_lens_model_list': self.index_lens_model_list}
        return lens_model_index_list


    def kwargs_flexion(self,bic_noflexion=1, bic_flexion=0):
        num_lens_model_list_in = len(self.lens_model_list_in)
        lens_remove_fixed_list = []
        lens_add_fixed_list = []
        for i in range(self.num_img):
            if i == self.img_index:
                print "lens model keep fixed in frame:", i + 1
            else:
                lens_remove_index = (i + 1) * num_lens_model_list_in - 1
                lens_remove_fixed_list.append([lens_remove_index, ['G1', 'G2', 'F1', 'F2'], [0, 0, 0, 0]])
                lens_add_fixed_list.append([lens_remove_index, ['G1', 'G2', 'F1', 'F2'], [0, 0, 0, 0]])
        flexion_remove_fixed = ['update_settings', {'lens_remove_fixed': lens_remove_fixed_list}]
        flexion_add_fixed = ['update_settings', {'lens_add_fixed': lens_add_fixed_list}]
        if bic_flexion > bic_noflexion:
            kwargs_fitting_flexion = flexion_add_fixed
            print ("no need to add flexion")
        else:
            kwargs_fitting_flexion = flexion_remove_fixed
        return kwargs_fitting_flexion



