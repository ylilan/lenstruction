from __future__ import print_function

class LenslightSpecify(object):
    """
     This class is used to do initialize kwargs of the source light models,
     the default light models are sersic and shapelets, ['SERSIC_ELLIPSE','SHAPELETS'].
    """
    def __init__(self, xcenter=[], ycenter=[]):
        """
        :param source_model_list: list of strings, name of lightmodel profiles.
        """
        self.x = xcenter
        self.y = ycenter

        if len(xcenter) > 0:
            self.source_model_list = ['SERSIC_ELLIPSE'] * len(xcenter)
        else:
            self.source_model_list = []


    def params(self, re=0.1,
               kwargs_init = None, kwargs_sigma = None, kwargs_fixed = None, kwargs_lower = None,
               kwargs_upper = None):
        """
        source parameters configuration in lenstronomy keywords arguments.
        :param re: float (arcsec unit), typical length scale in source plane
        :param betax: centerx in source plane
        :param betay: centery in source plane
        :param kwargs_init: list, initial keywords arguments
        :param kwargs_sigma: list, sigma keywords arguments
        :param kwargs_fixed: list, fixed keywords arguments
        :param kwargs_lower: list, lower boundary keywords arguments
        :param kwargs_upper: list, upper boundary keywords arguments
        :return: list of source parameters keywords arguments
        """
        kwargs_init_tmp = []
        kwargs_sigma_tmp = []
        kwargs_lower_tmp = []
        kwargs_upper_tmp = []
        kwargs_fixed_tmp =[]
        xcenter = self.x
        ycenter = self.y
        for i in range(len(xcenter)):
            kwargs_init_tmp.append({'R_sersic': re, 'center_x': xcenter[i], 'center_y': ycenter[i],
                                           'e1': 0., 'e2': 0., 'n_sersic': 4})
            kwargs_sigma_tmp.append({'R_sersic': re, 'n_sersic': 1, 'e1': 0.1, 'e2': 0.1, 'center_x': re, 'center_y': re})
            kwargs_lower_tmp.append({'R_sersic': 0, 'n_sersic': 0.5, 'e1': -1, 'e2': -1, 'center_x': xcenter[i] - 1,'center_y':  ycenter[i] - 1})
            kwargs_upper_tmp.append({'R_sersic': 100, 'n_sersic': 8, 'e1': 1, 'e2': 1,   'center_x': xcenter[i] + 1, 'center_y': ycenter[i] + 1})
            kwargs_fixed_tmp.append({'n_sersic': 4.4})

        if  kwargs_init is None:
            kwargs_source_init = kwargs_init_tmp
        else:
            kwargs_source_init = kwargs_init

        if  kwargs_sigma is None:
            kwargs_source_sigma = kwargs_sigma_tmp
        else:
            kwargs_source_sigma = kwargs_sigma

        if  kwargs_fixed is None:
            kwargs_fixed_source = kwargs_fixed_tmp
        else:
            kwargs_fixed_source = kwargs_fixed

        if  kwargs_lower is None:
            kwargs_lower_source = kwargs_lower_tmp
        else:
            kwargs_lower_source = kwargs_lower

        if kwargs_upper is None:
            kwargs_upper_source = kwargs_upper_tmp
        else:
            kwargs_upper_source = kwargs_upper
        source_params = [kwargs_source_init, kwargs_source_sigma, kwargs_fixed_source, kwargs_lower_source, kwargs_upper_source]
        return source_params



    def constrain(self, source_constrain=None):
        """
        :param source_constrain: dict, fitting constrains of source light models
        :return: fitting constrains of source light models
        """
        if source_constrain is None:
            source_constrain ={}
        return source_constrain

    def model_list(self):
        """
        :return: source model list
        """
        lens_light_model_list ={'lens_light_model_list': self.source_model_list}
        return lens_light_model_list