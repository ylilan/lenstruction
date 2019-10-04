from __future__ import print_function

class SourcePreparation(object):
    """
     This class is used to do initialize kwargs of the source light models,
     the default light models are sersic and shapelets, ['SERSIC_ELLIPSE','SHAPELETS'].
    """
    def __init__(self, source_model_list=['SERSIC_ELLIPSE','SHAPELETS']):
        """
        :param source_model_list: list of strings, name of lightmodel profiles.
        """
        self.source_model_list = source_model_list

    def params(self, deltaPix, source_params = None,betax=0, betay=0):
        """
        source parameters configuration in lenstronomy keywords arguments.
        :param betax: centerx in source plane
        :param betay:  centery in source plane
        :param deltaPix: pixel scale in image plane
        :param source_params: list of source parameters configuration,
        [kwargs_source_init, kwargs_source_sigma, kwargs_fixed_source, kwargs_lower_source,
                         kwargs_upper_source]
        :return:
        """
        kwargs_source_init = []
        kwargs_source_sigma = []
        kwargs_lower_source = []
        kwargs_upper_source = []
        kwargs_fixed_source =[]
        for source_type in self.source_model_list:
            if source_type == 'SERSIC_ELLIPSE':
                kwargs_source_init.append({'R_sersic': deltaPix, 'center_x': betax, 'center_y': betay,
                                           'e1': 0., 'e2': 0., 'n_sersic': 1})
                kwargs_source_sigma.append({'R_sersic': deltaPix, 'n_sersic': 0.1, 'e1': 0.1,
                                            'e2': 0.1, 'center_x': 0.03, 'center_y': 0.03})
                kwargs_lower_source.append({'R_sersic': 0, 'n_sersic': 0.1,
                                            'e1': -1, 'e2': -1, 'center_x': betax - 1,'center_y': betay - 1})
                kwargs_upper_source.append({'R_sersic': 100, 'n_sersic': 8, 'e1': 1, 'e2': 1,
                                            'center_x': betax + 1, 'center_y': betay + 1})
                kwargs_fixed_source.append({})

            elif source_type == 'SHAPELETS':
                kwargs_source_init.append({'n_max': -1, 'beta': 0.1, 'center_x':betax, 'center_y': betay})
                kwargs_source_sigma.append( {'beta': deltaPix,
                                             'center_x': deltaPix, 'center_y': deltaPix})
                kwargs_lower_source.append({'beta': 0,
                                            'center_x': betax - 1, 'center_y': betay - 1})
                kwargs_upper_source.append({'beta': 100, 'center_x':betax + 1,'center_y': betay + 1})
                kwargs_fixed_source.append({'n_max': -1, 'beta':deltaPix })
        if source_params is None:
            sourceparams = [kwargs_source_init, kwargs_source_sigma, kwargs_fixed_source, kwargs_lower_source,
                         kwargs_upper_source]
        else:
            sourceparams = source_params
        return sourceparams

    def constrain(self, source_constrain=None):
        """
        :param source_constrain: dict, fitting constrains of source light models
        :return: fitting constrains of source light models
        """
        if source_constrain is None:
            source_constrain ={'joint_source_with_source': [[0, 1, ['center_x', 'center_y']]]}
        return source_constrain

    def model_list(self):
        """
        :return: source model list
        """
        source_model_list ={'source_light_model_list': self.source_model_list}
        return source_model_list