from __future__ import print_function
from lenstronomy.Workflow.fitting_sequence import FittingSequence
import numpy as np


class ClsrWorkflow(object):
    def __init__(self, kwargs_data_joint, kwargs_model,lens_params,source_params, kwargs_constraints, kwargs_likelihood=None):
        """
        class to  manage cluster source reconstruction.
        This class inherited the FittingSequence class in Workflow module of lenstronomy.
        :param kwargs_data_joint: keywords arguments of [data, psf, numericals] in lenstronomy convention.
        :param kwargs_model: name of model list
        :param lens_params: lens model keywords arguments [kwargs_lens_init, kwargs_lens_sigma, kwargs_fixed_lens, kwargs_lower_lens, kwargs_upper_lens]
        :param source_params: source model keywords arguments [kwargs_source_init, kwargs_source_sigma, kwargs_fixed_source, kwargs_lower_source, kwargs_upper_source]
        :param kwargs_constraints: contraints on models
        :param kwargs_likelihood: options of calculating likelihood, see more: LikelihoodModule class in Sampling module of lenstronomy.
        """
        self.kwargs_data_joint =kwargs_data_joint
        self.multi_band_list = kwargs_data_joint.get('multi_band_list', [])
        self.kwargs_model =kwargs_model
        kwargs_params = {'lens_model': lens_params, 'source_model': source_params}
        self.kwargs_params= kwargs_params
        if kwargs_constraints is None:
            kwargs_constraints ={}
        if kwargs_likelihood is None:
            kwargs_likelihood = {'source_marg': False, 'check_positive_flux': True}
        self.fitting_seq_src = FittingSequence(kwargs_data_joint, kwargs_model, kwargs_constraints, kwargs_likelihood, kwargs_params)

    def run_fit_sequence(self,fitting_kwargs_list):
        """
        :param fitting_kwargs_list: list of [['string', {kwargs}], ..] with 'string being the specific fitting option and kwargs being the arguments passed to this option
        :return: fitting results
        """
        chain_list = self.fitting_seq_src.fit_sequence(fitting_kwargs_list)
        kwargs_result = self.fitting_seq_src.best_fit(bijective=False)
        bic_model = self.fitting_seq_src.bic
        return bic_model,chain_list, kwargs_result


    def lensmodel_comp(self, num_img, fixed_index, n_particles,n_iterations,sigma_scale, num_lens_model,flexion_option=True):
        """
        function to figure out the necessary of increasing lens model complexity. Currently, we only consider up to flexion term.
        :param n_particles: particles numbers of PSO
        :param n_iterations: iteration numbers of PSO
        :param sigma_scale: sigma scale for PSO
        :param num_img: int, numbers of lensed image
        :param fixed_index: int, index of fixed lensed image
        :param num_lens_model: numbers of strings contained in lens model list
        :param flexion_option: bool, default is taking flexion into consideration
        :return: pso results, fitting results of necessary lens model complexity, bic values
        """
        lens_remove_fixed_list = []
        lens_add_fixed_list = []
        for i in range(num_img):
            if i == fixed_index:
                print ("lens model keep fixed in frame:", i+1)
            else:
                lens_flexion_index = (i + 1) * num_lens_model - 1
                lens_remove_fixed_list.append([lens_flexion_index, ['G1', 'G2', 'F1', 'F2'], [0, 0, 0, 0]])
                lens_add_fixed_list.append([lens_flexion_index, ['G1', 'G2', 'F1', 'F2'], [0, 0, 0, 0]])
        flexion_add_fixed = [['update_settings', {'lens_add_fixed': lens_add_fixed_list}]]
        kwargs_pso = [['PSO', {'sigma_scale': sigma_scale, 'n_particles': n_particles, 'n_iterations': n_iterations}]]
        fitting_kwargs_fix =  flexion_add_fixed+ kwargs_pso
        bic_model_fix, chain_list_fix, kwargs_result_fix = self.run_fit_sequence(fitting_kwargs_fix)
        if flexion_option:
            flexion_remove_fixed = [['update_settings', {'lens_remove_fixed': lens_remove_fixed_list}]]
            fitting_kwargs_free = flexion_remove_fixed + kwargs_pso
            bic_model_free, chain_list_free, kwargs_result_free = self.run_fit_sequence(fitting_kwargs_free)
        else:
            bic_model_free = 10000000
        if bic_model_free > bic_model_fix:
            print ("No necessary to add flexion!")
            bic_list = [bic_model_fix]
            chain_list = [chain_list_fix]
            kwargs_result_list = [kwargs_result_fix]
            self._update_kwargs(kwargs_result_fix)
            _, _, _ = self.run_fit_sequence(flexion_add_fixed)
        elif bic_model_free < bic_model_fix:
            print("Flexion is needed!")
            bic_list = [bic_model_fix, bic_model_free]
            chain_list = [chain_list_fix, chain_list_free]
            kwargs_result_list = [kwargs_result_fix,kwargs_result_free]
        return chain_list, kwargs_result_list, bic_list


    def sourcemodel_comp(self, n_max_range=[0], sr = 0, n_particles=10, n_iterations=10,sigma_scale =1.0,
                         bic_model_in = [100000], chain_list_in = [0], kwargs_results_in = [0], bic_option=True) :
        """
        function to found the best fitting results among source models related to numbers of shapelets basis
        :param n_max_range: shapelets basis selection range
        :param sr: typical scale (") in source plane
        :param n_particles: particles numbers of PSO
        :param n_iterations: iteration numbers of PSO
        :param sigma_scale: sigma scale for PSO
        :param bic_model_in: BIC value of models before adding shapelets to source model
        :param chain_list_in: PSO chain results of models before adding shapelets to source model
        :param kwargs_results_in: fitting results of models before adding shapelets to source model
        :return: best-fit PSO results, best-fits results, PSO results for n_max_range, fitting results for n_max_range,
                bic values for all models
        """
        bic_model_list = bic_model_in
        chain_list_list = chain_list_in
        kwargs_result_list = kwargs_results_in
        bic_in_len = len(bic_model_in)
        bic_run = True
        beta0 = sr
        kwargs_pso = [['PSO', {'sigma_scale': sigma_scale, 'n_particles': n_particles, 'n_iterations': n_iterations}]]
        for nmax in n_max_range:
            if nmax < 0:
                raise ValueError("nmax can not be negative!",nmax)
            else:
                if nmax == n_max_range[0]:
                    start_kwargs_shapelet = [['update_settings', {'source_remove_fixed': [  [1, ['beta'], [beta0]] ]}]]
                else:
                    start_kwargs_shapelet = []
                beta_nmax = ((nmax + 1)) ** 0.5 * beta0
                fit_kwargs_shapelet = [['update_settings',
                                        {'source_add_fixed': [[1, ['n_max'], [nmax]]],
                                        'change_source_lower_limit': [[1, ['beta'], [beta_nmax]]]
                                         }
                                         ]]
                fitting_kwargs = start_kwargs_shapelet + fit_kwargs_shapelet + kwargs_pso
            if bic_run:
                print ("nmax",nmax,"fitting_kwargs",fitting_kwargs)
                bic_model,chain_list, kwargs_result = self.run_fit_sequence(fitting_kwargs)
                if bic_model >  bic_model_list[-1]:
                    if bic_option:
                        bic_run = False
                        if bic_model > bic_model_in[bic_in_len-1]:
                            print ("bic_model_in",bic_model_in)
                            print ("no necessary to add SHAPELETS !")
                            fix_kwargs_shapelet=[['update_settings', {'source_add_fixed': [[1, ['beta'], [sr]]]}]]
                            _, _, _ = self.run_fit_sequence(fix_kwargs_shapelet)
                    elif not bic_option:
                        chain_list_list.append(chain_list)
                        kwargs_result_list.append(kwargs_result)
                        bic_model_list.append(bic_model)
                    print ("no necessary to increase model complexity!")
                elif bic_model < bic_model_list[-1]:
                    chain_list_list.append(chain_list)
                    kwargs_result_list.append(kwargs_result)
                    bic_model_list.append(bic_model)
                    print (bic_model, "currently is the lowest BIC value in bic_model_list=", bic_model_list)
        bic_sourcemodel = bic_model_list[bic_in_len:]
        if bic_sourcemodel ==[]:
            chain_list_lowest = chain_list_in[-1]
            kwargs_result_lowest = kwargs_results_in[-1]
        else:
            index_bic_minima = np.where(bic_model_list == np.min(bic_model_list))[0][0]
            chain_list_lowest = chain_list_list[index_bic_minima]
            kwargs_result_lowest = kwargs_result_list[index_bic_minima]
        return  chain_list_lowest, kwargs_result_lowest, chain_list_list, kwargs_result_list, bic_model_list


    def _update_kwargs(self,kwargs_result):
        """
        :param kwargs_result: fitting results of a specific state
        :return: go back to a specific state
        """
        self.fitting_seq_src.update_state(kwargs_result)

