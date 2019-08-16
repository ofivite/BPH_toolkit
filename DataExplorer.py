#
#   TODO
#
# exception handling
# private/protected items
# is is_sum_w2 redundant?

import ROOT
from ROOT import RooFit as RF
from math import sqrt
from pandas import DataFrame
from misc import interactivity_yn
from StatTools import StatTools

class DataExplorer(StatTools):
    """Class exploring data-model relationships.
    """

    def __init__(self, label, data, model):
        super(DataExplorer, self).__init__()
        assert (type(label) is str), 'Label is not str'
        self.data = data
        self.model = model
        self.var = model.getObservables(data).iterator().Next()
        self.label = label
        self.is_fitted = False
        self.fit_status = -999
        self.chi2_test_status = -999

    def set_regions(self, num_of_sigma_window=3, num_of_sigma_to_sdb=2):
        """Set signal region (SR) window and distance to sidebands (SdR)
            SR=|m - mean| < window;
            SdR=|m - mean| > window + distance_to_sdb &
                |m - mean| < 2*window + distance_to_sdb
        NB: Assume that self.model is a double Gaussian

        Parameters
        ---------

        num_of_sigma_window: float, optional (default=3)
            number of effective sigmas in the window
        num_of_sigma_to_sdb: float, optional (default=2)
            number of effective sigmas in between SR and SdR

        Returns
        -------
        self, object
        """
        fr      = self.model.getVariables().find(f'fr_{self.label}').getVal()
        sigma_1 = self.model.getVariables().find(f'sigma_{self.label}_1').getVal()
        sigma_2 = self.model.getVariables().find(f'sigma_{self.label}_2').getVal()
        sigma_eff = sqrt(fr*sigma_1**2 + (1-fr)*sigma_2**2)  ### effective sigma of sum of two gaussians with common
        self.window = num_of_sigma_window*sigma_eff
        self.distance_to_sdb = num_of_sigma_to_sdb*sigma_eff
        return self

    def get_regions(self):
        """Reduce instance dataset with SR and SdR cuts

        Returns
        -------
        data_sig, data_sideband: tuple of RooDataSet
            datasets corresponding to events in SR and SdR
        """
        if self.is_fitted:
            raise Exception('Can\'t get regions: mean should be (if you mean it) MC-fixed value but not fitted to data.')
        mean = self.model.getParameters(self.data).find(f'mean_{self.label}').getVal()
        data_sig = self.data.reduce(f'TMath::Abs({self.var.GetName()} - {mean}) < {self.window}')
        data_sideband = self.data.reduce(f'TMath::Abs({self.var.GetName()} - {mean}) > {self.window + self.distance_to_sdb} && TMath::Abs({self.var.GetName()} -{mean}) < {2.*self.window + self.distance_to_sdb}')
        data_sig.SetName('SR')
        data_sideband.SetName('SdR')
        return data_sig, data_sideband

    def plot_regions(self, frame, y_sdb_left=0, y_sr=0, y_sdb_right=0, line_width=4):
        """Add vertical lines illustrating SR and SdR regions to the frame.
        NB: SR=|m-mean|<window;
            SdR=|m-mean|>window+distance_to_sdb &
                |m-mean|<2*window+distance_to_sdb

        Parameters
        ----------
        frame: RooPlot
            RooPlot frame to draw the lines on
        y_sdb_left: float, optional (default=0)
            y-coordinate for the left sideband region
        y_sr: float, optional (default=0)
            y-coordinate for the signal region
        y_sdb_right: float, optional (default=0)
            y-coordinate for the right sideband region

        Returns
        -------
        frame: RooPlot
        """
        mean = self.model.getParameters(self.data).find(f'mean_{self.label}').getVal()
        line_ll_sdb = (ROOT.TLine(mean - 2.*self.window - self.distance_to_sdb, 0, mean - 2.*self.window - self.distance_to_sdb, y_sdb_left),  ROOT.kBlue-8)
        line_lr_sdb = (ROOT.TLine(mean - self.window - self.distance_to_sdb,    0, mean - self.window - self.distance_to_sdb,    y_sdb_left),  ROOT.kBlue-8)
        line_rl_sdb = (ROOT.TLine(mean + 2.*self.window + self.distance_to_sdb, 0, mean + 2.*self.window + self.distance_to_sdb, y_sdb_right), ROOT.kBlue-8)
        line_rr_sdb = (ROOT.TLine(mean + self.window + self.distance_to_sdb   , 0, mean + self.window + self.distance_to_sdb,    y_sdb_right), ROOT.kBlue-8)
        line_l_sig  = (ROOT.TLine(mean - self.window,                           0, mean - self.window,                           y_sr)        , 47)
        line_r_sig  = (ROOT.TLine(mean + self.window,                           0, mean + self.window,                           y_sr)        , 47)
        lines = (line_ll_sdb, line_lr_sdb, line_rl_sdb, line_rr_sdb, line_l_sig, line_r_sig)
        for line, color in lines:
            line.SetLineColor(color)
            line.SetLineWidth(line_width)
            frame.addObject(line)
        return frame

    def fit(self, is_sum_w2, fix = [], fix_float=[]):
        """Fit instance data with instance model using fitTo() method. Extended or not is infered from the model. Set is_fitted=True.
        NB: the corresponding model parameters will be updated outside of the class instance after executing!

        Parameters
        ----------
        is_sum_w2: bool
            correct Hessian with data weights matrix to get correct errors, see RooFit tutorial rf403__weightedevts
        fix: list of RooRealVar, optional (default=[])
            variables from this list will be fixed in the fit and then released
        fix_float: list of RooRealVar, optional (default=[])
            variables from this list will be firstly setConstant(1) in the fit and then setConstant(0)

        Returns
        -------
        fit_results: RooFitResult
        """
        is_extended = self.model.canBeExtended()
        for param in fix:
            param.setConstant(1)
        self.model.fitTo(self.data, RF.Extended(is_extended), RF.SumW2Error(is_sum_w2))
        for param in fix_float:
            param.setConstant(1)
        self.model.fitTo(self.data, RF.Extended(is_extended), RF.SumW2Error(is_sum_w2))
        for param in fix_float:
            param.setConstant(0)
        self.model.fitTo(self.data, RF.Extended(is_extended), RF.SumW2Error(is_sum_w2))
        fit_results = self.model.fitTo(self.data, RF.Extended(is_extended), RF.SumW2Error(is_sum_w2), RF.Save())
        fit_results.Print()
        for param in fix:
            param.setConstant(0)
        self.is_fitted = True
        self.fit_status = fit_results.status()
        if is_sum_w2:
            print('\n\n' + 70*'~' + '\n' + ' '*30 + 'BEWARE!\n\nErrors might differ between two printed tables!\nThe last one from RooFitResult.Print() should be correct.\nYou might also want to consider chi2_fit() method as a cross-check,\nas in principle, that should give correct and more reliable results\n(but the normalization in this case will likely not be preserved \nand results might be unstable)\n' + 70*'~' + '\n\n')
        return fit_results

    def chi2_fit(self, fix_float=[], minos = False, poi = None):
        """Fit the instance data with binned chi2 method using Minuit2. Set is_fitted=True
        NB: weights presence is taken care of automatically

        Parameters
        ----------

        fix_float: list of RooRealVar, optional (default=[])
        variables from this list will be firstly setConstant(1) in the fit and then setConstant(0)

        minos: bool
            whether to calculate MINOS errors for POI

        poi: RooRealVar
            parameter of interest for which to calculate MINOS errors

        Returns
        -------
        self, object
        """
        hist_to_fit = ROOT.RooDataHist('hist_to_fit', 'hist_to_fit', ROOT.RooArgSet(self.var), self.data) ### binning is taken from the var's definition
        is_extended = self.model.canBeExtended()
        chi = ROOT.RooChi2Var("chi","chi", self.model, hist_to_fit, RF.Extended(is_extended), RF.DataError(ROOT.RooAbsData.Auto))
        m = ROOT.RooMinimizer(chi)
        m.setMinimizerType("Minuit2")
        m.setPrintLevel(3)
        m.minimize("Minuit2","minimize")
        for param in fix_float:
            param.setConstant(1)
        m.minimize("Minuit2","minimize")
        for param in fix_float:
            param.setConstant(0)
        self.fit_status = m.minimize("Minuit2","minimize")
        self.is_fitted = True
        if minos:
            if poi is None:
                raise TypeError('Poi is None by default: set it to a proper variable to run MINOS.')
            if self.data.isWeighted():
                interactivity_yn('The data is weighted and MINOS should not be used. Sure you want to proceed?')
            m.minos(ROOT.RooArgSet(poi))
        return m.save()

    def plot_on_frame(self, title=' ', plot_params=ROOT.RooArgSet(), **kwargs):
        """Plot the instance model with all its components and data on the RooPlot frame
        NB: signal component's name should starts with 'sig_', background - with 'bkgr_'
        Parameters
        ----------
        title: str, optional (default=' ')
            title for a RooPlot frame
        plot_params: RooArgSet, optional (default=RooArgSet)
            Set of parameters to be shown on the legend

        Returns
        -------
        frame: RooPlot
        """
        var_left  = self.var.getMin();
        var_right = self.var.getMax();
        var_nbins = self.var.numBins()

        frame = ROOT.RooPlot(" ", title, self.var, var_left, var_right, var_nbins)  # frame.getAttText().SetTextSize(0.053)
        self.data.plotOn(frame, RF.DataError(ROOT.RooAbsData.Auto))
        self.model.plotOn(frame, RF.LineColor(ROOT.kRed-6), RF.LineWidth(5)) #, RF.NormRange('full'), RF.Range('full')
        if plot_params.getSize() == 0:
            self.model.paramOn(frame, RF.Layout(0.55, 0.96, 0.9))
        else:
            self.model.paramOn(frame, RF.Layout(0.55, 0.96, 0.9), RF.Parameters(plot_params))

        iter = self.model.getComponents().iterator()
        iter_comp = iter.Next()
        while iter_comp:
            if iter_comp.GetName().startswith('sig_'):
                self.model.plotOn(frame, RF.Components(iter_comp.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4))
            if iter_comp.GetName().startswith('bkgr_'):
                self.model.plotOn(frame, RF.Components(iter_comp.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kBlue-8), RF.LineWidth(4))
            iter_comp = iter.Next()

        for component_name, range in kwargs.items():
            assert (type(component_name) is str), 'component type is not str'
            assert (type(range) is list), 'range is not list'
            assert (len(range) == 2), 'N elements in range != 2'
            self.model.plotOn(frame, RF.Components(component_name), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kGreen-5), RF.LineWidth(4), RF.Normalization(1.0), RF.Name(component_name), RF.Range(*range))

        frame.GetYaxis().SetTitle(f'Candidates / {round((var_right - var_left) * 1000. / var_nbins, 1)} MeV')
        frame.GetXaxis().SetTitleSize(0.04)
        frame.GetYaxis().SetTitleSize(0.04)
        frame.GetXaxis().SetLabelSize(0.033)
        frame.GetYaxis().SetLabelSize(0.033)
        frame.GetXaxis().SetTitleOffset(1.05)
        frame.GetYaxis().SetTitleOffset(1.3)
        return frame

    def write_to_workspace(self, poi, nuisances):
        """Create a workspace with the fitted to the data model, poi and nuisance parameters.

        Parameters
        ----------

        nuisances: list of RooRealVar
            nuisance parameters in statistical inference

        poi: RooRealVar
            parameter of interest in statistical inference

        Returns
        -------
        w: RooWorkspace
        """
        if not self.is_fitted:
            interactivity_yn('Model was not fitted to data, fit it first. Sure you want to proceed?')
        w = ROOT.RooWorkspace("w", True)
        Import = getattr(ROOT.RooWorkspace, 'import') # special trick to make things not crush
        Import(w, self.model)
        mc = ROOT.RooStats.ModelConfig("ModelConfig", w)
        mc.SetPdf(w.pdf(self.model.GetName()))
        mc.SetParametersOfInterest(ROOT.RooArgSet(w.var(poi.GetName())))
        mc.SetObservables(ROOT.RooArgSet(w.var(self.var.GetName())))
        mc.SetNuisanceParameters(ROOT.RooArgSet(*[w.var(nui.GetName()) for nui in nuisances]))
        mc.SetSnapshot(ROOT.RooArgSet(w.var(poi.GetName())))
        Import(w, mc, 'ModelConfig')
        Import(w, self.data, 'data')
        return w

    @staticmethod
    def extract_from_workspace(w):
        """Extract data, signal+background and background-only models from a given RooWorkspace.
        Background-only model is taken from the s+b model by setting the parameter of interest to be 0.

        NB: naming conventions assume that data's name is 'data', and that the parameter of interest is the first one and corresponds to the signal yield.

        Parameters:
        -----------

        w, RooWorkspace
            workspace with saved data and s+b model to unpack

        Returns:
        --------
        data, mc_sb, mc_b: RooAbsData, ModelConfig, ModelConfig
            Tuple with data, s+b and b-only models
        """
        data = w.obj("data")
        mc_sb = w.obj("ModelConfig")
        mc_sb.LoadSnapshot() # not sure whether it is useful
        mc_sb.SetName("mc_sb")
        poi = mc_sb.GetParametersOfInterest().first()

        mc_b = mc_sb.Clone()
        mc_b.SetName("B_only_model")
        oldval = poi.getVal()
        poi.setVal(0)
        mc_b.SetSnapshot(ROOT.RooArgSet(poi))
        poi.setVal(oldval)

        return data, mc_sb, mc_b

    @staticmethod
    def fix_shapes(workspaces_dict, models_dict, var_ignore_list):
        """Recursively fix model's parameters according to the values for the model in the corresponding workspace.
        Workspace with 'fix from' model must share in the dictionary the same 'binding' key with 'to fix' model.

        Parameters
        ----------

        workspaces_dict: dictionary
            dictionary with binding labels and workspaces which carry 'fix from' models

        models_dict: dictionary
            dictionary with binding labels and 'to fix' models

        var_ignore_list: list of RooRealVar
            parameters of models (e.g. means) which will not be setConstant (but the values will be set to a workspace one)

        Returns
        -------
        None
        """
        if len(workspaces_dict) < len(models_dict):
            raise Exception('There is more models than corresponding workspaces.')
        for key, s in models_dict.items():
            iter = s.getVariables().iterator()
            iter_comp = iter.Next()
            while iter_comp:
                if key not in workspaces_dict.keys():
                    raise Exception(f'Can\'t find the \'{key}\' key in the workspaces dictionary.')
                try:
                    val = workspaces_dict[key].var(iter_comp.GetName()).getVal()
                except:
                    if iter_comp in var_ignore_list:
                        print(f'Can\'t find {iter_comp.GetName()} variable in the workspace: skipping it as ignore_listed.')
                        pass
                    else:
                        raise Exception(f'Can\'t find {iter_comp.GetName()} variable in the workspace: check that models match.')
                else:
                    iter_comp.setVal(val)
                    if iter_comp.GetName() not in [v.GetName() for v in var_ignore_list]:
                        iter_comp.setConstant(1)
                finally:
                    iter_comp = iter.Next()

################################################################################################################################

print('\n\n         ~~~' + '\n\n')
