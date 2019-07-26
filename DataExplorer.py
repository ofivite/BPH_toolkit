import ROOT
from ROOT import RooFit as RF
from scipy.stats import chi2
from math import sqrt
from pandas import DataFrame

class DataExplorer(object):
    """Base class exploring data-model relationship in Bs->X(3872)phi study"""

    def __init__(self, label, data, model):
        super(DataExplorer, self).__init__()
        assert (type(label) is str), 'Label is not str'
        self.data = data
        self.model = model
        self.var = model.getObservables(data).iterator().Next()
        self.label = label
        self.is_fitted = False

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
        sigma_eff = sqrt(fr*sigma_1**2 + (1-fr)*sigma_2**2)  ### effective sigma of sum of two gaussians with common mean
        #
        self.window = num_of_sigma_window*sigma_eff #
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
            return
        mean = self.model.getParameters(self.data).find(f'mean_{self.label}').getVal()
        data_sig = self.data.reduce(f'TMath::Abs({self.var.GetName()} - {mean}) < {self.window}')
        data_sideband = self.data.reduce(f'TMath::Abs({self.var.GetName()} - {mean}) > {self.window + self.distance_to_sdb} && TMath::Abs({self.var.GetName()} -{mean}) < {2.*self.window + self.distance_to_sdb}')
        data_sig.SetName('sig')
        data_sideband.SetName('sideband')
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
        #
        for line, color in lines:
            line.SetLineColor(color)
            line.SetLineWidth(line_width)
            frame.addObject(line)
        return frame

    def fit(self, is_sum_w2, fix_float=[]):
        """Fit instance data with instance model using fitTo() method. Extended or not is infered from the model. Set is_fitted=True.
        NB: the corresponding model parameters will be updated outside of the class instance after executing!

        Parameters
        ----------
        is_sum_w2: bool
            correct Hessian with data weights matrix to get correct errors, see RooFit tutorial rf403__weightedevts
        fix_float: list of RooRealVar, optional (default=[])
            variables from this list will be firstly setConstant(1) in the fit and then setConstant(0)

        Returns
        -------
        fit_results: RooFitResult
        """
        is_extended = self.model.canBeExtended()
        self.model.fitTo(self.data, RF.Extended(is_extended), RF.SumW2Error(is_sum_w2))
        for param in fix_float:
            param.setConstant(1)
        self.model.fitTo(self.data, RF.Extended(is_extended), RF.SumW2Error(is_sum_w2))
        for param in fix_float:
            param.setConstant(0)
        self.model.fitTo(self.data, RF.Extended(is_extended), RF.SumW2Error(is_sum_w2))
        fit_results = self.model.fitTo(self.data, RF.Extended(is_extended), RF.SumW2Error(is_sum_w2), RF.Save())
        fit_results.Print()
        self.is_fitted = True
        if is_sum_w2:
            print('\n\n' + 65*'~' + '\n' + ' '*30 + 'BEWARE!\nErrors might differ between two printed tables!\nThe last one from RooFitResult.Print() should be correct.\nIf you want the errors to be reliable, opt for chi2_fit() method.\n' + 65*'~' + '\n\n')
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
        #
        m.minimize("Minuit2","minimize")
        for param in fix_float:
            param.setConstant(1)
        m.minimize("Minuit2","minimize")
        for param in fix_float:
            param.setConstant(0)
        m.minimize("Minuit2","minimize")
        m.minimize("Minuit2","minimize")
        self.is_fitted = True
        #
        if minos:
            if poi is None:
                raise TypeError('Poi is None by default: set it to a proper variable to run MINOS.')
            m.minos(ROOT.RooArgSet(poi))
        return self

    def plot_on_frame(self, title=' ', plot_params=ROOT.RooArgSet()):
        """Plot the instance model with all its components and data on the RooPlot frame

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
        #
        frame = ROOT.RooPlot(" ", title, self.var, var_left, var_right, var_nbins)  # frame.getAttText().SetTextSize(0.053)
        self.data.plotOn(frame, RF.DataError(ROOT.RooAbsData.Auto))
        self.model.plotOn(frame, RF.LineColor(ROOT.kRed-6), RF.LineWidth(5)) #, RF.NormRange('full'), RF.Range('full')
        self.model.paramOn(frame, RF.Layout(0.55, 0.96, 0.9), RF.Parameters(plot_params))
        #
        iter = self.model.getComponents().iterator()
        iter_comp = iter.Next()
        while iter_comp:
            if 'sig' in iter_comp.GetName().split('_'):
                self.model.plotOn(frame, RF.Components(iter_comp.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(47), RF.LineWidth(4))
            if 'bkgr' in iter_comp.GetName().split('_'):
                self.model.plotOn(frame, RF.Components(iter_comp.GetName()), RF.LineStyle(ROOT.kDashed), RF.LineColor(ROOT.kBlue-8), RF.LineWidth(4))
            iter_comp = iter.Next()
        #
        frame.GetYaxis().SetTitle(f'Candidates / {round((var_right - var_left) * 1000. / var_nbins, 1)} MeV')
        frame.GetXaxis().SetTitleSize(0.04)
        frame.GetYaxis().SetTitleSize(0.04)
        frame.GetXaxis().SetLabelSize(0.033)
        frame.GetYaxis().SetLabelSize(0.033)
        frame.GetXaxis().SetTitleOffset(1.05)
        frame.GetYaxis().SetTitleOffset(1.3)
        return frame

    def prepare_workspace(self, poi, nuisances):
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
            raise Exception('Model was not fitted to data, fit it first.')

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
    def extract(w):
        """Extract data, signal+background and background-only models from a given RooWorkspace.
        Background-only model is taken from the s+b model by setting the parameter of interest to be 0.

        NB: naming conventions assume that data's name is 'data', and that the parameter of interest is the first one and corresponds to the signal yield.

        Parameters:
        -----------

        w, RooWorkspace
            workspace with saved data and s+b model to unpack

        Returns:
        --------
        data, model_sb, model_b: RooAbsData, RooAbsPdf, RooAbsPdf
            Tuple with data, s+b and b-only models
        """
        data = w.obj("data")
        #
        model_sb = w.obj("ModelConfig")
        model_sb.SetName("model_sb")
        poi = model_sb.GetParametersOfInterest().first()
        #
        model_b = model_sb.Clone()
        model_b.SetName("B_only_model")
        oldval = poi.getVal()
        poi.setVal(0)
        model_b.SetSnapshot(ROOT.RooArgSet(poi))
        poi.setVal(oldval)
        return data, model_sb, model_b
