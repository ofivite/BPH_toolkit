import ROOT
from ROOT import RooFit as RF
from misc import interactivity_yn
from scipy.stats import chi2
from pandas import DataFrame
from math import sqrt

class StatTools:
    """Additional to DataExplorer module with methods for performing statistical inference.
    """
    def chi2_test(self, pvalue_threshold=0.05, nbins=-1):
        """Make goodness-of-fit chi2 test between the instance's data and model.
        NB: by default, binning is taken from the variable's definition. Otherwise, it is temporarily set to nbins value.

        Parameters
        ----------
        pvalue_threshold: float, optional (default=0.05)
            threshold for setting boolean flag self.chi2_test_status (pass/fail chi2 test)
        nbins: int/float, optional (default=-1: take the number of bins from the variable's definition)
            number of bins in calculating chi2

        Returns
        -------
        dict: Python dictionary with chi2, ndf and p-value of the test.
        """
        init_nbins = self.var.numBins()
        if nbins != -1:
            assert (nbins % 1 == 0 and nbins >= 0), 'nbins must be a positive integer'
            self.var.setBins(nbins)
        if not self.is_fitted:
            raise Exception('Model was not fitted to data, fit it first.')
        is_extended = self.model.canBeExtended()
        data_hist = ROOT.RooDataHist('data_hist', 'data_hist', ROOT.RooArgSet(self.var), self.data) # binning is taken from self.var definition
        nfloat = self.model.getParameters(self.data).selectByAttrib("Constant", ROOT.kFALSE).getSize()
        ndf = self.var.numBins() - nfloat
        chi2_var = ROOT.RooChi2Var("chi2_var","chi2_var", self.model, data_hist, RF.Extended(is_extended), RF.DataError(ROOT.RooAbsData.Auto))
        chi2_value = chi2_var.getVal()
        pvalue = 1 - chi2.cdf(chi2_value, ndf)
        self.chi2_test_status = 0 if pvalue > pvalue_threshold else 1
        self.var.setBins(init_nbins)
        return {f'{self.label}_{self.data.GetName()}': [chi2_value, ndf, pvalue]}

    @classmethod
    def asympt_signif(cls, w, print_level=-1):
        """Analytically calculate one-sided signal significance for a given in the workspace s+b model using RooStats.AsymptoticCalculator.

        Parameters
        ----------
        w: RooWorkspace
            workspace containing data and model to be opened. Usage of the method extract_from_workspace is assumed, so the workspace should contain data, s+b and b modelconfigs.
        print_level: int, optional (default=-1)
            verbosity of the output

        Returns
        -------
        as_result: HypoTestResult
            results of the test (for printing use Print() method)
        """
        data, mc_sb, mc_b = cls.extract_from_workspace(w)
        if data.isWeighted():
            interactivity_yn('It\'s not a good idea to do asymptotic significance calculation with weighted data. Sure you want to proceed?')
        ac = ROOT.RooStats.AsymptoticCalculator(data, mc_sb, mc_b)
        ac.SetPrintLevel(print_level)
        ac.SetOneSidedDiscovery(True)
        as_result = ac.GetHypoTest()
        as_result.Print()
        return as_result

    @classmethod
    def asympt_signif_ll(cls, w):
        """Analytically calculate one-sided significance for a given in the workspace s+b model by bare hands (through likelihoods). Might be useful as a cross-check to asympt_signif().
        NB: this gives more control on fitting procedure than asympt_signif() method

        Parameters
        ----------
        w: RooWorkspace
            workspace containing data and model to be opened. Usage of the method extract_from_workspace is assumed, so the workspace should contain data, s+b and b modelconfigs.

        Returns
        -------
        dict: dictionary
            results of the test: p-value, Z-value (signal significance), negative loglikelihood for s+b and b hypotheses respectively
        """
        data, mc_sb, mc_b = cls.extract_from_workspace(w)
        if data.isWeighted():
            interactivity_yn('It\'s not a good idea to do asymptotic significance calculation with weighted data. Sure you want to proceed?')
        num_poi = mc_sb.GetParametersOfInterest().getSize()
        if num_poi != 1:
            print(f'This implementation works only for one parameter of interest (you have {num_poi}). Either change pois in the workspace or see 1007.1727 paper.')
            return

        mc_sb.LoadSnapshot();
        model_sb = mc_sb.GetPdf()
        DE_sb = cls(label='sb', data=data, model=model_sb)
        rrr_sig = DE_sb.fit(is_sum_w2=data.isWeighted())

        mc_b.LoadSnapshot()
        mc_b.GetParametersOfInterest().first().setConstant()
        model_b = mc_b.GetPdf()
        DE_b = cls(label='b', data=data, model=model_b)
        rrr_null = DE_b.fit(is_sum_w2=data.isWeighted())

        nll_sig  = rrr_sig.minNll()
        nll_null = rrr_null.minNll()
        P = ROOT.TMath.Prob(2*(nll_null - nll_sig), 1) ## !!! should be always ndf = 1 = number of poi for this formula to work
        S = ROOT.TMath.ErfcInverse(P)*sqrt(2) # this yields same result as AsymptoticCalculator
        # S = ROOT.Math.gaussian_quantile_c(P, 1) # this is slightly different, might be python precision issues
        return {'P': P, 'S': S, 'nll_sig': nll_sig, 'nll_null': nll_null}

    @classmethod
    def toy_signif(cls, w, n_toys_null=1000, n_toys_alt=100, seed=333, save=False, save_folder='.', save_prefix='toy_signif'):
        """Calculate signal significance by generating toy samples to get the test statistic distribution under null and alt hypotheses. Use one-sided ProfileLikelihoodTestStat().

        Parameters
        ----------
        w: RooWorkspace
            workspace containing data and model to be opened. Usage of the method extract_from_workspace is assumed, so the workspace should contain data, s+b and b ModelConfigs.
        n_toys_null: int, optional (default=1000)
            number of toy samples to be generated for the null (b) hypothesis
        n_toys_alt: int, optional (default=100)
            number of toy samples to be generated for the alternative (s+b) hypothesis
        seed: int, optional (default=333)
            random seed for the generator
        save: bool, optional (default=False)
            whether save the plot or not
        save_folder: str, optional (default='.')
            path for saving
        save_prefix: str, optional (default='toy_signif')
            prefix to the file name

        Returns
        -------
        frame: RooPlot, HypoTestResult
            hypothesis test results (use Print() for printing them) and frame containing the test statistic plots (requires further drawing on the canvas)
        """
        ROOT.RooRandom.randomGenerator().SetSeed(seed)
        data, mc_sb, mc_b = cls.extract_from_workspace(w)
        fc = ROOT.RooStats.FrequentistCalculator(data, mc_sb, mc_b)
        fc.SetToys(n_toys_null, n_toys_alt); # fc.SetNToysInTails(500, 100)
        profll = ROOT.RooStats.ProfileLikelihoodTestStat(mc_sb.GetPdf())
        profll.SetOneSidedDiscovery(True)
        toymcs = ROOT.RooStats.ToyMCSampler(fc.GetTestStatSampler())
        toymcs.SetTestStatistic(profll)
        if not mc_sb.GetPdf().canBeExtended():
            toymcs.SetNEventsPerToy(1)
            print ('adjusting for non-extended formalism')
        fqResult = fc.GetHypoTest()
        fqResult.Print()
        plot_hypo_test = ROOT.RooStats.HypoTestPlot(fqResult)
        plot_hypo_test.SetLogYaxis(True)
        if save:
            c = ROOT.TCanvas()
            plot_hypo_test.Draw()
            c.Draw()
            c.SaveAs(f'{save_folder}/{save_prefix}.pdf')
        return plot_hypo_test, fqResult

    def toy_tstat(self, n_toys=1000, seed=333, save=False):
        ROOT.RooRandom.randomGenerator().SetSeed(seed)
        t_list = []
        #
        self.poi.setVal(0); self.poi.setConstant(1);
        self.chi2_fit() ### assuming the data to be RooDataHist()
        data_TH1 = self.data.createHistogram(self.var.GetName())
        max_error = max([data_TH1.GetBinError(i) for i in range(1, data_TH1.GetNbinsX() + 1)])

        for i_toy in range(n_toys):
            self.poi.setVal(0); self.poi.setConstant(1);
            toy_data = self.model.generate(ROOT.RooArgSet(self.var), self.data.sumEntries())
            toy_roohist = ROOT.RooDataHist('toy_TH1', 'toy_TH1', ROOT.RooArgSet(self.var), toy_data)
            toy_TH1 = toy_roohist.createHistogram(self.var.GetName())
            for i in range(1, toy_TH1.GetNbinsX()+1):
                toy_TH1.SetBinError(i, max_error)
            toy_hist = ROOT.RooDataHist('toy_hist', 'toy_hist', ROOT.RooArgList(self.var), RF.Import(toy_TH1))
            #
            self.chi2_fit(data=toy_hist)
            chi2_null = ROOT.RooChi2Var("chi2_null","chi2_null", self.model, toy_hist, RF.Extended(False), RF.DataError(ROOT.RooAbsData.Auto))
            chi2_null = chi2_null.getVal()
            #
            if save:
                c_null = ROOT.TCanvas("c_null", "c_null", 800, 600);
                frame_null = self.var.frame()
                toy_hist.plotOn(frame_null, RF.DataError(ROOT.RooAbsData.SumW2))
                self.model.plotOn(frame_null)
                frame_null.Draw()
                c_null.SaveAs(f'./tnull_toys/null_{i_toy}.pdf')
            #
            self.poi.setConstant(0);
            self.chi2_fit(data=toy_hist)
            chi2_sb = ROOT.RooChi2Var("chi2_sb","chi2_sb", self.model, toy_hist, RF.Extended(False), RF.DataError(ROOT.RooAbsData.Auto))
            chi2_sb = chi2_sb.getVal()
            t_list.append([chi2_null - chi2_sb, chi2_null, chi2_sb])

        if save:
            c_sb = ROOT.TCanvas("c_sb", "c_sb", 800, 600);
            frame_sb = self.var.frame()
            toy_hist.plotOn(frame_sb, RF.DataError(ROOT.RooAbsData.SumW2))
            self.model.plotOn(frame_sb)
            frame_sb.Draw()
            c_sb.SaveAs(f'./tnull_toys/sb_{i_toy}.pdf')

        df = DataFrame(t_list, columns=['t', 'chi2_null', 'chi2_sb'])
        df.to_pickle('t_.pkl')
        return df

    def plot_ll(self, poi, nbins=100, poi_min=-1, poi_max=-1, save=False, save_folder='.', save_prefix='pll'):
        """Plot the nominal and profiled likelihoods for the instance's data and model for the provided parameter of interest

        Parameters
        ----------
        poi: RooRealVar
            parameter of interest for which the likelihoods will be plotted
        nbins: int, optional (default=100)
            number of bins for plotting
        poi_min: float, optional (default=-1, set to be -5*poi_fit_error)
            left range for plotting
        poi_max: float, optional (default=-1, set to be +5*poi_fit_error)
            right range for plotting
        save: bool, optional (default=False)
            whether save the plot or not
        save_folder: str, optional (default='.')
            path for saving
        save_prefix: str, optional (default='pll')
            prefix to the file name

        Returns
        -------
        frame: RooPlot
            frame containing the likelihood plots (requires further drawing on the canvas)
        """
        poi_from_model = self.model.getVariables().find(poi.GetName())
        if poi_min == -1: poi_min = poi_from_model.getVal() - 5*poi_from_model.getError();
        if poi_max == -1: poi_max = poi_from_model.getVal() + 5*poi_from_model.getError();
        assert (nbins % 1 == 0 and nbins >= 0), 'nbins must be a positive integer'
        assert (poi_min < poi_max), 'left range value must be lower than right range one'

        nll = self.model.createNLL(self.data)
        pll = nll.createProfile(ROOT.RooArgSet(poi))
        frame_nll = poi.frame(RF.Bins(nbins), RF.Range(poi_min, poi_max))
        frame_nll.SetTitle('')
        nll.plotOn(frame_nll, RF.ShiftToZero(), RF.LineColor(ROOT.kGreen))
        pll.plotOn(frame_nll, RF.LineColor(ROOT.kRed))
        frame_nll.SetMaximum(25.)
        frame_nll.SetMinimum(0.)
        frame_nll.SetXTitle(poi.GetName())
        if save:
            c_ll = ROOT.TCanvas("c_ll", "c_ll", 800, 600)
            frame_nll.Draw()
            c_ll.SaveAs(f'{save_folder}/{save_prefix}.pdf')
        return frame_nll

    def plot_pull(self, save=False, save_folder='.', save_prefix='pull'):
        """Plot the pull between data and model.

        Parameters
        ----------
        save: bool, optional (default=False)
            whether save the plot or not
        save_folder: str, optional (default='.')
            path for saving
        save_prefix: str, optional (default='pull')
            prefix to the file name

        Returns
        -------
        frame: RooPlot
            frame containing the pull plot (requires further drawing on the canvas)
        """
        frame = self.var.frame()
        self.data.plotOn(frame)
        self.model.plotOn(frame)
        pull_hist = frame.pullHist()
        frame_pull = self.var.frame()
        frame_pull.addPlotable(pull_hist, 'P')

        line_null = ROOT.TLine(self.var.getMin(), 0, self.var.getMax(), 0)
        line_null.SetLineColor(ROOT.kRed-6);
        line_null.SetLineWidth(4);
        frame_pull.addObject(line_null)
        if frame_pull.GetMaximum() > 4:
            line_plus_3sigma = ROOT.TLine(self.var.getMin(), 3, self.var.getMax(), 3)
            line_plus_3sigma.SetLineColor(ROOT.kRed-6);
            line_plus_3sigma.SetLineWidth(4);
            line_plus_3sigma.SetLineStyle(2);
            frame_pull.addObject(line_plus_3sigma)
        if frame_pull.GetMinimum() < -4:
            line_minus_3sigma = ROOT.TLine(self.var.getMin(), -3, self.var.getMax(), -3)
            line_minus_3sigma.SetLineWidth(4)
            line_minus_3sigma.SetLineColor(ROOT.kRed-6)
            line_minus_3sigma.SetLineStyle(2)
            frame_pull.addObject(line_minus_3sigma)
        if save:
            c_pull = ROOT.TCanvas("c_pull", "c_pull", 800, 600)
            frame_pull.Draw()
            c_pull.SaveAs(f'{save_folder}/{save_prefix}_{self.var.GetName()}.pdf')
        return frame_pull

    def check_fit_bias(self, param_to_study, N_toys=1000, verbose=False, save=False, save_folder='.', save_prefix='bias_check'):
        """Using RooMCStudy() class make bias checks in fitted model's parameter param_to_study by repitative sampling of toys from the model and then fitting them.
        Normally, if mean or sigma of pull plot (should look like Gaussian) are within 3sigma from 0 and 1 respectively, the fit is not biased.
        NB: Fit extendability and weights presence are taken into account automatically.

        Parameters
        ----------
        param_to_study: RooRealVar
            variable for which the fit results will be accumulated
        N_toys: int, optional (default=1000)
            number of toy samples to be generated
        verbose: bool, optional (default=False)
            verbosity of the output
        save: bool, optional (default=False)
            whether save the plot or not
        save_folder: str, optional (default='.')
            path for saving
        save_prefix: str, optional (default='bias_check')
            prefix to the file name

        Returns
        -------
        frame: RooPlot, RooPlot, RooPlot
            frames with distributions for variable's fitted value, error and pull respectively (requires further drawing on the canvas)
        """
        if not self.is_fitted:
            raise Exception('Model was not fitted to data, fit it first.')
        is_extended = self.model.canBeExtended()
        is_sum_w2 = self.data.isWeighted()
        MC_manager = ROOT.RooMCStudy(self.model, ROOT.RooArgSet(self.var), RF.Extended(is_extended), RF.FitOptions(RF.SumW2Error(is_sum_w2), RF.Verbose(verbose)))
        MC_manager.generateAndFit(N_toys)

        frame_var = param_to_study.frame()
        MC_manager.plotParamOn(frame_var)
        frame_err = MC_manager.plotError(param_to_study)
        frame_pull = MC_manager.plotPull(param_to_study, -3, 3, 60, ROOT.kTRUE)

        if save:
            c_var = ROOT.TCanvas("c_var", "c_var", 800, 600)
            frame_var.Draw()
            c_err = ROOT.TCanvas("c_err", "c_err", 800, 600)
            frame_err.Draw()
            c_pull = ROOT.TCanvas("c_pull", "c_pull", 800, 600)
            frame_pull.Draw()
            c_var. SaveAs(f'{save_folder}/{save_prefix}_{self.var.GetName()}.pdf')
            c_err. SaveAs(f'{save_folder}/{save_prefix}_{self.var.GetName()}err.pdf')
            c_pull.SaveAs(f'{save_folder}/{save_prefix}_{self.var.GetName()}pull.pdf')
        return frame_var, frame_err, frame_pull
