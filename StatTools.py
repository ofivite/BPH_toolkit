import ROOT
from ROOT import RooFit as RF
from misc import interactivity_yn
from scipy.stats import chi2
from pandas import DataFrame
from math import sqrt

class StatTools:
    """Additional to DataExplorer module for performing statistical inference.
    """

    def chi2_test(self, pvalue_threshold = 0.05):
        """Make goodness-of-fit chi2 test between the instance's data and model.
        NB: binning is taken from the variable's definition.

        Parameters
        ----------

        pvalue_threshold: float, optional (default=0.05)
            threshold for setting boolean flag self.chi2_test_status (pass/fail chi2 test)

        Returns
        -------
        dict: Python dictionary with chi2, ndf and p-value of the test.
        """
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
        return {f'{self.label}_{self.data.GetName()}': [chi2_value, ndf, pvalue]}

    @classmethod
    def asympt_signif(cls, w):
        """Function to calculate one-sided significance for a given in the workspace s+b model using RooStats.AsymptoticCalculator.

        Parameters
        ----------

        w: RooWorkspace
            Workspace containing data and model to be opened. Usage of the method extract_from_workspace is assumed.

        Returns
        -------
        as_result, HypoTestResult
            Results of the test (for printing use Print() method)
        """
        data, mc_sb, mc_b = cls.extract_from_workspace(w)
        if data.isWeighted():
            interactivity_yn('It\'s not a good idea to do asymptotic significance calculation with weighted data. Sure you want to proceed?')
        ac = ROOT.RooStats.AsymptoticCalculator(data, mc_sb, mc_b)
        ac.SetOneSidedDiscovery(True)
        as_result = ac.GetHypoTest()
        as_result.Print()
        return as_result

    @classmethod
    def asympt_signif_ll(cls, w):
        """Function to calculate one-sided significance for a given in the workspace s+b model by bare hands (through likelihoods). Might be useful as a cross-check to asympt_signif().
        NB: this gives more control on fitting than AsymptoticCalculator
        Parameters
        ----------

        Returns
        -------
        as_result, HypoTestResult
            Results of the test (for printing use Print() method)
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
        print ('P=', P, ' nll_sig=', nll_sig, ' nll_null=', nll_null, '\n', 'S=', S)

    @staticmethod
    def toy_signif(w, n_toys = 1000, seed = 333):
        """
        // to be completed //
        """
        data, mc_sb, mc_b = self.extract_from_workspace(w)
        fc = ROOT.RooStats.FrequentistCalculator(data, mc_sb, mc_b)
        fc.SetToys(n_toys, n_toys/10); # fc.SetNToysInTails(500, 100)
        profll = ROOT.RooStats.ProfileLikelihoodTestStat(mc_sb.GetPdf())
        profll.SetOneSidedDiscovery(True)
        toymcs = ROOT.RooStats.ToyMCSampler(fc.GetTestStatSampler())
        toymcs.SetTestStatistic(profll)
        #
        if not mc_sb.GetPdf().canBeExtended():
            toymcs.SetNEventsPerToy(1)
            print ('adjusting for non-extended formalism')
        #
        fqResult = fc.GetHypoTest()
        fqResult.Print()
        #
        c = ROOT.TCanvas()
        plot = ROOT.RooStats.HypoTestPlot(fqResult)
        plot.SetLogYaxis(True)
        plot.Draw()
        c.Draw()
        c.SaveAs(f'./toy_signif_{mc_sb.GetPdf().GetName()}.pdf')
        return fqResult

    def toy_tstat(self, n_toys = 1000, seed = 333, save=False):
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

    def plot_ll(self, save=False, save_prefix = ''):
        nll = self.model.createNLL(self.data)
        pll = nll.createProfile(ROOT.RooArgSet(self.poi))
        #
        frame_nll = self.poi.frame(RF.Bins(100), RF.Range(ll_left, ll_right))
        frame_nll.SetTitle('')
        #
        nll.plotOn(frame_nll, RF.ShiftToZero(), RF.LineColor(ROOT.kGreen))
        pll.plotOn(frame_nll, RF.LineColor(ROOT.kRed))
        #
        frame_nll.SetMaximum(25.)
        frame_nll.SetMinimum(0.)
        if save:
            c_ll = ROOT.TCanvas("c_ll", "c_ll", 800, 600); ll_left = 0; ll_right = 200
            frame_nll.Draw()
            #
            line_width = 4
            line_5sigma = ROOT.TLine(ll_left, 12.5, ll_right, 12.5)
            line_5sigma.SetLineWidth(line_width); line_5sigma.SetLineColor(47)
            line_5sigma.Draw();
            #
            c_ll.SaveAs(f'{save_prefix}1_pll.pdf')
        return frame_nll

    def plot_pull(self, save=False, save_path='./fit_validation/', save_prefix = ''):
        frame = self.var.frame()
        self.data.plotOn(frame)
        self.model.plotOn(frame)
        pull_hist = frame.pullHist()
        #
        frame_pull = self.var.frame()
        frame_pull.addPlotable(pull_hist, 'P')
        if save:
            c_pull = ROOT.TCanvas("c_pull", "c_pull", 800, 600)
            frame_pull.Draw()
            c_pull.SaveAs(f'{save_path}{save_prefix}_{self.data.GetName()}.pdf')
        return frame_pull

    def plot_toys_pull(self, var_to_study, N_toys=100, N_gen=1, save=False, save_path='./fit_validation/', save_prefix = ''):
        """Make bias checks in fitted model parameter var_to_study by generating toys with RooMCStudy()
        """
        if not self.is_fitted:
            raise Exception('Model was not fitted to data, fit it first.')

        # width_N = 80 if self.mode == 'X' else 250
        # err_upper = 30 if self.mode == 'X' else 400; err_nbins = 30
        # var_lower = var_to_study.getVal() - width_N; var_upper = var_to_study.getVal() + width_N; var_nbins = 50
        #
        MC_manager = ROOT.RooMCStudy(self.model, ROOT.RooArgSet(self.var), RF.Extended(True), RF.FitOptions('mvl'))
        MC_manager.generateAndFit(N_toys, N_gen)
        #
        frame_var = var_to_study.frame() #var_lower, var_upper, var_nbins
        MC_manager.plotParamOn(frame_var)
        frame_err = MC_manager.plotError(var_to_study)
        frame_pull = MC_manager.plotPull(var_to_study, -3, 3, 60, ROOT.kTRUE)
        #
        if save:
            c_var = ROOT.TCanvas("c_var", "c_var", 800, 600)
            frame_var.Draw()
            c_err = ROOT.TCanvas("c_err", "c_err", 800, 600)
            frame_err.Draw()
            c_pull = ROOT.TCanvas("c_pull", "c_pull", 800, 600)
            frame_pull.Draw()
            #
            c_var. SaveAs(f'{save_path}{save_prefix}_{self.var.GetName()}.pdf')
            c_err. SaveAs(f'{save_path}{save_prefix}_{self.var.GetName()}err.pdf')
            c_pull.SaveAs(f'{save_path}{save_prefix}_{self.var.GetName()}pull.pdf')
        return MC_manager
