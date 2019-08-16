from DataExplorer import DataExplorer
import ROOT
from ROOT import RooFit as RF

N_GEN = 1000 # number of events to generate
SIG_FRACTION = 0.05
chi2_results = {}

# Observable
m = ROOT.RooRealVar("m", "m", 5.2, 5.3)
m.setBins(50)

# Parameters (you can tune them)
exp_par = ROOT.RooRealVar("exp_par", "exp_par", -30, -50, -10)
mean = ROOT.RooRealVar("mean", "", 5.25, 5.2, 5.3)
sigma = ROOT.RooRealVar("sigma", "", 0.005, 0.001, 0.05)

fraction = ROOT.RooRealVar('fraction', 'fraction', SIG_FRACTION)  # used only for data generation here
N_sig = ROOT.RooRealVar('N_sig', '', 100, 0, N_GEN)
N_bkgr = ROOT.RooRealVar('N_bkgr', '', 1000., 0, N_GEN)

# PDFs (note labelling for correct plotting)
sig = ROOT.RooGaussian("sig", "sig", m, mean, sigma)
bkgr = ROOT.RooExponential("bkgr", "bkgr", m, exp_par)
model_gen = ROOT.RooAddPdf("model_gen", "model_gen", ROOT.RooArgList(sig, bkgr), ROOT.RooArgList(fraction))
model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(sig, bkgr), ROOT.RooArgList(N_sig, N_bkgr))

# Sample N_GEN events
data = model_gen.generate(ROOT.RooArgSet(m), N_GEN)
data.reduce(f'{m.GetName()} > {m.getMin()} && {m.GetName()} < {m.getMax()}')

# Study and plot'em all
DE = DataExplorer(label='test', data=data, model=model)
fit_results = DE.chi2_fit(nbins = 50)
c = ROOT.TCanvas()
frame = DE.plot_on_frame(nbins=-1)
frame.Draw()

# Calculate statistical significance of signal observation
w = DE.write_to_workspace(poi=N_sig, nuisances= [exp_par, mean, sigma, N_bkgr])
# asympt_rrr = DE.asympt_signif_ll(w=w)
# DE.asympt_signif_ll(w=w) # another method
chi2_results = list(DE.chi2_test(pvalue_threshold=0.05, nbins=-1).values())[0]
print(f'\n\nchi2: {chi2_results[0]}\nndf: {chi2_results[1]}\np-value of chi2 test: {chi2_results[2]}\n')
print(f'fit status: {DE.fit_status}, chi2_test status: {DE.chi2_test_status}')

# c_ll = ROOT.TCanvas()
# frame_ll = DE.plot_ll(N_sig)
# frame_ll.Draw()
