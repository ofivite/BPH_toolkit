from DataExplorer import DataExplorer
import ROOT
from ROOT import RooFit as RF

N_GEN = 1000 # number of events to generate
SIG_FRACTION = 0.1 # fraction of the signal component
chi2_results = {}

# Observable
m = ROOT.RooRealVar("m", "mass [GeV]", 5.2, 5.3)
m.setBins(50)

# Parameters
exp_par = ROOT.RooRealVar("exp_par", "#lambda", -5., -50, 1.)
mean = ROOT.RooRealVar("mean", "mean [GeV]", 5.25, 5.2, 5.3)
sigma = ROOT.RooRealVar("sigma", "#sigma [GeV]", 0.005, 0.001, 0.05)

fraction = ROOT.RooRealVar('fraction', 'fraction', SIG_FRACTION)  # used only for data generation here
N_sig = ROOT.RooRealVar('N_sig', 'N_{sig}', 100, 0, N_GEN)
N_bkgr = ROOT.RooRealVar('N_bkgr', 'N_{bkgr}', 1000., 0, N_GEN)

# PDFs (note labelling for correct plotting)
sig = ROOT.RooGaussian("sig", "sig", m, mean, sigma)
bkgr = ROOT.RooExponential("bkgr", "bkgr", m, exp_par)
model_gen = ROOT.RooAddPdf("model_gen", "model_gen", ROOT.RooArgList(sig, bkgr), ROOT.RooArgList(fraction))
model = ROOT.RooAddPdf('model', 'model', ROOT.RooArgList(sig, bkgr), ROOT.RooArgList(N_sig, N_bkgr))

# Sample N_GEN events
data = model_gen.generate(ROOT.RooArgSet(m), N_GEN)
data = data.reduce(f'{m.GetName()} > {m.getMin()} && {m.GetName()} < {m.getMax()}')

# Fit and plot'em all
DE = DataExplorer(label='test', data=data, model=model)
fit_results = DE.fit(minos=True)
c = ROOT.TCanvas()
frame = DE.plot_on_frame()
frame.Draw()

# Make pull distribution for the fit
c_pull = ROOT.TCanvas()
frame_pull = DE.plot_pull()
frame_pull.Draw()

# Plot likelihood profiles
c_ll = ROOT.TCanvas()
frame_ll = DE.plot_ll(poi=N_sig)
frame_ll.Draw()

# Check whether there is bias in the fit
frame_var, frame_err, frame_pull = DE.check_fit_bias(param_to_study=N_sig, N_toys=100)
c_var = ROOT.TCanvas()
frame_var.Draw()
c_error = ROOT.TCanvas()
frame_err.Draw()
c_pull = ROOT.TCanvas()
frame_pull.Draw()

# Export data and model to workspace
w = DE.write_to_workspace(poi=N_sig, nuisances= [exp_par, mean, sigma, N_bkgr])

# Calculate statistical significance of signal observation with asymptotic approximation
asympt_rrr = DE.asympt_signif(w=w)
# DE.asympt_signif_ll(w=w) # another method

# Do chi^2 goodness-of-fit test and print fit and test status
chi2_results = list(DE.chi2_test(pvalue_threshold=0.05, nbins=-1).values())[0]
print(f'\n\nchi2: {chi2_results[0]}\nndf: {chi2_results[1]}\np-value of chi2 test: {chi2_results[2]}\n')
print(f'fit status: {DE.fit_status}, chi2_test status: {DE.chi2_test_status}')

# Calculate signal significance with toys
# Note that p-value might be zero if bkgr hypo is very unlikely since it might require lots of toys to get to the observed value of test statistic.
c_toy_signif = ROOT.TCanvas()
frame_toy_signif, _ = DE.toy_signif(w, n_toys_null = 100, n_toys_alt=10)
frame_toy_signif.Draw()
