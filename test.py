from DataExplorer import DataExplorer
import ROOT


# Observable
m = ROOT.RooRealVar("m", "m", 5.20, 5.30)

# Parameters
m0 = ROOT.RooRealVar("m0", "m0", 5.291, 5.20, 5.30)
k = ROOT.RooRealVar("k", "k", -30, -50, -10)

# PDF
argus = ROOT.RooArgusBG("argus", "argus", m, m0, k)

# Sample 1000 events
data = argus.generate(ROOT.RooArgSet(m), 1000)

c = ROOT.TCanvas()
DE = DataExplorer(label='test', data=data, model=argus)
fit_results = DE.fit(is_sum_w2=False)
frame = DE.plot_on_frame()
frame.Draw()
