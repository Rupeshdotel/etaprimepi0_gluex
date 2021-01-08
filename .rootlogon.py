# Content of .rootlogon.py
import ROOT
myStyle = ROOT.TStyle('MyStyle','My graphics style')
myStyle.SetCanvasColor(ROOT.kBlue) # My canvases are blue!
ROOT.gROOT.SetStyle('MyStyle')
