import os,sys
import ROOT as rt

rt.gStyle.SetOptStat(0)

anafile = "output_crossingpt_ana_test.root"
if len(sys.argv)>=2:
    anafile = sys.argv[1]

rfile = rt.TFile(anafile, "OPEN")

treenames = [("prefilter","mcxingptana_prefilter"),("postfilter","mcxingptana_postfilter")]
hists = {}
endptnames = ["top","bottom","upstream","downstream","anode","cathode","imgends","tot"]

# ------------------------------------------------------------------
# Canvases

cpos = rt.TCanvas("cpos","Truth End point position", 1400, 900)
cpos.Divide(3,2)

ctype = rt.TCanvas("ctype","Truth End point type", 1400, 400)
ctype.Divide(2,1)

ctypeev = rt.TCanvas("ctypeperevent","Truth End point type per event", 1400, 600)
ctypeev.Divide(3,2)

cdist = rt.TCanvas("cdist","Truth Dist", 1800, 1800)
cdist.Divide(3,3)

ctypepurity = rt.TCanvas("ctypepurity","Reco End point type purity", 1400, 400)
ctypepurity.Divide(2,1)

# TRUTH CROSSING POINT POINTS
for name,treename in treenames:

    tree = rfile.Get(treename)

    # POS
    cpos.cd(1)
    hpos_matched = {}
    hpos_totaled = {}
    hpos_ratio   = {}
    posx = ["x","y","z"]
    hpos_matched[0] = rt.TH1D("hposx_matched_"+name, "",15,0,260)
    hpos_matched[1] = rt.TH1D("hposy_matched_"+name, "",15,-120,120)
    hpos_matched[2] = rt.TH1D("hposz_matched_"+name, "",15,0,1050)
    hpos_totaled[0] = rt.TH1D("hposx_totaled_"+name, "",15,0,260)
    hpos_totaled[1] = rt.TH1D("hposy_totaled_"+name, "",15,-120,120)
    hpos_totaled[2] = rt.TH1D("hposz_totaled_"+name, "",15,0,1050)
    hpos_ratio[0] = rt.TH1D("hposx_eff_"+name, "",15,0,260)
    hpos_ratio[1] = rt.TH1D("hposy_eff_"+name, "",15,-120,120)
    hpos_ratio[2] = rt.TH1D("hposz_eff_"+name, "",15,0,1050)
    for n,x in enumerate(posx):
        cpos.cd(n+1)        
        tree.Draw("pos[%d]>>hpos%s_matched_%s"%(n,x,name),"matched==1 && flashmatched==1")
        tree.Draw("pos[%d]>>hpos%s_totaled_%s"%(n,x,name),"matched>=0 && flashmatched==1")
        cpos.cd(3+n+1)
        tree.Draw("pos[%d]>>hpos%s_eff_%s"%(n,x,name),"matched==1 && flashmatched==1")
        hpos_ratio[n].Divide( hpos_totaled[n] )
    hists[(name,"pos")] = {}
    hists[(name,"pos")]["matched"] = hpos_matched
    hists[(name,"pos")]["total"] = hpos_totaled
    hists[(name,"pos")]["ratio"] = hpos_ratio

    # TYPE
    ctype.cd(1)
    htype_matched = rt.TH1D("htype_matched_"+name, "",7,0,7)
    htype_totaled = rt.TH1D("htype_totaled_"+name, "",7,0,7)
    htype_ratio   = rt.TH1D("htype_ratio_"+name,   "",7,0,7)
    tree.Draw("truth_type>>htype_matched_"+name,"matched==1 && flashmatched==1")
    tree.Draw("truth_type>>htype_totaled_"+name,"matched>=0 && flashmatched==1")
    tree.Draw("truth_type>>htype_ratio_"+name,"matched==1 && flashmatched==1")
    htype_ratio.Divide( htype_totaled )
    hists[(name,"type")] = {}
    hists[(name,"type")]["matched"] = htype_matched
    hists[(name,"type")]["totaled"] = htype_totaled
    hists[(name,"type")]["ratio"] = htype_ratio    
    
    # DIST
    hdist = {}
    for n in range(0,7):
        cdist.cd(n+1)
        hdist[n] = rt.TH1D("hdist_type%d_%s"%(n,name),"",50,0,200)
        tree.Draw("dist>>hdist_type%d_%s"%(n,name),"truth_type==%d && flashmatched==1"%(n))
    hdist[7] = rt.TH1D("hdist_tot_"+name,"",50,0,200)
    tree.Draw("dist>>hdist_tot_"+name,"flashmatched==1")
    hists[(name,"dist")] = hdist

# PLOT

# distance of true point to nearest reco
for n in range(0,8):
    cdist.cd(n+1)
    hpre = hists[("prefilter","dist")][n]
    hpre.SetTitle(endptnames[n])
    hpre.SetLineColor(rt.kBlack)
    hpre.SetLineWidth(2)
    hpre.Draw()
    hpost = hists[("postfilter","dist")][n]
    hpost.SetLineColor(rt.kRed)
    hpost.SetLineWidth(2)    
    hpost.Draw("same")
cdist.Update()
    


# position
cpos.Draw()
for n,x in enumerate(["x","y","z"]):
    cpos.cd(1+n)
    hists[("prefilter","pos")]["total"][n].Draw()
    hists[("prefilter","pos")]["matched"][n].SetLineColor(rt.kBlack)    
    hists[("prefilter","pos")]["matched"][n].Draw("same")
    hists[("postfilter","pos")]["matched"][n].SetLineColor(rt.kRed)
    hists[("postfilter","pos")]["matched"][n].Draw("same")
    cpos.cd(4+n)
    hists[("prefilter","pos")]["ratio"][n].GetYaxis().SetRangeUser(0,1.0)    
    hists[("prefilter","pos")]["ratio"][n].Draw()
    hists[("postfilter","pos")]["ratio"][n].SetLineColor(rt.kRed)
    hists[("postfilter","pos")]["ratio"][n].Draw("same")
cpos.Update()

# type
ctype.Draw()
ctype.cd(1)
hists[("prefilter","type")]["totaled"].SetLineColor(rt.kBlack)
hists[("prefilter","type")]["totaled"].Draw()
#hists[("prefilter","type")]["matched"].SetLineColor(rt.kBlue)
hists[("prefilter","type")]["matched"].Draw("same")
hists[("postfilter","type")]["matched"].SetLineColor(rt.kRed)
hists[("postfilter","type")]["matched"].Draw("same")
ctype.cd(2)
hists[("prefilter","type")]["ratio"].GetYaxis().SetRangeUser(0,1.0)
hists[("prefilter","type")]["ratio"].Draw()
hists[("postfilter","type")]["ratio"].SetLineColor(rt.kRed)
hists[("postfilter","type")]["ratio"].Draw("same")
ctype.Update()

print "Filled Truth-level crossings"
    
# RECO
rtree = rfile.Get( "recoxingptana" )

# RECO TYPE
crecotype = rt.TCanvas("crecotype","Purity for Reco End Points", 1400, 400)
crecotype.Draw()
crecotype.Divide(2,1)
# fill hists
crecotype.cd(1)
hrecotype_totaled  = rt.TH1D("hrecotype_totaled", "",7,0,7)
hrecotype_matched  = rt.TH1D("hrecotype_matched", "",7,0,7)
hrecotype_filtered = rt.TH1D("hrecotype_filtered", "",7,0,7)
hrecotype_filtpass = rt.TH1D("hrecotype_filtpass", "",7,0,7)
hrecotype_ratio    = rt.TH1D("hrecotype_ratio",   "",7,0,7)
hrecotype_ratiofil = rt.TH1D("hrecotype_ratiofil","",7,0,7)
rtree.Draw("type>>hrecotype_totaled")
rtree.Draw("type>>hrecotype_matched", "matched==1") # matched to truth xingpt
rtree.Draw("type>>hrecotype_filtered", "passes==1") # all that passes
rtree.Draw("type>>hrecotype_filtpass", "matched==1 && passes==1") # matched and passes
rtree.Draw("type>>hrecotype_ratio",   "matched==1")
rtree.Draw("type>>hrecotype_ratiofil", "matched==1 && passes==1")

hrecotype_totaled.SetMinimum(0)
hrecotype_totaled.Draw()
hrecotype_matched.Draw("same")
hrecotype_matched.SetLineColor(rt.kRed)
hrecotype_filtered.SetLineColor(rt.kMagenta)
hrecotype_filtpass.SetLineColor(rt.kBlue)
hrecotype_ratiofil.SetLineColor(rt.kRed )
hrecotype_ratio.Divide( hrecotype_totaled )
hrecotype_ratiofil.Divide( hrecotype_filtered )

trecolen = rt.TLegend(0.1,0.8, 0.5, 0.6)
hrecotype_totaled.SetTitle("Number of Reco. Points")
hrecotype_totaled.Draw()
hrecotype_matched.Draw("same")
hrecotype_filtered.Draw("same")
hrecotype_filtpass.Draw("same")
trecolen.AddEntry( hrecotype_totaled,"All Reco.","L")
trecolen.AddEntry( hrecotype_matched,"Truth Matched","L")
trecolen.AddEntry( hrecotype_filtered,"Filtered","L")
trecolen.AddEntry( hrecotype_filtpass,"Filtered+Matched","L")
trecolen.Draw()

# purity
crecotype.cd(2)
hrecotype_ratio.GetYaxis().SetRangeUser(0,1)
hrecotype_ratio.SetTitle("Reco. End Pt. Purity")
hrecotype_ratio.Draw()
hrecotype_ratiofil.Draw("same")
crecotype.Update()
crecotype.SaveAs("recoxingpt_types.png")    

# DIST
# cdist = rt.TCanvas("cdist","Truth Dist", 1800, 1800)
# cdist.Draw()
# cdist.Divide(3,3)
# hdist = {}
# for n in range(0,7):
#     cdist.cd(n+1)
#     hdist[n] = rt.TH1D("hdist_type%d"%(n),"",50,0,200)
#     tree.Draw("dist>>hdist_type%d"%(n),"truth_type==%d"%(n))
#     hdist[7] = rt.TH1D("hdist_tot","",50,0,200)
#     tree.Draw("dist>>hdist_tot")

# for n in range(0,8):
#     cdist.cd(n+1)
#     hdist[n].Draw()
#     cdist.Update()
#     cdist.SaveAs("truthxingpt_%s_dist.png"%(name))    
 


raw_input()
