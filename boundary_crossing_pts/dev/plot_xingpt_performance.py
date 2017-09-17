import os,sys
import ROOT as rt

anafile = "output_crossingpt_ana.root"
if len(sys.argv)>=2:
    anafile = sys.argv[1]

rfile = rt.TFile(anafile, "OPEN")

treenames = ["mcxingptana_prefilter","mcxingptana_postfilter"]

# TRUTH CROSSING POINT POINTS
for name in treenames:

    tree = rfile.Get(name)

    # POS
    cpos = rt.TCanvas("cpos","[%s] Truth End point position"%(name), 1400, 900)
    cpos.Divide(3,2)
    cpos.cd(1)
    cpos.Draw()
    hpos_matched = {}
    hpos_totaled = {}
    hpos_ratio   = {}
    posx = ["x","y","z"]
    hpos_matched[0] = rt.TH1D("hposx_matched", "",15,0,260)
    hpos_matched[1] = rt.TH1D("hposy_matched", "",15,-120,120)
    hpos_matched[2] = rt.TH1D("hposz_matched", "",15,0,1050)
    hpos_totaled[0] = rt.TH1D("hposx_totaled", "",15,0,260)
    hpos_totaled[1] = rt.TH1D("hposy_totaled", "",15,-120,120)
    hpos_totaled[2] = rt.TH1D("hposz_totaled", "",15,0,1050)
    for n,x in enumerate(posx):    
        tree.Draw("pos[%d]>>hpos%s_matched"%(n,x),"matched==1")
        tree.Draw("pos[%d]>>hpos%s_totaled"%(n,x),"matched>=0")
    for n,x in enumerate(posx):
        cpos.cd(n+1)
        hpos_totaled[n].SetMinimum(0)
        hpos_totaled[n].Draw()
        hpos_matched[n].Draw("same")
        hpos_matched[n].SetLineColor(rt.kRed)
        hpos_ratio[n] = hpos_matched[n].Clone("hpos%s_ratio"%(x))
        hpos_ratio[n].Divide( hpos_totaled[n] )
        cpos.cd(4+n)
        hpos_ratio[n].GetYaxis().SetRangeUser(0,1)    
        hpos_ratio[n].Draw()
        cpos.Update()
        cpos.SaveAs("truthxingpt_%s_pos.png"%(name))


    # TYPE
    ctype = rt.TCanvas("ctype","[%s] Truth End point type"%(name), 1400, 400)
    ctype.Draw()
    ctype.Divide(2,1)
    ctype.cd(1)
    htype_matched = rt.TH1D("htype_matched", "",7,0,7)
    htype_totaled = rt.TH1D("htype_totaled", "",7,0,7)
    htype_ratio   = rt.TH1D("htype_ratio",   "",7,0,7)
    tree.Draw("truth_type>>htype_matched","matched==1 && flashmatched==1")
    tree.Draw("truth_type>>htype_totaled","matched>=0 && flashmatched==1")
    tree.Draw("truth_type>>htype_ratio","matched==1 && flashmatched==1")
    htype_totaled.SetMinimum(0)
    htype_totaled.Draw()
    htype_matched.Draw("same")
    htype_matched.SetLineColor(rt.kRed)
    htype_ratio.Divide( htype_totaled )

    htype_totaled.Draw()
    htype_matched.Draw("same")
    
    ctype.cd(2)
    htype_ratio.GetYaxis().SetRangeUser(0,1)
    htype_ratio.Draw()
    ctype.Update()
    ctype.SaveAs("truthxingpt_%s_truthtypes.png"%(name))    

    # TYPE: Top per event
    ctype = rt.TCanvas("ctypeperevent","[%s] Truth End point type per event"%(name), 1400, 600)
    ctype.Draw()
    ctype.Divide(3,2)
    ctype.cd(1)
    htype_top = rt.TH1D("htype_top", "",30,0,30)
    tree.Draw("truth_type>>htype_matched","matched==1 && flashmatched==1")
    tree.Draw("truth_type>>htype_totaled","matched>=0 && flashmatched==1")
    tree.Draw("truth_type>>htype_ratio","matched==1 && flashmatched==1")
    htype_totaled.SetMinimum(0)
    htype_totaled.Draw()
    htype_matched.Draw("same")
    htype_matched.SetLineColor(rt.kRed)
    htype_ratio.Divide( htype_totaled )

    htype_totaled.Draw()
    htype_matched.Draw("same")
    
    ctype.cd(2)
    htype_ratio.GetYaxis().SetRangeUser(0,1)
    htype_ratio.Draw()
    ctype.Update()
    ctype.SaveAs("truthxingpt_%s_truthtypes.png"%(name))    
    
    # DIST
    cdist = rt.TCanvas("cdist","Truth Dist", 1800, 1800)
    cdist.Draw()
    cdist.Divide(3,3)
    hdist = {}
    for n in range(0,7):
        cdist.cd(n+1)
        hdist[n] = rt.TH1D("hdist_type%d"%(n),"",50,0,200)
        tree.Draw("dist>>hdist_type%d"%(n),"truth_type==%d"%(n))
        hdist[7] = rt.TH1D("hdist_tot","",50,0,200)
        tree.Draw("dist>>hdist_tot")

    for n in range(0,8):
        cdist.cd(n+1)
        hdist[n].Draw()
        cdist.Update()
        cdist.SaveAs("truthxingpt_%s_dist.png"%(name))    

    # NPLANES
    # cplane = rt.TCanvas("cplane","Truth Plane", 1400, 400)
    # cplane.Draw()
    # cplane.Divide(2,1)
    # cplane.cd(1)
    # hplane_matched = rt.TH1D("hplane_matched", "",4,0,4)
    # hplane_totaled = rt.TH1D("hplane_totaled", "",4,0,4)
    # hplane_ratio   = rt.TH1D("hplane_ratio",   "",4,0,4)
    # tree.Draw("nplaneswcharge>>hplane_matched","matched==1")
    # tree.Draw("nplaneswcharge>>hplane_totaled","matched>=0")
    # tree.Draw("nplaneswcharge>>hplane_ratio","matched==1")
    # hplane_totaled.SetMinimum(0)
    # hplane_totaled.Draw()
    # hplane_matched.Draw("same")
    # hplane_matched.SetLineColor(rt.kRed)
    # hplane_ratio.Divide( hplane_totaled )

    # hplane_totaled.Draw()
    # hplane_matched.Draw("same")

    # cplane.cd(2)
    # hplane_ratio.GetYaxis().SetRangeUser(0,1)
    # hplane_ratio.Draw()
    # cplane.Update()

# RECO
rtree = rfile.Get( "recoxingptana" )

# RECO TYPE
crecotype = rt.TCanvas("crecotype","Reco End point type", 1400, 400)
crecotype.Draw()
crecotype.Divide(2,1)
crecotype.cd(1)
hrecotype_matched  = rt.TH1D("hrecotype_matched", "",7,0,7)
hrecotype_totaled  = rt.TH1D("hrecotype_totaled", "",7,0,7)
hrecotype_filtered = rt.TH1D("hrecotype_filtered", "",7,0,7)
hrecotype_filtpass = rt.TH1D("hrecotype_filtpass", "",7,0,7)
hrecotype_ratio    = rt.TH1D("hrecotype_ratio",   "",7,0,7)
hrecotype_ratiofil = rt.TH1D("hrecotype_ratiofil","",7,0,7)
rtree.Draw("type>>hrecotype_totaled", "matched>=0")
rtree.Draw("type>>hrecotype_filtered", "matched>=0 && passes==1")
rtree.Draw("type>>hrecotype_matched", "matched==1")
rtree.Draw("type>>hrecotype_filtpass", "matched==1 && passes==1")
rtree.Draw("type>>hrecotype_ratio",   "matched==1")
rtree.Draw("type>>hrecotype_ratiofil", "matched==1 && passes==1")
hrecotype_totaled.SetMinimum(0)
hrecotype_totaled.Draw()
hrecotype_matched.Draw("same")
hrecotype_matched.SetLineColor(rt.kRed)
hrecotype_filtered.SetLineColor(rt.kMagenta)
hrecotype_filtpass.SetLineColor(rt.kOrange+5)
hrecotype_ratiofil.SetLineColor(rt.kOrange+5)
hrecotype_ratio.Divide( hrecotype_totaled )
hrecotype_ratiofil.Divide( hrecotype_totaled )


hrecotype_totaled.Draw()
hrecotype_matched.Draw("same")
hrecotype_filtered.Draw("same")
hrecotype_filtpass.Draw("same")

crecotype.cd(2)
hrecotype_ratio.GetYaxis().SetRangeUser(0,1)
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
