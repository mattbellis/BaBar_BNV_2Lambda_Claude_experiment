#from ROOT.RooT import *
import ROOT

###############################################################
# Background and signal definition
###############################################################

################################################################################
################################################################################
def build_xy(var_ranges=[[5.2,5.3],[-0.2,0.2]]):
    #################################
    # Build two PDFs
    #################################
    xlo = var_ranges[0][0]
    xhi = var_ranges[0][1]
    ylo = var_ranges[1][0]
    yhi = var_ranges[1][1]

    x = ROOT.RooRealVar("x","m_{ES}",xlo, xhi, "GeV/c^{2}")
    y = ROOT.RooRealVar("y","#Delta E",ylo,yhi, "GeV")

    return x,y

################################################################################
################################################################################

################################################################################
# Crystal Barrel function: mES
################################################################################
def crystal_barrel_x(x):
    meanCB = ROOT.RooRealVar("meanCB","Gaussian #mu (CB) m_{ES}", 5.279)
    sigmaCB = ROOT.RooRealVar("sigmaCB"," Gaussian #sigma (CB) m_{ES}", 0.0028)
    alphaCB = ROOT.RooRealVar("alphaCB", "#alpha (CB) m_{ES}", 2.0)
    nCB = ROOT.RooRealVar("nCB","n of CB", 1.0)

    cb =     ROOT.RooCBShape("CB", "Crystal Barrel Shape PDF", x, meanCB, sigmaCB, alphaCB, nCB)

    pars = [meanCB,  sigmaCB, alphaCB, nCB]
    return pars, cb

################################################################################
# Crystal Barrel function: DeltaE
################################################################################
def crystal_barrel_y(y):
    meanCBdE = ROOT.RooRealVar("meanCBdE","Gaussian #mu (CB) #Delta E", 0.00)
    sigmaCBdE = ROOT.RooRealVar("sigmaCBdE","Gaussian #mu (CB) #Delta E", 0.020)
    alphaCBdE = ROOT.RooRealVar("alphaCBdE", "#alpha (CB) #Delta E", 2.0)
    nCBdE = ROOT.RooRealVar("nCBdE","n of CBdE", 1.0)

    # Second CB function: DeltaE
    sigmaCBdE_2 = ROOT.RooRealVar("sigmaCBdE_2","Gaussian #mu (CB) #Delta E (2)", 0.020)
    alphaCBdE_2 = ROOT.RooRealVar("alphaCBdE_2", "#alpha (CB) #Delta E (2)", 2.0)
    nCBdE_2 = ROOT.RooRealVar("nCBdE_2","n of CBdE_2", 1.0)

    cbdE =   ROOT.RooCBShape("CBdE", "Crystal Barrel Shape PDF: DeltaE", y, meanCBdE, sigmaCBdE, alphaCBdE, nCBdE)
    cbdE_2 = ROOT.RooCBShape("CBdE_2", "Crystal Barrel Shape PDF (2)", y, meanCBdE, sigmaCBdE_2, alphaCBdE_2, nCBdE_2)

    pars = [meanCBdE, sigmaCBdE, alphaCBdE, nCBdE, sigmaCBdE_2, alphaCBdE_2, nCBdE_2]
    return pars, cbdE, cbdE_2


################################################################################
# Double CB in dE
################################################################################
def double_cb_in_dE(cbdE, cbdE_2):

    ncbde1 = ROOT.RooRealVar("ncbde1","# cbde1 events,",500, 0, 1000000)
    ncbde2 = ROOT.RooRealVar("ncbde2","# cbde2 events",  50, 0, 1000000)
    double_cbdE = ROOT.RooAddPdf("double_cbdE","CBdE + CBdE_2",ROOT.RooArgList(cbdE, cbdE_2), ROOT.RooArgList(ncbde2))

    pars = [ncbde2]

    return pars, double_cbdE 

################################################################################
################################################################################
################################################################################
# Background PDF
################################################################################
################################################################################
################################################################################
# Linear in y (background)
################################################################################
def linear_in_y(y):
    p1 = ROOT.RooRealVar("poly1","Linear coefficient",-0.5) 
    rarglist = ROOT.RooArgList(p1)
    polyy = ROOT.RooPolynomial("polyy","Polynomial PDF",y, rarglist);

    pars = [p1]
    return pars, polyy 

################################################################################
# Argus background PDF
################################################################################
def argus_in_x(x):
    argpar = ROOT.RooRealVar("argpar","Argus shape par",-20.0)
    cutoff = ROOT.RooRealVar("cutoff","Argus cutoff",5.29)

    argus = ROOT.RooArgusBG("argus","Argus PDF",x,cutoff,argpar)

    pars = [argpar, cutoff]
    return pars, argus 
################################################################################

################################################################################
# Argus in NN
################################################################################
def argus_in_z(z):
    argpar_NN = ROOT.RooRealVar("argpar_NN","Argus shape par in NN",-7.0)
    cutoff_NN = ROOT.RooRealVar("cutoff_NN","Argus cutoff in NN",0.995)

    argpar_NN.setConstant(kFALSE)
    cutoff_NN.setConstant(kFALSE)

    argus_NN = ROOT.RooArgusBG("argus_NN","Argus NN PDF",z,cutoff_NN,argpar_NN)

    pars = [argpar_NN, cutoff_NN]
    return pars, argus_NN 
################################################################################

################################################################################
# BifurGaus in NN
################################################################################
def bifurgauss_in_z(z):
    mean_bfg = ROOT.RooRealVar("mean_bfg","Mean of bfg",0.975)
    sigma_bfg_L = ROOT.RooRealVar("sigma_bfg_L","Sigma L of bfg",0.50)
    sigma_bfg_R = ROOT.RooRealVar("sigma_bfg_R","Sigma R of bfg",0.01)

    mean_bfg.setConstant(kFALSE)
    sigma_bfg_L.setConstant(kFALSE)
    sigma_bfg_R.setConstant(kFALSE)

    bfg = ROOT.RooBifurGauss("bfg","BiFurGauss",z,mean_bfg,sigma_bfg_L,sigma_bfg_R)

    pars = [mean_bfg,sigma_bfg_L,sigma_bfg_R]
    return pars, bfg 
################################################################################

################################################################################
# Crystal Barrel function: NN (z)
################################################################################
def crystal_barrel_z(z):
    meanCB_NN = ROOT.RooRealVar("meanCB_NN","Gaussian #mu (CB) NN", 0.98)
    sigmaCB_NN = ROOT.RooRealVar("sigmaCB_NN"," Gaussian #sigma (CB) NN", 0.0028)
    alphaCB_NN = ROOT.RooRealVar("alphaCB_NN", "#alpha (CB) NN", 2.0)
    nCB_NN = ROOT.RooRealVar("nCB_NN","n of CB NN", 1.0)

    meanCB_NN.setConstant(kFALSE)
    sigmaCB_NN.setConstant(kFALSE)
    alphaCB_NN.setConstant(kFALSE)
    nCB_NN.setConstant(kTRUE)

    cb_NN =     ROOT.RooCBShape("CB_NN", "Crystal Barrel Shape PDF NN", z, meanCB_NN, sigmaCB_NN, alphaCB_NN, nCB_NN)

    pars = [meanCB_NN,  sigmaCB_NN, alphaCB_NN, nCB_NN]
    return pars, cb_NN


###############################################################
# ROOT.RooParametricStepFunction
###############################################################
def myRooKeys(z,data1):
    #############################

    kest1 = ROOT.RooKeysPdf("kest1","kest1",z,data1,ROOT.RooKeysPdf.NoMirror, 1.0)
    z.setBins(200, "cache")

    kc = ROOT.RooCachedPdf("kc","kc",kest1)

    # ROOT.RooDataHist
    rdh = kc.getCacheHist(ROOT.RooArgSet(z))

    # ROOT.RooHistPdf
    rhp = ROOT.RooHistPdf("nn_sig", "nn_sig", ROOT.RooArgSet(z), rdh)


    return [kest1, kc, rdh], rhp



################################################################################
#########################################
# Signal
#########################################
################################################################################
def sig_PDF(x, use_double_CB=False, workspace=None):
    pars = []

    cb = None
    cbdE = None
    cbdE_2 = None
    double_cbdE = None
    nn_sig = None
    sig_prod = None

    funcs = []

    funcs0 = []

    pars_0, cb = crystal_barrel_x(x)

    sig_prod = cb
    sig_prod.SetName("sig_pdf")

    pars += pars_0

    funcs = [cb, sig_prod]

    funcs += funcs0

    return pars, funcs, sig_prod


#########################################
# Background
#########################################
# Multiply the components
def bkg_PDF(x):

    nn_bkg = None
    bkg_prod = None
    argus = None
    polyy = None
    pars = []

    pars_a, argus = argus_in_x(x)
    pars_b = []

    bkg_prod = argus
    bkg_prod.SetName("bkg_pdf")

    pars += pars_a

    # Return all the sub-functions so that they stay active in the main program.
    funcs = [argus, bkg_prod]

    return pars, funcs, bkg_prod


#############################################################
#############################################################
def tot_PDF(x, use_double_CB=False, workspace=None):
    funcs = []

    pars_s, funcs_s, sig_pdf = sig_PDF(x, use_double_CB, workspace)
    pars_b, funcs_b, bkg_pdf = bkg_PDF(x)

    conv_factor_calc = ROOT.RooRealVar("conv_factor_calc","Conversion factor (calculated)",13.272) # Conversion factor, calculated
    conv_factor_fit  = ROOT.RooRealVar("conv_factor_fit", "Conversion factor (fit)",13.272) # Conversion factor, fit
    conv_factor_err  = ROOT.RooRealVar("conv_factor_err", "Error on conversion factor",0.32) # Conversion factor, fit

    branching_fraction = ROOT.RooRealVar("branching_fraction","Branching fraction",12.2)

    nbkg = ROOT.RooRealVar("nbkg","# bkg events,",150)
    nsig = ROOT.RooFormulaVar("nsig","conv_factor_fit*branching_fraction",ROOT.RooArgList(conv_factor_fit,branching_fraction))
    #nsig = ROOT.RooRealVar ("nsig","# sig events",150)

    # Gaussian constraint
    #gc = ROOT.RooFormulaVar("gc","exp(-(conv_factor_fit-conv_factor_calc)*(conv_factor_fit-conv_factor_calc)/(2.0*conv_factor_err*conv_factor_err))", \
    #ROOT.RooArgList(conv_factor_fit,conv_factor_calc,conv_factor_fit,conv_factor_calc,conv_factor_err,conv_factor_err))
    log_gc = ROOT.RooFormulaVar("log_gc","(conv_factor_calc-conv_factor_fit)*(conv_factor_calc-conv_factor_fit)/(2.0*conv_factor_err*conv_factor_err)", \
    ROOT.RooArgList(conv_factor_calc,conv_factor_fit,conv_factor_calc,conv_factor_fit,conv_factor_err,conv_factor_err))
            #ROOT.RooArgList(conv_factor_calc,conv_factor_fit,conv_factor_err))


    #sig_temp = ROOT.RooGenericPdf("sig_temp","gc*sig_pdf", ROOT.RooArgList(gc,sig_pdf))
    #bkg_temp = ROOT.RooGenericPdf("bkg_temp","gc*bkg_pdf", ROOT.RooArgList(gc,bkg_pdf))

    #total = ROOT.RooAddPdf("total","sig_temp + bkg_temp", ROOT.RooArgList(sig_temp, bkg_temp), ROOT.RooArgList(nsig, nbkg))
    #sub_total = ROOT.RooAddPdf("sub_total","sig_pdf + bkg_pdf", ROOT.RooArgList(sig_pdf, bkg_pdf), ROOT.RooArgList(nsig, nbkg))
    #sub_total = ROOT.RooAddPdf("sub_total","sig_pdf + bkg_pdf", ROOT.RooArgList(sig_pdf, bkg_pdf), ROOT.RooArgList(nsig, nbkg))

    #total = ROOT.RooGenericPdf("total","gc*sub_total", ROOT.RooArgList(gc,sub_total))

    total = ROOT.RooAddPdf("total","sig_pdf + bkg_pdf", ROOT.RooArgList(sig_pdf, bkg_pdf), ROOT.RooArgList(nsig, nbkg))

    pars = [nbkg, nsig, log_gc, conv_factor_calc, conv_factor_fit, conv_factor_err, branching_fraction]
    pars += pars_s
    pars += pars_b

    # Return all the sub-functions so that they stay active in the main program.
    funcs = [sig_pdf, bkg_pdf]
    funcs += funcs_s
    funcs += funcs_b

    return pars, funcs, total
#############################################################
#############################################################

#############################################################
#############################################################
