region_definitions = {}
region_definitions['Area'] = ['Fitting area', 'Signal area']
# Low and then high
region_definitions['signal MES'] = [5.27, 5.3]
region_definitions['fitting MES'] = [5.2, 5.3]

region_definitions['signal DeltaE'] = [-0.07, 0.07]
region_definitions['fitting DeltaE'] = [-0.2, 0.2]

region_definitions['sideband 1 DeltaE'] = [0.07, 0.14]
region_definitions['sideband 2 DeltaE'] = [-0.14, -0.07]
region_definitions['sideband MES'] = [5.27, 5.3]

region_definitions['inference'] = []
region_definitions['inference'].append([5.27, 5.3, -0.07, 0.07 ])
region_definitions['inference'].append([5.27, 5.3, -0.2, -0.07])
region_definitions['inference'].append([5.27, 5.3,  0.07, 0.20])
region_definitions['inference'].append([5.25, 5.27, -0.2, 0.0])
region_definitions['inference'].append([5.25, 5.27,  0.0, 0.20])
region_definitions['inference'].append([5.23, 5.25, -0.2, 0.0])
region_definitions['inference'].append([5.23, 5.25,  0.0, 0.20])
region_definitions['inference'].append([5.20, 5.23, -0.2, 0.0])
region_definitions['inference'].append([5.20, 5.23,  0.0, 0.20])

# For the Lambda0 identification
# Define the mass cuts around the lambda
lammass_world_average = 1.115683
width = 0.003 # GeV/c^2

lo = lammass_world_average - width
hi = lammass_world_average + width

region_definitions['Lambda0 mass'] = [lo, hi]
# Which variable is this for?
# FlightLen or postFitFlightLen
#region_definitions['Lambda0 flightlen'] = 1.00
# Lambda0postFitFlightSignificance
region_definitions['Lambda0 flightlen'] = 25.0



# Histogram definitions
hist_defs = {}

### Entries followed by lots of hashtags have outliers
hist_defs['BSphr'] =                           {"nbins":100, "lo":0,     "hi":0.2,   "label":"Sphericity"}
hist_defs["BpostFitDeltaE"]=                   {"nbins":100, "lo":-1.0,  "hi":1.0,   "label":"$\Delta$ E"}
hist_defs["BpostFitMes"]=                      {"nbins":100, "lo":5.2,   "hi":5.3,   "label":"$M_{ES}$ [GeV/c$^2$]"}
hist_defs["BThrust"]=                          {"nbins":100, "lo":0.9,   "hi":1.05,  "label":"B thrust"}
hist_defs["BCosThetaS"]=                       {"nbins":150, "lo":-1.05, "hi":1.05,  "label":"B Cos Theta S"}
hist_defs["Lambda0_unc_Mass"]=                 {"nbins":100, "lo":1.105, "hi":1.125, "label":"$\Lambda^0$ mass [GeV/c$^2$]"}
hist_defs["nTracks"]=                          {"nbins":100, "lo":0,     "hi":18,    "label":"number of tracks"}
hist_defs["nGoodTrkLoose"]=                    {"nbins":100, "lo":0,     "hi":14,    "label":"number of Good Tracks- Loose"}
hist_defs["nChargedTracks"]=                   {"nbins":100, "lo":0,     "hi":0.2,   "label":"number of Charged Tracks"}
hist_defs["R2"]=                               {"nbins":100, "lo":0,     "hi":1.05,  "label":"R2"}
hist_defs["R2All"]=                            {"nbins":100, "lo":0,     "hi":1,     "label":"R2All"}
hist_defs["thrustMag"]=                        {"nbins":100, "lo":.6,    "hi":1,     "label":"thrustMag"}
hist_defs["thrustMagAll"]=                     {"nbins":100, "lo":.6,    "hi":1,     "label":"thrustMagAll"}
hist_defs["thrustCosTh"]=                      {"nbins":100, "lo":0,     "hi":1,     "label":"thrustCosTh"}
hist_defs["thrustCosThAll"]=                   {"nbins":100, "lo":0,     "hi":1,     "label":"thrustCosThAll"}
hist_defs["sphericityAll"]=                    {"nbins":100, "lo":0,     "hi":0.75,  "label":"Sphericity"}
hist_defs["BCosSphr"]=                         {"nbins":100, "lo":-.8,   "hi":1,     "label":"BCosSphr"}
hist_defs["BCosThetaT"]=                       {"nbins":100, "lo":-1,    "hi":1,     "label":"BCosThetaT"}
hist_defs["BCosThrust"]=                       {"nbins":100, "lo":0,     "hi":1,     "label":"BCosThrust"}
hist_defs["BLegendreP2"]=                      {"nbins":100, "lo":0,     "hi":7,     "label":"BLegendreP2"}#################
hist_defs["BR2ROE"]=                           {"nbins":100, "lo":0,     "hi":1,     "label":"BR2ROE"}
hist_defs["BSphrROE"]=                         {"nbins":100, "lo":0,     "hi":1,     "label":"BSphrROE"}
hist_defs["BThrustROE"]=                       {"nbins":100, "lo":0.5,   "hi":1,     "label":"BThrustROE"}
hist_defs["Lambda0postFitFlight"]=             {"nbins":100, "lo":-1,    "hi":40,    "label":"Lambda0postFitFlight"} ######################
hist_defs["Lambda0postFitFlightSignificance"]= {"nbins":100, "lo":-20,   "hi":300,   "label":"Lambda0postFitFlightSignificance"} ##################
hist_defs["nTRK"]=                             {"nbins":20, "lo":0,      "hi":20,    "label":"# of charged tracks"}

hist_defs["BtagSideMes"]=                      {"nbins":100, "lo":4,     "hi":5,     "label":"tag-B $M_{ES}$ [GeV/c$^2$]"}


brfr_info= {}

brfr_info["p_lambda_conserved"]= {"branching_fraction":2.4e-7,"uncertainty_plus":1,"uncertainty_minus":.9}
