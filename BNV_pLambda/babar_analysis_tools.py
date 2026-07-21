import uproot
import awkward as ak

import pandas as pd

import matplotlib.pylab as plt
import numpy as np

import hist
from hist import Hist

import time

import myPIDselector

import math

import os

################################################################################
def load_datasets(topdir=None, sp_file_tag='Background_and_signal_SP_modes', \
                  collision_file_tag='Data', BNC=False, \
                  subset='Run1', sp_or_data=None, UNBLINDED=False):

    BLINDED_STRING = 'BLINDED'
    if UNBLINDED is True:
        BLINDED_STRING = 'UNBLINDED'


    subset_tag = ""
    if subset=='Run1' or subset=='run1':
        subset_tag = 'Only_Run_1'
    elif subset=='all' or subset=='All':
        subset_tag = 'All_runs'

    BNC_tag = ""
    if BNC is True:
        BNC_tag = "_BNC"

    data_sp, data_collision = None, None

    if sp_or_data=='sp' or sp_or_data=='SP' or sp_or_data==None:
        start= time.time()
        filename= f"{topdir}/{sp_file_tag}{BNC_tag}_{subset_tag}.parquet"
        print(f"Opening {filename}...")
        data_sp = ak.from_parquet(filename)
        print(f"Took {time.time()-start:.3f} seconds\n")


    if sp_or_data=='col' or sp_or_data=='collision' or sp_or_data==None:
        start= time.time()
        filename= f"{topdir}/{collision_file_tag}{BNC_tag}_{subset_tag}_{BLINDED_STRING}.parquet"
        if BNC is True:
            filename= f"{topdir}/{collision_file_tag}{BNC_tag}_{subset_tag}.parquet"
        print(f"Opening {filename}...")
        data_collision = ak.from_parquet(filename)
        print(f"Took {time.time()-start:.3f} seconds\n")

    return data_sp, data_collision 
################################################################################

################################################################################
def read_in_dataset_statistics(infilename='dataset_statistics.csv'):
    df = pd.read_csv(infilename)

    return df
################################################################################

################################################################################
def get_SP_cross_sections_and_labels(infilename='SP_cross_sections_and_labels.csv'):
    df = pd.read_csv(infilename)

    return df
################################################################################

################################################################################
def create_table_of_data_skims_statistics():
    df = read_in_dataset_statistics()
    
    dfspinfo = get_SP_cross_sections_and_labels()

    mask = dfspinfo['SP Mode']==1235
    bbbar_xsec = dfspinfo[mask]['Cross section [nb]'].values[0]

    mask = dfspinfo['SP Mode']==1237
    bbbar_xsec += dfspinfo[mask]['Cross section [nb]'].values[0]

    print(f"The total BBbar cross section is {bbbar_xsec} nb")

    mask = df['Data or MC']=='Data'
    #df[mask]

    mask = (df['Data or MC']=='Data') & (df['Skim']=='LambdaVeryVeryLoose')
    dftmp = df[mask][['Run', 'Luminosity (Data only) 1/pb','# of events (Data or MC)', '# of events (Data or MC) NOT SURE WHICH NUMBER TO USE']]

    dftmp['# of BBbar pairs'] = dftmp['Luminosity (Data only) 1/pb']*bbbar_xsec*1000

    dftmp['Run'] = dftmp['Run'].astype(int).astype(str)
    dftmp.loc['Total'] = dftmp.sum(numeric_only=True)

    dftmp.at['Total','Run'] = 'Total'

    header = []
    header.append('Run')
    header.append('Luminosity (1/pb)')
    header.append('\# skimmed events')
    header.append('\# org. events')
    header.append('\# BB pairs')

    caption = "Details of the numbers of events and luminosity from the {\\tt LambdaVeryVeryLoose} skim used in this analysis."
    label = 'tab:dataskims'

    df.style.to_latex(position_float='centering')

    output = dftmp.to_latex(index=False, header=header, float_format="%.1f", caption=caption, label=label, position='h')

    # Add in centering by replacing the first EOL with "EOL + \centering + EOL"
    output = output.replace('\n','\n\centering\n', 1)

    # Add an hline learn the bottom above the total
    output = output.replace('Total','\hline\nTotal', 1)

    return output


################################################################################

###################################################################################
def calculate_bits_for_PID_selector(trkidx, trk_selector_map, verbose=0):
    
    bits = None

    # If there is no trk index passed in, just calculate the bits for
    # all of the tracks
    if trkidx is None:
        trkidx = ak.local_index(trk_selector_map)

    # Grab the tracks that map on to the particle/collection we are interested in 
    subset_of_trk_selector_map = trk_selector_map[trkidx]
    if verbose:
        print("values in the subset of the trk selector map")
        print(subset_of_trk_selector_map)
        print()
        
    # We need this to convert the numbers in the selector to binary
    binary_repr_vec = np.vectorize(np.binary_repr)

    # Grab the number of entries in each so we can unflatten this later
    counts = ak.num(subset_of_trk_selector_map)
    
    # Now get the binary representation (as a string) for the flattened subset
    binrep = binary_repr_vec(ak.flatten(subset_of_trk_selector_map), width=32)

    if verbose:
        print("binary representation of selector map")
        print(binrep)
        print()

    # Convert the string to integers
    tempbits = np.array(binrep).astype(int)
    bits = ak.unflatten(tempbits,counts)

    if verbose:
        print("flattened integer representation of selectors as binary (int)")
        print(tempbits)
        print()
        print("unflattened integer representation of selectors as binary (int)")
        print(bits)
        print()

    return bits

##########################################################################

##########################################################################

def mask_PID_selection(bits, selector, pid_map_object):

    bit_to_look_for = pid_map_object.selectors.index(selector)

    place = int(math.pow(10,bit_to_look_for))
    #print(place)

    mask = bits // place % 10

    mask_bool = ak.values_astype(mask,bool)

    return mask_bool


##########################################################################
def scaling_value(spmode, dataset_information=None, cs_data=None, plot= False, verbose= False):
    mode= spmode 

    if spmode==0 or spmode==-999:
        return 1 # This is data
    
    mc_mask= (dataset_information["SP mode"]== mode) & (dataset_information["Skim"] != "LambdaVeryVeryLoose")
    # Pulls out all unskimmed MC with desired SP mode
    
    nevents_mc= dataset_information[mc_mask]["# of events (Data or MC)"].sum()
    # Sums the number of unskimmed MC events of the desired SP mode
                
    cs_mask= cs_data["SP Mode"]== mode

    cs= cs_data[cs_mask]["Cross section [nb]"]

    cs= cs.values[0]

    mask = (dataset_information['Data or MC'] == 'Data') & (dataset_information['Skim'] != 'LambdaVeryVeryLoose')
    int_lumi = dataset_information[mask]['Luminosity (Data only) 1/pb'].sum()
    
    n_exp_in_data= cs* int_lumi*1000 
    # Factor of 1000 arises from multiplying nanobarnes by 1/picobarnes
    
    scaling= n_exp_in_data/nevents_mc
    # MC is intentionally overgenerated, so to figure out how to weight the MC, divide the number generated by the total number of MC events

    
    if plot== True:
        spmask = data['spmode']== str(mode)
        x = data[spmask]['BpostFitMes'][:,0]
        plt.hist(x,bins=100, range=(5., 5.3), weights=scaling*np.ones(len(x)));
        plt.title(f"Scaling value for SP-{mode}: {scaling:.4f}")

    if verbose== True: 
        print(f"- Cross section for this SP mode is    {cs} nb")
        print(f"- # of events generated for SP-{mode}: {nevents_mc:13d}")
        print(f"- Number expected in data:             {n_exp_in_data:.1f}")
        print(f"- Integrated Luminosity:               {int_lumi:.1f} 1/pb")
        print("The scaling value for this SP mode is: ")

    return scaling

##########################################################################

def table_from_df(df, outfilename, caption="DEFAULT", label='tab:DEFAULT'):
    output = df.to_latex(index=False,
                  float_format="{:.4f}".format,
    )  # converts dataframe into latex readable text
    full_table = "\\begin{table}[h]\n" # initializes the table before the beginning of the tabular 
    full_table += "\\caption{" + f"{caption}" +"}\n" 
    full_table += "\\label{" + f"{label}" +"}\n" 
    full_table += "\\centering\n" 
    full_table += output #includes the converted dataframe in the table
    full_table += "\\end{table}" # ends the table, same purpose as begin{table} 
    filename= f"table_{outfilename}.tex"
    path= f"tables/{filename}"
    outfile = open(path, 'w')
    outfile.write(full_table)
    outfile.close()    
    return print(full_table) #make sure to return the print() of the full_table, otherwise it'll be one big string that latex can't handle


##########################################################################

def indices_to_booleans(indices, array_to_slice):
    whole_set, in_set = ak.unzip(ak.cartesian([ 
                                  ak.local_index(array_to_slice), indices ], nested=True))

    return ak.any(whole_set == in_set, axis=-1)

##########################################################################

# Only works when there is one candidate for the B and the Lambda0 right now
def build_antiproton_antimask(data, pps, selector = 'SuperLooseKMProtonSelection', IS_MC=True, verbose=0):

    if verbose:
        print("The MC and tracks for the first entry")
        idx = data['TRKMCIdx'][0]
        if IS_MC:
            mclund = data['mcLund'][0]

            for i,id in enumerate(idx):
                print(f"{i:2d}  {id:4d}   {mclund[id]}")
            print()

    
    lamd1idx = data['Lambda0d1Idx']
    lamd1Lund = data['Lambda0d1Lund']

    if verbose:
        print(f'lamd1idx\n{lamd1idx}')
        print(f'lamd1Lund\n{lamd1Lund}')
    
    d2idx = data['Bd2Idx']
    d2Lund = data['Bd2Lund']

    if verbose:
        print()
        print(f'B d2idx\n{d2idx}')
        print(f'B d2Lund\n{d2Lund}')
        print()

    qBd2 = (data['Bd2Lund'])/np.abs(data['Bd2Lund'])
    #print(qBd2)
    #print(qBd2[:,0])
    # Assuming we have only one B candidate here
    qBd2 = qBd2[:,0]

    qlamd1 = (data['Lambda0d1Lund'])/np.abs(data['Lambda0d1Lund'])
    qlamd1 = qlamd1[:,0]

    if verbose:
        print("Proton charges ----------------")
        print(f"qBd2\n{qBd2}\n")
        print(f"qlamd1\n{qlamd1}\n")
        print("Proton charges ----------------\n\n")

    qtrk = (data['TRKLund'])/np.abs(data['TRKLund'])

    if verbose:
        print("Track charges ----------------")
        print(f"qtrk\n{qtrk}\n")
        print("Track charges ----------------\n\n")
    
    
    trkidx_proton = data['pTrkIdx']

    if verbose:
        print(f"# of protons: {data['np']}")
        print(f"trkidx_proton (the track index for labeled protons) \n{trkidx_proton}")
        print()
    
    lamd1_trkidx = trkidx_proton[lamd1idx]
    d2_trkidx = trkidx_proton[d2idx]

    if verbose:
        print(f"lamd1_trkidx\n{lamd1_trkidx}\n")    
        print(f"B d2_trkidx\n{d2_trkidx}\n")
    

    trk_selector_map = data['pSelectorsMap']
    
    if verbose:
        print(f"qtrk\n{qtrk}\n")
        print(f"proton trk_selector_map\n{trk_selector_map}\n")
        print(f"pion   trk_selector_map\n{data['piSelectorsMap']}\n")
        print(f"kaon   trk_selector_map\n{data['KSelectorsMap']}\n")

    lamproton_selector_map = trk_selector_map[lamd1_trkidx]

    if verbose:
        print(lamproton_selector_map)
        print()
    
    bool_idx1 = indices_to_booleans(lamd1_trkidx, trk_selector_map)
    bool_idx2 = indices_to_booleans(d2_trkidx, trk_selector_map)

    if verbose:
        print(f"proton track selectors for the index of the proton from the lambda (and not that)")
        print(f"The boolean mask bool_idx1\n{bool_idx1}")
        print(trk_selector_map[bool_idx1])
        print(trk_selector_map[~bool_idx1])
        print()

        print(f"proton track selectors for the index of the proton from the B (and not that)")
        print(f"The boolean mask bool_idx1\n{bool_idx2}")
        print(trk_selector_map[bool_idx2])
        print(trk_selector_map[~bool_idx2])
        print()
        
    
        # Both protons
        print(f"Boolean for both protons (or) \n{(bool_idx1 | bool_idx2)}")
        print(f"trk selectors for both protons (or) \n{trk_selector_map[bool_idx1 | bool_idx2]}")
        print(f"trk selectors for other protons     \n{trk_selector_map[~(bool_idx1 | bool_idx2)]}")
        print()
    
    
    
    pbits = calculate_bits_for_PID_selector(None, trk_selector_map)

    if verbose:
        print(f"The pbits -----------------")
        print(f'All the pbits\n{pbits}')
        print(f"bool for protons from Lambda\n{bool_idx1}")
        print(f'pbits for protons from Lambda     pbits[bool_idx1]\n{pbits[bool_idx1]}')
        print(f'pbits for protons not from Lambda pbits[~bool_idx1]\n{pbits[~bool_idx1]}')
        print()

        print(f"bool for protons from B\n{bool_idx2}")
        print(f'pbits for protons from B     pbits[bool_idx2]\n{pbits[bool_idx2]}')
        print(f'pbits for protons not from B pbits[~bool_idx2]\n{pbits[~bool_idx2]}')
        print()

        
    mask_bool = mask_PID_selection(pbits, selector, pps)

    if verbose:
        print(f"mask_bool, (mask_PID_selection) the tracks that pass {selector}")
        print(mask_bool)
        print(pbits[mask_bool])
        print(pbits[~mask_bool])
        print()
        
        print(f"pbits that are not the final state lambda but are protons")
        print(pbits[~bool_idx1 & mask_bool])
        print()
        
        print(f"TRKcharge that are not the final state lambda but are protons")
        print(qtrk[~bool_idx1 & mask_bool])
        print()
    
    charge_test_idx1 = qtrk[~bool_idx1 & mask_bool] == -qlamd1
    charge_test_idx2 = qtrk[~bool_idx2 & mask_bool] == -qBd2

    if verbose:
        print(f"charge test idx1\n{charge_test_idx1}")
        print(f"charge test idx2\n{charge_test_idx2}")
    

    # This is True if there is an antiproton in the other final state particles
    # This is likely (I think) with background and less likely with signal
    charge_test = qtrk[~(bool_idx1 | bool_idx2) & mask_bool] == -qlamd1

    if verbose:
        print(f"charge test both: is neither Lam/B dau, passes selectors, opposite charge\n{charge_test}")

    if verbose:
        print("\nCounting..............")
        print(f'{charge_test[charge_test] = }')
        print(f'{ak.num(charge_test[charge_test]) = }')
        print(f'{ak.sum(ak.num(charge_test[charge_test])) = }')
        print()
    
    #nhave_opposite = ak.sum(ak.num(charge_test[charge_test]))
    #n = ak.num(charge_test, axis=0)
    
    have_opp = ak.num(charge_test[charge_test])
    
    nhave_opp = len(have_opp[have_opp>0])

    ndont_have_opp = len(have_opp[have_opp==0])

    n = len(data)

    if verbose:
        print(f"# of events:                                  {n}")
        print(f"# of events that don't have opposite protons: {ndont_have_opp}")
        print(f"# of events that       have opposite protons: {nhave_opp}")

    # Select events that *don't* have opposite protons
    # We need the opposite of this when we return it
    mask_no_antiprotons = ak.any(charge_test, axis=-1)

    return ~mask_no_antiprotons, charge_test

################################################################################
def create_empty_histograms(hist_defs): 
    ### Creates empty Hist object histograms based on the information in the dictionary above
    ### Then overwrites the information in the dictionary to be the hist object. All_hists goes from containing unconnected
    ### info about each variable to a single object containing all the same info 
    
    all_hists={}
    for var in hist_defs.keys():
        h = Hist.new.Reg(hist_defs[var]["nbins"], hist_defs[var]["lo"], hist_defs[var]["hi"], name='var', label=f"{hist_defs[var]['label']}") \
                 .StrCat([], name="SP", label="SP modes", growth=True)\
                 .StrCat([], name="cuts", label="Cuts", growth=True)\
                 .Weight()
    
        all_hists[var] = h

    return all_hists



################################################################################
def fill_histograms_v3(data, subset=None, empty_hists=None, spmodes=['998'], weights=[1.0], cuts=None):
    ### Takes the dictionary of objects we made before and fills them 
    ### with the correct information, based on SP mode and Cut. 
    ### Each cut pares down the background and hopefully makes the signal more apparent

    if empty_hists is None:
        print("Must pass in histograms to be filled")
        return None

    if cuts is None:
        print("Must pass in lists of cuts")
        return None

    if subset is None:
        subset = list(empty_hists.keys())

    # Save the result of the cuts in a dataframe
    df_dict = {}
    df_dict['var'] = []
    df_dict['cut'] = []
    df_dict['spmode'] = []
    df_dict['n'] = []

    for key in subset: 
        print(key)

        # Check to see if the variable is in the array
        if key not in list(data.fields):
            continue

        for spmode in spmodes:
            #print(spmode)
            weight = 1
            '''
            if spmode=='-999':
                weight = .005
            else:
                weight = weights[spmode]
            '''
            # Assuming this is passed in OK
            weight = weights[spmode]

            for cutname in cuts.keys():
                #print(f"filling with {cutname}  for {spmode}   for {key}")
                cut = cuts[cutname]['event']
                mask_sp = data['spmode']==spmode
                n = -1
                # Apply the cuts and fill the histograms
                #if key[0]=='B' or key.find('Lambda0')==0:
                #    x = ak.flatten(data[mask_ev][key][mask_part[mask_ev]])
                #else:
                #     x = data[mask_ev][key]
                #x = ak.flatten(data[key])
                x = data[key][cut & mask_sp]
                #print(type(x), type(x[0]))
                #if type(x[0])==int or type(x[0])==float:
                # Make sure array is not 0-length
                if len(x)<=0:
                    continue

                try:
                    float(x[0])
                        #1#x = ak.flatten(data[key][cut])
                except:
                    x = ak.flatten(x)

                n = len(x)
                if n>0:
                    #print(f"cuts = {cutname}")
                    empty_hists[key].fill(var=x, SP= spmode, cuts= f"{cutname}", weight= weight)

                # Fill the dataframe dictionary
                df_dict['var'].append(key)
                df_dict['cut'].append(cutname)
                df_dict['spmode'].append(spmode)
                df_dict['n'].append(n)

    df = pd.DataFrame.from_dict(df_dict)

    return df

################################################################################
def fill_histograms_v2(ak_arr, empty_hists, spmodes=['998'], weights=[1.0]):
    ### Takes the dictionary of objects we made before and fills them 
    ### with the correct information, based on SP mode and Cut. 
    ### Each cut pares down the background and hopefully makes the signal more apparent

    # Save the result of the cuts in a dataframe
    df_dict = {}
    df_dict['var'] = []
    df_dict['cut'] = []
    df_dict['spmode'] = []
    df_dict['n'] = []

    for key in empty_hists.keys(): 
        print(key)

        for spmode in spmodes:
            #print(spmode)
            weight = 1
            if spmode=='-999':
                weight = .005
            else:
                weight = weights[spmode]
            
            all_cuts = {}
            mask_sp= (ak_arr.spmode== spmode)            
            mask_fl0 = (ak_arr['Lambda0FlightLen']>=0)
            mask_fl1 = (ak_arr['Lambda0FlightLen']>=1)
            mask_ntrk = (ak_arr['nTRK'] > 4)

            # Blinding
            mes = ak_arr['BpostFitMes']
            de  = ak_arr['BpostFitDeltaE']           
            # Stuff
            blinding_mask=   (mes>5.27) & ((de>-.07) & (de<.07)) 
            fitarea_mask=  (mes>5.2) & ((de>-.2) & (de<.2)) 

            #var_mask = mask_fl & ~blinding_mask & fitarea_mask
            
            all_cuts['0'] = {'mask_ev':mask_sp, 'mask_part':mask_fl0}
            all_cuts['1'] = {'mask_ev':mask_sp & mask_ntrk, 'mask_part':mask_fl0}
            all_cuts['2'] = {'mask_ev':mask_sp & mask_ntrk, 'mask_part':mask_fl1}
            all_cuts['3'] = {'mask_ev':mask_sp & mask_ntrk, 'mask_part':mask_fl1 & ~blinding_mask & fitarea_mask}
            all_cuts['4'] = {'mask_ev':mask_sp & mask_ntrk, 'mask_part':mask_fl1 & fitarea_mask}

            for cutname,cutmasks in all_cuts.items():

                #print(f"cutname: {cutname}")

                mask_ev =   cutmasks['mask_ev']
                mask_part = cutmasks['mask_part']

                n = -1
                # Apply the cuts and fill the histograms
                if key[0]=='B' or key.find('Lambda0')==0:
                    x = ak.flatten(ak_arr[mask_ev][key][mask_part[mask_ev]])
                else:
                    x = ak_arr[mask_ev][key]

                n = len(x)
                empty_hists[key].fill(var=x, SP= spmode, cuts= f"{cutname}", weight= weight)

                # Fill the dataframe dictionary
                df_dict['var'].append(key)
                df_dict['cut'].append(cutname)
                df_dict['spmode'].append(spmode)
                df_dict['n'].append(n)

    df = pd.DataFrame.from_dict(df_dict)
    
    return df
    
################################################################################
def plot_histograms(all_hists, vars=[], bkg_spmodes=['998'], datamodes=['0'], sig_spmodes=['-999'], sig_weight=1, cut='0', save= True, overlay_data=True, only_stacked=False, fixed_grid=None, path='BNV_pLambda_plots', extra_tag=""):
    
    ### makes a directory (if it doesn't already exist) for these plots.
    ### plots will be saved to this dictionary if save= true
    current_dir= os.getcwd()
    directory = "BNV_pLambda_plots"
    #path= os.path.join(current_dir,directory)
    #if os.path.isdir(path)== False:
    #    os.mkdir(path)
    
    if len(vars) == 0:
        vars = list(all_hists.keys())

    ### color scheme dictionary
    cd= {}
    cd["998"]=  {"tab:blue"}
    cd["1005"]= {"tab:orange"}
    cd["-999"]= {"tab:brown"}
    cd["1235"]= {"tab:green"}
    cd["1237"]= {"tab:red"}
    cd["3981"]= {"tab:purple"}
    cd["3429"]= {"tab:pink"}
    cd["991"]=  {"tab:blue"}
    cd["0"]= {"tab:cyan"}

    #print("here: ", bkg_spmodes)

    if only_stacked and fixed_grid:
        # For 4 rows and 5 columns
        #width = fixed_grid[0] * 5
        #height = fixed_grid[1] * 3

        # For 1 row x 3 columns
        width = fixed_grid[0] * 12
        height = fixed_grid[1] * 1.5

        # For 1 row x 3 columns, not as high
        #width = fixed_grid[0] * 12
        #height = fixed_grid[1] * 0.75


        plt.figure(figsize=(width,height))           

    
    for axes_idx,var in enumerate(vars):
        
        h = all_hists[var]
                
        if only_stacked and not fixed_grid:   # if you only want the stacked sp and not the breakdown for individual modes 
            plt.figure(figsize=(5,3))
        elif only_stacked and fixed_grid:     # fixed grid is a grid of subplots 
            plt.subplot(fixed_grid[0], fixed_grid[1], axes_idx+1)
        else:
            plt.figure(figsize=(18,12))
            plt.subplot(3,3,1)
            
        #print("there: ", var, bkg_spmodes, cut)
        h[:,bkg_spmodes,cut].stack('SP')[:].project('var').plot(stack=True, histtype="fill")
        h[:,sig_spmodes,cut].project('var').plot(histtype="step", color='yellow', label= "signal")

        if overlay_data:
            h[:,datamodes,cut].project('var').plot(histtype="errorbar", color='black', label='Data')

        plt.legend(loc='upper left')
        plt.xlabel(plt.gca().get_xlabel(), fontsize=18)


        # If we are only plotting the stacked histograms, then we can skip over plotting them individually
        if not only_stacked:
            
            # Plot all the others 
            all_modes = bkg_spmodes + sig_spmodes + datamodes
            
            for idx,spmode in enumerate(all_modes):
                plt.subplot(3,3,idx+2)
                h[:,spmode,cut].project('var').plot(histtype="fill", label=spmode, color= cd[str(spmode)])
                plt.legend()
                plt.xlabel(plt.gca().get_xlabel(), fontsize=18)

        plt.tight_layout()
        
        '''
        if save== True:
            
            outfilename=f"plot_hist_cut{cut}_{var}.png" 
            if only_stacked and not fixed_grid:
                outfilename=f"plot_hist_cut{cut}_ONLY_STACKED_{var}{extra_tag}.png" 
                
            plt.savefig(f"{path}/{outfilename}")
        '''

    if save== True and fixed_grid and only_stacked:
        # name of .png saved to computer based on fields specified on function call 
        varnames = "_".join(vars)
        outfilename=f"plot_hist_cut{cut}_ONLY_STACKED_FIXED_GRID_{varnames}{extra_tag}.png" 
            
        plt.savefig(f"{path}/{outfilename}")
################################################################################

################################################################################
def get_signal_mask(data, region_definitions):
    mes= data["BpostFitMes"]
    de= data["BpostFitDeltaE"]
    signal_mask = (mes > region_definitions["signal MES"][0]) & ((de>region_definitions["signal DeltaE"][0]) & (de<region_definitions["signal DeltaE"][1])) 
    #fitarea_mask   =  (mes>5.2) & ((de>-.2) & (de<.2))
    return signal_mask 

################################################################################

################################################################################
def get_fit_mask(data, region_definitions):
    mes= data["BpostFitMes"]
    de= data["BpostFitDeltaE"]
    fit_mask = (mes > region_definitions["fitting MES"][0]) & ((de>region_definitions["fitting DeltaE"][0]) & (de<region_definitions["fitting DeltaE"][1])) 
    return fit_mask 

################################################################################
################################################################################
def get_lambda0_mask(data, region_definitions, flightlenvar='Lambda0FlightLen'):
    
    # Which variable are we using to cut upon
    cutvariable = data[flightlenvar]

    # Get the cut value from our region_definitions file
    # Cut on the flight length
    mask_fl = cutvariable>region_definitions['Lambda0 flightlen']

    # Cut on the mass
    lo = region_definitions['Lambda0 mass'][0]
    hi = region_definitions['Lambda0 mass'][1]

    m = data['Lambda0_unc_Mass']
    mask_lambda0 = (m>lo) & (m<hi) & mask_fl

    # We also only want to keep events with one B candidate and one Lambda0 candidate, 
    # so let's do that here
    nlambda0 = ak.num(m[mask_lambda0])
    mB = data['BpostFitMes']
    nB = ak.num(mB)
    #print(nlambda0)

    mask_event_nlambda0_and_nB = (nlambda0==1) & (nB==1)

    return mask_lambda0, mask_event_nlambda0_and_nB

################################################################################


################################################################################
def get_duplicates_mask(data, final_state='pLambda'):

    # Keep events with only 1 B candidate
    nB = ak.num(data['BMass'])
    nlambda0 = ak.num(data['Lambda0_unc_Mass'])
    
    mask_event_nlambda0_and_nB = None

    if final_state.lower()=='plambda' or  \
       final_state.lower()=='plam' or  \
       final_state.lower()=='plam0':
        mask_event_nlambda0_and_nB = (nlambda0==1) & (nB==1)
    elif final_state.lower()=='lam0lam0' or  \
       final_state.lower()=='lambdalambda' or  \
       final_state.lower()=='lamlam':
        mask_event_nlambda0_and_nB = (nlambda0==2) & (nB==1)

    return mask_event_nlambda0_and_nB

################################################################################

def PID_masks(data, \
              lamp_selector='VeryLooseKMProtonSelection', \
              lampi_selector='SuperLooseKMPionMicroSelection', \
              Bp_selector='VeryLooseKMProtonSelection', \
             verbosity=0):
    # Get these maps first
    pps = myPIDselector.PIDselector("p")
    pips = myPIDselector.PIDselector("pi")

    if verbosity:
        print("Names of selectors:\npions")
        print(pips.selectors)
        print("\nprotons\n")
        print(pps.selectors) 

    # Proton and pion information from the Lambda decay
    # These are the index of the proton (d1) and pion (d2) in those lists
    d1idx = data['Lambda0d1Idx']
    d2idx = data['Lambda0d2Idx']
    
    d1lund = data['Lambda0d1Lund']
    d2lund = data['Lambda0d2Lund']
    
    Bd2idx = data['Bd2Idx']
    Bd2lund = data['Bd2Lund']

    if verbosity==1:
        print(d1lund)
        print(d2lund)
        print(Bd2lund)
        print()
        
        print(d1idx)
        print(d2idx)
        print(Bd2idx)
        print()
    
    trkidx_proton = data['pTrkIdx']
    trk_selector_map_proton = data['pSelectorsMap']
    
    trkidx_pion = data['piTrkIdx']
    trk_selector_map_pion = data['piSelectorsMap']
    
    # Proton
    pbits = calculate_bits_for_PID_selector(trkidx_proton, trk_selector_map_proton, verbose=verbosity)
    # Pion
    pibits = calculate_bits_for_PID_selector(trkidx_pion, trk_selector_map_pion, verbose=verbosity)
    
    
    #selector_proton = 'TightKMProtonSelection'
    #selector_pion = 'TightKMPionMicroSelection'
    #print(f"Now trying to create a mask with {selector_proton}")
    #print(f"Now trying to create a mask with {selector_pion}")
    
    
    mask_bool_proton = mask_PID_selection(pbits[d1idx], lamp_selector, pps)
    mask_bool_protonB = mask_PID_selection(pbits[Bd2idx], Bp_selector, pps)
        
    mask_bool_pion = mask_PID_selection(pibits[d2idx], lampi_selector, pips)

    return mask_bool_proton, mask_bool_pion, mask_bool_protonB


##########################################################################
##########################################################################
def PID_masks_tiny_hydrogen(data, \
              p_selector='VeryLooseKMProtonSelection', \
              e_selector='SuperLooseKMElectronMicroSelection', \
             verbosity=0):
    # Get these maps first
    pps = myPIDselector.PIDselector("p")
    eps = myPIDselector.PIDselector("e")

    if verbosity:
        print("Names of selectors:\nelectrons")
        print(eps.selectors)
        print("\nprotons\n")
        print(pps.selectors) 

    # Proton and electron information from the neutron decay
    # These are the index of the proton (d1) and electron (d2) in those lists
    d1idx = data['nd1Idx']
    d2idx = data['nd2Idx']
    
    d1lund = data['nd1Lund']
    d2lund = data['nd2Lund']
    
    if verbosity==1:
        print(d1lund)
        print(d2lund)
        print(Bd2lund)
        print()
        
        print(d1idx)
        print(d2idx)
        print()
    
    trkidx_proton = data['pTrkIdx']
    trk_selector_map_proton = data['pSelectorsMap']
    
    trkidx_electron = data['eTrkIdx']
    trk_selector_map_electron = data['eSelectorsMap']
    
    # Proton
    pbits = calculate_bits_for_PID_selector(trkidx_proton, trk_selector_map_proton, verbose=verbosity)
    # electron
    ebits = calculate_bits_for_PID_selector(trkidx_electron, trk_selector_map_electron, verbose=verbosity)
    
    
    #selector_proton = 'TightKMProtonSelection'
    #selector_pion = 'TightKMPionMicroSelection'
    #print(f"Now trying to create a mask with {selector_proton}")
    #print(f"Now trying to create a mask with {selector_pion}")
    
    
    mask_bool_proton = mask_PID_selection(pbits[d1idx], p_selector, pps)
        
    mask_bool_electron = mask_PID_selection(ebits[d2idx], e_selector, eps)

    return mask_bool_proton, mask_bool_electron


##########################################################################


##########################################################################
def count_events_in_inference_regions(mes, DeltaE, region_definitions=None, tag="DEFAULT"):
    dict_temp = {'meslo':[], 'meshi':[], 'dElo':[], 'dEhi':[],'N':[]}
    for xlo,xhi,ylo,yhi in region_definitions['inference']:
        mask_temp = (mes>xlo) & (mes<xhi) & (DeltaE>ylo) & (DeltaE<yhi)
        N = len(mask_temp[mask_temp])
        print(f"{xlo:.2f} {xhi:.2f} {ylo:6.2f} {yhi:6.2f} {N:6d}")
        dict_temp['meslo'].append(xlo)
        dict_temp['meshi'].append(xhi)
        dict_temp['dElo'].append(ylo)
        dict_temp['dEhi'].append(yhi)
        dict_temp['N'].append(N)

    df = pd.DataFrame.from_dict(dict_temp)

    outfilename = f"inference_region_counts_{tag}.parquet"
    df.to_parquet(outfilename)


    return df,outfilename

##########################################################################



##########################################################################
def plot_mes_vs_DeltaE(mes, DeltaE, draw_signal_region=False, draw_sidebands=False, draw_inference_bins=False, tag=None, region_definitions=None, bins=100, zoom=False, plot_full=True, draw_fit_region=False):

    meslo = region_definitions['fitting MES'][0]
    meshi = region_definitions['fitting MES'][1]
    DeltaElo = -0.5
    DeltaEhi =  0.5

    sigmeslo = region_definitions['signal MES'][0]
    sigmeshi = region_definitions['signal MES'][1]
    sigDeltaElo = region_definitions['signal DeltaE'][0]
    sigDeltaEhi = region_definitions['signal DeltaE'][1]

    sbDE1lo = region_definitions['sideband 1 DeltaE'][0]
    sbDE1hi = region_definitions['sideband 1 DeltaE'][1]
    sbDE2lo = region_definitions['sideband 2 DeltaE'][0]
    sbDE2hi = region_definitions['sideband 2 DeltaE'][1]

    sbmeslo = region_definitions['sideband MES'][0]
    sbmeshi = region_definitions['sideband MES'][1]

    fitDeltaElo = region_definitions['fitting DeltaE'][0]
    fitDeltaEhi = region_definitions['fitting DeltaE'][1]

    if zoom==True and region_definitions is not None:
        DeltaElo = region_definitions['fitting DeltaE'][0]
        DeltaEhi = region_definitions['fitting DeltaE'][1]

    h= Hist(
        hist.axis.Regular(bins,meslo,meshi,name= "sig_BPFM", label= "M$_{ES}$ [GeV/c$^2$]", flow= True),
        hist.axis.Regular(bins,DeltaElo,DeltaEhi,name= "bkg_BPFMDE", label= "$\Delta$E [GeV]", flow= True),
    )

    # normal fill
    h.fill(mes, DeltaE)

    #plt.gca()

    if plot_full:
        h.plot2d_full(
                #main_cmap="coolwarm",
            main_cmap="plasma",
            top_ls="--",
            top_color="orange",
            top_lw=2,
            side_ls=":",
            side_lw=2,
            side_color="steelblue",
        )
    else:
        h.plot2d(
                #main_cmap="coolwarm",
            cmap="plasma",
        )

    #plt.xlim(5.1,5.3)
    #plt.ylim(-.5,.5)
    #plt.show()

    def draw_box(xlo,xhi, ylo,yhi,fmt='w--'):
        plt.plot([xlo,xlo,xhi,xhi,xlo],[ylo,yhi,yhi,ylo,ylo], fmt, linewidth= 4)
        mask_temp = (mes>xlo) & (mes<xhi) & (DeltaE>ylo) & (DeltaE<yhi)
        #print(mask_temp)
        print(ylo,yhi)
        print(f"{xlo:.2f} {xhi:.2f} {ylo:6.2f} {yhi:6.2f} {len(mask_temp[mask_temp]):6d}")



    #'''
    if draw_signal_region==True:
        print("Signal region")
        #plt.plot([meslo,meslo,meshi,meshi,meslo],[DeltaElo,DeltaEhi,DeltaEhi,DeltaElo,DeltaElo], "w-", linewidth= 4)
        #plt.plot([sigmeslo,sigmeslo,sigmeshi,sigmeshi,sigmeslo],[sigDeltaElo,sigDeltaEhi,sigDeltaEhi,sigDeltaElo,sigDeltaElo], "w--", linewidth= 4)
        draw_box(sigmeslo,sigmeshi,sigDeltaElo,sigDeltaEhi)
    #'''
    if draw_sidebands==True:
        print("Sidebands")
        draw_box(sbmeslo,sbmeshi,sbDE1lo,sbDE1hi)
        draw_box(sbmeslo,sbmeshi,sbDE2lo,sbDE2hi)

    if draw_fit_region==True:
        print("Fitting region")
        plt.plot([meslo,meshi],[fitDeltaElo, fitDeltaElo], 'w--')
        plt.plot([meslo,meshi],[fitDeltaEhi, fitDeltaEhi], 'w--')

    if draw_inference_bins==True:
        print("Inference bins")
        for xlo,xhi,ylo,yhi in region_definitions['inference']:
            draw_box(xlo,xhi,ylo,yhi,fmt='y--')

    plt.xlabel(plt.gca().get_xlabel(), fontsize=18)
    plt.ylabel(plt.gca().get_ylabel(), fontsize=18)
    #plt.tight_layout()

    if tag is not None:
        plt.savefig(f'BNV_pLambda_plots/plot_{tag}_de_vs_mes.png')
        plt.tight_layout()

    signal_mask = (mes > region_definitions["signal MES"][0]) & ((DeltaE>region_definitions["signal DeltaE"][0]) & (DeltaE<region_definitions["signal DeltaE"][1])) 

    fit_mask = (mes > region_definitions["fitting MES"][0]) & ((DeltaE>region_definitions["fitting DeltaE"][0]) & (DeltaE<region_definitions["fitting DeltaE"][1])) 

    sideband1_mask = (mes > region_definitions["sideband MES"][0]) & ((DeltaE>region_definitions["sideband 1 DeltaE"][0]) & (DeltaE<region_definitions["sideband 1 DeltaE"][1])) 
    sideband2_mask = (mes > region_definitions["sideband MES"][0]) & ((DeltaE<region_definitions["sideband 2 DeltaE"][0]) & (DeltaE>region_definitions["sideband 2 DeltaE"][1])) 

    nsig = len(mes[signal_mask])
    nfit = len(mes[fit_mask])
    nside1 = len(mes[sideband1_mask])
    nside2 = len(mes[sideband2_mask])

    print(f'signal: {nsig}   fit: {nfit}  s1: {nside1}  s2: {nside2}  sum(s1,s2): {nside1 + nside2}')

    d = {}
    d['nsig'] = nsig
    d['nfit'] = nfit
    d['nside1'] = nside1
    d['nside2'] = nside2

    return d


##########################################################################
##########################################################################
def spherical_to_cartesian(r, costh, phi):

    #print(p3_spher)
    #print(costh)
    theta = np.arccos(costh)

    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    pmag = np.sqrt(x**2 + y**2 + z**2)
    #print(pmag, r)

    return x,y,z, pmag

##########################################################################
##########################################################################
def mass_from_spherical(p4s_spherical):

    etot,px,py,pz = 0,0,0,0

    for p4s_sph in p4s_spherical:
        #print("Here")
        #print(p4s_sph, p4s_sph[1:])
        x,y,z,pmag = spherical_to_cartesian(p4s_sph[1:])
        mass = p4s_sph[0]
        e = np.sqrt(mass**2 + pmag**2)

        etot += e
        px += x
        py += y
        pz += z

    #print("p4: ", etot, px, py, pz)
    m2 = etot**2 - (px**2 + py**2 + pz**2)
    if m2>=0:
        return np.sqrt(m2)
    else:
        return -np.sqrt(-m2)


##########################################################################
##########################################################################
def fill_new_entry_with_tag_side_B(data, final_state='pLambda'):
    '''
    # This is when I was trying to get the B-tag-side mass to have
    # the same shape as the other event shape variables, which
    # depend on the B candidate
    # Create a dummy set of data the size of our primary array
    x = -999*ak.ones_like(data['BSphr'], dtype=float)

    # We're going to flatten and then unflatten it so that it has the same shape
    # as the other event-shape variables
    xvals = ak.flatten(x) # This is when we thought we would have a jagged array
    n = ak.num(x)
    idx = ak.local_index(xvals)

    # Calculate the tag side B mass for entries with only one B and one Lambda
    mask_event_duplicates= get_duplicates_mask(data)
    m = tag_side_B(data[mask_event_duplicates], verbose=0, reverse_for_testing=False)

    # Numpy will be easier for what we want to do
    xvals = xvals.to_numpy()
    xvals[idx[mask_event_duplicates]] = m

    # Roll it back up to the original shape
    x = ak.unflatten(xvals, n)
    '''

    # Calculate the tag side B mass for entries with only one B and one Lambda
    mask_event_duplicates= get_duplicates_mask(data, final_state=final_state)
    m = tag_side_B(data[mask_event_duplicates], verbose=0, reverse_for_testing=False)

    x= -999*ak.ones_like(data['spmode'], dtype=float)
    idx = ak.local_index(x)

    # Numpy will be easier for what we want to do
    x = x.to_numpy()
    x[idx[mask_event_duplicates]] = m

    data['BtagSideMes'] = x

##########################################################################
##########################################################################
def tag_side_B(data, verbose=0, reverse_for_testing=False):

    if verbose:
        idx = data['TRKMCIdx'][0]
        mclund = data['mcLund'][0]
        nTRK = data['nTRK'][0]

        print("The MC and tracks for the first entry")

        print(f"nTRK: {nTRK}")
        for i,id in enumerate(idx):
            print(f"{i:2d}  {id:4d}   {mclund[id]}")
        print()

    lamd1idx = data['Lambda0d1Idx']
    lamd1Lund = data['Lambda0d1Lund']
    lamd2idx = data['Lambda0d2Idx']
    lamd2Lund = data['Lambda0d2Lund']

    TRKp3CM = data['TRKp3CM']
    TRKcosthCM = data['TRKcosthCM']
    TRKphiCM = data['TRKphiCM']

    gammap3CM = data['gammap3CM']
    gammacosthCM = data['gammacosthCM']
    gammaphiCM = data['gammaphiCM']

    #gammap3CM = data['pi0p3CM']
    #gammacosthCM = data['pi0costhCM']
    #gammaphiCM = data['pi0phiCM']

    if verbose:
        print(f'lamd1idx\n{lamd1idx}')
        print(f'lamd1Lund\n{lamd1Lund}')
        print(f'lamd2idx\n{lamd2idx}')
        print(f'lamd2Lund\n{lamd2Lund}')

    d2idx = data['Bd2Idx']
    d2Lund = data['Bd2Lund']

    if verbose:
        print()
        print(f'B d2idx\n{d2idx}')
        print(f'B d2Lund\n{d2Lund}')
        print()

    trkidx_proton = data['pTrkIdx']
    trkidx_pion = data['piTrkIdx']

    if verbose:
        print(f"# of protons: {data['np']}")
        print(f"trkidx_proton (the track index for labeled protons) \n{trkidx_proton}")
        print(f"# of pions: {data['npi']}")
        print(f"trkidx_pion (the track index for labeled pions) \n{trkidx_pion}")

        print()

    lamd1_trkidx = trkidx_proton[lamd1idx]
    lamd2_trkidx = trkidx_pion[lamd2idx]

    d2_trkidx = trkidx_proton[d2idx]

    if verbose:
        print(f"lamd1_trkidx\n{lamd1_trkidx}\n")
        print(f"lamd2_trkidx\n{lamd2_trkidx}\n")
        print(f"B d2_trkidx\n{d2_trkidx}\n")

    # This does not do things correctly when we have *two* candidates for either 
    # the lambda or B because it cuts out all of the tracks for both candidates but 
    # is assuming they are both part of one candidate
    bool_proton1 = indices_to_booleans(lamd1_trkidx, TRKp3CM)
    bool_pi = indices_to_booleans(lamd2_trkidx, TRKp3CM)

    bool_proton2 = indices_to_booleans(d2_trkidx, TRKp3CM)

    mask_tag_side_track = None
    if not reverse_for_testing:
        mask_tag_side_track = ~(bool_proton1 | bool_proton2 | bool_pi)
    else:
        # Use the candidates on the signal side
        mask_tag_side_track = (bool_proton1 | bool_proton2 | bool_pi)

    if verbose:
        print(f"bool_proton1: {bool_proton1}")
        print(f"bool_proton2: {bool_proton2}")
        print(f"bool_pi:      {bool_pi}")
        print(f"mask all:     {mask_tag_side_track}")

    TRKxCM, TRKyCM, TRKzCM, TRKpmagCM = spherical_to_cartesian(TRKp3CM[mask_tag_side_track], TRKcosthCM[mask_tag_side_track], TRKphiCM[mask_tag_side_track])
    gammaxCM, gammayCM, gammazCM, gammapmagCM = spherical_to_cartesian(gammap3CM, gammacosthCM, gammaphiCM)

    if verbose:
        print(f"TRKxCM: {TRKxCM}")
        print(f"TRKyCM: {TRKyCM}")
        print(f"TRKzCM: {TRKzCM}")
        print()
        print(f"gammaxCM: {gammaxCM}")
        print(f"gammayCM: {gammayCM}")
        print(f"gammazCM: {gammazCM}")
        print()

    TRKxCMtot = ak.sum(TRKxCM, axis=-1)
    TRKyCMtot = ak.sum(TRKyCM, axis=-1)
    TRKzCMtot = ak.sum(TRKzCM, axis=-1)

    gammaxCMtot = ak.sum(gammaxCM, axis=-1)
    gammayCMtot = ak.sum(gammayCM, axis=-1)
    gammazCMtot = ak.sum(gammazCM, axis=-1)

    pxCMtot = TRKxCMtot
    pyCMtot = TRKyCMtot
    pzCMtot = TRKzCMtot

    # Add in the photons
    if not reverse_for_testing:
        pxCMtot = pxCMtot + gammaxCMtot
        pyCMtot = pyCMtot + gammayCMtot
        pzCMtot = pxCMtot + gammazCMtot


    beamE = data['eeE']
    beamx = data['eePx']
    beamy = data['eePy']
    beamz = data['eePz']

    m = np.sqrt(beamE**2 - (beamx**2 + beamy**2 + beamz**2))

    if verbose:
        print(f'beam mass: {m}')
        print(f'beam E:    {beamE}')
        print(f'beam Px:   {beamx}')
        print(f'beam Py:   {beamy}')
        print(f'beam Pz:   {beamz}')

    Bpseudomass2 = (m/2)**2 - (pxCMtot**2 + pyCMtot**2 + pzCMtot**2)
    Bpseudomass = np.sqrt(Bpseudomass2)
    #mask_pos = Bpseudomass2>=0
    #Bpseudomass = -999*ak.ones_like(Bpseudomass2)
    #Bpseudomass[mask_pos] = np.sqrt(Bpseudomass2[mask_pos])
    #Bpseudomass[~mask_pos] = -np.sqrt(-Bpseudomass2[mask_pos])

    return Bpseudomass

##########################################################################
##########################################################################
##########################################################################################
def munge_mask_shapes(mask_larger, mask_smaller):

    # mask_larger will be bigger and have some Trues and Falses
    #
    # We want a mask the same size as mask_larger but where the Trues
    # have been replaces by the values in mask_smaller

    # Get the original indices
    idx = ak.local_index(mask_larger)
    
    # Make a mask of Falses. This will be the replacement mask
    mask = ak.zeros_like(mask_larger, dtype=bool)

    # It's a bit easier to work with numpy for this step
    mask = mask.to_numpy()

    # Get the indices where the larger mask is true and set it equal to the 
    # values of the smaller mask
    mask[idx[mask_larger]] = mask_smaller

    return mask

##########################################################################################
def get_final_masks(data_temp, region_definitions=None, tag="DEFAULT", IS_MC=True, is_BNC=False, final_state='pLambda'):

    if region_definitions is None:
        print("Need to pass in the region definitions")
        print("Exiting...")
        return None

    #dcuts = {"cut":[], "name":[], "event":[], "candidates":[]}
    dcuts = {}
    
    #########################################################################
    # First select events without duplicate candidates
    #########################################################################
    mask_event_duplicates= get_duplicates_mask(data_temp, final_state=final_state)

    ##########################################################################
    # Redefine the data array
    # tO only use the events with one candidate for Lambda0 and B
    ##########################################################################
    data = data_temp[mask_event_duplicates]

    #dcuts['cut'].append(1)
    dcuts[1] = {}
    dcuts[1]["name"] = "cut duplicates"
    dcuts[1]["event"] = mask_event_duplicates
    dcuts[1]["candidates"] = None
    
    
    #########################################################################
    # Events with Mes and DeltaE
    #########################################################################
    mask_candidates_fit = get_fit_mask(data, region_definitions=region_definitions)

    mask_event_fit_region = ak.any(mask_candidates_fit, axis=-1)

    # Munge the shape to ensure it is the same shape as the original data array
    mask_event_fit_region = munge_mask_shapes(mask_event_duplicates, mask_event_fit_region)

    # We need to make sure it is the same size as the original file
    
    dcuts[2] = {}
    dcuts[2]["name"] = "fitting region"
    dcuts[2]["event"] = mask_event_fit_region
    dcuts[2]["candidates"] = mask_candidates_fit
    
    
    #########################################################################
    # Cut on Lambda and number of candidates
    #########################################################################
    lambda_var = 'Lambda0postFitFlightSignificance'
    #lambda_var = 'Lambda0FlightLen'
    #lambda_var = 'Lambda0postFitFlightLen'
    mask_candidates_lambda0, mask_event_nlambda0_and_nB = get_lambda0_mask(data, \
                                                                    region_definitions=region_definitions, \
                                                                    flightlenvar=lambda_var)

    # Munge the shape to ensure it is the same shape as the original data array
    mask_event_nlambda0_and_nB = munge_mask_shapes(mask_event_duplicates, mask_event_nlambda0_and_nB)

    dcuts[3] = {} 
    dcuts[3]["name"] = "Lambda0 cuts / nB / nLambda"
    dcuts[3]["event"] = mask_event_nlambda0_and_nB
    dcuts[3]["candidates"] = mask_candidates_lambda0
    
    #########################################################################
    # PID cuts
    #########################################################################
    
    mask_bool_proton, mask_bool_pion, mask_bool_protonB = PID_masks(data, \
                  lamp_selector='SuperLooseKMProtonSelection', \
                  lampi_selector='VeryTightKMPionMicroSelection', \
                  Bp_selector='SuperTightKMProtonSelection', \
                  verbosity=0)

    #mask_event = mask_event_nlambda0_and_nB & mask_event_fit_region

    # We can do this because we selected the data where there is only one candidate for th
    # B and Lambda0
    mask_pid =      mask_bool_proton & \
                    mask_bool_pion & \
                    mask_bool_protonB

    mask_event_pid = ak.any(mask_pid, axis=-1)

    # Munge the shape to ensure it is the same shape as the original data array
    mask_event_pid = munge_mask_shapes(mask_event_duplicates, mask_event_pid)
    
    dcuts[4] = {} 
    dcuts[4]["name"] = "PID cuts"
    dcuts[4]["event"] = mask_event_pid
    dcuts[4]["candidates"] = mask_pid
    

    #########################################################################
    # Antiproton
    #########################################################################

    pps = myPIDselector.PIDselector("p")
    selector_to_test = "TightKMProtonSelection"
    
    # Is MC?
    # Note that we use the previous cuts on PID and fitting region

    mask_event_no_antiprotons, ct = build_antiproton_antimask(data, pps, selector_to_test, IS_MC=IS_MC, verbose=0)
    
    #mask_event = mask_no_antiprotons & mask_event_pid & mask_event_nlambda0_and_nB & mask_event_fit_region

    # Munge the shape to ensure it is the same shape as the original data array
    mask_event_no_antiprotons = munge_mask_shapes(mask_event_duplicates, mask_event_no_antiprotons)

    
    dcuts[6] = {}
    dcuts[6]["name"] = "antiproton cuts"
    dcuts[6]["event"] = mask_event_no_antiprotons
    dcuts[6]["candidates"] = None

    #########################################################################
    # Everything
    #########################################################################

    mask_event = mask_event_no_antiprotons & mask_event_pid & mask_event_nlambda0_and_nB & mask_event_fit_region

    dcuts[-1] = {}
    dcuts[-1]["name"] = "all"
    dcuts[-1]["event"] = mask_event
    dcuts[-1]["candidates"] = None
    
    return dcuts

############################################################################
############################################################################
##########################################################################################
def get_numbers_for_cut_flow(data, region_definitions=None,tag="DEFAULT", spmodes=None, verbose=False, dcuts=None):

    if region_definitions is None:
        print("Need to pass in the region definitions")
        print("Exiting...")
        return None

    # Check to make sure the user passed in some SP modes
    if spmodes is None:
        print("Will default to using all the spmodes in the file")
        spmodes = np.unique(data['spmode'].to_list())
        print(spmodes)
        print()

    if dcuts is None:
        # Need to get the original duplicates mask for any other cuts we might generate outside the function
        dcuts = get_final_masks(data, region_definitions=region_definitions)

    df_dict = {"cut":[], "name":[], "nevents":[], "pct":[], "tag":[], "spmode":[]}

    #########################################################################
    # Go through all the cuts
    #########################################################################
    for spmode in spmodes:
        mask_sp = data['spmode']==spmode

        #########################################################################
        # Org
        #########################################################################
        norg = len(data[mask_sp])

        df_dict["cut"].append(0)
        df_dict["name"].append("org")
        df_dict["nevents"].append(norg)
        df_dict["tag"].append(tag)
        df_dict["spmode"].append(spmode)

        if verbose:
            print(f"{norg} events in original file")

        df_dict["pct"].append(100*norg/norg)

        for key in dcuts.keys():
            if verbose:
                print(f'{key:3d} {dcuts[key]["name"]}')

            # Get the event mask
            mask_event = dcuts[key]['event']

            mask = mask_sp & mask_event
            n = len(mask[mask])

            df_dict["cut"].append(key)
            df_dict["name"].append(dcuts[key]["name"])
            df_dict["nevents"].append(n)
            df_dict["tag"].append(tag)
            df_dict["spmode"].append(spmode)
            df_dict["pct"].append(100*n/norg)


    df = pd.DataFrame.from_dict(df_dict)

    return df

############################################################################
############################################################################

##########################################################################
##########################################################################
def dump_awkward_to_dataframe(arr, fields_to_dump=None, write_to_filename=None):
    '''
    subset = ['spmode', 'BpostFitMes', 'BpostFitDeltaE', 'Lambda0_unc_Mass', \
          'BtagSideMes', 'BSphr', 'BThrust', 'BCosThetaS', \
          'R2', 'R2All', \
          'thrustMag', 'thrustMagAll', 'thrustCosTh', 'thrustCosThAll', 'sphericityAll', \
          'BCosSphr', 'BCosThetaT', 'BCosThrust', 'BLegendreP2', 'BR2ROE', 'BSphrROE', \
          'BThrustROE']
    '''
    if fields_to_dump is None:
        fields_to_dump = arr.fields

    # Get the type of array of a single dim array.
    # This should be more codified, but right now I am just hoping that
    # the first entry is something that should be flattened
    ak_array_type= type(arr[fields_to_dump[0]])

    df_dict={}
    nentries = -1
    for i,var in enumerate(fields_to_dump):
        x= arr[var] ##in each event, cut on the above cuts and pull out the info from each of the variables listed above
        if type(x[0]) == ak_array_type:
            x= ak.flatten(arr[var])
        df_dict[var] = x
        if i==0:
            nentries = len(x)
        else:
            if nentries != len(x):
                print("Irregular number of entries!")
                print("Not filling the dictionary!\n")
                for key,vals in df_dict.items():
                    print(f"{key:16s}  {len(vals)}")
                print(f"{val:16d} {len(x)}")

    df_out= pd.DataFrame.from_dict(df_dict)

    #if dropna==True:
    #    df_out.dropna(inplace=True)

    if write_to_filename is not None:
        df_out.to_parquet(write_to_filename)

    return df_out



##########################################################################
##########################################################################
#def punzi_fom_nn(model_aft_train, sp_data, threshold, sp_998_df, sp_999_df, sig_disc= 4, scaling= 0.3):
def punzi_fom_nn(df_sp, df_col, sig_sp_mode='-999', region_definitions = None, sigma = 4.0, BNC=False):

    # Collision data
    mask = (df_col['cut_-1'] == True) 
    if BNC is True:
        mask = (df_col['cut_2'] == True) 
        mask = mask & (df_col['cut_3'] == True) 
        mask = mask & (df_col['cut_4'] == True) 
        
    
    df_col_tmp = df_col[mask]

    # SP
    mask = (df_sp['cut_-1'] == True) 
    if BNC is True:
        mask = (df_sp['cut_2'] == True) 
        mask = mask & (df_sp['cut_3'] == True) 
        mask = mask & (df_sp['cut_4'] == True) 

    mask = mask & (df_sp['spmode'] == sig_sp_mode)
    mask = mask & (df_sp['used_in_sig_train'] == False)
    df_sp_tmp = df_sp[mask]

    meslo = region_definitions['signal MES'][0]
    meshi = region_definitions['signal MES'][1]
    
    delo = region_definitions['signal DeltaE'][0]
    dehi = region_definitions['signal DeltaE'][1]

    messidelo = region_definitions['sideband MES'][0]
    messidehi = region_definitions['sideband MES'][1]
    
    desidelo1 = region_definitions['sideband 1 DeltaE'][0]
    desidehi1 = region_definitions['sideband 1 DeltaE'][1]
    
    desidelo2 = region_definitions['sideband 2 DeltaE'][0]
    desidehi2 = region_definitions['sideband 2 DeltaE'][1]

    # Print statements
    print(f'{meslo = }        {meshi = }')
    print(f'{messidelo = }    {messidehi = }')
    print(f'{delo = }         {dehi = }')
    print(f'{desidelo1 = }     {desidehi1 = }')
    print(f'{desidelo2 = }     {desidehi2 = }')
    
    fom_dict = {}
    fom_dict['thresh'] = []
    fom_dict['nbkg_sb1'] = []
    fom_dict['nbkg_sb2'] = []
    fom_dict['nbkg'] = []
    fom_dict['nsig'] = []

    # Collision data
    mes_col = df_col_tmp['BpostFitMes']
    de_col = df_col_tmp['BpostFitDeltaE']

    mask1_col = (mes_col>messidelo) & (mes_col<messidehi) & (de_col>desidelo1) & (de_col<desidehi1)    
    mask2_col = (mes_col>messidelo) & (mes_col<messidehi) & (de_col>desidelo2) & (de_col<desidehi2)

    # SP
    mes_sp = df_sp_tmp['BpostFitMes']
    de_sp = df_sp_tmp['BpostFitDeltaE']

    mask_sp = (mes_sp>meslo) & (mes_sp<meshi) & (de_sp>delo) & (de_sp<dehi) 

    for thresh in np.arange(0,1,0.01):
        
        # Collision data
        mask_thresh_col = df_col_tmp['proba'] > thresh

        nsb1 = len(df_col_tmp[mask1_col & mask_thresh_col])        
        nsb2 = len(df_col_tmp[mask2_col & mask_thresh_col])
    
        # Collision data
        mask_thresh_sp = df_sp_tmp['proba'] > thresh

        nsig = len(df_sp_tmp[mask_sp & mask_thresh_sp])        
    
        #print(nsb1, nsb2, nsig)
        
        fom_dict['thresh'].append(thresh)
        fom_dict['nbkg_sb1'].append(nsb1)
        fom_dict['nbkg_sb2'].append(nsb2)
        #fom_dict['nbkg'].append((nsb1 + nsb2)/2)
        fom_dict['nbkg'].append(nsb1 + nsb2)

        fom_dict['nsig'].append(nsig)

    df_fom = pd.DataFrame.from_dict(fom_dict)
    df_fom['sig_pct'] = df_fom['nsig'] / df_fom['nsig'].iloc[0]

    # Number of signal estimation for S / sqrt(S+B)
    N_S0 = 50
    df_fom['N_S'] = N_S0*df_fom['sig_pct']


    #sigma = 4.0
    
    df_fom['fom'] = df_fom['sig_pct'] / (np.sqrt(df_fom['nbkg']) + sigma/2.0)
    df_fom['fom_std'] = df_fom['N_S'] / np.sqrt(df_fom['N_S'] + df_fom['nbkg'])

    return df_fom
##########################################################################
##########################################################################
def calculate_conversion_factor(df_sp, df_col, decay='BNV', region_definitions=None, sig_eff_after_ML=1.0):

    nBpairs = 470.89e6
    # Assuming 0.28%
    #nBpairs_err =  1.32
    # Assuming 0.6%
    nBpairs_err =  2.83e6

    #nB_bf = [0.484, 0.484, 0.516, 0.516, 0.516, 0.516, 0.516]
    #nB_bf_err = [0.006,0.006,0.006,0.006,0.006,0.006]

    # From PDG
    nB_bf = [0.514]
    nB_bf_err = [0.006]

    # initial numbers for signal SP
    #skim_eff = [0.467, 0.504, 0.553, 0.569, 0.553, 0.569]
    #skim_eff_err = [0.001, 0.001, 0.001, 0.001, 0.001, 0.001]

    # BNC
    #skim_eff = [0.16] # Is this the efficiency after the skim?
    skim_eff_due_to_ML = [sig_eff_after_ML] # Is this the efficiency after the skim?
    skim_eff_due_to_ML_err = [0.0015] # Assuming that we go from ~60000 events to about 20000 after the ML selection


    ######### HOW ARE THESE DIFFERENT FROM ABOVE?
    # initial numbers for signal SP
    #initial_numbers = [22000, 22000, 28000, 28000, 28000, 28000]
    #final_numbers =   [12387, 11223, 14536, 13371, 15867, 14760]

    # These are the numbers before and after BtaTupleMaker
    # BNC
    initial_numbers = [95243]
    final_numbers_after_BTM =   [44905]
    final_numbers =   [44905]
    if decay=='BNV':
        initial_numbers = [2*99966] # I put 2 runs of ~100k 
        # After BtaTupleMaker
        final_numbers_after_BTM =   [119215]
        # After rectangular cuts
        final_numbers =   [59937]

    print("At the skim stage...")
    print(f"initial numbers (before BTM):       {initial_numbers[0]}")
    print(f"final   numbers (after BTM):        {final_numbers_after_BTM[0]}")
    print(f"final   numbers (after rect cuts):  {final_numbers[0]}")
    print(f"eff (after BTM)                  :  {final_numbers_after_BTM[0]/initial_numbers[0]:.4f}%")
    print(f"eff (due to rectcuts)            :  {final_numbers[0]/final_numbers_after_BTM[0]:.4f}%")
    print(f"eff (after rectcuts)             :  {final_numbers[0]/initial_numbers[0]:.4f}%")
    print()

    nmodes = len(initial_numbers)

    # Baryon branching fractions
    baryon_bf =     [0.639]
    baryon_bf_err = [0.005]

    # Tracking errors
    # http://www.slac.stanford.edu/BFROOT/www/Physics/TrackEfficTaskForce/TauEff/R24/TauEff.html
    trk_err_per_trk = 0.128/100.0 # percent error
    trk_err_pct_l0 = 3.0 * trk_err_per_trk # 3 tracks
    trk_err_pct_lc = 4.0 * trk_err_per_trk # 4 tracks

    #trk_pct_err = [trk_err_pct_lc, trk_err_pct_lc, trk_err_pct_l0, trk_err_pct_l0, trk_err_pct_l0, trk_err_pct_l0]
    trk_pct_err = [trk_err_pct_l0]

    # PID errors
    pid_err_p  = 0.010
    pid_err_pi = 0.010
    pid_err_k  = 0.012
    #pid_err_e  = 0.004
    #pid_err_mu = 0.007
    pid_err_e  = 0.01
    pid_err_mu = 0.025

    pid_pct_err = []
    # Precise
    #pid_pct_err.append(math.sqrt(pid_err_p*pid_err_p + pid_err_pi*pid_err_pi + pid_err_k*pid_err_k + pid_err_mu*pid_err_mu))
    #pid_pct_err.append(math.sqrt(pid_err_p*pid_err_p + pid_err_pi*pid_err_pi + pid_err_k*pid_err_k + pid_err_e*pid_err_e))
    #pid_pct_err.append(math.sqrt(pid_err_p*pid_err_p + pid_err_pi*pid_err_pi + pid_err_mu*pid_err_mu))
    #pid_pct_err.append(math.sqrt(pid_err_p*pid_err_p + pid_err_pi*pid_err_pi + pid_err_e*pid_err_e))
    #pid_pct_err.append(math.sqrt(pid_err_p*pid_err_p + pid_err_pi*pid_err_pi + pid_err_mu*pid_err_mu))
    #pid_pct_err.append(math.sqrt(pid_err_p*pid_err_p + pid_err_pi*pid_err_pi + pid_err_e*pid_err_e))

    # BNC or BNV
    pid_pct_err.append(math.sqrt(pid_err_p*pid_err_p + pid_err_p*pid_err_p + pid_err_pi*pid_err_pi))


    # Estimate
    #pid_pct_err.append(0.025)
    #pid_pct_err.append(0.025)
    #pid_pct_err.append(0.020)
    #pid_pct_err.append(0.020)
    #pid_pct_err.append(0.020)
    #pid_pct_err.append(0.020)

    # Eff calculations
    conv_factors = []
    conv_factor_errs = []
    total_efficiencies = []
    total_efficiency_errs = []
    for i in range(0,nmodes):
        n0 = initial_numbers[i]
        n =  final_numbers[i]
        eff = n/float(n0)
        # https://lss.fnal.gov/archive/test-tm/2000/fermilab-tm-2286-cd.pdf
        #eff_err = math.sqrt((eff*(1.0-eff))/n0)
        #eff_err = (1/n0) * math.sqrt(n*(1.0-n/n0))
        eff_err = (1/n0) * math.sqrt(n*(1.0-n/n0))

        pre_ML_eff = eff
        pre_ML_eff_err = eff_err

        # Total efficiency due to 
        # BTM
        # Rectangular cuts
        # ML cut
        eff *= skim_eff_due_to_ML[i]

        conv_factor = (nBpairs*2.0*nB_bf[i]) * eff * baryon_bf[i]

        # Calculate all the percent errors. 
        pct_errs = []

        # number of Bs
        pct_errs.append(nBpairs_err/float(nBpairs))
        # B branching fraction
        pct_errs.append(nB_bf_err[i]/float(nB_bf[i]))

        # Efficiency
        pct_errs.append(skim_eff_due_to_ML_err[i]/skim_eff_due_to_ML[i])
        pct_errs.append(pre_ML_eff_err/pre_ML_eff)

        # Branching fractions
        pct_errs.append(baryon_bf_err[i]/baryon_bf[i])

        # Tracking
        pct_errs.append(trk_pct_err[i])
        # PID?
        pct_errs.append(pid_pct_err[i])

        eff_tot_err =  (eff_err/eff)*(eff_err/eff)
        eff_tot_err += trk_pct_err[i]*trk_pct_err[i]
        eff_tot_err += pid_pct_err[i]*pid_pct_err[i]

        tot_pct_err = 0.0
        for pe in pct_errs:
            tot_pct_err += pe*pe
            #print "%f %f %f %f" % (tot_pct_err, math.sqrt(tot_pct_err), pe*pe, pe)

        conv_factor_err = math.sqrt(tot_pct_err)
        #print "conv_factor_err: %f" % (conv_factor_err)

        # Convert back to a number, rather than a percentage
        conv_factor_err *= conv_factor

        total_eff = eff
        #total_eff_err = eff*math.sqrt(eff_tot_err)
        # Going to use the total fractional uncertainties to calculate the signal efficiency
        total_eff_err = total_eff * (conv_factor_err/conv_factor)

        output = "%d\ttrk_pct_err: %6.4f\n" % (i, trk_pct_err[i])
        output += " \tpid_pct_err: %6.4f\n" % (pid_pct_err[i])
        output += " \tnBpairs_err: %6.4f\n" % (nBpairs_err/float(nBpairs))
        output += " \tnB_bf_err: %6.4f\n" % (nB_bf_err[i]/float(nB_bf[i]))
        output += " \tbaryon_bf_err: %6.4f\n\n" % (baryon_bf_err[i]/baryon_bf[i])

        output += f"\tpre_ML_eff: {pre_ML_eff:6.4f} +/- {pre_ML_eff_err:6.4f} ({pre_ML_eff_err/pre_ML_eff:6.4f} frac uncert) (we'll use this for the MC stat error)\n"
        output += f"\ttotal eff:  {eff:6.4f} +/- {total_eff_err:6.4f}\n\n" 
                #(pre_skim_eff,eff_err, eff,eff*math.sqrt(eff_tot_err))
        output += f"\tconv_factor: {conv_factor:6.2f} +/- {conv_factor_err:6.3f} (pct_err: {100*conv_factor_err/conv_factor:6.3f})\n"
                #(conv_factor,conv_factor_err,100*conv_factor_err/conv_factor)

        print(output)
        print()
        output = f"conv_factor: {conv_factor:6.4e} +/- {conv_factor_err:6.3e} (pct_err: {100*conv_factor_err/conv_factor:6.3f})"
        print(output)

        conv_factors.append(conv_factor)
        conv_factor_errs.append(conv_factor_err)
        total_efficiency_errs.append(total_eff_err)
        total_efficiencies.append(total_eff)

    return total_efficiencies, total_efficiency_errs, conv_factors, conv_factor_errs

##########################################################################
##########################################################################

#""""
def add_probas_to_dfs(workspace, df_col, df_sp):

    model = workspace['model']
    x_train = workspace['x_train']
    y_train = workspace['y_train']
    x_test = workspace['x_test']
    y_test = workspace['y_test']

    #idx_x_train     = workspace['idx_sig_train']
    #idx_x_not_train = workspace['idx_sig_not_train']

    # Get the training vars
    training_vars = model.feature_names

    print("Training vars initially: ")
    print(training_vars)
    print(len(training_vars))

    remove_these = ['proba', 'used_in_sig_train', 'used_in_bkg_train']
    #print(type(training_vars))
    for r in remove_these:
        print(f"Checking for {r}...")
        if sum(training_vars.isin([r]).astype(int))>0:
            idx = training_vars.get_loc(r)
            print(idx, training_vars)
            training_vars = training_vars.delete(idx)
            print(training_vars)
            print(f"removed {r}")

    print("Training vars finally: ")
    print(training_vars)
    print(len(training_vars))


    #x_train = df_col[training_vars]#[idx_x_train]

    #### FOR SP
    #for i in threshold:
    #    output_df= see_stuff(sig_samp= sig_samp,bkg_samp= bkg_samp, thresh= i, verbose= False, df=MC_data, df_col=coll_data)

    # 3. Scale the test data using the same scaler
    #scaler = StandardScaler()

    # SHOULD THIS BE THE ACTUAL TRAINING SP
    #x_dummy = scaler.fit_transform(x_train)

    #print(x_dummy.T[-1])
    #print(len(x_dummy), len(x_dummy[0]))

    #print("training vars after" )
    #print(training_vars)

    # Collision
    x_test  = df_col[training_vars].values#[idx_x_not_train]
    print(x_test.T[-1])

    #x_test = scaler.transform(x_test)

    proba = model.predict_proba(x_test)

    #print("proba: ")
    #print(proba)

    df_col['proba'] = proba[:,0]

    #print(df_col.columns)

    # SP
    # Collision
    x_test  = df_sp[training_vars].values#[idx_x_not_train]
    #x_test = scaler.transform(x_test)

    #print(x_test)

    proba = model.predict_proba(x_test)

    #print("proba: ")
    #print(proba)

    df_sp['proba'] = proba[:,0]

    #print(df_sp.columns)

    idx_bkg_train = workspace['idx_bkg_train']
    idx_sig_train = workspace['idx_sig_train']

    nentries = len(df_sp)
    df_sp['used_in_sig_train'] = np.zeros(nentries, dtype=bool)
    df_sp.loc[idx_sig_train, "used_in_sig_train"] = True

    df_sp['used_in_bkg_train'] = np.zeros(nentries, dtype=bool)
    df_sp.loc[idx_bkg_train, 'used_in_bkg_train'] = True


    #x_col_proba = model.predict_proba(df_col_tmp)

    #### FOR COLLISION
    #y_proba_col_sig = model.predict_proba()

    #sp998= sp_data["spmode"]== "998"
    #N_bkg = len(sp_998_df[sp998]) ## total number of background events (sp 998)
    #signal_before= len(sp_999_df)
    #signal_after= len(sp_999_df)
    #efficiency = signal_after/signal_before ## the accuracy of the model after training with the SP
    #fom = efficiency(threshold)/(np.sqrt(N_bkg(threshold)+sig_disc/2))
    #return fom
    return 0
#"""
