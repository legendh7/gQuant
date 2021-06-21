# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 01:28:41 2021

@author: jmh
"""

import pandas as pd
import math
import numpy as np
import re
import os
from datetime import date
import time

#used for deconvolution and deisotope
PROTON = 1.00727646677

elements = ''
glycan_units = {'HexNAc':{'name': 'HexNAc', 'mass':203.07937251951, 'element_composition': 'C8H13O5N0' },
                'Hex':{'name': 'Hex', 'mass':162.0528234185, 'element_composition': 'C6H10O5N1' },
                'Sia':{'name': 'Sia', 'mass':291.09541650647, 'element_composition': 'C11H17O8N1' },
                'dHex':{'name': 'dHex', 'mass':146.05790879894, 'element_composition': 'C6H10O4N0' } }

glycan_derivatization = {'2_AB':{'name': '2_AB', 'mass':119.06037, 'element_composition': 'C7H7N2O' },
                         '2_AA':{'name': '2_AA', 'mass':120.04439, 'element_composition': 'C7H6N1O1' }}

Na_mz = 22.98922
K_mz = 38.96316
H_mz = 1.00728
H2O_mz = 18.01002

def calc_iso_pattern():
    '''
    designed to calculate distinct glycan isotope patterns with given
    '''
    pass

def gen_glycan_elements(g):
    '''
    to generate element composition for given glycan
    '''
    pass

def gen_glycan_db(max_glycan_mass):
    '''
    '''
    smallest_n_glycan = 'N2H3'
    glycan_pattern = ''
    # rules, least glycan structure = N2H3
    # 

human_db = 'Human_glycan_344.qdb'
mamalian_db = 'Mammalian_glycan_419.qdb'

def load_glycan_db(db_to_choose):
    if db_to_choose == 'Human':
        return pd.read_csv(human_db)
    elif db_to_choose == 'Mammalian':
        return pd.read_csv(mamalian_db)
    else:
        raise Exception('Invalid glycan database')

def glycan_db_with_mods(initial_gdb, derivatization = {}, mode = 'positive', adducts = {'H':1.0078,'Na':22.9898},is_full_de = False):
    '''
    To generate glycan database with given adducts(charge carriers), charges, chemical modifications(reduce end/free end, derivatization)
    @param: 
    initial_gdb: gdb load from qdb file
    adducts: adduct list specify by user
    max_charge: charges to consider, also specify by user, default value is 1
    mods: derivatization list specified by user, 2-AB, free, 2-AA, etc
    
    @return: custom gdb with above considerations
    note that in current version, alteration above was not considered in the isotope ratio calculation as it was manually generated and as the smaller size of known modifications.
    '''
    
    l = len(initial_gdb)
    new_gdb = pd.DataFrame(columns = ['Hex', 'HexNAc', 'NeuAc', 'Xyl', 'dHex', 'composition',
       'mw_residue', 'mw_free_end', 'element_compos_free_end',
       'derivative', 'adduct', 'mw', 'first', 'second', 'third',
       'fourth', 'fifth', 'sixth', '1st', '2nd', '3rd', '4th', '5th',
       '6th'])

    if mode == 'positive':
        # then m/z set as MH+ or MNa+ etc
        if derivatization:
            # derivatization is not empty, only one derivatization allowed，not permethylation or per etherlation allowed
            if is_full_de:
                #consider only full deri
                for i in range(l):
                    for adduct in adducts:
                        #should make sure derivatization is one
                        origin_row = initial_gdb.iloc[i,:].copy()
                        # assume derivatization not 100% derivative set as "-"

                        origin_row['adduct'] = adduct
                        
                        origin_row['derivative'] = list(derivatization.keys())[0]
                        mw = origin_row['mw_free_end'] + list(derivatization.values())[0] + adducts[adduct]
                        origin_row['mw'] = mw
                        
                        #add one copy of derivatized glycan
                        new_gdb = new_gdb.append(origin_row, ignore_index=True)
            else:        
                for i in range(l):
                    for adduct in adducts:
                        #should make sure derivatization is one
                        origin_row = initial_gdb.iloc[i,:].copy()
                        # assume derivatization not 100% derivative set as "-"
                        origin_row['derivative'] = 'FreeEnd'
                        origin_row['adduct'] = adduct
                        mw = origin_row['mw_free_end'] + adducts[adduct]
                        origin_row['mw'] = mw
                        
                        #add one copy of underivatized glycan 
                        new_gdb = new_gdb.append(origin_row, ignore_index=True)
                        
                        
                        origin_row['derivative'] = list(derivatization.keys())[0]
                        mw = origin_row['mw_free_end'] + list(derivatization.values())[0] + adducts[adduct]
                        origin_row['mw'] = mw
                        
                        #add one copy of derivatized glycan
                        new_gdb = new_gdb.append(origin_row, ignore_index=True)
                    
                    
        else:
            # no derivatization is performed on glycans only adducts(charge carrier considered)
            for i in range(l):
                for adduct in adducts:
                    origin_row = initial_gdb.iloc[i,:].copy()
                    # assume derivatization not 100% derivative set as "-"
                    origin_row['derivative'] = 'FreeEnd'
                    origin_row['adduct'] = adduct
                    mw = origin_row['mw_free_end'] + adducts[adduct]
                    origin_row['mw'] = mw
                    
                    #add one copy of underivatized glycan 
                    new_gdb = new_gdb.append(origin_row, ignore_index=True)


    if mode == 'negative':
        # then m/z set as M-H or M+Na-2H etc to be negatively charged
        if derivatization:
            # derivatization is not empty, only one derivatization allowed，not permethylation or per etherlation allowed 
                if is_full_de:
                    for i in range(l):
                        for adduct in adducts:
                            #should make sure derivatization is one
                            origin_row = initial_gdb.iloc[i,:].copy()
                            # assume derivatization not 100% derivative set as "-"
    
                            origin_row['adduct'] = adduct
                                                    
                            origin_row['derivative'] = list(derivatization.keys())[0]
                            mw = origin_row['mw_free_end'] + list(derivatization.values())[0] + adducts[adduct] - 2*PROTON
                            origin_row['mw'] = mw
                            
                            #add one copy of derivatized glycan
                            new_gdb = new_gdb.append(origin_row, ignore_index=True)

                else:
                    for i in range(l):
                        for adduct in adducts:
                            #should make sure derivatization is one
                            origin_row = initial_gdb.iloc[i,:].copy()
                            # assume derivatization not 100% derivative set as "-"
                            origin_row['derivative'] = '-'
                            origin_row['adduct'] = adduct
                            mw = origin_row['mw_free_end'] + adducts[adduct] - 2*PROTON
                            origin_row['mw'] = mw
                            
                            #add one copy of underivatized glycan 
                            new_gdb = new_gdb.append(origin_row, ignore_index=True)
                            
                            
                            origin_row['derivative'] = list(derivatization.keys())[0]
                            mw = origin_row['mw_free_end'] + list(derivatization.values())[0] + adducts[adduct] - 2*PROTON
                            origin_row['mw'] = mw
                            
                            #add one copy of derivatized glycan
                            new_gdb = new_gdb.append(origin_row, ignore_index=True)
                    
                    
        else:
            # no derivatization is performed on glycans only adducts(charge carrier considered)
            for i in range(l):
                for adduct in adducts:
                    origin_row = initial_gdb.iloc[i,:].copy()
                    # assume derivatization not 100% derivative set as "-"
                    origin_row['derivative'] = '-'
                    origin_row['adduct'] = adduct
                    mw = origin_row['mw_free_end'] + adducts[adduct] - 2*PROTON
                    origin_row['mw'] = mw
                    
                    #add one copy of underivatized glycan 
                    new_gdb = new_gdb.append(origin_row, ignore_index=True)
                    
    return new_gdb

def load_ms_data(ms_data_directory, note_sign = '#', header = None, sep = ',', skip_line= None):
    '''
    '''    
    signal = []
    with open(ms_data_directory) as f:
        fo = f.readlines()
        signal = [] #can also set as dict or nparray
        
        
        for line in fo:
            # remove '\n'
            line = line[:-1]
            # skip blanck line
            if len(line.split(sep)) == 0:
                continue
            # skip lines starts with a-z A-Z
            if line.startswith(r'T'):
                continue
            # skip identifier
            elif line.split()[0] == note_sign:
                continue
            
            # to capture the m/z vs. signal pattern, careaobut sep type
            elif len(line.split(sep)) == 2 and isinstance(float(line.split(sep)[0]), float) and isinstance(float(line.split(sep)[1]), float):
                    mz, intensity = line.split(sep)
                    signal.append(((float(mz)), float(intensity)))
        
    if len(signal) == 0:
        raise Exception('Data not found!')
    
    return signal
    
    
def deisotope(spect, max_charge, max_next_peaks = 8, min_next_peaks = 3, mz_tol_ppm = 30.0, iso_shift = 1.00287 ):
    '''
    This function was developed to figure out the isotope envelop patterns within a spectrum
    
    @spect: spectrum containing peak and intensity pairs in a dictionary form {'peak1':intensity1,'peak2': intensity2....}
    @max_charge: precursor charge as reported from vendor instrument
    @max_next_peak: for seaking the possible isotope envelop and summing intensity of final reported monoisotope peaks.
    @min_next_peak: (smallest number - 1) of peaks required for an envelop
    @mz_tol_ppm: peak m/z tolerence factor in ppm unit
    @iso_shift: equal to a H+ add_up, [MH]+ --> [MHn]n+
    
    @return: a list with tuples containing monoisotope peaks, summed intensities and coresponding charges
             peaks without any isotope was disgarded
             
    @note: single peak without isotope found was also included in finally list with intensity retained and charge set as 'undefined'
    
    '''
    
    # make sure that spectrum was not empty
    assert(len(spect) > 0)
    
    # two tuples
#    mz_list, intens_list = zip(*sorted(spect.items()))
    # currently a list
    mz_list, intens_list = list(zip(*sorted(spect)))
    
    # initiated monoisotope list
    monoisotope_peak_list = []
    full_isotope_envelope_list = []
    temp_isotope_envelope_list = []    
    
    length = len(spect)
    
    #negetive mode if necessary
    if max_charge < 0:
        charges = [-x for x in range(int(abs(max_charge)), 0, -1)]        
    # positve mode
    else:
        charges = [ x for x in range(int(max_charge), 0, -1) ]
    
    # Sequential search through the peaklist to find the most possible isotope envelop
    # for peak deisotope
    
    for i in range(1, length):
        temp_isotope_envelope_list = []
        
        for z in charges:
            
            found = False
            
            # i = 1.2.....
            j = i -1
            suspect_peak = float(mz_list[i]) - float(iso_shift / abs(z))
            mz_tol = suspect_peak * mz_tol_ppm / 1000000 
                            
            while j >= 0:
                if suspect_peak - mz_tol <= mz_list[j] <= suspect_peak + mz_tol and intens_list[i]/3 < intens_list[j] < intens_list[i]*3:
                    found = True
                    break
                j = j - 1
                    
            # if preceding peak found, jump out of charge loop and to next peak
            if found:
                break
            
            found = 1
            last_intensity = intens_list[i]
            m = 1
                        
            local_max = False
            
            for i_envelop in range(1, max_next_peaks + 1):
                if i + i_envelop >= length:
                    break
                
                target_mz = float(mz_list[i]) + iso_shift * i_envelop / abs(z)
                
                # how to calculate charge based tolerance in this case, solved, treated as peak * 20 ppm since it charge independent.
                has_peak_result = has_peak_list(spect, target_mz, mz_tol_ppm)
                
                if has_peak_result == False:
                    break
                else:
                    
                    # check if the isotope fit an envelop with only one highest peak
                    m += 1
                    mz, intensity = has_peak_result
                    #original is 3, due to the fact that relative quant 1:10 lead to lower former and later peak ,consider enlarge it to higher 2021.04.07
                    if  last_intensity / 3 < intensity < last_intensity * 3:
                        local_max = True
                    elif local_max == True and intensity > last_intensity:
                        break
                    
                    found += 1
                    temp_isotope_envelope_list.append(tuple((mz_list[i + i_envelop], intens_list[i + i_envelop], float(z), found)))
                    last_intensity = intensity
                    
            if found >=3:                
                # return a tuple containing three dimension information
                monoisotope_peak_list.append(tuple((mz_list[i], intens_list[i], float(z), found)))
                #not sure why it always missing, add the monoisotope here
                full_isotope_envelope_list.append((mz_list[i], intens_list[i], float(z), 1, found))
                for single_isotope in temp_isotope_envelope_list:
                    full_isotope_envelope_list.append(single_isotope)
                #initiate temp envelope peak list
                
                break
            # loop go to last charge and still no isotope found. set it as undefine
            # Note that z == 1 can not be absent otherwise lower charge peaks will be missed
            elif z ==1 and i_envelop == 1:
#                monoisotope_peak_list.append(tuple((mz_list[i], intensity_sum, 'undefine' , found)))
                break
            
    return monoisotope_peak_list, full_isotope_envelope_list



def deisotope_profile(spect, max_charge, mz_calibrate = 0.0, max_next_peaks = 8, min_next_peaks = 3, mz_tol_ppm = 30.0, iso_shift = 1.00287 ):
    '''
    This function was developed to figure out the isotope envelop patterns within a spectrum
    
    @spect: spectrum containing peak and intensity pairs in a dictionary form {'peak1':intensity1,'peak2': intensity2....}
    @max_charge: precursor charge as reported from vendor instrument
    @max_next_peak: for seaking the possible isotope envelop and summing intensity of final reported monoisotope peaks.
    @min_next_peak: (smallest number - 1) of peaks required for an envelop
    @mz_tol_ppm: peak m/z tolerence factor in ppm unit
    @iso_shift: equal to a H+ add_up, [MH]+ --> [MHn]n+
    
    @return: a list with tuples containing monoisotope peaks, summed intensities and coresponding charges
             peaks without any isotope was disgarded
             
    @note: single peak without isotope found was also included in finally list with intensity retained and charge set as 'undefined'
    
    '''
    
    # make sure that spectrum was not empty
    assert(len(spect) > 0)
    
    # two tuples
#    mz_list, intens_list = zip(*sorted(spect.items()))
    # currently a list
    mz_list, intens_list = list(zip(*sorted(spect)))
    
    if mz_calibrate:
        mz_list = [m_to_cal + mz_calibrate for m_to_cal in mz_list]
        
    # initiated monoisotope list
    monoisotope_peak_list = []
    full_isotope_envelope_list = []
    temp_isotope_envelope_list = []    
    
    length = len(spect)
    
    #negetive mode if necessary
    if max_charge < 0:
        charges = [-x for x in range(int(abs(max_charge)), 0, -1)]        
    # positve mode
    else:
        charges = [ x for x in range(int(max_charge), 0, -1) ]
    
    # Sequential search through the peaklist to find the most possible isotope envelop
    # for peak deisotope
    
    for i in range(0, length):
        
        
        
        for z in charges:
            temp_isotope_envelope_list = []
            found = False
            
            # i = 1.2.....
            j = i -1
            suspect_peak = float(mz_list[i]) - float(iso_shift / abs(z))
            mz_tol = suspect_peak * mz_tol_ppm / 1000000 
                            
            while j >= 0:
                if suspect_peak - mz_tol <= mz_list[j] <= suspect_peak + mz_tol and intens_list[i]/3 < intens_list[j] < intens_list[i]*3:
                    found = True
                    break
                j = j - 1
                    
            # if preceding peak found, current peak was not monoisotopic peak, jump out of charge loop and to next peak
            if found:
                break
            
            found = 1
            last_intensity = intens_list[i]

            
            for i_envelop in range(1, max_next_peaks + 1):
                if i + i_envelop >= length:
                    break
                
                target_mz = float(mz_list[i]) + iso_shift * i_envelop / abs(z)
                
                # how to calculate charge based tolerance in this case, solved, treated as peak * 20 ppm since it charge independent.
                has_peak_result = has_peak_list(spect, target_mz, mz_tol_ppm)
                
                
                if has_peak_result == False:
                    break
                else:
                    
                    # check if the isotope fit an envelop with only one highest peak

                    curr_mz, curr_intensity = has_peak_result
                    #original is 3, due to the fact that relative quant 1:10 lead to lower former and later peak ,consider enlarge it to higher 2021.04.07
                    if last_intensity / 3 < curr_intensity < last_intensity * 3:
                        found += 1
                        temp_isotope_envelope_list.append(tuple((curr_mz, curr_intensity, float(z), found)))
                        last_intensity = curr_intensity
                    else:
                        break
                    

                           
            if found >=3:                
                # return a tuple containing three dimension information

                monoisotope_peak_list.append(tuple((mz_list[i], intens_list[i], float(z), found)))
                #not sure why it always missing, add the monoisotope here
                full_isotope_envelope_list.append((mz_list[i], intens_list[i], float(z), 1, found))
                
                for single_isotope in temp_isotope_envelope_list:
                    full_isotope_envelope_list.append(single_isotope)
                #initiate temp envelope peak list
                
                break
            # loop go to last charge and still no isotope found. set it as undefine
            # Note that z == 1 can not be absent otherwise lower charge peaks will be missed
            elif z ==1 and i_envelop == 1:
#                monoisotope_peak_list.append(tuple((mz_list[i], intensity_sum, 'undefine' , found)))
                break
            
    return monoisotope_peak_list, full_isotope_envelope_list


def has_peak_list(spect_list, mz, mz_tol_ppm = 20):
    
    '''
    this function was developed to confirm whether a peak is pesented
    in a spectrum
    @spect: spect in dict form
    @mz: mz to search
    @mz_tol_ppm: peak m/z tolerence factor in ppm unit
    
    @return tuple containing peak and intensity i.e. (peak, intensity), if not found ,reuturn False
    
    '''
    
    # no peaks contained
    if len(spect_list) == 0:
        return False
    
#==============================================================================
#     # only one peak contained
#     elif len(spect_list) == 1:
#         return spect[0][0]
#==============================================================================
    
    # binary search
    peak_lst = sorted(spect_list)    
    low = 0
    high = len(spect_list) - 1
    mz_tol = mz * mz_tol_ppm / 1000000
              
    while low <= high:
        mid = int( (low + high) / 2 )
        if mz < peak_lst[mid][0] and peak_lst[mid][0] - mz > mz_tol:
            high = mid - 1
            
        elif mz > peak_lst[mid][0] and mz - peak_lst[mid][0] > mz_tol:
            low = mid + 1
            
        else:
            return peak_lst[mid]
    # peaks not found
    return False    

def has_glycan(spect_list, mz, mz_tol_ppm = 20):
    
    '''
    this function was developed to confirm whether a peak is pesented
    in a spectrum
    @spect: spect in dict form
    @mz: mz to search
    @mz_tol_ppm: peak m/z tolerence factor in ppm unit
    
    @return tuple containing peak and intensity i.e. (peak, intensity), if not found ,reuturn False
    
    '''
    
    # no peaks contained
    if len(spect_list) == 0:
        return False
    
#==============================================================================
#     # only one peak contained
#     elif len(spect_list) == 1:
#         return spect[0][0]
#==============================================================================
    
    # binary search
    peak_lst = sorted(spect_list)    
    low = 0
    high = len(spect_list) - 1
    mz_tol = mz * mz_tol_ppm / 1000000
              
    while low <= high:
        mid = int( (low + high) / 2 )
        if mz < peak_lst[mid][0] and peak_lst[mid][0] - mz > mz_tol:
            high = mid - 1
            
        elif mz > peak_lst[mid][0] and mz - peak_lst[mid][0] > mz_tol:
            low = mid + 1
            
        else:
            return peak_lst[mid]
    # peaks not found
    return False    

    
def match_glycan_ppm(peaks, gdb_with_adduct, mz_tol, delta_mass =0.0000, mz_calibr = 0.0):
    '''
    design to search gdb for glycan mapping
    search for both low- and heavy-(if delta_mass > 0) labeled glycan
    
    @ parma: mz_calibr: used to calibr mass in case of mass shift. set as 0.0 if 
    @return: peak list with glycan assignment
    
    '''
    
    match_results_ppm = pd.DataFrame(columns = ['m/z', 'charge', 'intensity','glycan', 'derivatization', 'adduct', 'composition','experimental_mass',
   'theoretical_mass', 'mz_diff', 'mz_diff_ppm', 'paired_m/z', 'paired_intensity', 'L/H ratio', 'label_type', 'quant_H_mz', 'quant_L_mz', 'quant_H_intensity', 'quant_L_intensity', 'calibr_ratio','coment' ])

    assert(len(peaks)>0)
#    charge = 1.0

    i = 0
    for peak in peaks:
        #peak in the form of tuple (m/z, intensity, charge, num_isotopes)
        #assumed that H+addcuts, should not be a big problem with maldi with singaly charged
        peak_mass = peak[0] * peak[2]
        
        
        # peak should excceed 910.33248
        if peak_mass< 910.0:
            continue
        
        # may change to binary search in the future to accelerate match speed
        if not delta_mass:
            # no need to search heavy labeled
            for i in range(len(gdb_with_adduct)):
                
                mz_diff = peak_mass - gdb_with_adduct.iloc[i,:]['mw']
                
                mz_diff_ppm = abs(peak_mass - gdb_with_adduct.iloc[i,:]['mw'])/gdb_with_adduct.iloc[i,:]['mw'] * 1000000.0
                
                # For maldi normaly one charge form found
                glycan = gdb_with_adduct.iloc[i,:]['composition']+ '+' + gdb_with_adduct.iloc[i,:]['derivative']+ '+' + str(int(peak[2])) + '*'+ gdb_with_adduct.iloc[i,:]['adduct']  +'+light'
                if mz_diff_ppm <= mz_tol:
                    match_results_ppm = match_results_ppm.append(pd.DataFrame({'m/z':[peak[0]], "charge":[peak[2]],'intensity':[peak[1]], 'glycan':[glycan],'composition':[gdb_with_adduct.iloc[i,:]['composition']],
                                                                       'derivatization':[gdb_with_adduct.iloc[i,:]['derivative']], 'adduct':[gdb_with_adduct.iloc[i,:]['adduct']],'theoretical_mass': [gdb_with_adduct.iloc[i,:]['mw']],
                                                                       'experimental_mass':[peak_mass], "mz_diff":[mz_diff],"mz_diff_ppm":[mz_diff_ppm],'label_type':['L']}), ignore_index=True)
        else:
            # need to search heavy labeled glycan
            for i in range(len(gdb_with_adduct)):
                
                mz_diff_light = peak_mass - gdb_with_adduct.iloc[i,:]['mw']
                mz_diff_heavy = peak_mass - gdb_with_adduct.iloc[i,:]['mw'] - delta_mass
                
                mz_diff_light_ppm = abs(mz_diff_light)/gdb_with_adduct.iloc[i,:]['mw'] * 1000000.0
                mz_diff_heavy_ppm = abs(mz_diff_heavy)/gdb_with_adduct.iloc[i,:]['mw'] * 1000000.0
                
                # For maldi normaly one charge form found
                glycan_light = gdb_with_adduct.iloc[i,:]['composition']+ '+' + gdb_with_adduct.iloc[i,:]['derivative']+ '+' + str(int(peak[2])) + '*'+ gdb_with_adduct.iloc[i,:]['adduct'] +'+light'
                glycan_heavy = gdb_with_adduct.iloc[i,:]['composition']+ '+' + gdb_with_adduct.iloc[i,:]['derivative']+ '+' + str(int(peak[2])) + '*'+ gdb_with_adduct.iloc[i,:]['adduct'] +'+heavy'
                
                if mz_diff_light_ppm <= mz_tol:
                    match_results_ppm = match_results_ppm.append(pd.DataFrame({'m/z':[peak[0]], "charge":[peak[2]],'intensity':[peak[1]], 'glycan':[glycan_light],'composition':[gdb_with_adduct.iloc[i,:]['composition']],
                                                                       'derivatization':[gdb_with_adduct.iloc[i,:]['derivative']], 'adduct':[gdb_with_adduct.iloc[i,:]['adduct']],'theoretical_mass': [gdb_with_adduct.iloc[i,:]['mw']],
                                                                       'experimental_mass':[peak_mass], "mz_diff":[mz_diff_light],"mz_diff_ppm":[mz_diff_light_ppm],'label_type':['L']}), ignore_index=True)
                if mz_diff_heavy_ppm <= mz_tol:
                    match_results_ppm = match_results_ppm.append(pd.DataFrame({'m/z':[peak[0]], "charge":[peak[2]],'intensity':[peak[1]], 'glycan':[glycan_heavy],'composition':[gdb_with_adduct.iloc[i,:]['composition']],
                                                                       'derivatization':[gdb_with_adduct.iloc[i,:]['derivative']], 'adduct':[gdb_with_adduct.iloc[i,:]['adduct']],'theoretical_mass': [gdb_with_adduct.iloc[i,:]['mw']],
                                                                       'experimental_mass':[peak_mass], "mz_diff":[mz_diff_heavy],"mz_diff_ppm":[mz_diff_heavy_ppm],'label_type':['H']}), ignore_index=True)
                                                                      
    return match_results_ppm


def match_glycan_dalton(peaks, gdb_with_adduct, mz_tol, mz_calibr = 0.0):
    '''
    design to serarch gdb for glycan mapping
    
    @ parma: mz_calibr: used to calibr mass in case of mass shift. set as 0.0 if 
    @return: peak list with glycan assignment
    
    '''
    
    match_results = pd.DataFrame(columns = ['m/z', 'charge', 'intensity','glycan', 'derivatization', 'adduct', 'composition','experimental_mass',
   'theoretical_mass', 'delta_mass', 'paired_m/z', 'paired_intensity', 'ratio' ])

    assert(len(peaks)>0)

    i = 0
    for peak in peaks:
        #peak in the form of tuple (m/z, intensity, charge, num_isotopes)
        #assumed that H+addcuts, should not be a big problem with maldi with singaly charged
        peak_mass = peak[0] * peak[2]
        
        # may change to binary search in the future to accelerate match speed
        for i in range(len(gdb_with_adduct)):
            
            delta_mass = peak_mass - gdb_with_adduct.iloc[i,:]['mw'] 

            glycan = gdb_with_adduct.iloc[i,:]['composition']+gdb_with_adduct.iloc[i,:]['derivative']+str(int(peak[2])) + gdb_with_adduct.iloc[i,:]['adduct']
            if delta_mass <= mz_tol:
                match_results = match_results.append(pd.DataFrame({'m/z':[peak[0]], "charge":[peak[2]],'intensity':[peak[1]], 'glycan':[glycan],'composition':[gdb_with_adduct.iloc[i,:]['composition']],
                                                                   'derivatization':[gdb_with_adduct.iloc[i,:]['derivative']], 'adduct':[gdb_with_adduct.iloc[i,:]['adduct']],'theoretical_mass': [gdb_with_adduct.iloc[i,:]['mw']],
                                                                   'experimental_mass':[peak_mass], "delta_mass":[delta_mass]}), ignore_index=True)

                
    return match_results

def find_peak_list(spect_list, mz, mz_tol_ppm = 20):
    
    '''
    this function was developed to confirm whether a peak is pesented
    in a spectrum
    @spect: spect in dict form
    @mz: mz to search
    @mz_tol_ppm: peak m/z tolerence factor in ppm unit
    
    @return tuple containing peak and intensity i.e. (peak, intensity), if not found ,reuturn False
    
    '''
    
    # no peaks contained
    if len(spect_list) == 0:
        return False
    
    # binary search
    peak_lst = sorted(spect_list)    
    low = 0
    high = len(spect_list) - 1
    mz_tol = mz * mz_tol_ppm / 1000000
              
    while low <= high:
        mid = int( (low + high) / 2 )
        if mz < peak_lst[mid][0] and peak_lst[mid][0] - mz > mz_tol:
            high = mid - 1
            
        elif mz > peak_lst[mid][0] and mz - peak_lst[mid][0] > mz_tol:
            low = mid + 1
            
        else:
            return peak_lst[mid]
    # peaks not found
    return False    

def gquant(match_results, full_envelop_list, gdb_with_adduct, spect_peak, delta_mass, quant_mode = 'highest' ):
    '''
    quantiation function of gquant
    >params:
    @matche_results: matched glycans and peak from
    @spect: input spect in centroding form (list of peak and intensity)
    @gdb_with_adduct: used to fetch glycan isotope pattern
    @delta_mass: important for MS based relative quantitation (2H 2.0151 3H 3.02293 Da)
    @mode: by monoisotopic or highest
    
    >return: dataframe containing quantitation results and function.    
    '''
    
    if delta_mass == 0:
        # no heavy or light condition, will add new function in the future
        match_results['L/H ratio'] = 'N/A'
        return match_results
    
    if delta_mass < 0:
        raise Exception('Delta mass for relative quantitation must be positve')
    
    for i in match_results.index:
        matched_mz = match_results.loc[i,'m/z']
        glycan = match_results.loc[i,'composition']
        qlabel = match_results.loc[i,'label_type']
        if qlabel == 'L':
            
            low_mz = matched_mz
            #matched glycan equals to 
            #locate and return match paired peak
            gnv = locate_envelope(matched_mz, full_envelop_list, delta_mass)
            
            if gnv:
                # i.e. peak with delta_mass interval found
                high_mz = gnv[0]
                
                # assined values
                match_results.loc[i,'paired_m/z'] = high_mz
                match_results.loc[i,'paired_intensity'] = gnv[1]

                #for ratio calculation and calibration
                glycan_info =  gdb_with_adduct.loc[gdb_with_adduct.composition == glycan]
                glycan_info_list = np.array(glycan_info)
                
                if quant_mode == 'highest':
                    #quant with hihgest isotope
                    l_to_h_ratio, ql_mz, ql_intens, qh_mz, qh_intens, calibr_r = gquant_highest(full_envelop_list,glycan_info_list ,low_mz, high_mz)

                    #in case no enough isotopes to support quantitation
                    match_results.loc[i,'L/H ratio'] = l_to_h_ratio                    
                    match_results.loc[i,'quant_L_mz'] = ql_mz
                    match_results.loc[i,'quant_L_intensity'] = ql_intens
                    
                    match_results.loc[i,'quant_H_mz'] = qh_mz
                    match_results.loc[i,'quant_H_intensity'] = qh_intens
                    
                    match_results.loc[i,'calibr_ratio'] = calibr_r
                    
                if quant_mode == 'mono':
                    #quant with hihgest isotope
                    l_to_h_ratio, ql_mz, ql_intens, qh_mz, qh_intens, calibr_r = gquant_normal(full_envelop_list,glycan_info_list ,low_mz, high_mz)

                    #in case no enough isotopes to support quantitation
                    match_results.loc[i,'L/H ratio'] = l_to_h_ratio                    
                    match_results.loc[i,'quant_L_mz'] = ql_mz
                    match_results.loc[i,'quant_L_intensity'] = ql_intens
                    
                    match_results.loc[i,'quant_H_mz'] = qh_mz
                    match_results.loc[i,'quant_H_intensity'] = qh_intens
                    
                    match_results.loc[i,'calibr_ratio'] = calibr_r
            else:
                # only light peak found
                match_results.loc[i,'L/H ratio'] = 0
                    
        if qlabel == 'H':
                #matched glycan equals to 
                #locate and return match paired peak

                high_mz = matched_mz
                
                #heavy labeled glycan matched, search for light glycan

                low_mz = high_mz - abs(delta_mass)
                #find low_mz and calculate L/H ratio                
                gnv = has_peak_list_for_H( spect_peak, low_mz)
                
                if gnv:
                    # i.e. peak with delta_mass interval found
                    low_mz = gnv[0]
                    
                    # assined values
                    match_results.loc[i,'paired_m/z'] = low_mz
                    match_results.loc[i,'paired_intensity'] = gnv[1]
    
                    #for ratio calculation and calibration
                    glycan_info =  gdb_with_adduct.loc[gdb_with_adduct.composition == glycan]
                    glycan_info_list = np.array(glycan_info)
                    
                    if quant_mode == 'highest':
                        #quant with hihgest isotope
                        l_to_h_ratio, ql_mz, ql_intens, qh_mz, qh_intens, calibr_r = gquant_highest_H(spect_peak, glycan_info_list ,low_mz, high_mz)
    
                        #in case no enough isotopes to support quantitation
                        match_results.loc[i,'L/H ratio'] = l_to_h_ratio 
                        match_results.loc[i,'quant_L_mz'] = ql_mz
                        match_results.loc[i,'quant_L_intensity'] = ql_intens
                        
                        match_results.loc[i,'quant_H_mz'] = qh_mz
                        match_results.loc[i,'quant_H_intensity'] = qh_intens
                        
                        match_results.loc[i,'calibr_ratio'] = calibr_r
                        
                    if quant_mode == 'mono':
                        l_to_h_ratio, ql_mz, ql_intens, qh_mz, qh_intens, calibr_r = gquant_normal_H(spect_peak, glycan_info_list ,low_mz, high_mz)
    
                        #in case no enough isotopes to support quantitation
                        match_results.loc[i,'L/H ratio'] = l_to_h_ratio 
                        match_results.loc[i,'quant_L_mz'] = ql_mz
                        match_results.loc[i,'quant_L_intensity'] = ql_intens
                        
                        match_results.loc[i,'quant_H_mz'] = qh_mz
                        match_results.loc[i,'quant_H_intensity'] = qh_intens
                        
                        match_results.loc[i,'calibr_ratio'] = calibr_r
                        

                else:
                    # only heavy peak found
                    match_results.loc[i,'L/H ratio'] = np.inf
    

    order = ['m/z', 'charge', 'intensity','composition','L/H ratio','glycan', 'derivatization', 'adduct', 'experimental_mass', 'theoretical_mass', 'mz_diff', 'mz_diff_ppm', 'paired_m/z', 'paired_intensity', 'label_type', 'quant_H_mz', 'quant_L_mz', 'quant_H_intensity', 'quant_L_intensity', 'calibr_ratio','coment' ]
    match_results = match_results[order]      
    
    match_results.sort_values(by = ['m/z', 'mz_diff_ppm', 'derivatization'])
    
    return match_results

def has_peak_list_for_H(sp, mz_to_found, mz_tol_ppm = 20):
    
    '''
    this function was developed to confirm whether a peak is pesented
    in a spectrum
    @low_mz_to_found: low mz to found
    @sp: spect in list form
    @mz_tol_ppm: peak m/z tolerence factor in ppm unit
    
    @return tuple containing peak and intensity i.e. (peak, intensity), if not found ,reuturn False
    
    '''
    
    # no peaks contained
    if len(sp) == 0:
        return False
    
#==============================================================================
#     # only one peak contained
#     elif len(spect_list) == 1:
#         return spect[0][0]
#==============================================================================
    
    # binary search
    peak_lst_to_find = sorted(sp)    
    low = 0
    high = len(sp) - 1

    mz_tol = mz_to_found * mz_tol_ppm / 1000000
              
    while low <= high:
        mid = int( (low + high) / 2 )
        if mz_to_found < peak_lst_to_find[mid][0] and peak_lst_to_find[mid][0] - mz_to_found > mz_tol:
            high = mid - 1
            
        elif mz_to_found > peak_lst_to_find[mid][0] and mz_to_found - peak_lst_to_find[mid][0] > mz_tol:
            low = mid + 1
            
        else:
            return peak_lst_to_find[mid]
    # peaks not found
    return False

def gquant_highest_H(sp,glycan_info_list, low_mz, high_mz):
    '''
    
    in current case, delta mass was integral folds of 1.0 Da
    @return: calibr Heavy to light ratio
    
    '''

#    glycan_info =  glycan_db.loc[glycan_db.composition == 'glycan_structure']
#    
    glycan_mz = glycan_info_list[0][7]
    
    next_count = 0    
    current_next_isotope = has_peak_list_for_H(sp, low_mz)
    paried_next_isotope = has_peak_list_for_H(sp, high_mz)
    
    # introduce new parameter of isotope interval
    # isotope interval equal was used to check whether isotope inteference works if it was small, calibr_ratio calculation was peform to calibrate final ratio
    isotope_interval = abs(int(high_mz - low_mz))
    
    if isotope_interval >= 15:
        #Then isotope inteference can be ignore
        
        # choose the highest epitope
        if glycan_mz < 2250:
            #monoisotope chosen as quant peak
            calibr_ratio = 0.0
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio 
            
    
        if 2250 < glycan_mz < 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 1
            current_next_isotope_1 = has_peak_list_for_H(sp, low_mz+PROTON)
            paried_next_isotope_1 = has_peak_list_for_H(sp, high_mz+PROTON) 
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_1 and paried_next_isotope_1:
                calibr_ratio = 0.0
                return (paried_next_isotope_1[1] - calibr_ratio*current_next_isotope_1[1]) / current_next_isotope_1[1], current_next_isotope_1[0], current_next_isotope_1[1], paried_next_isotope_1[0], paried_next_isotope_1[1], calibr_ratio
                
            
            # The second isotopic peaks chosen for quantitation
    
        elif glycan_mz >= 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 2
            
            current_next_isotope_2 = has_peak_list_for_H(sp, low_mz+PROTON*2)
            paried_next_isotope_2 = has_peak_list_for_H(sp, high_mz+PROTON*2)
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_2 and paried_next_isotope_2:
                calibr_ratio = 0.0                
                return (paried_next_isotope_2[1] - calibr_ratio*current_next_isotope_2[1]) / current_next_isotope_2[1], current_next_isotope_2[0], current_next_isotope_2[1], paried_next_isotope_2[0], paried_next_isotope_2[1], calibr_ratio
                
    elif isotope_interval == 14:
        if glycan_mz < 2250:
            # could still have isotope value from gdb
            calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio
            #monoisotope chosen as quant peak
    
        if 2250 < glycan_mz < 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 1
            current_next_isotope_1 = has_peak_list_for_H(sp, low_mz+PROTON)
            paried_next_isotope_1 = has_peak_list_for_H(sp, high_mz+PROTON) 
            
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_1 and paried_next_isotope_1:
                calibr_ratio = 0.0
                return (paried_next_isotope_1[1] - calibr_ratio*current_next_isotope_1[1]) / current_next_isotope_1[1], current_next_isotope_1[0], current_next_isotope_1[1], paried_next_isotope_1[0], paried_next_isotope_1[1], calibr_ratio

            # The second isotopic peaks chosen for quantitation
    
        elif glycan_mz >= 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 2
            current_next_isotope_2 = has_peak_list_for_H(sp, low_mz+PROTON*2)
            paried_next_isotope_2 = has_peak_list_for_H(sp, high_mz+PROTON*2)
            
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_2 and paried_next_isotope_2:
                #drawback of this algorithum that the database only contain 6 isotopes can report i2 and i7 might have some quantitaion bias for rel labeling interval higher than 4
                calibr_ratio = 0.0
                return (paried_next_isotope_2[1] - calibr_ratio*current_next_isotope_2[1]) / current_next_isotope_2[1], current_next_isotope_2[0], current_next_isotope_2[1], paried_next_isotope_2[0], paried_next_isotope_2[1], calibr_ratio
                
    
    elif isotope_interval ==13:
        if glycan_mz < 2250:
            # could still have isotope value from gdb
            calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio
            #monoisotope chosen as quant peak
    
        if 2250 < glycan_mz < 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 1
            current_next_isotope_1 = has_peak_list_for_H(sp, low_mz+PROTON)
            paried_next_isotope_1 = has_peak_list_for_H(sp, high_mz+PROTON) 
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_1 and paried_next_isotope_1:
                calibr_ratio = glycan_info_list[0][19 + isotope_interval] /glycan_info_list[0][19]
                return (paried_next_isotope_1[1] - calibr_ratio*current_next_isotope_1[1]) / current_next_isotope_1[1], current_next_isotope_1[0], current_next_isotope_1[1], paried_next_isotope_1[0], paried_next_isotope_1[1], calibr_ratio

            # The second isotopic peaks chosen for quantitation
    
        elif glycan_mz >= 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 2
            current_next_isotope_2 = has_peak_list_for_H(sp, low_mz+PROTON*2)
            paried_next_isotope_2 = has_peak_list_for_H(sp, high_mz+PROTON*2)
            calibr_ratio = 0.0
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_2 and paried_next_isotope_2:
                
                return (paried_next_isotope_2[1] - calibr_ratio*current_next_isotope_2[1]) / current_next_isotope_2[1], current_next_isotope_2[0], current_next_isotope_2[1], paried_next_isotope_2[0], paried_next_isotope_2[1], calibr_ratio

    elif 0 < isotope_interval < 13:
        if glycan_mz < 2250:
            # could still have isotope value from gdb
            
            calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
            final_ratio = (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1]
            return final_ratio, current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio
            #monoisotope chosen as quant peak
    
        if 2250 < glycan_mz < 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 1

            current_next_isotope_1 = has_peak_list_for_H(sp, low_mz+PROTON)
            paried_next_isotope_1 = has_peak_list_for_H(sp, high_mz+PROTON)            
#            current_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
#            paried_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_1 and paried_next_isotope_1:

                calibr_ratio = glycan_info_list[0][19 + isotope_interval] /glycan_info_list[0][19]
                return (paried_next_isotope_1[1] - calibr_ratio*current_next_isotope_1[1]) / current_next_isotope_1[1], current_next_isotope_1[0], current_next_isotope_1[1], paried_next_isotope_1[0], paried_next_isotope_1[1], calibr_ratio

            # The second isotopic peaks chosen for quantitation
    
        elif glycan_mz >= 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 2
            current_next_isotope_2 = has_peak_list_for_H(sp, low_mz+PROTON*2)
            paried_next_isotope_2 = has_peak_list_for_H(sp, high_mz+PROTON*2)  
#            current_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
#            paried_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_2 and paried_next_isotope_2:

                calibr_ratio = glycan_info_list[0][20 + isotope_interval] /glycan_info_list[0][20]
                return (paried_next_isotope_2[1] - calibr_ratio*current_next_isotope_2[1]) / current_next_isotope_2[1], current_next_isotope_2[0], current_next_isotope_2[1], paried_next_isotope_2[0], paried_next_isotope_2[1], calibr_ratio

    
    #if used l0 and h0 for calculation, calibr_ratio is fixed
    if isotope_interval <= 15:

        calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
    else:
        calibr_ratio = 0.0
    final_ratio = (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1]
    # return by default 
    return final_ratio,current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio     

def gquant_normal_H(sp,glycan_info_list ,low_mz, high_mz):
    '''
    
    in current case, delta mass was integral folds of 1.0 Da
    @return: calibr Heavy to light ratio
    
    '''

#    glycan_info =  glycan_db.loc[glycan_db.composition == 'glycan_structure']
#    
    glycan_mz = glycan_info_list[0][7]
     
    current_next_isotope = has_peak_list_for_H(sp, low_mz)
    paried_next_isotope = has_peak_list_for_H(sp, high_mz)
    
    # introduce new parameter of isotope interval
    # isotope interval equal was used to check whether isotope inteference works if it was small, calibr_ratio calculation was peform to calibrate final ratio
    isotope_interval = abs(int(high_mz - low_mz))
    
    if isotope_interval >= 15:
        #Then isotope inteference can be ignore
        
        # choose the highest epitope
        if glycan_mz < 2250:
            #monoisotope chosen as quant peak
            calibr_ratio = 0.0
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio 

    elif 0< isotope_interval <= 14:

            calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio
            #monoisotope chosen as quant peak

    #if used l0 and h0 for calculation, calibr_ratio is fixed
    if isotope_interval <= 15:

        calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
    else:
        calibr_ratio = 0.0
    final_ratio = (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1]
    # return by default 
    return final_ratio,current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio

         
def locate_envelope(mz, full_envelope, delta_mass):
    '''
    to locate matched mz and its paired mz
    @return: tuple list of isotopes
    '''    
    if delta_mass > 0:
        
        for i in range(len(full_envelope)):
        # equal float might cause some problem and should be replaced in futre
        
            #search right
            if full_envelope[i][0] < mz:
                continue
            # Need to consider the error mapping ,QIT report >6ppm and 5800 ~3ppm for paired
            elif abs(full_envelope[i][0] -  mz - delta_mass) / mz * 1000000.0 < 10.0 :
                return full_envelope[i]
                
    if delta_mass < 0:
        for i in range(len(full_envelope)):
        # equal float might cause some problem and should be replaced in futre
        
            #search left
            if full_envelope[i][0] > mz:
                continue
            # Need to consider the error mapping ,QIT report >6ppm and 5800 ~3ppm for paired
            elif abs(full_envelope[i][0] -  mz - delta_mass) / mz * 1000000.0 < 10.0 :
                return full_envelope[i]
                
    return None


def next_isotope_in_the_same_envelope(full_envelope, current_mz, next_count):
    '''
    return next envelop for quantitation
    @return: if found return next, if not return None
    '''    
    
    for i in range(len(full_envelope)-next_count):
        
        if full_envelope[i][0] == current_mz:
            # Need to return next_count i+1 i+2 i+3 etc
            if full_envelope[i+next_count][3] - full_envelope[i][3]  == next_count:
                return full_envelope[i+next_count]            
    return None        

def gquant_highest(full_envelop_list,glycan_info_list ,low_mz, high_mz):
    '''
    
    in current case, delta mass was integral folds of 1.0 Da
    @return: calibr Heavy to light ratio
    
    '''

#    glycan_info =  glycan_db.loc[glycan_db.composition == 'glycan_structure']
#    
    glycan_mz = glycan_info_list[0][7]
    
    next_count = 0    
    current_next_isotope = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
    paried_next_isotope = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
    
    # introduce new parameter of isotope interval
    # isotope interval equal was used to check whether isotope inteference works if it was small, calibr_ratio calculation was peform to calibrate final ratio
    isotope_interval = abs(int(high_mz - low_mz))
    
    if isotope_interval >= 15:
        #Then isotope inteference can be ignore
        
        # choose the highest epitope
        if glycan_mz < 2250:
            #monoisotope chosen as quant peak
            calibr_ratio = 0.0
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio 
            
    
        if 2250 < glycan_mz < 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 1
            current_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
            paried_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_1 and paried_next_isotope_1:
                calibr_ratio = 0.0
                return (paried_next_isotope_1[1] - calibr_ratio*current_next_isotope_1[1]) / current_next_isotope_1[1], current_next_isotope_1[0], current_next_isotope_1[1], paried_next_isotope_1[0], paried_next_isotope_1[1], calibr_ratio
                
            
            # The second isotopic peaks chosen for quantitation
    
        elif glycan_mz >= 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 2
            current_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
            paried_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_2 and paried_next_isotope_2:
                calibr_ratio = 0.0                
                return (paried_next_isotope_2[1] - calibr_ratio*current_next_isotope_2[1]) / current_next_isotope_2[1], current_next_isotope_2[0], current_next_isotope_2[1], paried_next_isotope_2[0], paried_next_isotope_2[1], calibr_ratio
                
    elif isotope_interval == 14:
        if glycan_mz < 2250:
            # could still have isotope value from gdb
            calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio
            #monoisotope chosen as quant peak
    
        if 2250 < glycan_mz < 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 1
            current_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
            paried_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_1 and paried_next_isotope_1:
                calibr_ratio = 0.0
                return (paried_next_isotope_1[1] - calibr_ratio*current_next_isotope_1[1]) / current_next_isotope_1[1], current_next_isotope_1[0], current_next_isotope_1[1], paried_next_isotope_1[0], paried_next_isotope_1[1], calibr_ratio

            # The second isotopic peaks chosen for quantitation
    
        elif glycan_mz >= 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 2
            current_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
            paried_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_2 and paried_next_isotope_2:
                #drawback of this algorithum that the database only contain 6 isotopes can report i2 and i7 might have some quantitaion bias for rel labeling interval higher than 4
                calibr_ratio = 0.0
                return (paried_next_isotope_2[1] - calibr_ratio*current_next_isotope_2[1]) / current_next_isotope_2[1], current_next_isotope_2[0], current_next_isotope_2[1], paried_next_isotope_2[0], paried_next_isotope_2[1], calibr_ratio
                
    
    elif isotope_interval ==13:
        if glycan_mz < 2250:
            # could still have isotope value from gdb
            calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio
            #monoisotope chosen as quant peak
    
        if 2250 < glycan_mz < 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 1
            current_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
            paried_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_1 and paried_next_isotope_1:
                calibr_ratio = glycan_info_list[0][19 + isotope_interval] /glycan_info_list[0][19]
                return (paried_next_isotope_1[1] - calibr_ratio*current_next_isotope_1[1]) / current_next_isotope_1[1], current_next_isotope_1[0], current_next_isotope_1[1], paried_next_isotope_1[0], paried_next_isotope_1[1], calibr_ratio

            # The second isotopic peaks chosen for quantitation
    
        elif glycan_mz >= 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 2
            current_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
            paried_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            calibr_ratio = 0.0
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_2 and paried_next_isotope_2:
                
                return (paried_next_isotope_2[1] - calibr_ratio*current_next_isotope_2[1]) / current_next_isotope_2[1], current_next_isotope_2[0], current_next_isotope_2[1], paried_next_isotope_2[0], paried_next_isotope_2[1], calibr_ratio

    elif 0 < isotope_interval < 13:
        if glycan_mz < 2250:
            # could still have isotope value from gdb

            calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
            final_ratio = (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1]
            return final_ratio, current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio
            #monoisotope chosen as quant peak
    
        if 2250 < glycan_mz < 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 1
            
            current_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
            paried_next_isotope_1 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_1 and paried_next_isotope_1:

                calibr_ratio = glycan_info_list[0][19 + isotope_interval] /glycan_info_list[0][19]
                return (paried_next_isotope_1[1] - calibr_ratio*current_next_isotope_1[1]) / current_next_isotope_1[1], current_next_isotope_1[0], current_next_isotope_1[1], paried_next_isotope_1[0], paried_next_isotope_1[1], calibr_ratio

            # The second isotopic peaks chosen for quantitation
    
        elif glycan_mz >= 3850:
            # if evelope not sufficent enough to performed return None
            next_count = 2
            current_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
            paried_next_isotope_2 = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
            
            #make sure that it found both next isotope, if not, return original ratio
            if current_next_isotope_2 and paried_next_isotope_2:

                calibr_ratio = glycan_info_list[0][20 + isotope_interval] /glycan_info_list[0][20]
                return (paried_next_isotope_2[1] - calibr_ratio*current_next_isotope_2[1]) / current_next_isotope_2[1], current_next_isotope_2[0], current_next_isotope_2[1], paried_next_isotope_2[0], paried_next_isotope_2[1], calibr_ratio

    
    #if used l0 and h0 for calculation, calibr_ratio is fixed
    if isotope_interval <= 15:

        calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
    else:
        calibr_ratio = 0.0
    final_ratio = (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1]
    # return by default 
    return final_ratio,current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio     


def gquant_normal(full_envelop_list,glycan_info_list ,low_mz, high_mz):
    '''
    
    in current case, delta mass was integral folds of 1.0 Da
    @return: calibr Heavy to light ratio
    
    '''

#    glycan_info =  glycan_db.loc[glycan_db.composition == 'glycan_structure']
#    
    glycan_mz = glycan_info_list[0][7]
    
    next_count = 0    
    current_next_isotope = next_isotope_in_the_same_envelope(full_envelop_list, low_mz, next_count)
    paried_next_isotope = next_isotope_in_the_same_envelope(full_envelop_list, high_mz, next_count)
    
    # introduce new parameter of isotope interval
    # isotope interval equal was used to check whether isotope inteference works if it was small, calibr_ratio calculation was peform to calibrate final ratio
    isotope_interval = abs(int(high_mz - low_mz))
    
    if isotope_interval >= 15:
        #Then isotope inteference can be ignore
        
        # choose the highest epitope
        if glycan_mz < 2250:
            #monoisotope chosen as quant peak
            calibr_ratio = 0.0
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio 

    elif 0< isotope_interval <= 14:

            calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
            return (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1], current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio
            #monoisotope chosen as quant peak

    #if used l0 and h0 for calculation, calibr_ratio is fixed
    if isotope_interval <= 15:

        calibr_ratio = glycan_info_list[0][18 + isotope_interval] /glycan_info_list[0][18]
    else:
        calibr_ratio = 0.0
    final_ratio = (paried_next_isotope[1] - calibr_ratio*current_next_isotope[1]) / current_next_isotope[1]
    # return by default 
    return final_ratio,current_next_isotope[0], current_next_isotope[1], paried_next_isotope[0], paried_next_isotope[1], calibr_ratio 
    
def merge_quant_results():
    pass

def filter_quant_results(spect, quant_res, percent_intensity_threshold):
    pass

def cal_signal_to_noise_ratio(spect, mode = 'mean'):
    
    '''
    signal to noise ratio = signal / noise_level
    '''
    
    
    
    
    
    pass

def peak_integration():
    pass

def integrate(a, b, n, func):
    
    '''
    '''
    
    h = (b-a) / float(n)
    xk = [ a + i*h for i in range(1,n)]
    return h/2 * (func(a) + 2 * sum_fun_xk(xk, func)+ func(b))
    
def base_peak(peaklist):
    '''
    return the maximum peak of a given peak list
    @peaklist, narray (mz, intensity)
    @return: mz, intensity pair with max intensity
    '''
    return max(sorted(peaklist)[:-1, :])

def cal_centroid(data):
    

    """
    cal_centroid function was adapted from pymzml(https://pymzml.readthedocs.io/en/latest/_modules/pymzml/spec.html)
    
    Follows are original discription from pyzmzl as quoted：
    Perform a Gauss fit to centroid the peaks for the property
    centroided_peaks.

    Returns:
        centroided_peaks (list): list of centroided m/z, i tuples
        
    """
    #specify m/z and intensity
    intens = [ i[1] for i in data]
    mzs = [ i[0] for i in data]
    tmp = []
    
    for pos, i in enumerate(intens[:-1]):
        if pos <= 1:
            continue
        if 0 < intens[pos - 1] < i > intens[pos + 1] > 0:
            
            #to pick peaks from profiles
            x1 = float(mzs[pos - 1])
            y1 = float(intens[pos - 1])
            x2 = float(mzs[pos])
            y2 = float(intens[pos])
            x3 = float(mzs[pos + 1])
            y3 = float(intens[pos + 1])
            
            #to avoid severe bias to gaussian dist with 10time threshold
            if x2 - x1 > (x3 - x2) * 10 or (x2 - x1) * 10 < x3 - x2:
                # no gauss fit if distance between mz values is too large
                continue
#            if y3 == y1:
#                # i.e. a reprofiledSpec
#                x1  = mzs[pos-5]
#                y1  = intens[pos-5]
#                x3  = mzs[pos+7]
#                y3  = intens[pos+7]

            try:
                
                # to calculate the normal distribution coefficients miu and sigma 
                double_log = math.log(y2 / y1) / math.log(y3 / y1)
                mue = (double_log * (x1 * x1 - x3 * x3) - x1 * x1 + x2 * x2) / (2 * (x2 - x1) - 2 * double_log * (x3 - x1))
                
                c_squarred = (x2 * x2 - x1 * x1 - 2 * x2 * mue + 2 * x1 * mue) / (2 * math.log(y1 / y2))

                # return sumed intensity of peak-0.4Da - peak +0.4Da 
                A = cal_intensity(data, mue)
                #A = y1 * math.exp((x1 - mue) * (x1 - mue) / (2 * c_squarred))
            except ZeroDivisionError:
                continue
            
            tmp.append((mue, A))
            
    return tmp

def cal_intensity(profile, p, m = 'full'):
    
    '''
    to calculate spectrum based peak intensity
    profile: curve
    p: peaks
    m: methods, and full equals to a 1 Da interval to be sumed up
    '''
    sumed_peak_intensity = 0.0
    #sum the intensity between -0.4 < target mz < 0.4Da
    for point in profile:
        if abs(point[0] - p) < 0.4:
            sumed_peak_intensity += point[1]
            
    return sumed_peak_intensity
            
            
        
    

def cal_noise_level(spect_lst, mode = 'median'):
    '''
    simply calculate noise level by mean or median intensity for given profiled spect_lst
    spect_lst: profiled
    '''    
    if mode == 'mean':
        return 0.4*np.mean([i[1] for i in spect_lst])
    if mode == 'median':
        return 0.4*np.median([i[1] for i in spect_lst])
    


def remove_noise(peak_lst, noise_level, snr):
    
    if noise_level == 0:
        return peak_lst
        
    noise_removed_peaks = [(p, i) for p, i in peak_lst if i / noise_level >= snr]
    return noise_removed_peaks

def cal_profile():
    pass


def save_res(result_buffer, fin, fout):
    '''
    create directory and save results 
    @result_buffer: final quantitation results
    @fin: filename of input data
    @fout
    
    @return:
    '''
    
    f_i = 1
    today = date.today()
    str_date = today.isoformat()
    
    if not fout:
        #ff0 was unspecified
        while True:
#            fout_name = fin.split('/')[-1].split('.')[0] +'_out_' + str(f_i) + '.txt'
            res_path = os.getcwd() + '\\' + fin.split('/')[-1].split('.')[0] + str_date + '.results.'+ str(f_i)
            if os.path.exists(res_path):
                f_i += 1
            else:
                break
    
    os.makedirs(res_path)
    
    csv_name_total = res_path+'\\' + fin.split('/')[-1].split('.')[0] +'_out_' + str(f_i) +  '_total' + '.csv'
    result_buffer.to_csv(csv_name_total, index = False, header = True)

