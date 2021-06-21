# -*- coding: utf-8 -*-
"""
Created on Sat Dec 29 13:02:52 2018

@author: hjm

# Copyright 2019 Jiangming Huang
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

"""

#from glyc_quant_main import *
from Tkinter import *
import ttk
import tkFileDialog as tf
from gquant_main_calibr_py2 import *
import os
import re

global ms_data_files


#Main Panel
root = Tk()
root.title("gQuant-Data processing tool for glycan stable isotope labeling")
root.geometry('550x650')

#Section 1 for data loading and database setting
g1 = LabelFrame(root, text="File IO setting", padx=5, pady=5)
g1.place( x = 10, y = 25, width=520, height=120)

# input settings, accept mass data output as plaintext or comma seperated variables (CSVs)
l1 = Label(g1, text="Select MS data: ")
l1.place(x = 2, y = 10)
ms_data_text = StringVar()
ms_data = Entry(g1, textvariable = ms_data_text, width = 45)
ms_data_text.set("")
ms_data.place(x = 100, y = 12)

# simply check if valid
def is_valid_data(fnames):
    '''
    '''
    pattern1 = r'.+/.txt'
    pattern2 = r'.+/.csv'
    result1 = re.findall(pattern1, fnames)
    result2 = re.findall(pattern2, fnames)

    if result1 == [] and result2 == []:
        return False
    elif result1 == []:
        return result2
    else:
        return result1
        
def data_set():
    global ms_data_files
    ms_data_files = tf.askopenfilenames(initialdir = os.getcwd(),title = "Select MS Data files",filetypes = (("MSMS files","*.txt"),("all files","*.*")))
    ms_data_text.set(ms_data_files)

b1 = Button(g1, text = "Browse", command = data_set)
b1.place(x = 440, y = 8)

l2 = Label(g1, text="Output folder: ")
l2.place(x = 2, y = 50)
output_file_name = StringVar()
output_file = Entry(g1, textvariable = output_file_name, width = 45)
output_file_name.set("")
output_file.place(x = 100, y = 52)

# To specify output directory
def output_set():
    output = tf.askopenfilename(initialdir = os.getcwd(),title = "Select Output Directory",filetypes = (("Comma Separated file","*.csv"),("all files","*.*")))
    output_file_name.set(output)
    
b2 = Button(g1, text = "Browse", command = output_set)
b2.place(x = 440, y = 48)

#Section 2 for general mass setting including exported mass spec data type()
g2 = LabelFrame(root, text="General Mass Setting", padx=5, pady=5)
g2.place( x = 10, y = 150, width=520, height=250)

l3_1 = Label(g2, text="MS mode:")
l3_1.place(x = 2, y = 10)
ms_mode = StringVar()
d_c1 = ttk.Combobox(g2, width = 16, textvariable = ms_mode, state = 'readonly')
d_c1['value'] = ('positive', 'negative')
d_c1.current(0)
d_c1.place(x = 100, y = 12)

l3_2 = Label(g2, text="MS tolerance:")
l3_2.place(x = 255, y = 10)
d_value2 = DoubleVar()
dv2 = Entry(g2, textvariable = d_value2, width = 10)
d_value2.set(100.0)
dv2.place(x = 350, y = 12)
unit2 = StringVar()
d_c2 = ttk.Combobox(g2, width = 5, textvariable = unit2, state = 'readonly')
# =currently only support ppm
d_c2['value'] = ('ppm')
d_c2.current(0)
d_c2.place(x =430, y = 12)

# data type setting for further processing of data point(profile) or direct count(centroid), as supported by vendor software.
l4_1 = Label(g2, text="Data type:")
l4_1.place(x = 2, y = 50)
unit3 = StringVar()
d_c3 = ttk.Combobox(g2, width = 16, textvariable = unit3, state = 'readonly',)
d_c3['value'] = ('profile (curve)', 'centroid (peak)')
d_c3.current(1)
d_c3.place(x = 100, y = 52)


l4_1_2 = Label(g2, text="Delimiter:")
l4_1_2.place(x = 255, y = 50)
delimi = StringVar()
d_l1 = ttk.Combobox(g2, width = 16, textvariable = delimi, state = 'readonly',)
d_l1['value'] = ('tab (\\t)', 'comma (,)', 'space ( )', 'semicolon (;)')
d_l1.current(0)
d_l1.place(x = 350, y = 52)


l4_1_3_1 = Label(g2, text="mz_calibr:")
l4_1_3_1.place(x = 2, y = 88)
mz_calibr = DoubleVar()
mz_c = Entry(g2, textvariable = mz_calibr, width = 10)
mz_calibr.set(0.0)
mz_c.place(x = 100, y = 90)
l4_1_3_2 = Label(g2, text="Da")
l4_1_3_2.place(x = 180, y = 88)


l4_2_1 = Label(g2, text="Min intensity:")
l4_2_1.place(x = 255, y = 127)
quant_percent_threshold = DoubleVar()
qpt = Entry(g2, textvariable = quant_percent_threshold, width = 10)
quant_percent_threshold.set(0.0)
qpt.place(x = 350, y = 127)
l4_2_2 = Label(g2, text="%")
l4_2_2.place(x = 420, y = 125)

l4_3 = Label(g2, text="S/N value:")
l4_3.place(x = 2, y = 125)
signal_to_noise = DoubleVar()
snv = Entry(g2, textvariable = signal_to_noise, width = 10)
signal_to_noise.set(20.0)
snv.place(x = 100, y = 127)

l4_4 = Label(g2, text="Max charges:")
l4_4.place(x = 255, y = 88)
max_charge = IntVar()
max_charge_value = Spinbox(g2,from_=1,to=8,increment=1, textvariable = max_charge, width = 10)
max_charge.set(1)
max_charge_value.place(x = 350, y = 90)

#***************************************************
# include common adducts
#
#***************************************************

l4_5 = Label(g2, text="Adduct(s):")
l4_5.place(x = 2, y = 160)

s2_1 = IntVar()
cb2_1 = Checkbutton(g2, text="H+", variable = s2_1)
cb2_1.place(x = 100, y = 160)
cb2_1.select()

s2_2 = IntVar()
cb2_2 = Checkbutton(g2, text="Na+", variable = s2_2)
cb2_2.place(x = 200, y = 160)
cb2_2.select()

s2_3 = IntVar()
cb2_3 = Checkbutton(g2, text="K+", variable = s2_3)
cb2_3.place(x = 300, y = 160)
cb2_3.select()

s2_4 = IntVar()
cb2_4 = Checkbutton(g2, text="Li+", variable = s2_4)
cb2_4.place(x = 400, y = 160)


l4_6 = Label(g2, text="New adduct")
new_adduct_name = StringVar()
new_adduct_name_input = Entry(g2, textvariable = new_adduct_name, width = 4)
new_adduct_name.set("X")


l4_7 = Label(g2, text="Adduct mass")
new_adduct_mass = DoubleVar()
new_adduct_mass_input = Entry(g2, textvariable = new_adduct_mass, width = 8)
new_adduct_mass.set("0.0000")

def new_adduct():
    '''
    to add one user specify adduct ions with name and mass information
    @effect click to show and hide
    '''
    
    if s2_5.get():
        # checked, expect s2_5 == 1. show function and entry
        l4_6.place(x = 200, y = 192)
        new_adduct_name_input.place(x = 280, y = 192)
        l4_7.place(x = 335, y = 190)
        new_adduct_mass_input.place(x = 420, y = 192)
        
    else:
        # unchecked, will hide automatically
    
        l4_6.place_forget()
        new_adduct_name_input.place_forget()
        l4_7.place_forget()
        new_adduct_mass_input.place_forget()
        
   
s2_5 = IntVar()
cb2_5 = Checkbutton(g2, text="Other", command = new_adduct,variable = s2_5)
cb2_5.place(x = 100, y = 190)


#Section 3 contain simple setting for glycan de novo sequencing
g3 = LabelFrame(root, text="Glycan parameters", padx=5, pady=5)
g3.place( x = 10, y = 405, width=520, height=130)

l5_1 = Label(g3, text="Derivatization \nglycan mass:")
l5_1.place(x = 2, y = 4)
#label_mass_value = DoubleVar()
#lv = Entry(g3, textvariable = label_mass_value, width = 9)
#label_mass_value.set(0.0)
#lv.place(x = 100, y = 12)

l5_1_1 = Label(g3, text="Derivatization \ntype & mass:")

label_name = StringVar()
ln = Entry(g3, textvariable = label_name, width = 7)
label_name.set('PFBHA')
label_mass_value = DoubleVar()
lv = Entry(g3, textvariable = label_mass_value, width = 10)
label_mass_value.set(195.01019)

label_unit = Label(g3, text="Da")

label_mass_unit = StringVar()
l_c = ttk.Combobox(g3, width = 16, textvariable = label_mass_unit, state = 'readonly')
l_c['value'] = ('FreeEnd(+0.0 Da)', 'RedEnd(+2.0 Da)','2AB(+119.06 Da)', '2AA(+120.06 Da)')
l_c.current(0)
l_c.place(x = 100, y = 14)


def custom_derivatization():
    
    # hide extra widget
    l5_1.place_forget()
    l_c.place_forget()
    l5_1_1.place(x = 2, y = 4)
    ln.place(x = 100, y = 14)
    lv.place(x = 155, y = 14)
    label_unit.place(x = 220, y = 14)
    # Switch Function
    b6.config(text="Default Derivatization", command = revert_derivatization)

def revert_derivatization():
    '''
    
    #partially adapted from https://www.cnblogs.com/brightyuxl/archive/2018/10/22/9832248.html
    '''
    l5_1_1.place_forget()
    ln.place_forget()
    lv.place_forget()
    label_unit.place_forget()
    l5_1.place(x = 2, y = 4)
    l_c.place(x = 100, y = 12)
    

    # buttom switch
    b6.config(text='Custom Derivatization', command = custom_derivatization)

b6 = Button(g3, text="Custom Derivatization", command = custom_derivatization, state = 'normal')
b6.place(x = 245, y = 10)   


is_full_deri = IntVar()
cb_is_full_deri = Checkbutton(g3, text="Fully derivated", variable = is_full_deri)
cb_is_full_deri.place(x = 390, y = 10)

#l5_2 = Label(g3, text="Delta mass for\nrelative quant:")
#l5_2.place(x = 250, y = 4)
#glycan_rel_quant_value = DoubleVar()
#dv = Entry(g3, textvariable = glycan_rel_quant_value, width = 9)
#glycan_rel_quant_value.set(2.0)
#dv.place(x = 350, y = 12)
#delta_mass_value = StringVar()
#d_c = ttk.Combobox(g3, width = 5, textvariable = delta_mass_value, state = 'readonly')
#d_c['value'] = ('Da')
#d_c.current(0)
#d_c.place(x = 425, y = 12)

l5_2 = Label(g3, text="Delta mass for\nrelative quant:")
l5_2.place(x = 2, y = 54)
glycan_rel_quant_value = DoubleVar()
dv = Entry(g3, textvariable = glycan_rel_quant_value, width = 10)
glycan_rel_quant_value.set(2.0151)
dv.place(x = 100, y = 62)
delta_mass_value = StringVar()
d_c = ttk.Combobox(g3, width = 4, textvariable = delta_mass_value, state = 'readonly')
d_c['value'] = ('Da')
d_c.current(0)
d_c.place(x = 185, y = 62)

#l5_3 = Label(g3, text="Largest glycan \nmass:")
#l5_3.place(x = 2, y = 60)
#g_mass_text = DoubleVar()
#g_mass = Entry(g3,textvariable = g_mass_text, width = 18)
##Entry(g3, textvariable = g_node_text, width = 10)
#g_mass_text.set(4500.0)
#g_mass.place(x = 100, y = 62)


l5_4 = Label(g3, text="GlycanDatabase:")
l5_4.place(x = 250, y = 60)
gdb = StringVar()
gdb_cbb = ttk.Combobox(g3, width = 16, textvariable = gdb, state = 'readonly')
gdb_cbb['value'] = ('Human', 'Mammalian')
gdb_cbb.current(0)
gdb_cbb.place(x =350, y = 62)


#Simple indicator for running status can be seen here
l6 = Label(root, text="Running state")
l6.place(x = 10, y = 545)
d_state = StringVar()
ds = Entry(root, textvariable = d_state, bg = 'light grey', state = 'disabled', width = 25)
d_state.set("IDLE")
ds.place(x = 110, y = 547)

global res 

def match_n_quant():
    '''
    Initialization
    '''
    if ms_data_text.get() == '' or ms_data_text.get() == "Select correct data file!":
        ms_data_text.set("Select correct data file!")
        ms_data.config(fg = 'red')
        return None
    if is_valid_data(ms_data_text.get()):
        ms_data.config(fg = 'black')
    # Initial file input and output

    global ms_data_files
    if not ms_data_files or ms_data_text == "Select correct data file!":
        ms_data_text.set("Select correct data file!")
        ms_data.config(fg = 'red')
        return None
    else:
        ms_data.config(fg = 'black')
    
    file_out = output_file_name.get()
    
    
    # initial general mass setting
    ms_mode = d_c1.get()
    mz_tol = d_value2.get()
    mz_tol_unit = d_c2.get()
    data_type = d_c3.get()
    
    dl_type = '\t'
    if delimi.get() == 'tab (\\t)':
        dl_type = '\t'
    if delimi.get() == 'comma (,)':
        dl_type= ','
    if delimi.get() == 'space ( )':
        dl_type = ' '
    if delimi.get() == 'semicolon (;)':
        dl_type = ';'
    
    peak_intensity_threshold = quant_percent_threshold.get()    
    snr = signal_to_noise.get()
    
    consider_max_charge = max_charge.get()
    
    check_adducts = {}
    is_H = s2_1.get()
    is_Na = s2_2.get()
    is_K = s2_3.get()
    is_Li = s2_4.get()

    is_full_deri_st = is_full_deri.get()
    
    # add to list
    if is_H:
        check_adducts.update({'H':1.0078})
    if is_Na:
        check_adducts.update({'Na':22.9898})
    if is_K:
        check_adducts.update({'K':38.9637})
    if is_Li:
        check_adducts.update({'Li':6.941})
    if s2_5.get():
        X_name = new_adduct_name.get()
        adduct_X = new_adduct_mass.get()
        if X_name and adduct_X:
            check_adducts.update({X_name:float(adduct_X)})
    # no mass set to new adduct        
    
    # initiate glycan match and quantitation parameters
    glycan_deravatization_mass =label_mass_value.get()
    glycan_deravatization_name = label_name.get()
    
    check_deriva = {}
    if l_c.get() == 'FreeEnd(+0.0 Da)' :
        check_deriva = {}
    if l_c.get() == 'RedEnd(+2.0 Da)':
        check_deriva= {'Reduce End':2.0151}
    if l_c.get() == '2AB(+119.06 Da)':
        check_deriva = {'2AB': 119.06 }
    if l_c.get() == '2AA(+120.06 Da)':
        check_deriva = { '2AA':120.06}

    if b6["text"] == "Default Derivatization" and glycan_deravatization_name and glycan_deravatization_mass:
        check_deriva = {glycan_deravatization_name:glycan_deravatization_mass}    
        
    
    glycan_rel_quant_delta_mass = glycan_rel_quant_value.get()
    glycan_database_type =  gdb_cbb.get()
    
    
    print "*********************************************************************"
    print "gQuant runed with following parameters:"
    print "====================================================================="
    print "mass spect related parameters:"     
    print "Input file: %s;" 
    print ms_data_files
    print "Output directory: %s;" %file_out
    print "mass spectrometry mode: %s" %ms_mode
    print "m/z tolerance: %s %s" %(mz_tol,mz_tol_unit)
    print "Exported data form: %s" %data_type
    print "Peak intensity threshold to quant: %s%%" %peak_intensity_threshold
    print "signal to noise ratio threshold: %s" %snr
    print "max charge considered: %s" %consider_max_charge
    print "adduct list(s): "
    print check_adducts
    print "====================================================================="
    print "glycan match and quantitation parameters:"    
    print "glycan derivatization mass: %s" %glycan_deravatization_mass
    print "glycan relative quantitation delta mass: %s" %glycan_rel_quant_delta_mass
    print "max charge consier: %s" %consider_max_charge
    print "glycan database selected: %s" %glycan_database_type
    print "*********************************************************************"

    
    quant_it(ms_data_files, file_out, data_type, dl_type, snr, consider_max_charge, glycan_database_type,check_deriva, check_adducts,is_full_deri_st,mz_tol,glycan_rel_quant_delta_mass)
    
    global res
    res = os.getcwd()
#    d_state.set("Task finished in %s min!" % str(t_time))
    ds.config(bg = 'green')
    b5.config(state = 'normal')
        
def cancel_n_exit():
    root.destroy()


def open_res_path():

    if res and os.path.exists(res):
        if os.path.isfile(res):
            import win32process
            try:
                win32process.CreateProcess(res, '',None , None , 0 ,win32process. CREATE_NO_WINDOW , None , None ,win32process.STARTUPINFO())
            except Exception, e:
                print(e)
        else:
            os.startfile(str(res))
    
    else:
        print('directory not found')

b3 = Button(root, text="Quant It!", command = match_n_quant, bg = 'wheat')
b3.place(x = 390, y = 545)
b4 = Button(root, text=" Exit! ", command = cancel_n_exit, bg = 'Grey')
b4.place(x = 475, y = 545)
b5 = Button(root, text=" Open result folder  ", command = open_res_path, state = 'disabled')
b5.place(x = 390, y = 585)


def quant_it(param1, param2, param3, param4, param5, param6, param7,param8, param9,param10,param11,param12):
    try:
        if param3 == 'profile (curve)':
            # exported datapoints
            stat_data = pd.DataFrame(columns = ['File', 'data_point', 'num_peak_with_envelope','total_envelope', 'cost_time'])
            
            for files in param1:
                #auto process
                if files.split('.')[-1] == 'txt':
                    start_time = time.time()
                    a3_profile = load_ms_data(files,sep = param5)
                    nl = cal_noise_level(a3_profile)
                    peak_a = cal_centroid(a3_profile)
                    peak_remove_noise = remove_noise(peak_a, nl, param5)
                    peak_remove_noise_calibr = [(a[0]+param6_c, a[1]) for a in peak_remove_noise]
                    peaksn, full_envelopen = deisotope_profile(peak_remove_noise_calibr)
                    
                    
                    gdb2 = load_glycan_db(param7)
                    
                    gdb_with_adduct = glycan_db_with_mods(initial_gdb = gdb2, derivatization = param8, adducts = param9, is_full_de = param10)
                    
                    mg7 = match_glycan_ppm(peaksn,gdb_with_adduct, mz_tol = param11, delta_mass=param12 )
                    gg5 = gquant(mg7,full_envelopen,gdb_with_adduct,delta_mass =param12 )
                    
                    num_point = len(peak_a)
                    num_peak = len(peaksn)
                    num_env = len(full_envelopen)
    
                    if param2:
                        while True:
                            if os.path.exists(param2):
                                param2 = param2.split('.')[0] + '01'+'.csv'
                            else:
                                file_out2 = param2
                                break
                    else:
                        fout_name_pre = files.split('.')[0] + 'out'
                        while True:
                            if os.path.exists(fout_name_pre + '.csv'):
                                fout_name_pre = fout_name_pre + '01'
                            else:
                                file_out2 = fout_name_pre+'.csv'
                                break

                    gg5.to_csv(file_out2, index = False)
                    
                    end_time = time.time()
                    total_time = round((end_time - start_time) / 60.0, 2)
                    stat_data = stat_data.append(pd.DataFrame({'File':[files], "data_point":[num_point],'num_peak_with_envelope':[num_peak], 'total_envelope':[num_env],'cost_time':[total_time]}), ignore_index=True)
                    total_time = 0
                    num_point = 0
                    num_peak = 0
                    num_env = 0            
                    print 'Done'
            
            fstat = 'stat.csv'
            while True:
                if os.path.exists(fstat):
                    fstat = fstat.split('.')[0] + '01'+'.csv'
                else:
                    break                
            stat_data.to_csv(fstat, index = False)
    
        elif param3 == 'centroid (peak)':
            # exported peaks
            stat_data = pd.DataFrame(columns = ['File', 'data_point', 'num_peak_with_envelope','total_envelope', 'cost_time'])

            print param1
            print param2
            print param3
            print param4
            print param5
            print param6
            print param7
            print param8
            print param9
            print param10
            print param11
            print param12

                        
            for files in param1:
                print files
                #auto process
                if files.split('.')[-1] == 'txt':
                    start_time = time.time()
            
                    a3 = load_ms_data(files, sep = param4)
                    
                    peaks, full_envelope = deisotope_profile(a3,param6)
                    
                    gdb2 = load_glycan_db(param7)
    
                    gdb_with_adduct = glycan_db_with_mods(initial_gdb = gdb2, derivatization = param8, adducts = param9, is_full_de = param10)
                    
                    print '1'
                    mg7 = match_glycan_ppm(peaks,gdb_with_adduct, float(param11), float(param12))
                    #mg7 = match_glycan_ppm(peaks,gdb_with_adduct, mz_tol = 150.0, delta_mass= 2.0151, mz_calibr = 0.0)
                    print '2'
                    gg5 = gquant(mg7,full_envelope,gdb_with_adduct, a3, delta_mass= float(param12))
                    
                    num_point = len(a3)
                    num_peak = len(peaks)
                    num_env = len(full_envelope)
                            
                    if param2:
                        while True:
                            if os.path.exists(param2):
                                param2 = param2.split('.')[0] + '01'+'.csv'
                            else:
                                file_out2 = param2
                                break
                    else:
                        fout_name_pre = files.split('.')[0] + 'out'
                        while True:
                            if os.path.exists(fout_name_pre + '.csv'):
                                fout_name_pre = fout_name_pre + '01'
                            else:
                                file_out2 = fout_name_pre+'.csv'
                                break

                    gg5.to_csv(file_out2, index = False)

                    end_time = time.time()
                    total_time = round((end_time - start_time) / 60.0, 2)
                    stat_data = stat_data.append(pd.DataFrame({'File':[files], "data_point":[num_point],'num_peak_with_envelope':[num_peak], 'total_envelope':[num_env],'cost_time':[total_time]}), ignore_index=True)
                    print files + 'cost' + str(total_time)      
                    total_time = 0
                    num_point = 0
                    num_peak = 0
                    num_env = 0      
                    print 'Done'
                    
            fstat = 'stat.csv'
            while True:
                if os.path.exists(fstat):
                    fstat = fstat.split('.')[0] + '01'+'.csv'
                else:
                    break                
            stat_data.to_csv(fstat, index = False)
            
    except Exception as e:
        print 'check carefully', e



root.mainloop()
