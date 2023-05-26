import sys
sys.path.insert(1,"./")
import numpy as np
import math
from rocketcea.cea_obj_w_units import CEA_Obj
import Toolbox.RListGenerator as RListGenerator
import Toolbox.RocketCEAAssister as RA
import Toolbox.IsentropicEquations as IE
import Toolbox.RocketEquation as RE
import difflib
import re as regex
from rocketprops.rocket_prop import get_prop
import Toolbox.Constant as const
import matplotlib.pyplot as plt
from matplotlib import cm
import Analysis.FirstOrderCalcs as FAC
import Components.ThrustChamber as ThrustChamber
import Components.CoolingSystem as CS
import Components.StructuralApproximation as SA
import Toolbox.CADAssistant as CAD
from scipy.optimize import minimize_scalar
import scipy.optimize
import Main
import GUIFuncs
import os
import PySimpleGUI as sg
import json
import pandas as pd

def listSelectPopup(listnames):
    layout = [
        [sg.Text("Select A List to Plot Against (list of x values)")],
        [sg.Listbox(listnames, enable_events=True, size=(40,20), key='SELECTED')],
        [sg.Button("Plot List Against it's Own Range", key = "NOLIST")]
    ]
    
    window = sg.Window('POPUP', layout).Finalize()
    
    while True:
        event, values = window.read()

        if event == sg.WINDOW_CLOSED:
            break
        elif event == 'SELECTED':
            break
        elif event == 'NOLIST':
            window.close()
            return None
        else:
            print('OVER')

    window.close()
    if values and values['SELECTED']:
        return values['SELECTED'][0]
def updateDisplayNames(maininputs):
            maindisplaynames = []
            for name in maininputs.keys():
                try:
                    try:
                        try:
                            maindisplaynames.append(name + " =   " + str(int(10000*maininputs[name])/10000) + "   " + unitdict[name])
                        except:
                            maindisplaynames.append(name + " =   " + str(maininputs[name]) + "   " + unitdict[name])
                    except:
                        try:
                            maindisplaynames.append(name + " =   " + str(int(10000*maininputs[name])/10000) + "   " + 'no units yet')
                        except:
                            maindisplaynames.append( name + " =   " + str(maininputs[name]) + "   " + 'no units yet')
                except:
                    maindisplaynames.append(name + " =   Null   " + 'no units yet')
            return maindisplaynames
def get_input_layout(canDelete = False, canAdd = True):
    if canDelete:
        return [
            [sg.Text('input'), sg.InputText(default_text='0', key='input')],
            [sg.Button('Change Input', key='change')],
            [sg.Button('Delete', key='delete')],
            [sg.Text('Optional Name of New Variable'), sg.InputText(default_text='0', key='addname')],
            [sg.Button('Add New Variable', key='add')]]
    else:
        return [
            [sg.Text('input'), sg.InputText(default_text='0', key='input')],
            [sg.Button('Change Input', key='change')],
            [sg.Text('Optional Name of New Variable'), sg.InputText(default_text='0', key='addname')],
            [sg.Button('Add New Variable', key='add')]]
                    
        
def run_selected_function(function_name, inputs):
    if function_name == "main":
        return Main.main(inputs["args"], inputs["config_title"], output=inputs["output"])
    if function_name == "Everything_MassFixed": #needs params['drymass']
        return GUIFuncs.Everything_MassFixed(inputs['args'], inputs['chanelgeometry'])
    if function_name == "FirstOrderCalcs":
        return GUIFuncs.FirstOrderSolver(inputs)

def open_function_window(function_name):
    if function_name == 'main':

        
        try:
            maininputs = {'thrust' : params["thrust"],
                            'fuelname' : params["fuelname"],
                            'oxname' : params["oxname"],
                            'pc' : params["pc"],
                            'rc' : params['rc'],
                            'dp' : params['dp'],
                            'time' : params['time'],
                            'lstar' : params['lstar'],
                            'throat_radius_curvature' : params['throat_radius_curvature'],
                            'TWR' : params['TWR']}
        except:
            maininputs = {'thrust' : 1000,
                            'fuelname' : "null",
                            'oxname' : 'null',
                            'pc' : 0}
        maindisplaynames = updateDisplayNames(maininputs)
        
        config_column = [
            [sg.Text("\/ Available Arguments \/")],
            [
                sg.Listbox(
                    values=displaynames, enable_events=True, size=(40, 30), key="-AVAIL-"
                )
            ],
            [
                sg.Text("Config Name"),
                sg.In(size=(15, 1), enable_events=True, key="-CONFIG NAME-")
            ]
            
        ]

        # For now will only show the name of the file that was chosen
        function_caller_column = [
            [sg.Text("\/ Selected Arguments \/")],
            [
                sg.Listbox(
                    values=maindisplaynames, enable_events=True, size=(40, 20), key="-SELECTED-"
                )
            ],
            
            [
                sg.Text("Output"),
                sg.Combo(["True - Set Config Name!","False"], size = 20, default_value='False', k="-OUTPUT-", enable_events=True, readonly=True)
            ],
            [sg.Button("Run Function")]
                
        ]

        # ----- Full layout -----
        main_layout = [
            [
                sg.Column(config_column),
                sg.Column(function_caller_column),
            ]
        ]
        window = sg.Window("Main Inputs", main_layout)

# Run the Event Loop
        while True:
            event, values = window.read()
            if event == "Exit" or event == sg.WIN_CLOSED:
                break
            if event == "-AVAIL-":
                paramname = values["-AVAIL-"][0].split(" ")[0]
                try:
                    maininputs[paramname] = params[paramname]
                except:
                    maininputs[paramname] = 'Null'
                maindisplaynames = updateDisplayNames(maininputs)
                window['-SELECTED-'].update(maindisplaynames)
            if event == "-SELECTED-":
                paramname = values["-SELECTED-"][0].split(" ")[0]
                inputs_layout = get_input_layout(canDelete=True)
                inputs_window = sg.Window(function_name + ' inputs', inputs_layout)
                inputs_event, inputs_values = inputs_window.read()
                while inputs_event != 'change' and inputs_event != 'delete' and inputs_event != 'add':
                    inputs_event, inputs_values = inputs_window.read()
                    if inputs_event == "Exit" or inputs_event == sg.WIN_CLOSED:
                        break
                inputs_window.close()
                if inputs_event == 'change':
                    try:
                        try:
                            maininputs[paramname] = float(inputs_values["input"])
                        except:
                            maininputs[paramname] = inputs_values["input"]
                    except:
                        maininputs[paramname] = 'Null'
                    maindisplaynames = updateDisplayNames(maininputs)
                    window['-SELECTED-'].update(maindisplaynames)
                if inputs_event == 'delete':
                    del maininputs[paramname]
                    maindisplaynames = updateDisplayNames(maininputs)
                    window['-SELECTED-'].update(maindisplaynames)
                if inputs_event == 'add':
                    try:
                        try:
                            maininputs[inputs_values['addname']]= float(inputs_values["input"])
                        except:
                            maininputs[inputs_values['addname']]= inputs_values["input"]
                    except:
                        pass
                    maindisplaynames = updateDisplayNames(maininputs)
                    window['-SELECTED-'].update(maindisplaynames)
            if event == "Run Function":
                inputs={}
                inputs["args"] = maininputs
                print(inputs)
                try:
                    inputs["config_title"] = values['-CONFIG NAME-']
                    inputs["output"] = values['-OUTPUT-'] == "True - Set Config Name!"
                except:
                    inputs["config_title"] = "YOU FORGOT A TITLE IDIOT"
                    inputs["output"] = False
                returnparams = run_selected_function("main", inputs)
                break
        window.close()
        return returnparams
    elif function_name == "RocketEquation_MassFixed":
        layout = [
            [sg.Text("Input Dry mass (kg)")],
            [sg.InputText(default_text='200', key='drymass')],
            #[sg.Text("Input Alittiude Target (m)")],
            #[sg.InputText(default_text='100000', key='target')],
            [sg.Button("Run")]
        ]
        
        window = sg.Window('POPUP', layout).Finalize()
        
        while True:
            event, values = window.read()

            if event == sg.WINDOW_CLOSED:
                break
            elif event == 'Run':
                break

        window.close()
        params["drymass"] = float(values["drymass"])
        #params["target"] = float(values["target"])
        return GUIFuncs.RocketEquation_MassFixed(params)
    elif function_name == "Everything_MassFixed":
        pass
    elif function_name == 'ChanelSolver':
        chanelGeometry = {'chlist': None, 'twlist': None, 'nlist': None, 'ewlist': None, 'helicitylist': None, 'c2l': None}
        print(chanelGeometry)
        
        layout = [
            [sg.Text('LoadDefaults'), sg.Button("DEFAULTS")], # Thinking there should be a standard set, perhaps from a cvs file
            [sg.Text('Input chlist (MATH EQUATION)'), sg.Input(key='chlist'), sg.Checkbox('Constant', default=False, key='ConstantChlist')], 
            [sg.Text('Input twlist (MATH EQUATION)'), sg.Input(key='twlist')],
            [sg.Text('Input nlist (MATH EQUATION)'), sg.Input(key='nlist')],
            [sg.Text('Input ewlist (MATH EQUATION)'), sg.Input(key='ewlist')],
            [sg.Text('Input helicitylist (MATH EQUATION)'), sg.Input(key='helicitylist')],
            [sg.Text('Input dxlist (MATH EQUATION)'), sg.Input(key='dxlist')],
            [sg.Text('Input c2l (MATH EQUATION)'), sg.Input(key='c2l')],
            [sg.Button("Run")]
        ]
        
        window = sg.Window('ChanelGeom', layout)
        df = pd.read_csv(r'Defaults/DefaultChanelGeom.csv')
        while True:
            event, values = window.read()

            if event == sg.WINDOW_CLOSED:
                break
            elif event == 'Run':
                if values['chlist']:
                    if values['ConstantChlist'] == True:
                        l = 200 #default len 200 good?
                        chanelGeometry['chlist'] = [values['chlist']]*l
                        print("set chlist with default")
                        print(chanelGeometry['chlist'])
                    else:
                        # function for chlist, probably best if it is called from GUI Funcs w/ parameter (input or inputs(list format), 0)
                        print("set chlist with function")
                else: 
                    print("chlist set to default")
                    chanelGeometry['chlist'] = list(df['chlist'])
                
                if values['twlist']:
                    # function for twlist, probably best if it is called from GUI Funcs w/ parameter (input or inputs, 1)
                    print("set twlist")
                else: 
                    print("twlist set to default")
                    chanelGeometry['twlist'] = list(df['twlist'])
                
                if values['nlist']:
                    # function for nlist, probably best if it is called from GUI Funcs w/ parameter (input or inputs, 2)
                    print("set nlist")
                else: 
                    print("nlist set to default")
                    chanelGeometry['nlist'] = list(df['nlist'])
                
                if values['ewlist']:
                    # function for ewlist, probably best if it is called from GUI Funcs w/ parameter (input or inputs, 3)
                    print("set ewlist")
                else: 
                    print("ewlist set to default")
                    chanelGeometry['ewlist'] = list(df['ewlist'])
                
                if values['helicitylist']:
                    # function for helicitylist, probably best if it is called from GUI Funcs w/ parameter (input or inputs, 4)
                    print("set helicitylist")
                else: 
                    print("helicitylist set to default")
                    chanelGeometry['helicitylist'] = list(df['helicitylist'])
                
                if values['dxlist']:
                    # function for dxlist, probably best if it is called from GUI Funcs w/ parameter (input or inputs, 5)
                    print("set dxlist")
                else: 
                    print("dxlist set to default")
                    chanelGeometry['dxlist'] = list(df['dxlist'])
                
                if values['c2l']:
                    # function for c2l, probably best if it is called from GUI Funcs w/ parameter (input or inputs, 6)
                    print("set c2l")
                else: 
                    print("c2l set to default")
                    chanelGeometry['c2l'] = list(df['c2l'])
                window.close()
                
            elif event == 'DEFAULTS':
                chanelGeometry['chlist'] = list(df['chlist'])
                chanelGeometry['twlist'] = list(df['twlist'])
                chanelGeometry['nlist'] = list(df['nlist'])
                chanelGeometry['ewlist'] = list(df['ewlist'])
                chanelGeometry['helicitylist'] = list(df['helicitylist'])
                chanelGeometry['dxlist'] = list(df['dxlist'])
                chanelGeometry['c2l'] = list(df['c2l'])
                window.close()

        window.close()
        print(chanelGeometry)
        #Trivial to add whatever function now that layout is done
        GUIFuncs.ChanelSolver(params, chanelGeometry)


        


# First the window layout in 2 columns
workingfuncnames = ['main','coolingsystem','rocketequation',"Everything_MassFixed","FirstOrderCalcs","RocketEquation_MassFixed", "ChanelSolver"]
sg.theme('Dark Purple 3')
paramnames = [
    'thrust', #Newtons
    'time', #s
    'impulse', #N*s
    'rho_ox', #Kg/M^3
    'rho_fuel', #Kg/M^3
    'pc', #Pa, if cr is specified, this is Pressure at end of combustor
    'pinj',#Pa, only useful if you specify CR, otherwise assumed to be pc
    'pe', #Pa
    'g', #m/s^2
    'rm', #o/f by weight
    'phi', #ratio from stoich (1 is stoich, >1 is fuel rich)
    'at', # m^2, area of throat
    'rt', # m, radius of throat 
    'cr', # contraction ratio
    'rc', # m, combustion chamber radius
    'ac', # m^2, area combustion chamber
    'l_star', # m, volume cc/area throat
    'mol_weight', # kg/mol
    'gamma', # in cc
    'gamma_exit', # in exit
    'gamma_throat', # in throat
    'isp', # s
    'temp_c', # K, chamber temp
    'rg', # specific gas constant (SI units if what they are)
    'pr_throat',
    'rho_throat',
    'temp_e',
    'v_exit',
    'a_exit',
    'mach_exit',
    'temp_throat',
    'p_throat',
    'v_throat',
    'mdot',
    'mdot_ox',
    'mdot_fuel',
    'er',
    'cstar',
    'cf',
    'c_eff',
    'rho_av',
    'vc',
    'theta_con',
    'lc',
    'theta_div',
    'ln_conical',
    'ln_bell',
    'throat_radius_curvature',
    'ae',
    're',
    'nv',
    'nvstar',
    'nf',
    'nw',
    'fuelname',
    'oxname',
    'CEA',
    'pambient',
    'cf_efficiency', # Huzel and Huang page 16
    'isp_efficiency',
    'thetac', # converging section angle
    'thetai', # diverging angle at throat
    'thetae'] # diverign angle at exit
unitdict = {
    'thrust' : "N",
    'time' : "s",
    'impulse' : "N*s",
    'rho_ox' : 'kg/M^3',
    'rho_fuel': 'kg/M^3',
    'pc' : "Pa, if cr is specified, this is Pressure at end of combustor",
    'pinj' : "Pa, only useful if you specify CR, otherwise assumed to be pc",
    'pe' : 'Pa',
    'g' : 'm/s^2',
    'rm' : 'o/f by mass',
    'phi' : 'ratio from stoich (1 is stoich, >1 is fuel rich)',
    'at' :  'm^2, area of throat',
    'rt' :  'm, radius of throat' ,
    'cr' :  'contraction area ratio',
    'rc' :  'm, combustion chamber radius',
    'ac' :  'm^2, area combustion chamber',
    'l_star' :  'm, volume cc/area throat',
    'mol_weight' :  'kg/mol',
    'gamma' :  'in cc',
    'gamma_exit' :  'in exit',
    'gamma_throat' :  'in throat',
    'isp' :  's',
    'temp_c' :  'K, chamber temp',
    'rg' :  'specific gas constant (SI units whatever they are)',
    'pr_throat' : "dimensionless?",
    'rho_throat' : 'kg/M^3',
    'temp_e' : 'K',
    'v_exit' : 'm/s',
    'a_exit' : 'm/s, speed of sound',
    'mach_exit' : 'Mach',
    'temp_throat' : 'K',
    'p_throat' : 'Pa',
    'v_throat' : 'm/s',
    'mdot' : 'kg/s',
    'mdot_ox' : 'kg/s',
    'mdot_fuel' : 'kg/s',
    'er' : 'Expansion Area Ratio',
    'cstar' : 'm/s',
    'cf' : 'dimensionless?',
    'c_eff' : 'm/s',
    'rho_av' : 'kg/s',
    'vc' : 'm^3',
    'theta_con' : 'radians' ,
    'lc' : 'm, combustor length',
    'theta_div' : 'radians',
    'ln_conical' : 'm, 15 degree conical',
    'ln_bell' : 'm',
    'throat_radius_curvature' : 'm',
    'ae' : 'm^2',
    're' : 'm^2',
    'nv' : 'dimenionless',
    'nvstar' : 'dimenionless',
    'nf' : 'dimenionless',
    'nw' : 'dimenionless',
    'fuelname' : '',
    'oxname' : '',
    'CEA' : '<- this is an object',
    'pambient' : 'Pa',
    'cf_efficiency' :  'Huzel and Huang page 16',
    'isp_efficiency' :  'dimensionless',
    'thetac' :  'converging section angle',
    'thetai' :  'diverging angle at throat',
    'thetae' : 'diverign angle at exit',# diverign angle at exit
    'kin_visc_fuel' : "Pa s m3/kg",
    'kin_visc_ox' : "Pa s m3/kg",
    'dyn_visc_fuel' : "Pa s",
    'dyn_visc_ox' : "Pa s"}
params={
        'thrust': 4000 * const.lbToN,  # Newtons
        'time': 40,  # s
        # 'rho_ox' : 1141, #Kg/M^3
        # 'rho_fuel' : 842,
        'pc': 300 * const.psiToPa,
        'pe': 14.7 * const.psiToPa,
       # 'phi':1,
        'cr': None,
        'TWR' : 4,
        'lstar': 1.24,
        'fuelname': 'Ethanol_75',
        'oxname': 'N2O',
        'throat_radius_curvature': .0254 *2,
        'dp': 150 * const.psiToPa,
        'rc' : .11,
        'thetac' : (35*math.pi/180)}
storedlists = {}
displaynames = updateDisplayNames(params)
listdisplaynames = list(storedlists.keys())
config_column = [
    [
        sg.Button("Load Config"),
        sg.In(size=(15, 1), enable_events=True, key="-LOAD FOLDER-"),
        sg.FileBrowse(),
    ],
    [
        sg.Listbox(
            values=displaynames, enable_events=True, size=(60, 20), key="-PARAMS-"
        )
    ],
    [
        sg.Button("Save Config"),
        sg.InputText(size=(25, 1),default_text='config name - RENAME THIS', key='-SAVE FILE NAME-')
    ],
    [
        sg.Text("Save Folder Location"),
        sg.In(size=(15, 1), enable_events=True, key="-SAVE FOLDER-"),
        sg.FolderBrowse(),
        
    ],
    [sg.Text("Save File Type"),
        sg.Combo(["csv","txt",'JSON'], size = 10, 
            default_value='JSON', k="-SAVE STYLE-", enable_events=True, readonly=True)
    ],
    [sg.Button("Clear Params", key = "-CLEAR-"), sg.Button("Clear Params to basics for FAC", key = "-CLEAR TO BASICS-")],
    [sg.Text("Search Params"), sg.Input(size=(20,1), enable_events=True, key="SEARCHBAR")],
    [sg.Listbox(values=[], size=(60,5), key ="-SEARCH RESULTS-")]
]

# For now will only show the name of the file that was chosen
function_caller_column = [
    [sg.Text("Saved Lists (Click to Plot)")],
    [sg.Listbox(
            values=listdisplaynames, enable_events=True, size=(40, 25), key="-LISTS-"
        )],
    [sg.Button("Run Function"),
        sg.Combo(workingfuncnames, size = 20, k="-FUNC NAME-", enable_events=True, readonly=True)],
    [sg.Button("Update Params using First Order Calcs",key='-FAC-')],
]

# ----- Full layout -----
layout = [
    [
        sg.Column(config_column),
        sg.VSeperator(),
        sg.Column(function_caller_column),
    ]
]

window = sg.Window("ROOPERT", layout)

# Run the Event Loop
while True:
    event, values = window.read()
    if event == "Exit" or event == sg.WIN_CLOSED:
        break
    if event == "-LISTS-":
        ylistname = values["-LISTS-"][0]
        ylist = storedlists[ylistname]
        listnames = []
        for key in storedlists.keys():
            if len(storedlists[key]) == len(ylist):
                listnames.append(key)
        xlistname = listSelectPopup(listnames)
        if xlistname is None:
            xlist = np.arange(len(ylist))
            xlistname = "Range of " + ylistname
        else:
            xlist = storedlists[xlistname]
        plt.plot(xlist,ylist)
        #plt.title('Unemployment Rate Vs Year', fontsize=14)
        plt.xlabel(xlistname, fontsize=14)
        plt.ylabel(ylistname, fontsize=14)
        plt.grid(True)
        plt.show(block=False)

    if event == "Run Function":
        function_name = values["-FUNC NAME-"]
        funcoutput = open_function_window(function_name)
        if function_name == 'ChanelSolver':
            pass

        if function_name == 'main':
            params = funcoutput
            displaynames = updateDisplayNames(params)
        if function_name == "RocketEquation_MassFixed":
            for key, value in funcoutput.items():
                storedlists[key] = value
            listdisplaynames = list(storedlists.keys())
        window['-PARAMS-'].update(displaynames)
        window['-LISTS-'].update(listdisplaynames)
    if event == "-PARAMS-":
                paramname = values["-PARAMS-"][0].split(" ")[0]
                inputs_layout = get_input_layout(canDelete=True)
                inputs_window = sg.Window(f"Change {paramname} or add new param", inputs_layout)
                inputs_event, inputs_values = inputs_window.read()
                while inputs_event != 'change' and inputs_event != 'delete' and inputs_event != 'add':
                    inputs_event, inputs_values = inputs_window.read()
                    if inputs_event == "Exit" or inputs_event == sg.WIN_CLOSED:
                        break
                inputs_window.close()
                if inputs_event == 'change':
                    try:
                        try: #try to cast it as a float, if you cant then just leave it as a string
                            params[paramname] = float(inputs_values["input"])
                        except:
                            params[paramname] = inputs_values["input"]
                    except:
                        params[paramname] = None
                    displaynames = updateDisplayNames(params)
                    window['-PARAMS-'].update(displaynames)
                if inputs_event == 'delete':
                    del params[paramname]
                    displaynames = updateDisplayNames(params)
                    window['-PARAMS-'].update(displaynames)
                if inputs_event == 'add':
                    try:
                        try: #try to cast it as a float, if you cant then just leave it as a string
                            params[inputs_values['addname']] = float(inputs_values["input"])
                        except:
                            params[inputs_values['addname']] = inputs_values["input"]
                    except:
                        pass
                    displaynames = updateDisplayNames(params)
                    window['-PARAMS-'].update(displaynames)
    if event == "Save Config":
        if values["-SAVE STYLE-"] == 'JSON':
            folder = values["-SAVE FOLDER-"]
            if folder =="":
                path=os.path.join("Configs",values['-SAVE FILE NAME-']+"."+values["-SAVE STYLE-"])
            else:
                path=os.path.join(folder,values['-SAVE FILE NAME-']+"."+values["-SAVE STYLE-"])
            try:
                with open(path, 'w') as file:
                    file.write(json.dumps(params))
            except:
                params['CEA'] = None
                with open(path, 'w') as file:
                    file.write(json.dumps(params))
                displaynames = updateDisplayNames(params)
                window['-PARAMS-'].update(displaynames)
        if values["-SAVE STYLE-"] == 'csv':
            d = {"variable" : [], "value" : [], "unit/description" : []}
            print(d["variable"])
            for i in displaynames:
                v = i.split("  ")
                d["variable"].append(v[0][0:len(v[0])-2])
                d["value"].append(v[1])
                d["unit/description"].append(v[2])
            pd_params = pd.DataFrame.from_dict(d, orient="index")
            folder = values["-SAVE FOLDER-"]
            if folder == "":
                path=os.path.join("Configs",values['-SAVE FILE NAME-']+"."+values["-SAVE STYLE-"])
            else:
                path=os.path.join(folder, values['-SAVE FILE NAME-'] + "." + values["-SAVE STYLE-"])
            pd_params.to_csv(path)
            
    if event == "Load Config":
        filename = values["-LOAD FOLDER-"]
        with open(filename, 'r') as file:
            params = json.load(file)
            displaynames = updateDisplayNames(params)
        window['-PARAMS-'].update(displaynames)
    if event == "-FAC-":
        params = run_selected_function("FirstOrderCalcs", params)
        displaynames = updateDisplayNames(params)
        window['-PARAMS-'].update(displaynames)
    if event == "-CLEAR-":
        answer = sg.popup_yes_no('Are You Sure You Want to Clear all Params?')
        if answer == 'Yes':
            params = {}
            displaynames = updateDisplayNames(params)
            window['-PARAMS-'].update(displaynames)
    if event == "-CLEAR TO BASICS-":
        answer = sg.popup_yes_no('Are You Sure You Want to Clear to Basics?')
        if answer == 'Yes':
            newparams={}
            for basics in ['thrust' ,
                            'time',
                            'fuelname' ,
                            'oxname' ,
                            'pc' ,
                            'pe' ,
                            'throat_radius_curvature' ,
                            'dp' ]:
                try:
                    newparams[basics] = params[basics]
                except:
                    newparams[basics] = None
            
            params = newparams
            displaynames = updateDisplayNames(params)
            window['-PARAMS-'].update(displaynames)
    if values['SEARCHBAR'] != '':                         # if a keystroke entered in search field
        search = values['SEARCHBAR']
        new_values = [x for x in displaynames if search in x]  # do the filtering
        window['-SEARCH RESULTS-'].update(new_values)     # display in the listbox
    if values['SEARCHBAR'] == '':                         # if searchbox is empty
        window['-SEARCH RESULTS-'].update([])     # display in the listbox
window.close()