# -*- coding: utf-8 -*-
"""
Created on Mon Aug 13 18:30:43 2018

@author: jingwui
"""

### This module is to print the Gene circuit using the modified code from quickplot.py
#   based on DNAplotlib.py

#    import matplotlib
#    matplotlib.use('Agg')

# Other modules we require
#import argparse
import dnaplotlib as dpl
import matplotlib.pyplot as plt


def PlotCircuitfunc(input):
    # Types mapping
    types = {}
    types['p'] = 'Promoter'
    types['i'] = 'Ribozyme'
    types['r'] = 'RBS'
    types['c'] = 'CDS'
    types['t'] = 'Terminator'
    types['s'] = 'Spacer'
    types['='] = 'Scar'
    types['o'] = 'Origin'
    types['es'] = 'EmptySpace'

    # Colours mapping
    colors = {}
    colors['black']  = (0.00,0.00,0.00)
    colors['gray']   = (0.60,0.60,0.60)
    colors['red']    = (0.89,0.10,0.11)
    colors['orange'] = (1.00,0.50,0.00)
    colors['yellow'] = (1.00,1.00,0.00)
    colors['white']  = (1.00,1.00,1.00)
    colors['green']  = (0.20,0.63,0.17)
    colors['blue']   = (0.12,0.47,0.71)
    colors['purple'] = (0.42,0.24,0.60)
    colors['lightred']    = (0.98,0.60,0.60)
    colors['lightorange'] = (0.99,0.75,0.44)
    colors['lightyellow'] = (1.00,1.00,0.60)
    colors['lightgreen']  = (0.70,0.87,0.54)
    colors['lightblue']   = (0.65,0.81,0.89)
    colors['lightpurple'] = (0.79,0.70,0.84)

    # Generate the parts list from the arguments
    part_list = []
    part_idx = 1
    for el in input.split(' '):
        if el != '':
            part_parts = el.split('.')
            if len(part_parts) >= 2:
                part_short_type = part_parts[0]
                part_fwd = True
                if part_short_type[0] == '-':
                    part_fwd = False
                    part_short_type = el[1]
                if part_short_type in types.keys():
                    part_type = types[part_short_type]
                    part_color = part_parts[1]
                    part_rgb = (0,0,0)
                    if part_color in colors.keys():
                        part_rgb = colors[part_color]
                        if len(part_parts) >= 3:
                            part_label = part_parts[2]
                        else:
                            part_label = ''
                        if part_short_type == 'c':
                            part_label_yoffset = 0
                            part_label_color = 'white'
                            part_label_size = 10
                            part_label_style = 'italic'
                        elif part_short_type == 'o':
                            part_label_yoffset = 0; #-10
                            part_label_color = 'black'
                            part_label_size = 8
                            part_label_style = 'normal'
                        elif (part_short_type == 'p') or (part_short_type == 'r'):
                            if len(part_parts) > 3:
                                part_label_yoffset = -5
                                part_label_color = part_parts[3]
                                part_label_size = 9
                                part_label_style = 'normal'
                            else:
                                part_label_yoffset = -5
                                part_label_color = 'black'
                                part_label_size = 9
                                part_label_style = 'normal'
                        else:
                            part_label_yoffset = -5
                            part_label_color = 'black'
                            part_label_size = 9
                            part_label_style = 'normal'
                            if part_parts[0][0] == '-':
                                part_label_yoffset = 5
                                
                    part_list.append( {'name'  : str(part_idx),'type'  : part_type, 'fwd'   : part_fwd,
                                       'opts'  : {'color': part_rgb, 'label': part_label, 'label_y_offset': part_label_yoffset,
                                                  'label_color': part_label_color, 'label_size': part_label_size,
                                                  'label_style': part_label_style}} )
                    
        part_idx = part_idx+1

    #print(part_list) #to check the part_list to ensure the settings are correct 
    return part_list


### This function definition is run and plot the sbol-compliant gene circuit
def Run_PlotCircuit(Input, Regulations = None):

    dr = dpl.DNARenderer()

    ## Parse the command line inputs
    #parser = argparse.ArgumentParser(description="one line quick plot")
    #parser.add_argument("-input",  dest="input",  required=True,  help="\"p.orange p.lightblue i.lightred r.green c.orange t.purple -t.black -c.yellow -p.yellow\"", metavar="string")
    #parser.add_argument("-output", dest="output", required=False, help="output pdf filename")
    #args = parser.parse_args()


    # Input = "p.black.pBAD r.black.rbsD c.red.RFP t.black -t.black.BBaB0015 -c.green.GFP -r.black.rbs34 -i.black -p.black.pTet o.black.CoIE1"

    #print('Hello: ', args.input)

    # Process the arguments
    design = PlotCircuitfunc(Input)

    reg_renderers = dr.std_reg_renderers()
    part_renderers = dr.SBOL_part_renderers()

    #Regulations = None
    # For adding regulation control (Activation, Repression, Connection)
#    Regulations = [{'type': 'Repression', 'from_part': {'start': 40, 'end': 40},
#             'to_part': {'start': 150,'end': 150, 'fwd': True},
#             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]

    # Generate the figure
    fig = plt.figure(figsize=(1.0, 1.0), dpi = 100)
    ax = fig.add_subplot(1,1,1)

    # Plot the design
    dna_start, dna_end = dr.renderDNA(ax, design, part_renderers, Regulations, reg_renderers)
    max_dna_len = dna_end-dna_start

    print('Max Dna length: ', max_dna_len)

    # Format the axis
    ax.set_xticks([])
    ax.set_yticks([])

    # Set bounds
    ax.set_xlim([(-0.0*max_dna_len),
    		        max_dna_len+(0.0*max_dna_len)])
    ax.set_ylim([-25,25])
    ax.set_aspect('equal')
    ax.set_axis_off()

    # Update the size of the figure to fit the constructs drawn
    fig_x_dim = max_dna_len/30#/60.0
    print('x_dim: ', fig_x_dim)
    if fig_x_dim < 1.0:
    	fig_x_dim = 1.0
    fig_y_dim = 1.8 #fig_x_dim*0.4
    plt.gcf().set_size_inches( fig_x_dim, fig_y_dim, forward=True )


    # Save the figure
    plt.tight_layout()
    #fig.savefig(args.output, transparent=True)
    plt.show()


# This function definition serves to test this module
def Main_PlotCircuit():

    plt.close() # close all the figures

#    ReporterColor = 'green.'
#    Reporter = 'GFP'#'sfGFP' 

    ReporterColor = 'red.'
    Reporter = 'RFP'

    #dna length
    #origin = 17.0
    #promoter = 14.0
    #rbs = 14.0
    #c = 32.0
    #t = 12.0
    #= = 10.0
    #p+r+c+t = 72
    #p+r+c+t+o+= = 99
    #p+r+1/2 c = 43 (connector from gene)


    Origin = ['p15A','CoIE1', 'pSB1A3']
    #Gene = ['gRNA', 'dCas9']
    #Gene = ['EL222', 'gRNA']
    #Gene = ['AraC']
    #Gene = ['TetR']
    Gene = ['LacI']

    #Input = '=.black'

#    Input = "p.black r.black c.orange t.black " + \
#            "p.black r.black c.blue t.black " + \
#            "p.black r.black c." + ReporterColor + Reporter + " " + "t.black o.black." + Origin[1]

#    Input = "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[0] + " =.white "+ "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[1]
#
#    RepressorName = input('Please insert the name of repressor: ')
#    Input = "p.black r.black c.blue."+ RepressorName + " t.black "+ "p.black r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[0]
#    print('Input: ', Input)
#    print('Input Type: ', type(Input))
#    Regulations = [{'type': 'Repression', 'from_part': {'start': 40, 'end': 40},
#         'to_part': {'start': 78,'end': 78, 'fwd': True},
#         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]

#    Input = "p.black r.black c.orange" + " t.black o.black." + Origin[0] + " =.white "+ \
#                "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[1]
#
#    Regulations = [{'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
#     'to_part': {'start': 110,'end': 110, 'fwd': True},
#     'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
#     {'type': 'Activation', 'from_part': {'start': 102, 'end': 102},
#     'to_part': {'start': 102,'end': 102, 'fwd': True},
#     'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]


        ### ANDgate ###
#    Input = "p.black r.black c.orange." + Gene[0] + " t.black " + \
#    "p.black r.black c." + ReporterColor + Reporter + " " + "t.black o.black." + Origin[0]
#
#    Regulations = [{'type': 'Activation', 'from_part': {'start': 75, 'end': 75},
#     'to_part': {'start': 75,'end': 75, 'fwd': True},
#     'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}},
#    {'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
#     'to_part': {'start': 83,'end': 83, 'fwd': True},
#     'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}}]

#    Input = "p.black r.black c.orange."+ Gene[0] + " t.black o.black." + Origin[0] + " =.white " + \
#    "p.black r.black c.blue."+ Gene[1] + " t.black o.black." + Origin[1] + " =.white " + \
#    "p.black r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[2]
#
#    Regulations = [{'type': 'Connection', 'from_part': {'start': 43, 'end': 43},
#         'to_part': {'start': 142,'end': 142, 'fwd': True},
#         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
#        {'type': 'Repression', 'from_part': {'start': 142, 'end': 142},
#         'to_part': {'start': 205,'end': 205, 'fwd': True},
#         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
    
#    Input = "p.black r.black c.orange."+ Gene[0] + " t.black o.black." + Origin[1] + " =.white " + \
#                "p.black r.black c.blue."+ Gene[1] + " t.black " + \
#                "p.black r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[0]
#
#    Regulations = [{'type': 'Connection', 'from_part': {'start': 43, 'end': 43},
#         'to_part': {'start': 142,'end': 142, 'fwd': True},
#         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
#        {'type': 'Repression', 'from_part': {'start': 142, 'end': 142},
#         'to_part': {'start': 178,'end': 178, 'fwd': True},
#         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
        
#    ### OR gate
#    Input = "p.black r.black c."+ ReporterColor + Reporter + " t.black "+ "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[1]
#    
#    Regulations = None
    
    ### pBAD inducible (Two parts) single plasmid ###
    
    Input = "p.black r.black c.orange."+ Gene[0] + " t.black " + \
            "p.black.pLac r.black c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[1]

    Regulations = [{'type': 'Repression', 'from_part': {'start': 43, 'end': 43},
         'to_part': {'start': 78,'end': 78, 'fwd': True},
         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5, 'first_arc_y_offset': -3.5}}, 
        {'type': 'Repression', 'from_part': {'start': 66, 'end': 66},
         'to_part': {'start': 66,'end': 66, 'fwd': True},
         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5, 'second_arc_y_offset': 7}}]
    
    ### Inducible (Two parts) double plasmids
    
#    Input = "p.black r.black c.orange."+ Gene[0] + " t.black o.black." + Origin[0] +\
#        " =.white " + "p.black r.black c."+ ReporterColor + Reporter + " " + "t.black o.black." + Origin[1]
#
#    Regulations = [{'type': 'Repression', 'from_part': {'start': 43, 'end': 43},
#         'to_part': {'start': 105,'end': 105, 'fwd': True},
#         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5, 'first_arc_y_offset': -3.5}}, 
#        {'type': 'Repression', 'from_part': {'start': 80, 'end': 80},
#         'to_part': {'start': 80,'end': 80, 'fwd': True},
#         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5, 'second_arc_y_offset': 7}}]
    
    ### pTet inducible (Two parts) ###
    
#    Input = "p.black r.black c.orange."+ Gene[0] + " t.black " + \
#            "p.black r.black c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[1]
#
#    Regulations = [{'type': 'Repression', 'from_part': {'start': 49, 'end': 49},
#         'to_part': {'start': 79,'end': 79, 'fwd': True},
#         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}}, 
#        {'type': 'Repression', 'from_part': {'start': 40, 'end': 40},
#         'to_part': {'start': 40,'end': 40, 'fwd': True},
#         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
    
    ### Inducible Systems (Only Single Part) ###
    #Input = "p.black.pBAD r.black c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[1]
    
#    Input = "p.black.pLac r.black c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[1]
#
#    Regulations = [{'type': 'Activation', 'from_part': {'start': 7, 'end': 7},
#         'to_part': {'start': 7,'end': 7, 'fwd': True},
#         'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}}]
    
#    Regulations = [{'type': 'Activation', 'from_part': {'start': 7, 'end': 7},
#         'to_part': {'start': 7,'end': 7, 'fwd': True},
#         'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]

    ### Constitutive Promoter ###
#    Input = "p.black r.black.rbs64.red c." + ReporterColor + Reporter+ " " + "t.black o.black." + Origin[1]
#    
#    Regulations = None
    
#    Input = "p.black r.black c.orange.JingWui" + " t.black o.black." + Origin[0] + " =.white " \
#    "p.black r.black c.blue" + " t.black o.black." + Origin[1] + " =.white " \
#    "p.black r.black c." + ReporterColor + Reporter + " t.black o.black." + Origin[2]
#
#    Regulations = [{'type': 'Activation', 'from_part': {'start': 44, 'end': 44},
#     'to_part': {'start': 209,'end': 209, 'fwd': True},
#     'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
#     {'type': 'Activation', 'from_part': {'start': 143, 'end': 143},
#     'to_part': {'start': 201,'end': 201, 'fwd': True},
#     'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]

#    Input = "p.black r.black c.orange" + " t.black p.black r.black c.blue t.black o.black." + Origin[0] + " =.white "+ \
#                "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[1]
#
#    Regulations = [{'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
#     'to_part': {'start': 182,'end': 182, 'fwd': True},
#     'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
#     {'type': 'Activation', 'from_part': {'start': 115, 'end': 115},
#     'to_part': {'start': 174,'end': 174, 'fwd': True},
#     'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]

#    Input = "p.black r.black c.orange" + " t.black p.black r.black c.blue t.black o.black." + Origin[0] + " =.white "+ \
#    "p.black r.black c."+ ReporterColor + Reporter + " t.black o.black." + Origin[1]
#
#    Regulations = [{'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
#             'to_part': {'start': 182,'end': 182, 'fwd': True},
#             'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
#             {'type': 'Activation', 'from_part': {'start': 115, 'end': 115},
#             'to_part': {'start': 174,'end': 174, 'fwd': True},
#             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]

#    Input = "p.black r.black c.orange" + " t.black o.black." + Origin[0] + " =.white " \
#            "p.black r.black c.blue" + " t.black o.black." + Origin[1] + " =.white " \
#            "p.black r.black c." + ReporterColor + Reporter + " t.black o.black." + Origin[2]

#    Regulations = [{'type': 'Activation', 'from_part': {'start': 43, 'end': 43},
#             'to_part': {'start': 155,'end': 155, 'fwd': True},
#             'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
#             {'type': 'Activation', 'from_part': {'start': 115, 'end': 115},
#             'to_part': {'start': 147,'end': 147, 'fwd': True},
#             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]

#    Regulations = [{'type': 'Activation', 'from_part': {'start': 44, 'end': 44},
#             'to_part': {'start': 209,'end': 209, 'fwd': True},
#             'opts': {'color': (1.00,0.50,0.00), 'linewidth': 1.5}},
#             {'type': 'Activation', 'from_part': {'start': 143, 'end': 143},
#             'to_part': {'start': 201,'end': 201, 'fwd': True},
#             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
    #Input = "p.black.pBAD r.black.rbsD c.red.RFP t.black -t.black.BBaB0015 -c.green.GFP -r.black.rbs34 -i.black -p.black.pTet o.black.CoIE1"
    #Input = "p.black.pBAD r.black.rbsD c.red.RFP t.black o.blue.E8K =.white -t.black.BBaB0015 -c.green.GFP -r.black.rbs34 -i.black -p.black.pTet o.black.CoIE1"
    #Input = "p.black.pBAD r.black.rbsD c.red.RFP t.black o.blue.E8K"
    print('Input Type: ',type(Input))
    #Regulations = None
    # For adding regulation control (Activation, Repression, Connection)
#    Regulations = [{'type': 'Repression', 'from_part': {'start': 8, 'end': 8},
#             'to_part': {'start': 8,'end': 8, 'fwd': True},
#             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]
#    Regulations = [{'type': 'Repression', 'from_part': {'start': 40, 'end': 40},
#             'to_part': {'start': 105,'end': 105, 'fwd': True},
#             'opts': {'color': (0.12,0.47,0.71), 'linewidth': 1.5}}]

    Run_PlotCircuit(Input, Regulations)


### Run the main function to check if the module functions properly
if __name__ == '__main__':
    Main_PlotCircuit()