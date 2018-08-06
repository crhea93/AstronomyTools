'''
DEFINE INPUTS IN MAIN FUNCTION AT END OF FILE
This file will create a WVT diagram created from the WVT bin information using
paraview as a renderer (a word??). You can either color the cones randomly or
based off of an external field variable such as temperature
---------------------------------------------------
---------------------------------------------------
Inputs:
    input_dir - Full path of input directory (e.g. '/home/user/Documents/TestData/')
    input_cones - Full path of WVT bin-specific Location information containing centroid values (e.g. input_dir+'WVT_data_bins.txt')
    input_field - Full path of WVT bin-specific Field information (e.g. input_dir+'WVT_temperature_bins.txt')
        If you are using random colors, then it doesn't matter what you put
    output - Full path of output directory and file name (e.g. input_dir+'WVT_Diagram')
    color - Full path to external field variable file (e.g. input_dir+'WVT_bin_field.txt')
        If you wish for the colors to be random just input "Random"
    colormap - Name of colormap for cones (e.g. 'viridis')
---------------------------------------------------
List of Functions (In order of appearance):
    Render_Cones --> Script to create cones
---------------------------------------------------
ADDITIONAL INFORMATION:
    -- List of colormaps found here:
        https://matplotlib.org/users/colormaps.html
    -- Algorithm description found here:
        https://spacesymmetrystructure.wordpress.com/2009/03/16/is-there-anything-new-to-say-about-voronoi-diagrams/
---------------------------------------------------
---------------------------------------------------
Carter Rhea
https://carterrhea.com
carterrhea93@gmail.com
'''

from paraview.simple import *
import random
import numpy as np
from matplotlib import cm

#---------------------INPUTS---------------------#
input_dir = '/home/crhea/Desktop/12036/'
input_cones = input_dir+'WVT_data_bins.txt'
input_field = input_dir+'WVT_temperature_bins.txt'
output = input_dir+'WVT_Diagram'
color = 'Random'
colormap = 'plasma'
#-----------------------------------------------#


# Create a cone and assign it as the active object

def Render_Cones(input_cones,input_field,output,color,colormap):
    x = []
    y = []
    weight = []
    with open(input_cones, 'r') as f:
        lines = f.read().splitlines()
    for line in lines[2:]:
        x.append(float(line.split(" ")[0]))
        y.append(float(line.split(" ")[1]))
        weight.append(float(line.split(" ")[2]))
    n_cones = len(lines[2:])
    cmap = cm.get_cmap(colormap)
    if color == 'Random':
        print(color)
        RGB = [cmap(random.randint(0,256)) for i in range(n_cones)]
    if color != 'Random':
        print(color)
        field_vals = []#np.linspace(0,256,n_cones)/256 #Have to normalize since There are only 256 color choices....
        with open(input_field, 'r') as f:
            lines_field = f.read().splitlines()
        for line in lines_field[2:]:
            field_vals.append(line.split(" ")[1])
        RGB = [cmap(norm_fVal) for norm_fVal in field_vals]
    renderView1 = GetActiveViewOrCreate('RenderView')
    renderView1.ViewSize = [4000, 4000]
    renderView1.OrientationAxesVisibility = 0
    renderView1.UseGradientBackground = 1
    for i in range(0,n_cones):
        if i%10 == 0:
            print("Plotting Cone Number "+str(i))
        c1 = Cone()
        c1.Resolution=100
        c1.Center = [x[i],y[i],-1]
        c1.Height = weight[i]
        c1.Radius = 5.0
        c1.Direction = [0,0,1]
        cone = GetActiveSource()
        coneDisplay = GetDisplayProperties(cone, view=renderView1)
        coneDisplay.DiffuseColor = [RGB[i][0], RGB[i][1], RGB[i][2]]
        Show()
    renderView1.ResetCamera()
    Render()
    SaveScreenshot(output+'.png', magnification=1, quality=100)
    return None

def main():
    Render_Cones(input_cones,input_field,output,color,colormap)
    return None
main()
