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
input_dir = '%%%'
input_cones = input_dir+'WVT_data_bins.txt'
input_field = input_dir+'Temperature_Bins.txt'
output = input_dir+'WVT_Diagram'
color = input_dir+'Temperature_Bins.txt'
colormap = 'coolwarm'
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
        RGB = [cmap(random.randint(0,256)) for i in range(n_cones)]
    if color != 'Random':
        field_vals = []#np.linspace(0,256,n_cones)/256 #Have to normalize since There are only 256 color choices....
        with open(input_field, 'r') as f:
            lines_field = f.read().splitlines()
        for line in lines_field[1:]:
            field_vals.append(float(line.split(" ")[1]))
        #field_vals = (np.max(field_vals)-field_vals)/(np.max(field_vals)-np.min(field_vals))
        RGB = [cmap(norm_fVal) for norm_fVal in field_vals]
    with open('Bin_with_temp.txt','w') as f:
        f.write("xposition yposition weight temperature \n")
        for line_num in range(len(x)):
            f.write(str(x[line_num])+" "+str(y[line_num])+" "+str(weight[line_num])+" "+str(field_vals[line_num])+" \n")
    renderView1 = GetActiveViewOrCreate('RenderView')
    #renderView1.ViewSize = [4000, 4000]
    renderView1.OrientationAxesVisibility = 0
    renderView1.UseGradientBackground = 1
    '''for i in range(0,n_cones):
        if i%10 == 0:
            print("Plotting Cone Number "+str(i))
        c1 = Cone()
        c1.Resolution=100
        c1.Center = [x[i],y[i],-1]
        c1.Height = weight[i]
        c1.Radius = 50.0
        c1.Direction = [0,0,1]
        cone = GetActiveSource()
        coneDisplay = GetDisplayProperties(cone, view=renderView1)
        #coneDisplay.DiffuseColor = [RGB[i][0], RGB[i][1], RGB[i][2]]
        Show()
    renderView1.ResetCamera()
    Render()'''
    bin_with_temptxt = CSVReader(FileName=['Bin_with_temp.txt'])

    # Properties modified on bin_with_temptxt
    bin_with_temptxt.FieldDelimiterCharacters = ' '

    # Create a new 'SpreadSheet View'
    spreadSheetView1 = CreateView('SpreadSheetView')
    spreadSheetView1.ColumnToSort = ''
    spreadSheetView1.BlockSize = 1024L
    # uncomment following to set a specific view size
    # spreadSheetView1.ViewSize = [400, 400]

    # get layout
    layout1 = GetLayout()

    # place view in the layout
    layout1.AssignView(2, spreadSheetView1)

    # show data in view
    bin_with_temptxtDisplay = Show(bin_with_temptxt, spreadSheetView1)
    # trace defaults for the display properties.
    bin_with_temptxtDisplay.FieldAssociation = 'Row Data'
    bin_with_temptxtDisplay.CompositeDataSetIndex = [0]

    # destroy spreadSheetView1
    Delete(spreadSheetView1)
    del spreadSheetView1

    # close an empty frame
    layout1.Collapse(2)

    # find view
    renderView1 = FindViewOrCreate('RenderView1', viewtype='RenderView')
    # uncomment following to set a specific view size
    # renderView1.ViewSize = [1092, 788]

    # set active view
    SetActiveView(renderView1)

    # create a new 'Table To Points'
    tableToPoints1 = TableToPoints(Input=bin_with_temptxt)
    tableToPoints1.XColumn = ''
    tableToPoints1.YColumn = ''
    tableToPoints1.ZColumn = ''

    # Properties modified on tableToPoints1
    tableToPoints1.XColumn = 'xposition'
    tableToPoints1.YColumn = 'yposition'
    tableToPoints1.a2DPoints = 1
    tableToPoints1.KeepAllDataArrays = 1

    # show data in view
    tableToPoints1Display = Show(tableToPoints1, renderView1)
    # trace defaults for the display properties.
    tableToPoints1Display.ColorArrayName = [None, '']
    tableToPoints1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    tableToPoints1Display.GlyphType = 'Cone'

    # reset view to fit data
    #renderView1.ResetCamera()

    # create a new 'Glyph'
    glyph1 = Glyph(Input=tableToPoints1,
        GlyphType='Cone')


    # Properties modified on glyph1
    glyph1.GlyphType = 'Cone'
    glyph1.Scalars = ['POINTS', 'weight']
    glyph1.ScaleMode = 'scalar'
    glyph1.ScaleFactor = 1.0
    glyph1.GlyphMode = 'All Points'

    # Properties modified on glyph1.GlyphType
    glyph1.GlyphType.Resolution = 50
    glyph1.GlyphType.Radius = 2.0
    glyph1.GlyphType.Direction = [0.0, 0.0, 1.0]

    # show data in view
    glyph1Display = Show(glyph1, renderView1)
    # trace defaults for the display properties.
    glyph1Display.ColorArrayName = [None, '']
    glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
    glyph1Display.GlyphType = 'Cone'

    # set scalar coloring
    ColorBy(glyph1Display, ('POINTS', 'temperature'))

    # rescale color and/or opacity maps used to include current data range
    glyph1Display.RescaleTransferFunctionToDataRange(True)

    # show color bar/color legend
    glyph1Display.SetScalarBarVisibility(renderView1, True)

    # get color transfer function/color map for 'temperature'
    temperatureLUT = GetColorTransferFunction('temperature')

    # get opacity transfer function/opacity map for 'temperature'
    temperaturePWF = GetOpacityTransferFunction('temperature')

    #### saving camera placements for all active views

    # current camera placement for renderView1
    renderView1.InteractionMode = '2D'
    renderView1.CameraPosition = [435.228145985, 386.645663727, 10000.0]
    renderView1.CameraFocalPoint = [435.228145985, 386.645663727, 0.0]
    renderView1.CameraParallelScale = 116.22463957194203

    #### uncomment the following to render all views
    # RenderAllViews()
    # alternatively, if you want to write images, you can use SaveScreenshot(...).
    SaveScreenshot(output+'.png', magnification=1, quality=100)
    return None

def main():
    Render_Cones(input_cones,input_field,output,color,colormap)
    return None
main()
