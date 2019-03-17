#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'Legacy VTK Reader'
PartVTK = LegacyVTKReader(FileNames=['PART.vtk'])

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1221, 697]

# get color transfer function/color map for 'Body'
bodyLUT = GetColorTransferFunction('Body')

# show data in view
PartVTKDisplay = Show(PartVTK, renderView1)
# trace defaults for the display properties.
PartVTKDisplay.ColorArrayName = ['CELLS', 'Body']
PartVTKDisplay.LookupTable = bodyLUT
PartVTKDisplay.ScalarOpacityUnitDistance = 105.85501319308176

# show color bar/color legend
PartVTKDisplay.SetScalarBarVisibility(renderView1, True)

# get opacity transfer function/opacity map for 'Body'
bodyPWF = GetOpacityTransferFunction('Body')

## Properties modified on PartVTKDisplay
#PartVTKDisplay.Opacity = 0.99

# current camera placement for renderView1
renderView1.CameraPosition = [-436.04734004418157, -12146.694241507623, 6942.194098078725]
renderView1.CameraFocalPoint = [486.8189234974118, 13169.176910129594, -7147.388758230402]
renderView1.CameraViewUp = [-0.0006240608698340998, 0.4863249784298808, 0.8737777897744979]
renderView1.CameraParallelScale = 6200.376426408011

## reset view to fit data
renderView1.ResetCamera()

## current camera placement for renderView1
#renderView1.CameraPosition = [-1459.4317894886253, 321.54710564598145, -2344.1436059056227]
#renderView1.CameraFocalPoint = [2781.218505859375, -0.3740234375, -1591.9608764648438]
#renderView1.CameraViewUp = [-0.17399172352685593, 0.00882810469225902, 0.9847075427311891]
#renderView1.CameraParallelScale = 1117.8024810146403

# save screenshot
SaveScreenshot('PART_clup.png', magnification=1, quality=100, view=renderView1)

quit()
