&miscparameters
 new_controls=.TRUE.
/
&vtktfmparameters
 make_same=.FALSE.
 max_number_of_panels=50
 angle_units='degrees'
 option='panel'
 number_of_transforms=50,
 preserve_input_structure=.TRUE.
/
&vtkfiles
vtk_input_file='RES.vtk'
vtk_output_file='RES_tfm.vtk'
/
&panelarrayparameters
panel_bodies=1,5,7,12,10,18
panel_transform=1,2,3,4,5,6
/
&positionparameters ! 1
position_transform=2, ! offsca - scale and offset in Cartesians
position_offset=0., 0., 500
position_scale=1.,1.,1.
/
&positionparameters ! 2
position_transform=6, ! poloidal tilt
position_offset=-5.
/
&positionparameters ! 3
position_transform=7, ! toroidal tilt
position_offset=15
/
&positionparameters ! 4
position_transform=12, ! toroidal rotation
position_offset=0., 0., 30
/
&positionparameters ! 5
position_transform=22, ! poloidal rotation
position_offset=0., 0., -10
/
&positionparameters ! 6
position_transform=42, ! displace by first towards (R,Z) given by last two
position_offset=-10., 1000., 0
/
