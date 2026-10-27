import os
import shutil
import openmc

shutil.copy('geometry/materials.xml', 'materials.xml')
shutil.copy('geometry/geometry.xml', 'geometry.xml')

geometry = openmc.Geometry.from_xml('geometry.xml')
cells = geometry.get_all_cells()

fuel_cell = None
clad_cell = None
water_cell = None
control_cell = None
reflector_cell = None

for cell in cells.values():
    if cell.fill and hasattr(cell.fill, 'name'):
        if cell.fill.name == 'fuel':
            fuel_cell = cell
        elif cell.fill.name == 'cladding':
            clad_cell = cell
        elif cell.fill.name == 'moderator':
            if water_cell is None:
                water_cell = cell
        elif cell.fill.name == 'control':
            control_cell = cell
        elif cell.fill.name == 'reflector':
            reflector_cell = cell

plot_xy = openmc.Plot()
plot_xy.filename = 'geometry_xy'
plot_xy.basis = 'xy'
plot_xy.origin = [0.0, 0.0, 0.0]
plot_xy.width = [160.0, 160.0]
plot_xy.pixels = [800, 800]
plot_xy.color_by = 'cell'
plot_xy.colors = {
    fuel_cell: (255, 0, 0),
    clad_cell: (128, 128, 128),
    water_cell: (0, 0, 255),
    control_cell: (0, 255, 0),
    reflector_cell: (255, 255, 0)
}

plot_xz = openmc.Plot()
plot_xz.filename = 'geometry_xz'
plot_xz.basis = 'xz'
plot_xz.origin = [0.0, 0.0, 0.0]
plot_xz.width = [160.0, 450.0]
plot_xz.pixels = [800, 800]
plot_xz.color_by = 'cell'
plot_xz.colors = {
    fuel_cell: (255, 0, 0),
    clad_cell: (128, 128, 128),
    water_cell: (0, 0, 255),
    control_cell: (0, 255, 0),
    reflector_cell: (255, 255, 0)
}

plots = openmc.Plots([plot_xy, plot_xz])
plots.export_to_xml('plots.xml')

openmc.plot_geometry()