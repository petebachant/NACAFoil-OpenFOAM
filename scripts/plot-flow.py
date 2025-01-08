"""Plot the flowfield."""

import argparse
import glob

import pyvista


def plot_flow(case: str):
    vtk_files = glob.glob(
        f"cases/{case}/postProcessing/surfaces/*/midPlane.vtp"
    )
    assert len(vtk_files) == 1, "Could not find exactly one file"
    data = pyvista.read(vtk_files[0])
    plotter = pyvista.Plotter(off_screen=True)
    plotter.add_mesh(data, scalars="U")
    plotter.view_xz()
    plotter.camera.focal_point = (1, 0, 0)
    plotter.camera.zoom(6.0)
    plotter.screenshot(f"figures/{case}-umag.png")
    plotter.export_html(f"figures/{case}-umag.html")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("angle_of_attack")
    args = parser.parse_args()
    case = f"naca0012-re2e5-aoa-{args.angle_of_attack}"
    plot_flow(case)
