"""Open a generated C-grid in Gmsh for interactive inspection."""

import argparse
from pathlib import Path

import gmsh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mesh", type=Path, help="Path to the generated .msh file")
    args = parser.parse_args()
    if not args.mesh.is_file():
        parser.error(f"Mesh file not found: {args.mesh}")

    gmsh.initialize(readConfigFiles=False)
    try:
        gmsh.open(str(args.mesh.resolve()))

        # This 2D case is extruded one cell deep. Dimension 3 counts those
        # volume cells and leaves out the surface quads on the boundaries.
        _, volume_element_tags, _ = gmsh.model.mesh.getElements(3)
        cell_count = sum(len(tags) for tags in volume_element_tags)
        if cell_count == 0:
            raise SystemExit(f"{args.mesh.name} contains no volume cells.")
        print(f"{args.mesh.name}: {cell_count:,} volume cells", flush=True)

        # Start in the x-y plane and show the surface grid as a wireframe.
        gmsh.option.setNumber("General.Trackball", 0)
        gmsh.option.setNumber("Mesh.SurfaceEdges", 1)
        gmsh.option.setNumber("Mesh.VolumeEdges", 0)
        gmsh.fltk.run()
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()
