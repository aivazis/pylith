#!/usr/bin/env python3

import gmsh

DOMAIN_X = DOMAIN_Y = 8.0e+3
DX = 4.0e+3

def generate_mesh(cell="tri"):

    gmsh.initialize()
    gmsh.model.add("mesh")

    x0 = -0.5 * DOMAIN_X
    y0 = -0.5 * DOMAIN_Y
    v0 = gmsh.model.geo.addPoint(x0, y0, 0, DX)
    v1 = gmsh.model.geo.addPoint(x0+DOMAIN_X, y0, 0, DX)
    v2 = gmsh.model.geo.addPoint(x0+DOMAIN_X, y0+DOMAIN_Y, 0, DX)
    v3 = gmsh.model.geo.addPoint(x0, y0+DOMAIN_Y, 0, DX)

    e0 = gmsh.model.geo.addLine(v0, v1)
    e1 = gmsh.model.geo.addLine(v1, v2)
    e2 = gmsh.model.geo.addLine(v2, v3)
    e3 = gmsh.model.geo.addLine(v3, v0)

    c0 = gmsh.model.geo.addCurveLoop([e0, e1, e2, e3])
    s0 = gmsh.model.geo.addPlaneSurface([c0])

    gmsh.model.geo.synchronize()

    boundary = gmsh.model.addPhysicalGroup(0, [e0, e1, e2, e3])
    gmsh.model.setPhysicalName(1, boundary, "boundary")

    material = gmsh.model.addPhysicalGroup(2, [s0], 24)
    gmsh.model.setPhysicalName(2, material, "material-id")

    # right angle tris
    gmsh.option.setNumber("Mesh.Algorithm", 8)

    if cell == "quad":
        gmsh.option.setNumber("Mesh.RecombineAll", 1)

    gmsh.model.mesh.generate(2)
    gmsh.write(f"{cell}_gmsh.msh")

    gmsh.fltk.run()
    gmsh.finalize()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--cell", action="store", dest="cell", required=True, choices=["tri","quad"])
    args = parser.parse_args()

    generate_mesh(args.cell)
