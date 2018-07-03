// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TopologyOps.hh" // implementation of object methods

#include "pylith/utils/error.hh" // USES PYLITH_CHECK_ERROR

#include <iostream> // USES std::cout
#include <cassert> // USES assert()

// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::createFault(pylith::topology::Mesh* faultMesh,
                                         const pylith::topology::Mesh& mesh,
                                         PetscDMLabel groupField) {
    PYLITH_METHOD_BEGIN;

    assert(faultMesh);
    PetscErrorCode err;

    faultMesh->coordsys(mesh.coordsys());
    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);

    // Convert fault to a DM
    PetscDM subdm = NULL;
    const char *groupName = "";

    if (groupField) {err = DMLabelGetName(groupField, &groupName);PYLITH_CHECK_ERROR(err);}
    err = DMPlexCreateSubmesh(dmMesh, groupField, 1, PETSC_FALSE, &subdm);PYLITH_CHECK_ERROR(err);
    // Check that no cell have all vertices on the fault
    if (groupField) {
        IS subpointIS;
        const PetscInt *dmpoints;
        PetscInt defaultValue, cStart, cEnd, vStart, vEnd;

        err = DMLabelGetDefaultValue(groupField, &defaultValue);PYLITH_CHECK_ERROR(err);
        err = DMPlexCreateSubpointIS(subdm, &subpointIS);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetHeightStratum(subdm, 0, &cStart, &cEnd);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetDepthStratum(dmMesh, 0, &vStart, &vEnd);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(subpointIS, &dmpoints);PYLITH_CHECK_ERROR(err);
        for (PetscInt c = cStart; c < cEnd; ++c) {
            PetscBool invalidCell = PETSC_TRUE;
            PetscInt *closure = NULL;
            PetscInt closureSize;

            err = DMPlexGetTransitiveClosure(dmMesh, dmpoints[c], PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
            for (PetscInt cl = 0; cl < closureSize*2; cl += 2) {
                PetscInt value = 0;

                if ((closure[cl] < vStart) || (closure[cl] >= vEnd)) { continue;}
                err = DMLabelGetValue(groupField, closure[cl], &value);PYLITH_CHECK_ERROR(err);
                if (value == defaultValue) {invalidCell = PETSC_FALSE;break;}
            } // for
            err = DMPlexRestoreTransitiveClosure(dmMesh, dmpoints[c], PETSC_TRUE, &closureSize, &closure);PYLITH_CHECK_ERROR(err);
            if (invalidCell) {
                std::ostringstream msg;
                msg << "Ambiguous fault surface. Cell "<<dmpoints[c]<<" has all of its vertices on the fault.";
                err = ISRestoreIndices(subpointIS, &dmpoints);PYLITH_CHECK_ERROR(err);
                err = ISDestroy(&subpointIS);PYLITH_CHECK_ERROR(err);
                err = DMDestroy(&subdm);PYLITH_CHECK_ERROR(err);
                throw std::runtime_error(msg.str());
            } // if
        } // for
        err = ISRestoreIndices(subpointIS, &dmpoints);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&subpointIS);PYLITH_CHECK_ERROR(err);
    } // if
    err = DMPlexOrient(subdm);PYLITH_CHECK_ERROR(err);

    std::string submeshLabel = "fault_" + std::string(groupName);
    faultMesh->dmMesh(subdm, submeshLabel.c_str());

    PYLITH_METHOD_END;
} // createFault


// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::create(pylith::topology::Mesh* mesh,
                                    const pylith::topology::Mesh& faultMesh,
                                    PetscDMLabel faultBdLabel,
                                    const int materialId) {
    assert(mesh);
    PetscDM sdm = NULL;
    PetscDM dm = mesh->dmMesh();assert(dm);
    PetscDMLabel subpointMap = NULL, label = NULL, mlabel = NULL;
    PetscInt dim, cMax, cEnd, numCohesiveCellsOld;
    PetscErrorCode err;

    // Have to remember the old number of cohesive cells
    err = DMPlexGetHeightStratum(dm, 0, NULL, &cEnd);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetHybridBounds(dm, &cMax, NULL, NULL, NULL);PYLITH_CHECK_ERROR(err);
    numCohesiveCellsOld = cEnd - (cMax < 0 ? cEnd : cMax);
    // Create cohesive cells
    err = DMPlexGetSubpointMap(faultMesh.dmMesh(), &subpointMap);PYLITH_CHECK_ERROR(err);
    err = DMLabelDuplicate(subpointMap, &label);PYLITH_CHECK_ERROR(err);
    err = DMLabelClearStratum(label, mesh->dimension());PYLITH_CHECK_ERROR(err);
    // Fix over-aggressive completion of boundary label
    err = DMGetDimension(dm, &dim);PYLITH_CHECK_ERROR(err);
    if (faultBdLabel && (dim > 2)) {
        PetscIS bdIS;
        const PetscInt *bd;
        PetscInt fStart, fEnd, n, i;

        err = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);PYLITH_CHECK_ERROR(err);
        err = DMLabelGetStratumIS(faultBdLabel, 1, &bdIS);PYLITH_CHECK_ERROR(err);
        err = ISGetLocalSize(bdIS, &n);PYLITH_CHECK_ERROR(err);
        err = ISGetIndices(bdIS, &bd);PYLITH_CHECK_ERROR(err);
        for (i = 0; i < n; ++i) {
            const PetscInt p = bd[i];

            // Remove faces
            if ((p >= fStart) && (p < fEnd)) {
                const PetscInt *edges,   *verts, *supportA, *supportB;
                PetscInt numEdges, numVerts, supportSizeA, sA, supportSizeB, sB, val, bval, e, s;
                PetscBool found = PETSC_FALSE;

                err = DMLabelClearValue(faultBdLabel, p, 1);PYLITH_CHECK_ERROR(err);
                // Remove the cross edge
                err = DMPlexGetCone(dm, p, &edges);PYLITH_CHECK_ERROR(err);
                err = DMPlexGetConeSize(dm, p, &numEdges);PYLITH_CHECK_ERROR(err);
                if (numEdges != 3) {
                    std::ostringstream msg;
                    msg << "Internal error while creating fault mesh. Face "<<p<<" has "<<numEdges<<" edges != 3.";
                    throw std::logic_error(msg.str());
                }
                for (e = 0; e < numEdges; ++e) {
                    err = DMPlexGetCone(dm, edges[e], &verts);PYLITH_CHECK_ERROR(err);
                    err = DMPlexGetConeSize(dm, edges[e], &numVerts);PYLITH_CHECK_ERROR(err);
                    if (numVerts != 2) {
                        std::ostringstream msg;
                        msg << "Internal error while creating fault mesh. Edge "<<edges[e]<<" has "<<numVerts<<" vertices != 2.";
                        throw std::logic_error(msg.str());
                    }
                    err = DMPlexGetSupportSize(dm, verts[0], &supportSizeA);PYLITH_CHECK_ERROR(err);
                    err = DMPlexGetSupport(dm, verts[0], &supportA);PYLITH_CHECK_ERROR(err);
                    for (s = 0, sA = 0; s < supportSizeA; ++s) {
                        err = DMLabelGetValue(label, supportA[s], &val);PYLITH_CHECK_ERROR(err);
                        err = DMLabelGetValue(faultBdLabel, supportA[s], &bval);PYLITH_CHECK_ERROR(err);
                        if (( val >= 0) && ( bval >= 0) ) { ++sA;}
                    }
                    err = DMPlexGetSupportSize(dm, verts[1], &supportSizeB);PYLITH_CHECK_ERROR(err);
                    err = DMPlexGetSupport(dm, verts[1], &supportB);PYLITH_CHECK_ERROR(err);
                    for (s = 0, sB = 0; s < supportSizeB; ++s) {
                        err = DMLabelGetValue(label, supportB[s], &val);PYLITH_CHECK_ERROR(err);
                        err = DMLabelGetValue(faultBdLabel, supportB[s], &bval);PYLITH_CHECK_ERROR(err);
                        if (( val >= 0) && ( bval >= 0) ) { ++sB;}
                    }
                    if ((sA > 2) && (sB > 2)) {
                        err = DMLabelClearValue(faultBdLabel, edges[e], 1);PYLITH_CHECK_ERROR(err);
                        found = PETSC_TRUE;
                        break;
                    }
                }
                if (!found) {
                    std::ostringstream msg;
                    msg << "Internal error while creating fault mesh. Face "<<p<<" has no cross edge.";
                    throw std::logic_error(msg.str());
                }
            }
        }
        err = ISRestoreIndices(bdIS, &bd);PYLITH_CHECK_ERROR(err);
        err = ISDestroy(&bdIS);PYLITH_CHECK_ERROR(err);
    }
    // Completes the set of cells scheduled to be replaced
    err = DMPlexLabelCohesiveComplete(dm, label, faultBdLabel, PETSC_FALSE, faultMesh.dmMesh());PYLITH_CHECK_ERROR(err);
    err = DMPlexConstructCohesiveCells(dm, label, NULL, &sdm);PYLITH_CHECK_ERROR(err);

    err = DMGetLabel(sdm, "material-id", &mlabel);PYLITH_CHECK_ERROR(err);
    if (mlabel) {
        err = DMPlexGetHeightStratum(sdm, 0, NULL, &cEnd);PYLITH_CHECK_ERROR(err);
        err = DMPlexGetHybridBounds(sdm, &cMax, NULL, NULL, NULL);PYLITH_CHECK_ERROR(err);
        assert(cEnd > cMax + numCohesiveCellsOld);
        for (PetscInt cell = cMax; cell < cEnd - numCohesiveCellsOld; ++cell) {
            PetscInt onBd;

            /* Eliminate hybrid cells on the boundary of the split from cohesive label,
             * they are marked with -(cell number) since the hybrid cell number aliases vertices in the old mesh */
            err = DMLabelGetValue(label, -cell, &onBd);PYLITH_CHECK_ERROR(err);
            //if (onBd == dim) continue;
            err = DMLabelSetValue(mlabel, cell, materialId);PYLITH_CHECK_ERROR(err);
        }
    }
    err = DMLabelDestroy(&label);PYLITH_CHECK_ERROR(err);

    PetscReal lengthScale = 1.0;
    err = DMPlexGetScale(dm, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(sdm, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);
    mesh->dmMesh(sdm);
} // create


// ----------------------------------------------------------------------
// Form a parallel fault mesh using the cohesive cell information
void
pylith::faults::TopologyOps::createFaultParallel(pylith::topology::Mesh* faultMesh,
                                                 const pylith::topology::Mesh& mesh,
                                                 const int materialId,
                                                 const char* label) {
    PYLITH_METHOD_BEGIN;

    assert(faultMesh);
    const char    *labelname = "material-id";
    PetscErrorCode err;

    faultMesh->coordsys(mesh.coordsys());

    PetscDM dmMesh = mesh.dmMesh();assert(dmMesh);
    PetscDM dmFaultMesh;

    const PetscBool hasLagrangeConstraints = PETSC_TRUE;
    err = DMPlexCreateCohesiveSubmesh(dmMesh, hasLagrangeConstraints, labelname, materialId, &dmFaultMesh);PYLITH_CHECK_ERROR(err);
    err = DMViewFromOptions(dmFaultMesh, NULL, "-pylith_fault_dm_view");PYLITH_CHECK_ERROR(err);
    err = DMPlexOrient(dmFaultMesh);PYLITH_CHECK_ERROR(err);
    std::string meshLabel = "fault_" + std::string(label);

    PetscReal lengthScale = 1.0;
    err = DMPlexGetScale(dmMesh, PETSC_UNIT_LENGTH, &lengthScale);PYLITH_CHECK_ERROR(err);
    err = DMPlexSetScale(dmFaultMesh, PETSC_UNIT_LENGTH, lengthScale);PYLITH_CHECK_ERROR(err);

    faultMesh->dmMesh(dmFaultMesh, meshLabel.c_str());

    PYLITH_METHOD_END;
} // createFaultParallel


// ----------------------------------------------------------------------
void
pylith::faults::TopologyOps::classifyCellsDM(PetscDM dmMesh,
                                             PetscInt vertex,
                                             const int depth,
                                             const int faceSize,
                                             PetscInt firstCohesiveCell,
                                             PointSet& replaceCells,
                                             PointSet& noReplaceCells,
                                             const int debug) {
    // Replace all cells on a given side of the fault with a vertex on the fault
    PointSet vReplaceCells;
    PointSet vNoReplaceCells;
    const PetscInt *support;
    PetscInt supportSize, s, classifyTotal = 0;
    PetscBool modified = PETSC_FALSE;
    PetscErrorCode err;

    if (debug) {std::cout << "Checking fault vertex " << vertex << std::endl;}
    err = DMPlexGetSupportSize(dmMesh, vertex, &supportSize);PYLITH_CHECK_ERROR(err);
    err = DMPlexGetSupport(dmMesh, vertex, &support);PYLITH_CHECK_ERROR(err);
    for (s = 0; s < supportSize; ++s) {
        const PetscInt point = support[s];

        if (point >= firstCohesiveCell) { return;}
        if (replaceCells.find(point)   != replaceCells.end()) {vReplaceCells.insert(point);}
        if (noReplaceCells.find(point) != noReplaceCells.end()) { vNoReplaceCells.insert(point);}
        modified = PETSC_TRUE;
        ++classifyTotal;
    }
    PetscInt classifySize = vReplaceCells.size() + vNoReplaceCells.size();

    while (modified && (classifySize < classifyTotal)) {
        modified = PETSC_FALSE;
        for (s = 0; s < supportSize; ++s) {
            const PetscInt point = support[s];
            PetscBool classified = PETSC_FALSE;

            if (debug) {
                const PetscInt *cone;
                PetscInt coneSize;

                std::cout << "Checking neighbor " << vertex << std::endl;
                err = DMPlexGetConeSize(dmMesh, vertex, &coneSize);PYLITH_CHECK_ERROR(err);
                err = DMPlexGetCone(dmMesh, vertex, &cone);PYLITH_CHECK_ERROR(err);
                for (PetscInt c = 0; c < coneSize; ++c) {
                    std::cout << "  cone point " << cone[c] << std::endl;
                }
            }
            if (vReplaceCells.find(point) != vReplaceCells.end()) {
                if (debug) { std::cout << "  already in replaceCells" << std::endl;}
                continue;
            } // if
            if (vNoReplaceCells.find(point) != vNoReplaceCells.end()) {
                if (debug) { std::cout << "  already in noReplaceCells" << std::endl;}
                continue;
            } // if
            if (point >= firstCohesiveCell) {
                if (debug) { std::cout << "  already a cohesive cell" << std::endl;}
                continue;
            } // if
             // If neighbor shares a face with anyone in replaceCells, then add
            for (PointSet::const_iterator c_iter = vReplaceCells.begin(); c_iter != vReplaceCells.end(); ++c_iter) {
                const PetscInt *coveringPoints;
                PetscInt numCoveringPoints, points[2];

                points[0] = point;points[1] = *c_iter;
                err = DMPlexGetMeet(dmMesh, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
                err = DMPlexRestoreMeet(dmMesh, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
                if (numCoveringPoints == faceSize) {
                    if (debug) { std::cout << "    Scheduling " << point << " for replacement" << std::endl;}
                    vReplaceCells.insert(point);
                    modified = PETSC_TRUE;
                    classified = PETSC_TRUE;
                    break;
                } // if
            } // for
            if (classified) { continue;}
            // It is unclear whether taking out the noReplace cells will speed this up
            for (PointSet::const_iterator c_iter = vNoReplaceCells.begin(); c_iter != vNoReplaceCells.end(); ++c_iter) {
                const PetscInt *coveringPoints;
                PetscInt numCoveringPoints, points[2];

                points[0] = point;points[1] = *c_iter;
                err = DMPlexGetMeet(dmMesh, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
                err = DMPlexRestoreMeet(dmMesh, 2, points, &numCoveringPoints, &coveringPoints);PYLITH_CHECK_ERROR(err);
                if (numCoveringPoints == faceSize) {
                    if (debug) { std::cout << "    Scheduling " << point << " for no replacement" << std::endl;}
                    vNoReplaceCells.insert(point);
                    modified = PETSC_TRUE;
                    classified = PETSC_TRUE;
                    break;
                } // for
            } // for
        }
        if (debug) {
            std::cout << "classifySize: " << classifySize << std::endl;
            std::cout << "classifyTotal: " << classifyTotal << std::endl;
            std::cout << "vReplaceCells.size: " << vReplaceCells.size() << std::endl;
            std::cout << "vNoReplaceCells.size: " << vNoReplaceCells.size() << std::endl;
        }
        assert(size_t(classifySize) < vReplaceCells.size() + vNoReplaceCells.size());
        classifySize = vReplaceCells.size() + vNoReplaceCells.size();
        if (classifySize > classifyTotal) {
            std::ostringstream msg;
            msg << "Internal error classifying cells during creation of cohesive cells."
                << "  classifySize: " << classifySize << ", classifyTotal: " << classifyTotal;
            throw std::logic_error(msg.str());
        } // if
    }
    replaceCells.insert(vReplaceCells.begin(), vReplaceCells.end());
    // More checking
    noReplaceCells.insert(vNoReplaceCells.begin(), vNoReplaceCells.end());
} // classifyCellsDM


// End of file
