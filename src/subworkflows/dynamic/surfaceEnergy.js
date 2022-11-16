import { removeNewLinesAndExtraSpaces } from "@exabyte-io/code.js/dist/utils";

function getIOUnitEndpointOptions(query, projection = "{}") {
    return {
        params: {
            query,
            projection,
        },
    };
}

function getAssignmentUnitInput(unit, name) {
    return [
        {
            name, scope: unit.flowchartId,
        },
    ];
}

const getSurfaceEnergySubworkflowUnits = function ({ scfUnit, unitBuilders }) {
    const { IOUnitConfigBuilder, AssignmentUnitConfigBuilder, AssertionUnitConfigBuilder } = unitBuilders;

    let input,
        endpointOptions;
    endpointOptions = getIOUnitEndpointOptions("{'_id': MATERIAL_ID}");
    const getSlabUnit = new IOUnitConfigBuilder("io-slab", "materials", endpointOptions).build();

    input = getAssignmentUnitInput(getSlabUnit, "DATA");
    const setSlabUnit = new AssignmentUnitConfigBuilder("slab", "SLAB", "DATA[0]", input).build();

    endpointOptions = getIOUnitEndpointOptions("{'_id': SLAB.metadata.bulkId}");
    const getBulkUnit = new IOUnitConfigBuilder("io-bulk", "materials", endpointOptions).build();

    const BULKValue = "DATA[0] if DATA else None";
    input = getAssignmentUnitInput(getBulkUnit, "DATA");
    const setBulkUnit = new AssignmentUnitConfigBuilder("bulk", "BULK", BULKValue, input).build();

    const assertBulkUnit = new AssertionUnitConfigBuilder("assert-bulk", "BULK != None", "Bulk material does not exist!").build();

    const query = removeNewLinesAndExtraSpaces(`{
        'exabyteId': BULK.exabyteId,
        'data.name': 'total_energy',
        'group': {'$regex': ''.join((SUBWORKFLOW.application.shortName, ':'))}
    }`);
    // Do not confuse `sort` with `$sort`. `sort` is used in meteor collections and `$sort` is used in Mongo queries.
    const projection = "{'sort': {'precision.value': -1}, 'limit': 1}";
    endpointOptions = getIOUnitEndpointOptions(query, projection);
    const getEBulkUnit = new IOUnitConfigBuilder("io-e-bulk", "refined-properties", endpointOptions).build();

    input = getAssignmentUnitInput(getEBulkUnit, "DATA");
    const EBULKValue = "DATA[0].data.value if DATA else None";
    const setEBulkUnit = new AssignmentUnitConfigBuilder("e-bulk", "E_BULK", EBULKValue, input).build();

    const assertEBulkUnit = new AssertionUnitConfigBuilder("assert-e-bulk", "E_BULK != None", "E_BULK does not exist!").build();

    const AValue = "np.linalg.norm(np.cross(SLAB.lattice.vectors.a, SLAB.lattice.vectors.b))";
    const setSurfaceUnit = new AssignmentUnitConfigBuilder("surface", "A", AValue).build();

    const setNBulkUnit = new AssignmentUnitConfigBuilder("n-bulk", "N_BULK", "len(BULK.basis.elements)").build();
    const setNSlabUnit = new AssignmentUnitConfigBuilder("n-slab", "N_SLAB", "len(SLAB.basis.elements)").build();

    input = getAssignmentUnitInput(scfUnit, "total_energy");
    const setESlabUnit = new AssignmentUnitConfigBuilder("e-slab", "E_SLAB", "total_energy", input).build();

    const results = [{ name: "surface_energy" }];
    const SEValue = "1 / (2 * A) * (E_SLAB - E_BULK * (N_SLAB/N_BULK))";
    const surfaceEnergyUnit = new AssignmentUnitConfigBuilder("surface-energy", "SURFACE_ENERGY", SEValue, [], results).build();

    return [
        getSlabUnit,
        setSlabUnit,
        getBulkUnit,
        setBulkUnit,
        assertBulkUnit,
        getEBulkUnit,
        setEBulkUnit,
        assertEBulkUnit,
        setSurfaceUnit,
        setNBulkUnit,
        setNSlabUnit,
        scfUnit,
        setESlabUnit,
        surfaceEnergyUnit,
    ];
};

export { getSurfaceEnergySubworkflowUnits };
