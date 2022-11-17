import { UNIT_TYPES } from "../../../enums";

/**
 * @summary Get QptIrr units used in phonon map calculations
 * @param unitBuilders {Object} unit builders
 * @param unitFactoryCls {*} unit factory class
 * @param application {*} application instance
 * @returns {[{head: boolean, preProcessors: [], postProcessors: [], name: *, flowchartId: *, type: *, results: [], monitors: []},*]}
 */
function getQpointIrrep({ unitBuilders, unitFactoryCls, application }) {
    const { ExecutionUnitConfigBuilder } = unitBuilders;

    const pythonUnit = new ExecutionUnitConfigBuilder(
        "python",
        application,
        "python",
        "espresso_xml_get_qpt_irr",
    ).build();

    const assignmentUnit = unitFactoryCls.create({
        type: UNIT_TYPES.assignment,
        input: [
            {
                scope: pythonUnit.flowchartId,
                name: "STDOUT",
            },
        ],
        operand: "Q_POINTS",
        value: "json.loads(STDOUT)",
    });

    return [pythonUnit, assignmentUnit];
}

export { getQpointIrrep };
