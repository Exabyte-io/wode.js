import { Application } from "@exabyte-io/ade.js";
import { expect } from "chai";

import { createUnit } from "../src/subworkflows/create";
import { builders } from "../src/units/builders";
import { UnitFactory } from "../src/units/factory";
import { createWorkflows } from "../src/workflows";

describe("units", () => {
    it("can be cloned with new flowchartId", () => {
        const workflows = createWorkflows({});
        const exampleWorkflow = workflows[0];
        const exampleSubworkflow = exampleWorkflow.subworkflows[0];
        const exampleUnit = exampleSubworkflow.units[0];
        const exampleUnitClone = exampleUnit.clone();
        // eslint-disable-next-line no-unused-expressions
        expect(exampleUnitClone).to.exist;
        expect(exampleUnit.flowchartId).to.not.equal(exampleUnitClone.flowchartId);
    });

    it("can create execution unit", () => {
        const unit = createUnit({
            config: {
                type: "executionBuilder",
                config: {
                    name: "test",
                    execName: "pw.x",
                    flavorName: "pw_scf",
                    flowchartId: "test",
                },
            },
            application: new Application({ name: "espresso" }),
            unitBuilders: builders,
            unitFactoryCls: UnitFactory,
        });

        const expectedResults = [
            { name: "atomic_forces" },
            { name: "fermi_energy" },
            { name: "pressure" },
            { name: "stress_tensor" },
            { name: "total_energy" },
            { name: "total_energy_contributions" },
            { name: "total_force" },
        ];

        expect(unit.flavor.results).to.deep.equal(expectedResults);
        expect(unit.results).to.deep.equal(expectedResults);
    });
});
