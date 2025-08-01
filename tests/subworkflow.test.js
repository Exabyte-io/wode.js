import ApplicationRegistry from "@exabyte-io/ade.js/dist/js/ApplicationRegistry";
import { expect } from "chai";

import { createSubworkflowByName } from "../src/subworkflows";
import { AssignmentUnit, ConditionUnit } from "../src/units";

const assignmentUnitData = {
    type: "assignment",
    application: { name: "espresso", version: "6.3" },
};

const conditionUnitData = {
    type: "condition",
    application: { name: "espresso", version: "6.3" },
};

describe("subworkflows", () => {
    it("have updateContext function", () => {
        const subworkflow = createSubworkflowByName({
            appName: "espresso",
            swfName: "total_energy",
        });
        expect(typeof subworkflow.updateContext).to.be.equal("function");
    });
    it("can update context", () => {
        const subworkflow = createSubworkflowByName({
            appName: "espresso",
            swfName: "total_energy",
        });
        const newContext = { testKey: "testValue" };
        subworkflow.updateContext(newContext);
        expect(subworkflow.context).to.include(newContext);
    });
    it("add unit to list end", () => {
        const subworkflow = createSubworkflowByName({
            appName: "espresso",
            swfName: "total_energy",
        });

        expect(subworkflow.units.length).to.be.equal(1);
        expect(subworkflow.units[0]._json.type).to.be.equal("execution");

        const assignementUnit = new AssignmentUnit(assignmentUnitData);
        subworkflow.addUnit(assignementUnit, -1);

        expect(subworkflow.units.length).to.be.equal(2);
        expect(subworkflow.units[0]._json.type).to.be.equal("execution");
        expect(subworkflow.units[1]._json.type).to.be.equal("assignment");
    });
    it("add unit to list head", () => {
        const subworkflow = createSubworkflowByName({
            appName: "espresso",
            swfName: "total_energy",
        });

        expect(subworkflow.units.length).to.be.equal(1);
        expect(subworkflow.units[0]._json.type).to.be.equal("execution");

        const assignementUnit = new AssignmentUnit(assignmentUnitData);
        subworkflow.addUnit(assignementUnit, 0);

        expect(subworkflow.units.length).to.be.equal(2);
        expect(subworkflow.units[0]._json.type).to.be.equal("assignment");
        expect(subworkflow.units[1]._json.type).to.be.equal("execution");
    });
    it("add unit in the middle list of two", () => {
        const subworkflow = createSubworkflowByName({
            appName: "espresso",
            swfName: "total_energy",
        });
        expect(subworkflow.units.length).to.be.equal(1);
        expect(subworkflow.units[0]._json.type).to.be.equal("execution");

        const assignementUnit = new AssignmentUnit(assignmentUnitData);
        const conditionUnit = new ConditionUnit(conditionUnitData);
        subworkflow.addUnit(assignementUnit, -1);

        expect(subworkflow.units.length).to.be.equal(2);
        expect(subworkflow.units[0]._json.type).to.be.equal("execution");
        expect(subworkflow.units[1]._json.type).to.be.equal("assignment");

        subworkflow.addUnit(conditionUnit, 1);

        expect(subworkflow.units.length).to.be.equal(3);
        expect(subworkflow.units[0]._json.type).to.be.equal("execution");
        expect(subworkflow.units[1]._json.type).to.be.equal("condition");
        expect(subworkflow.units[2]._json.type).to.be.equal("assignment");
    });
    it("can update application", () => {
        const subworkflow = createSubworkflowByName({
            appName: "espresso",
            swfName: "total_energy",
        });

        const assignementUnit = new AssignmentUnit(assignmentUnitData);
        subworkflow.addUnit(assignementUnit, -1);

        expect(subworkflow.units.length).to.be.equal(2);
        expect(subworkflow.units[0]._json.type).to.be.equal("execution");
        expect(subworkflow.units[1]._json.type).to.be.equal("assignment");
        expect(subworkflow.units[0].application.version).to.be.equal("6.3");
        expect(subworkflow.units[1].application?.version).to.be.equal(undefined);

        const newApplication = ApplicationRegistry.createApplication({
            name: "espresso",
            version: "6.7.0",
        });

        expect(newApplication.version).to.be.equal("6.7.0");

        subworkflow.setApplication(newApplication);

        expect(subworkflow.application.version).to.be.equal("6.7.0");
        expect(subworkflow.units[0].application?.version).to.be.equal("6.7.0");
        expect(subworkflow.units[1].application?.version).to.be.equal(undefined);
    });
});
