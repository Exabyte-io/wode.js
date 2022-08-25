import { expect } from "chai";

import { createWorkflows } from "../src/workflows";

describe("units", () => {
    it("can be cloned with new flowchartId", () => {
        const workflows = createWorkflows({});
        const exampleWorkflow = workflows[0];
        const exampleSubworkflow = exampleWorkflow.subworkflows[0];
        const exampleUnit = exampleSubworkflow.units[0];
        const exampleUnitClone = exampleUnit.clone();
        expect(exampleUnitClone).to.exist;
        expect(exampleUnit.flowchartId).to.not.equal(exampleUnitClone.flowchartId);
    });
});
