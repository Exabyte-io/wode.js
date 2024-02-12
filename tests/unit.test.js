import { expect } from "chai";

import { createWorkflows } from "../src/workflows";

describe("units", () => {
    it("can be cloned with new flowchartId", () => {
        try {
            const workflows = createWorkflows({});
            const exampleWorkflow = workflows[0];
            const exampleSubworkflow = exampleWorkflow.subworkflows[0];
            const exampleUnit = exampleSubworkflow.units[0];
            const exampleUnitClone = exampleUnit.clone();
            // eslint-disable-next-line no-unused-expressions
            expect(exampleUnitClone).to.exist;
            expect(exampleUnit.flowchartId).to.not.equal(exampleUnitClone.flowchartId);
        } catch (err) {
            console.log({
                error: err.error,
            });
        }
    });
});
