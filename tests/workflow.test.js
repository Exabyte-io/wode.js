import { expect } from "chai";

import { createWorkflows } from "../src/workflows";
import { createWorkflow } from "../src/workflows/create";
import { workflowData as allWorkflowData } from "../src/workflows/workflows";

describe("workflows", () => {
    it("can all be created", () => {
        const workflows = createWorkflows({});
        workflows.map((wf) => {
            // eslint-disable-next-line no-unused-expressions
            expect(wf).to.exist;
            // eslint-disable-next-line no-unused-expressions
            expect(wf.isValid()).to.be.true;

            // expect(wf.validate()).to.be.true;
            return null;
        });
    });
});

describe("workflow property", () => {
    it("isMultiMaterial is read correctly", () => {
        // Nudged Elastic Band is multi-material
        const mmWorkflow = createWorkflow({
            appName: "espresso",
            workflowData: allWorkflowData.workflows.espresso.neb,
        });
        // eslint-disable-next-line no-unused-expressions
        expect(mmWorkflow.isMultiMaterial).to.be.true;
    });
});
