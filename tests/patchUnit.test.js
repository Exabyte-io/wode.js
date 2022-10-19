import { expect } from "chai";

import { Workflow } from "../src/workflows";
import { createWorkflow } from "../src/workflows/create";
import { workflowData as allWorkflowData } from "../src/workflows/workflows";

describe("assignement unit", () => {
    it("can be patched from workflow asset", () => {
        // add patch to workflow config
        const patches = [
            {
                index: 0,
                setProp: { operand: "someOtherOperand", value: 42 },
                type: "assignment",
            },
        ];
        // classification_workflow has assignment units
        const { classification_workflow } = allWorkflowData.workflows.exabyteml;
        classification_workflow.units[0].patchUnitConfig = patches;

        const workflow = createWorkflow({
            appName: "exabyteml",
            workflowData: classification_workflow,
            workflowCls: Workflow,
        });
        const assignmentUnit = workflow.subworkflows[0].units[0];
        expect(assignmentUnit.operand).to.be.equal(patches[0].setProp.operand);
        expect(assignmentUnit.value).to.be.equal(patches[0].setProp.value);
    });
});
