import { expect } from "chai";
import { _ } from "lodash";

import { Workflow } from "../src/workflows";
import { createWorkflow } from "../src/workflows/create";
import { workflowData as allWorkflowData } from "../src/workflows/workflows";

describe("assignement unit", () => {
    it("can be patched from workflow asset", () => {
        // add patch to workflow config
        const patches = [
            {
                index: 0,
                type: "assignment",
                config: {
                    attributes: { operand: "someOtherOperand", value: 42 },
                },
            },
        ];
        // classification_workflow has assignment units
        const { classification_workflow } = allWorkflowData.workflows.exabyteml;
        const wfConfig = _.cloneDeep(classification_workflow);
        wfConfig.units[0].unitConfigs = patches;

        const workflow = createWorkflow({
            appName: "exabyteml",
            workflowData: wfConfig,
            workflowCls: Workflow,
        });
        const assignmentUnit = workflow.subworkflows[0].units[0];
        expect(assignmentUnit.operand).to.be.equal(patches[0].config.attributes.operand);
        expect(assignmentUnit.value).to.be.equal(patches[0].config.attributes.value);
    });
});
