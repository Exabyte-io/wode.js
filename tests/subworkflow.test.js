import { expect } from "chai";

import { createSubworkflowByName } from "../src/subworkflows";

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
});
