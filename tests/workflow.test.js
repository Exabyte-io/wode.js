import { Template } from "@exabyte-io/ade.js";
import { expect } from "chai";

import { ContextProviderRegistry } from "../src/context/registry";
import { createWorkflows } from "../src/workflows";

// patch to make all context providers available
Template.providerRegistry = ContextProviderRegistry;

describe("workflows", () => {
    it("can all be created", () => {
        const workflows = createWorkflows({});
        workflows.map((wf) => {
            // eslint-disable-next-line no-unused-expressions
            expect(wf).to.exist;
            return null;
        });
    });
});
