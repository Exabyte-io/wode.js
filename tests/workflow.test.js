import { expect } from "chai";

import { createWorkflows } from "../src/workflows";

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
