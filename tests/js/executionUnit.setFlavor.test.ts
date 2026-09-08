import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import esseSchemas from "@mat3ra/esse/dist/js/schemas.json";
import type { ApplicationSchema } from "@mat3ra/esse/dist/js/types";
import { ApplicationRegistry } from "@mat3ra/standata";
import StandataDriver from "@mat3ra/standata/dist/js/StandataDriver";
import { expect } from "chai";
import type { JSONSchema7 } from "json-schema";

import { ExecutionUnit } from "../../src/js";

const HELLO_WORLD = "hello_world";
const MATTERSIM = "mlff:mattersim:cell_relaxation";
const LASSO = "pyml:model:lasso_regression:sklearn";

const inputNames = (unit: ExecutionUnit) => unit.input.map((item) => item.template.name);
const inputContents = (unit: ExecutionUnit) => unit.input.map((item) => item.template.content);

describe("ExecutionUnit input files", () => {
    let pythonApplication: ApplicationSchema;

    before(() => {
        JSONSchemasInterface.setSchemas(esseSchemas as JSONSchema7[]);
    });

    beforeEach(() => {
        ApplicationRegistry.setDriver(new StandataDriver());
        pythonApplication = new ApplicationRegistry().findApplication({ name: "python" });
    });

    const createHelloWorldUnit = () => {
        const unit = new ExecutionUnit({ application: pythonApplication });
        expect(unit.flavor.name).to.equal(HELLO_WORLD);
        return unit;
    };

    it("takes the input files of the flavor it is switched to", () => {
        const unit = createHelloWorldUnit();
        expect(inputNames(unit)).to.deep.equal(["script.py", "requirements.txt"]);

        unit.setFlavor(MATTERSIM);

        expect(inputNames(unit)).to.deep.equal(["script.py", "utils.py", "requirements.txt"]);
        expect(inputContents(unit)).to.deep.equal(
            new ApplicationRegistry()
                .getInput(pythonApplication, unit.flavor)
                .map((template) => template.content),
        );
    });

    it("replaces same-named input files carried over from the previous flavor", () => {
        const unit = createHelloWorldUnit();
        const helloWorldContents = inputContents(unit);

        unit.setFlavor(MATTERSIM);

        // "script.py" and "requirements.txt" exist under both flavors, so a stale row would look
        // like a match by name and keep serving hello_world content.
        expect(unit.input[0].template.content).to.not.equal(helloWorldContents[0]);
        expect(unit.input[2].template.content).to.not.equal(helloWorldContents[1]);
    });

    it("renames input files a flavor declares under a different name", () => {
        const unit = createHelloWorldUnit();

        unit.setFlavor(LASSO);

        expect(inputNames(unit)).to.deep.equal([
            "model_lasso_regression_sklearn.py",
            "requirements.txt",
        ]);
    });

    it("switches back to the input files of the original flavor", () => {
        const unit = createHelloWorldUnit();
        const helloWorldInput = unit.toJSON().input;

        unit.setFlavor(MATTERSIM);
        unit.setFlavor(HELLO_WORLD);

        expect(unit.input).to.deep.equal(helloWorldInput);
    });

    it("drops manual edits made under the previous flavor", () => {
        const unit = createHelloWorldUnit();
        unit.input = unit.input.map((item) => ({
            ...item,
            template: { ...item.template, content: "print('edited')" },
            isManuallyChanged: true,
        }));

        unit.setFlavor(MATTERSIM);

        expect(inputContents(unit)).to.not.include("print('edited')");
        expect(unit.input.every((item) => !item.isManuallyChanged)).to.equal(true);
    });

    it("keeps manual edits when reconstructed from saved JSON", () => {
        const unit = createHelloWorldUnit();
        unit.setFlavor(MATTERSIM);
        unit.input = unit.input.map((item, index) =>
            index === 1
                ? {
                      ...item,
                      template: { ...item.template, content: "# edited utils" },
                      isManuallyChanged: true,
                  }
                : item,
        );

        const reconstructed = new ExecutionUnit(unit.toJSON());

        expect(inputNames(reconstructed)).to.deep.equal([
            "script.py",
            "utils.py",
            "requirements.txt",
        ]);
        expect(reconstructed.input[1].template.content).to.equal("# edited utils");
        expect(reconstructed.input[1].isManuallyChanged).to.equal(true);
    });

    it("keeps input files renamed by the user when reconstructed from saved JSON", () => {
        const unit = createHelloWorldUnit();
        unit.input = unit.input.map((item, index) =>
            index === 0 ? { ...item, template: { ...item.template, name: "my_script.py" } } : item,
        );

        const reconstructed = new ExecutionUnit(unit.toJSON());

        expect(inputNames(reconstructed)).to.deep.equal(["my_script.py", "requirements.txt"]);
    });

    it("never fills two input slots from the same saved row", () => {
        const unit = createHelloWorldUnit();
        unit.setFlavor(MATTERSIM);
        // A saved unit whose "utils.py" row is missing: the remaining rows must not be reused to
        // fill the gap, or "requirements.txt" would show up twice.
        const json = unit.toJSON();
        json.input = json.input.filter((item) => item.template.name !== "utils.py");

        const reconstructed = new ExecutionUnit(json);

        expect(inputNames(reconstructed)).to.deep.equal([
            "script.py",
            "utils.py",
            "requirements.txt",
        ]);
    });
});
