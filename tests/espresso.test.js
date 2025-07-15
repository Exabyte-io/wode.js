import { Constraint } from "@mat3ra/made/dist/js/constraints/constraints";
import { Material } from "@mat3ra/made/dist/js/material";
import { expect } from "chai";

import QEPWXContextProvider from "../src/context/providers/by_application/espresso/QEPWXContextProvider";

describe("QEPWXContextProvider.atomicPositionsWithConstraints", () => {
    const expectedOutput = `Si     0.000000000    0.000000000    0.000000000 1 0 1
Si     0.250000000    0.250000000    0.250000000 0 0 0`;

    it("returns each atom on its own line when input is array", () => {
        const material = Material.createDefault();
        const constraints = [
            new Constraint({ id: 0, value: [true, false, true] }),
            new Constraint({ id: 1, value: [false, false, false] }),
        ];
        material.setBasisConstraints(constraints);
        const result = QEPWXContextProvider.atomicPositionsWithConstraints(material);
        expect(result).to.equal(expectedOutput);
    });
});
