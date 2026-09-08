import JSONSchemasInterface from "@mat3ra/esse/dist/js/esse/JSONSchemasInterface";
import esseSchemas from "@mat3ra/esse/dist/js/schemas.json";
import { Material } from "@mat3ra/made";
import { expect } from "chai";
import type { JSONSchema7 } from "json-schema";

import type { OrderedMaterial } from "../../src/js/context/mixins/MaterialContextMixin";
import NWChemInputDataManager from "../../src/js/context/providers/by_application/nwchem/NWChemInputDataManager";

/**
 * NWChem keeps its geometry inside the input file, so the relaxed geometry reaches the next unit
 * only if that unit renders ` restart nwchem` and omits its own geometry block. Mirrors the
 * RESTART_MODE line in QEPWXInputDataManager.
 */
describe("NWChemInputDataManager restart flag", () => {
    before(() => {
        JSONSchemasInterface.setSchemas(esseSchemas as JSONSchema7[]);
    });

    const buildData = (jobHasParent: boolean, workflowHasRelaxation: boolean) => {
        const material = Material.createDefault() as OrderedMaterial;
        const provider = new NWChemInputDataManager(
            {},
            { material, jobHasParent, workflowHasRelaxation },
        );
        return provider.getDefaultData();
    };

    it("does not restart for a standalone workflow with no relaxation", () => {
        expect(buildData(false, false).RESTART).to.be.false;
    });

    it("restarts when the workflow has a relaxation subworkflow", () => {
        expect(buildData(false, true).RESTART).to.be.true;
    });

    // Espresso restarts on jobHasParent too, because prepare_restart symlinks the parent's whole
    // outdir. NWChem's symlinks only perm/, and the RTDB is at the work-dir root, so a child job
    // would get `restart nwchem` with nothing to restart from.
    it("does not restart merely because the job continues from a parent job", () => {
        expect(buildData(true, false).RESTART).to.be.false;
    });
});
