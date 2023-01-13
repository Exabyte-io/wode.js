import { JSONSchemaFormDataProvider, MaterialContextMixin } from "@exabyte-io/code.js/dist/context";
import { deepClone } from "@exabyte-io/code.js/dist/utils";
import { Made } from "@exabyte-io/made.js";
import { mix } from "mixwith";

export class BoundaryConditionsFormDataProvider extends mix(JSONSchemaFormDataProvider).with(
    MaterialContextMixin,
) {
    static Material = Made.Material;

    get boundaryConditions() {
        return this.material.metadata.boundaryConditions || {};
    }

    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return {
            type: this.boundaryConditions.type || "pbc",
            offset: this.boundaryConditions.offset || 0,
            electricField: 0,
            targetFermiEnergy: 0,
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {
            type: { "ui:disabled": true },
            offset: { "ui:disabled": true },
            electricField: {},
            targetFermiEnergy: {},
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get humanName() {
        return "Boundary Conditions";
    }

    yieldDataForRendering() {
        const data = deepClone(this.yieldData());
        data.boundaryConditions.offset *= Made.coefficients.ANGSTROM_TO_BOHR;
        data.boundaryConditions.targetFermiEnergy *= Made.coefficients.EV_TO_RY;
        data.boundaryConditions.electricField *= Made.coefficients.EV_A_TO_RY_BOHR;
        return data;
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            type: "object",
            properties: {
                type: {
                    type: "string",
                    title: "Type",
                    default: this.defaultData.type,
                },
                offset: {
                    type: "number",
                    title: "Offset (A)",
                    default: this.defaultData.offset,
                },
                electricField: {
                    type: "number",
                    title: "Electric Field (eV/A)",
                    default: this.defaultData.electricField,
                },
                targetFermiEnergy: {
                    type: "number",
                    title: "Target Fermi Energy (eV)",
                    default: this.defaultData.targetFermiEnergy,
                },
            },
        };
    }
}
