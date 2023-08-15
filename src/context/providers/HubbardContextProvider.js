import { Application } from "@exabyte-io/ade.js";
import {
    ApplicationContextMixin,
    JSONSchemaFormDataProvider,
    MaterialContextMixin,
} from "@exabyte-io/code.js/dist/context";
import { Made } from "@exabyte-io/made.js";
import { mix } from "mixwith";

const defaultHubbardConfig = {
    atomicSpecies: "",
    atomicOrbital: "",
    hubbardUValue: 0.01,
};

export class HubbardContextProvider extends mix(JSONSchemaFormDataProvider).with(
    ApplicationContextMixin,
    MaterialContextMixin,
) {
    static Application = Application;

    static Material = Made.Material;

    constructor(config) {
        super(config);
        this.uniqueElements = this.material.Basis.uniqueElements;
    }

    // eslint-disable-next-line class-methods-use-this
    get uiSchema() {
        return {
            atomicSpecies: {},
            atomicOrbital: {},
            hubbardUValue: {},
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get defaultData() {
        return defaultHubbardConfig;
    }

    // eslint-disable-next-line class-methods-use-this
    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U (Quantum Espresso) calculation.",
            type: "object",
            properties: {
                atomicSpecies: {
                    type: "string",
                    enum: this.uniqueElements,
                    default: defaultHubbardConfig.atomicSpecies,
                },
                atomicOrbital: {
                    type: "string",
                    enum: ["2p", "3s", "3p"],
                    default: defaultHubbardConfig.atomicOrbital,
                },
                hubbardUValue: {
                    type: "number",
                    default: defaultHubbardConfig.hubbardUValue,
                },
            },
        };
    }
}
