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
            "ui:options": {
                addable: true,
                orderable: false,
                removable: true,
            },
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
            type: "array",
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Atomic species",
                        enum: this.uniqueElements,
                        default: defaultHubbardConfig.atomicSpecies,
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Atomic orbital (Hubbard manifold)",
                        enum: ["2p", "3s", "3p"],
                        default: defaultHubbardConfig.atomicOrbital,
                    },
                    hubbardUValue: {
                        type: "number",
                        title: "Hubbard U value",
                        default: defaultHubbardConfig.hubbardUValue,
                    },
                },
            },
        };
    }
}
