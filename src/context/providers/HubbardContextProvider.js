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
    atomicOrbital: "2p",
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

    get defaultData() {
        return [{ ...defaultHubbardConfig, atomicSpecies: this.uniqueElements }];
    }

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: false,
                removable: true,
            },
            title: {
                classNames: "col-xs-12",
            },
            items: {
                atomicSpecies: this.defaultFieldStyles,
                atomicOrbital: this.defaultFieldStyles,
                hubbardUValue: this.defaultFieldStyles,
            },
        };
    }

    // eslint-disable-next-line class-methods-use-this
    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U calculation.",
            type: "array",
            items: {
                type: "object",
                properties: {
                    atomicSpecies: {
                        type: "string",
                        title: "Atomic species",
                        enum: this.uniqueElements,
                        default: this.uniqueElements[0],
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Atomic orbital",
                        enum: ["2p", "3s", "3p", "3d", "4s", "4p", "4d"],
                        default: defaultHubbardConfig.atomicOrbital,
                    },
                    hubbardUValue: {
                        type: "number",
                        title: "Hubbard U value",
                        default: defaultHubbardConfig.hubbardUValue,
                    },
                },
            },
            minItems: 1,
        };
    }
}
