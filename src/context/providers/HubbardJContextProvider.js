import { HubbardUContextProvider } from "./HubbardUContextProvider";

const defaultHubbardConfig = {
    paramType: "U",
    atomicSpecies: "",
    atomicOrbital: "2p",
    value: 1.0,
};

export class HubbardJContextProvider extends HubbardUContextProvider {
    get defaultData() {
        return [
            {
                ...defaultHubbardConfig,
                atomicSpecies: this.uniqueElements?.length > 0 ? this.uniqueElements[0] : "",
            },
        ];
    }

    get uiSchemaStyled() {
        return {
            "ui:options": {
                addable: true,
                orderable: true,
                removable: true,
            },
            title: {
                "ui:classNames": "col-xs-12",
            },
            items: {
                paramType: this.defaultFieldStyles,
                atomicSpecies: this.defaultFieldStyles,
                atomicOrbital: this.defaultFieldStyles,
                value: this.defaultFieldStyles,
            },
        };
    }

    get jsonSchema() {
        return {
            $schema: "http://json-schema.org/draft-04/schema#",
            title: "",
            description: "Hubbard parameters for DFT+U+J calculation.",
            type: "array",
            items: {
                type: "object",
                properties: {
                    paramType: {
                        type: "string",
                        title: "Species",
                        enum: ["U", "J", "B", "E2", "E3"],
                        default: defaultHubbardConfig.paramType,
                    },
                    atomicSpecies: {
                        type: "string",
                        title: "Species",
                        enum: this.uniqueElements,
                        default: this.uniqueElements?.length > 0 ? this.uniqueElements[0] : "",
                    },
                    atomicOrbital: {
                        type: "string",
                        title: "Orbital",
                        enum: [
                            "2p",
                            "3s",
                            "3p",
                            "3d",
                            "4s",
                            "4p",
                            "4d",
                            "4f",
                            "5s",
                            "5p",
                            "5d",
                            "5f",
                            "6s",
                            "6p",
                            "6d",
                            "7s",
                            "7p",
                            "7d",
                        ],
                        default: defaultHubbardConfig.atomicOrbital,
                    },
                    value: {
                        type: "number",
                        title: "Value (eV)",
                        default: defaultHubbardConfig.value,
                    },
                },
            },
            minItems: 1,
        };
    }
}
